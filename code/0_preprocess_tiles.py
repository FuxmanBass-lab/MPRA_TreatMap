#!/usr/bin/env python3
"""
Build a wide activity matrix for MPRA heatmap from control activity files,
multiple treatment-effect files, and a tile-id strand map.

Inputs
------
1) Control activity files (TSV/CSV) in a directory (e.g. data/control_activity/)
   Expected columns (case-insensitive, best-effort autodetect):
     - ID column: one of {IID, ID, tile_id}
     - log2FC:    one of {log2FoldChange, log2fc, log2FC}
     - padj:      one of {padj, qvalue, q_value}
     - coverage:   require ctrl_mean >= 50 (configurable via --dna-threshold)
   Thresholds for **controls** (AND): padj < 0.05, log2FC > 2.0
   If multiple files correspond to the same base label (e.g. DMSO1/2/3),
   they are collapsed to a single column by taking the maximum log2FC
   among files that pass thresholds (done separately for sense/antisense).

2) Treatment effect files (TSV/CSV) in a directory (e.g. data/treatments_effect/)
   Expected columns (case-insensitive):
     - ID column: {ID, IID, tile_id}
     - log2FC:    {log2FoldChange, log2fc, log2FC}
     - padj:      {padj, qvalue, q_value}
     - coverage:   require dna_mean  >= 50 (configurable via --dna-threshold)
   Thresholds for **treatments**:
     - primary keep: padj < 0.25 and |log2FC| >= 0.5
     - complement rule: if a tile passes primary keep, also include its
       complementary strand **if** its padj < 0.25 (log2FC can be anything)
     - replicate activity gate (per treatment file):
         • compute mean of control replicates and mean of treatment replicates from columns like B2_DMSO_r1_activity, ...
         • take max(mean_control, mean_treatment) and require it > 2 (rows failing this are discarded)
         • Choose control label as the first present in KNOWN_CONTROLS (canonical, case‑insensitive; aliases supported: heat→HEATCTRL, water→H2O, ctrl/control→CONTROL).
   If multiple files exist for one treatment label (e.g. DOX vs DMSO1, DOX vs DMSO2),
   per tile we keep the record with **largest |log2FC|** (tie -> smaller padj).

3) Strand map CSV (required): maps sense ↔ antisense and provides genomic bp.
   Required columns (case-sensitive):
     - sense, antisense, bp

4) Virus name (optional): only used if --filter-by-virus is provided; otherwise the
   entire strand map is used as-is. (Assumes your map is already virus-specific.)

Output
------
A wide CSV with columns:
  bp, log2fc_<LABEL>, log2fc_<LABEL>_antisense, ...
where <LABEL> include control labels (e.g. DMSO, H2O, HeatCtrl) and
all treatment labels (e.g. DOX, H2O2, Dex, ABT, ...).
Treatment labels are preserved exactly as derived from the filenames (no case/style changes).

By default, one row is written for **every tile/bp present in the strand map** (even if all values are empty). Use --drop-empty-rows to only write rows that have at least one non-empty value in any control/treatment column.

Usage
-----
python 0_preprocess_tiles.py \
  --controls-dir ../data/control_activity \
  --treatments-dir ../data/treatments_effect \
  --map ../data/strands/ebv_strands_match.csv \
  --virus EBV \
  --out ../results/EBV/ebv_act_R.csv \
  --verbose

# Optional: if your strands CSV includes multiple viruses, you can subset via:
#   --virus EBV --filter-by-virus

Optional thresholds / verbosity:
  --control-padj 0.05 --control-lfc 2.0 \
  --treat-padj 0.25 --treat-lfc 0.585 \
  --dna-threshold 50 \
  --verbose
"""
from __future__ import annotations
import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple, Set

# -------------------------
# Virus alias handling
# -------------------------
def _virus_patterns(q: str) -> list[str]:
    """Expand a user query like 'EBV' into a list of lowercase substrings to
    match in sense/antisense IDs. Allows common aliases.
    """
    if not q:
        return []
    ql = q.strip().lower()
    # built-in alias map (lowercased)
    ALIASES = {
        "ebv": ["epstein_barr_virus", "herpesvirus:epstein_barr_virus", "epstein-barr"],
        "hhv-8": ["kaposi_sarcoma_[hhv-8]", "herpesvirus:kaposi_sarcoma_[hhv-8]", "hhv8"],
        "adenovirus": ["adenovirus"],
    }
    pats = [ql]
    pats.extend(ALIASES.get(ql, []))
    # also split on commas so users can pass multiple tokens, e.g. "EBV, HHV-8"
    if "," in q:
        pats.extend([p.strip().lower() for p in q.split(",") if p.strip()])
    # dedupe while preserving order
    out = []
    for p in pats:
        if p and p not in out:
            out.append(p)
    return out

# -------------------------
# Small CSV/TSV utilities
# -------------------------

def _open(path: str, mode: str = "r"):
    if path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    return open(path, mode, encoding="utf-8", newline="")

def _detect_delimiter(path: str, default: str = ",") -> str:
    try:
        with _open(path, "r") as fh:
            sample = fh.read(4096)
            fh.seek(0)
            dialect = csv.Sniffer().sniff(sample, delimiters=[",", "\t", ";", "|"])
            return dialect.delimiter
    except Exception:
        ext = Path(path).suffix.lower()
        return "\t" if ext == ".tsv" else default

def read_table(path: str) -> Tuple[List[str], List[List[str]]]:
    delim = _detect_delimiter(path)
    with _open(path, "r") as fh:
        rdr = csv.reader(fh, delimiter=delim)
        rows = list(rdr)
    if not rows:
        raise ValueError(f"No rows in file: {path}")
    return rows[0], rows[1:]

def _norm(s: str) -> str:
    return s.strip().lower().replace(" ", "_")

def find_col(header: List[str], candidates: Iterable[str], fallback: Optional[int] = None) -> int:
    # exact
    for c in candidates:
        if c in header:
            return header.index(c)
    # case-insensitive normalized
    idx = {_norm(h): i for i, h in enumerate(header)}
    for c in candidates:
        i = idx.get(_norm(c))
        if i is not None:
            return i
    if fallback is not None:
        return fallback
    raise KeyError(f"None of the columns {list(candidates)} found in header: {header}")


def to_float(x: str) -> Optional[float]:
    try:
        if x is None:
            return None
        x = x.strip()
        if x == "" or x.upper() == "NA":
            return None
        return float(x)
    except Exception:
        return None

# -------------------------
# Activity replicate detection utilities (for treatment files)
# -------------------------
import re

# Flexible patterns:
#   * separators may be '_', '-', or '.'
#   * replicate token may be 'rN' or 'repN'
#   * trailing 'activity' (and anything after) is optional
#   * case-insensitive
ACTIVITY_PATTERNS = [
    re.compile(r"^.*?[_.-]([A-Za-z0-9.+-]+)[_.-](?:r|rep)\s*0*(\d+)(?:[_.-]activity.*)?$", re.IGNORECASE),
    # common variant without leading batch prefix
    re.compile(r"^([A-Za-z0-9.+-]+)[_.-](?:r|rep)\s*0*(\d+)(?:[_.-]activity.*)?$", re.IGNORECASE),
]

KNOWN_CONTROLS = ["DMSO","H2O","HeatCtrl","CONTROL"]

# Canonicalize condition tokens seen in replicate activity column names and filenames
# e.g., "heat" -> "HeatCtrl", "ctrl" -> "CONTROL", "water" -> "H2O"
COND_ALIASES = {
    "dmso": "DMSO",
    "vehicle": "DMSO",
    "h2o": "H2O",
    "water": "H2O",
    "h2o2": "H2O2",
    "heatctrl": "HeatCtrl",
    "heat": "Heat",
    "ctrl": "CONTROL",
    "control": "CONTROL",
}

def canon_cond(token: str) -> str:
    if token is None:
        return ""
    t = token.strip().lower()
    return COND_ALIASES.get(t, t.upper())

# Pretty label helper for output columns
def pretty_treatment_label(token: str) -> str:
    """Return a nice output label: controls stay uppercase; treatment names default to TitleCase.
    Chemistry-like tokens with digits (e.g., H2O2) stay as-is.
    """
    c = canon_cond(token)
    if c in KNOWN_CONTROLS:
        return c
    # keep tokens containing digits as-is (e.g., H2O2)
    if any(ch.isdigit() for ch in c):
        return c
    # TitleCase (first letter upper, rest lower), but preserve existing TitleCase like "Heat"
    return c[:1].upper() + c[1:].lower() if c else c

from typing import Dict, List, Optional

def _fallback_parse_activity(col: str) -> Optional[str]:
    """Very permissive fallback: split on separators and look for ... <COND> <r/repN> activity ..."""
    parts = re.split(r"[_.-]+", col)
    # scan triplets: cond, r/rep token, maybe 'activity'
    for j in range(len(parts) - 1):
        tok = parts[j]
        nxt = parts[j+1]
        if re.match(r"^(?:r|rep)\s*\d+$", nxt, flags=re.IGNORECASE):
            return tok
    return None

def detect_activity_replicates(header: List[str]) -> Dict[str, List[int]]:
    """Scan header for replicate activity columns; return {COND: [col_indices,...]} with COND canonicalized.
    Supports patterns like '..._COND_r1_activity', 'B2-Heat-rep2-activity', 'Heat.r3'.
    """
    cond2idx: Dict[str, List[int]] = {}
    for i, h in enumerate(header):
        cond_token = None
        for pat in ACTIVITY_PATTERNS:
            m = pat.match(h)
            if m:
                cond_token = m.group(1)
                break
        if cond_token is None:
            cond_token = _fallback_parse_activity(h)
        if not cond_token:
            continue
        cond = canon_cond(cond_token)
        cond2idx.setdefault(cond, []).append(i)
    return cond2idx

def mean_cols(row: List[str], idxs: List[int]) -> Optional[float]:
    if not idxs:
        return None
    vals: List[float] = []
    for i in idxs:
        if i < len(row):
            v = to_float(row[i])
            if v is not None:
                vals.append(v)
    if not vals:
        return None
    return sum(vals) / float(len(vals))

# -------------------------
# Label parsing
# -------------------------

def base_label_from_control_fname(fname: str) -> str:
    """Collapse replicate-like names into a base label (DMSO, H2O, HeatCtrl, H2O2)."""
    s = Path(fname).stem.upper()
    if "H2O2" in s:
        return "H2O2"
    if "H2O" in s:
        return "H2O"
    if "DMSO" in s:
        return "DMSO"
    if "HEATCTRL" in s:
        return "HeatCtrl"
    m = re.search(r"([A-Z]+)[^A-Z]*_ACTIVITY", s)
    return m.group(1).title() if m else s.title()

def label_from_treatment_fname(fname: str) -> str:
    s = Path(fname).stem
    m = re.search(r"comparison_([^_]+)_vs_", s, re.IGNORECASE)
    if m:
        return m.group(1)
    parts = re.split(r"[_-]", s)
    return parts[0]

# -------------------------
# Core processing
# -------------------------

def load_map(map_csv: str, virus: Optional[str] = None, filter_by_virus: bool = False):
    hdr, rows = read_table(map_csv)
    sense_i = find_col(hdr, ["sense"])
    anti_i  = find_col(hdr, ["antisense"])
    bp_i    = find_col(hdr, ["bp"])

    sense2anti: Dict[str, str] = {}
    anti2sense: Dict[str, str] = {}
    id2bp: Dict[str, int] = {}

    kept = 0
    patterns = _virus_patterns(virus) if (filter_by_virus and virus) else []
    for r in rows:
        if max(sense_i, anti_i, bp_i) >= len(r):
            continue
        s_id = r[sense_i].strip()
        a_id = r[anti_i].strip()
        bp_v = r[bp_i].strip()
        if not s_id or not a_id or not bp_v:
            continue
        if patterns:
            s_low = s_id.lower(); a_low = a_id.lower()
            if not any((p in s_low) or (p in a_low) for p in patterns):
                continue
        try:
            bp = int(float(bp_v))
        except Exception:
            continue
        sense2anti[s_id] = a_id
        anti2sense[a_id] = s_id
        id2bp[s_id] = bp
        id2bp[a_id] = bp
        kept += 1

    if kept == 0:
        raise RuntimeError(f"No rows found in map {map_csv} (after optional virus filter)")
    return sense2anti, anti2sense, id2bp

def merge_controls(controls_dir: str, sense2anti, anti2sense, id2bp,
                   padj_thr: float, lfc_thr: float, dna_thr: float, verbose: bool=False):
    """Return ((sense_by_label, anti_by_label), labels) where inner dicts are {label: {bp: log2fc}}."""
    sense_by_label: Dict[str, Dict[int, float]] = {}
    anti_by_label:  Dict[str, Dict[int, float]] = {}
    labels: List[str] = []

    files = sorted([str(p) for p in Path(controls_dir).glob("*.tsv")] +
                   [str(p) for p in Path(controls_dir).glob("*.csv")])
    if verbose:
        print(f"[controls] found {len(files)} files in {controls_dir}", file=sys.stderr)

    for f in files:
        try:
            hdr, rows = read_table(f)
        except Exception as e:
            print(f"[controls] skip {f}: {e}", file=sys.stderr)
            continue
        try:
            id_i  = find_col(hdr, ["IID","ID","tile_id"], fallback=0)
            lfc_i = find_col(hdr, ["log2FoldChange","log2fc","log2FC"])
            pad_i = find_col(hdr, ["padj","qvalue","q_value"])
            try:
                cov_i = find_col(hdr, ["ctrl_mean", "control_mean", "ctrl", "control"])
            except Exception:
                cov_i = None
                if verbose:
                    print(f"[controls] {f}: WARNING: no 'ctrl_mean' column found; skipping DNA coverage filter", file=sys.stderr)
        except Exception as e:
            print(f"[controls] {f}: {e}", file=sys.stderr)
            continue

        label = base_label_from_control_fname(f)
        if label not in sense_by_label:
            sense_by_label[label] = {}
            anti_by_label[label]  = {}
            labels.append(label)
        sense_acc = sense_by_label[label]
        anti_acc  = anti_by_label[label]

        kept_rows = 0
        for r in rows:
            if max(id_i, lfc_i, pad_i) >= len(r):
                continue
            tid  = r[id_i].strip()
            if tid not in id2bp:
                continue
            # Coverage filter
            cov = to_float(r[cov_i]) if cov_i is not None and cov_i < len(r) else None
            if cov is None or cov < dna_thr:
                continue
            lfc  = to_float(r[lfc_i])
            padj = to_float(r[pad_i])
            if lfc is None or padj is None:
                continue
            if not (padj < padj_thr and lfc > lfc_thr):
                continue
            bp = id2bp[tid]
            if tid in sense2anti:
                old = sense_acc.get(bp)
                if old is None or lfc > old:
                    sense_acc[bp] = lfc
            elif tid in anti2sense:
                old = anti_acc.get(bp)
                if old is None or lfc > old:
                    anti_acc[bp] = lfc
            kept_rows += 1
        if verbose:
            print(f"[controls] {label}: kept {kept_rows}", file=sys.stderr)

    return (sense_by_label, anti_by_label), labels

def best_of_records(a: Tuple[float,float], b: Tuple[float,float]) -> Tuple[float,float]:
    # choose by max |lfc|, tie-breaker min padj
    lfc_a, pad_a = a; lfc_b, pad_b = b
    abs_a, abs_b = abs(lfc_a), abs(lfc_b)
    if abs_b > abs_a:
        return b
    if abs_b < abs_a:
        return a
    if pad_b is not None and pad_a is not None and pad_b < pad_a:
        return b
    return a

def merge_treatments(treat_dir: str, sense2anti, anti2sense, id2bp,
                     padj_thr: float, lfc_thr: float, dna_thr: float, verbose: bool=False):
    # per-label best records
    label_best: Dict[str, Dict[str, Tuple[float, float]]] = {}

    files = sorted([str(p) for p in Path(treat_dir).glob("*.tsv")] +
                   [str(p) for p in Path(treat_dir).glob("*.csv")])
    if verbose:
        print(f"[treat] found {len(files)} files in {treat_dir}", file=sys.stderr)

    for f in files:
        try:
            hdr, rows = read_table(f)
        except Exception as e:
            print(f"[treat] skip {f}: {e}", file=sys.stderr)
            continue
        try:
            id_i  = find_col(hdr, ["ID","IID","tile_id"], fallback=0)
            lfc_i = find_col(hdr, ["log2FoldChange","log2fc","log2FC"])
            pad_i = find_col(hdr, ["padj","qvalue","q_value"])
            try:
                dna_i = find_col(hdr, ["dna_mean", "dna", "dna_coverage", "dna_cov"])  # primary is 'dna_mean'
            except Exception:
                dna_i = None
                if verbose:
                    print(f"[treat] {f}: WARNING: no 'dna_mean' column found; skipping DNA coverage filter", file=sys.stderr)
        except Exception as e:
            print(f"[treat] {f}: {e}", file=sys.stderr)
            continue

        raw_label = label_from_treatment_fname(f)
        # Detect replicate activity columns and pick control/treatment sets for this file
        cond2idx = detect_activity_replicates(hdr)  # keys are canonical already
        # --- VERBOSE DEBUGGING OF REPLICATE DETECTION ---
        if verbose:
            if cond2idx:
                summary = ", ".join(f"{k}:{len(v)}" for k, v in sorted(cond2idx.items()))
                print(f"[treat] {Path(f).name}: detected replicate conditions -> {summary}", file=sys.stderr)
            else:
                print(f"[treat] {Path(f).name}: no replicate-like activity columns detected in header", file=sys.stderr)
        # --- END DEBUGGING BLOCK ---
        ctrl_label = None
        if cond2idx:
            # choose control label by priority among known controls present
            for c in KNOWN_CONTROLS:
                if c in cond2idx:
                    ctrl_label = c
                    break
            # canonicalize the label from filename for matching
            label_c = canon_cond(raw_label)
            # if labeled treatment is not found, but there are exactly two conditions, pick the other as treatment
            if label_c not in cond2idx and len(cond2idx) == 2:
                other = [k for k in cond2idx.keys() if k != (ctrl_label or "")]
                treat_for_means = other[0] if other else label_c
            else:
                treat_for_means = label_c
            ctrl_cols  = cond2idx.get(ctrl_label or "", [])
            treat_cols = cond2idx.get(treat_for_means, [])
            gate_defined = bool(ctrl_cols or treat_cols)
            # Preserve the original label as taken from the filename
            label = raw_label
        else:
            ctrl_cols, treat_cols = [], []
            gate_defined = False
            # Preserve the original label as taken from the filename
            label = raw_label
        ctrl_name = ctrl_label if cond2idx else None
        treat_name = treat_for_means if cond2idx else None

        d = label_best.setdefault(label, {})
        scanned = 0
        passed_activity = 0
        for r in rows:
            if max(id_i, lfc_i, pad_i) >= len(r):
                continue
            tid  = r[id_i].strip()
            if tid not in id2bp:
                continue
            # DNA coverage filter
            if dna_i is not None:
                dna_v = to_float(r[dna_i]) if dna_i < len(r) else None
                if dna_v is None or dna_v < dna_thr:
                    continue
            else:
                # if column missing, we accept the row (already warned above)
                pass
            # Replicate activity gate: require max(mean_ctrl, mean_treat) > 2 (only if gate is defined)
            if gate_defined:
                mc = mean_cols(r, ctrl_cols) if ctrl_cols else None
                mt = mean_cols(r, treat_cols) if treat_cols else None
                # Require at least one of the means; if you want BOTH, uncomment next line and remove the line after
                # if mc is None or mt is None: continue
                candidates = [x for x in (mc, mt) if x is not None]
                if not candidates or max(candidates) <= 2.0:
                    continue
                passed_activity += 1
            lfc  = to_float(r[lfc_i])
            padj = to_float(r[pad_i])
            if lfc is None or padj is None:
                continue
            cur = (lfc, padj)
            prev = d.get(tid)
            d[tid] = cur if prev is None else best_of_records(prev, cur)
            scanned += 1
        if verbose:
            extra = f" (replicates: ctrl={ctrl_name}, treat={treat_name})" if (ctrl_name or treat_name) else ""
            if 'gate_defined' in locals() and gate_defined:
                print(f"[treat] {label}: scanned {scanned} rows; activity>2: {passed_activity}{extra}", file=sys.stderr)
            else:
                print(f"[treat] {label}: scanned {scanned} rows; no replicate gate{extra}", file=sys.stderr)

    # Apply thresholds + complement rule, then split to sense/antisense dicts keyed by bp
    sense_out: Dict[str, Dict[int, float]] = {}
    anti_out:  Dict[str, Dict[int, float]] = {}
    labels: List[str] = []
    hits_per_id: Dict[str, Set[str]] = {}

    for label, recs in label_best.items():
        primary_pass = {
            tid for tid, (lfc, padj) in recs.items()
            if padj is not None and padj < padj_thr and abs(lfc) >= lfc_thr
        }
        include = set(primary_pass)
        # complement rule
        for tid in list(primary_pass):
            comp = sense2anti.get(tid) or anti2sense.get(tid)
            if comp and comp in recs:
                lfc_c, pad_c = recs[comp]
                if pad_c is not None and pad_c < padj_thr:
                    include.add(comp)
        # record per-ID treatment hits
        for tid in include:
            s = hits_per_id.setdefault(tid, set())
            s.add(label)
        s_map: Dict[int, float] = {}
        a_map: Dict[int, float] = {}
        for tid in include:
            lfc, _ = recs[tid]
            bp = id2bp.get(tid)
            if bp is None:
                continue
            if tid in sense2anti:
                s_map[bp] = lfc
            elif tid in anti2sense:
                a_map[bp] = lfc
        sense_out[label] = s_map
        anti_out[label]  = a_map
        labels.append(label)

    return sense_out, anti_out, labels, hits_per_id

# -------------------------
# Assembly & write
# -------------------------

def assemble_and_write(out_csv: str,
                        controls, c_labels: List[str],
                        treats_s: Dict[str, Dict[int, float]],
                        treats_a: Dict[str, Dict[int, float]],
                        t_labels: List[str],
                        id2bp: Dict[str, int], verbose: bool=False, drop_empty: bool=False):
    sense_ctrl, anti_ctrl = controls

    # All bps that appear in the map (keeps full axis)
    all_bps = sorted({bp for bp in id2bp.values()})

    # Column order: controls (alpha) then treatments (alpha)
    c_labels = sorted(set(c_labels))
    t_labels = sorted(set(t_labels))

    cols: List[str] = ["bp"]
    for lab in c_labels:
        cols += [f"log2fc_{lab}", f"log2fc_{lab}_antisense"]
    for lab in t_labels:
        cols += [f"log2fc_{lab}", f"log2fc_{lab}_antisense"]

    with _open(out_csv, "w") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        written = 0
        for bp in all_bps:
            row_vals: List[Optional[float]] = []
            row = [bp]
            for lab in c_labels:
                v1 = sense_ctrl.get(lab, {}).get(bp, "")
                v2 = anti_ctrl.get(lab, {}).get(bp, "")
                row.extend([v1, v2])
                row_vals.extend([v1, v2])
            for lab in t_labels:
                v1 = treats_s.get(lab, {}).get(bp, "")
                v2 = treats_a.get(lab, {}).get(bp, "")
                row.extend([v1, v2])
                row_vals.extend([v1, v2])
            if drop_empty:
                # consider empty if all are '' or None
                any_val = any((v is not None and v != "") for v in row_vals)
                if not any_val:
                    continue
            w.writerow(row)
            written += 1
    if verbose:
        print(f"[out] wrote {out_csv} with {written} rows (of {len(all_bps)} bp in map) and {len(cols)} columns", file=sys.stderr)

def write_treatment_hits(out_csv: str, hits: Dict[str, Set[str]]):
    if not out_csv:
        return
    rows = []
    for tid in sorted(hits.keys()):
        labels = sorted(hits[tid])
        rows.append((tid, ",".join(labels)))
    with _open(out_csv, "w") as fh:
        w = csv.writer(fh)
        w.writerow(["ID", "treatments"])  # header
        for tid, labs in rows:
            w.writerow([tid, labs])

# -------------------------
# CLI
# -------------------------

def parse_args(argv: Optional[List[str]] = None):
    p = argparse.ArgumentParser(description="Preprocess control/treatment activity into a wide log2FC matrix per bp")
    p.add_argument("--controls-dir", required=True, help="Directory of control activity TSV/CSV files")
    p.add_argument("--treatments-dir", required=True, help="Directory of treatment-effect TSV/CSV files")
    p.add_argument("--map", required=True, help="CSV mapping sense/antisense to bp (columns: sense,antisense,bp)")
    p.add_argument("--virus", required=False, default=None, help="(Optional) Virus name substring used only if --filter-by-virus is set")
    p.add_argument("--filter-by-virus", action="store_true", help="If set, filter the map by --virus; otherwise include all rows of the map")
    p.add_argument("--out", required=True, help="Output CSV path")
    p.add_argument("--treat-hits-out", default="", help="Optional CSV listing tile IDs with treatment hits (columns: ID,treatments)")
    p.add_argument("--control-padj", type=float, default=0.05)
    p.add_argument("--control-lfc", type=float, default=2.0)
    p.add_argument("--treat-padj", type=float, default=0.25)
    p.add_argument("--treat-lfc", type=float, default=0.5)
    p.add_argument("--dna-threshold", type=float, default=50.0, help="Minimum DNA coverage (ctrl_mean for controls, dna_mean for treatments)")
    p.add_argument("--verbose", action="store_true")
    p.add_argument("--drop-empty-rows", action="store_true", help="If set, omit rows where all control/treatment values are empty")
    return p.parse_args(argv)

def main(argv: Optional[List[str]] = None) -> int:
    a = parse_args(argv)
    if not Path(a.controls_dir).is_dir():
        print(f"ERROR: controls-dir not found: {a.controls_dir}", file=sys.stderr)
        return 2
    if not Path(a.treatments_dir).is_dir():
        print(f"ERROR: treatments-dir not found: {a.treatments_dir}", file=sys.stderr)
        return 2
    if not Path(a.map).exists():
        print(f"ERROR: map CSV not found: {a.map}", file=sys.stderr)
        return 2

    sense2anti, anti2sense, id2bp = load_map(a.map, a.virus, a.filter_by_virus)

    (ctrl_s, ctrl_a), c_labels = merge_controls(
        a.controls_dir, sense2anti, anti2sense, id2bp,
        padj_thr=a.control_padj, lfc_thr=a.control_lfc, dna_thr=a.dna_threshold, verbose=a.verbose
    )

    tr_s, tr_a, t_labels, hits_per_id = merge_treatments(
        a.treatments_dir, sense2anti, anti2sense, id2bp,
        padj_thr=a.treat_padj, lfc_thr=a.treat_lfc, dna_thr=a.dna_threshold, verbose=a.verbose
    )

    assemble_and_write(a.out, (ctrl_s, ctrl_a), c_labels, tr_s, tr_a, t_labels, id2bp, verbose=a.verbose, drop_empty=a.drop_empty_rows)
    if a.treat_hits_out:
        write_treatment_hits(a.treat_hits_out, hits_per_id)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
