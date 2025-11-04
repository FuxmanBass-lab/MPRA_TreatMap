#!/usr/bin/env python3
"""
Generate three FASTA files per treatment condition by aggregating tiles across all viruses.

Classes for a given run (requires --control and --treatment):
  1) control: tiles that are ACTIVE in the control column but have NO activity in the treatment column (both treatment strands NaN)
  2) upregulated ("up"): tiles with positive log2fc in the treatment (chosen strand rule below)
  3) downregulated ("down"): tiles with negative log2fc in the treatment (chosen strand rule below)

Rules:
  - For each bp within a virus, choose *one* strand only (sense OR antisense) per treatment.
    * If both strands have values, choose the strand with larger absolute log2fc; if tie, prefer sense.
    * If only one strand has a value, use that strand.
    * If neither strand has a value, classify as control and take the sense strand by default.
  - Do not duplicate tiles across classes for a given treatment.
  - Output files are grouped across all viruses (one set of three FASTAs per treatment).

Input layout (assumed from your repo structure):
  - Strands map per virus in:   data/strands/*_strands_match.csv
      Columns: sense, number, antisense, bp
  - Sequences for all tiles in: data/all_tiles_seqs.tsv
      Columns: tile_id, Sequence
  - Activity per virus in:      results/<VIRUS>/*_act_R.csv
      Columns include: bp, log2fc_<TREAT>, log2fc_<TREAT>_antisense, ...

Usage example:
  python code/03_fasta_for_jammer.py \
    --control Kwater \
    --treatment TNF \
    --strands-dir data/strands \
    --seq-tsv     data/all_tiles_seqs.tsv \
    --results-dir results \
    --outdir      downstream/fasta

      python code/03_fasta_for_jammer.py \
    --control Kwater \
    --treatment INFg \
    --strands-dir data/strands \
    --seq-tsv     data/all_tiles_seqs.tsv \
    --results-dir results \
    --outdir      downstream/fasta

      python code/03_fasta_for_jammer.py \
    --control Kwater \
    --treatment ILB1 \
    --strands-dir data/strands \
    --seq-tsv     data/all_tiles_seqs.tsv \
    --results-dir results \
    --outdir      downstream/fasta

      python code/03_fasta_for_jammer.py \
    --control Awater \
    --treatment INFa \
    --strands-dir data/strands \
    --seq-tsv     data/all_tiles_seqs.tsv \
    --results-dir results \
    --outdir      downstream/fasta

      python code/03_fasta_for_jammer.py \
    --control Awater \
    --treatment POLYDADT \
    --strands-dir data/strands \
    --seq-tsv     data/all_tiles_seqs.tsv \
    --results-dir results \
    --outdir      downstream/fasta

      python code/03_fasta_for_jammer.py \
    --control Awater \
    --treatment SENDAI \
    --strands-dir data/strands \
    --seq-tsv     data/all_tiles_seqs.tsv \
    --results-dir results \
    --outdir      downstream/fasta

"""

from __future__ import annotations
import argparse
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

import pandas as pd

# ------------------------------
# Helpers
# ------------------------------

def read_sequences(seq_tsv: Path) -> Dict[str, str]:
    """Read tile sequences table (TSV with columns: tile_id, Sequence).

    Returns dict mapping tile_id -> sequence (uppercased, no spaces).
    """
    df = pd.read_csv(seq_tsv, sep="\t").rename(columns=str.strip)
    if not {"tile_id", "Sequence"}.issubset(df.columns):
        raise ValueError(
            f"{seq_tsv} must have columns: tile_id, Sequence; has {list(df.columns)}"
        )
    seqs = (
        df[["tile_id", "Sequence"]]
        .dropna()
        .assign(Sequence=lambda d: d["Sequence"].astype(str).str.replace(r"\s+", "", regex=True).str.upper())
        .set_index("tile_id")["Sequence"]
        .to_dict()
    )
    return seqs


def read_strands_map(csv_path: Path) -> pd.DataFrame:
    """Read strands mapping with columns (sense, number, antisense, bp).

    Returns a DataFrame indexed by bp with columns 'sense' and 'antisense'.
    """
    # auto-detect delimiter (comma or tab)
    df = pd.read_csv(csv_path, sep=None, engine="python")
    # strip surrounding whitespace and normalize for matching
    df = df.rename(columns=lambda c: str(c).strip())
    # Accept common variations and be case-insensitive
    cols = list(df.columns)
    cols_norm_map = {}
    for c in cols:
        norm = str(c).strip().lower().replace(" ", "").replace("_", "")
        cols_norm_map[norm] = c

    def resolve_with_aliases(preferred: str, aliases: list[str]) -> str:
        keys = [preferred, *aliases]
        for k in keys:
            k_norm = k.strip().lower().replace(" ", "").replace("_", "")
            if k_norm in cols_norm_map:
                return cols_norm_map[k_norm]
        # final attempt: look for startswith matches (e.g., "bpindex")
        for norm, orig in cols_norm_map.items():
            if norm.startswith(preferred.strip().lower().replace(" ", "").replace("_", "")):
                return orig
        raise ValueError(f"{csv_path} missing required column '{preferred}' (tried aliases: {aliases})")

    bp_col = resolve_with_aliases(
        "bp",
        aliases=["basepair", "basepairs", "base_pair", "position", "pos", "bpindex", "bp_id", "index", "tile_bp", "bpno", "bpnum"],
    )
    sense_col = resolve_with_aliases(
        "sense",
        aliases=["sense_id", "tile_sense", "plus", "plus_id", "forward", "forward_id", "strand_plus", "sense_tile"],
    )
    antisense_col = resolve_with_aliases(
        "antisense",
        aliases=["antisense_id", "anti_sense", "minus", "minus_id", "reverse", "reverse_id", "strand_minus", "antisense_tile"],
    )

    sub = df[[bp_col, sense_col, antisense_col]].rename(
        columns={bp_col: "bp", sense_col: "sense", antisense_col: "antisense"}
    )
    # ensure 'bp' numeric for merging
    sub["bp"] = pd.to_numeric(sub["bp"], errors="coerce")
    sub = sub.dropna(subset=["bp"]).astype({"bp": int})
    return sub


def discover_activity_files(results_dir: Path) -> List[Path]:
    """Find per-virus *_act_R.csv files under results/<virus>/.
    Returns a list of paths.
    """
    paths: List[Path] = []
    for vdir in sorted([p for p in results_dir.iterdir() if p.is_dir()]):
        for p in vdir.glob("*_act_R.csv"):
            paths.append(p)
    if not paths:
        raise FileNotFoundError(f"No *_act_R.csv found under {results_dir}")
    return paths


def pick_strand_and_class(
    s_val: Optional[float], a_val: Optional[float]
) -> Tuple[str, Optional[str]]:
    """Decide which strand to take and which class ('up'/'down'/None).

    Returns (strand, cls) where strand in {'sense','antisense'} and
    cls in {'up','down',None}. If cls is None, it's *control* (no value).

    Tie-breaking: if both present with equal |val|, prefer sense.
    """
    s = None if pd.isna(s_val) else float(s_val)
    a = None if pd.isna(a_val) else float(a_val)

    if s is None and a is None:
        return ("sense", None)  # control, default to sense strand

    if s is not None and a is None:
        return ("sense", "up" if s > 0 else ("down" if s < 0 else None))

    if s is None and a is not None:
        return ("antisense", "up" if a > 0 else ("down" if a < 0 else None))

    # both present: choose larger absolute value
    if abs(s) > abs(a):
        return ("sense", "up" if s > 0 else ("down" if s < 0 else None))
    if abs(a) > abs(s):
        return ("antisense", "up" if a > 0 else ("down" if a < 0 else None))
    # tie
    return ("sense", "up" if s > 0 else ("down" if s < 0 else None))


def to_fasta_lines(ids: Iterable[str], seqs: Dict[str, str], missing: Set[str]) -> List[str]:
    lines: List[str] = []
    for tid in ids:
        seq = seqs.get(tid)
        if not seq:
            missing.add(tid)
            continue
        lines.append(f">{tid}\n{seq}\n")
    return lines

def choose_strand(vs: Optional[float], va: Optional[float]) -> Tuple[str, Optional[float]]:
    s = None if pd.isna(vs) else float(vs)
    a = None if pd.isna(va) else float(va)
    if s is None and a is None:
        return ("sense", None)
    if s is not None and a is None:
        return ("sense", s)
    if s is None and a is not None:
        return ("antisense", a)
    if abs(s) > abs(a):
        return ("sense", s)
    if abs(a) > abs(s):
        return ("antisense", a)
    return ("sense", s)


# ------------------------------
# Main logic
# ------------------------------

def build_fastas(
    strands_dir: Path, seq_tsv: Path, results_dir: Path, outdir: Path,
    control_name: str, treatment_name: str,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    seqs = read_sequences(seq_tsv)

    act_files = discover_activity_files(results_dir)

    per_virus: List[Tuple[str, pd.DataFrame, pd.DataFrame]] = []
    for act_path in act_files:
        virus = act_path.parent.name
        guess = (strands_dir / f"{virus.lower()}_strands_match.csv")
        if not guess.exists():
            guess = strands_dir / f"{virus}_strands_match.csv"
        if not guess.exists():
            raise FileNotFoundError(f"Strands map not found for {virus}: tried {virus.lower()}_strands_match.csv and {virus}_strands_match.csv")
        strands_df = read_strands_map(guess)
        act_df = pd.read_csv(act_path)
        if "bp" not in act_df.columns:
            raise ValueError(f"{act_path} missing 'bp' column")
        act_df["bp"] = pd.to_numeric(act_df["bp"], errors="coerce").astype("Int64")
        act_df = act_df.dropna(subset=["bp"]).astype({"bp": int})
        per_virus.append((virus, strands_df, act_df))

    agg: Dict[str, Set[str]] = {"control": set(), "up": set(), "down": set()}

    for virus, strands_df, act_df in per_virus:
        merged = strands_df.merge(act_df, on="bp", how="left", validate="1:1")
        s_ctrl = f"log2fc_{control_name}"
        a_ctrl = f"log2fc_{control_name}_antisense"
        s_treat = f"log2fc_{treatment_name}"
        a_treat = f"log2fc_{treatment_name}_antisense"
        # If a column is missing for a virus, treat as all-NaN
        ctrl_s_vals = merged[s_ctrl] if s_ctrl in merged.columns else pd.Series([pd.NA] * len(merged))
        ctrl_a_vals = merged[a_ctrl] if a_ctrl in merged.columns else pd.Series([pd.NA] * len(merged))
        trt_s_vals  = merged[s_treat] if s_treat in merged.columns else pd.Series([pd.NA] * len(merged))
        trt_a_vals  = merged[a_treat] if a_treat in merged.columns else pd.Series([pd.NA] * len(merged))

        for sense_id, antisense_id, vs_c, va_c, vs_t, va_t in zip(
            merged["sense"], merged["antisense"], ctrl_s_vals, ctrl_a_vals, trt_s_vals, trt_a_vals
        ):
            # Determine treatment class/strand first
            trt_strand, trt_val = choose_strand(vs_t, va_t)
            ctrl_strand, ctrl_val = choose_strand(vs_c, va_c)

            # 1) Control class: active in control (has value on either strand) AND treatment has no value on either strand
            if pd.isna(vs_t) and pd.isna(va_t) and (not pd.isna(vs_c) or not pd.isna(va_c)):
                chosen_id = sense_id if ctrl_strand == "sense" else antisense_id
                agg["control"].add(chosen_id)
                continue

            # 2) Treatment up/down classes when treatment has a value
            if trt_val is not None:
                chosen_id = sense_id if trt_strand == "sense" else antisense_id
                if trt_val > 0:
                    agg["up"].add(chosen_id)
                elif trt_val < 0:
                    agg["down"].add(chosen_id)
                # exactly zero -> ignore

    label = f"{treatment_name}_vs_{control_name}"
    for cls in ("control", "up", "down"):
        ids = sorted(agg[cls])
        missing: Set[str] = set()
        lines = to_fasta_lines(ids, seqs, missing)
        fname = {
            "control": f"{label}_control.fa",
            "up": f"{label}_upregulated.fa",
            "down": f"{label}_downregulated.fa",
        }[cls]
        out_path = outdir / fname
        with open(out_path, "w") as fh:
            fh.writelines(lines)
        print(f"[wrote] {out_path} (n={len(ids)}; missing_seqs={len(missing)})")
        if missing:
            miss_path = outdir / f"{label}_{cls}_MISSING_IDS.txt"
            with open(miss_path, "w") as fh:
                for tid in sorted(missing):
                    fh.write(f"{tid}\n")
            print(f"  [warn] {len(missing)} IDs missing from sequence table -> {miss_path}")


# ------------------------------
# CLI
# ------------------------------

def main():
    ap = argparse.ArgumentParser(description="Generate FASTAs for a specified control vs treatment across viruses (strand-aware)")
    ap.add_argument("--control", required=True, help="Control column base name (e.g., DMSO, H2O, HeatCtrl)")
    ap.add_argument("--treatment", required=True, help="Treatment column base name (e.g., TNF, Dex, INFG)")
    ap.add_argument("--strands-dir", type=Path, default=Path("data/strands"), help="Directory with *_strands_match.csv files")
    ap.add_argument("--seq-tsv", type=Path, default=Path("data/all_tiles_seqs.tsv"), help="TSV with columns: tile_id, Sequence")
    ap.add_argument("--results-dir", type=Path, default=Path("results"), help="Directory containing per-virus subfolders with *_act_R.csv")
    ap.add_argument("--outdir", type=Path, default=Path("results/fasta"), help="Where to write FASTAs")
    args = ap.parse_args()

    build_fastas(
        strands_dir=args.strands_dir,
        seq_tsv=args.seq_tsv,
        results_dir=args.results_dir,
        outdir=args.outdir,
        control_name=args.control,
        treatment_name=args.treatment,
    )

if __name__ == "__main__":
    main()