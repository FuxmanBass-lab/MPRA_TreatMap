#!/usr/bin/env python3
"""
Pull sequences for a set of tile IDs.

Takes two tabular files:
  1) IDs file (e.g., Significant_tile_ids.csv) whose first column contains IDs
     (your example calls it "Strand").
  2) Sequences file (e.g., seq.tsv) with columns like: tile_id, Sequence

Outputs a CSV with columns: ID,Sequence for the IDs that are present in both.

Examples
--------
# CSV of IDs + TSV of sequences, write to file
python 0_pull_seq.py \
  --ids Significant_tile_ids.csv \
  --seqs seq.tsv \
  --out significant_id_sequences.csv

# Explicit column overrides if your headers are different
python 0_pull_seq.py \
  --ids Significant_tile_ids.csv --ids-col Strand \
  --seqs seq.tsv --seqs-id-col tile_id --seq-col Sequence \
  --out significant_id_sequences.csv
"""
from __future__ import annotations
import argparse
import csv
import sys
from pathlib import Path
from typing import Iterable, List, Dict, Tuple, Optional

# ------------------------------
# Helpers
# ------------------------------

def _open(path: str):
    """Open plain text file; supports '-' for stdin/stdout context (read only)."""
    if path == "-":
        return sys.stdin
    return open(path, "r", encoding="utf-8", newline="")


def _detect_delimiter(path: str, default: str = ",") -> str:
    """Heuristic delimiter detection using csv.Sniffer, fallback by extension."""
    try:
        with _open(path) as fh:
            sample = fh.read(2048)
            fh.seek(0)
            sniffer = csv.Sniffer()
            dialect = sniffer.sniff(sample, delimiters=[",", "\t", ";", "|"])
            return dialect.delimiter
    except Exception:
        # Fallback by extension
        ext = Path(path).suffix.lower()
        if ext == ".tsv":
            return "\t"
        return default


def _read_tabular(path: str, *, force_delim: Optional[str] = None) -> Tuple[List[str], List[List[str]]]:
    delim = force_delim or _detect_delimiter(path)
    with _open(path) as fh:
        rdr = csv.reader(fh, delimiter=delim)
        rows = list(rdr)
    if not rows:
        raise ValueError(f"No rows found in file: {path}")
    header = rows[0]
    data = rows[1:]
    return header, data


def _norm(s: str) -> str:
    return s.strip().lower().replace(" ", "_")


def _index_col(header: List[str], preferred: Optional[str], fallbacks: Iterable[str]) -> int:
    """Find a column index by preferred name (exact or case-insensitive),
    then by any fallback names; finally fall back to the first column.
    """
    # Exact match first
    if preferred and preferred in header:
        return header.index(preferred)
    # Case-insensitive match
    if preferred is not None:
        ln = _norm(preferred)
        for i, h in enumerate(header):
            if _norm(h) == ln:
                return i
    # Try fallbacks
    for fb in fallbacks:
        for i, h in enumerate(header):
            if _norm(h) == _norm(fb):
                return i
    # Default: first column
    return 0


# ------------------------------
# Core logic
# ------------------------------

def collect_ids(ids_path: str, ids_col: Optional[str]) -> List[str]:
    header, data = _read_tabular(ids_path)
    col_idx = _index_col(header, ids_col, fallbacks=("id", "tile_id", "strand"))
    ids: List[str] = []
    seen = set()
    for row in data:
        if col_idx >= len(row):
            continue
        tid = row[col_idx].strip()
        if tid and tid not in seen:
            ids.append(tid)
            seen.add(tid)
    return ids


def load_sequences(seqs_path: str, id_col: Optional[str], seq_col: Optional[str]) -> Dict[str, str]:
    header, data = _read_tabular(seqs_path)
    id_idx = _index_col(header, id_col, fallbacks=("tile_id", "id"))
    seq_idx = _index_col(header, seq_col, fallbacks=("sequence", "seq"))
    m: Dict[str, str] = {}
    for row in data:
        if id_idx >= len(row) or seq_idx >= len(row):
            continue
        key = row[id_idx].strip()
        val = row[seq_idx].strip()
        if key:
            m[key] = val
    return m


def write_output(pairs: List[Tuple[str, str]], out_path: str):
    if out_path == "-":
        wfh = sys.stdout
        close_later = False
    else:
        wfh = open(out_path, "w", encoding="utf-8", newline="")
        close_later = True
    try:
        w = csv.writer(wfh)
        w.writerow(["ID", "Sequence"])  # header
        for tid, seq in pairs:
            w.writerow([tid, seq])
    finally:
        if close_later:
            wfh.close()


# ------------------------------
# CLI
# ------------------------------

def parse_args(argv: Optional[List[str]] = None):
    p = argparse.ArgumentParser(description="Pull sequences for a list of tile IDs")
    p.add_argument("--ids", required=True, help="CSV/TSV containing IDs (first column by default)")
    p.add_argument("--seqs", required=True, help="CSV/TSV containing sequences (tile_id & Sequence columns)")
    p.add_argument("--out", default="-", help="Output CSV path (default: stdout)")
    # Optional explicit column names
    p.add_argument("--ids-col", default=None, help="Column in --ids that contains the IDs (default: first col; recognizes 'Strand')")
    p.add_argument("--seqs-id-col", default=None, help="ID column in --seqs (default: tile_id)")
    p.add_argument("--seq-col", default=None, help="Sequence column in --seqs (default: Sequence)")
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    # Collect IDs (preserve order, de-dupe)
    ids = collect_ids(args.ids, args.ids_col)
    if not ids:
        print(f"ERROR: No IDs found in {args.ids}", file=sys.stderr)
        return 2

    # Load sequences map
    seq_map = load_sequences(args.seqs, args.seqs_id_col, args.seq_col)
    if not seq_map:
        print(f"ERROR: No sequences found in {args.seqs}", file=sys.stderr)
        return 2

    # Join in the order of the ids file
    pairs: List[Tuple[str, str]] = []
    missing: List[str] = []
    for tid in ids:
        s = seq_map.get(tid)
        if s is None:
            missing.append(tid)
        else:
            pairs.append((tid, s))

    # Report missing matches to stderr (but still write matched output)
    if missing:
        print(f"WARNING: {len(missing)} IDs not found in sequences file (first 5 shown): {missing[:5]}", file=sys.stderr)

    write_output(pairs, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
