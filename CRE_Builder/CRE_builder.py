#!/usr/bin/env python3
"""
CRE_builder.py

Build cis-regulatory elements (CREs) from MPRA tile-level activities.

Definitions
-----------
- Tiles are 200 bp long, positioned by their middle point (e.g., first tile at 100 bp).
- Consecutive tiles are 50 bp apart (centers differ by 50).
- A tile spans [center - 100, center + 100].
- A CRE is formed by joining *active* tiles that are consecutive **by at least one strand**
  (i.e., consecutive tiles both have sense activity > 0 OR both have antisense
  activity > 0). Strand agreement can vary across joins; a run may switch
  strands only where the consecutive pair shares at least one active strand.
- The CRE activity for each strand is the maximum activity of that strand over the
  tiles included in the CRE.

Outputs
-------
- Per **condition** CSV combining all viruses: `<out_dir>/<condition>.csv`  (conditions include `<name>` and `<name>_down`)
  Columns: `virus,start,end,sense_activity,antisense_activity` (1-based inclusive coords)
- Per **virus & condition** BED (optional): `<bed_dir>/<virus>.<condition>.bed` (for condition `activity`, filename is `<virus>.bed`)
  Columns: `chrom  start0  end  sense_activity  antisense_activity` (0-based, half-open)

Usage
-----
# Single CSV (with one or many conditions)
python CRE_builder.py \
  --in ../../results/05_downstream/virus_tracks/csv/ad14_track.csv \
  --out-dir ../../results/05_downstream/cre_csv/ \
  --bed-dir ../../results/05_downstream/cre_bed/

# Batch a directory of track CSVs with a filename pattern (e.g., "*_track.csv" or "*.csv")
python CRE_builder.py \
  --in-dir ../../results/05_downstream/virus_tracks/csv/ \
  --pattern "*_track.csv" \
  --out-dir ../../results/05_downstream/cre_csv/ \
  --bed-dir ../../results/05_downstream/cre_bed/

# To disable the support filter: add --no-support-filter

Notes
-----
- Tiles are considered active for a *sign* only if the strand's activity magnitude ≥ 0.585 and matches that sign (> +0.585 for UP, < -0.585 for DOWN). Values with |activity| < 0.585 are treated as neutral (no sign).
- Conditions are inferred from column pairs: `<cond>` and `<cond>_antisense` (legacy `sense_strength`/`antisense_strength` → condition `activity`).
- Empty outputs are suppressed: we do not write CSVs or BEDs if there are no CREs for that condition/virus.
- Optional inclusion filter (default ON) can be disabled with --no-support-filter. Use --support-filter to turn it back on explicitly.
"""
from __future__ import annotations

import argparse
import os
import re
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple
from glob import glob

import pandas as pd

# Constants reflecting the MPRA tiling design
TILE_LEN = 200
SHIFT = 50
HALF = TILE_LEN // 2  # 100
SIGN_THRESH = 0.585  # absolute activity threshold to assign a sign


@dataclass
class Tile:
    center: int
    sense: float
    antisense: float

    @property
    def start(self) -> int:
        # 1-based inclusive
        return self.center - HALF

    @property
    def end(self) -> int:
        # 1-based inclusive
        return self.center + HALF


@dataclass
class CRE:
    virus: str
    first_center: int
    last_center: int
    max_sense: float
    max_antisense: float
    n_tiles: int = 1
    both_strands_hit: bool = False
    cross_tile_complement: bool = False  # True if an adjacent tile showed the complementary strand with same sign

    def extend_with(self, tile: Tile, s_on: bool, a_on: bool) -> None:
        self.last_center = tile.center
        self.n_tiles += 1
        # record if any tile in this CRE had both strands active with the current sign
        if s_on and a_on:
            self.both_strands_hit = True
        if tile.sense > self.max_sense:
            self.max_sense = tile.sense
        if tile.antisense > self.max_antisense:
            self.max_antisense = tile.antisense

    @property
    def start(self) -> int:
        return self.first_center - HALF

    @property
    def end(self) -> int:
        return self.last_center + HALF

    def to_row(self) -> Tuple[str, int, int, float, float]:
        return (self.virus, self.start, self.end, self.max_sense, self.max_antisense)


def _infer_virus_from_path(path: str) -> str:
    base = os.path.basename(path)
    # Examples: ad14_track.csv -> ad14, ebv_track.csv -> ebv, poly_bk_track.csv -> poly_bk
    m = re.match(r"(.+?)_track\.csv$", base)
    return m.group(1) if m else os.path.splitext(base)[0]


def _detect_conditions(df: pd.DataFrame) -> List[Tuple[str, str, str]]:
    """Return list of (sense_col, antisense_col, condition_name).

    - Pair columns like `<cond>` and `<cond>_antisense`.
    - Support legacy `sense_strength`/`antisense_strength` -> `activity`.
    - Ignore `bp` and any columns without a matching antisense partner.
    """
    cols = set(df.columns)
    conds: List[Tuple[str, str, str]] = []

    # Legacy support
    if {"sense_strength", "antisense_strength"}.issubset(cols):
        conds.append(("sense_strength", "antisense_strength", "activity"))

    # Generic pairs
    for c in cols:
        if c in {"bp", "sense_strength", "antisense_strength"}:
            continue
        a = f"{c}_antisense"
        if a in cols:
            conds.append((c, a, c))

    # Deduplicate while preserving order in dataframe
    seen = set()
    ordered: List[Tuple[str, str, str]] = []
    for col in df.columns:
        if col == "bp":
            continue
        if col == "sense_strength":
            key = ("sense_strength", "antisense_strength", "activity")
            if key not in seen and key in set(conds):
                seen.add(key)
                ordered.append(key)
        else:
            a = f"{col}_antisense"
            key = (col, a, col)
            if key in set(conds) and key not in seen:
                seen.add(key)
                ordered.append(key)
    return ordered


def _load_tiles_for_condition(df: pd.DataFrame, sense_col: str, antisense_col: str) -> List[Tile]:
    # Fill NaNs with 0 to treat missing as inactive
    df2 = df[["bp", sense_col, antisense_col]].copy()
    df2 = df2.sort_values("bp").reset_index(drop=True)
    df2[[sense_col, antisense_col]] = df2[[sense_col, antisense_col]].fillna(0.0)
    tiles = [Tile(int(row.bp), float(row._asdict()[sense_col]), float(row._asdict()[antisense_col])) for row in df2.itertuples(index=False)]
    return tiles


def build_cres_from_tiles(tiles: List[Tile], virus: str, sign: str = "pos", enforce_support_filter: bool = True) -> List[CRE]:
    """Collapse tile-level activities into CREs.

    Rules:
    - Consider only *sign-assigned* tiles (for sign="pos": sense ≥ SIGN_THRESH OR antisense ≥ SIGN_THRESH; for sign="neg": sense ≤ -SIGN_THRESH OR antisense ≤ -SIGN_THRESH).
    - Scan in order of `bp`.
    - Start a new CRE at the first active tile.
    - Extend the current CRE to include the next active tile **iff** its center is exactly
      previous_center + SHIFT (i.e., consecutive) **and** the pair shares at least one
      active strand (both sense-active, or both antisense-active). Otherwise, close the
      current CRE and start a new one at the new active tile.
    - A singleton active tile forms a 200 bp CRE.
    - CRE sense/antisense activities are the maxima across member tiles.
    - Include a CRE only if **either** (i) at least one tile within it has the complementary strand above the |SIGN_THRESH| with the **same sign**, **or** (ii) it contains at least two consecutive sign-matching tiles (i.e., not a singleton).
    - The above inclusion rule can be toggled with enforce_support_filter (True = apply rule; False = keep all CREs).
    """
    if not tiles:
        return []

    # Filter to active tiles only (keep bp order)
    active_tiles: List[Tuple[Tile, Tuple[bool, bool]]] = []
    for t in tiles:
        if sign == "neg":
            s_on = t.sense <= -SIGN_THRESH
            a_on = t.antisense <= -SIGN_THRESH
        else:
            s_on = t.sense >= SIGN_THRESH
            a_on = t.antisense >= SIGN_THRESH
        # If strands have opposite signs (e.g., sense ≥ +thr and antisense ≤ -thr),
        # the tile will be included in BOTH passes because we run pos and neg separately.
        if s_on or a_on:
            active_tiles.append((t, (s_on, a_on)))

    if not active_tiles:
        return []

    cres: List[CRE] = []

    # Initialize with the first active tile
    (t0, (s0, a0)) = active_tiles[0]
    current = CRE(virus=virus, first_center=t0.center, last_center=t0.center,
                  max_sense=t0.sense, max_antisense=t0.antisense,
                  n_tiles=1, both_strands_hit=(s0 and a0))

    prev_center = t0.center
    prev_flags = (s0, a0)

    for (t, (s_on, a_on)) in active_tiles[1:]:
        # Neighbor must share the same sign on at least one strand w.r.t the previous tile
        is_consecutive = (t.center - prev_center) == SHIFT
        shares_strand = (prev_flags[0] and s_on) or (prev_flags[1] and a_on)
        if is_consecutive and shares_strand:
            # record if complementary strand is active across the boundary
            # (e.g., previous tile sense active & current antisense active, or vice versa)
            if (prev_flags[0] and a_on) or (prev_flags[1] and s_on):
                current.cross_tile_complement = True
            current.extend_with(t, s_on, a_on)
        else:
            # close current
            if enforce_support_filter:
                # keep if any of the following holds:
                # (a) has ≥2 tiles (consecutive sign-matching tiles),
                # (b) any tile had both strands with matching sign (same tile complementary),
                # (c) complementary strand is active on an adjacent tile (cross-tile complementary)
                if current.n_tiles >= 2 or current.both_strands_hit or current.cross_tile_complement:
                    cres.append(current)
            else:
                cres.append(current)
            # start a new CRE at this tile
            current = CRE(virus=virus, first_center=t.center, last_center=t.center,
                          max_sense=t.sense, max_antisense=t.antisense,
                          n_tiles=1, both_strands_hit=(s_on and a_on))
        prev_center = t.center
        prev_flags = (s_on, a_on)

    # append the last CRE (respect filter toggle)
    if enforce_support_filter:
        if current.n_tiles >= 2 or current.both_strands_hit or current.cross_tile_complement:
            cres.append(current)
    else:
        cres.append(current)
    return cres


def _iter_track_csvs(in_dir: str, pattern: str) -> Iterable[str]:
    pattern_path = os.path.join(in_dir, pattern)
    for path in sorted(glob(pattern_path)):
        if os.path.isfile(path):
            yield path


def _write_bed(df: pd.DataFrame, path: str) -> None:
    # BED columns: chrom (virus), start0 (0-based), end (1-based), sense_activity, antisense_activity
    if df.empty:
        return  # skip creating empty BEDs
    bed_df = df.copy()
    bed_df["start0"] = bed_df["start"] - 1
    bed_df = bed_df[["virus", "start0", "end", "sense_activity", "antisense_activity"]]
    bed_df.rename(columns={"virus": "chrom"}, inplace=True)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    bed_df.to_csv(path, sep="\t", header=False, index=False)


def _process_csv_for_conditions(csv_path: str, virus: Optional[str], enforce_support_filter: bool) -> Dict[str, pd.DataFrame]:
    """Return mapping: condition -> CRE dataframe for that single virus file."""
    if virus is None:
        virus = _infer_virus_from_path(csv_path)
    df = pd.read_csv(csv_path)
    if "bp" not in df.columns:
        raise ValueError(f"{csv_path} missing 'bp' column")
    cond_specs = _detect_conditions(df)
    if not cond_specs:
        raise ValueError(f"{csv_path} has no detectable condition column pairs.")

    result: Dict[str, pd.DataFrame] = {}
    for sense_col, antisense_col, cond_name in cond_specs:
        tiles = _load_tiles_for_condition(df, sense_col, antisense_col)
        # Positive (up) condition
        cres_up = build_cres_from_tiles(tiles, virus=virus, sign="pos", enforce_support_filter=enforce_support_filter)
        rows_up = [c.to_row() for c in cres_up]
        result[cond_name] = pd.DataFrame(rows_up, columns=["virus", "start", "end", "sense_activity", "antisense_activity"])

        # Negative (down) condition
        cres_down = build_cres_from_tiles(tiles, virus=virus, sign="neg", enforce_support_filter=enforce_support_filter)
        rows_down = [c.to_row() for c in cres_down]
        result[f"{cond_name}_down"] = pd.DataFrame(rows_down, columns=["virus", "start", "end", "sense_activity", "antisense_activity"])
    return result


def main():
    p = argparse.ArgumentParser(description="Build CREs from MPRA tile-level activity tracks (supports multiple conditions).")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--in", dest="in_path", help="Input track CSV (columns: bp and pairs <cond>, <cond>_antisense)")
    g.add_argument("--in-dir", dest="in_dir", help="Directory containing track CSV files")
    p.add_argument("--pattern", default="*_track.csv", help="Filename pattern to match within --in-dir (e.g., '*_track.csv', '*.csv')")
    p.add_argument("--out-dir", required=True, help="Directory to write per-condition CRE CSVs")
    p.add_argument("--virus", default=None, help="Virus name (if not inferrable from filename)")
    p.add_argument("--no-header", action="store_true", help="Write CSVs without header rows")
    p.add_argument("--bed-dir", default=None, help="Directory for per-virus per-condition BED files (optional)")

    mx = p.add_mutually_exclusive_group()
    mx.add_argument("--support-filter", dest="support_filter", action="store_true", help="Require complementary strand sign or ≥2 consecutive sign-matching tiles to keep a CRE (default ON).")
    mx.add_argument("--no-support-filter", dest="support_filter", action="store_false", help="Disable the support filter; keep all sign-matched CREs (may include singleton CREs).")
    p.set_defaults(support_filter=True)

    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # condition -> list of dfs (across viruses)
    combined: Dict[str, List[pd.DataFrame]] = {}
    # For BEDs: condition -> virus -> df
    bed_map: Dict[str, Dict[str, pd.DataFrame]] = {}

    if args.in_path:
        vname = args.virus or _infer_virus_from_path(args.in_path)
        by_cond = _process_csv_for_conditions(args.in_path, vname, args.support_filter)
        for cond, dfv in by_cond.items():
            combined.setdefault(cond, []).append(dfv)
            bed_map.setdefault(cond, {})[vname] = dfv
    else:
        for csv_path in _iter_track_csvs(args.in_dir, args.pattern):
            vname = _infer_virus_from_path(csv_path)
            by_cond = _process_csv_for_conditions(csv_path, vname, args.support_filter)
            for cond, dfv in by_cond.items():
                combined.setdefault(cond, []).append(dfv)
                bed_map.setdefault(cond, {})[vname] = dfv

    # Write per-condition CSVs
    for cond, frames in combined.items():
        if frames:
            df = pd.concat(frames, ignore_index=True)
            if not df.empty:
                df["start"] = df["start"].astype(int)
                df["end"] = df["end"].astype(int)
                df.to_csv(out_path := os.path.join(args.out_dir, f"{cond}.csv"), index=False, header=not args.no_header)
            # else: skip writing empty CSV
        # else: no frames → skip writing

    # Write per-virus per-condition BEDs
    if args.bed_dir:
        for cond, virus_map in bed_map.items():
            for virus, dfv in virus_map.items():
                if dfv.empty:
                    continue  # do not write empty BEDs
                cond_clean = "activity" if cond == "activity" else cond
                if cond_clean == "activity":
                    bed_name = f"{virus}.bed"
                else:
                    bed_name = f"{virus}.{cond}.bed"
                bed_path = os.path.join(args.bed_dir, bed_name)
                _write_bed(dfv, bed_path)


if __name__ == "__main__":
    main()