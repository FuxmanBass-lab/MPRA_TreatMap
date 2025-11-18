#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import pathlib
import glob
from typing import Optional, List, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import gridspec
from matplotlib.colors import Normalize, ListedColormap

######################################################################
# Prefiltered/wide matrix workflow helpers
######################################################################

def detect_bp_col(df: pd.DataFrame) -> str:
    for col in ["bp", "pos", "position", "genomic_position"]:
        if col in df.columns:
            return col
    raise ValueError("Could not detect bp column (tried: bp, pos, position, genomic_position)")

def detect_paired_tracks(df: pd.DataFrame) -> list[tuple[str, str]]:
    pairs = []
    cols = set(df.columns)
    for c in sorted(df.columns):
        if c.endswith("(+)"):
            sense = c
            anti = c[:-3] + "(-)"
            if anti in cols:
                pairs.append((sense, anti))
    return sorted(pairs, key=lambda t: t[0])

def remaining_heatmap_tracks(df: pd.DataFrame, paired: list[tuple[str, str]]) -> list[str]:
    bp_col = detect_bp_col(df)
    paired_cols = set()
    for s, a in paired:
        paired_cols.add(s)
        paired_cols.add(a)
    return [c for c in df.columns if c != bp_col and c not in paired_cols]


# ---------------------------------------------------------------------
# Plotting: enhancer-style + heatmap panel
# ---------------------------------------------------------------------

def plot_prefiltered_matrix(
    df: pd.DataFrame,
    title: str,
    out_png: str,
    gene_df: pd.DataFrame = None,
    heatmap_vmin: float = 0.0,
    heatmap_vmax: float = 4.0,
    dpi: int = 300,
):
    """
    Plot prefiltered wide matrix: main panel for paired tracks (sense/antisense mirrored), heatmap for remaining tracks.
    """
    if df.empty:
        print(f"[warn] No data to plot for {title}; skipping")
        return

    bp_col = detect_bp_col(df)
    x_bp = df[bp_col].to_numpy()
    x_kb = x_bp / 1000.0
    genome_size = float(np.nanmax(x_bp) - np.nanmin(x_bp))
    width = max(12.0, min(28.0, 20.0 * (genome_size / 100_000.0)))

    paired = detect_paired_tracks(df)
    heatmap_tracks = remaining_heatmap_tracks(df, paired)
    want_main = len(paired) > 0
    want_heatmap = len(heatmap_tracks) > 0
    want_genes = gene_df is not None and not gene_df.empty

    # Layout: up to 3 panels
    if want_main and want_heatmap and want_genes:
        fig = plt.figure(figsize=(width, 3.0 + 0.7 + 1.4), dpi=dpi)
        gs = gridspec.GridSpec(3, 1, height_ratios=[3.0, 0.7, 1.4], hspace=0.25)
        ax = fig.add_subplot(gs[0])
        axh = fig.add_subplot(gs[1], sharex=ax)
        axg = fig.add_subplot(gs[2], sharex=ax)
    elif want_main and want_heatmap:
        fig = plt.figure(figsize=(width, 3.0 + 0.8), dpi=dpi)
        gs = gridspec.GridSpec(2, 1, height_ratios=[3.0, 0.8], hspace=0.25)
        ax = fig.add_subplot(gs[0])
        axh = fig.add_subplot(gs[1], sharex=ax)
        axg = None
    elif want_main and want_genes:
        fig = plt.figure(figsize=(width, 3.0 + 1.4), dpi=dpi)
        gs = gridspec.GridSpec(2, 1, height_ratios=[3.0, 1.4], hspace=0.25)
        ax = fig.add_subplot(gs[0])
        axh = None
        axg = fig.add_subplot(gs[1], sharex=ax)
    elif want_main:
        fig, ax = plt.subplots(figsize=(width, 3.0), dpi=dpi)
        axh = None
        axg = None
    elif want_heatmap and want_genes:
        fig = plt.figure(figsize=(width, 0.8 + 1.4), dpi=dpi)
        gs = gridspec.GridSpec(2, 1, height_ratios=[0.8, 1.4], hspace=0.25)
        ax = None
        axh = fig.add_subplot(gs[0])
        axg = fig.add_subplot(gs[1], sharex=axh)
    elif want_heatmap:
        fig, axh = plt.subplots(figsize=(width, 0.8), dpi=dpi)
        ax = None
        axg = None
    else:
        print(f"[warn] No paired or heatmap tracks to plot for {title}; skipping")
        return

    # ----- MAIN paired panel -----
    if want_main:
        ymax = 0.0
        for sense_col, anti_col in paired:
            y_sense = df[sense_col].to_numpy(dtype=float)
            y_anti = -df[anti_col].to_numpy(dtype=float)
            ax.plot(x_kb, y_sense, label=sense_col, linewidth=0.8)
            ax.plot(x_kb, y_anti, label=anti_col, linewidth=0.8)
            ymax = max(ymax, float(np.nanmax(np.abs(y_sense))), float(np.nanmax(np.abs(y_anti))))
        ax.axhline(0, color="black", linewidth=0.7)
        ax.set_ylabel("Paired track strength")
        ax.set_title(title)
        num_ticks = max(6, int(genome_size / 4000))
        ax.xaxis.set_major_locator(plt.MaxNLocator(num_ticks))
        if not np.isfinite(ymax) or ymax == 0:
            ymax = 1.0
        ax.set_ylim(-1.05 * ymax, 1.05 * ymax)
        ax.legend(frameon=False, loc="best", fontsize=7)
        if not want_heatmap and not want_genes:
            ax.set_xlabel("Genomic position (kb)")

    # ----- HEATMAP panel -----
    if want_heatmap:
        mat = []
        xh = df[bp_col].to_numpy()
        for t in heatmap_tracks:
            mat.append(df[t].to_numpy(dtype=float))
        mat = np.array(mat)  # shape (n_tracks, n_bp)
        cmap = plt.cm.Reds
        norm = Normalize(vmin=heatmap_vmin, vmax=heatmap_vmax, clip=True)
        extent = (xh.min() / 1000.0, xh.max() / 1000.0, 0, len(heatmap_tracks))
        axh.imshow(mat, aspect="auto", cmap=cmap, norm=norm, extent=extent, origin="lower")
        axh.set_yticks(np.arange(len(heatmap_tracks)) + 0.5)
        axh.set_yticklabels(heatmap_tracks, fontsize=6)
        axh.set_ylabel("non-control", fontsize=7)
        if want_main or want_genes:
            axh.tick_params(axis="x", labelbottom=False)
        else:
            axh.set_xlabel("Genomic position (kb)")
        for i in range(len(heatmap_tracks) + 1):
            axh.axhline(i, color="white", linewidth=0.2)

    # ----- GENE panel -----
    if want_genes:
        gene_edge = (0.20, 0.40, 0.70, 0.95)
        gene_fill = (0.20, 0.40, 0.70, 0.12)
        bar_h = 0.15
        def bp2kb(v):
            return np.asarray(v, dtype=float) / 1000.0
        genes = []
        for _, row in gene_df.iterrows():
            s_kb = bp2kb(row["start"]); e_kb = bp2kb(row["end"])
            if not np.isfinite(s_kb) or not np.isfinite(e_kb):
                continue
            if e_kb < s_kb:
                s_kb, e_kb = e_kb, s_kb
            strand_key = "+" if str(row.get("strand", "+")) == "+" else "-"
            genes.append((strand_key, float(s_kb), float(e_kb), str(row.get("name", ""))))
        y_base_plus = 0.1
        y_base_minus = -0.1
        if want_main or want_heatmap:
            ax_gene = axg
        else:
            ax_gene = axg
        for strand_key, s_kb, e_kb, name in sorted(genes, key=lambda t: (t[0], t[1], t[2])):
            if strand_key == "+":
                y0 = y_base_plus
            else:
                y0 = y_base_minus
            poly = mpatches.Rectangle((s_kb, y0), e_kb - s_kb, bar_h,
                                      edgecolor=gene_edge, facecolor=gene_fill, linewidth=0.5)
            ax_gene.add_patch(poly)
            ax_gene.text((s_kb + e_kb) / 2.0, y0 + bar_h + 0.02, name,
                         ha="center", va="bottom", fontsize=5.5, color=gene_edge, clip_on=False)
        ax_gene.set_ylim(-0.4, 0.5)
        ax_gene.set_yticks([])
        ax_gene.set_xlabel("Genomic position (kb)")
        ax_gene.spines["top"].set_visible(False)
        ax_gene.spines["right"].set_visible(False)
        ax_gene.spines["left"].set_visible(False)

    fig.tight_layout()
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    fig.savefig(out_png, dpi=dpi)
    plt.close(fig)
    print(f"[ok] wrote {out_png}")



# ---------------------------------------------------------------------
# activity helpers (same contract as your other project)
# ---------------------------------------------------------------------

def _activity_sign(activity: str) -> int:
    a = (activity or "silencer").strip().lower()
    if a not in {"silencer", "enhancer"}:
        raise ValueError(f"Unsupported activity: {activity}")
    return +1 if a == "enhancer" else -1


def _activity_label(activity: str) -> tuple[str, str]:
    sgn = _activity_sign(activity)
    if sgn == +1:
        return ("Activation", "log2FC")
    else:
        return ("Repression", "−log2FC")


# ---------------------------------------------------------------------
# IO helpers (copy of your style)
# ---------------------------------------------------------------------

def read_mapping(path: str | os.PathLike) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python")
    rename_map = {}
    for c in df.columns:
        lc = c.strip().lower()
        if lc in {"bp", "pos", "position", "genomic_position"}:
            rename_map[c] = "bp"
        elif lc.startswith("sense"):
            rename_map[c] = "sense"
        elif lc.startswith("antisense") or lc.startswith("anti"):
            rename_map[c] = "antisense"
        elif lc in {"number", "idx", "index"}:
            rename_map[c] = c
    df = df.rename(columns=rename_map)

    required = {"sense", "antisense", "bp"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing required columns: {sorted(missing)}")

    df = df[list(required | {"number"} & set(df.columns))].copy()
    df["bp"] = pd.to_numeric(df["bp"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["bp"]).copy()
    df["bp"] = df["bp"].astype(int)
    return df


def find_annotation_file(genomes_dir: str, virus_key: str) -> str | None:
    vdir = pathlib.Path(genomes_dir) / virus_key
    candidates: list[str] = []
    for ext in ("*.gtf", "*.gff3", "*.GTF", "*.GFF3"):
        candidates.extend(glob.glob(str(vdir / ext)))
    return candidates[0] if candidates else None


def _attr_get(attrs: str, keys=("gene_name", "Name", "gene", "gene_id", "locus_tag", "product")) -> str:
    for k in keys:
        for token in str(attrs).split(";"):
            token = token.strip()
            if not token:
                continue
            if token.startswith(k + "="):
                return token.split("=", 1)[1].strip().strip('"')
            if token.startswith(k + " "):
                return token.split(" ", 1)[1].strip().strip('"')
    return ""


def read_annotation_tidy(path: str) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"],
        dtype={"seqid": str, "start": int, "end": int, "strand": str, "type": str},
    )
    if df.empty:
        return pd.DataFrame(columns=["seqid", "start", "end", "strand", "name"])
    order = ["gene", "CDS", "mRNA", "exon", "ORF"]
    for t in order:
        sub = df.loc[df["type"].str.lower() == t.lower()].copy()
        if not sub.empty:
            df = sub
            break
    df["name"] = df["attributes"].astype(str).map(_attr_get)
    df["name"] = df["name"].replace({"": np.nan})
    if df["name"].isna().any():
        df.loc[df["name"].isna(), "name"] = df.loc[df["name"].isna()] \
            .apply(lambda r: f"{r['type']}:{int(r['start'])}-{int(r['end'])}", axis=1)
    keep = df[["seqid", "start", "end", "strand", "name"]].copy()
    agg = keep.groupby(["seqid", "name", "strand"], as_index=False).agg({"start": "min", "end": "max"})
    return agg[["seqid", "start", "end", "strand", "name"]]


# ---------------------------------------------------------------------
# MPRA reading (enhancer/silencer exactly like you do)
# ---------------------------------------------------------------------

def read_mpra(path: str | os.PathLike,
              padj_thr: float,
              min_lfc: float,
              min_ctrl_mean: float,
              activity: str) -> pd.DataFrame:
    sgn = _activity_sign(activity)
    df = pd.read_csv(path, sep=None, engine="python")
    df = df.rename(columns={
        "Id": "ID",
        "id": "ID",
        "log2fc": "log2FoldChange",
        "log2fcchange": "log2FoldChange",
        "padj": "padj",
    })
    required = {"ID", "log2FoldChange", "padj"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"MPRA table missing required columns: {sorted(missing)}")
    has_ctrl = "ctrl_mean" in df.columns
    if has_ctrl:
        df["ctrl_mean"] = pd.to_numeric(df["ctrl_mean"], errors="coerce")

    df["log2FoldChange"] = pd.to_numeric(df["log2FoldChange"], errors="coerce")
    df["padj"] = pd.to_numeric(df["padj"], errors="coerce")

    if activity.strip().lower() == "enhancer":
        sign_ok = df["log2FoldChange"] > 0
    else:
        sign_ok = df["log2FoldChange"] < 0

    base_filt = (
        df["log2FoldChange"].notna()
        & df["padj"].notna()
        & sign_ok
        & (df["padj"] <= padj_thr)
        & (df["log2FoldChange"].abs() >= min_lfc)
    )
    if has_ctrl:
        filt = base_filt & (df["ctrl_mean"].fillna(0) >= min_ctrl_mean)
    else:
        filt = base_filt

    out = df.loc[filt, ["ID", "log2FoldChange"]].copy()
    out["strength"] = (sgn * out["log2FoldChange"]).clip(lower=0.0)
    out = out.drop(columns=["log2FoldChange"])
    return out


def read_mpra_all_strengths(path: str | os.PathLike, activity: str) -> dict[str, float]:
    sgn = _activity_sign(activity)
    df = pd.read_csv(path, sep=None, engine="python")
    df = df.rename(columns={"Id": "ID", "id": "ID", "log2fc": "log2FoldChange"})
    if "ID" not in df.columns or "log2FoldChange" not in df.columns:
        raise ValueError("MPRA table must contain ID and log2FoldChange")
    df["log2FoldChange"] = pd.to_numeric(df["log2FoldChange"], errors="coerce")
    df = df.dropna(subset=["ID", "log2FoldChange"])
    df["strength_any"] = (sgn * df["log2FoldChange"]).clip(lower=0.0)
    return df.set_index("ID")["strength_any"].to_dict()


# ---------------------------------------------------------------------
# Track building (same as your project)
# ---------------------------------------------------------------------

def build_tracks(map_df: pd.DataFrame,
                 mpra_df: pd.DataFrame,
                 any_strength: dict[str, float],
                 include_complement_anyway: bool,
                 include_neighbor_anyway: bool,
                 enforce_support: bool) -> pd.DataFrame:
    id_to_strength = mpra_df.set_index("ID")["strength"].to_dict()
    keep_set = set(id_to_strength.keys())

    m = map_df.sort_values("bp").reset_index(drop=True).copy()
    sense_ids = m["sense"].astype(object).to_numpy()
    anti_ids = m["antisense"].astype(object).to_numpy()

    base_sense = np.array([sid in keep_set for sid in sense_ids], dtype=bool)
    base_anti = np.array([aid in keep_set for aid in anti_ids], dtype=bool)

    n = len(m)
    if n == 0:
        return pd.DataFrame({"bp": [], "sense_strength": [], "antisense_strength": []})

    if enforce_support:
        neigh_sense = np.zeros(n, dtype=bool)
        neigh_anti = np.zeros(n, dtype=bool)
        if n > 1:
            neigh_sense[1:] |= base_sense[:-1]
            neigh_sense[:-1] |= base_sense[1:]
            neigh_anti[1:] |= base_anti[:-1]
            neigh_anti[:-1] |= base_anti[1:]
        comp_for_sense = base_anti.copy()
        comp_for_anti = base_sense.copy()
        comp_neigh_sense = np.zeros(n, dtype=bool)
        comp_neigh_anti = np.zeros(n, dtype=bool)
        if n > 1:
            comp_neigh_sense[1:] |= base_anti[:-1]
            comp_neigh_sense[:-1] |= base_anti[1:]
            comp_neigh_anti[1:] |= base_sense[:-1]
            comp_neigh_anti[:-1] |= base_sense[1:]
        include_sense = base_sense & (neigh_sense | comp_for_sense | comp_neigh_sense)
        include_anti = base_anti & (neigh_anti | comp_for_anti | comp_neigh_anti)
    else:
        include_sense = base_sense
        include_anti = base_anti

    sense_sig = pd.Series(sense_ids).map(id_to_strength).fillna(0.0).to_numpy(float)
    anti_sig = pd.Series(anti_ids).map(id_to_strength).fillna(0.0).to_numpy(float)

    sense_strength = np.where(include_sense, sense_sig, 0.0)
    anti_strength = np.where(include_anti, anti_sig, 0.0)

    if include_complement_anyway:
        comp_anti_vals = pd.Series(anti_ids).map(any_strength).fillna(0.0).to_numpy(float)
        comp_sense_vals = pd.Series(sense_ids).map(any_strength).fillna(0.0).to_numpy(float)
        anti_strength = np.where(include_sense, np.maximum(anti_strength, comp_anti_vals), anti_strength)
        sense_strength = np.where(include_anti, np.maximum(sense_strength, comp_sense_vals), sense_strength)

    if include_neighbor_anyway:
        sense_any = pd.Series(sense_ids).map(any_strength).fillna(0.0).to_numpy(float)
        anti_any = pd.Series(anti_ids).map(any_strength).fillna(0.0).to_numpy(float)
        inc_idx = np.flatnonzero(include_sense | include_anti)
        for idx in inc_idx:
            left = idx - 1
            right = idx + 1
            if left >= 0:
                sense_strength[left] = max(sense_strength[left], sense_any[left])
                anti_strength[left] = max(anti_strength[left], anti_any[left])
            if right < n:
                sense_strength[right] = max(sense_strength[right], sense_any[right])
                anti_strength[right] = max(anti_strength[right], anti_any[right])

    out = pd.DataFrame({
        "bp": m["bp"].astype(int),
        "sense_strength": sense_strength.astype(float),
        "antisense_strength": anti_strength.astype(float),
    })
    out = out.drop_duplicates(subset=["bp"])
    return out


# ---------------------------------------------------------------------
# Heatmap helper: read wide matrix and keep selected columns
# ---------------------------------------------------------------------

def read_heatmap_matrix(path: str,
                        tracks: list[str] | None,
                        bp_colnames=("bp", "pos", "position", "genomic_position")) -> tuple[pd.DataFrame, str, list[str]]:
    df = pd.read_csv(path)
    # detect bp col
    bp_col = None
    for c in bp_colnames:
        if c in df.columns:
            bp_col = c
            break
    if bp_col is None:
        raise ValueError("heatmap matrix must have a 'bp' (or pos/position) column")

    if tracks:
        missing = [t for t in tracks if t not in df.columns]
        if missing:
            raise ValueError(f"heatmap tracks not found in {path}: {missing}")
        cols = [bp_col] + tracks
        df = df[cols].copy()
    else:
        # everything except bp
        tracks = [c for c in df.columns if c != bp_col]

    df = df.sort_values(bp_col)
    return df, bp_col, tracks


# ---------------------------------------------------------------------
# Plotting: enhancer-style + heatmap panel
# ---------------------------------------------------------------------

def plot_enhancer_with_heatmap(
    track_df: pd.DataFrame,
    title: str,
    out_png: str,
    y_label: str,
    gene_df: pd.DataFrame | None,
    heatmap_df: pd.DataFrame | None,
    heatmap_bp: str | None,
    heatmap_tracks: list[str] | None,
    heatmap_vmin: float,
    heatmap_vmax: float,
    dpi: int = 300,
):
    """
    Top: your enhancer-style mirrored line plot.
    Bottom (optional): stacked heatmap rows, one per non-control track.
    """
    if track_df.empty:
        print(f"[warn] No data to plot for {title}; skipping")
        return

    x_bp = track_df["bp"].to_numpy()
    x_kb = x_bp / 1000.0
    y_sense = track_df["sense_strength"].to_numpy()
    y_anti = -track_df["antisense_strength"].to_numpy()

    genome_size = float(np.nanmax(x_bp) - np.nanmin(x_bp))
    width = max(12.0, min(28.0, 20.0 * (genome_size / 100_000.0)))

    # how many panels?
    want_heatmap = heatmap_df is not None and heatmap_tracks
    want_genes = gene_df is not None and not gene_df.empty

    if want_heatmap and want_genes:
        # 3 rows: main, heatmap, genes
        fig = plt.figure(figsize=(width, 3.0 + 0.7 + 1.4), dpi=dpi)
        gs = gridspec.GridSpec(3, 1, height_ratios=[3.0, 0.7, 1.4], hspace=0.25)
        ax = fig.add_subplot(gs[0])
        axh = fig.add_subplot(gs[1], sharex=ax)
        axg = fig.add_subplot(gs[2], sharex=ax)
    elif want_heatmap:
        fig = plt.figure(figsize=(width, 3.0 + 0.8), dpi=dpi)
        gs = gridspec.GridSpec(2, 1, height_ratios=[3.0, 0.8], hspace=0.25)
        ax = fig.add_subplot(gs[0])
        axh = fig.add_subplot(gs[1], sharex=ax)
        axg = None
    elif want_genes:
        fig = plt.figure(figsize=(width, 3.0 + 1.4), dpi=dpi)
        gs = gridspec.GridSpec(2, 1, height_ratios=[3.0, 1.4], hspace=0.25)
        ax = fig.add_subplot(gs[0])
        axh = None
        axg = fig.add_subplot(gs[1], sharex=ax)
    else:
        fig, ax = plt.subplots(figsize=(width, 3.0), dpi=dpi)
        axh = None
        axg = None

    # ----- MAIN enhancer-style plot -----
    ax.plot(x_kb, y_sense, label="+ strand (sense)", linewidth=0.9)
    ax.plot(x_kb, y_anti, label="- strand (antisense)", linewidth=0.9)
    ax.axhline(0, color="black", linewidth=0.7)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    num_ticks = max(6, int(genome_size / 4000))
    ax.xaxis.set_major_locator(plt.MaxNLocator(num_ticks))
    ymax = max(float(np.nanmax(np.abs(y_sense))), float(np.nanmax(np.abs(y_anti))))
    if not np.isfinite(ymax) or ymax == 0:
        ymax = 1.0
    ax.set_ylim(-1.05 * ymax, 1.05 * ymax)
    ax.legend(frameon=False, loc="best")

    # ----- HEATMAP PANEL -----
    if axh is not None:
        # we want: rows = heatmap_tracks, cols = bp
        mat = []
        xh = heatmap_df[heatmap_bp].to_numpy()
        for t in heatmap_tracks:
            mat.append(heatmap_df[t].to_numpy(dtype=float))
        mat = np.array(mat)  # shape (n_tracks, n_bp)

        # use a simple colormap; you can tweak to match R
        cmap = plt.cm.Reds
        norm = Normalize(vmin=heatmap_vmin, vmax=heatmap_vmax, clip=True)

        # imshow: x is bp (in kb), y is track index
        extent = (xh.min() / 1000.0, xh.max() / 1000.0, 0, len(heatmap_tracks))
        axh.imshow(mat, aspect="auto", cmap=cmap, norm=norm, extent=extent, origin="lower")
        axh.set_yticks(np.arange(len(heatmap_tracks)) + 0.5)
        axh.set_yticklabels(heatmap_tracks, fontsize=6)
        axh.set_ylabel("non-control", fontsize=7)
        axh.tick_params(axis="x", labelbottom=False)
        # thin grid lines
        for i in range(len(heatmap_tracks) + 1):
            axh.axhline(i, color="white", linewidth=0.2)

    # ----- GENE PANEL -----
    if axg is not None:
        # basic rectangles like your other script
        gene_edge = (0.20, 0.40, 0.70, 0.95)
        gene_fill = (0.20, 0.40, 0.70, 0.12)
        bar_h = 0.15

        def bp2kb(v):
            return np.asarray(v, dtype=float) / 1000.0

        genes = []
        for _, row in gene_df.iterrows():
            s_kb = bp2kb(row["start"]); e_kb = bp2kb(row["end"])
            if not np.isfinite(s_kb) or not np.isfinite(e_kb):
                continue
            if e_kb < s_kb:
                s_kb, e_kb = e_kb, s_kb
            strand_key = "+" if str(row.get("strand", "+")) == "+" else "-"
            genes.append((strand_key, float(s_kb), float(e_kb), str(row.get("name", ""))))

        # simple lanes
        y_base_plus = 0.1
        y_base_minus = -0.1
        for strand_key, s_kb, e_kb, name in sorted(genes, key=lambda t: (t[0], t[1], t[2])):
            if strand_key == "+":
                y0 = y_base_plus
            else:
                y0 = y_base_minus
            poly = mpatches.Rectangle((s_kb, y0), e_kb - s_kb, bar_h,
                                      edgecolor=gene_edge, facecolor=gene_fill, linewidth=0.5)
            axg.add_patch(poly)
            axg.text((s_kb + e_kb) / 2.0, y0 + bar_h + 0.02, name,
                     ha="center", va="bottom", fontsize=5.5, color=gene_edge, clip_on=False)
        axg.set_ylim(-0.4, 0.5)
        axg.set_yticks([])
        axg.set_xlabel("Genomic position (kb)")
        axg.spines["top"].set_visible(False)
        axg.spines["right"].set_visible(False)
        axg.spines["left"].set_visible(False)
    else:
        ax.set_xlabel("Genomic position (kb)")

    fig.tight_layout()
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    fig.savefig(out_png, dpi=dpi)
    plt.close(fig)
    print(f"[ok] wrote {out_png}")


# ---------------------------------------------------------------------
# main
# ---------------------------------------------------------------------

def virus_title_from_filename(fp: str | os.PathLike) -> str:
    name = pathlib.Path(fp).stem
    if name.endswith("_strands_match"):
        name = name[: -len("_strands_match")]
    return name


def main():
    p = argparse.ArgumentParser(
        description="Plot genome-wide enhancer/silencer tracks (sense vs antisense) + optional heatmap panel."
    )
    p.add_argument("--map-dir", default=None, help="Directory with per-virus *_strands_match.csv files (required unless --matrix is used)")
    p.add_argument("--mpra", default=None, help="MPRA results table (required unless --matrix is used)")
    p.add_argument("--outdir", required=True, help="Where to write per-virus PNGs")
    p.add_argument("--csv-outdir", default=None, help="Where to write per-virus track CSVs")
    p.add_argument("--padj", type=float, default=0.05)
    p.add_argument("--min-lfc", type=float, default=1.0)
    p.add_argument("--min-ctrl-mean", type=float, default=50.0)
    p.add_argument("--activity", choices=["enhancer", "silencer"], default="enhancer")
    p.add_argument("--include-complement-anyway", action="store_true", default=False)
    p.add_argument("--include-neighbor-anyway", action="store_true", default=False)
    p.add_argument("--no-support-filter", action="store_true", default=False)
    p.add_argument("--genomes-dir", default=None)
    p.add_argument("--plot-genes", action="store_true", default=False)

    # NEW: heatmap extras
    p.add_argument("--heatmap-matrix", default=None,
                   help="Wide CSV: bp + non-control columns to show as stacked heatmap")
    p.add_argument("--heatmap-tracks", default=None,
                   help="Comma-separated list of columns from --heatmap-matrix to show (in order)")
    p.add_argument("--heatmap-vmin", type=float, default=0.0)
    p.add_argument("--heatmap-vmax", type=float, default=4.0)

    # NEW: prefiltered matrix workflow
    group = p.add_argument_group("Prefiltered matrix workflow")
    group.add_argument("--matrix", default=None, type=str,
                       help="Prefiltered wide CSV (from R pipeline) for direct plotting")
    group.add_argument("--matrix-title", default=None, type=str,
                       help="Plot title for --matrix workflow (default: filename stem)")

    args = p.parse_args()
    enforce_support = not args.no_support_filter

    # Prefiltered matrix workflow: branch early if --matrix is provided
    if args.matrix is not None:
        import pathlib as _pathlib
        import os as _os
        df = pd.read_csv(args.matrix)
        plot_prefiltered_matrix(
            df,
            title=(args.matrix_title or _pathlib.Path(args.matrix).stem),
            out_png=_os.path.join(args.outdir, (args.matrix_title or _pathlib.Path(args.matrix).stem) + ".png"),
            gene_df=None,
            heatmap_vmin=args.heatmap_vmin,
            heatmap_vmax=args.heatmap_vmax,
            dpi=300
        )
        return

    # Validate required arguments for original workflow
    if args.map_dir is None or args.mpra is None:
        p.error("--map-dir and --mpra are required unless --matrix is provided")

    # Else: original workflow
    mpra = read_mpra(
        args.mpra,
        padj_thr=args.padj,
        min_lfc=args.min_lfc,
        min_ctrl_mean=args.min_ctrl_mean,
        activity=args.activity,
    )
    any_strength = read_mpra_all_strengths(args.mpra, activity=args.activity)
    noun, y_suffix = _activity_label(args.activity)
    y_label = f"{noun} strength ({y_suffix})"

    map_dir = pathlib.Path(args.map_dir)
    mapping_paths = sorted(map_dir.glob("*_strands_match.csv"))
    if not mapping_paths:
        raise FileNotFoundError(f"No mapping CSVs found in {args.map_dir}")

    os.makedirs(args.outdir, exist_ok=True)
    csv_outdir = args.csv_outdir or args.outdir
    os.makedirs(csv_outdir, exist_ok=True)

    # heatmap matrix (optional)
    if args.heatmap_matrix:
        tracks = [t.strip() for t in args.heatmap_tracks.split(",")] if args.heatmap_tracks else None
        heatmap_df, heatmap_bp, heatmap_cols = read_heatmap_matrix(args.heatmap_matrix, tracks)
    else:
        heatmap_df = None
        heatmap_bp = None
        heatmap_cols = None

    for mp in mapping_paths:
        virus = virus_title_from_filename(mp)
        mdf = read_mapping(mp)
        track = build_tracks(mdf, mpra, any_strength,
                             include_complement_anyway=args.include_complement_anyway,
                             include_neighbor_anyway=args.include_neighbor_anyway,
                             enforce_support=enforce_support)
        # write CSV
        track_csv = os.path.join(csv_outdir, f"{virus}_track.csv")
        track.to_csv(track_csv, index=False)

        # gene selection
        gene_df = None
        if args.plot_genes and args.genomes_dir:
            ann_path = find_annotation_file(args.genomes_dir, virus)
            if ann_path and os.path.exists(ann_path):
                try:
                    gfull = read_annotation_tidy(ann_path)
                    min_bp = int(mdf["bp"].min()); max_bp = int(mdf["bp"].max())
                    best_sid = None; best_hits = -1; best_span = -1
                    for sid, sub in gfull.groupby("seqid"):
                        hits = ((sub["end"] >= min_bp) & (sub["start"] <= max_bp)).sum()
                        span = (sub["end"].max() - sub["start"].min()) if not sub.empty else -1
                        if (hits > best_hits) or (hits == best_hits and span > best_span):
                            best_sid, best_hits, best_span = sid, hits, span
                    if best_sid is not None:
                        gene_df = gfull.loc[gfull["seqid"] == best_sid, ["start", "end", "strand", "name"]].copy()
                except Exception as e:
                    print(f"[warn] could not parse genes for {virus}: {e}")

        out_png = os.path.join(args.outdir, f"{virus}.png")
        plot_enhancer_with_heatmap(
            track_df=track,
            title=virus,
            out_png=out_png,
            y_label=y_label,
            gene_df=gene_df,
            heatmap_df=heatmap_df,
            heatmap_bp=heatmap_bp,
            heatmap_tracks=heatmap_cols,
            heatmap_vmin=args.heatmap_vmin,
            heatmap_vmax=args.heatmap_vmax,
        )


if __name__ == "__main__":
    main()