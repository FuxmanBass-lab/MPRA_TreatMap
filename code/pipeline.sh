#!/usr/bin/env bash
set -euo pipefail

# ---------------------------
# MPRA Heatmap Pipeline (DRY)
# ---------------------------
# This script replaces repetitive per-virus blocks with a configurable loop.
# Each virus is defined once below, and the same 3-step pipeline runs:
#   1) 0_preprocess_tiles.py
#   2) 0_pull_seq.py
#   3) 2_plot_heatmap.r (standard and combined tracks)
#
# You can add a new virus by adding entries to the config section only.

# ------- Defaults -------
CONTROL_LFC=2
TREAT_LFC=0.585
DNA_THRESHOLD=50
VERBOSE="--verbose"

# Default condition set (11-condition panel)
DEF_CONDS="Awater,Kwater,ILB1,INFa,INFg,POLYDADT,SENDAI,TNF"
DEF_CTLS="Awater,Kwater"
DEF_ORDER="Kwater,ILB1,INFg,TNF,Awater,INFa,POLYDADT,SENDAI"
DEF_KB_STEP=2000

# ---------------------------
# Per‑virus configuration
# Key names must match the directory names used under ../results/
# ---------------------------
# Required fields per key:
#   MAP[key]          : path to strands map CSV
#   ACCESSION[key]    : NCBI accession to embed in BED
#   PREFIX[key]       : file prefix for BED
#   GENES_BED[key]    : path to nearby genes BED (omit or leave empty if none)
# Optional fields:
#   GENES_ANNO[key]   : genes annotation TSV (omit if none)
#   CONDS[key]        : CSV list of conditions
#   CTLS[key]         : CSV list of control conditions
#   ORDER[key]        : CSV order for plotting
#   TITLE[key]        : Human-friendly figure title
#   OUT_EXT[key]      : png or pdf (default png)
#   GENE_SCALE[key]   : numeric scale for gene track (e.g., 1.5). If set, enables --gene-track-scale
#   RUN_COMBINED[key] : if "no", skip the combined track figure (default runs both)
#   EXTRA_STD[key]    : extra args to pass to 2_plot_heatmap.r for standard plot
#   EXTRA_COMB[key]   : extra args to pass to 2_plot_heatmap.r for combined plot
#   TITLE_COMB[key]   : Human-friendly figure title for combined plot

declare -A MAP ACCESSION PREFIX GENES_BED GENES_ANNO CONDS CTLS ORDER TITLE OUT_EXT GENE_SCALE RUN_COMBINED EXTRA_STD EXTRA_COMB TITLE_COMB

# --- EBV ---
MAP[EBV]=../data/strands/ebv_strands_match.csv
ACCESSION[EBV]=NC_007605.1
PREFIX[EBV]=EBV
GENES_BED[EBV]=../results/EBV/EBV_nearby_genes.bed
GENES_ANNO[EBV]=../results/EBV/gene_annotations.tsv
TITLE[EBV]="MPRA Differential Activity – EBV"
GENE_SCALE[EBV]=1.5

# --- HCMV ---
MAP[HCMV]=../data/strands/hcmv_strands_match.csv
ACCESSION[HCMV]=FJ616285.1
PREFIX[HCMV]=HCMV
GENES_BED[HCMV]=../results/HCMV/HCMV_nearby_genes.bed
GENES_ANNO[HCMV]=../results/HCMV/gene_annotations.tsv
TITLE[HCMV]="MPRA Differential Activity – HCMV"
EXTRA_STD[HCMV]="--height=18 --width=22 --gene-gap 0.5 --gene-label-nudge 6 --gene-label-padding 1.5 --gene-label-force 1.5"
EXTRA_COMB[HCMV]="--height=18 --width=22 --gene-gap 0.5 --gene-label-nudge 1 --gene-label-padding 1.5 --gene-label-force 1.5"

# --- KHSV ---
MAP[KHSV]=../data/strands/khsv_strands_match.csv
ACCESSION[KHSV]=MZ712172.1
PREFIX[KHSV]=KHSV
GENES_BED[KHSV]=../results/KHSV/KHSV_nearby_genes.bed
GENES_ANNO[KHSV]=../results/KHSV/gene_annotations.tsv
TITLE[KHSV]="MPRA Differential Activity – KHSV"
GENE_SCALE[KHSV]=1.5

# --- VZV ---
MAP[VZV]=../data/strands/vzv_strands_match.csv
ACCESSION[VZV]=KU926311.1
PREFIX[VZV]=VZV
GENES_BED[VZV]=../results/VZV/VZV_nearby_genes.bed
GENES_ANNO[VZV]=../results/VZV/gene_annotations.tsv
TITLE[VZV]="MPRA Differential Activity – VZV"
GENE_SCALE[VZV]=1.5

# --- HSV1 ---
MAP[HSV1]=../data/strands/hsv1_strands_match.csv
ACCESSION[HSV1]=KT899744.1
PREFIX[HSV1]=HSV1
GENES_BED[HSV1]=../results/HSV1/HSV1_nearby_genes.bed
GENES_ANNO[HSV1]=../results/HSV1/gene_annotations.tsv
TITLE[HSV1]="MPRA Differential Activity – HSV1"
GENE_SCALE[HSV1]=1.5
TITLE_COMB[HSV1]="MPRA Differential Activity – HSV1 KOS (combined)"

# --- HSV2 ---
MAP[HSV2]=../data/strands/hsv2_strands_match.csv
ACCESSION[HSV2]=OM370995.1
PREFIX[HSV2]=HSV2
GENES_BED[HSV2]=../results/HSV2/HSV2_nearby_genes.bed
TITLE[HSV2]="MPRA Differential Activity – HSV2"
GENE_SCALE[HSV2]=1.5

# --- HHV6B ---
MAP[HHV6B]=../data/strands/hhv6b_strands_match.csv
ACCESSION[HHV6B]=MW536483.1
PREFIX[HHV6B]=HHV6B
GENES_BED[HHV6B]=../results/HHV6B/HHV6B_nearby_genes.bed
TITLE[HHV6B]="MPRA Differential Activity – HHV6B"

# --- HHV7 ---
MAP[HHV7]=../data/strands/hhv7_strands_match.csv
ACCESSION[HHV7]=AF037218.1
PREFIX[HHV7]=HHV7
GENES_BED[HHV7]=../results/HHV7/HHV7_nearby_genes.bed
TITLE[HHV7]="MPRA Differential Activity – HHV7"

# --- HIV1 CH058 ---
MAP[HIV1_CH058]=../data/strands/hiv1_ch058_strands_match.csv
ACCESSION[HIV1_CH058]=JN944907
PREFIX[HIV1_CH058]=HIV1_CH058
GENES_BED[HIV1_CH058]=../results/HIV1_CH058/HIV1_CH058_nearby_genes.bed
TITLE[HIV1_CH058]="MPRA Differential Activity – HIV1 CH058"

# --- HIV2 ROD ---
MAP[HIV2_ROD]=../data/strands/hiv2_rod_strands_match.csv
ACCESSION[HIV2_ROD]=M15390.1
PREFIX[HIV2_ROD]=HIV2_ROD
GENES_BED[HIV2_ROD]=../results/HIV2_ROD/HIV2_ROD_nearby_genes.bed
TITLE[HIV2_ROD]="MPRA Differential Activity – HIV2 ROD"

# --- HIV2 GH-1 ---
MAP[HIV2_GH]=../data/strands/hiv2_gh_strands_match.csv
ACCESSION[HIV2_GH]=325709
PREFIX[HIV2_GH]=HIV2_GH
GENES_BED[HIV2_GH]=../results/HIV2_GH/HIV2_GH_nearby_genes.bed
TITLE[HIV2_GH]="MPRA Differential Activity – HIV2 GH-1"

# --- HTLV1 ---
MAP[HTLV1]=../data/strands/htlv1_strands_match.csv
ACCESSION[HTLV1]=MH392265.1
PREFIX[HTLV1]=HTLV1
GENES_BED[HTLV1]=../results/HTLV1/HTLV1_nearby_genes.bed
TITLE[HTLV1]="MPRA Differential Activity – HTLV1"

# --- HTLV2 ---
MAP[HTLV2]=../data/strands/htlv2_strands_match.csv
ACCESSION[HTLV2]=M10060
PREFIX[HTLV2]=HTLV2
GENES_BED[HTLV2]=../results/HTLV2/HTLV2_nearby_genes.bed
TITLE[HTLV2]="MPRA Differential Activity – HTLV2"

# --- HTLV3 ---
MAP[HTLV3]=../data/strands/htlv3_strands_match.csv
ACCESSION[HTLV3]=DQ093792.1
PREFIX[HTLV3]=HTLV3
GENES_BED[HTLV3]=../results/HTLV3/HTLV3_nearby_genes.bed
TITLE[HTLV3]="MPRA Differential Activity – HTLV3"

# --- Polyomavirus BK ---
MAP[POLY_BK]=../data/strands/poly_bk_strands_match.csv
ACCESSION[POLY_BK]=J02039.1
PREFIX[POLY_BK]=POLY_BK
GENES_BED[POLY_BK]=../results/POLY_BK/POLY_BK_nearby_genes.bed
TITLE[POLY_BK]="MPRA Differential Activity – POLY_BK"

# --- Polyomavirus JC ---
MAP[POLY_JC]=../data/strands/poly_jc_strands_match.csv
ACCESSION[POLY_JC]=NC_001699.1
PREFIX[POLY_JC]=POLY_JC
GENES_BED[POLY_JC]=../results/POLY_JC/POLY_JC_nearby_genes.bed
TITLE[POLY_JC]="MPRA Differential Activity – POLY_JC"

# --- Merkel Cell Polyomavirus ---
MAP[POLY_MERK]=../data/strands/poly_merk_strands_match.csv
ACCESSION[POLY_MERK]=207705714
PREFIX[POLY_MERK]=POLY_MERK
GENES_BED[POLY_MERK]=../results/POLY_MERK/POLY_MERK_nearby_genes.bed
TITLE[POLY_MERK]="MPRA Differential Activity – POLY_MERK"

# --- HPV1 ---
MAP[HPV1]=../data/strands/hpv1_strands_match.csv
ACCESSION[HPV1]=NC_001356.1
PREFIX[HPV1]=HPV1
GENES_BED[HPV1]=../results/HPV1/HPV1_nearby_genes.bed
TITLE[HPV1]="MPRA Differential Activity – HPV1"

# --- HPV2 ---
MAP[HPV2]=../data/strands/hpv2_strands_match.csv
ACCESSION[HPV2]=MN605988.1
PREFIX[HPV2]=HPV2
GENES_BED[HPV2]=../results/HPV2/HPV2_nearby_genes.bed
TITLE[HPV2]="MPRA Differential Activity – HPV2"

# --- HPV5 ---
MAP[HPV5]=../data/strands/hpv5_strands_match.csv
ACCESSION[HPV5]=JN211194.1
PREFIX[HPV5]=HPV5
GENES_BED[HPV5]=../results/HPV5/HPV5_nearby_genes.bed
TITLE[HPV5]="MPRA Differential Activity – HPV5"

# --- HPV6B ---
MAP[HPV6B]=../data/strands/hpv6b_strands_match.csv
ACCESSION[HPV6B]=X00203
PREFIX[HPV6B]=HPV6B
GENES_BED[HPV6B]=../results/HPV6B/HPV6B_nearby_genes.bed
TITLE[HPV6B]="MPRA Differential Activity – HPV6B"

# --- HPV11 ---
MAP[HPV11]=../data/strands/hpv11_strands_match.csv
ACCESSION[HPV11]=M14119
PREFIX[HPV11]=HPV11
GENES_BED[HPV11]=../results/HPV11/HPV11_nearby_genes.bed
TITLE[HPV11]="MPRA Differential Activity – HPV11"

# --- HPV16 ---
MAP[HPV16]=../data/strands/hpv16_strands_match.csv
ACCESSION[HPV16]=NC_001526.4
PREFIX[HPV16]=HPV16
GENES_BED[HPV16]=../results/HPV16/HPV16_nearby_genes.bed
TITLE[HPV16]="MPRA Differential Activity – HPV16"

# --- HPV18 ---
MAP[HPV18]=../data/strands/hpv18_strands_match.csv
ACCESSION[HPV18]=X05015
PREFIX[HPV18]=HPV18
GENES_BED[HPV18]=../results/HPV18/HPV18_nearby_genes.bed
TITLE[HPV18]="MPRA Differential Activity – HPV18"

# --- HIV1 Rejo ---
MAP[HIV1_rejo]=../data/strands/hiv1_rejo_strands_match.csv
ACCESSION[HIV1_rejo]=JN944911
PREFIX[HIV1_rejo]=HIV1_rejo
GENES_BED[HIV1_rejo]=../results/HIV1_rejo/HIV1_rejo_nearby_genes.bed
TITLE[HIV1_rejo]="MPRA Differential Activity – HIV1 Rejo"
EXTRA_STD[HIV1_rejo]="--height=8 --width=16"
EXTRA_COMB[HIV1_rejo]="--height=8 --width=16"

# --- HPV52 ---
MAP[HPV52]=../data/strands/hpv52_strands_match.csv
ACCESSION[HPV52]=X74481.1
PREFIX[HPV52]=HPV52
GENES_BED[HPV52]=../results/HPV52/HPV52_nearby_genes.bed
TITLE[HPV52]="MPRA Differential Activity – HPV52"
TITLE_COMB[HPV52]="MPRA Differential Activity – HPV52 (combined)"

# --- Adenoviruses ---
MAP[AD1]=../data/strands/ad1_strands_match.csv
ACCESSION[AD1]=OQ518269.1
PREFIX[AD1]=AD1
GENES_BED[AD1]=../results/AD1/AD1_nearby_genes.bed
TITLE[AD1]="MPRA Differential Activity – AD1"
TITLE_COMB[AD1]="MPRA Differential Activity – AD1 (combined)"

MAP[AD3]=../data/strands/ad3_strands_match.csv
ACCESSION[AD3]=AY599834.1
PREFIX[AD3]=AD3
GENES_BED[AD3]=../results/AD3/AD3_nearby_genes.bed
TITLE[AD3]="MPRA Differential Activity – AD3"
TITLE_COMB[AD3]="MPRA Differential Activity – AD3 (combined)"

MAP[AD4]=../data/strands/ad4_strands_match.csv
ACCESSION[AD4]=KX384949.1
PREFIX[AD4]=AD4
GENES_BED[AD4]=../results/AD4/AD4_nearby_genes.bed
TITLE[AD4]="MPRA Differential Activity – AD4"
TITLE_COMB[AD4]="MPRA Differential Activity – AD4 (combined)"

MAP[AD5]=../data/strands/ad5_strands_match.csv
ACCESSION[AD5]=AY601635.1
PREFIX[AD5]=AD5
GENES_BED[AD5]=../results/AD5/AD5_nearby_genes.bed
TITLE[AD5]="MPRA Differential Activity – AD5"
TITLE_COMB[AD5]="MPRA Differential Activity – AD5 (combined)"

MAP[AD7]=../data/strands/ad7_strands_match.csv
ACCESSION[AD7]=AY594255
PREFIX[AD7]=AD7
GENES_BED[AD7]=../results/AD7/AD7_nearby_genes.bed
TITLE[AD7]="MPRA Differential Activity – AD7"
TITLE_COMB[AD7]="MPRA Differential Activity – AD7 (combined)"

MAP[AD11]=../data/strands/ad11_strands_match.csv
ACCESSION[AD11]=AF532578.1
PREFIX[AD11]=AD11
GENES_BED[AD11]=../results/AD11/AD11_nearby_genes.bed
TITLE[AD11]="MPRA Differential Activity – AD11"
TITLE_COMB[AD11]="MPRA Differential Activity – AD11 (combined)"

MAP[AD14]=../data/strands/ad14_strands_match.csv
ACCESSION[AD14]=AY803294.1
PREFIX[AD14]=AD14
GENES_BED[AD14]=../results/AD14/AD14_nearby_genes.bed
TITLE[AD14]="MPRA Differential Activity – AD14"
TITLE_COMB[AD14]="MPRA Differential Activity – AD14 (combined)"

MAP[AD37]=../data/strands/ad37_strands_match.csv
ACCESSION[AD37]=AB448775.1
PREFIX[AD37]=AD37
GENES_BED[AD37]=../results/AD37/AD37_nearby_genes.bed
TITLE[AD37]="MPRA Differential Activity – AD37"
TITLE_COMB[AD37]="MPRA Differential Activity – AD3 (combined)"

# --- "Full" small panel variants (PDF outputs) ---
MAP[EBV_BB]=../data/strands/ebv_strands_match.csv
ACCESSION[EBV_BB]=NC_007605.1
PREFIX[EBV_BB]=EBV
GENES_BED[EBV_BB]=../results/EBV_BB/EBV_nearby_genes.bed
GENES_ANNO[EBV_BB]=../results/EBV_BB/gene_annotations.tsv
CONDS[EBV_BB]="H2O,H2O2,Oligo,TNF,INFG"
CTLS[EBV_BB]="H2O"
ORDER[EBV_BB]="H2O,H2O2,Oligo,TNF,INFG"
TITLE[EBV_BB]="MPRA Differential Activity – EBV"
OUT_EXT[EBV_BB]=pdf

MAP[HCMV_BB]=../data/strands/hcmv_strands_match.csv
ACCESSION[HCMV_BB]=FJ616285.1
PREFIX[HCMV_BB]=HCMV
GENES_BED[HCMV_BB]=../results/HCMV_BB/HCMV_nearby_genes.bed
GENES_ANNO[HCMV_BB]=../results/HCMV_BB/gene_annotations.tsv
CONDS[HCMV_BB]="H2O,H2O2,Oligo,TNF,INFG"
CTLS[HCMV_BB]="H2O"
ORDER[HCMV_BB]="H2O,H2O2,Oligo,TNF,INFG"
TITLE[HCMV_BB]="MPRA Differential Activity – HCMV"
OUT_EXT[HCMV_BB]=pdf

# ---------------------------
# Helpers
# ---------------------------
function cfg() { # get config with default
  local arr="$1" key="$2" def="$3"
  local val=""
  eval "val=\"\${$arr[$key]:-}\""
  if [[ -z "${val}" ]]; then
    echo "${def}"
  else
    echo "${val}"
  fi
}

run_pipeline_for() {
  local key="$1"
  local outdir="../results/${key}"

  local map_file="${MAP[$key]}"
  local accession="${ACCESSION[$key]}"
  local prefix="${PREFIX[$key]}"
  local genes_bed="$(cfg GENES_BED "$key" "")"
  local genes_anno="$(cfg GENES_ANNO "$key" "")"

  local conds="$(cfg CONDS "$key" "$DEF_CONDS")"
  local ctls="$(cfg CTLS "$key" "$DEF_CTLS")"
  local order="$(cfg ORDER "$key" "$DEF_ORDER")"
  local title="$(cfg TITLE "$key" "$key")"
  local out_ext="$(cfg OUT_EXT "$key" "pdf")"
  local gene_scale="$(cfg GENE_SCALE "$key" "")"
  local run_combined="$(cfg RUN_COMBINED "$key" "yes")"

  local title_comb
  if [[ -n "${TITLE_COMB[$key]:-}" ]]; then
    title_comb="${TITLE_COMB[$key]}"
  else
    title_comb="${title}"
  fi

  mkdir -p "${outdir}"

  echo "[i] ${key}: preprocess"
  python 0_preprocess_tiles.py \
    --controls-dir ../data/control_activity \
    --treatments-dir ../data/treatments_effect \
    --map "${map_file}" \
    --out "${outdir}/${key,,}_act_R.csv" \
    --treat-hits-out "${outdir}/treatment_hits.csv" \
    --control-lfc "${CONTROL_LFC}" \
    --treat-lfc "${TREAT_LFC}" \
    --dna-threshold "${DNA_THRESHOLD}" \
    --require-neighbor-active \
    ${VERBOSE}

  echo "[i] ${key}: pull sequences"
  python 0_pull_seq.py \
    --ids "${outdir}/treatment_hits.csv" \
    --seqs ../data/all_tiles_seqs.tsv \
    --out "${outdir}/significant_tiles_sequences.csv"

  echo "[i] ${key}: make BED"
  Rscript 1_generate_bed.r \
    --excel="${outdir}/significant_tiles_sequences.csv" \
    --seq-col=sequence \
    --id-col=ID \
    --accession="${accession}" \
    --outdir="${outdir}/" \
    --prefix="${prefix}"

  # Build common plot args
  local matrix_file="${outdir}/${key,,}_act_R.csv"
  local out_base="${outdir}/${key}_mpra_heatmap"
  local plot_args=( --matrix "${matrix_file}"
                    --conditions "${conds}"
                    --controls "${ctls}"
                    --order "${order}"
                    --kb-step "${DEF_KB_STEP}"
                    --title "${title}" )

  # genes bed
  if [[ -n "${genes_bed}" ]]; then
    plot_args+=( --genes-bed "${genes_bed}" )
  fi
  if [[ -n "${genes_anno}" ]]; then
    plot_args+=( --genes-anno "${genes_anno}" )
  fi
  if [[ -n "${gene_scale}" ]]; then
    plot_args+=( --gene-track-scale "${gene_scale}" )
  fi

  local extra_std="${EXTRA_STD[$key]:-}"
  local extra_comb="${EXTRA_COMB[$key]:-}"

  echo "[i] ${key}: plot (standard)"
  if [[ -n "${extra_std}" ]]; then
    # shellcheck disable=SC2206
    extra_split_std=( ${extra_std} )
    Rscript 2_plot_heatmap.r \
      "${plot_args[@]}" \
      --out "${out_base}.${out_ext}" \
      "${extra_split_std[@]}"
  else
    Rscript 2_plot_heatmap.r \
      "${plot_args[@]}" \
      --out "${out_base}.${out_ext}"
  fi

  if [[ "${run_combined}" == "yes" ]]; then
    echo "[i] ${key}: plot (combined)"
    # Use combined title override if provided
    plot_args_comb=( "${plot_args[@]}" )
    # Replace last occurrence of --title value by rebuilding pair
    plot_args_comb=( "--matrix" "${matrix_file}" "--conditions" "${conds}" "--controls" "${ctls}" "--order" "${order}" "--kb-step" "${DEF_KB_STEP}" "--title" "${title_comb}" )

    # Choose extra args for combined
    local combined_extras
    if [[ -n "${extra_comb}" ]]; then
      combined_extras="${extra_comb}"
    else
      combined_extras="${extra_std}"
    fi

    if [[ -n "${genes_bed}" ]]; then
      plot_args_comb+=( --genes-bed "${genes_bed}" )
    fi
    if [[ -n "${genes_anno}" ]]; then
      plot_args_comb+=( --genes-anno "${genes_anno}" )
    fi
    if [[ -n "${gene_scale}" ]]; then
      plot_args_comb+=( --gene-track-scale "${gene_scale}" )
    fi

    if [[ -n "${combined_extras}" ]]; then
      # shellcheck disable=SC2206
      extra_split_comb=( ${combined_extras} )
      Rscript 2_plot_heatmap.r \
        "${plot_args_comb[@]}" \
        --combined \
        --out "${out_base}_combined.${out_ext}" \
        "${extra_split_comb[@]}"
    else
      Rscript 2_plot_heatmap.r \
        "${plot_args_comb[@]}" \
        --combined \
        --out "${out_base}_combined.${out_ext}"
    fi
  fi
}

# ---------------------------
# Execution
# ---------------------------
# Order matters only for readability of logs
VIRUSES=(
  EBV
  HCMV
  KHSV
  VZV
  HSV1
  HSV2
  HHV6B
  HHV7
  HIV1_CH058
  HIV2_ROD
  HIV2_GH
  HTLV1
  HTLV2
  HTLV3
  POLY_BK
  POLY_JC
  POLY_MERK
  HPV1
  HPV2
  HPV5
  HPV6B
  HPV11
  HPV16
  HPV18
  HIV1_rejo
  HPV52
  AD1
  AD3
  AD4
  AD5
  AD7
  AD11
  AD14
  AD37
)

for v in "${VIRUSES[@]}"; do
  run_pipeline_for "${v}"
done
