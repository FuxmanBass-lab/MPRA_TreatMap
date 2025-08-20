#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Usage Example:
#   Rscript code/2_plot_heatmap.r \
#     --matrix=results/matrix.csv \
#     --conditions=WT,KO,Mut \
#     --controls=WT \
#     --order=WT,KO,Mut \
#     --genes-bed=annotations/genes.bed \
#     --genes-anno=annotations/genes_anno.csv \
#     --out=figures/heatmap_ebv.pdf \
#     --title="EBV MPRA Differential Activity" \
#     --kb-step=25000
#     --combine
# ------------------------------------------------------------------------------
# Script: 2_plot_heatmap.r
# Description:
#   EXACT style match to legacy EBV MPRA plot provided by Berkay.
#   - No additional filtering; uses a pre-filtered wide matrix.
#   - Conditions are user-specified; controls are a subset.
#   - Order is user-specified; each condition renders two tracks (sense/antisense).
#   - Controls use the green gradient; treatments use the blue↔red diverging gradient.
#   - Bands, gaps, separators, gene arrows, labels, axes, and theme match legacy code.
#
# Inputs:
#   --matrix        Wide CSV from preprocessing with columns:
#                     bp, log2fc_<COND>, log2fc_<COND>_antisense
#   --conditions    Comma-separated condition labels (in matrix)
#   --controls      Comma-separated subset of conditions (for style only)
#   --order         Comma-separated display order of conditions (default: --conditions)
#   --genes-bed     BED file: chrom,start,end,name,score,strand  (name must match gene_id in anno)
#   --genes-anno    CSV/TSV: either (A) gene_id, IE, E, E.L, L, Latent (0/1)
#                                  or (B) Gene, Immediate Early, Early, Early-Late, Late, Latent (0/1)
#                    (If B, Gene is treated as gene_id; Immediate Early→IE, Early-Late→E.L)
#
# Options:
#   --out           Output image path (png/pdf/svg). If empty, display only.
#   --title         Plot title (default: "MPRA Differential Activity – EBV")
#   --width         Width in inches (default: 14)
#   --height        Height in inches (default: 8)
#   --dpi           Resolution when saving (default: 300)
#   --kb-step       X-axis major tick step in bp (default: 25000)
#   --seed          Random seed (default: 42)
#   --gene-gap      Extra vertical gap between heatmap bands and gene tracks (default: 1.0 units)
#   --combine       If set, combine sense/antisense per condition using max(log2fc)
#   Note: If --kb-step is 0 or too large for the genome span, the script auto-selects pretty breaks.
# ------------------------------------------------------------------------------

# -------------------------
# Minimal flag parser
# -------------------------
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  kv <- list()
  i <- 1
  n <- length(args)
  while (i <= n) {
    a <- args[[i]]
    if (!grepl('^--', a)) { i <- i + 1; next }  # ignore stray tokens

    # Case 1: --key=value
    if (grepl('^--[A-Za-z0-9_.-]+=.*', a)) {
      sp <- strsplit(sub('^--', '', a), '=', fixed = TRUE)[[1]]
      key <- sp[1]
      val <- paste(sp[-1], collapse='=')  # allow '=' in value
      kv[[key]] <- val
      i <- i + 1
      next
    }

    # Strip leading dashes to get the key
    key <- sub('^--', '', a)

    # Case 2: --flag (boolean) or --key value
    if (i == n) {
      # last token, treat as boolean TRUE
      kv[[key]] <- TRUE
      i <- i + 1
    } else {
      nxt <- args[[i+1]]
      if (!grepl('^--', nxt)) {
        # --key value
        kv[[key]] <- nxt
        i <- i + 2
      } else {
        # boolean flag
        kv[[key]] <- TRUE
        i <- i + 1
      }
    }
  }

  # required
  if (is.null(kv$matrix))      stop('--matrix is required')
  if (is.null(kv$conditions))  stop('--conditions is required (comma-separated)')

  # defaults
  if (is.null(kv$controls))    kv$controls <- ''
  if (is.null(kv$order))       kv$order <- kv$conditions
  if (is.null(kv$`genes-bed`)) kv$`genes-bed` <- ''
  if (is.null(kv$`genes-anno`))kv$`genes-anno` <- ''
  if (is.null(kv$out))         kv$out <- ''
  if (is.null(kv$title))       kv$title <- 'MPRA Differential Activity – EBV'
  if (is.null(kv$width))       kv$width <- '14'
  if (is.null(kv$height))      kv$height <- '8'
  if (is.null(kv$dpi))         kv$dpi <- '300'
  if (is.null(kv$`kb-step`))   kv$`kb-step` <- '25000'
  if (is.null(kv$seed))        kv$seed <- '42'
  if (is.null(kv$`gene-gap`))   kv$`gene-gap` <- '2.0'
  if (is.null(kv$combine))     kv$combine <- FALSE

  # coerce types
  kv$width   <- as.numeric(kv$width)
  kv$height  <- as.numeric(kv$height)
  kv$dpi     <- as.numeric(kv$dpi)
  kv$kb_step <- as.numeric(kv$`kb-step`)
  kv$seed    <- as.integer(kv$seed)
  kv$gene_gap <- as.numeric(kv$`gene-gap`)
  kv$combine <- tolower(as.character(kv$combine)) %in% c('1','true','t','yes','y')

  # vectors
  kv$conditions_vec <- trimws(strsplit(as.character(kv$conditions), ',', fixed = TRUE)[[1]])
  kv$order_vec      <- trimws(strsplit(as.character(kv$order), ',', fixed = TRUE)[[1]])
  kv$controls_vec   <- if (nzchar(kv$controls)) trimws(strsplit(as.character(kv$controls), ',', fixed = TRUE)[[1]]) else character(0)

  return(kv)
}

opt <- parse_args()
set.seed(opt$seed)

# -------------------------
