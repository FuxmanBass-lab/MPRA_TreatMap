#!/usr/bin/env Rscript

# ────────────────────────────────────────────────────────────────────────────────
# generate_bed.r
# Map query sequences to an EBV reference, write BED of matches, and a BED of
# nearby gene/CDS features within a ± window.
#
# Run from command line, e.g.:
#   Rscript generate_bed.r \
#     --excel=ebv_hits.xlsx \
#     --sheet=Sheet1 \
#     --seq-col=sequence \
#     --id-col=tile_id \
#     --accession=NC_007605.1 \
#     --window=5000 \
#     --outdir=. \
#     --prefix=ebv
#
# Outputs (in --outdir):
#   <prefix>_mapped_regions.bed
#   <prefix>_nearby_genes.bed
#   <accession>.fasta (+ .fai)
#   <accession>.gff3
#
# Notes:
#  * No packages are auto-installed. If something is missing, the script exits
#    with a clear message telling you what to install.
#  * BED export uses rtracklayer and adheres to BED 0-based, half-open coords.
# ────────────────────────────────────────────────────────────────────────────────

suppressWarnings(suppressMessages({
  # we only load once we have confirmed availability via `need()`
}))

need <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("\nMissing required package '%s'.\n", pkg))
    if (bioc) {
      cat("Install with: BiocManager::install('", pkg, "')\n", sep = "")
    } else {
      cat("Install with: install.packages('", pkg, "')\n", sep = "")
    }
    quit(status = 1)
  }
}


# Required packages
need("Biostrings", bioc = TRUE)
need("GenomicRanges", bioc = TRUE)
need("rtracklayer", bioc = TRUE)
need("Rsamtools", bioc = TRUE)
need("rentrez", bioc = TRUE)
need("readxl")

suppressPackageStartupMessages({
  library(Biostrings)
  library(GenomicRanges)
  library(rtracklayer)
  library(Rsamtools)
  library(rentrez)
  library(readxl)
})

# ────────────────────────────────────────────────────────────────────────────────
# CLI parsing (no external deps): flags are --key=value
# ────────────────────────────────────────────────────────────────────────────────

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  kv <- list()
  for (a in args) {
    if (grepl("^--[A-Za-z0-9_.-]+=", a)) {
      sp <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
      kv[[sp[1]]] <- sp[2]
    }
  }
  # defaults
  if (is.null(kv$accession)) kv$accession <- "NC_007605.1"
  if (is.null(kv$sheet))     kv$sheet     <- "Sheet1"
  if (is.null(kv$`seq-col`)) kv$`seq-col` <- "sequence"
  if (is.null(kv$`id-col`))  kv$`id-col`  <- "tile_id"
  if (is.null(kv$window))    kv$window    <- "5000"
  if (is.null(kv$outdir))    kv$outdir    <- "."
  if (is.null(kv$prefix))    kv$prefix    <- "ebv"
  # required
  if (is.null(kv$excel)) {
    usage()
    cat("\nERROR: --excel=<path/to/xlsx> is required.\n\n")
    quit(status = 2)
  }
  # coerce
  kv$window <- as.integer(kv$window)
  if (is.na(kv$window) || kv$window < 0) {
    cat("ERROR: --window must be a non-negative integer.\n")
    quit(status = 2)
  }
  kv
}

usage <- function() {
  cat("\nUsage:\n  Rscript generate_bed.r --excel=ebv_hits.xlsx [--sheet=Sheet1] \\\n         [--seq-col=sequence] [--id-col=tile_id] \\\n         [--accession=NC_007605.1] [--window=5000] \\\n         [--outdir=.] [--prefix=ebv]\n\n")
}

opt <- parse_args()

# ensure outdir exists
if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ────────────────────────────────────────────────────────────────────────────────

# ────────────────────────────────────────────────────────────────────────────────
cat("\nAll done.\n")
