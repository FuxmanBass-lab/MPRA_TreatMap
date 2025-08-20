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
# Fetch/reference files (cached to --outdir)
# ────────────────────────────────────────────────────────────────────────────────

fa_path   <- file.path(opt$outdir, paste0(opt$accession, ".fasta"))
fa_fai    <- paste0(fa_path, ".fai")
gff3_path <- file.path(opt$outdir, paste0(opt$accession, ".gff3"))

fetch_if_missing <- function(path, fetch_fun, label) {
  if (!file.exists(path)) {
    cat(sprintf("[%s] Downloading to %s ...\n", label, path))
    txt <- fetch_fun()
    writeLines(txt, path)
  } else {
    cat(sprintf("[%s] Using cached %s\n", label, path))
  }
}

# FASTA
fetch_if_missing(
  fa_path,
  function() rentrez::entrez_fetch(db = "nuccore", id = opt$accession, rettype = "fasta", retmode = "text"),
  "FASTA"
)

# index FASTA if missing
if (!file.exists(fa_fai)) {
  cat(sprintf("[FASTA] Indexing %s ...\n", fa_path))
  Rsamtools::indexFa(fa_path)
}

# GFF3 (for annotation overlap)
fetch_if_missing(
  gff3_path,
  function() rentrez::entrez_fetch(db = "nuccore", id = opt$accession, rettype = "gff3", retmode = "text"),
  "GFF3"
)

# ────────────────────────────────────────────────────────────────────────────────
# Load genome and queries
# ────────────────────────────────────────────────────────────────────────────────

ebv_dna <- Biostrings::readDNAStringSet(fa_path)
if (length(ebv_dna) < 1) {
  cat("ERROR: FASTA contains no sequences.\n")
  quit(status = 1)
}
contig_id <- strsplit(names(ebv_dna)[1], " ")[[1]][1]
ebv_seq   <- ebv_dna[[1]]
cat(sprintf("Loaded contig %s (length = %d bp)\n", contig_id, length(ebv_seq)))

# Read Excel or CSV
if (!file.exists(opt$excel)) {
  cat(sprintf("ERROR: Input file not found: %s\n", opt$excel))
  quit(status = 1)
}

file_ext <- tolower(tools::file_ext(opt$excel))
if (file_ext %in% c("xlsx", "xls")) {
  df <- readxl::read_excel(opt$excel, sheet = opt$sheet)
} else if (file_ext == "csv") {
  df <- read.csv(opt$excel, stringsAsFactors = FALSE)
} else {
  cat(sprintf("ERROR: Unsupported input file extension '%s'. Supported extensions are .xlsx, .xls, and .csv\n", file_ext))
  quit(status = 1)
}

# tolerate different casings
colnames(df) <- tolower(colnames(df))
seq_col <- tolower(opt$`seq-col`)
id_col  <- tolower(opt$`id-col`)

if (!all(c(seq_col, id_col) %in% colnames(df))) {
  cat(sprintf("ERROR: Columns '%s' and '%s' must exist in %s (have: %s)\n",
              seq_col, id_col, opt$excel, paste(colnames(df), collapse = ", ")))
  quit(status = 1)
}

df <- df[!is.na(df[[seq_col]]) & df[[seq_col]] != "", , drop = FALSE]
if (nrow(df) == 0) {
  cat("ERROR: No non-empty sequences found in the input file.\n")
  quit(status = 1)
}

queries <- as.character(df[[seq_col]])
ids     <- as.character(df[[id_col]])
names(queries) <- ids

# ────────────────────────────────────────────────────────────────────────────────
# Exact-match mapping (forward + reverse complement) with a progress bar
# ────────────────────────────────────────────────────────────────────────────────

hits_list <- vector("list", length(queries))

pb_total <- length(queries)
cat(sprintf("Mapping %d query sequences ...\n", pb_total))
pb <- txtProgressBar(min = 0, max = pb_total, style = 3)

for (i in seq_along(queries)) {
  query_s <- queries[i]
  query_id <- names(queries)[i]
  gr_all <- GenomicRanges::GRanges()

  if (nchar(query_s) == 0) {
    setTxtProgressBar(pb, i); next
  }

  q <- Biostrings::DNAString(query_s)

  # Forward
  fwd <- Biostrings::matchPattern(q, ebv_seq, max.mismatch = 0)
  if (length(fwd) > 0) {
    gr_all <- c(gr_all, GenomicRanges::GRanges(
      seqnames = contig_id,
      ranges   = IRanges::ranges(fwd),
      strand   = "+",
      query_id = query_id
    ))
  }

  # Reverse complement
  rev <- Biostrings::matchPattern(Biostrings::reverseComplement(q), ebv_seq, max.mismatch = 0)
  if (length(rev) > 0) {
    gr_all <- c(gr_all, GenomicRanges::GRanges(
      seqnames = contig_id,
      ranges   = IRanges::ranges(rev),
      strand   = "-",
      query_id = query_id
    ))
  }

  hits_list[[i]] <- gr_all
  setTxtProgressBar(pb, i)
}
close(pb)

all_hits <- do.call(c, hits_list)
if (length(all_hits) == 0) {
  cat("No exact matches were found in the reference.\n")
  quit(status = 0)
}

# name + score for BED
mcols(all_hits)$name  <- all_hits$query_id
mcols(all_hits)$score <- as.integer(1000L)
# keep only name/score
mcols(all_hits) <- mcols(all_hits)[, c("name", "score"), drop = FALSE]

# write BED of mapped regions
mapped_bed <- file.path(opt$outdir, paste0(opt$prefix, "_mapped_regions.bed"))
rtracklayer::export(all_hits, mapped_bed, format = "bed")
cat(sprintf("Wrote mapped regions to %s\n", normalizePath(mapped_bed)))

# ────────────────────────────────────────────────────────────────────────────────
# Nearby gene/CDS features within ±window
# ────────────────────────────────────────────────────────────────────────────────

genome_len <- length(ebv_seq)
# clamp to genome
ext <- GenomicRanges::GRanges(
  seqnames = GenomicRanges::seqnames(all_hits),
  ranges = IRanges::IRanges(
    start = pmax(GenomicRanges::start(all_hits) - opt$window, 1L),
    end   = pmin(GenomicRanges::end(all_hits)   + opt$window, genome_len)
  ),
  strand = "*"
)

# import GFF3
cat(sprintf("Importing annotations from %s ...\n", gff3_path))
gff_gr <- rtracklayer::import(gff3_path, format = "gff3")

# keep gene/CDS entries
feature_types <- c("gene", "CDS")
if (!"type" %in% colnames(mcols(gff_gr))) {
  # Some NCBI GFF3s use 'type' in the GRanges top-level slot already
}

gene_feats <- gff_gr[gff_gr$type %in% feature_types]

# overlaps (robust to older GenomicRanges versions)
ov <- GenomicRanges::findOverlaps(gene_feats, ext, ignore.strand = TRUE)
if (length(ov) > 0) {
  nearby_genes <- gene_feats[unique(S4Vectors::queryHits(ov))]
} else {
  nearby_genes <- gene_feats[0]
}

# Give a friendly name column if possible
nm <- NULL
if ("Name" %in% colnames(mcols(nearby_genes))) nm <- mcols(nearby_genes)$Name
if (is.null(nm) && "gene" %in% colnames(mcols(nearby_genes))) nm <- mcols(nearby_genes)$gene
if (is.null(nm)) nm <- paste0(nearby_genes$type, "_", seq_along(nearby_genes))

mcols(nearby_genes)$name  <- nm
mcols(nearby_genes)$score <- as.integer(0L)
# keep only name/score
mcols(nearby_genes) <- mcols(nearby_genes)[, c("name", "score"), drop = FALSE]

nearby_bed <- file.path(opt$outdir, paste0(opt$prefix, "_nearby_genes.bed"))
rtracklayer::export(nearby_genes, nearby_bed, format = "bed")
cat(sprintf("Wrote nearby gene/CDS features to %s\n", normalizePath(nearby_bed)))

# ────────────────────────────────────────────────────────────────────────────────
# Done
# ────────────────────────────────────────────────────────────────────────────────
cat("\nAll done.\n")
