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
# Libraries
# -------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(ggnewscale)
  library(gggenes)
  library(ggrepel)
  library(grid)   # unit()
  library(readr)
  library(stringr)
  library(tools)
})

# -------------------------
# Load matrix, derive columns
# -------------------------
mat <- read.csv(opt$matrix, stringsAsFactors = FALSE, check.names = FALSE)
if (!"bp" %in% names(mat)) stop("Input matrix must have a 'bp' column")

# Respect explicit order
if (!setequal(opt$order_vec, opt$conditions_vec)) {
  warning("--order and --conditions differ; intersecting and using --order sequence")
}
ord <- intersect(opt$order_vec, opt$conditions_vec)

# Build column order: sense/antisense or combined
if (!opt$combine) {
  # Build column order: each condition → sense, antisense
  column_order <- as.vector(rbind(paste0("log2fc_", ord), paste0("log2fc_", ord, "_antisense")))
  missing <- setdiff(column_order, names(mat))
  if (length(missing) > 0) stop("Missing columns in --matrix: ", paste(missing, collapse=", "))
  # Controls columns (for styling only)
  ctrl_set <- intersect(ord, opt$controls_vec)
  control_columns <- as.vector(rbind(paste0("log2fc_", ctrl_set), paste0("log2fc_", ctrl_set, "_antisense")))
} else {
  # Combine sense and antisense using max(log2fc) per bp
  combined_cols <- character(0)
  for (c in ord) {
    s_col <- paste0("log2fc_", c)
    a_col <- paste0("log2fc_", c, "_antisense")
    miss <- setdiff(c(s_col, a_col), names(mat))
    if (length(miss) > 0) stop("Missing columns in --matrix for combine: ", paste(miss, collapse=", "))
    new_col <- paste0("log2fc_", c, "_combined")
    v1 <- suppressWarnings(as.numeric(mat[[s_col]]))
    v2 <- suppressWarnings(as.numeric(mat[[a_col]]))
    comb <- ifelse(is.na(v1) & is.na(v2), NA_real_, pmax(v1, v2, na.rm = TRUE))
    mat[[new_col]] <- comb
    combined_cols <- c(combined_cols, new_col)
  }
  column_order <- combined_cols
  ctrl_set <- intersect(ord, opt$controls_vec)
  control_columns <- paste0("log2fc_", ctrl_set, "_combined")
}

# Long format and classify type
data_long <- mat %>%
  select(bp, all_of(column_order)) %>%
  pivot_longer(-bp, names_to = "condition", values_to = "log2fc") %>%
  mutate(
    condition = factor(condition, levels = column_order),
    type      = ifelse(condition %in% control_columns, "control", "treatment"),
    log2fc    = as.numeric(log2fc)
  ) %>%
  filter(!is.na(log2fc))

# -------------------------
# Build band layout (exact legacy-style)
# -------------------------
band_height <- 1
gap_within  <- 0.15   # between sense/antisense pairs
gap_between <- 0.5    # between control block and treatment block

# Place the big gap *after the last control band in the actual layout order*,
# even if controls are interleaved or ordered arbitrarily.
ctrl_band_indices <- match(control_columns, column_order)
ctrl_band_indices <- ctrl_band_indices[!is.na(ctrl_band_indices)]
last_ctrl_idx <- if (length(ctrl_band_indices) > 0) max(ctrl_band_indices) else 0

pos <- 0; ymid <- c()
for (i in seq_along(column_order)) {
  ymid <- c(ymid, pos)
  pos  <- pos - band_height

  # Decide the gap after this band
  if (i < length(column_order)) {
    if (last_ctrl_idx > 0 && i == last_ctrl_idx) {
      # Big gap between control block and treatment block
      pos <- pos - gap_between
    } else {
      # Within-block gap
      if (!opt$combine) {
        # In split mode, apply the within-gap after every sense/antisense pair
        if (i %% 2 == 0) pos <- pos - gap_within
      } else {
        # In combined mode, apply the same within-gap after every band
        pos <- pos - gap_within
      }
    }
  }
}
layout <- tibble(
  condition = factor(column_order, levels = column_order),
  ymid      = ymid
) %>% mutate(
  ymin  = ymid - band_height/2,
  ymax  = ymid + band_height/2,
  start = min(mat$bp, na.rm=TRUE),
  end   = max(mat$bp, na.rm=TRUE),
  # Pretty labels matching legacy text, updated for combine
  label = if (!opt$combine) {
    ifelse(grepl("_antisense$", condition),
           paste(sub("^log2fc_", "", sub("_antisense$", "", condition)),
                ifelse(condition %in% control_columns, "antisense", "+ treatment antisense")),
           paste(sub("^log2fc_", "", condition),
                ifelse(condition %in% control_columns, "sense", "+ treatment sense")))
  } else {
    base <- sub("^log2fc_", "", sub("_combined$", "", condition))
    ifelse(condition %in% control_columns, base, paste(base, "+ treatment"))
  }
)
lowest_y <- min(layout$ymin)

# Separator line y between control block and first treatment band
midGap <- NA_real_
if (last_ctrl_idx > 0 && last_ctrl_idx < length(column_order)) {
  y_top  <- layout$ymax[last_ctrl_idx]
  y_bot  <- layout$ymin[last_ctrl_idx + 1]
  midGap <- (y_top + y_bot)/2
}

# -------------------------
# Genes: from BED + annotations (stage & latent)
# -------------------------
manual <- NULL
if (nzchar(opt$`genes-bed`)) {
  bed <- tryCatch({ read_tsv(opt$`genes-bed`, col_names = FALSE, show_col_types = FALSE, progress = FALSE) }, error=function(e) NULL)
  if (!is.null(bed) && ncol(bed) >= 6) {
    names(bed)[1:6] <- c("chrom","start","end","name","score","strand")
    manual <- bed %>% transmute(
      gene_id = name,
      strand  = strand,
      start   = as.numeric(start),
      end     = as.numeric(end)
    )
  } else if (!is.null(bed) && ncol(bed) >= 4) { # fallback if no score/strand
    names(bed)[1:4] <- c("chrom","start","end","name")
    manual <- bed %>% mutate(strand = "+") %>% transmute(
      gene_id = name, strand = strand,
      start = as.numeric(start), end = as.numeric(end)
    )
  }
}

# -------------------------
# Stage annotations (robust to column name variants and x/0/1)
# -------------------------
annots <- NULL
has_annots <- FALSE
if (nzchar(opt$`genes-anno`)) {
  ext <- tolower(file_ext(opt$`genes-anno`))
  ga <- tryCatch({
    if (ext %in% c("tsv","tab","txt")) readr::read_tsv(opt$`genes-anno`, show_col_types = FALSE, progress = FALSE)
    else readr::read_csv(opt$`genes-anno`, show_col_types = FALSE, progress = FALSE)
  }, error = function(e) NULL)

  if (!is.null(ga)) {
    # Normalize column names: trim, lower, collapse spaces/underscores/hyphens
    norm <- function(x) {
      x <- gsub("[\t\r\n]", "", x)
      x <- trimws(x)
      x <- tolower(x)
      x <- gsub("[ _]+", " ", x)
      x
    }
    nms_raw <- names(ga)
    nms <- norm(nms_raw)

    # Map known variants to canonical keys
    key_map <- c(
      "gene" = "gene_id",
      "gene id" = "gene_id",
      "gene_id" = "gene_id",
      "name" = "gene_id",
      "immediate early" = "ie",
      "ie" = "ie",
      "early" = "e",
      "e" = "e",
      "early-late" = "e.l",
      "early late" = "e.l",
      "e-l" = "e.l",
      "e.l" = "e.l",
      "late" = "l",
      "l" = "l",
      "latent" = "latent"
    )
    canon <- vapply(nms, function(x) if (x %in% names(key_map)) key_map[[x]] else x, character(1))
    names(ga) <- canon

    # Ensure required columns exist (create if missing)
    if (!"gene_id" %in% names(ga)) {
      stop("--genes-anno is missing a gene identifier column (accepted: Gene, gene_id, name)")
    }
    for (cn in c("ie","e","e.l","l","latent")) if (!cn %in% names(ga)) ga[[cn]] <- NA

    # Coerce flags to 0/1 from 1/0/x/blank
    coerce01 <- function(v) {
      if (is.null(v)) return(0L)
      v <- as.character(v)
      v[is.na(v)] <- "0"
      v <- trimws(tolower(v))
      v <- ifelse(v %in% c("x","1","true","t","yes","y"), "1", v)
      suppressWarnings(as.integer(ifelse(v == "1", 1L, 0L)))
    }
    ga$IE     <- coerce01(ga$ie)
    ga$E      <- coerce01(ga$e)
    ga$`E.L`  <- coerce01(ga$`e.l`)
    ga$L      <- coerce01(ga$l)
    ga$Latent <- coerce01(ga$latent)

    annots <- ga %>% dplyr::transmute(
      gene_id,
      stage = dplyr::case_when(
        IE == 1 ~ "IE",
        E  == 1 ~ "E",
        `E.L` == 1 ~ "E-L",
        L  == 1 ~ "L",
        Latent == 1 ~ NA_character_,  # latent handled via is_latent
        TRUE ~ NA_character_
      ),
      is_latent = Latent == 1
    )
    # Flag used later to control legends/labels when annotations are absent
    has_annots <- !is.null(annots)
  }
}

# Build sense/antisense gene rows & compute positions like legacy
genes <- NULL
if (!is.null(manual)) {
  gene_top_base <- 1.0 + opt$gene_gap
  gene_bottom_base <- 3.0 + opt$gene_gap
  sense_df <- manual %>% filter(strand == "+") %>%
    mutate(xmid = (start + end)/2,
           subrow = (floor(xmid/1000) %% 3) + 1,
           track  = lowest_y - gene_top_base - (subrow-1)*1.5,
           strand = "+")
  antisense_df <- manual %>% filter(strand == "-") %>%
    mutate(xmid = (start + end)/2,
           subrow = (floor(xmid/1000) %% 4) + 1,
           track  = lowest_y - gene_bottom_base - (subrow-1)*1.5,
           strand = "-")
  genes <- bind_rows(sense_df, antisense_df)
  if (!is.null(annots)) {
    # Keep ONLY genes that appear in the annotations table
    genes <- genes %>% inner_join(annots, by = "gene_id")
  } else {
    # No annotations provided → keep genes but with NA stage/latent
    genes <- genes %>% mutate(stage = NA_character_, is_latent = NA)
  }
  genes <- genes %>% mutate(
    nudge_offset = ifelse(strand == "+", 0.5, -0.5),
    hjust        = ifelse(strand == "+", 0, 1)
  )
}

# -------------------------
# Plot (exact styling)
# -------------------------
set.seed(opt$seed)

p <- ggplot() +
  # (A) grey bands
  geom_rect(data = layout, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax),
            fill = "#F7F7F7") +

  # (B) control bars (green gradient)
  geom_segment(
    data = left_join(filter(data_long, type == "control"), layout, by = "condition"),
    aes(x = bp, xend = bp, y = ymin, yend = ymax, color = log2fc),
    linewidth = 0.9, lineend = "butt"
  ) +
  scale_color_gradientn(
    colors = c("#edf8e9",
"#bae4b3","#74c476","#31a354","#006d2c"),
    values = scales::rescale(c(3,3.5,4.5,5)),
    limits = c(1,6), oob = scales::squish,
    na.value = "white",
    name = expression(log[2]*"FC (control)")
  ) +
  ggnewscale::new_scale_color() +

  # (C) treatment bars (diverging)
  geom_segment(
    data = left_join(filter(data_long, type == "treatment"), layout, by = "condition"),
    aes(x = bp, xend = bp, y = ymin, yend = ymax, color = log2fc),
    linewidth = 0.9, lineend = "butt"
  ) +
  scale_color_gradientn(
    colors = c("#053061","#2166ac","#4393c3","white",
               "#fddbc7","#f4a582","#d6604d","#b2182b"),
    values = scales::rescale(c(-3,-2,-1,0,0.5,1,2,3)),
    breaks = c(-3,-2,-1,0,1,2,3),
    limits = c(-4,4), oob = scales::squish,
    name   = expression(log[2]*"FC (treatment)")
  ) +
  ggnewscale::new_scale_color()

# (D) gene arrows + labels (if available)
if (!is.null(genes)) {
  p <- p +
    gggenes::geom_gene_arrow(
      data = genes,
      aes(xmin = start, xmax = end, y = track,
          forward = strand == "+", fill = stage, linetype = is_latent),
      arrowhead_height = unit(3,"mm"),
      color = "black", size = 0.3
    )

  if (has_annots) {
    p <- p +
      scale_fill_manual(
        values   = c(IE="#E41A1C", E="darkorange", `E-L`="#4DAF4A", L="#984EA3"),
        na.value = "gray80", name="Stage"
      ) +
      scale_linetype_manual(
        values   = c(`FALSE`="solid", `TRUE`="dashed"),
        na.value = "solid", name="Latent?"
      )
  } else {
    p <- p +
      scale_fill_manual(
        values   = c(IE="#E41A1C", E="darkorange", `E-L`="#4DAF4A", L="#984EA3"),
        na.value = "gray80", guide = "none"
      ) +
      scale_linetype_manual(
        values   = c(`FALSE`="solid", `TRUE`="dashed"),
        na.value = "solid", guide = "none"
      )
  }

  p <- p +
    ggrepel::geom_text_repel(
      data = genes,
      aes(x = xmid, y = track, label = gene_id, color = stage, hjust = hjust),
      direction     = "both",
      nudge_y       = genes$nudge_offset,
      segment.size  = 0.2,
      segment.color = "grey50",
      box.padding   = 0.4,
      point.padding = 0.2,
      force         = 1,
      max.overlaps  = Inf,
      size          = 3,
      inherit.aes   = FALSE
    )

  if (has_annots) {
    p <- p + scale_color_manual(
      name   = "Stage",
      values = c(IE="#E41A1C", E="darkorange", `E-L`="#4DAF4A", L="#984EA3"),
      na.value = "black"
    )
  } else {
    p <- p + scale_color_manual(
      values = c(IE="#E41A1C", E="darkorange", `E-L`="#4DAF4A", L="#984EA3"),
      na.value = "black", guide = "none"
    )
  }
}

# (E) separator & boundaries
if (!is.na(midGap)) {
  p <- p + geom_hline(yintercept = midGap, color = "white", linewidth = 0.2)
}
p <- p +
  geom_vline(xintercept = c(min(mat$bp, na.rm=TRUE), max(mat$bp, na.rm=TRUE)),
             color = "gray90", linewidth = 0.3)

# (F) axes
# Compute robust x breaks: if kb_step is too large for the span (or <=0),
# fall back to pretty() to ensure visible tick marks.
x_min <- floor(min(mat$bp, na.rm=TRUE)/1000)*1000
x_max <- ceiling(max(mat$bp, na.rm=TRUE)/1000)*1000
if (is.na(opt$kb_step) || opt$kb_step <= 0) {
  x_breaks <- scales::pretty_breaks(n = 6)(c(x_min, x_max))
} else {
  x_breaks <- seq(x_min, x_max, by = opt$kb_step)
  if (length(x_breaks) < 3) {
    x_breaks <- scales::pretty_breaks(n = 6)(c(x_min, x_max))
  }
}

p <- p +
  scale_y_continuous(
    breaks = layout$ymid, labels = layout$label,
    expand = c(0,0),
    limits = c(if (!is.null(genes)) min(genes$track + genes$nudge_offset) - 0.2 else min(layout$ymin) - 0.5,
               max(layout$ymax) + 0.5)
  ) +
  scale_x_continuous(
    name = "Genome Position (kb)",
    breaks = x_breaks,
    labels = function(x) x/1000
  )

# (G) theme & title (exact match)
p <- p +
  labs(title = opt$title, y = NULL) +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y     = element_text(size = 10),
    axis.title.x    = element_text(margin = margin(t = 10)),
    plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.grid      = element_blank()
  )

print(p)

if (nzchar(opt$out)) {
  ggsave(filename = opt$out, plot = p, width = opt$width, height = opt$height, dpi = opt$dpi, units = "in")
}
