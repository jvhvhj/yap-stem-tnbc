#!/usr/bin/env Rscript
# Purpose: Figure 4 orthogonal program transfer
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 5; orthogonal transfer in Fig. 4.

options(stringsAsFactors = FALSE, scipen = 999, width = 240, encoding = "UTF-8")
set.seed(20260810)

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "Figure4F_GSE180286_program_transfer")

.libPaths(c(
  file.path(root, "0703_rebuild/R_library"),
  file.path(root, "0710_final_evidence_rebuild/vendor"),
  ".software/r-library",
  .libPaths()
))

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggalluvial)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(svglite)
  library(Cairo)
  library(digest)
})

read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA", fileEncoding = "UTF-8")
sha256 <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)

paths <- list(
  transfer = file.path(out_dir, "Figure4F_program_gene_transfer.tsv"),
  waterfall = file.path(out_dir, "Figure4F_gene_effect_waterfall.tsv"),
  enrichment = file.path(out_dir, "Figure4F_signature_enrichment.tsv"),
  audit = file.path(out_dir, "Figure4F_transfer_audit.tsv"),
  provenance = file.path(out_dir, "GSE180286_sample_patient_provenance_audit.tsv"),
  palette = file.path(root, "project_semantic_palette.tsv"),
  style_mapping = file.path(out_dir, "panel_style_mapping.tsv"),
  reference_mapping = file.path(out_dir, "Figure4F_reference_mapping.md"),
  legend = file.path(out_dir, "Figure4F_updated_legend.md"),
  data_script = file.path(out_dir, "scripts/01_prepare_Figure4F_sources.R"),
  plot_script = file.path(out_dir, "scripts/02_plot_Figure4F.R")
)
missing <- names(paths)[!vapply(paths, file.exists, logical(1))]
if (length(missing)) stop("Missing required Figure 4F input(s): ", paste(missing, collapse = ", "))
input_hash_before <- vapply(paths, sha256, character(1))

transfer <- read_tsv(paths$transfer)
waterfall <- read_tsv(paths$waterfall)
enrich <- read_tsv(paths$enrichment)
audit <- read_tsv(paths$audit)
provenance <- read_tsv(paths$provenance)
palette <- read_tsv(paths$palette)

audit_value <- function(field) {
  z <- audit$value[audit$audit_field == field]
  if (length(z) != 1L) stop("Audit field missing or duplicated: ", field)
  z
}
get_col <- function(token, fallback) {
  z <- palette$hex[match(token, palette$token)]
  ifelse(length(z) == 1L && !is.na(z), z, fallback)
}

COL_WARM <- get_col("Positive", "#D7605C")
COL_COOL <- get_col("Negative", "#2F6DB3")
COL_NEUTRAL <- get_col("Neutral", "#F7F7F7")
COL_UNAVAILABLE <- "#AEB5BA"
COL_TEXT <- "#3D4650"
COL_MUTED <- "#6D7680"
COL_GRID <- "#DDE1E4"
FONT <- "Arial"

stopifnot(
  nrow(transfer) == 146L,
  sum(transfer$detectable_in_GSE180286) == 141L,
  nrow(waterfall) == 141L,
  sum(waterfall$expected_direction_preserved) == 132L,
  sum(!waterfall$expected_direction_preserved) == 9L,
  sum(waterfall$GSE180286_High_vs_Other_effect > 0) == 132L,
  sum(waterfall$GSE180286_High_vs_Other_effect == 0) == 1L,
  sum(waterfall$GSE180286_High_vs_Other_effect < 0) == 8L,
  nrow(enrich) == as.integer(audit_value("GSE180286_ranking_genes")),
  length(unique(enrich$enrichment_NES)) == 1L,
  length(unique(enrich$enrichment_BH_FDR)) == 1L,
  nrow(provenance[provenance$included_in_transfer_effect, , drop = FALSE]) == 2L,
  all(provenance$patient[provenance$included_in_transfer_effect] == "P5")
)

theme_pub <- theme_classic(base_family = FONT, base_size = 7.2) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(size = 7.7, face = "bold", hjust = 0, colour = COL_TEXT, margin = margin(b = 3)),
    axis.title = element_text(size = 6.7, colour = COL_TEXT),
    axis.text = element_text(size = 6.3, colour = COL_TEXT),
    axis.line = element_line(colour = COL_TEXT, linewidth = 0.35),
    axis.ticks = element_line(colour = COL_TEXT, linewidth = 0.3),
    axis.ticks.length = grid::unit(1.1, "mm"),
    legend.title = element_text(size = 6.4, face = "bold", colour = COL_TEXT),
    legend.text = element_text(size = 6.2, colour = COL_TEXT),
    plot.margin = margin(3, 4, 3, 4)
  )

## ------------------------------------------------------------------
## F1: count-proportional horizontal Sankey.
## ------------------------------------------------------------------

make_ribbon <- function(id, x0, x1, lo0, hi0, lo1, hi1, fill, n = 120L) {
  tt <- seq(0, 1, length.out = n)
  ee <- 3 * tt^2 - 2 * tt^3
  xx <- x0 + (x1 - x0) * tt
  top <- hi0 + (hi1 - hi0) * ee
  bottom <- lo0 + (lo1 - lo0) * ee
  data.frame(
    id = id,
    x = c(xx, rev(xx)),
    y = c(top, rev(bottom)),
    fill = fill,
    stringsAsFactors = FALSE
  )
}

flow_ribbons <- rbind(
  make_ribbon("measurable", 0.06, 1.00, -68, 73, -70.5, 70.5, COL_WARM),
  make_ribbon("unavailable", 0.06, 1.00, -73, -68, -87, -82, COL_UNAVAILABLE),
  make_ribbon("preserved", 1.06, 2.10, -61.5, 70.5, -66, 66, COL_WARM),
  make_ribbon("not_preserved", 1.06, 2.10, -70.5, -61.5, -83, -74, COL_COOL)
)
flow_nodes <- data.frame(
  xmin = c(0.00, 1.00, 1.00, 2.10, 2.10),
  xmax = c(0.06, 1.06, 1.06, 2.16, 2.16),
  ymin = c(-73, -70.5, -87, -66, -83),
  ymax = c(73, 70.5, -82, 66, -74),
  border = c(COL_TEXT, COL_TEXT, COL_UNAVAILABLE, COL_WARM, COL_COOL),
  stringsAsFactors = FALSE
)

p_flow <- ggplot() +
  geom_polygon(
    data = flow_ribbons,
    aes(x = x, y = y, group = id, fill = I(fill)),
    alpha = 0.76, colour = NA
  ) +
  geom_rect(
    data = flow_nodes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, colour = I(border)),
    fill = "white", linewidth = 0.35
  ) +
  annotate("text", x = 0.03, y = 87, label = "GSE176078\nprogram (146)",
           family = FONT, fontface = "bold", size = 5.7 / 2.845, lineheight = 0.92, colour = COL_TEXT) +
  annotate("text", x = 1.03, y = 86, label = "Measurable in\nGSE180286\n141 (96.6%)",
           family = FONT, fontface = "bold", size = 5.4 / 2.845, lineheight = 0.92, colour = COL_TEXT) +
  annotate("text", x = 2.03, y = 46, label = "Expected direction\npreserved\n132/141 (93.6%)",
           hjust = 1, family = FONT, fontface = "bold", size = 5.4 / 2.845,
           lineheight = 0.91, colour = "white") +
  annotate("text", x = 2.20, y = -78.5, label = "Not preserved /\nnear-zero  9",
           hjust = 0, family = FONT, fontface = "bold", size = 5.2 / 2.845,
           lineheight = 0.91, colour = COL_COOL) +
  annotate("text", x = 1.09, y = -84.5, label = "Unavailable  5",
           hjust = 0, family = FONT, size = 5.1 / 2.845, colour = COL_MUTED) +
  coord_cartesian(xlim = c(-0.10, 3.12), ylim = c(-95, 96), clip = "off") +
  labs(title = "Program retention") +
  theme_void(base_family = FONT) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(size = 7.7, face = "bold", hjust = 0, colour = COL_TEXT, margin = margin(b = 3)),
    plot.margin = margin(3, 4, 3, 4)
  )

## ------------------------------------------------------------------
## F2: signed effect stem-and-endpoint waterfall.
## ------------------------------------------------------------------

waterfall$plot_class <- factor(
  ifelse(waterfall$expected_direction_preserved, "Expected direction preserved (132)", "Not preserved / near-zero (9)"),
  levels = c("Expected direction preserved (132)", "Not preserved / near-zero (9)")
)
effect_cols <- c("Expected direction preserved (132)" = COL_WARM, "Not preserved / near-zero (9)" = COL_COOL)
label_pos <- waterfall[!is.na(waterfall$label) & waterfall$GSE180286_High_vs_Other_effect > 0, ]
label_neg <- waterfall[!is.na(waterfall$label) & waterfall$GSE180286_High_vs_Other_effect <= 0, ]
y_min <- min(waterfall$GSE180286_High_vs_Other_effect)
y_max <- max(waterfall$GSE180286_High_vs_Other_effect)
median_effect <- median(waterfall$GSE180286_High_vs_Other_effect)

p_water <- ggplot(waterfall, aes(x = effect_rank, y = GSE180286_High_vs_Other_effect, colour = plot_class)) +
  geom_hline(yintercept = 0, colour = COL_TEXT, linewidth = 0.42) +
  geom_segment(aes(xend = effect_rank, y = 0, yend = GSE180286_High_vs_Other_effect), linewidth = 0.30, alpha = 0.78) +
  geom_point(size = 0.72, alpha = 0.96) +
  geom_text_repel(
    data = label_pos, aes(label = label), family = FONT, size = 5.9 / 2.845,
    nudge_x = 7, nudge_y = 0.04, direction = "both", box.padding = 0.25, point.padding = 0.1,
    min.segment.length = 0, segment.size = 0.22, seed = 20260810, max.overlaps = Inf,
    show.legend = FALSE
  ) +
  annotate(
    "text", x = 140, y = y_max * 1.18,
    label = sprintf("132/141 positive\nMedian effect = %.3f", median_effect),
    hjust = 1, vjust = 1, family = FONT, fontface = "bold",
    size = 5.7 / 2.845, lineheight = 0.92, colour = COL_TEXT
  ) +
  geom_text_repel(
    data = label_neg, aes(label = label), family = FONT, size = 5.9 / 2.845,
    nudge_x = -7, nudge_y = -0.022, direction = "both", box.padding = 0.25, point.padding = 0.1,
    min.segment.length = 0, segment.size = 0.22, seed = 20260810, max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = effect_cols, name = NULL) +
  scale_x_continuous(breaks = c(1, 35, 70, 105, 141), limits = c(0, 142), expand = c(0, 0)) +
  scale_y_continuous(limits = c(y_min - 0.055, y_max * 1.24), expand = c(0, 0)) +
  labs(
    title = "Transferred-gene effects",
    x = "Transferred genes ordered by GSE180286 effect",
    y = "GSE180286 High-minus-Other effect"
  ) +
  theme_pub +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.margin = margin(0, 0, 1, 0),
    legend.key.width = grid::unit(3.5, "mm")
  )

## ------------------------------------------------------------------
## F3: predefined-signature weighted running-ES fingerprint.
## ------------------------------------------------------------------

nes <- unique(enrich$enrichment_NES)
fdr <- unique(enrich$enrichment_BH_FDR)
es <- unique(enrich$enrichment_ES)
hits_df <- enrich[enrich$predefined_signature_member, ]
curve_min <- min(enrich$running_ES)
curve_max <- max(enrich$running_ES)
barcode_y0 <- min(-0.085, curve_min - 0.025)
barcode_y1 <- barcode_y0 + 0.060

p_gsea <- ggplot(enrich, aes(x = rank, y = running_ES)) +
  geom_hline(yintercept = 0, colour = COL_MUTED, linewidth = 0.32) +
  geom_area(fill = COL_WARM, alpha = 0.12) +
  geom_line(colour = COL_WARM, linewidth = 0.78) +
  geom_segment(
    data = hits_df,
    aes(x = rank, xend = rank, y = barcode_y0, yend = barcode_y1),
    inherit.aes = FALSE, colour = COL_TEXT, linewidth = 0.30, alpha = 0.96
  ) +
  annotate(
    "text", x = 0.97 * max(enrich$rank), y = curve_max * 1.08,
    label = sprintf("ES = %.2f   NES = %.2f\nBH FDR = %.2g", es, nes, fdr),
    hjust = 1, vjust = 1, family = FONT, fontface = "bold", size = 6.2 / 2.845, colour = COL_TEXT
  ) +
  annotate(
    "text", x = 0.03 * max(enrich$rank), y = barcode_y0 - 0.012,
    label = "High-associated", hjust = 0, vjust = 1,
    family = FONT, size = 5.9 / 2.845, colour = COL_WARM
  ) +
  annotate(
    "text", x = 0.97 * max(enrich$rank), y = barcode_y0 - 0.012,
    label = "Other-associated", hjust = 1, vjust = 1,
    family = FONT, size = 5.9 / 2.845, colour = COL_COOL
  ) +
  scale_x_continuous(breaks = c(1, 10000, 20000, 30000), labels = comma, expand = c(0, 0)) +
  scale_y_continuous(limits = c(barcode_y0 - 0.05, curve_max * 1.13), expand = c(0, 0)) +
  labs(
    title = "Predefined-signature enrichment",
    x = "GSE180286 transcriptome rank",
    y = "Running enrichment score"
  ) +
  theme_pub

## ------------------------------------------------------------------
## F4: concise provenance and limitation strip.
## ------------------------------------------------------------------

footer_text <- "GSE180286 orthogonal dataset transfer; not patient-level replication."
p_footer <- ggplot() +
  annotate("segment", x = 0, xend = 1, y = 0.82, yend = 0.82, colour = COL_GRID, linewidth = 0.45) +
  annotate("text", x = 0.5, y = 0.36, label = footer_text, hjust = 0.5, vjust = 0.5,
           family = FONT, size = 6.4 / 2.845, colour = COL_MUTED) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void(base_family = FONT) +
  theme(plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(0, 4, 0, 4))

main_row <- p_flow + p_water + p_gsea + plot_layout(widths = c(1.48, 1.42, 1.38))
combined <- main_row / p_footer +
  plot_layout(heights = c(1, 0.09)) +
  plot_annotation(
    title = "F  Orthogonal transfer of the YAP-Stem transcriptional program",
    theme = theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(family = FONT, size = 9.4, face = "bold", hjust = 0, colour = COL_TEXT, margin = margin(b = 4)),
      plot.margin = margin(4, 5, 3, 5, unit = "mm")
    )
  )

pdf_file <- file.path(out_dir, "Figure4F_GSE180286_program_transfer.pdf")
svg_file <- file.path(out_dir, "Figure4F_GSE180286_program_transfer.svg")
png_file <- file.path(out_dir, "Figure4F_GSE180286_program_transfer.png")

ggsave(pdf_file, combined, width = 180, height = 88, units = "mm", device = grDevices::cairo_pdf, bg = "white")
ggsave(svg_file, combined, width = 180, height = 88, units = "mm", device = svglite::svglite, bg = "white")
ggsave(png_file, combined, width = 180, height = 88, units = "mm", dpi = 600, bg = "white")
stopifnot(
  file.exists(pdf_file), file.info(pdf_file)$size > 10000,
  file.exists(svg_file), file.info(svg_file)$size > 10000,
  file.exists(png_file), file.info(png_file)$size > 10000
)

input_hash_after <- vapply(paths, sha256, character(1))
stopifnot(identical(input_hash_before, input_hash_after))

source_hashes <- data.frame(
  file_role = c(names(paths), "figure_pdf", "figure_svg", "figure_png"),
  file = c(unname(unlist(paths)), pdf_file, svg_file, png_file),
  sha256 = c(unname(input_hash_before), sha256(pdf_file), sha256(svg_file), sha256(png_file)),
  stringsAsFactors = FALSE
)
write_tsv(source_hashes, file.path(out_dir, "Figure4F_source_hashes.tsv"))

qc_lines <- c(
  "Figure 4F QC - orthogonal transfer of the GSE176078 transcriptional program",
  "==========================================================================",
  "DATA_GATE: PASS",
  "Discovery program source: GSE176078 score-independent frozen 146-gene program",
  "Discovery program modified: NO",
  "Validation-dataset gene reselection: NO",
  "Source genes: 146",
  "Detectable in GSE180286: 141",
  "Unavailable in GSE180286: 5",
  paste0("Unavailable genes: ", audit_value("unavailable_genes")),
  "Expected direction preserved: 132/141 (93.6%)",
  "Not preserved / near-zero: 9/141 (8 negative; 1 exactly zero)",
  paste0("Median transferred effect: ", format(median_effect, digits = 12)),
  paste0("GSE180286 ranked genes: ", audit_value("GSE180286_ranking_genes")),
  paste0("Enrichment ES: ", audit_value("enrichment_ES")),
  paste0("Enrichment NES: ", audit_value("enrichment_NES")),
  paste0("Enrichment BH FDR: ", audit_value("enrichment_BH_FDR")),
  "Predefined complete measurable signature used for enrichment: YES",
  "Exact transfer source uses pooled P5 Primary/LN High-minus-Other mean effects: YES",
  "Formal sample-blocked model in the exact 132/141 and NES=1.628 source: NO",
  "GSM5457211_E2020-1 mapping: P5 Primary, 5297 TNBC malignant cells",
  "GSM5457212_E2020-2 mapping: P5 LN+, 162 TNBC malignant cells",
  "Both included samples map to the same patient P5: YES",
  "Earlier two-patient note reconciled: P4 contributed one cell to the old object but was excluded from transfer effect",
  "Full sample counts moved to legend/provenance audit: YES",
  "Concise limitation note stated in figure: YES",
  "P4 single-cell sample excluded: YES",
  "Independent patient-level replication claim: NO",
  "Figure 4H pathway metrics used: NO",
  "Large-scale external patient-validation data used: NO",
  "Dataset label Wu in figure or required TSV outputs: NO",
  "Yan label or data in figure or required TSV outputs: NO",
  "Panel width: 180 mm",
  "Panel height: 88 mm",
  "PDF: vector",
  "SVG: vector",
  "PNG: 600 dpi, white background",
  "PLOT_STATUS: GENERATED_PENDING_VISUAL_QC"
)
writeLines(qc_lines, file.path(out_dir, "Figure4F_QC.txt"), useBytes = TRUE)

cat("Figure 4F rendered: PDF, SVG, PNG\n")
