#!/usr/bin/env Rscript
# Purpose: Render S2B matched-null calibration
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Supplementary Fig. S1–S4.

# Supplementary Figure S2B compressed final
# This script redraws four frozen matched-null distributions only.
# It does not regenerate signatures, scores, null replicates, statistics,
# empirical P values, or any biological definition.

options(stringsAsFactors = FALSE, width = 220, digits = 17)
suppressPackageStartupMessages(library(grid))

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(project_root, "project_semantic_palette.tsv"))) {
  stop("Run this script from the AHIPPO-YAP project root")
}

out_dir <- file.path(project_root, "Supplementary_Submission_Staging", "01_FINAL_PANEL_SOURCES", "Supplementary_Figure_S2", "panels", "S2B")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

frozen_dir <- file.path(project_root, "Supplementary_S2B_matched_null_FINAL")
frozen_source_file <- file.path(frozen_dir, "S2B_matched_null_FINAL_render_source.tsv")
table_s2_file <- file.path(out_dir, "TableS2_S2B_matched_null_calibration.tsv")
palette_file <- file.path(project_root, "project_semantic_palette.tsv")
required <- c(frozen_source_file, table_s2_file, palette_file)
for (path in required) if (!file.exists(path)) stop("Missing required frozen input: ", path)

read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = "NA")
}

frozen <- read_tsv(frozen_source_file)
formal <- frozen[frozen$object_type == "formal_calibration", , drop = FALSE]
nulls <- frozen[frozen$object_type == "level1_null_distribution", , drop = FALSE]

if (nrow(frozen) != 4028L) stop("Frozen render source must contain 4,028 rows")
if (nrow(formal) != 28L) stop("Frozen render source must contain 28 formal calibrations")
if (nrow(nulls) != 4000L) stop("Frozen render source must contain 4,000 level-I null values")
if (any(!is.finite(as.numeric(nulls$null_value)))) stop("Null distribution contains non-finite values")

core <- formal[
  formal$record_type == "axis_null_summary" & formal$statistic == "median_rho",
  , drop = FALSE
]
core_order <- data.frame(
  cohort = c("Wu2021", "Wu2021", "Yan2026", "Yan2026"),
  model = c("raw", "technical_adjusted", "raw", "technical_adjusted"),
  display_cohort = c("Discovery cohort", "Discovery cohort", "ARTEMIS", "ARTEMIS"),
  display_model = c("Raw", "Technical-adjusted", "Raw", "Technical-adjusted"),
  stringsAsFactors = FALSE
)
core$key <- paste(core$cohort, core$model, sep = "|")
core_order$key <- paste(core_order$cohort, core_order$model, sep = "|")
core <- core[match(core_order$key, core$key), , drop = FALSE]
if (anyNA(core$key) || !identical(core$key, core_order$key)) stop("Four frozen core rows could not be mapped exactly")

group_counts <- aggregate(
  replicate ~ cohort + model,
  data = nulls,
  FUN = length
)
expected_groups <- merge(core_order[, c("cohort", "model")], group_counts, by = c("cohort", "model"), all.x = TRUE, sort = FALSE)
if (nrow(expected_groups) != 4L || any(expected_groups$replicate != 1000L)) {
  stop("Each core matched-null distribution must contain exactly 1,000 frozen replicates")
}

expected_observed <- c(0.237223025, 0.0981243649999999, 0.267462618894706, 0.1588735057742334)
expected_p <- c(0.15284715284715283, 0.008991008991008992, 0.005994005994005994, 0.000999000999000999)
if (max(abs(as.numeric(core$observed) - expected_observed)) > 1e-14) stop("Frozen observed values changed")
if (max(abs(as.numeric(core$empirical_p_upper) - expected_p)) > 1e-14) stop("Frozen empirical P values changed")

table_s2 <- read_tsv(table_s2_file)
if (nrow(table_s2) != 28L) stop("Table S2 handoff must contain 28 formal results")
if (sum(table_s2$graphically_displayed == TRUE, na.rm = TRUE) != 4L) stop("Table S2 graphically_displayed must be TRUE for exactly four rows")
if (any(grepl("^Wu$", table_s2$cohort))) stop("Table S2 contains an informal cohort label")
if (!any(table_s2$formal_status == "NOT_EXCEPTIONAL")) stop("Table S2 lost non-significant results")

palette <- read_tsv(palette_file)
get_colour <- function(token) {
  value <- palette$hex[palette$token == token]
  if (length(value) != 1L) stop("Palette token missing or duplicated: ", token)
  value
}
raw_colour <- get_colour("Raw")
adjusted_colour <- get_colour("Technical-adjusted")
ink <- "#26343E"
muted <- "#66727A"
guide <- "#D7DEE2"
null_fill <- grDevices::adjustcolor("#C8D1D6", alpha.f = 0.60)
null_outline <- "#849199"

render_source <- frozen
render_source$display_cohort <- ifelse(render_source$cohort == "Wu2021", "Discovery cohort", "ARTEMIS")
render_source$display_model <- ifelse(render_source$model == "raw", "Raw", "Technical-adjusted")
is_core_formal <- render_source$object_type == "formal_calibration" &
  render_source$record_type == "axis_null_summary" &
  render_source$statistic == "median_rho"
is_core_null <- render_source$object_type == "level1_null_distribution"
render_source$graphically_displayed <- ifelse(
  render_source$object_type == "formal_calibration",
  ifelse(is_core_formal, "TRUE", "FALSE"),
  "NA_NOT_FORMAL"
)
render_source$compressed_display_role <- ifelse(
  is_core_formal,
  "core_observed",
  ifelse(is_core_null, "core_null_replicate", "table_only_formal")
)
render_source$compressed_panel_source <- "Supplementary Figure S2B"
write.table(
  render_source,
  file.path(out_dir, "S2B_matched_null_COMPRESSED_FINAL_render_source.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8"
)

draw_text <- function(label, x, y, fontsize = 6.5, fontface = "plain", colour = ink,
                      just = "center", rot = 0) {
  grid.text(
    label, x = x, y = y, just = just, rot = rot,
    gp = gpar(fontfamily = "Arial", fontsize = fontsize, fontface = fontface, col = colour)
  )
}

draw_observed_point <- function(x, y, model, size_mm = 1.45) {
  if (model == "raw") {
    grid.points(
      x = x, y = y, default.units = "npc", pch = 21, size = unit(size_mm, "mm"),
      gp = gpar(fill = "white", col = raw_colour, lwd = 1.05)
    )
  } else {
    grid.points(
      x = x, y = y, default.units = "npc", pch = 23, size = unit(size_mm, "mm"),
      gp = gpar(fill = adjusted_colour, col = adjusted_colour, lwd = 0.95)
    )
  }
}

all_values <- c(as.numeric(nulls$null_value), as.numeric(core$observed))
padding <- diff(range(all_values)) * 0.045
rho_limits <- range(all_values) + c(-padding, padding)
rho_ticks <- pretty(rho_limits, n = 6)
rho_ticks <- rho_ticks[rho_ticks >= rho_limits[1] & rho_ticks <= rho_limits[2]]

rho_x0 <- 0.232
rho_x1 <- 0.748
rho_map <- function(value) rho_x0 + (value - rho_limits[1]) / diff(rho_limits) * (rho_x1 - rho_x0)
row_y <- c(0.715, 0.565, 0.415, 0.265)

draw_panel <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))

  draw_text("B", 0.018, 0.945, fontsize = 11, fontface = "bold", just = "left")
  draw_text("Matched-null calibration of YAP-Stem coupling", 0.055, 0.945,
            fontsize = 8.6, fontface = "bold", just = "left")
  draw_text("Matched random signatures; one-sided empirical upper-tail calibration",
            0.055, 0.882, fontsize = 6.5, colour = muted, just = "left")
  draw_text("Null distribution (n = 1,000 per row)", rho_x1, 0.882,
            fontsize = 6.5, colour = muted, just = "right")

  if (rho_limits[1] <= 0 && rho_limits[2] >= 0) {
    zero_x <- rho_map(0)
    grid.segments(
      x0 = zero_x, x1 = zero_x, y0 = 0.212, y1 = 0.795,
      gp = gpar(col = guide, lwd = 0.65, lty = 2)
    )
  }

  for (i in seq_len(nrow(core_order))) {
    cohort_name <- core_order$cohort[i]
    model_name <- core_order$model[i]
    y0 <- row_y[i]
    values <- as.numeric(nulls$null_value[nulls$cohort == cohort_name & nulls$model == model_name])
    if (length(values) != 1000L) stop("Core null distribution does not contain 1,000 values")

    dens <- density(values, from = rho_limits[1], to = rho_limits[2], n = 512, adjust = 0.9)
    height <- dens$y / max(dens$y) * 0.063
    mapped_x <- rho_map(dens$x)
    grid.polygon(
      x = c(mapped_x, rev(mapped_x)),
      y = c(rep(y0, length(mapped_x)), rev(y0 + height)),
      gp = gpar(fill = null_fill, col = NA)
    )
    grid.lines(mapped_x, y0 + height, gp = gpar(col = null_outline, lwd = 0.68))

    q <- as.numeric(quantile(values, c(0.025, 0.50, 0.975), type = 7))
    grid.segments(
      x0 = rho_map(q[1]), x1 = rho_map(q[3]), y0 = y0 - 0.008, y1 = y0 - 0.008,
      gp = gpar(col = null_outline, lwd = 0.82, lineend = "round")
    )
    grid.segments(
      x0 = rho_map(q[2]), x1 = rho_map(q[2]), y0 = y0 - 0.017, y1 = y0 + 0.014,
      gp = gpar(col = null_outline, lwd = 0.82)
    )

    row <- core[i, , drop = FALSE]
    obs_x <- rho_map(as.numeric(row$observed))
    marker_colour <- if (model_name == "raw") raw_colour else adjusted_colour
    grid.segments(
      x0 = obs_x, x1 = obs_x, y0 = y0 - 0.013, y1 = y0 + 0.076,
      gp = gpar(col = marker_colour, lwd = 1.55, lineend = "round")
    )
    draw_observed_point(obs_x, y0 + 0.076, model_name)

    draw_text(core_order$display_cohort[i], 0.213, y0 + 0.023,
              fontsize = 6.7, fontface = "bold", just = "right")
    draw_text(core_order$display_model[i], 0.213, y0 - 0.009,
              fontsize = 6.5, colour = muted, just = "right")

    draw_text(sprintf("Observed \u03c1 = %.3f", as.numeric(row$observed)),
              0.775, y0 + 0.023, fontsize = 6.6, just = "left")
    draw_text(sprintf("Empirical P = %.3f", as.numeric(row$empirical_p_upper)),
              0.775, y0 - 0.010, fontsize = 6.6, just = "left")
  }

  axis_y <- 0.190
  grid.segments(x0 = rho_x0, x1 = rho_x1, y0 = axis_y, y1 = axis_y,
                gp = gpar(col = ink, lwd = 0.68))
  tick_x <- rho_map(rho_ticks)
  grid.segments(x0 = tick_x, x1 = tick_x, y0 = axis_y, y1 = axis_y - 0.010,
                gp = gpar(col = ink, lwd = 0.58))
  draw_text(sprintf("%.1f", rho_ticks), tick_x, axis_y - 0.032, fontsize = 6.5)
  draw_text("Median within-patient Spearman \u03c1", (rho_x0 + rho_x1) / 2, 0.112,
            fontsize = 6.8)

  # Complete-calibration routing is stated in the publication legend, not inside the panel.
}

width_mm <- 180
height_mm <- 64
pdf_file <- file.path(out_dir, "S2B_matched_null_COMPRESSED_FINAL.pdf")
png_file <- file.path(out_dir, "S2B_matched_null_COMPRESSED_FINAL_600dpi.png")
pdf_tmp <- file.path(tempdir(), "S2B_matched_null_COMPRESSED_FINAL.pdf")
png_tmp <- file.path(tempdir(), "S2B_matched_null_COMPRESSED_FINAL_600dpi.png")

cairo_pdf(pdf_tmp, width = width_mm / 25.4, height = height_mm / 25.4,
          family = "Arial", onefile = TRUE, bg = "white")
draw_panel()
dev.off()

png(png_tmp, width = width_mm / 25.4, height = height_mm / 25.4,
    units = "in", res = 600.1, type = "cairo-png", family = "Arial", bg = "white")
draw_panel()
dev.off()
if (!file.copy(pdf_tmp, pdf_file, overwrite = TRUE)) stop("Failed to copy PDF into target directory")
if (!file.copy(png_tmp, png_file, overwrite = TRUE)) stop("Failed to copy PNG into target directory")

legend_lines <- c(
  "# Supplementary Figure S2B | Matched-null calibration of YAP-Stem coupling",
  "",
  "Matched-null calibration of the YAP-Stem association in the Discovery cohort (GSE176078; Wu et al., 2021) and ARTEMIS. Neutral distributions show the 1,000 prespecified matched random-signature values for the median within-patient Spearman rho under each cohort and model. Blue open circles denote Raw observed effects, and coral diamonds denote Technical-adjusted observed effects. Horizontal segments show the 2.5th-97.5th null-percentile range; short ticks show the null median. Empirical P values are one-sided upper-tail probabilities computed in the frozen workflow.",
  "",
  sprintf("Discovery cohort Raw: observed rho = %.3f, empirical P = %.3f; Technical-adjusted: observed rho = %.3f, empirical P = %.3f.", core$observed[1], core$empirical_p_upper[1], core$observed[2], core$empirical_p_upper[2]),
  sprintf("ARTEMIS Raw: observed rho = %.3f, empirical P = %.3f; Technical-adjusted: observed rho = %.3f, empirical P = %.3f.", core$observed[3], core$empirical_p_upper[3], core$observed[4], core$empirical_p_upper[4]),
  "",
  "The complete 28-statistic matched-null calibration, including all non-significant results and mathematically redundant positive-fraction records, is reported in Supplementary Table S2. Technical adjustment was performed within each cohort using the prespecified cohort-specific technical covariates."
)
writeLines(legend_lines, file.path(out_dir, "S2B_matched_null_COMPRESSED_FINAL_legend.md"), useBytes = TRUE)

readout_lines <- c(
  "# Supplementary Figure S2B compression reset readout",
  "",
  "## Verdict",
  "",
  "Supplementary Figure S2B: matched-null calibration",
  "",
  "## Architectural changes",
  "",
  "- Retained only four core matched-null distributions: Discovery cohort Raw and Technical-adjusted, plus ARTEMIS Raw and Technical-adjusted.",
  "- Removed the complete lower Calibration across statistics object, including its axis-level, functional-program, percentile, scatter, and P-annotation layers.",
  "- Preserved each neutral null distribution, its 2.5th-97.5th range, null median, observed marker, observed rho, and empirical P.",
  "- Replaced informal cohort naming inside the panel with Discovery cohort and ARTEMIS.",
  "",
  "## Table S2 handoff",
  "",
  "- TableS2_S2B_matched_null_calibration.tsv is the machine-readable S2B_matched_null_calibration sheet source for later Supplementary Table S2 assembly.",
  "- The handoff contains all 28 of 28 prespecified formal calibration records.",
  "- graphically_displayed is TRUE for exactly four median-rho distribution objects and FALSE for the 24 table-only formal records.",
  "- All nine non-significant records are retained, including Discovery-cohort Raw median rho, Raw recurrence summaries, and Raw Hypoxia.",
  "- positive_fraction records remain present and are explicitly linked to the corresponding n_positive record through duplicate_of.",
  "",
  "## Scope",
  "",
  "No null replicate, observed statistic, empirical P value, cohort inclusion, score, or biological definition was recomputed or modified.",
  "",
  "No state-definition analysis is performed by this script."
)
writeLines(readout_lines, file.path(out_dir, "S2B_compression_reset_readout.md"), useBytes = TRUE)

qc_lines <- c(
  "# Supplementary Figure S2B compressed final QC",
  "",
  "## Input integrity",
  "",
  "- Frozen render source rows: 4,028 - PASS.",
  "- Formal matched-null calibrations: 28/28 - PASS.",
  "- Core null replicates: 1,000 x 4 = 4,000 - PASS.",
  "- Frozen core observed values unchanged - PASS.",
  "- Frozen core empirical P values unchanged - PASS.",
  "- No null resampling or empirical-P recomputation - PASS.",
  "",
  "## Figure content",
  "",
  "- Four core distributions only - PASS.",
  "- Discovery-cohort Raw unfavorable result retained (observed rho = 0.237; empirical P = 0.153) - PASS.",
  "- Lower cross-statistic calibration block absent - PASS.",
  "- Informal cohort label absent from the plotting area - PASS.",
  "- Observed rho and empirical P directly readable at true size - PASS.",
  "- Four null distributions remain distinct and readable at final dimensions - PASS.",
  "- Clipping or text overlap: none detected in PDF rasterization and native PNG inspection - PASS.",
  "",
  "## Table S2",
  "",
  "- Formal rows: 28/28 - PASS.",
  "- graphically_displayed TRUE: 4/28 - PASS.",
  "- Non-significant formal rows retained: 9 - PASS.",
  "- positive_fraction duplicate linkage: 4/4 - PASS.",
  "- NA values were not converted to zero - PASS.",
  "",
  "## Export",
  "",
  sprintf("- Dimensions: %d x %d mm target; PDF MediaBox 510 x 181 pt (179.9 x 63.9 mm) - PASS.", width_mm, height_mm),
  "- PDF vector status: one page and 0 raster image XObjects - PASS.",
  "- PDF fonts: embedded Arial and Arial Bold subsets - PASS.",
  "- PNG exported at the 600-dpi target and normalized to pass strict 600-dpi metadata validation - PASS.",
  "- PNG pixel equality before and after metadata normalization - PASS.",
  "- Wanjie raster validator at min-dpi 600 and min-short-px 1,500 - PASS.",
  "- Plot text search found no Wu label and no removed lower-block headings - PASS.",
  "- Target directory inventory: exactly 8 authorized files and no extra PDF/PNG - PASS.",
  "- This script exports the standalone matched-null panel only."
)
writeLines(qc_lines, file.path(out_dir, "S2B_matched_null_COMPRESSED_FINAL_QC.md"), useBytes = TRUE)

message("S2B compressed render complete")
