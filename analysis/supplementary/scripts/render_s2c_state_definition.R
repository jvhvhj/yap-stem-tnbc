#!/usr/bin/env Rscript
# Purpose: Render S2C state-definition robustness
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Supplementary Fig. S1–S4.

# Supplementary Figure S2C candidate
# Mechanical visualization of frozen patient-level Program146 effects.
# No state, score, effect, or statistical test is recomputed here.

options(stringsAsFactors = FALSE)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_dir_override <- Sys.getenv("S2C_OUTPUT_DIR", unset = "")
out_dir <- if (nzchar(out_dir_override)) {
  normalizePath(out_dir_override, winslash = "/", mustWork = FALSE)
} else {
  file.path(project_root, "Supplementary_S2C_state_definition_candidate")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source_dir <- file.path(project_root, "Supplementary_S2C_state_definition_candidate")
render_source_file <- file.path(source_dir, "S2C_state_definition_render_source.tsv")
palette_file <- file.path(project_root, "project_semantic_palette.tsv")

if (!file.exists(render_source_file)) stop("Missing frozen render source: ", render_source_file)
if (!file.exists(palette_file)) stop("Missing project semantic palette: ", palette_file)

d <- read.delim(render_source_file, check.names = FALSE, stringsAsFactors = FALSE,
                na.strings = c("", "NA"), fileEncoding = "UTF-8")
palette <- read.delim(palette_file, check.names = FALSE, stringsAsFactors = FALSE)

get_colour <- function(token) {
  z <- palette$hex[palette$token == token]
  if (length(z) != 1L) stop("Palette token missing or duplicated: ", token)
  z
}

raw_colour <- get_colour("Raw")
adjusted_colour <- get_colour("Technical-adjusted")
ink <- "#26343E"
muted <- "#66727A"
guide <- "#D9E0E4"
zero_colour <- "#98A3AA"

required_columns <- c(
  "source_record_type", "source_row", "patient_id", "state_definition_source",
  "state_definition_display", "definition_qualifier", "definition_order",
  "model_source", "model_display", "model_order", "testable_status", "high_n",
  "other_n", "Program146_effect", "measurable_genes", "requested_genes",
  "graphically_displayed_status", "graphical_reason", "n_cells", "lane_q1",
  "lane_median", "lane_q3", "positive_n", "testable_n",
  "definition_centre_npc", "lane_centre_npc", "jitter_offset_npc", "point_y_npc",
  "x_display_min", "x_display_max", "source_file", "source_sha256"
)
if (!identical(names(d), required_columns)) stop("Render-source schema changed")
if (nrow(d) != 474L) stop("Render source must retain all 474 authoritative records")
if (sum(d$source_record_type == "patient") != 468L) stop("Expected 468 patient records")
if (sum(d$source_record_type == "summary") != 6L) stop("Expected six formal summary records")

patient_all <- d[d$source_record_type == "patient", , drop = FALSE]
patient_plot <- patient_all[patient_all$graphically_displayed_status == "PLOTTED_POINT", , drop = FALSE]
summary_data <- d[d$source_record_type == "summary", , drop = FALSE]

if (nrow(patient_plot) != 466L) stop("Expected 466 testable plotted patient effects")
if (any(!is.finite(patient_plot$Program146_effect))) stop("A plotted effect is non-finite")
if (any(patient_plot$Program146_effect <= 0)) stop("A plotted effect is zero or negative")
if (length(unique(patient_plot$source_sha256)) != 1L) stop("Source hash inconsistency")
if (!all(patient_plot$measurable_genes == 141L & patient_plot$requested_genes == 146L)) {
  stop("Program146 coverage changed")
}

p59 <- patient_all[
  patient_all$patient_id == "P59" &
    patient_all$state_definition_source == "Primary pooled tertile",
  , drop = FALSE
]
if (nrow(p59) != 2L || any(p59$testable_status != "NON_TESTABLE") ||
    any(p59$high_n != 50) || any(p59$other_n != 0) ||
    any(!is.na(p59$Program146_effect)) ||
    any(p59$graphically_displayed_status != "NON_TESTABLE_NOT_PLOTTED")) {
  stop("P59 pooled non-testable rule failed")
}

summary_data <- summary_data[order(summary_data$definition_order, summary_data$model_order), ]
expected <- data.frame(
  definition_order = rep(1:3, each = 2),
  model_order = rep(1:2, 3),
  testable_n = c(77, 77, 78, 78, 78, 78),
  positive_n = c(77, 77, 78, 78, 78, 78),
  median = c(
    0.117553979158401, 0.0791866546683832,
    0.094327911734581, 0.0635475032759759,
    0.0967041105031967, 0.0661845241857039
  )
)
if (!identical(as.integer(summary_data$definition_order), expected$definition_order) ||
    !identical(as.integer(summary_data$model_order), expected$model_order) ||
    any(as.integer(summary_data$testable_n) != expected$testable_n) ||
    any(as.integer(summary_data$positive_n) != expected$positive_n) ||
    any(abs(summary_data$lane_median - expected$median) > 1e-12)) {
  stop("Formal summary rows no longer reproduce the approved six results")
}

x_limits <- unique(c(d$x_display_min, d$x_display_max))
x_limits <- x_limits[is.finite(x_limits)]
if (length(x_limits) != 2L) stop("Display limits are not uniquely defined")
x_limits <- range(x_limits)
if (min(patient_plot$Program146_effect) <= x_limits[1] ||
    max(patient_plot$Program146_effect) >= x_limits[2]) stop("Display limits would clip points")
if (!(x_limits[1] < 0 && x_limits[2] > 0)) stop("Shared axis must include effect=0")

x_ticks <- seq(0, floor(x_limits[2] / 0.1) * 0.1, by = 0.1)
if (tail(x_ticks, 1) < x_limits[2] - 0.04) x_ticks <- c(x_ticks, tail(x_ticks, 1) + 0.1)
x_ticks <- x_ticks[x_ticks <= x_limits[2] + 1e-12]

plot_x0 <- 0.300
plot_x1 <- 0.775
count_x <- 0.900
axis_y <- 0.095

x_map <- function(x) {
  plot_x0 + (x - x_limits[1]) / diff(x_limits) * (plot_x1 - plot_x0)
}

draw_text <- function(label, x, y, fontsize = 6.5, fontface = "plain",
                      colour = ink, just = "center", rot = 0) {
  grid::grid.text(
    label, x = x, y = y, just = just, rot = rot,
    gp = grid::gpar(fontfamily = "Arial", fontsize = fontsize,
                    fontface = fontface, col = colour)
  )
}

draw_shape_legend <- function() {
  grid::grid.points(
    x = 0.385, y = 0.830, default.units = "npc",
    pch = 21, size = grid::unit(1.18, "mm"),
    gp = grid::gpar(fill = "white", col = raw_colour, lwd = 0.70)
  )
  draw_text("Raw", 0.405, 0.830, fontsize = 6.4, just = "left")
  grid::grid.points(
    x = 0.525, y = 0.830, default.units = "npc",
    pch = 23, size = grid::unit(1.20, "mm"),
    gp = grid::gpar(fill = adjusted_colour, col = ink, lwd = 0.58)
  )
  draw_text("Technical-adjusted", 0.548, 0.830, fontsize = 6.4, just = "left")
}

draw_lane <- function(definition_order, model_order) {
  lane <- patient_plot[
    patient_plot$definition_order == definition_order &
      patient_plot$model_order == model_order,
    , drop = FALSE
  ]
  sm <- summary_data[
    summary_data$definition_order == definition_order &
      summary_data$model_order == model_order,
    , drop = FALSE
  ]
  if (nrow(sm) != 1L) stop("Lane summary mapping failed")
  if (nrow(lane) != sm$testable_n) stop("Lane point count differs from formal testable n")

  y_lane <- unique(lane$lane_centre_npc)
  if (length(y_lane) != 1L) stop("Lane centre inconsistency")
  series_colour <- if (model_order == 1L) raw_colour else adjusted_colour
  iqr_colour <- grDevices::adjustcolor(series_colour, alpha.f = 0.62)

  grid::grid.segments(
    x0 = x_map(sm$lane_q1), x1 = x_map(sm$lane_q3), y0 = y_lane, y1 = y_lane,
    gp = grid::gpar(col = iqr_colour, lwd = 1.80, lineend = "round")
  )

  if (model_order == 1L) {
    grid::grid.points(
      x = x_map(lane$Program146_effect), y = lane$point_y_npc,
      default.units = "npc", pch = 21, size = grid::unit(0.95, "mm"),
      gp = grid::gpar(
        fill = grDevices::adjustcolor("white", alpha.f = 0.94),
        col = grDevices::adjustcolor(raw_colour, alpha.f = 0.78), lwd = 0.55
      )
    )
  } else {
    grid::grid.points(
      x = x_map(lane$Program146_effect), y = lane$point_y_npc,
      default.units = "npc", pch = 23, size = grid::unit(0.95, "mm"),
      gp = grid::gpar(
        fill = grDevices::adjustcolor(adjusted_colour, alpha.f = 0.72),
        col = grDevices::adjustcolor(ink, alpha.f = 0.58), lwd = 0.48
      )
    )
  }

  med_x <- x_map(sm$lane_median)
  grid::grid.segments(
    x0 = med_x, x1 = med_x, y0 = y_lane - 0.014, y1 = y_lane + 0.014,
    gp = grid::gpar(col = ink, lwd = 1.08, lineend = "round")
  )
}

draw_panel <- function() {
  grid::grid.newpage()

  draw_text("C", 0.022, 0.965, fontsize = 11.0, fontface = "bold", just = c("left", "top"))
  draw_text("Program146 effects across alternative", 0.087, 0.958,
            fontsize = 7.6, fontface = "bold", just = c("left", "top"))
  draw_text("YAP-Stem state definitions", 0.087, 0.908,
            fontsize = 7.6, fontface = "bold", just = c("left", "top"))
  draw_shape_legend()

  draw_text("Positive n/N", count_x, 0.835, fontsize = 6.4,
            fontface = "bold", colour = muted)
  draw_text("Raw → adjusted", count_x, 0.805, fontsize = 6.2,
            fontface = "bold", colour = muted)

  x_zero <- x_map(0)
  grid::grid.segments(
    x0 = x_zero, x1 = x_zero, y0 = 0.145, y1 = 0.750,
    gp = grid::gpar(col = zero_colour, lwd = 0.65, lty = 2)
  )

  for (sep_y in c(0.540, 0.320)) {
    grid::grid.segments(
      x0 = 0.020, x1 = 0.975, y0 = sep_y, y1 = sep_y,
      gp = grid::gpar(col = guide, lwd = 0.42)
    )
  }

  label_lines <- list(
    c("Pooled tertile", "(primary)"),
    c("Within-patient", "tertile"),
    c("Patient-z", "tertile")
  )
  definition_centres <- c(0.650, 0.430, 0.210)
  count_labels <- c("77/77 → 77/77", "78/78 → 78/78", "78/78 → 78/78")

  for (i in 1:3) {
    centre <- definition_centres[i]
    draw_text(label_lines[[i]][1], 0.020, centre + 0.058,
              fontsize = 6.7, fontface = "bold", just = "left")
    draw_text(label_lines[[i]][2], 0.020, centre + 0.026,
              fontsize = 6.3, colour = if (i == 1) muted else ink, just = "left")
    draw_text("Raw", 0.280, centre + 0.028, fontsize = 6.2, colour = muted, just = "right")
    draw_text("Adjusted", 0.280, centre - 0.028, fontsize = 6.2, colour = muted, just = "right")
    draw_text(count_labels[i], count_x, centre, fontsize = 6.4, colour = ink)
    draw_lane(i, 1L)
    draw_lane(i, 2L)
  }

  grid::grid.segments(
    x0 = plot_x0, x1 = plot_x1, y0 = axis_y, y1 = axis_y,
    gp = grid::gpar(col = ink, lwd = 0.62)
  )
  tick_x <- x_map(x_ticks)
  grid::grid.segments(
    x0 = tick_x, x1 = tick_x, y0 = axis_y, y1 = axis_y - 0.010,
    gp = grid::gpar(col = ink, lwd = 0.52)
  )
  draw_text(sprintf("%.1f", x_ticks), tick_x, axis_y - 0.032, fontsize = 6.2)
  draw_text("Program146 High − Other effect", (plot_x0 + plot_x1) / 2, 0.018,
            fontsize = 6.7)
}

output_stem <- if (nzchar(out_dir_override)) "S2C_state_definition_FINAL" else "S2C_state_definition_candidate"
pdf_file <- file.path(out_dir, paste0(output_stem, ".pdf"))
png_file <- file.path(out_dir, paste0(output_stem, "_600dpi.png"))
pdf_tmp <- tempfile(pattern = "S2C_candidate_", fileext = ".pdf")
png_tmp <- tempfile(pattern = "S2C_candidate_", fileext = ".png")

grDevices::cairo_pdf(
  filename = pdf_tmp,
  width = 88 / 25.4,
  height = 78 / 25.4,
  family = "Arial",
  bg = "white",
  onefile = TRUE
)
draw_panel()
grDevices::dev.off()

grDevices::png(
  filename = png_tmp,
  width = 88,
  height = 78,
  units = "mm",
  res = 601,
  type = "windows",
  bg = "white"
)
draw_panel()
grDevices::dev.off()

if (!file.copy(pdf_tmp, pdf_file, overwrite = TRUE)) stop("Failed to copy PDF to output directory")
if (!file.copy(png_tmp, png_file, overwrite = TRUE)) stop("Failed to copy PNG to output directory")
unlink(c(pdf_tmp, png_tmp))

write_utf8_copy <- function(lines, destination) {
  tmp <- tempfile(pattern = "S2C_text_", fileext = ".md")
  writeLines(lines, tmp, useBytes = TRUE)
  ok <- file.copy(tmp, destination, overwrite = TRUE)
  unlink(tmp)
  if (!ok) stop("Failed to copy text output: ", destination)
}

legend_lines <- c(
  "# Supplementary Figure S2C | Program146 effects across alternative YAP-Stem state definitions",
  "",
  "Patient-level Program146 High − Other effects in ARTEMIS under three prespecified YAP-Stem state definitions: the primary pooled upper tertile, a within-patient upper tertile using the frozen cohort-z joint axis, and a patient-z upper tertile after within-patient standardization of the YAP and Stemness component scores. Each point is one testable patient; open blue circles show Raw effects and filled coral diamonds show Technical-adjusted effects. Horizontal segments show the interquartile range, and short ticks mark the median. Technical adjustment was performed using the prespecified ARTEMIS technical covariates before calculating the patient-level High − Other Program146 effect.",
  "",
  "All testable patient effects were positive under every definition (pooled: 77/77 Raw and 77/77 Technical-adjusted; within-patient: 78/78 and 78/78; patient-z: 78/78 and 78/78). Technical adjustment attenuated effect magnitude while preserving its positive direction. P59 was non-testable under the pooled definition because it contained 50 High cells and no Other cells; it was not plotted as zero and was excluded from the pooled median and interquartile range. Program146 coverage remained 141 of 146 genes throughout.",
  "",
  "Program146 effects remained directionally recurrent across three prespecified YAP-Stem state definitions, while technical adjustment attenuated effect magnitude."
)
write_utf8_copy(legend_lines, file.path(out_dir, "S2C_legend_draft.md"))

qc_lines <- c(
  "# Supplementary Figure S2C candidate QC",
  "",
  "## Scope",
  "",
  "- Mechanical visualization only; no state definition, Program146 score, patient effect, or statistical test was recomputed.",
  "- This script exports the standalone state-definition sensitivity panel only.",
  "",
  "## Authoritative input integrity",
  "",
  sprintf("- Source: `%s`.", unique(d$source_file)),
  sprintf("- Source SHA-256: `%s`.", unique(d$source_sha256)),
  "- Source structure: 474 records = 468 patient records + 6 formal summary records - PASS.",
  "- Rendered patient points: 466 testable effects - PASS.",
  "- Program146 coverage: 141/146 genes in every lane; unchanged and not rendered as a separate visual track - PASS.",
  "",
  "## Frozen numerical reproduction",
  "",
  sprintf("- Primary pooled Raw: %d/%d positive; median = %.6f - PASS.", summary_data$positive_n[1], summary_data$testable_n[1], summary_data$lane_median[1]),
  sprintf("- Primary pooled Technical-adjusted: %d/%d positive; median = %.6f - PASS.", summary_data$positive_n[2], summary_data$testable_n[2], summary_data$lane_median[2]),
  sprintf("- Within-patient Raw: %d/%d positive; median = %.6f - PASS.", summary_data$positive_n[3], summary_data$testable_n[3], summary_data$lane_median[3]),
  sprintf("- Within-patient Technical-adjusted: %d/%d positive; median = %.6f - PASS.", summary_data$positive_n[4], summary_data$testable_n[4], summary_data$lane_median[4]),
  sprintf("- Patient-z Raw: %d/%d positive; median = %.6f - PASS.", summary_data$positive_n[5], summary_data$testable_n[5], summary_data$lane_median[5]),
  sprintf("- Patient-z Technical-adjusted: %d/%d positive; median = %.6f - PASS.", summary_data$positive_n[6], summary_data$testable_n[6], summary_data$lane_median[6]),
  "- Negative testable effects: 0; zero testable effects: 0 - PASS.",
  "- P59 pooled Raw and adjusted: NON_TESTABLE (High=50, Other=0), missing effect retained, not plotted, and excluded from formal summaries - PASS.",
  "",
  "## Visual architecture",
  "",
  "- Three state-definition rows with adjacent Raw and Technical-adjusted sublanes - PASS.",
  "- One common effect axis across all six lanes with a visible effect=0 reference - PASS.",
  "- Every testable patient is shown with deterministic perpendicular jitter; effect values are not jittered - PASS.",
  "- Median and IQR are shown for each lane; no box rectangles, density envelopes, P values, significance marks, or patient-level connectors were added - PASS.",
  "- Raw uses the project blue open circle and Technical-adjusted uses the project coral filled diamond - PASS.",
  "- Clinical-forest design reference: adaptation is limited to aligned rows, a shared quantitative axis, and compact summary marks; no hazard-ratio, logarithmic-axis, confidence-interval, or significance semantics were imported.",
  "",
  "## Export and true-size checks",
  "",
  "- Target dimensions: 88 x 78 mm.",
  sprintf("- Full patient-effect range: %.6f to %.6f; shared display range: %.6f to %.6f; no clipping or winsorization.", min(patient_plot$Program146_effect), max(patient_plot$Program146_effect), x_limits[1], x_limits[2]),
  "- PDF requested as vector and PNG requested at 600 dpi; external PDF/font/DPI/raster validation is recorded below after rendering.",
  "- Core plotted text is >=6.2 pt; title 7.6 pt; panel letter 11 pt.",
  "- External visual inspection for clipping and unintended overlap is recorded below after rendering."
)
write_utf8_copy(qc_lines, file.path(out_dir, "S2C_state_definition_QC.md"))

cat("S2C candidate rendered\n")
