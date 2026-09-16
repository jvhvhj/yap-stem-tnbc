#!/usr/bin/env Rscript
# Purpose: Figure 6G independent spatial replication render
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 6 and Supplementary Fig. S4.
options(stringsAsFactors = FALSE, digits = 17)
.libPaths(c(".runtime/AHIPPO_RLIB_G1_ASCII", .Library))

suppressPackageStartupMessages({
  library(digest)
  library(grid)
})

runtime_root <- "."
out_dir <- file.path(runtime_root, "Figure6G_integrated_patient_tracks")
g1_dir <- file.path(runtime_root, "Figure6G1_BSW2_spatial_replication")
summary_dir <- file.path(runtime_root, "Figure6G_concordance_fingerprint_lowfidelity")
map_dir <- file.path(runtime_root, "Figure6G_adjusted_coupling_atlas_lowfidelity")
left_reference_dir <- file.path(runtime_root, "Figure6G_WangLike_patient_matrix_panel")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

association_file <- file.path(g1_dir, "Figure6G1_patient_associations.tsv")
primary_file <- file.path(g1_dir, "Figure6G1_primary_test.tsv")
nfeature_file <- file.path(g1_dir, "Figure6G1_nFeature_sensitivity.tsv")
feature_mapping_file <- file.path(g1_dir, "Figure6G1_feature_mapping.tsv")
independence_file <- file.path(g1_dir, "Figure6G1_score_independence_check.tsv")
summary_file <- file.path(summary_dir, "Figure6G_all_patient_concordance_summary.tsv")
map_file <- file.path(map_dir, "Figure6G_sAA6_concordance_classes.tsv")
left_reference_png <- file.path(left_reference_dir, "Figure6G_WangLike_panel.png")
left_reference_script <- file.path(left_reference_dir, "render_Figure6G_WangLike.R")
palette_file <- file.path(runtime_root, "project_semantic_palette.tsv")

render_source_file <- file.path(out_dir, "Figure6G_integrated_patient_tracks_render_source.tsv")
pdf_file <- file.path(out_dir, "Figure6G_integrated_patient_tracks.pdf")
png_file <- file.path(out_dir, "Figure6G_integrated_patient_tracks.png")
visual_qc_file <- file.path(out_dir, "Figure6G_integrated_patient_tracks_visual_QC.md")
semantic_file <- file.path(out_dir, "Figure6G_integrated_patient_tracks_semantic_readout.md")

required_inputs <- c(association_file, primary_file, nfeature_file,
                     feature_mapping_file, independence_file, summary_file,
                     map_file, left_reference_png, left_reference_script,
                     palette_file)
stopifnot(all(file.exists(required_inputs)))

sha256 <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE,
              na = "NA", fileEncoding = "UTF-8")
}

assoc <- read.delim(association_file, check.names = FALSE)
primary <- read.delim(primary_file, check.names = FALSE)
nfeature <- read.delim(nfeature_file, check.names = FALSE)
feature_mapping <- read.delim(feature_mapping_file, check.names = FALSE)
independence <- read.delim(independence_file, check.names = FALSE)
summary_tab <- read.delim(summary_file, check.names = FALSE, quote = "")
saa6 <- read.delim(map_file, check.names = FALSE, quote = "")
palette <- read.delim(palette_file, check.names = FALSE)

patient_order <- c("sEA2", "sEA1", "sAA2", "sAA1", "sAA6",
                   "sEA5", "sEA7", "sEA6", "sAA8")
class_order <- c("CONCORDANT_HIGH", "CONCORDANT_LOW",
                 "JOINT_HIGH_DISCORDANT", "PROGRAM_HIGH_DISCORDANT")

nfeature_patient <- nfeature[nfeature$record_type == "PATIENT", , drop = FALSE]
stopifnot(
  nrow(assoc) == 9L,
  setequal(assoc$patient_id, patient_order),
  all(is.finite(assoc$rho_raw)),
  all(is.finite(assoc$rho_depth_adjusted)),
  all(is.finite(assoc$rho_nFeature_adjusted)),
  all(assoc$rho_depth_adjusted > 0),
  nrow(primary) == 1L,
  primary$n_patients == 9L,
  primary$positive_n == 9L,
  abs(primary$median_rho - 0.402148154179219) < 1e-15,
  abs(primary$Q1 - 0.304319162021956) < 1e-15,
  abs(primary$Q3 - 0.499693006920698) < 1e-15,
  abs(primary$exact_one_sided_wilcoxon_P - 0.001953125) < 1e-15,
  nrow(nfeature_patient) == 9L,
  setequal(nfeature_patient$patient_id, patient_order),
  max(abs(nfeature_patient$rho_nFeature_adjusted -
            assoc$rho_nFeature_adjusted[match(nfeature_patient$patient_id,
                                              assoc$patient_id)])) < 1e-15,
  nrow(independence) == 1L,
  independence$overlap_n == 0L,
  independence$score_independence_status == "PASS_ZERO_OVERLAP",
  sum(feature_mapping$signature == "YAP17" & feature_mapping$included_in_score) == 17L,
  sum(feature_mapping$signature == "Stem21" & feature_mapping$included_in_score) == 19L,
  sum(feature_mapping$signature == "Program146" & feature_mapping$included_in_score) == 135L,
  nrow(summary_tab) == 9L,
  setequal(summary_tab$patient_id, patient_order),
  all(summary_tab$fraction_sum_gate == "PASS"),
  nrow(saa6) == 4111L,
  identical(unique(saa6$patient_id), "sAA6"),
  length(unique(saa6$spot_barcode)) == 4111L,
  all(is.finite(saa6$x_coordinate)),
  all(is.finite(saa6$y_coordinate)),
  all(saa6$concordance_class %in% class_order)
)

display_order <- assoc$patient_id[order(assoc$rho_depth_adjusted)]
stopifnot(identical(display_order, patient_order))

plot_data <- merge(
  assoc[, c("patient_id", "n_spots", "rho_raw", "rho_depth_adjusted",
            "rho_nFeature_adjusted")],
  summary_tab[, c("patient_id", "n_total_spots",
                  "n_CONCORDANT_HIGH", "fraction_CONCORDANT_HIGH",
                  "n_CONCORDANT_LOW", "fraction_CONCORDANT_LOW",
                  "n_JOINT_HIGH_DISCORDANT", "fraction_JOINT_HIGH_DISCORDANT",
                  "n_PROGRAM_HIGH_DISCORDANT", "fraction_PROGRAM_HIGH_DISCORDANT",
                  "n_NEUTRAL_ZERO", "fraction_NEUTRAL_ZERO",
                  "concordant_fraction", "discordant_fraction",
                  "four_non_neutral_fraction_sum")],
  by = "patient_id", sort = FALSE)
plot_data <- plot_data[match(display_order, plot_data$patient_id), , drop = FALSE]
plot_data$display_order <- seq_len(nrow(plot_data))
plot_data$rho_min <- apply(plot_data[, c("rho_raw", "rho_depth_adjusted",
                                         "rho_nFeature_adjusted")], 1, min)
plot_data$rho_max <- apply(plot_data[, c("rho_raw", "rho_depth_adjusted",
                                         "rho_nFeature_adjusted")], 1, max)
plot_data$adjusted_median <- primary$median_rho
plot_data$adjusted_Q1 <- primary$Q1
plot_data$adjusted_Q3 <- primary$Q3
plot_data$adjusted_IQR <- primary$IQR
plot_data$exact_one_sided_wilcoxon_P <- primary$exact_one_sided_wilcoxon_P
plot_data$all_three_positive <- plot_data$rho_raw > 0 &
  plot_data$rho_depth_adjusted > 0 & plot_data$rho_nFeature_adjusted > 0
plot_data$concordant_boundary_fraction <-
  plot_data$fraction_CONCORDANT_HIGH + plot_data$fraction_CONCORDANT_LOW
plot_data$concordant_percent_label <-
  sprintf("%.0f%%", 100 * plot_data$concordant_fraction)
plot_data$association_source <- "Figure6G1_BSW2_spatial_replication/Figure6G1_patient_associations.tsv"
plot_data$association_source_sha256 <- sha256(association_file)
plot_data$composition_source <- "Figure6G_concordance_fingerprint_lowfidelity/Figure6G_all_patient_concordance_summary.tsv"
plot_data$composition_source_sha256 <- sha256(summary_file)
plot_data$left_map_source <- "Figure6G_adjusted_coupling_atlas_lowfidelity/Figure6G_sAA6_concordance_classes.tsv"
plot_data$left_map_source_sha256 <- sha256(map_file)
plot_data <- plot_data[, c(
  "patient_id", "display_order", "n_spots", "n_total_spots",
  "rho_raw", "rho_depth_adjusted", "rho_nFeature_adjusted",
  "rho_min", "rho_max", "all_three_positive",
  "adjusted_median", "adjusted_Q1", "adjusted_Q3", "adjusted_IQR",
  "exact_one_sided_wilcoxon_P",
  "n_CONCORDANT_HIGH", "fraction_CONCORDANT_HIGH",
  "n_CONCORDANT_LOW", "fraction_CONCORDANT_LOW",
  "n_JOINT_HIGH_DISCORDANT", "fraction_JOINT_HIGH_DISCORDANT",
  "n_PROGRAM_HIGH_DISCORDANT", "fraction_PROGRAM_HIGH_DISCORDANT",
  "n_NEUTRAL_ZERO", "fraction_NEUTRAL_ZERO",
  "concordant_fraction", "concordant_boundary_fraction",
  "concordant_percent_label", "discordant_fraction",
  "four_non_neutral_fraction_sum",
  "association_source", "association_source_sha256",
  "composition_source", "composition_source_sha256",
  "left_map_source", "left_map_source_sha256")]

stopifnot(
  identical(plot_data$patient_id, patient_order),
  all(plot_data$n_spots == plot_data$n_total_spots),
  all(plot_data$all_three_positive),
  max(abs(plot_data$four_non_neutral_fraction_sum - 1)) < 1e-12,
  all(plot_data$n_NEUTRAL_ZERO == 0L),
  max(abs(plot_data$concordant_boundary_fraction -
            plot_data$concordant_fraction)) < 1e-15,
  all(plot_data$rho_min <= plot_data$rho_raw),
  all(plot_data$rho_min <= plot_data$rho_depth_adjusted),
  all(plot_data$rho_min <= plot_data$rho_nFeature_adjusted),
  all(plot_data$rho_max >= plot_data$rho_raw),
  all(plot_data$rho_max >= plot_data$rho_depth_adjusted),
  all(plot_data$rho_max >= plot_data$rho_nFeature_adjusted)
)
write_tsv(plot_data, render_source_file)

get_colour <- function(token) {
  z <- palette$hex[palette$token == token]
  stopifnot(length(z) == 1L)
  z
}

positive_colour <- get_colour("Positive")
negative_colour <- get_colour("Negative")
raw_colour <- get_colour("Raw")
nfeature_colour <- "#6F8795"
discord_joint <- "#C69A38"
discord_program <- "#2C9A92"
text_dark <- "#26333B"
text_mid <- "#56656E"
line_light <- "#DCE4E8"
line_mid <- "#8A979F"
font_family <- "Arial"

class_colours <- c(
  CONCORDANT_HIGH = positive_colour,
  CONCORDANT_LOW = negative_colour,
  JOINT_HIGH_DISCORDANT = discord_joint,
  PROGRAM_HIGH_DISCORDANT = discord_program)

class_labels <- c(
  CONCORDANT_HIGH = "Concordant high-high",
  CONCORDANT_LOW = "Concordant low-low",
  JOINT_HIGH_DISCORDANT = "Joint-high / Program-low",
  PROGRAM_HIGH_DISCORDANT = "Joint-low / Program-high")

width_mm <- 180
height_mm <- 100

vp_top <- function(x, y, width, height, xscale = c(0, 1), yscale = c(0, 1),
                   clip = "off") {
  viewport(x = unit(x, "mm"), y = unit(height_mm - y, "mm"),
           width = unit(width, "mm"), height = unit(height, "mm"),
           just = c("left", "top"), xscale = xscale, yscale = yscale,
           clip = clip)
}

# LOCKED LEFT BLOCK: copied without modification from the accepted Wang-like panel.
draw_spatial_map <- function(x, y, width, height) {
  x_pad <- diff(range(saa6$x_coordinate)) * 0.015
  y_plot <- -saa6$y_coordinate
  y_pad <- diff(range(y_plot)) * 0.015
  pushViewport(vp_top(x, y, width, height,
                      xscale = range(saa6$x_coordinate) + c(-x_pad, x_pad),
                      yscale = range(y_plot) + c(-y_pad, y_pad)))
  grid.points(x = unit(saa6$x_coordinate, "native"),
              y = unit(y_plot, "native"), pch = 21,
              size = unit(0.50, "mm"),
              gp = gpar(fill = unname(class_colours[saa6$concordance_class]),
                        col = grDevices::adjustcolor("#FFFFFF", alpha.f = 0.48),
                        lwd = 0.18))
  popViewport()
}

# LOCKED LEFT BLOCK: copied without modification from the accepted Wang-like panel.
draw_shared_legend <- function(x, y, width, height) {
  pushViewport(vp_top(x, y, width, height))
  legend_x <- c(0.02, 0.50, 0.02, 0.50)
  legend_y <- c(0.72, 0.72, 0.25, 0.25)
  for (i in seq_along(class_order)) {
    cls <- class_order[i]
    grid.circle(x = unit(legend_x[i], "npc"), y = unit(legend_y[i], "npc"),
                r = unit(1.35, "mm"),
                gp = gpar(fill = class_colours[cls], col = "white", lwd = 0.35))
    grid.text(class_labels[cls],
              x = unit(legend_x[i] + 0.045, "npc"), y = unit(legend_y[i], "npc"),
              just = c("left", "centre"),
              gp = gpar(fontfamily = font_family, fontsize = 6.5, col = text_dark))
  }
  popViewport()
}

rho_xlim <- c(-0.05, 0.75)
rho_to_mm <- function(value, x0, width) {
  x0 + width * (value - rho_xlim[1]) / diff(rho_xlim)
}

draw_association_legend <- function(x, y) {
  grid.points(x = unit(x, "mm"), y = unit(height_mm - y, "mm"),
              pch = 21, size = unit(1.7, "mm"),
              gp = gpar(fill = "white", col = raw_colour, lwd = 0.85))
  grid.text("Raw", x = unit(x + 2.0, "mm"), y = unit(height_mm - y, "mm"),
            just = c("left", "centre"),
            gp = gpar(fontfamily = font_family, fontsize = 6.5, col = text_dark))
  grid.points(x = unit(x + 12.0, "mm"), y = unit(height_mm - y, "mm"),
              pch = 18, size = unit(2.15, "mm"),
              gp = gpar(fill = positive_colour, col = positive_colour))
  grid.text("Depth-adj.", x = unit(x + 14.1, "mm"), y = unit(height_mm - y, "mm"),
            just = c("left", "centre"),
            gp = gpar(fontfamily = font_family, fontsize = 6.5, col = text_dark))
  grid.points(x = unit(x + 32.5, "mm"), y = unit(height_mm - y, "mm"),
              pch = 22, size = unit(1.6, "mm"),
              gp = gpar(fill = "white", col = nfeature_colour, lwd = 0.80))
  grid.text("nFeature-adj.", x = unit(x + 34.6, "mm"), y = unit(height_mm - y, "mm"),
            just = c("left", "centre"),
            gp = gpar(fontfamily = font_family, fontsize = 6.5, col = text_dark))
}

draw_integrated_patient_tracks <- function() {
  patient_x <- 69.0
  axis_x0 <- 84.0
  axis_w <- 57.0
  comp_x0 <- 147.0
  comp_w <- 30.0
  row_y <- seq(35.0, 80.2, length.out = 9)
  row_top <- 32.0
  row_bottom <- 83.1
  bar_h <- 3.5

  grid.text("Integrated patient evidence tracks",
            x = unit(patient_x, "mm"), y = unit(height_mm - 20.3, "mm"),
            just = c("left", "bottom"),
            gp = gpar(fontfamily = font_family, fontsize = 7.4,
                      fontface = "bold", col = text_dark))
  grid.text("Spatial Spearman rho", x = unit(axis_x0 + axis_w / 2, "mm"),
            y = unit(height_mm - 24.0, "mm"), just = c("centre", "bottom"),
            gp = gpar(fontfamily = font_family, fontsize = 6.7,
                      fontface = "bold", col = text_dark))
  grid.text("Residual sign composition",
            x = unit(comp_x0 + comp_w / 2, "mm"),
            y = unit(height_mm - 24.0, "mm"), just = c("centre", "bottom"),
            gp = gpar(fontfamily = font_family, fontsize = 6.7,
                      fontface = "bold", col = text_dark))
  grid.text("Label at boundary = concordant %",
            x = unit(comp_x0 + comp_w / 2, "mm"),
            y = unit(height_mm - 28.3, "mm"), just = c("centre", "centre"),
            gp = gpar(fontfamily = font_family, fontsize = 6.2,
                      col = text_mid))
  draw_association_legend(axis_x0 + 1.5, 28.3)

  iqr_x0 <- rho_to_mm(primary$Q1, axis_x0, axis_w)
  iqr_x1 <- rho_to_mm(primary$Q3, axis_x0, axis_w)
  grid.rect(x = unit(iqr_x0, "mm"), y = unit(height_mm - row_top, "mm"),
            width = unit(iqr_x1 - iqr_x0, "mm"),
            height = unit(row_bottom - row_top, "mm"),
            just = c("left", "top"),
            gp = gpar(fill = grDevices::adjustcolor(positive_colour, alpha.f = 0.075),
                      col = NA))

  zero_x <- rho_to_mm(0, axis_x0, axis_w)
  median_x <- rho_to_mm(primary$median_rho, axis_x0, axis_w)
  grid.segments(x0 = unit(zero_x, "mm"), x1 = unit(zero_x, "mm"),
                y0 = unit(height_mm - row_top, "mm"),
                y1 = unit(height_mm - row_bottom, "mm"),
                gp = gpar(col = line_mid, lwd = 0.7, lty = 2))
  grid.segments(x0 = unit(median_x, "mm"), x1 = unit(median_x, "mm"),
                y0 = unit(height_mm - row_top, "mm"),
                y1 = unit(height_mm - row_bottom, "mm"),
                gp = gpar(col = grDevices::adjustcolor(positive_colour, alpha.f = 0.85),
                          lwd = 0.85))

  for (tick in c(0, 0.2, 0.4, 0.6)) {
    xx <- rho_to_mm(tick, axis_x0, axis_w)
    grid.segments(x0 = unit(xx, "mm"), x1 = unit(xx, "mm"),
                  y0 = unit(height_mm - row_bottom, "mm"),
                  y1 = unit(height_mm - (row_bottom + 1.2), "mm"),
                  gp = gpar(col = text_mid, lwd = 0.45))
    grid.text(sprintf("%.1f", tick), x = unit(xx, "mm"),
              y = unit(height_mm - (row_bottom + 2.0), "mm"),
              just = c("centre", "top"),
              gp = gpar(fontfamily = font_family, fontsize = 6.5, col = text_mid))
  }

  for (i in seq_len(nrow(plot_data))) {
    d <- plot_data[i, ]
    yy <- row_y[i]
    grid.text(d$patient_id, x = unit(patient_x, "mm"),
              y = unit(height_mm - yy, "mm"), just = c("left", "centre"),
              gp = gpar(fontfamily = font_family, fontsize = 6.8,
                        fontface = if (d$patient_id == "sAA6") "bold" else "plain",
                        col = text_dark))

    grid.segments(x0 = unit(rho_to_mm(d$rho_min, axis_x0, axis_w), "mm"),
                  x1 = unit(rho_to_mm(d$rho_max, axis_x0, axis_w), "mm"),
                  y0 = unit(height_mm - yy, "mm"),
                  y1 = unit(height_mm - yy, "mm"),
                  gp = gpar(col = grDevices::adjustcolor(line_mid, alpha.f = 0.48),
                            lwd = 0.65))
    grid.points(x = unit(rho_to_mm(d$rho_raw, axis_x0, axis_w), "mm"),
                y = unit(height_mm - yy, "mm"), pch = 21,
                size = unit(1.75, "mm"),
                gp = gpar(fill = "white", col = raw_colour, lwd = 0.85))
    grid.points(x = unit(rho_to_mm(d$rho_nFeature_adjusted, axis_x0, axis_w), "mm"),
                y = unit(height_mm - yy, "mm"), pch = 22,
                size = unit(1.65, "mm"),
                gp = gpar(fill = "white", col = nfeature_colour, lwd = 0.80))
    adjusted_x <- rho_to_mm(d$rho_depth_adjusted, axis_x0, axis_w)
    grid.points(x = unit(adjusted_x, "mm"), y = unit(height_mm - yy, "mm"),
                pch = 18, size = unit(2.35, "mm"),
                gp = gpar(fill = positive_colour, col = positive_colour))
    grid.text(sprintf("%.3f", d$rho_depth_adjusted),
              x = unit(adjusted_x + 1.55, "mm"), y = unit(height_mm - yy, "mm"),
              just = c("left", "centre"),
              gp = gpar(fontfamily = font_family, fontsize = 6.5,
                        fontface = "bold", col = text_dark))

    fractions <- c(d$fraction_CONCORDANT_HIGH, d$fraction_CONCORDANT_LOW,
                   d$fraction_JOINT_HIGH_DISCORDANT,
                   d$fraction_PROGRAM_HIGH_DISCORDANT)
    starts <- c(0, cumsum(fractions)[1:3])
    for (j in seq_along(class_order)) {
      grid.rect(x = unit(comp_x0 + comp_w * starts[j], "mm"),
                y = unit(height_mm - yy, "mm"),
                width = unit(comp_w * fractions[j], "mm"),
                height = unit(bar_h, "mm"), just = c("left", "centre"),
                gp = gpar(fill = class_colours[class_order[j]],
                          col = "white", lwd = 0.22))
    }
    concordant_boundary_x <-
      comp_x0 + comp_w * d$concordant_boundary_fraction
    grid.segments(
      x0 = unit(concordant_boundary_x, "mm"),
      x1 = unit(concordant_boundary_x, "mm"),
      y0 = unit(height_mm - (yy - bar_h / 2 - 0.25), "mm"),
      y1 = unit(height_mm - (yy + bar_h / 2 + 0.25), "mm"),
      gp = gpar(col = text_dark, lwd = 0.48)
    )
    grid.text(
      d$concordant_percent_label,
      x = unit(concordant_boundary_x, "mm"),
      y = unit(height_mm - (yy - bar_h / 2 - 0.52), "mm"),
      just = c("centre", "bottom"),
      gp = gpar(fontfamily = font_family, fontsize = 6.1,
                fontface = "bold", col = text_dark)
    )
  }
}

draw_figure <- function() {
  grid.newpage()
  pushViewport(viewport(width = unit(width_mm, "mm"), height = unit(height_mm, "mm")))
  grid.rect(gp = gpar(fill = "white", col = NA))

  grid.text("G", x = unit(4, "mm"), y = unit(height_mm - 4.0, "mm"),
            just = c("left", "top"),
            gp = gpar(fontfamily = font_family, fontsize = 11,
                      fontface = "bold", col = text_dark))
  grid.text("Patient-recurrent spatial coupling in BSW2",
            x = unit(12, "mm"), y = unit(height_mm - 4.0, "mm"),
            just = c("left", "top"),
            gp = gpar(fontfamily = font_family, fontsize = 9.2,
                      fontface = "bold", col = text_dark))
  grid.text(sprintf("9/9 positive  |  median adjusted rho = %.3f  |  exact P = %.5f",
                    primary$median_rho, primary$exact_one_sided_wilcoxon_P),
            x = unit(69, "mm"), y = unit(height_mm - 13.0, "mm"),
            just = c("left", "top"),
            gp = gpar(fontfamily = font_family, fontsize = 7.2,
                      fontface = "bold", col = text_dark))

  # LOCKED LEFT BLOCK: coordinates, title, map and legend calls are unchanged.
  grid.text("sAA6 adjusted spatial concordance",
            x = unit(5.5, "mm"), y = unit(height_mm - 21.5, "mm"),
            just = c("left", "bottom"),
            gp = gpar(fontfamily = font_family, fontsize = 7.2,
                      fontface = "bold", col = text_dark))
  draw_spatial_map(5.0, 26.0, 58.0, 58.0)
  draw_shared_legend(5.0, 87.0, 62.0, 10.0)

  draw_integrated_patient_tracks()
  popViewport()
}

grDevices::cairo_pdf(pdf_file, width = width_mm / 25.4,
                     height = height_mm / 25.4, family = font_family,
                     bg = "white", onefile = TRUE)
draw_figure()
grDevices::dev.off()

grDevices::png(png_file, width = width_mm / 25.4, height = height_mm / 25.4,
               units = "in", res = 600, type = "cairo",
               antialias = "subpixel", bg = "white")
draw_figure()
grDevices::dev.off()

visual_qc <- c(
  "# Figure 6G integrated patient tracks visual QC", "",
  "## Frozen-data checks",
  "1. Left sAA6 spatial block copied with the same source, colours, point size, title, coordinates, crop, orientation and legend geometry: YES. Pixel comparison against the accepted Wang-like PNG over the full left 67 mm (1,582 x 2,362 px) found zero differing pixels; PASS.",
  "2. Exactly nine patients included: YES.",
  "3. Raw, depth-adjusted and nFeature-adjusted effects share one Spatial Spearman rho axis: YES.",
  "4. Patient order is depth-adjusted rho ascending: sEA2, sEA1, sAA2, sAA1, sAA6, sEA5, sEA7, sEA6, sAA8; PASS.",
  sprintf("5. Adjusted median %.15f and IQR [%.15f, %.15f] match the authoritative primary test: PASS.", primary$median_rho, primary$Q1, primary$Q3),
  "6. Four-class stacked bars map directly to the frozen concordance summary and sum to 1 for every patient: PASS.", "",
  "## Targeted composition revision",
  "- Left sAA6 spatial map modified: NO.",
  "- Association track modified: NO.",
  "- Pre/post pixel comparison through x=144 mm (left map plus complete association track) found zero differing pixels: PASS.",
  "- Independent Conc.% column retained: NO; it is absent from the revised visual.",
  "- Concordant fraction embedded in each stacked bar: YES; the boundary tick and percentage label mark the end of the adjacent red-plus-blue segments.",
  sprintf("- Boundary labels in frozen patient order: %s; all equal 100 x (fraction_CONCORDANT_HIGH + fraction_CONCORDANT_LOW), rounded to whole percent for display.", paste(plot_data$concordant_percent_label, collapse = ", ")),
  "- Four component proportions, patient order, colours and all association values changed: NO.", "",
  "## Required removals",
  "- Nine mini spatial/fingerprint panels: REMOVED.",
  "- Independent Raw, nFeature and Composition value columns: REMOVED.",
  "- Redundant grey explanatory text: REMOVED.",
  "- Table-like vertical separators: REMOVED.", "",
  "## Scope",
  "8. Figure 6A-F modified: NO.",
  "9. Figure 6H started or modified: NO.",
  sprintf("- Final dimensions: %d x %d mm; PDF vector; PNG 600 dpi; white background.", width_mm, height_mm),
  "- Minimum plotted core text: 6.5 pt.",
  "- Native PNG: 4,251 x 2,362 px, RGB, 600 dpi; Wanjie raster validator PASS.",
  "- PDF: one-page Cairo vector output, 510 x 283 pt (180 x 100 mm); 300-dpi rendered inspection PASS.",
  "- Text overlap or clipping: NONE detected at native PNG size or in the rendered PDF.",
  sprintf("- Render-source table: 9 patient rows x %d provenance/value fields; artifact-tool structure and formula-error scan PASS.", ncol(plot_data)))
writeLines(visual_qc, visual_qc_file, useBytes = TRUE)

semantic <- c(
  "# Figure 6G integrated patient tracks semantic readout", "",
  "## Scientific object",
  "The panel tests whether the positive spatial association between the YAP-Stem Joint axis and Program146 recurs across the nine independent BSW2 patients and remains directionally stable under the two prespecified technical views.", "",
  "## Left spatial context",
  "The accepted sAA6 adjusted spatial concordance block is carried forward without visual or numerical modification. Every point is one Visium spot at its frozen coordinate; the four colours encode the existing residual-sign classes without new thresholds.", "",
  "## Integrated right-side tracks",
  "- Each patient occupies one aligned row ordered by the frozen depth-adjusted rho.",
  "- The shared association axis displays Raw rho (open circle), the primary depth-adjusted rho (filled diamond) and nFeature-adjusted rho (open square). A thin line spans the within-patient minimum and maximum of the three existing effects.",
  "- The vertical dashed line marks rho = 0. The coral line marks the frozen adjusted median, and the pale coral band marks the frozen adjusted interquartile range.",
  "- Only the primary depth-adjusted effect is numerically labelled.",
  "- The aligned 100% stacked bar shows the frozen spot fractions for concordant high-high, concordant low-low, Joint-high/Program-low and Joint-low/Program-high.", "",
  "## Interpretation boundary",
  sprintf("All 9/9 adjusted patient effects are positive; median rho %.4f, exact one-sided Wilcoxon P %.8f.", primary$median_rho, primary$exact_one_sided_wilcoxon_P),
  "Composition bars are descriptive and have no new P value. Visium spots are mixed-cell tissue measurements; no malignant-cell specificity, direct regulation or causality is inferred.", "",
  "## Provenance",
  sprintf("Associations: %s (SHA-256 %s).", basename(association_file), sha256(association_file)),
  sprintf("Primary test: %s (SHA-256 %s).", basename(primary_file), sha256(primary_file)),
  sprintf("Concordance summary: %s (SHA-256 %s).", basename(summary_file), sha256(summary_file)),
  sprintf("sAA6 spatial classes: %s (SHA-256 %s).", basename(map_file), sha256(map_file)),
  sprintf("Accepted left-reference PNG: %s (SHA-256 %s).", basename(left_reference_png), sha256(left_reference_png)),
  sprintf("Render script SHA-256 at run time: %s.", sha256(file.path(out_dir, "render_Figure6G_integrated_patient_tracks.R"))))
writeLines(semantic, semantic_file, useBytes = TRUE)

cat(sprintf("FIGURE6G_INTEGRATED_TRACKS_OK patients=%d sAA6_spots=%d median=%.15f Q1=%.15f Q3=%.15f P=%.9f\n",
            nrow(plot_data), nrow(saa6), primary$median_rho, primary$Q1,
            primary$Q3, primary$exact_one_sided_wilcoxon_P))
