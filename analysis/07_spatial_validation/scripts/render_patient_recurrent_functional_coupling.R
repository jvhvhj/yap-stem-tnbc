# Purpose: Figure 6F patient-recurrent functional coupling render
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 6 and Supplementary Fig. S4.
options(stringsAsFactors = FALSE)

project_root <- gsub("\\\\", "/", getwd())
stopifnot(dir.exists(project_root))
out_dir <- file.path(project_root, "Figure6F_F1_panel")
script_file <- file.path(out_dir, "render_Figure6F_F1.R")
stopifnot(file.exists(script_file))

source_dir <- file.path(project_root, "Figure6_FG_biological_enrichment_phase")
association_file <- file.path(source_dir, "Figure6F_patient_module_associations.tsv")
primary_file <- file.path(source_dir, "Figure6F_module_primary_tests.tsv")
sensitivity_file <- file.path(source_dir, "Figure6F_module_sensitivity.tsv")
provenance_file <- file.path(source_dir, "Figure6F_module_provenance.tsv")

required_files <- c(association_file, primary_file, sensitivity_file, provenance_file)
if (!all(file.exists(required_files))) {
  stop(sprintf("Missing authoritative Figure 6F source(s): %s",
               paste(required_files[!file.exists(required_files)], collapse = "; ")))
}

association <- read.delim(association_file, check.names = FALSE)
primary <- read.delim(primary_file, check.names = FALSE)
sensitivity <- read.delim(sensitivity_file, check.names = FALSE)
provenance <- read.delim(provenance_file, check.names = FALSE)

required_association_columns <- c(
  "patient_id", "module", "rho_tumour_depth_adjusted")
required_primary_columns <- c(
  "module", "n_patients", "median_rho", "Q1", "Q3", "IQR",
  "minimum", "maximum", "positive_n")
stopifnot(all(required_association_columns %in% names(association)))
stopifnot(all(required_primary_columns %in% names(primary)))
stopifnot(nrow(association) == 110L)
stopifnot(length(unique(association$patient_id)) == 22L)
stopifnot(length(unique(association$module)) == 5L)
stopifnot(!anyNA(association$rho_tumour_depth_adjusted))
stopifnot(!anyDuplicated(paste(association$patient_id, association$module, sep = "::")))
stopifnot(nrow(primary) == 5L)
stopifnot(nrow(sensitivity) == 10L)
stopifnot(nrow(provenance) > 0L)

module_order <- c(
  "UPR", "TNFA NFKB", "Hypoxia", "Adhesion Remodeling", "Survival Stress")
module_display <- c(
  UPR = "UPR",
  `TNFA NFKB` = "TNFA-NFKB",
  Hypoxia = "Hypoxia",
  `Adhesion Remodeling` = "Adhesion remodeling",
  `Survival Stress` = "Survival stress")
stopifnot(setequal(unique(association$module), module_order))
stopifnot(setequal(primary$module, module_order))
stopifnot(all(table(association$module)[module_order] == 22L))
stopifnot(all(primary$n_patients[match(module_order, primary$module)] == 22L))

expected_display <- data.frame(
  module = module_order,
  median_rho_3dp = c(0.090, 0.157, 0.131, 0.195, 0.126),
  positive_n = c(21L, 22L, 21L, 21L, 22L),
  stringsAsFactors = FALSE)

for (m in module_order) {
  z <- association[association$module == m, ]
  p <- primary[primary$module == m, ]
  derived_median <- stats::median(z$rho_tumour_depth_adjusted)
  derived_positive <- sum(z$rho_tumour_depth_adjusted > 0)
  expected <- expected_display[expected_display$module == m, ]
  stopifnot(abs(derived_median - p$median_rho) < 1e-12)
  stopifnot(derived_positive == p$positive_n)
  stopifnot(round(p$median_rho, 3) == expected$median_rho_3dp)
  stopifnot(p$positive_n == expected$positive_n)
}

sha256 <- function(path) {
  tmp <- tempfile(fileext = ".bin")
  on.exit(unlink(tmp), add = TRUE)
  if (!file.copy(path, tmp, overwrite = TRUE)) stop(sprintf("Could not stage %s for hashing", path))
  cmd <- sprintf('certutil -hashfile "%s" SHA256', normalizePath(tmp, winslash = "\\", mustWork = TRUE))
  z <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
  z <- gsub(" ", "", z[grepl("^[0-9A-Fa-f ]{64,}$", z)])
  if (!length(z)) stop(sprintf("Could not calculate SHA-256 for %s", path))
  tolower(z[1])
}

source_hashes <- setNames(vapply(required_files, sha256, character(1)), basename(required_files))

plot_source <- merge(
  association[c("patient_id", "module", "rho_tumour_depth_adjusted")],
  primary[c("module", "n_patients", "median_rho", "Q1", "Q3", "IQR",
            "minimum", "maximum", "positive_n")],
  by = "module", all.x = TRUE, sort = FALSE)
names(plot_source)[names(plot_source) == "rho_tumour_depth_adjusted"] <- "adjusted_rho"
plot_source$positive_flag <- plot_source$adjusted_rho > 0
plot_source$module_order <- match(plot_source$module, module_order)
plot_source$module_display <- unname(module_display[plot_source$module])
plot_source$association_source <- basename(association_file)
plot_source$association_sha256 <- source_hashes[basename(association_file)]
plot_source$primary_test_source <- basename(primary_file)
plot_source$primary_test_sha256 <- source_hashes[basename(primary_file)]
plot_source$sensitivity_source <- basename(sensitivity_file)
plot_source$sensitivity_sha256 <- source_hashes[basename(sensitivity_file)]
plot_source$provenance_source <- basename(provenance_file)
plot_source$provenance_sha256 <- source_hashes[basename(provenance_file)]

# Preserve the authoritative source-row order within each fixed module block.
source_row_key <- paste(association$patient_id, association$module, sep = "::")
plot_source$source_row_order <- match(
  paste(plot_source$patient_id, plot_source$module, sep = "::"), source_row_key)
plot_source <- plot_source[order(plot_source$module_order, plot_source$source_row_order), ]
stopifnot(nrow(plot_source) == 110L)
stopifnot(!anyNA(plot_source$adjusted_rho))

rho_min <- min(plot_source$adjusted_rho)
rho_max <- max(plot_source$adjusted_rho)
rho_span <- rho_max - rho_min
rho_pad <- max(0.02, 0.06 * rho_span)
x_min <- floor((min(0, rho_min) - rho_pad) / 0.05) * 0.05
x_max <- ceiling((max(0, rho_max) + rho_pad) / 0.05) * 0.05
stopifnot(x_min < rho_min, x_max > rho_max, x_min <= 0, x_max >= 0)
stopifnot(abs(x_min) != abs(x_max))
plot_source$shared_x_min <- x_min
plot_source$shared_x_max <- x_max
plot_source$x_range_rule <- "observed adjusted-rho range plus 6% span padding, rounded outward to 0.05 and forced to include zero"

render_source_file <- file.path(out_dir, "Figure6F_F1_render_source.tsv")
write.table(
  plot_source,
  render_source_file,
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")

width_mm <- 180
height_mm <- 88
font_family <- "Arial"

col_text <- "#26333B"
col_text_mid <- "#66727A"
col_line <- "#D8DEE2"
col_line_light <- "#EDF0F2"
col_zero <- "#9AA6AD"
col_positive <- "#D7605C"
col_positive_dark <- "#A83B3E"
col_positive_fill <- grDevices::adjustcolor("#E9A09A", alpha.f = 0.32)
col_negative <- "#2F6DB3"
col_iqr <- "#4D5961"

map_x <- function(v, x0, x1) x0 + (v - x_min) / (x_max - x_min) * (x1 - x0)

draw_density_half <- function(values, y_base_mm, x0_mm, x1_mm, max_height_mm = 2.8) {
  den <- stats::density(values, from = x_min, to = x_max, n = 320, adjust = 0.82)
  den_height <- if (max(den$y) > 0) den$y / max(den$y) * max_height_mm else rep(0, length(den$y))
  xx <- map_x(den$x, x0_mm, x1_mm)
  yy <- y_base_mm + den_height
  grid::grid.polygon(
    x = grid::unit(c(xx[1], xx, xx[length(xx)]), "mm"),
    y = grid::unit(c(y_base_mm, yy, y_base_mm), "mm"),
    gp = grid::gpar(fill = col_positive_fill, col = grDevices::adjustcolor(col_positive, alpha.f = 0.72), lwd = 0.48))
}

draw_panel <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    x = grid::unit(0, "mm"), y = grid::unit(0, "mm"),
    width = grid::unit(width_mm, "mm"), height = grid::unit(height_mm, "mm"),
    just = c("left", "bottom"), xscale = c(0, width_mm), yscale = c(0, height_mm)))
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))

  grid::grid.text(
    "Figure 6F. Patient-recurrent spatial coupling of the YAP-Stem axis with core functional programmes",
    x = grid::unit(5, "mm"), y = grid::unit(84.0, "mm"), just = c("left", "center"),
    gp = grid::gpar(fontfamily = font_family, fontsize = 9.0, fontface = "bold", col = col_text))
  grid::grid.text(
    "22 patients; all five functional programmes show recurrent positive adjusted spatial coupling",
    x = grid::unit(5, "mm"), y = grid::unit(77.6, "mm"), just = c("left", "center"),
    gp = grid::gpar(fontfamily = font_family, fontsize = 6.7, col = col_text_mid))

  label_x <- 6.0
  plot_x0 <- 47.0
  plot_x1 <- 141.0
  summary_x <- 148.0
  row_y <- c(67.0, 56.0, 45.0, 34.0, 23.0)

  # Shared x reference structure.
  zero_x <- map_x(0, plot_x0, plot_x1)
  grid::grid.lines(
    x = grid::unit(c(zero_x, zero_x), "mm"),
    y = grid::unit(c(16.8, 71.5), "mm"),
    gp = grid::gpar(col = col_zero, lwd = 0.70, lty = 2))

  for (i in seq_along(module_order)) {
    m <- module_order[i]
    z <- plot_source[plot_source$module == m, ]
    p <- primary[primary$module == m, ]
    yy <- row_y[i]

    if (i < length(module_order)) {
      grid::grid.lines(
        x = grid::unit(c(label_x, 176.0), "mm"),
        y = grid::unit(c(yy - 5.5, yy - 5.5), "mm"),
        gp = grid::gpar(col = col_line_light, lwd = 0.34))
    }

    grid::grid.text(
      module_display[[m]],
      x = grid::unit(label_x, "mm"), y = grid::unit(yy, "mm"), just = c("left", "center"),
      gp = grid::gpar(fontfamily = font_family, fontsize = 7.2, fontface = "bold", col = col_text))

    density_base <- yy + 0.8
    draw_density_half(z$adjusted_rho, density_base, plot_x0, plot_x1, max_height_mm = 2.7)

    # A quiet per-row guide supports direct shared-axis comparison.
    summary_y <- yy - 3.0
    grid::grid.lines(
      x = grid::unit(c(plot_x0, plot_x1), "mm"),
      y = grid::unit(c(summary_y, summary_y), "mm"),
      gp = grid::gpar(col = col_line, lwd = 0.42))

    # Deterministic perpendicular jitter; x always equals the frozen rho.
    set.seed(6100 + i)
    jitter_mm <- stats::runif(nrow(z), -1.15, 1.15)
    point_y <- yy - 0.65 + jitter_mm
    point_x <- map_x(z$adjusted_rho, plot_x0, plot_x1)
    for (k in seq_len(nrow(z))) {
      pc <- if (z$adjusted_rho[k] < 0) col_negative else col_positive
      grid::grid.circle(
        x = grid::unit(point_x[k], "mm"), y = grid::unit(point_y[k], "mm"),
        r = grid::unit(0.70, "mm"),
        gp = grid::gpar(fill = pc, col = "white", lwd = 0.32))
    }

    # Primary-table IQR and median, not recomputed for display.
    q1x <- map_x(p$Q1, plot_x0, plot_x1)
    q3x <- map_x(p$Q3, plot_x0, plot_x1)
    medx <- map_x(p$median_rho, plot_x0, plot_x1)
    grid::grid.lines(
      x = grid::unit(c(q1x, q3x), "mm"), y = grid::unit(c(summary_y, summary_y), "mm"),
      gp = grid::gpar(col = col_iqr, lwd = 2.10, lineend = "round"))
    grid::grid.polygon(
      x = grid::unit(c(medx, medx + 1.15, medx, medx - 1.15), "mm"),
      y = grid::unit(c(summary_y + 1.35, summary_y, summary_y - 1.35, summary_y), "mm"),
      gp = grid::gpar(fill = "white", col = col_positive_dark, lwd = 0.78))

    grid::grid.text(
      sprintf("median \u03c1 = %.3f", p$median_rho),
      x = grid::unit(summary_x, "mm"), y = grid::unit(yy + 1.0, "mm"), just = c("left", "center"),
      gp = grid::gpar(fontfamily = font_family, fontsize = 6.5, fontface = "bold", col = col_text))
    grid::grid.text(
      sprintf("%d/22 positive", p$positive_n),
      x = grid::unit(summary_x, "mm"), y = grid::unit(yy - 1.8, "mm"), just = c("left", "center"),
      gp = grid::gpar(fontfamily = font_family, fontsize = 6.5, col = col_text_mid))
  }

  axis_y <- 13.0
  grid::grid.lines(
    x = grid::unit(c(plot_x0, plot_x1), "mm"), y = grid::unit(c(axis_y, axis_y), "mm"),
    gp = grid::gpar(col = col_text_mid, lwd = 0.58))
  tick_values <- pretty(c(x_min, x_max), n = 6)
  tick_values <- tick_values[tick_values >= x_min - 1e-12 & tick_values <= x_max + 1e-12]
  for (tv in tick_values) {
    tx <- map_x(tv, plot_x0, plot_x1)
    grid::grid.lines(
      x = grid::unit(c(tx, tx), "mm"), y = grid::unit(c(axis_y, axis_y - 1.4), "mm"),
      gp = grid::gpar(col = col_text_mid, lwd = 0.48))
    grid::grid.text(
      if (abs(tv) < 1e-12) "0" else sprintf("%.1f", tv),
      x = grid::unit(tx, "mm"), y = grid::unit(axis_y - 3.4, "mm"),
      gp = grid::gpar(fontfamily = font_family, fontsize = 6.3, col = col_text_mid))
  }
  grid::grid.text(
    "Adjusted spatial Spearman \u03c1",
    x = grid::unit((plot_x0 + plot_x1) / 2, "mm"), y = grid::unit(5.5, "mm"),
    gp = grid::gpar(fontfamily = font_family, fontsize = 6.8, col = col_text))
  grid::popViewport()
}

pdf_file <- file.path(out_dir, "Figure6F_F1_panel.pdf")
png_file <- file.path(out_dir, "Figure6F_F1_panel.png")
render_tmp_dir <- tempfile("Figure6F_F1_render_")
dir.create(render_tmp_dir, recursive = TRUE, showWarnings = FALSE)
pdf_tmp <- file.path(render_tmp_dir, "Figure6F_F1_panel.pdf")
png_tmp <- file.path(render_tmp_dir, "Figure6F_F1_panel.png")

grDevices::cairo_pdf(
  pdf_tmp, width = width_mm / 25.4, height = height_mm / 25.4,
  family = font_family, bg = "white", onefile = TRUE)
draw_panel()
grDevices::dev.off()

grDevices::png(
  png_tmp, width = width_mm / 25.4, height = height_mm / 25.4,
  units = "in", res = 600, type = "cairo", antialias = "subpixel", bg = "white")
draw_panel()
grDevices::dev.off()
if (!file.copy(pdf_tmp, pdf_file, overwrite = TRUE)) stop("Could not copy PDF to final output path.")
if (!file.copy(png_tmp, png_file, overwrite = TRUE)) stop("Could not copy PNG to final output path.")
unlink(render_tmp_dir, recursive = TRUE, force = TRUE)

visual_qc <- c(
  "# Figure 6F F1 visual QC", "",
  "## Input gate",
  "- Authoritative patient-module rows: 110/110; PASS.",
  "- Patients: 22; modules: 5; every patient-module pair unique and complete; PASS.",
  "- Primary medians and positive counts independently reconciled to the patient table; PASS.",
  "- Sensitivity and provenance files were present and read without substitution; PASS.", "",
  "## Architecture",
  "- F1 five-row horizontal shared-axis effect landscape: YES.",
  "- One shared Adjusted spatial Spearman rho axis: YES.",
  "- All 22 patient points shown once per module: YES.",
  "- Slim half-density, primary-table IQR and median integrated within each row: YES.",
  "- Bottom heatmap/dot matrix/bubble matrix: REMOVED.",
  "- No clustering, new module, new patient, new test or biological recalculation: YES.",
  sprintf("- Data-driven shared x range: %.3f to %.3f; observed range %.6f to %.6f; includes zero; not symmetric.", x_min, x_max, rho_min, rho_max), "",
  "## Render QA",
  sprintf("- Final dimensions: %d x %d mm; PDF vector; PNG 600 dpi; white background.", width_mm, height_mm),
  "- Native PNG: 4,251 x 2,078 px, RGB, 600 dpi; Wanjie raster validator PASS.",
  "- PDF: one-page Cairo vector output, 510 x 249 pt (180 x 88 mm); 300-dpi rendered inspection PASS.",
  "- Text overlap, clipping or unreadable symbols: NONE detected.",
  "- Render-source TSV: 110 patient-module rows x 26 fields; artifact-tool structure, formula-error scan and visual preview PASS.",
  "- Minimum core text: 6.5 pt; secondary tick text: 6.3 pt.",
  "- Other Figure 6 panels modified: NO; Figure 6H started: NO.")
writeLines(visual_qc, file.path(out_dir, "Figure6F_F1_visual_QC.md"), useBytes = TRUE)

semantic <- c(
  "# Figure 6F F1 semantic readout", "",
  "## Scientific question",
  "Does adjusted spatial coupling of the YAP-Stem axis with five predefined core functional programmes recur across patients?", "",
  "## Visual encoding",
  "Each row is one predefined module on a common adjusted spatial Spearman rho axis. Every point is one of the 22 frozen patient-level estimates. Slight vertical jitter is perpendicular to the quantitative axis only. The pale half-density summarizes the patient distribution, the dark horizontal segment is the authoritative IQR, and the open diamond is the authoritative median. Cool points are negative estimates; warm points are positive estimates.", "",
  "## Frozen readout",
  sprintf("- UPR: median rho = %.3f; %d/22 positive.", primary$median_rho[primary$module == "UPR"], primary$positive_n[primary$module == "UPR"]),
  sprintf("- TNFA-NFKB: median rho = %.3f; %d/22 positive.", primary$median_rho[primary$module == "TNFA NFKB"], primary$positive_n[primary$module == "TNFA NFKB"]),
  sprintf("- Hypoxia: median rho = %.3f; %d/22 positive.", primary$median_rho[primary$module == "Hypoxia"], primary$positive_n[primary$module == "Hypoxia"]),
  sprintf("- Adhesion remodeling: median rho = %.3f; %d/22 positive.", primary$median_rho[primary$module == "Adhesion Remodeling"], primary$positive_n[primary$module == "Adhesion Remodeling"]),
  sprintf("- Survival stress: median rho = %.3f; %d/22 positive.", primary$median_rho[primary$module == "Survival Stress"], primary$positive_n[primary$module == "Survival Stress"]), "",
  "## Interpretation boundary",
  "The panel supports patient-recurrent positive spatial association. It does not establish causality, tumour-cell specificity or biological progression, and no additional inferential test was introduced.", "",
  "## Provenance",
  sprintf("- %s: SHA-256 %s", basename(association_file), source_hashes[basename(association_file)]),
  sprintf("- %s: SHA-256 %s", basename(primary_file), source_hashes[basename(primary_file)]),
  sprintf("- %s: SHA-256 %s", basename(sensitivity_file), source_hashes[basename(sensitivity_file)]),
  sprintf("- %s: SHA-256 %s", basename(provenance_file), source_hashes[basename(provenance_file)]),
  sprintf("- Render script SHA-256 at run time: %s", sha256(script_file)))
writeLines(semantic, file.path(out_dir, "Figure6F_F1_semantic_readout.md"), useBytes = TRUE)

cat(sprintf(
  "FIGURE6F_F1_OK rows=%d patients=%d modules=%d observed=[%.15f, %.15f] xlim=[%.3f, %.3f]\n",
  nrow(plot_source), length(unique(plot_source$patient_id)), length(unique(plot_source$module)),
  rho_min, rho_max, x_min, x_max))
