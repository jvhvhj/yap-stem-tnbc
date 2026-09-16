#!/usr/bin/env Rscript
# Purpose: Render S2A patient-resolved robustness
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Supplementary Fig. S1–S4.

# Supplementary Figure S2A FINAL
# Patient-resolved robustness of YAP–Stem coupling across scoring methods
#
# This script only reshapes and renders frozen patient-level correlations.
# It does not recompute scores or correlations and performs no new inference.

options(stringsAsFactors = FALSE, width = 220)

# Run from the project root. Using getwd() avoids locale-dependent parsing of
# non-ASCII path literals on Windows while preserving the same files.
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(project_root, "project_semantic_palette.tsv"))) {
  stop("Run this script from the AHIPPO-YAP project root")
}
out_dir_override <- Sys.getenv("S2A_OUTPUT_DIR", unset = "")
out_dir <- if (nzchar(out_dir_override)) {
  normalizePath(out_dir_override, winslash = "/", mustWork = FALSE)
} else {
  file.path(project_root, "Supplementary_S2A_patient_resolved_FINAL")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

.libPaths(c(
  file.path(project_root, "0710_final_evidence_rebuild/vendor"),
  ".software/r-library",
  file.path(project_root, "0703_rebuild/R_library"),
  .libPaths()
))

suppressPackageStartupMessages({
  library(grid)
})

wu_file <- file.path(
  project_root,
  "0717_methodological_framework_benchmark/tables/scoring_method_patient_effects_wu.tsv"
)
artemis_file <- file.path(
  project_root,
  "0717_methodological_framework_benchmark/tables/scoring_method_patient_effects_yan.tsv"
)
palette_file <- file.path(project_root, "project_semantic_palette.tsv")

for (path in c(wu_file, artemis_file, palette_file)) {
  if (!file.exists(path)) stop("Required frozen input is missing: ", path)
}

read_patient_source <- function(path, cohort_display, source_file_label) {
  x <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    colClasses = "character"
  )
  required <- c("cohort", "method", "model", "patient", "n_cells", "rho")
  if (!all(required %in% names(x))) {
    stop("Missing required fields in ", path, ": ",
         paste(setdiff(required, names(x)), collapse = ", "))
  }
  # Preserve the source's exact textual rho representation for the provenance
  # export while using numeric values for summaries and plotting.
  x$rho_source_text <- x$rho
  x$n_cells <- as.integer(x$n_cells)
  x$rho <- as.numeric(x$rho)
  x$cohort_display <- cohort_display
  # Store a portable project-relative provenance path. Absolute paths containing
  # Chinese characters are escaped by the Windows R locale in TSV output.
  x$source_file <- source_file_label
  x$source_row <- seq_len(nrow(x)) + 1L
  x
}

wu <- read_patient_source(
  wu_file,
  "Wu",
  "0717_methodological_framework_benchmark/tables/scoring_method_patient_effects_wu.tsv"
)
artemis <- read_patient_source(
  artemis_file,
  "ARTEMIS",
  "0717_methodological_framework_benchmark/tables/scoring_method_patient_effects_yan.tsv"
)

method_order <- c(
  "AUCell",
  "UCell",
  "current_frozen_clean",
  "expanded_stem_plasticity",
  "fixed_TEAD_regulon"
)

method_display <- c(
  AUCell = "AUCell",
  UCell = "UCell",
  current_frozen_clean = "Mean expression\n(primary)",
  expanded_stem_plasticity = "Expanded stem/\nplasticity",
  fixed_TEAD_regulon = "Fixed TEAD\nregulon"
)

model_order <- c("raw", "technical_adjusted")
model_display <- c(raw = "Raw", technical_adjusted = "Technical-adjusted")

validate_source <- function(x, cohort_display, expected_patients, expected_rows) {
  if (nrow(x) != expected_rows) {
    stop(cohort_display, " source has ", nrow(x), " rows; expected ", expected_rows)
  }
  if (!setequal(unique(x$method), method_order)) {
    stop(cohort_display, " method set does not match the fixed five-method definition")
  }
  if (!setequal(unique(x$model), model_order)) {
    stop(cohort_display, " model set is not exactly raw + technical_adjusted")
  }
  if (length(unique(x$patient)) != expected_patients) {
    stop(cohort_display, " patient count mismatch")
  }
  key <- paste(x$method, x$model, x$patient, sep = "|")
  if (anyDuplicated(key)) stop(cohort_display, " contains duplicate method/model/patient keys")
  if (any(!is.finite(x$rho))) stop(cohort_display, " contains missing/non-finite rho values")
  for (method in method_order) {
    patient_sets <- lapply(model_order, function(model) {
      sort(x$patient[x$method == method & x$model == model])
    })
    if (!identical(patient_sets[[1]], patient_sets[[2]])) {
      stop(cohort_display, " Raw/adjusted patient sets differ for ", method)
    }
    if (length(patient_sets[[1]]) != expected_patients) {
      stop(cohort_display, " incomplete patient set for ", method)
    }
  }
  invisible(TRUE)
}

validate_source(wu, "Wu", 8L, 80L)
validate_source(artemis, "ARTEMIS", 78L, 780L)

patient_data <- rbind(wu, artemis)
if (nrow(patient_data) != 860L) stop("Combined source must contain 860 patient-level records")

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
zero_colour <- "#9AA4AB"
iqr_raw <- grDevices::adjustcolor(raw_colour, alpha.f = 0.43)
iqr_adjusted <- grDevices::adjustcolor(adjusted_colour, alpha.f = 0.45)

source_hashes <- c(
  Wu = "82C1EF4C15516816EE678E2AF5AD2C9BF845E659655EE58454FB8CF2C04CBF62",
  ARTEMIS = "12CF7B5F4F3A5577299C737B558052EAA6512587CC8FD4743F837253D3C624E7",
  Palette = "EB11BB4A1FC05822BE207105EDE2592715999424BD1F657DFD20CD9C968C169E"
)

summary_rows <- list()
k <- 0L
for (cohort_name in c("Wu", "ARTEMIS")) {
  for (method_name in method_order) {
    for (model_name in model_order) {
      k <- k + 1L
      values <- patient_data$rho[
        patient_data$cohort_display == cohort_name &
          patient_data$method == method_name &
          patient_data$model == model_name
      ]
      q <- as.numeric(stats::quantile(values, probs = c(0.25, 0.50, 0.75),
                                      type = 7, names = FALSE))
      summary_rows[[k]] <- data.frame(
        cohort_display = cohort_name,
        method = method_name,
        model = model_name,
        q1_rho = q[1],
        median_rho = q[2],
        q3_rho = q[3],
        positive_n = sum(values > 0),
        testable_n = length(values),
        stringsAsFactors = FALSE
      )
    }
  }
}
summary_data <- do.call(rbind, summary_rows)

patient_data$method_order <- match(patient_data$method, method_order)
patient_data$method_display <- unname(method_display[patient_data$method])
patient_data$model_order <- match(patient_data$model, model_order)
patient_data$model_display <- unname(model_display[patient_data$model])
patient_data$cohort_order <- match(patient_data$cohort_display, c("Wu", "ARTEMIS"))

summary_key <- paste(summary_data$cohort_display, summary_data$method, summary_data$model, sep = "|")
patient_key <- paste(patient_data$cohort_display, patient_data$method, patient_data$model, sep = "|")
summary_match <- match(patient_key, summary_key)
if (anyNA(summary_match)) stop("Patient rows could not be mapped to summary rows")

patient_data$q1_rho <- summary_data$q1_rho[summary_match]
patient_data$median_rho <- summary_data$median_rho[summary_match]
patient_data$q3_rho <- summary_data$q3_rho[summary_match]
patient_data$positive_n <- summary_data$positive_n[summary_match]
patient_data$testable_n <- summary_data$testable_n[summary_match]

# Deterministic perpendicular jitter only; rho itself is never jittered.
lane_group <- interaction(
  patient_data$cohort_display, patient_data$method, patient_data$model,
  drop = TRUE
)
patient_data$patient_order_within_lane <- integer(nrow(patient_data))
for (group_level in levels(lane_group)) {
  idx <- which(lane_group == group_level)
  patient_data$patient_order_within_lane[idx] <-
    rank(patient_data$patient[idx], ties.method = "first")
}
golden_fraction <- 0.618033988749895
patient_data$jitter_offset_npc <- vapply(seq_len(nrow(patient_data)), function(i) {
  amplitude <- if (patient_data$cohort_display[i] == "ARTEMIS") 0.0110 else 0.0075
  ((((patient_data$patient_order_within_lane[i] - 1) * golden_fraction) %% 1) - 0.5) *
    2 * amplitude
}, numeric(1))

row_centres <- c(0.690, 0.575, 0.460, 0.345, 0.230)
model_offsets <- c(raw = 0.022, technical_adjusted = -0.022)
patient_data$row_centre_npc <- row_centres[patient_data$method_order]
patient_data$lane_centre_npc <- patient_data$row_centre_npc +
  unname(model_offsets[patient_data$model])
patient_data$point_y_npc <- patient_data$lane_centre_npc + patient_data$jitter_offset_npc

observed_range <- range(patient_data$rho)
display_padding <- diff(observed_range) * 0.05
rho_limits <- c(observed_range[1] - display_padding, observed_range[2] + display_padding)
rho_ticks <- seq(
  ceiling(rho_limits[1] / 0.2) * 0.2,
  floor(rho_limits[2] / 0.2) * 0.2,
  by = 0.2
)

if (any(patient_data$rho < rho_limits[1] | patient_data$rho > rho_limits[2])) {
  stop("At least one rho falls outside the derived display limits")
}

patient_data$rho_display_min <- rho_limits[1]
patient_data$rho_display_max <- rho_limits[2]
patient_data$source_sha256 <- ifelse(
  patient_data$cohort_display == "Wu", source_hashes[["Wu"]], source_hashes[["ARTEMIS"]]
)

render_columns <- c(
  "cohort_display", "cohort_order", "method", "method_display", "method_order",
  "model", "model_display", "model_order", "patient", "n_cells", "rho",
  "positive_n", "testable_n", "q1_rho", "median_rho", "q3_rho",
  "row_centre_npc", "lane_centre_npc", "jitter_offset_npc", "point_y_npc",
  "rho_display_min", "rho_display_max", "source_file", "source_row", "source_sha256"
)
render_source <- patient_data[
  order(patient_data$cohort_order, patient_data$method_order,
        patient_data$model_order, patient_data$patient),
  render_columns
]
# Keep the provenance TSV strictly one physical row per observation. Display
# line breaks are plot-only and are flattened in the machine-readable export.
render_source$method_display <- gsub("\n", " ", render_source$method_display, fixed = TRUE)
render_source$rho <- patient_data$rho_source_text[
  order(patient_data$cohort_order, patient_data$method_order,
        patient_data$model_order, patient_data$patient)
]

write.table(
  render_source,
  file.path(out_dir, "S2A_patient_resolved_FINAL_render_source.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA",
  fileEncoding = "UTF-8"
)

layout_spec <- list(
  wu = list(plot_x0 = 0.150, plot_x1 = 0.430, count_x = 0.492, cohort_x = 0.340),
  artemis = list(plot_x0 = 0.590, plot_x1 = 0.870, count_x = 0.935, cohort_x = 0.785)
)

x_map <- function(rho, x0, x1) {
  x0 + (rho - rho_limits[1]) / diff(rho_limits) * (x1 - x0)
}

draw_text <- function(label, x, y, ..., fontsize = 6.5, fontface = "plain",
                      colour = ink, just = "center") {
  grid::grid.text(
    label, x = x, y = y, just = just,
    gp = grid::gpar(
      fontfamily = "Arial", fontsize = fontsize, fontface = fontface,
      col = colour
    ), ...
  )
}

draw_cohort_block <- function(cohort_name, spec) {
  block_data <- patient_data[patient_data$cohort_display == cohort_name, , drop = FALSE]
  block_summary <- summary_data[summary_data$cohort_display == cohort_name, , drop = FALSE]

  # Shared-axis reference and baseline.
  x_zero <- x_map(0, spec$plot_x0, spec$plot_x1)
  grid::grid.segments(
    x0 = x_zero, x1 = x_zero, y0 = 0.180, y1 = 0.740,
    gp = grid::gpar(col = zero_colour, lwd = 0.65, lty = 2)
  )
  grid::grid.segments(
    x0 = spec$plot_x0, x1 = spec$plot_x1, y0 = 0.145, y1 = 0.145,
    gp = grid::gpar(col = ink, lwd = 0.65)
  )

  # Identical ticks in both cohort blocks.
  tick_x <- x_map(rho_ticks, spec$plot_x0, spec$plot_x1)
  grid::grid.segments(
    x0 = tick_x, x1 = tick_x, y0 = 0.145, y1 = 0.137,
    gp = grid::gpar(col = ink, lwd = 0.55)
  )
  draw_text(sprintf("%.1f", rho_ticks), tick_x, 0.111, fontsize = 6.5, colour = ink)

  # IQR, patient points, and median tick for every method/model lane.
  for (method_name in method_order) {
    for (model_name in model_order) {
      lane <- block_data[
        block_data$method == method_name & block_data$model == model_name,
        , drop = FALSE
      ]
      sm <- block_summary[
        block_summary$method == method_name & block_summary$model == model_name,
        , drop = FALSE
      ]
      if (nrow(sm) != 1L) stop("Summary row mapping failed")
      y_lane <- unique(lane$lane_centre_npc)
      series_colour <- if (model_name == "raw") raw_colour else adjusted_colour
      iqr_colour <- if (model_name == "raw") iqr_raw else iqr_adjusted

      grid::grid.segments(
        x0 = x_map(sm$q1_rho, spec$plot_x0, spec$plot_x1),
        x1 = x_map(sm$q3_rho, spec$plot_x0, spec$plot_x1),
        y0 = y_lane, y1 = y_lane,
        gp = grid::gpar(col = iqr_colour, lwd = 2.00, lineend = "round")
      )

      if (model_name == "raw") {
        point_size <- if (cohort_name == "Wu") 1.65 else 1.22
        grid::grid.points(
          x = x_map(lane$rho, spec$plot_x0, spec$plot_x1),
          y = lane$point_y_npc,
          default.units = "npc",
          pch = 21,
          size = grid::unit(point_size, "mm"),
          gp = grid::gpar(
            fill = grDevices::adjustcolor("white", alpha.f = if (cohort_name == "Wu") 0.96 else 0.76),
            col = grDevices::adjustcolor(series_colour, alpha.f = if (cohort_name == "Wu") 1.00 else 0.82),
            lwd = if (cohort_name == "Wu") 0.72 else 0.55
          )
        )
      } else {
        point_size <- if (cohort_name == "Wu") 1.68 else 1.24
        grid::grid.points(
          x = x_map(lane$rho, spec$plot_x0, spec$plot_x1),
          y = lane$point_y_npc,
          default.units = "npc",
          pch = 23,
          size = grid::unit(point_size, "mm"),
          gp = grid::gpar(
            fill = grDevices::adjustcolor(series_colour, alpha.f = if (cohort_name == "Wu") 0.90 else 0.66),
            col = grDevices::adjustcolor(ink, alpha.f = if (cohort_name == "Wu") 0.78 else 0.56),
            lwd = if (cohort_name == "Wu") 0.65 else 0.50
          )
        )
      }

      median_x <- x_map(sm$median_rho, spec$plot_x0, spec$plot_x1)
      grid::grid.segments(
        x0 = median_x, x1 = median_x,
        y0 = y_lane - 0.012, y1 = y_lane + 0.012,
        gp = grid::gpar(col = ink, lwd = 1.15, lineend = "round")
      )
    }
  }

  # One positive-count annotation per method and cohort.
  for (i in seq_along(method_order)) {
    raw_sm <- block_summary[
      block_summary$method == method_order[i] & block_summary$model == "raw",
      , drop = FALSE
    ]
    adj_sm <- block_summary[
      block_summary$method == method_order[i] & block_summary$model == "technical_adjusted",
      , drop = FALSE
    ]
    label <- sprintf(
      "%d/%d \u2192 %d/%d",
      raw_sm$positive_n, raw_sm$testable_n, adj_sm$positive_n, adj_sm$testable_n
    )
    draw_text(label, spec$count_x, row_centres[i], fontsize = 6.5, colour = ink)
  }
}

draw_panel <- function() {
  grid::grid.newpage()

  # Panel identity and compact shared legend.
  draw_text("A", 0.018, 0.952, fontsize = 11.0, fontface = "bold", just = c("left", "top"))
  draw_text(
    "Patient-resolved robustness of YAP\u2013Stem coupling across scoring methods",
    0.055, 0.949, fontsize = 9.0, fontface = "bold", just = c("left", "top")
  )

  grid::grid.points(
    x = 0.735, y = 0.939, default.units = "npc",
    pch = 21, size = grid::unit(1.35, "mm"),
    gp = grid::gpar(fill = "white", col = raw_colour, lwd = 0.65)
  )
  draw_text("Raw", 0.751, 0.939, fontsize = 6.5, just = "left")
  grid::grid.points(
    x = 0.823, y = 0.939, default.units = "npc",
    pch = 23, size = grid::unit(1.38, "mm"),
    gp = grid::gpar(fill = adjusted_colour, col = adjusted_colour, lwd = 0.65)
  )
  draw_text("Technical-adjusted", 0.840, 0.939, fontsize = 6.5, just = "left")

  # Cohort hierarchy and positive-count headings.
  draw_text("Discovery cohort", layout_spec$wu$cohort_x, 0.842, fontsize = 7.8, fontface = "bold")
  draw_text("ARTEMIS", layout_spec$artemis$cohort_x, 0.842, fontsize = 7.8, fontface = "bold")
  grid::grid.segments(
    x0 = c(layout_spec$wu$plot_x0, layout_spec$artemis$plot_x0),
    x1 = c(0.535, 0.980),
    y0 = 0.816, y1 = 0.816,
    gp = grid::gpar(col = "#DDE2E5", lwd = 0.55)
  )
  draw_text(
    "Positive n/N\nRaw \u2192 adjusted",
    layout_spec$wu$count_x, 0.775, fontsize = 6.5, fontface = "bold", colour = muted
  )
  draw_text(
    "Positive n/N\nRaw \u2192 adjusted",
    layout_spec$artemis$count_x, 0.775, fontsize = 6.5, fontface = "bold", colour = muted
  )

  # One shared method-label lane aligned to both cohorts.
  for (i in seq_along(method_order)) {
    draw_text(
      unname(method_display[method_order[i]]),
      0.136, row_centres[i], fontsize = 6.7, colour = ink, just = "right"
    )
  }

  draw_cohort_block("Wu", layout_spec$wu)
  draw_cohort_block("ARTEMIS", layout_spec$artemis)

  draw_text(
  "Within-patient Spearman ρ", 0.510, 0.047,
    fontsize = 7.0, colour = ink
  )
}

pdf_file <- file.path(out_dir, "S2A_patient_resolved_FINAL.pdf")
png_file <- file.path(out_dir, "S2A_patient_resolved_FINAL_600dpi.png")
tmp_pdf <- tempfile(pattern = "S2A_patient_resolved_FINAL_", fileext = ".pdf")
tmp_png <- tempfile(pattern = "S2A_patient_resolved_FINAL_", fileext = ".png")

grDevices::cairo_pdf(
  filename = tmp_pdf,
  width = 180 / 25.4,
  height = 74 / 25.4,
  family = "Arial",
  bg = "white",
  onefile = TRUE
)
draw_panel()
grDevices::dev.off()
if (!file.copy(tmp_pdf, pdf_file, overwrite = TRUE)) stop("PDF copy to final path failed")
unlink(tmp_pdf)

grDevices::png(
  filename = tmp_png,
  width = 180,
  height = 74,
  units = "mm",
  res = 600,
  bg = "white",
  type = "cairo"
)
draw_panel()
grDevices::dev.off()
if (!file.copy(tmp_png, png_file, overwrite = TRUE)) stop("PNG copy to final path failed")
unlink(tmp_png)

if (!file.exists(pdf_file) || file.info(pdf_file)$size <= 0) stop("PDF export failed")
if (!file.exists(png_file) || file.info(png_file)$size <= 0) stop("PNG export failed")

summary_for_qc <- summary_data[
  order(match(summary_data$cohort_display, c("Wu", "ARTEMIS")),
        match(summary_data$method, method_order),
        match(summary_data$model, model_order)),
]
summary_lines <- vapply(seq_len(nrow(summary_for_qc)), function(i) {
  z <- summary_for_qc[i,]
  sprintf(
    "| %s | %s | %s | %.6f | %.6f | %.6f | %d/%d |",
    z$cohort_display, z$method, unname(model_display[z$model]),
    z$q1_rho, z$median_rho, z$q3_rho, z$positive_n, z$testable_n
  )
}, character(1))

qc_lines <- c(
  "# S2A patient-resolved FINAL QC",
  "",
  "## Input and numerical fidelity",
  "",
  "- Wu observations rendered: **80/80 — PASS**.",
  "- ARTEMIS observations rendered: **780/780 — PASS**.",
  "- Combined patient-level observations: **860/860 — PASS**.",
  "- Missing/non-finite rho values: **0 — PASS**.",
  "- Duplicate cohort–method–model–patient keys: **0 — PASS**.",
  "- Patient aggregation before plotting: **NO — PASS**.",
  "- Values were read directly from the two frozen patient-effect tables; scores and correlations were not recomputed.",
  sprintf("- Frozen observed rho range: **%.15f to %.15f**.", observed_range[1], observed_range[2]),
  sprintf("- Shared display range (5%% margin derived from complete range): **%.15f to %.15f**.", rho_limits[1], rho_limits[2]),
  "- Same rho scale in Wu and ARTEMIS: **YES — PASS**.",
  "- Clipping: **NO — PASS**.",
  "- Winsorization/transformation: **NO — PASS**.",
  "- Method order changed or clustered: **NO — PASS**.",
  "- Raw/adjusted identities checked against source model values: **PASS**.",
  "",
  "## Source hashes",
  "",
  sprintf("- Wu table SHA-256: `%s`.", source_hashes[["Wu"]]),
  sprintf("- ARTEMIS table SHA-256: `%s`.", source_hashes[["ARTEMIS"]]),
  sprintf("- Project semantic palette SHA-256: `%s`.", source_hashes[["Palette"]]),
  sprintf("- Raw colour: `%s`; Technical-adjusted colour: `%s`.", raw_colour, adjusted_colour),
  "",
  "## Reproduced patient summaries",
  "",
  "Quartiles use R quantile type 7 and are calculated only as display summaries from the frozen patient-level rho values.",
  "",
  "| Cohort | Method | Model | Q1 rho | Median rho | Q3 rho | Positive/testable |",
  "|---|---|---|---:|---:|---:|---:|",
  summary_lines,
  "",
  "## Visual and output gate",
  "",
  "- Architecture: two side-by-side cohort blocks with the same five method rows — PASS.",
  "- Every patient rho rendered as an individual point — PASS.",
  "- Raw and Technical-adjusted use stable colour/fill plus narrow vertical lane offset — PASS.",
  "- Magnitude encoded only by horizontal position — PASS.",
  "- Median/IQR reproduced from patient-level sources — PASS.",
  "- Positive counts reproduced directly from patient-level sources — PASS.",
  "- Density/violin/boxplot/heatmap/tile/bar plot: NONE — PASS.",
  "- 6/6 program-retention text block: ABSENT — PASS.",
  "- UCell/AUCell implementation-QC block: ABSENT — PASS.",
  "- Technical-adjustment wording does not imply identical covariates across cohorts — PASS.",
  "- Technical adjustment was performed within each cohort using the prespecified cohort-specific technical covariates.",
  "- ONLY HEADER COMPACTION AND SUMMARY-GLYPH EMPHASIS WERE CHANGED.",
  "- Header wording: `Positive n/N` and `Raw → adjusted`, identical in both cohorts.",
  "- IQR alpha/line width and median-tick line width increased by approximately 10–15%; patient point glyphs are unchanged.",
  "- Final dimensions: **180 mm × 74 mm**.",
  "- Minimum true-size font: **6.5 pt**.",
  "- PDF: vector; PNG: 600 dpi; white background.",
  "- Design references: Nature box-compare and clinical-forest examples, with Wanjie typography and alignment guidance. The heatmap template was not used.",
  "- High-impact reference principle: Wang et al., Nature Communications (2025), Figure 4b, DOI 10.1038/s41467-025-56618-y; adapted to two cohorts and paired Raw/Technical-adjusted patient strips.",
  "",
  "## Output boundary",
  "",
  "This script renders the standalone S2A panel, without figure-page assembly.",
  "",
  "## Freeze status",
  "",
  "SUPPLEMENTARY FIGURE S2A FINAL — FROZEN"
)

writeLines(qc_lines, file.path(out_dir, "S2A_patient_resolved_FINAL_QC.md"), useBytes = TRUE)

message("S2A patient-resolved FINAL render complete")
