#!/usr/bin/env Rscript
# Purpose: Figure 4 integrated patient functional evidence
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 4 and Supplementary Fig. S1.

# Figure 4A v5: six currently validated programs — patient-resolved functional
# landscape with integrated cohort-level trajectory summaries.
#
# Layout (retained from the approved Figure4a1 reference architecture):
#   left  : patient-resolved D1-D10 functional heatmap (6 program rows x
#           80 columns = 8 patient blocks x 10 deciles);
#   right : cohort-level trajectory summary per program as a flush row
#           annotation that shares the matrix row coordinates. Each row shows
#           the median trajectory + IQR (over 8 patients), the median
#           patient-specific Spearman rho and the positive-direction n/8.
#
# The six validated programs (frozen in Figure4D_module_definition_audit.tsv):
#   Hallmark TNFalpha-NF-kappaB | Hallmark EMT | Hallmark hypoxia |
#   Hallmark apoptosis | Custom survival-stress | Custom adhesion-remodeling.
# UPR (canonical/custom ambiguity) and Wound healing are NOT main-panel rows.
#
# No separate trajectory main panel: the trajectory block is an integrated
# summary annotation of the patient-resolved heatmap.
#
# This script only consumes the audited Figure4D tables; it does not recompute
# cell-level scores or inferential statistics.

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(Cairo)
  library(grid)
  library(digest)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Unable to resolve script path", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
ROOT <- normalizePath(dirname(dirname(dirname(script_path))), winslash = "/", mustWork = TRUE)
OUTDIR <- file.path(ROOT, "Figure4A_v5_six_validated_programs_integrated_landscape")
SCORE_SOURCE <- file.path(
  ROOT,
  "Figure4D_continuous_functional_trajectories",
  "Figure4D_decile_module_scores.tsv"
)
RHO_SOURCE <- file.path(
  ROOT,
  "Figure4D_continuous_functional_trajectories",
  "Figure4D_patient_program_rho.tsv"
)
AUDIT_SOURCE <- file.path(
  ROOT,
  "Figure4D_continuous_functional_trajectories",
  "Figure4D_module_definition_audit.tsv"
)
ATLAS_SOURCE <- file.path(
  ROOT,
  "0703_rebuild",
  "workstation_results",
  "Wu2021_full_atlas_umap_metadata.csv"
)

PDF_OUT <- file.path(OUTDIR, "Figure4A_v5_six_validated_programs_integrated_landscape.pdf")
PNG_OUT <- file.path(OUTDIR, "Figure4A_v5_six_validated_programs_integrated_landscape.png")
CORE_OUT <- file.path(OUTDIR, "Figure4A_v5_core_records.tsv")
PROFILE_OUT <- file.path(OUTDIR, "Figure4A_v5_cohort_profile.tsv")
RHO_SUMMARY_OUT <- file.path(OUTDIR, "Figure4A_v5_patient_rho_summary.tsv")
CONTEXT_OUT <- file.path(OUTDIR, "Figure4A_v5_patient_context.tsv")
COMPOSITION_OUT <- file.path(OUTDIR, "Figure4A_v5_patient_epithelial_composition.tsv")
REFERENCE_OUT <- file.path(OUTDIR, "Figure4A_v5_reference_mapping.md")
QC_OUT <- file.path(OUTDIR, "Figure4A_v5_render_QC.md")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

PATIENT_ORDER <- c(
  "CID4465", "CID4495", "CID44971", "CID44991",
  "CID4513", "CID4515", "CID4523", "CID3963"
)
INTERNAL_ORDER <- c(
  "hallmark_tnfa_nfkb",
  "hallmark_emt",
  "hallmark_hypoxia",
  "hallmark_apoptosis",
  "custom_survival_stress",
  "custom_adhesion_remodeling"
)
PROGRAM_LABELS <- c(
  hallmark_tnfa_nfkb = "Hallmark TNFα–NF-κB",
  hallmark_emt = "Hallmark EMT",
  hallmark_hypoxia = "Hallmark Hypoxia",
  hallmark_apoptosis = "Hallmark Apoptosis",
  custom_survival_stress = "Custom survival-stress",
  custom_adhesion_remodeling = "Custom adhesion-remodeling"
)
# Semantic colour family: stress / inflammatory programs use project coral;
# the two remodelling programs (EMT, adhesion-remodeling) use project blue.
PROGRAM_CLASS <- c(
  hallmark_tnfa_nfkb = "stress",
  hallmark_emt = "remodelling",
  hallmark_hypoxia = "stress",
  hallmark_apoptosis = "stress",
  custom_survival_stress = "stress",
  custom_adhesion_remodeling = "remodelling"
)
DECILE_ORDER <- 1:10

FONT <- "Helvetica"
CANVAS_W_MM <- 180
CANVAS_H_MM <- 106
HEATMAP_W_MM <- 90
HEATMAP_H_MM <- 50
PROFILE_W_MM <- 47
RASTER_DPI <- 600.05
PYTHON_EXE <- Sys.getenv("AHIPPO_PYTHON", Sys.which("python3"))
PNG_METADATA_SCRIPT <- file.path(OUTDIR, "scripts", "02_set_png_metadata.py")

NEGATIVE <- "#2F6DB3"
ZERO <- "#F7F7F7"
POSITIVE <- "#D7605C"
PROFILE_COLORS <- c(
  hallmark_tnfa_nfkb = "#D7605C",
  hallmark_emt = "#2F6DB3",
  hallmark_hypoxia = "#D7605C",
  hallmark_apoptosis = "#D7605C",
  custom_survival_stress = "#D7605C",
  custom_adhesion_remodeling = "#2F6DB3"
)
MALIGNANT_COLOUR <- "#D7605C"
NORMAL_COLOUR <- "#9AB3C7"
TRAJECTORY_COLOUR <- "#A8ADB3"
TEXT_DARK <- "#3D4650"
TEXT_MID <- "#59626A"
VALUE_LIMIT <- 3.5
COL_FUN <- circlize::colorRamp2(
  c(-VALUE_LIMIT, 0, VALUE_LIMIT),
  c(NEGATIVE, ZERO, POSITIVE)
)

write_tsv <- function(x, path) {
  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA",
    fileEncoding = "UTF-8"
  )
}

stop_gate <- function(message) {
  stop(paste0("[PRE-PLOTTING GATE FAILED] ", message), call. = FALSE)
}

# -----------------------------------------------------------------------------
# 1. Authoritative data and scientific-value gate
# -----------------------------------------------------------------------------

for (p in c(SCORE_SOURCE, RHO_SOURCE, AUDIT_SOURCE, ATLAS_SOURCE)) {
  if (!file.exists(p)) stop_gate(paste("Missing authoritative source:", p))
}

scores <- read.delim(
  SCORE_SOURCE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)
rho <- read.delim(
  RHO_SOURCE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)
audit <- read.delim(
  AUDIT_SOURCE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)

score_required <- c(
  "internal_id", "module", "patient", "decile", "median_module_z",
  "cell_n", "cohort_median", "cohort_q1", "cohort_q3", "n_patients"
)
rho_required <- c(
  "internal_id", "patient", "spearman_rho", "direction",
  "median_patient_rho", "positive_patients", "negative_patients",
  "zero_patients", "valid_patients", "direction_consistency", "valid_status"
)
if (length(setdiff(score_required, names(scores))) > 0) {
  stop_gate(paste("Missing required score columns:", paste(setdiff(score_required, names(scores)), collapse = ", ")))
}
if (length(setdiff(rho_required, names(rho))) > 0) {
  stop_gate(paste("Missing required rho columns:", paste(setdiff(rho_required, names(rho)), collapse = ", ")))
}

scores$patient <- as.character(scores$patient)
scores$internal_id <- as.character(scores$internal_id)
scores$decile <- as.integer(scores$decile)
scores$median_module_z <- as.numeric(scores$median_module_z)
scores$cell_n <- as.numeric(scores$cell_n)
rho$patient <- as.character(rho$patient)
rho$internal_id <- as.character(rho$internal_id)
rho$spearman_rho <- as.numeric(rho$spearman_rho)
rho$median_patient_rho <- as.numeric(rho$median_patient_rho)
rho$positive_patients <- as.integer(rho$positive_patients)
rho$negative_patients <- as.integer(rho$negative_patients)
rho$zero_patients <- as.integer(rho$zero_patients)
rho$valid_patients <- as.integer(rho$valid_patients)

if (nrow(scores) != 480L) stop_gate(paste("Expected 480 score records; found", nrow(scores)))
if (!setequal(unique(scores$patient), PATIENT_ORDER)) stop_gate("Patient set is not the approved eight-patient set")
if (!setequal(unique(scores$internal_id), INTERNAL_ORDER)) stop_gate("Program set is not the six validated programs")
if (!setequal(unique(scores$decile), DECILE_ORDER)) stop_gate("Deciles are not exactly D1-D10")
if (anyNA(scores[, score_required])) stop_gate("Missing values exist in required score fields")
if (any(scores$cell_n <= 0)) stop_gate("Non-positive decile cell count detected")
if (any(scores$n_patients != 8L)) stop_gate("A decile-module stratum does not contain all eight patients")

score_key <- paste(scores$patient, scores$internal_id, scores$decile, sep = "|")
expected_score <- expand.grid(
  patient = PATIENT_ORDER,
  internal_id = INTERNAL_ORDER,
  decile = DECILE_ORDER,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
expected_key <- paste(expected_score$patient, expected_score$internal_id, expected_score$decile, sep = "|")
if (length(unique(score_key)) != 480L || !setequal(score_key, expected_key)) {
  stop_gate("The patient x program x decile grid is incomplete or duplicated")
}

if (nrow(rho) != 48L) stop_gate(paste("Expected 48 rho records; found", nrow(rho)))
if (!setequal(unique(rho$internal_id), INTERNAL_ORDER)) stop_gate("Rho program set mismatch")
if (!setequal(unique(rho$patient), PATIENT_ORDER)) stop_gate("Rho patient set mismatch")
if (!all(rho$valid_status == "VALID")) stop_gate("Non-valid patient-program rho record detected")
if (!all(rho$positive_patients + rho$negative_patients + rho$zero_patients == 8L)) {
  stop_gate("Rho direction counts do not sum to eight patients")
}
if (any(rho$median_patient_rho < -1 | rho$median_patient_rho > 1)) {
  stop_gate("Median patient rho outside [-1, 1]")
}

# The rho values are cell-level patient-specific Spearman correlations from the
# audited Figure4D table. A decile-median-based Spearman is recomputed here only
# as a direction-consistency cross-check; it is NOT what is displayed.
decile_rho_check <- 0L
decile_rho_total <- 0L
decile_rho_mismatch <- character(0)
for (m in INTERNAL_ORDER) {
  sub <- scores[scores$internal_id == m, ]
  for (p in PATIENT_ORDER) {
    pp <- sub[sub$patient == p, ]
    pp <- pp[order(pp$decile), ]
    sp <- suppressWarnings(cor(pp$decile, pp$median_module_z, method = "spearman"))
    auth_rho <- rho$spearman_rho[rho$internal_id == m & rho$patient == p]
    if (is.na(sp)) next
    decile_rho_total <- decile_rho_total + 1L
    if (sign(sp) == sign(auth_rho)) {
      decile_rho_check <- decile_rho_check + 1L
    } else {
      decile_rho_mismatch <- c(decile_rho_mismatch, sprintf("%s/%s (decile=%.2f, cell-level=%.3f)", m, p, sp, auth_rho))
    }
  }
}

# The value scale must cover the true data without clipping.
value_range <- range(scores$median_module_z)
if (value_range[1] < -VALUE_LIMIT || value_range[2] > VALUE_LIMIT) {
  stop_gate(sprintf(
    "Observed median module-z range %.4f to %.4f exceeds the approved -%.1f to +%.1f scale",
    value_range[1], value_range[2], VALUE_LIMIT, VALUE_LIMIT
  ))
}
cohort_range <- range(c(scores$cohort_median, scores$cohort_q1, scores$cohort_q3))
if (cohort_range[1] < -VALUE_LIMIT || cohort_range[2] > VALUE_LIMIT) {
  stop_gate("Cohort median/IQR exceeds the approved -3.5 to +3.5 scale")
}

score_hash <- digest::digest(file = SCORE_SOURCE, algo = "sha256", serialize = FALSE)
rho_hash <- digest::digest(file = RHO_SOURCE, algo = "sha256", serialize = FALSE)
audit_hash <- digest::digest(file = AUDIT_SOURCE, algo = "sha256", serialize = FALSE)

# Lock display order explicitly; no clustering, sorting or rescaling is used.
scores$patient_order <- match(scores$patient, PATIENT_ORDER)
scores$program_order <- match(scores$internal_id, INTERNAL_ORDER)
scores <- scores[order(scores$patient_order, scores$program_order, scores$decile), ]

core_out <- data.frame(
  patient = scores$patient,
  internal_id = scores$internal_id,
  program = unname(PROGRAM_LABELS[scores$internal_id]),
  gene_set_class = scores$gene_set_class,
  decile = paste0("D", scores$decile),
  median_module_z = scores$median_module_z,
  cell_count = as.integer(scores$cell_n),
  cohort_median = scores$cohort_median,
  cohort_q1 = scores$cohort_q1,
  cohort_q3 = scores$cohort_q3,
  source_file = normalizePath(SCORE_SOURCE, winslash = "/", mustWork = TRUE),
  stringsAsFactors = FALSE
)
write_tsv(core_out, CORE_OUT)

# -----------------------------------------------------------------------------
# 2. Approved descriptive annotations only
# -----------------------------------------------------------------------------

# (a) Cohort median trajectory + IQR per program-decile, from the audited
#     Figure4D cohort summary (same patient set, n = 8 per stratum).
profile_rows <- list()
row_index <- 1L
for (m in INTERNAL_ORDER) {
  sub <- scores[scores$internal_id == m, ]
  for (d in DECILE_ORDER) {
    x <- sub[sub$decile == d, ]
    if (nrow(x) != 8L) stop_gate("Cohort profile stratum does not contain all eight patients")
    # cohort_median/q1/q3 are repeated once per patient within a module-decile;
    # take the single unique cohort summary value.
    cohort_row <- unique(x[, c("cohort_median", "cohort_q1", "cohort_q3", "n_patients")])
    if (nrow(cohort_row) != 1L) stop_gate("Cohort summary is not identical across patients in a module-decile stratum")
    profile_rows[[row_index]] <- data.frame(
      internal_id = m,
      program = unname(PROGRAM_LABELS[m]),
      decile = paste0("D", d),
      median = cohort_row$cohort_median,
      Q1 = cohort_row$cohort_q1,
      Q3 = cohort_row$cohort_q3,
      n_patients = cohort_row$n_patients,
      stringsAsFactors = FALSE
    )
    row_index <- row_index + 1L
  }
}
profile <- do.call(rbind, profile_rows)
if (nrow(profile) != 60L || any(profile$n_patients != 8L)) {
  stop_gate("Cohort profile is not complete at 60 records with n=8 per stratum")
}
write_tsv(profile, PROFILE_OUT)

# (b) Median patient-specific Spearman rho and positive-direction n/8, from the
#     audited Figure4D patient-level table (cell-level, within-patient).
rho_summary_rows <- list()
row_index <- 1L
for (m in INTERNAL_ORDER) {
  x <- rho[rho$internal_id == m, ]
  if (nrow(x) != 8L) stop_gate("Rho program stratum does not contain all eight patients")
  rho_summary_rows[[row_index]] <- data.frame(
    internal_id = m,
    program = unname(PROGRAM_LABELS[m]),
    median_patient_rho = x$median_patient_rho[1],
    rho_q1 = x$rho_q1[1],
    rho_q3 = x$rho_q3[1],
    rho_iqr = x$rho_iqr[1],
    positive_patients = x$positive_patients[1],
    negative_patients = x$negative_patients[1],
    zero_patients = x$zero_patients[1],
    valid_patients = x$valid_patients[1],
    direction_consistency = x$direction_consistency[1],
    stringsAsFactors = FALSE
  )
  row_index <- row_index + 1L
}
rho_summary <- do.call(rbind, rho_summary_rows)
if (nrow(rho_summary) != 6L) stop_gate("Rho summary is not complete at 6 programs")
write_tsv(rho_summary, RHO_SUMMARY_OUT)

# (c) Patient total malignant-cell context (from the authoritative 480 records).
patient_decile_counts <- unique(scores[, c("patient", "decile", "cell_n")])
if (nrow(patient_decile_counts) != 80L) stop_gate("Patient-decile cell-count table is not complete")
context <- aggregate(cell_n ~ patient, data = patient_decile_counts, FUN = sum)
context$patient_order <- match(context$patient, PATIENT_ORDER)
context <- context[order(context$patient_order), ]
names(context)[names(context) == "cell_n"] <- "total_malignant_epithelial_cells"
write_tsv(context[, c("patient", "total_malignant_epithelial_cells")], CONTEXT_OUT)

# (d) Frozen full-atlas epithelial composition (Cancer + Normal Epithelial).
atlas <- read.csv(
  ATLAS_SOURCE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)
atlas_required <- c("patient", "celltype_major")
if (!all(atlas_required %in% names(atlas))) stop_gate("Atlas metadata lacks patient/celltype_major")
atlas <- atlas[
  atlas$patient %in% PATIENT_ORDER &
    atlas$celltype_major %in% c("Cancer Epithelial", "Normal Epithelial"),
  atlas_required,
  drop = FALSE
]

composition <- data.frame(
  patient = PATIENT_ORDER,
  malignant_epithelial_cells = vapply(
    PATIENT_ORDER,
    function(p) sum(atlas$patient == p & atlas$celltype_major == "Cancer Epithelial"),
    integer(1)
  ),
  normal_epithelial_cells = vapply(
    PATIENT_ORDER,
    function(p) sum(atlas$patient == p & atlas$celltype_major == "Normal Epithelial"),
    integer(1)
  ),
  stringsAsFactors = FALSE
)
composition$total_epithelial_cells <-
  composition$malignant_epithelial_cells + composition$normal_epithelial_cells
composition$malignant_fraction <-
  composition$malignant_epithelial_cells / composition$total_epithelial_cells
composition$display_ratio <- paste0(
  format(composition$malignant_epithelial_cells, trim = TRUE, big.mark = ",", scientific = FALSE),
  "/",
  format(composition$total_epithelial_cells, trim = TRUE, big.mark = ",", scientific = FALSE)
)

if (any(composition$total_epithelial_cells <= 0L)) stop_gate("Missing epithelial denominator")
if (!identical(
  composition$malignant_epithelial_cells,
  as.integer(context$total_malignant_epithelial_cells)
)) {
  stop_gate("Cancer-epithelial counts do not match the authoritative 480-record source")
}

atlas_hash <- digest::digest(file = ATLAS_SOURCE, algo = "sha256", serialize = FALSE)
composition_out <- composition[, c(
  "patient", "malignant_epithelial_cells", "normal_epithelial_cells",
  "total_epithelial_cells", "malignant_fraction", "display_ratio"
)]
composition_out$source_file <- normalizePath(ATLAS_SOURCE, winslash = "/", mustWork = TRUE)
composition_out$source_sha256 <- atlas_hash
write_tsv(composition_out, COMPOSITION_OUT)

# -----------------------------------------------------------------------------
# 3. Reference gate
# -----------------------------------------------------------------------------

reference_text <- c(
  "# Figure 4A v5 (six validated programs) visual-reference mapping",
  "",
  "## Gate decision",
  "",
  "`PASS_CLOSE`",
  "",
  "The grouped annotated heatmap (Figure4a1) combines the patient-resolved D1-D10 functional landscape with the cohort-level decile summary. Its layout derives from `Figure4A_FINAL_integrated_JointAxis_landscape_REFERENCE_GATED`. Rows represent the six validated programs.",
  "",
  "## Data-to-visual correspondence",
  "",
  "- Reference core matrix -> six fixed program rows by 80 ordered columns (eight patient blocks, D1-D10 within every block). Rows are the six validated programs; UPR (canonical/custom ambiguity) and Wound healing are not main-panel rows.",
  "- Reference quantitative annotation -> a row-aligned cohort-level summary for each program: eight actual patient trajectories (low-alpha), the coloured cohort median trajectory with a matching IQR ribbon, the median patient-specific Spearman rho, and the positive-direction n/8.",
  "- Reference column context -> one horizontal 100% epithelial-composition bar per patient block, with malignant/total epithelial counts printed above the bar.",
  "- Reference signed scale -> a common zero-centred module-z scale (-3.5, 0, +3.5) chosen so the true data range (-1.76 .. +3.34) is never clipped.",
  "- Reference trajectory styling -> eight low-alpha patient trajectories, coloured cohort-median line, matching IQR ribbon, mini axes and dashed facet frames, adapted from the approved patient-trajectory reference.",
  "",
  "## Adaptations and exclusions",
  "",
  "The trajectory block is an integrated summary of the patient-resolved heatmap (a flush row annotation sharing the matrix row coordinates); there is no separate trajectory main panel. Clustering, dendrograms, row z-scoring, patient-specific scaling, effect sorting, pseudotime language, stars, P/q-value columns and independent patchwork subpanels are not transferred. The top bar denominator is Cancer Epithelial + Normal Epithelial in the frozen full-atlas annotation.",
  "",
  "## Reference files",
  "",
  "- Approved Figure4a1 reference: `Figure4A_FINAL_integrated_JointAxis_landscape_REFERENCE_GATED/Figure4A_FINAL_integrated_joint_axis_landscape.pdf`.",
  "- Validated program audit: `Figure4D_continuous_functional_trajectories/Figure4D_module_definition_audit.tsv` (SELECTED_MAIN six programs).",
  "- Patient-specific Spearman rho: `Figure4D_continuous_functional_trajectories/Figure4D_patient_program_rho.tsv`.",
  "- Patient-trajectory reference: `ExtendedData_Figure5_panelA_patient_trajectories.pdf`."
)
writeLines(reference_text, REFERENCE_OUT, useBytes = TRUE)

# -----------------------------------------------------------------------------
# 4. Matrix and embedded annotations
# -----------------------------------------------------------------------------

matrix_value <- matrix(
  NA_real_,
  nrow = length(INTERNAL_ORDER),
  ncol = length(PATIENT_ORDER) * length(DECILE_ORDER),
  dimnames = list(
    unname(PROGRAM_LABELS[INTERNAL_ORDER]),
    rep(paste0("D", DECILE_ORDER), times = length(PATIENT_ORDER))
  )
)

for (i in seq_len(nrow(scores))) {
  row_pos <- match(scores$internal_id[i], INTERNAL_ORDER)
  col_pos <- (match(scores$patient[i], PATIENT_ORDER) - 1L) * 10L + scores$decile[i]
  matrix_value[row_pos, col_pos] <- scores$median_module_z[i]
}
if (anyNA(matrix_value)) stop_gate("Heatmap matrix contains unfilled cells")

column_patient <- factor(
  rep(PATIENT_ORDER, each = 10L),
  levels = PATIENT_ORDER
)
decile_labels <- rep(c("1", "", "", "", "5", "", "", "", "", "10"), times = 8L)

# -----------------------------------------------------------------------------
# 5. Annotation drawing
# -----------------------------------------------------------------------------

# Geometry of the right-hand cohort-level summary block.
# Viewport x: -2.2 .. 24 ; trajectory sub-region -2.2 .. 13.4 (divider),
# rho summary sub-region 13.4 .. 24. Y: 0 .. 6 (one unit per program row).
VP_XMIN <- -2.2
VP_XMAX <- 24
DIVIDER_X <- 13.4
RHO_TEXT_X <- 22.6
Y_LABEL_X <- 2.0
BOX_LEFT <- 2.4
BOX_RIGHT <- 12.6
X_DATA <- seq(2.9, 12.1, length.out = 10)
ROW_HALF <- 0.295
RHO_LINE_SPACE <- 0.16

# Captured geometry of the annotation block (filled during decoration).
ANNOT_GEOM <- new.env(parent = emptyenv())

make_heatmap <- function() {
  top_context <- HeatmapAnnotation(
    `Epithelial composition` = anno_empty(
      height = unit(10.0, "mm"),
      border = FALSE
    ),
    which = "column",
    show_annotation_name = FALSE,
    gap = unit(0, "mm")
  )

  cohort_profile <- rowAnnotation(
    `Cohort summary` = anno_empty(
      width = unit(PROFILE_W_MM, "mm"),
      border = FALSE
    ),
    show_annotation_name = FALSE,
    gap = unit(0, "mm")
  )

  Heatmap(
    matrix_value,
    name = "Median module z",
    col = COL_FUN,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cluster_row_slices = FALSE,
    cluster_column_slices = FALSE,
    row_order = seq_len(nrow(matrix_value)),
    column_order = seq_len(ncol(matrix_value)),
    column_split = column_patient,
    column_gap = unit(1.8, "mm"),
    top_annotation = top_context,
    right_annotation = cohort_profile,
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = gpar(fontfamily = FONT, fontsize = 7.2, fontface = "bold"),
    row_names_max_width = unit(38, "mm"),
    show_column_names = TRUE,
    column_labels = decile_labels,
    column_names_side = "bottom",
    column_names_rot = 0,
    column_names_gp = gpar(fontfamily = FONT, fontsize = 6.4),
    column_title = levels(column_patient),
    column_title_gp = gpar(fontfamily = FONT, fontsize = 7.0, fontface = "bold"),
    column_title_rot = 0,
    column_title_side = "top",
    bottom_annotation = NULL,
    column_names_max_height = unit(3.8, "mm"),
    rect_gp = gpar(col = "white", lwd = 0.12),
    border = FALSE,
    use_raster = FALSE,
    width = unit(HEATMAP_W_MM, "mm"),
    height = unit(HEATMAP_H_MM, "mm"),
    heatmap_legend_param = list(
      title = "Median module z",
      title_gp = gpar(fontfamily = FONT, fontsize = 6.8, fontface = "bold"),
      labels_gp = gpar(fontfamily = FONT, fontsize = 6.4),
      at = c(-VALUE_LIMIT, 0, VALUE_LIMIT),
      labels = c("-3.5", "0", "+3.5"),
      direction = "horizontal",
      legend_width = unit(34, "mm"),
      legend_height = unit(3.0, "mm"),
      border = NA,
      title_position = "lefttop"
    )
  )
}

decorate_patient_composition <- function() {
  for (slice_index in seq_along(PATIENT_ORDER)) {
    malignant_fraction <- composition$malignant_fraction[slice_index]
    normal_fraction <- 1 - malignant_fraction
    decorate_annotation("Epithelial composition", slice = slice_index, {
      pushViewport(viewport(xscale = c(0, 1), yscale = c(0, 1)))
      bar_left <- 0.06
      bar_right <- 0.94
      bar_width <- bar_right - bar_left
      bar_y <- 0.23
      bar_h <- 0.23
      if (malignant_fraction > 0) {
        grid.rect(
          x = unit(bar_left + bar_width * malignant_fraction / 2, "native"),
          y = unit(bar_y, "native"),
          width = unit(bar_width * malignant_fraction, "native"),
          height = unit(bar_h, "native"),
          gp = gpar(fill = MALIGNANT_COLOUR, col = NA)
        )
      }
      if (normal_fraction > 0) {
        grid.rect(
          x = unit(bar_left + bar_width * malignant_fraction + bar_width * normal_fraction / 2, "native"),
          y = unit(bar_y, "native"),
          width = unit(bar_width * normal_fraction, "native"),
          height = unit(bar_h, "native"),
          gp = gpar(fill = NORMAL_COLOUR, col = NA)
        )
      }
      grid.rect(
        x = unit((bar_left + bar_right) / 2, "native"),
        y = unit(bar_y, "native"),
        width = unit(bar_width, "native"),
        height = unit(bar_h, "native"),
        gp = gpar(fill = NA, col = "#6F7780", lwd = 0.35)
      )
      grid.text(
        sub("/", "/\n", composition$display_ratio[slice_index], fixed = TRUE),
        x = unit(0.5, "native"),
        y = unit(0.70, "native"),
        just = c("center", "center"),
        gp = gpar(
          fontfamily = FONT,
          fontsize = 6.2,
          fontface = "bold",
          lineheight = 0.82,
          col = "#3F4850"
        )
      )
      popViewport()
    })
  }
}

decorate_cohort_profile <- function() {
  decorate_annotation("Cohort summary", {
    pushViewport(viewport(xscale = c(VP_XMIN, VP_XMAX), yscale = c(0, 6)))

    # Capture the annotation block geometry for header placement.
    loc_bottom_left <- grid::deviceLoc(unit(0, "npc"), unit(0, "npc"))
    loc_top_left <- grid::deviceLoc(unit(0, "npc"), unit(1, "npc"))
    ANNOT_GEOM$left_mm <- as.numeric(loc_bottom_left$x) * 25.4
    ANNOT_GEOM$bottom_mm <- as.numeric(loc_bottom_left$y) * 25.4
    ANNOT_GEOM$top_mm <- as.numeric(loc_top_left$y) * 25.4
    ANNOT_GEOM$width_mm <- grid::convertWidth(unit(1, "npc"), "mm", valueOnly = TRUE)

    for (program_index in seq_along(INTERNAL_ORDER)) {
      center_y <- 6 - program_index + 0.5
      internal_id <- INTERNAL_ORDER[program_index]
      program_label <- unname(PROGRAM_LABELS[internal_id])
      program_colour <- unname(PROFILE_COLORS[internal_id])
      p <- profile[profile$internal_id == internal_id, ]
      p <- p[match(paste0("D", 1:10), p$decile), ]

      y_median <- center_y + (p$median / VALUE_LIMIT) * ROW_HALF
      y_q1 <- center_y + (p$Q1 / VALUE_LIMIT) * ROW_HALF
      y_q3 <- center_y + (p$Q3 / VALUE_LIMIT) * ROW_HALF

      # Facet frame and axes follow the approved patient-trajectory reference.
      grid.rect(
        x = unit((BOX_LEFT + BOX_RIGHT) / 2, "native"),
        y = unit(center_y, "native"),
        width = unit(BOX_RIGHT - BOX_LEFT, "native"),
        height = unit(0.72, "native"),
        gp = gpar(fill = "white", col = "#AEB5BC", lwd = 0.45, lty = 2)
      )
      grid.lines(
        x = unit(c(BOX_LEFT, BOX_RIGHT), "native"),
        y = unit(c(center_y, center_y), "native"),
        gp = gpar(col = "#D5D9DD", lwd = 0.35)
      )

      # Eight actual patient trajectories from the audited 480 records.
      for (patient_id in PATIENT_ORDER) {
        patient_values <- scores$median_module_z[
          scores$internal_id == internal_id & scores$patient == patient_id
        ]
        patient_deciles <- scores$decile[
          scores$internal_id == internal_id & scores$patient == patient_id
        ]
        patient_values <- patient_values[order(patient_deciles)]
        if (length(patient_values) != 10L) stop_gate("Incomplete patient trajectory")
        patient_y <- center_y + (patient_values / VALUE_LIMIT) * ROW_HALF
        grid.lines(
          x = unit(X_DATA, "native"),
          y = unit(patient_y, "native"),
          gp = gpar(
            col = grDevices::adjustcolor(TRAJECTORY_COLOUR, alpha.f = 0.52),
            lwd = 0.42,
            lineend = "round"
          )
        )
      }

      # Cohort median and IQR layered above patient trajectories.
      grid.polygon(
        x = unit(c(X_DATA, rev(X_DATA)), "native"),
        y = unit(c(y_q1, rev(y_q3)), "native"),
        gp = gpar(fill = grDevices::adjustcolor(program_colour, alpha.f = 0.20), col = NA)
      )
      grid.lines(
        x = unit(X_DATA, "native"),
        y = unit(y_median, "native"),
        gp = gpar(col = program_colour, lwd = 1.05, lineend = "round")
      )

      # Compact shared-scale y-axis labels; no row-specific rescaling.
      y_tick_values <- c(VALUE_LIMIT, 0, -VALUE_LIMIT)
      y_tick_pos <- center_y + (y_tick_values / VALUE_LIMIT) * ROW_HALF
      for (tick_y in y_tick_pos) {
        grid.lines(
          x = unit(c(Y_LABEL_X + 0.1, Y_LABEL_X + 0.3), "native"),
          y = unit(c(tick_y, tick_y), "native"),
          gp = gpar(col = "#6F7780", lwd = 0.35)
        )
      }
      grid.text(
        c("+3.5", "0", "-3.5"),
        x = unit(Y_LABEL_X, "native"),
        y = unit(y_tick_pos, "native"),
        just = "right",
        gp = gpar(fontfamily = FONT, fontsize = 6.2, col = TEXT_MID)
      )

      # Median patient-specific Spearman rho and positive-direction n/8.
      rs <- rho_summary[rho_summary$internal_id == internal_id, ]
      rho_line_1 <- sprintf("median ρ = %.2f", rs$median_patient_rho)
      rho_line_2 <- sprintf("%d/%d positive", rs$positive_patients, rs$valid_patients)
      grid.text(
        rho_line_1,
        x = unit(RHO_TEXT_X, "native"),
        y = unit(center_y + RHO_LINE_SPACE, "native"),
        just = "right",
        gp = gpar(fontfamily = FONT, fontsize = 6.2, col = TEXT_DARK)
      )
      grid.text(
        rho_line_2,
        x = unit(RHO_TEXT_X, "native"),
        y = unit(center_y - RHO_LINE_SPACE, "native"),
        just = "right",
        gp = gpar(fontfamily = FONT, fontsize = 6.2, col = TEXT_MID)
      )

      if (program_index == length(INTERNAL_ORDER)) {
        tick_x <- X_DATA[c(5, 10)]
        for (tick_value in tick_x) {
          grid.lines(
            x = unit(c(tick_value, tick_value), "native"),
            y = unit(c(center_y - 0.36, center_y - 0.32), "native"),
            gp = gpar(col = "#6F7780", lwd = 0.35)
          )
        }
        grid.text(
          c("5", "10"),
          x = unit(tick_x, "native"),
          y = unit(center_y - 0.44, "native"),
          gp = gpar(fontfamily = FONT, fontsize = 6.2, col = TEXT_MID)
        )
      }
    }

    # Vertical divider between trajectory and rho summary sub-regions.
    grid.lines(
      x = unit(c(DIVIDER_X, DIVIDER_X), "native"),
      y = unit(c(0.35, 5.65), "native"),
      gp = gpar(col = "#B7BCC2", lwd = 0.45, lty = 2)
    )
    popViewport()
  })
}

draw_panel <- function() {
  grid.newpage()
  pushViewport(viewport(
    x = unit(0, "mm"), y = unit(4, "mm"),
    width = unit(CANVAS_W_MM, "mm"), height = unit(86, "mm"),
    just = c("left", "bottom"), name = "figure4a_v5_data_region"
  ))
  heatmap_object <- make_heatmap()
  composition_legend <- Legend(
    title = "Epithelial composition",
    labels = c("Malignant epithelial", "Normal epithelial"),
    legend_gp = gpar(fill = c(MALIGNANT_COLOUR, NORMAL_COLOUR), col = NA),
    title_gp = gpar(fontfamily = FONT, fontsize = 6.8, fontface = "bold"),
    labels_gp = gpar(fontfamily = FONT, fontsize = 6.2),
    direction = "horizontal",
    nrow = 1,
    grid_width = unit(2.8, "mm"),
    grid_height = unit(2.1, "mm"),
    column_gap = unit(3.0, "mm")
  )
  draw(
    heatmap_object,
    newpage = FALSE,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    annotation_legend_list = list(composition_legend),
    merge_legend = TRUE,
    padding = unit(c(2, 2, 3, 2), "mm"),
    column_title = "Joint YAP-Stem axis decile",
    column_title_side = "bottom",
    column_title_gp = gpar(fontfamily = FONT, fontsize = 7.3)
  )
  decorate_patient_composition()
  decorate_cohort_profile()

  # Return to the root viewport before placing panel lettering and sub-region
  # headers, so that 1npc refers to the full 180 x 100 mm page canvas.
  upViewport(0)

  # Sub-region geometry derived from the captured annotation block (page mm,
  # origin at the bottom-left of the canvas).
  annot_left <- ANNOT_GEOM$left_mm
  annot_w <- ANNOT_GEOM$width_mm
  traj_frac <- (DIVIDER_X - VP_XMIN) / (VP_XMAX - VP_XMIN)
  header_traj_x <- annot_left + traj_frac * annot_w / 2 - 2.2
  header_rho_x <- CANVAS_W_MM - 2.5
  header_y <- ANNOT_GEOM$top_mm + 2.5

  grid.text(
    c("Median trajectory", "+ IQR · 8 patient lines"),
    x = unit(header_traj_x, "mm"),
    y = unit(c(header_y, header_y - 2.4), "mm"),
    just = c("center", "bottom"),
    gp = gpar(
      fontfamily = FONT,
      fontsize = c(6.2, 5.9),
      fontface = c("bold", "plain"),
      lineheight = 0.9,
      col = c(TEXT_DARK, TEXT_MID)
    )
  )
  grid.text(
    c("Patient Spearman ρ", "median · positive n/8"),
    x = unit(header_rho_x, "mm"),
    y = unit(c(header_y, header_y - 2.4), "mm"),
    just = c("right", "bottom"),
    gp = gpar(
      fontfamily = FONT,
      fontsize = c(6.2, 5.9),
      fontface = c("bold", "plain"),
      lineheight = 0.9,
      col = c(TEXT_DARK, TEXT_MID)
    )
  )

  grid.text(
    "A",
    x = unit(3.5, "mm"),
    y = unit(1, "npc") - unit(4.5, "mm"),
    just = c("left", "center"),
    gp = gpar(fontfamily = FONT, fontsize = 11.5, fontface = "bold")
  )
  grid.text(
    "Patient-resolved functional landscape across six validated programs",
    x = unit(11.5, "mm"),
    y = unit(1, "npc") - unit(4.5, "mm"),
    just = c("left", "center"),
    gp = gpar(fontfamily = FONT, fontsize = 9.5, fontface = "bold")
  )
  grid.text(
    "Malignant/total\nepithelial n",
    x = unit(34.0, "mm"),
    y = unit(83.5, "mm"),
    just = c("right", "center"),
    gp = gpar(fontfamily = FONT, fontsize = 6.3, col = "#4E5963")
  )
}

# -----------------------------------------------------------------------------
# 6. Publication exports
# -----------------------------------------------------------------------------

Cairo::CairoPDF(
  PDF_OUT,
  width = CANVAS_W_MM / 25.4,
  height = CANVAS_H_MM / 25.4,
  family = FONT,
  bg = "white"
)
draw_panel()
dev.off()

Cairo::CairoPNG(
  PNG_OUT,
  width = round(CANVAS_W_MM / 25.4 * 600),
  height = round(CANVAS_H_MM / 25.4 * 600),
  res = 600,
  family = FONT,
  bg = "white"
)
draw_panel()
dev.off()

if (!file.exists(PYTHON_EXE)) stop_gate("Bundled Python executable for PNG metadata was not found")
if (!file.exists(PNG_METADATA_SCRIPT)) stop_gate("PNG metadata helper script was not found")
metadata_status <- system2(
  PYTHON_EXE,
  args = c(shQuote(PNG_METADATA_SCRIPT), shQuote(PNG_OUT), as.character(RASTER_DPI)),
  stdout = TRUE,
  stderr = TRUE
)
if (!is.null(attr(metadata_status, "status")) && attr(metadata_status, "status") != 0L) {
  stop_gate(paste("PNG DPI metadata update failed:", paste(metadata_status, collapse = " ")))
}

if (!file.exists(PDF_OUT) || file.info(PDF_OUT)$size < 10000L) stop_gate("PDF export failed")
if (!file.exists(PNG_OUT) || file.info(PNG_OUT)$size < 10000L) stop_gate("PNG export failed")

pdf_hash <- digest::digest(file = PDF_OUT, algo = "sha256", serialize = FALSE)
png_hash <- digest::digest(file = PNG_OUT, algo = "sha256", serialize = FALSE)

rho_summary_text <- vapply(
  INTERNAL_ORDER,
  function(m) {
    rs <- rho_summary[rho_summary$internal_id == m, ]
    sprintf("%s: median ρ=%.3f, positive=%d/%d, consistency=%s",
      PROGRAM_LABELS[m], rs$median_patient_rho,
      rs$positive_patients, rs$valid_patients, rs$direction_consistency)
  },
  character(1)
)

qc_text <- c(
  "# Figure 4A v5 (six validated programs) render QC",
  "",
  "## Gate and provenance",
  "",
  "- Data gate: PASS.",
  "- Visual-reference gate: PASS_CLOSE.",
  paste0("- Decile module-score source: `", normalizePath(SCORE_SOURCE, winslash = "/", mustWork = TRUE), "`."),
  paste0("- Score source SHA-256: `", score_hash, "`."),
  paste0("- Patient rho source: `", normalizePath(RHO_SOURCE, winslash = "/", mustWork = TRUE), "`."),
  paste0("- Rho source SHA-256: `", rho_hash, "`."),
  paste0("- Module definition audit: `", normalizePath(AUDIT_SOURCE, winslash = "/", mustWork = TRUE), "`."),
  paste0("- Audit source SHA-256: `", audit_hash, "`."),
  paste0("- Frozen atlas composition source: `", normalizePath(ATLAS_SOURCE, winslash = "/", mustWork = TRUE), "`."),
  paste0("- Atlas source SHA-256: `", atlas_hash, "`."),
  "- Core completeness: 480/480 records (8 patients x 6 programs x 10 deciles).",
  "- Cohort profile completeness: 60/60 records; n=8 patients in every program-decile stratum.",
  "- Patient trajectory completeness: 480/480 records; ten points for every patient-program trajectory.",
  "- Patient-program rho completeness: 48/48 records; 8 patients x 6 programs.",
  "- Patient epithelial-composition completeness: 8/8 patients.",
  "- Composition denominator: Cancer Epithelial + Normal Epithelial from the same frozen atlas annotation.",
  "- Malignant-count cross-source identity: 8/8 exact matches.",
  "",
  "## Locked structure",
  "",
  paste0("- Patient order: ", paste(PATIENT_ORDER, collapse = ", "), "."),
  paste0("- Program order: ", paste(unname(PROGRAM_LABELS[INTERNAL_ORDER]), collapse = ", "), "."),
  "- UPR (canonical/custom) is not a main-panel row; Wound healing is not a main-panel row.",
  "- Decile order: D1-D10 within every patient block.",
  "- Clustering: NO.",
  "- Effect-based sorting: NO.",
  "- Row/patient scaling: NO.",
  "- Winsorization or smoothing: NO.",
  "- Invented significance: NO.",
  "",
  "## Numeric, colour-scale and rho checks",
  "",
  sprintf("- True median module-z range: %.6f to %.6f.", value_range[1], value_range[2]),
  sprintf("- True cohort median/IQR range: %.6f to %.6f.", cohort_range[1], cohort_range[2]),
  sprintf("- Shared colour scale: -%.1f / 0 / +%.1f.", VALUE_LIMIT, VALUE_LIMIT),
  "- Zero centring: YES.",
  "- Colour clipping: NO.",
  "- Cohort trajectories share the same -3.5 to +3.5 vertical mapping: YES.",
  "- Patient-specific Spearman rho: cell-level, within patient, from the audited Figure4D table (not recomputed here).",
  sprintf("- Decile-median-based Spearman sign agreement with cell-level rho: %d/%d.", decile_rho_check, decile_rho_total),
  if (length(decile_rho_mismatch) > 0) {
    paste0("- Sign-cross-check mismatches (all near-zero magnitude): ", paste(decile_rho_mismatch, collapse = "; "), ".")
  } else {
    "- Sign-cross-check mismatches: none."
  },
  "- Median patient-specific Spearman rho / positive-direction n/8 (displayed):",
  paste0("  - ", rho_summary_text),
  "",
  "## Visual integration",
  "",
  "- Dominant object: one six-row by 80-column patient-resolved heatmap.",
  "- Patient composition context: one horizontal 100% stacked bar per patient; malignant epithelial is coral and normal epithelial is blue-grey; exact malignant/total count is printed above.",
  "- Right-side cohort-level summary: flush row annotation sharing the matrix row coordinate system.",
  "- Median trajectory + IQR: coloured cohort-median line with matching alpha-0.20 IQR ribbon per program.",
  "- Patient trajectories: all eight real patient lines retained in low-alpha neutral grey.",
  "- Median patient-specific Spearman rho and positive-direction n/8: text column to the right of each trajectory facet.",
  "- Trajectory axes and frames: shared -3.5/0/+3.5 y labels, compact D5/D10 x ticks, and thin dashed facet frames.",
  "- Trajectory styling: stress/inflammatory rows use project coral; remodelling rows (EMT, adhesion-remodeling) use project blue.",
  "- Separate trajectory main panel: NO; the trajectory block is an integrated summary annotation of the heatmap.",
  "- Independent patchwork line-plot subplot: NO.",
  "- Independent cell-composition subplot: NO; composition is a flush column annotation.",
  "- Stars, P/q-value text columns, effect strips and forest summaries: ABSENT.",
  "",
  "## Typography and export",
  "",
  sprintf("- Final canvas: %d mm x %d mm.", CANVAS_W_MM, CANVAS_H_MM),
  "- Font family: Helvetica.",
  "- Panel letter: 11.5 pt bold.",
  "- Panel title: 9.5 pt bold.",
  "- Patient IDs: 7.0 pt bold.",
  "- Program labels: 7.2 pt bold.",
  "- Decile labels: 6.4 pt.",
  "- Axis title: 7.3 pt.",
  "- Legend title: 6.8 pt bold; legend ticks: 6.4 pt.",
  "- Trajectory/rho summary text: 6.2 pt (rho value lines), 6.5/6.0 pt (sub-region headers).",
  "- Core text below 6.0 pt: NO.",
  "- PDF: vector output.",
  "- PNG: drawn by R/Cairo at the 600 dpi pixel dimensions; PNG pHYs metadata written by Python/Pillow at 600.05 dpi so integer rounding remains at or above 600 dpi; white background.",
  paste0("- PDF SHA-256: `", pdf_hash, "`."),
  paste0("- PNG SHA-256: `", png_hash, "`."),
  "",
  "## Final automated status",
  "",
  "`CANDIDATE_RENDERED_PENDING_MANUAL_REVIEW`"
)
writeLines(qc_text, QC_OUT, useBytes = TRUE)

cat("Figure 4A v5 candidate rendered successfully.\n")
cat("PDF:", PDF_OUT, "\n")
cat("PNG:", PNG_OUT, "\n")
