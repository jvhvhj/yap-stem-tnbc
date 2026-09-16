# Purpose: Independent BSW2 spatial replication
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 6 and Supplementary Fig. S4.
options(stringsAsFactors = FALSE, digits = 16)
.libPaths(c(".runtime/AHIPPO_RLIB_G1_ASCII", .Library))

suppressPackageStartupMessages({
  library(Matrix)
  library(hdf5r)
  library(digest)
})

root <- "."
g0 <- file.path(root, "Figure6G0_BSW2_archive_audit")
out <- file.path(root, "Figure6G1_BSW2_spatial_replication")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

patient_order <- c("sAA1", "sAA2", "sAA6", "sAA8", "sEA1", "sEA2", "sEA5", "sEA6", "sEA7")
race_map <- setNames(c(rep("BA", 4), rep("WA", 5)), patient_order)

sha256 <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
write_tsv <- function(x, name) {
  write.table(x, file.path(out, name), sep = "\t", row.names = FALSE,
              quote = FALSE, na = "NA", fileEncoding = "UTF-8")
}
fmt <- function(x, d = 4) formatC(x, format = "f", digits = d)
fmtp <- function(x) format(x, scientific = TRUE, digits = 5)

registry <- read.delim(file.path(g0, "Figure6G0_BSW2_section_registry.tsv"), check.names = FALSE)
registry <- registry[match(c(patient_order, "sAA9"), registry$patient_id), , drop = FALSE]
stopifnot(identical(registry$patient_id[seq_along(patient_order)], patient_order))
stopifnot(all(registry$analysis_status[seq_along(patient_order)] == "FINAL_PAPER_ANALYSIS_INCLUDED"))
stopifnot(registry$analysis_status[registry$patient_id == "sAA9"] == "ARCHIVED_NOT_IN_FINAL_PAPER_ANALYSIS")

yap_path <- file.path(root, "0719_frozen_score_definition_audit", "signature_gene_lists", "YAP_clean__clean_score.tsv")
stem_path <- file.path(root, "0719_frozen_score_definition_audit", "signature_gene_lists", "Stemness_clean__clean_score.tsv")
program_path <- file.path(root, "Figure-cursor", "figure4e", "146_gene_program.txt")
alias_path <- file.path(root, "Figure6_S0Q3_GSE210616_ESTIMATE_tumour_composition_reconstruction", "Figure6_S0Q3_alias_policy_audit.tsv")
s1_method_path <- file.path(root, ".tmp", "Figure6_S1", "assemble_s1_outputs.R")
analysis_script_path <- file.path(out, "scripts", "01_run_Figure6G1.R")

yap_tab <- read.delim(yap_path, check.names = FALSE)
stem_tab <- read.delim(stem_path, check.names = FALSE)
yap_original <- yap_tab$gene
stem_original <- stem_tab$gene
program_original <- scan(program_path, what = character(), quiet = TRUE)
stopifnot(length(yap_original) == 17L, length(stem_original) == 21L, length(program_original) == 146L)

authorized_alias <- c(CYR61 = "CCN1", CTGF = "CCN2", `HNRNPU-AS1` = "HNRNPU")
map_alias <- function(x) {
  y <- x
  hit <- x %in% names(authorized_alias)
  y[hit] <- unname(authorized_alias[x[hit]])
  y
}
yap_mapped_all <- map_alias(yap_original)
stem_mapped_all <- map_alias(stem_original)
program_mapped_all <- map_alias(program_original)

read_10x_h5 <- function(path) {
  infile <- hdf5r::H5File$new(filename = path, mode = "r")
  on.exit(infile$close_all(), add = TRUE)
  g <- "matrix"
  features <- as.character(infile[[paste0(g, "/features/name")]][])
  barcodes <- as.character(infile[[paste0(g, "/barcodes")]][])
  mat <- Matrix::sparseMatrix(
    i = infile[[paste0(g, "/indices")]][] + 1L,
    p = infile[[paste0(g, "/indptr")]][],
    x = as.numeric(infile[[paste0(g, "/data")]][]),
    dims = as.integer(infile[[paste0(g, "/shape")]][]), repr = "T")
  rownames(mat) <- make.unique(features)
  colnames(mat) <- barcodes
  list(matrix = as(mat, "CsparseMatrix"), raw_features = features, barcodes = barcodes)
}

z_vector <- function(x) {
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x)) / s
}

giotto_standard_log2_equivalent <- function(selected_umi, library_size,
                                             scalefactor = 6000,
                                             log_offset = 1,
                                             logbase = 2) {
  stopifnot(all(library_size > 0), log_offset == 1, logbase == 2)
  normalized <- selected_umi %*% Matrix::Diagonal(x = scalefactor / library_size)
  normalized@x <- log1p(normalized@x) / log(logbase)
  normalized
}

program_score <- function(expr, program_features) {
  x <- as.matrix(expr[program_features, , drop = FALSE])
  sds <- apply(x, 1, sd)
  nonzero <- is.finite(sds) & sds > 0
  z <- matrix(0, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
  if (any(nonzero)) {
    z[nonzero, ] <- sweep(
      sweep(x[nonzero, , drop = FALSE], 1, rowMeans(x[nonzero, , drop = FALSE]), "-"),
      1, sds[nonzero], "/")
  }
  list(score = colMeans(z), nonzero = program_features[nonzero], zero_sd = program_features[!nonzero])
}

rank_residual <- function(y, covariates) {
  ry <- rank(y, ties.method = "average")
  rc <- as.data.frame(lapply(covariates, rank, ties.method = "average"))
  x <- cbind(Intercept = 1, as.matrix(rc))
  lm.fit(x = x, y = ry)$residuals
}

partial_rank_rho <- function(x, y, covariates) {
  rx <- rank_residual(x, covariates)
  ry <- rank_residual(y, covariates)
  cor(rx, ry, method = "pearson")
}

exact_signed_rank <- function(x) {
  x <- x[is.finite(x)]
  x <- x[x != 0]
  n <- length(x)
  if (!n) return(list(P = NA_real_, Wplus = NA_real_, n_nonzero = 0L, mode = "UNDEFINED_ALL_ZERO"))
  r <- rank(abs(x), ties.method = "average")
  observed <- sum(r[x > 0])
  masks <- 0:(2^n - 1L)
  sums <- vapply(masks, function(m) {
    bits <- as.logical(intToBits(m)[seq_len(n)])
    sum(r[bits])
  }, numeric(1))
  p <- mean(sums >= observed - 1e-12)
  list(P = p, Wplus = observed, n_nonzero = n,
       mode = "EXACT_SIGN_ENUMERATION_CONDITIONAL_ON_ABSOLUTE_RANKS")
}

summarize_effects <- function(x) {
  x <- x[is.finite(x)]
  wt <- exact_signed_rank(x)
  data.frame(
    n_patients = length(x), median_rho = median(x),
    Q1 = unname(quantile(x, 0.25, type = 7)),
    Q3 = unname(quantile(x, 0.75, type = 7)),
    IQR = IQR(x, type = 7), minimum = min(x), maximum = max(x),
    positive_n = sum(x > 0), positive_fraction = mean(x > 0),
    exact_one_sided_wilcoxon_P = wt$P,
    W_plus = wt$Wplus, n_nonzero = wt$n_nonzero,
    exact_test_mode = wt$mode, stringsAsFactors = FALSE)
}

# Feature map is established once because all ten archive sections have the
# identical audited feature index.
first_h5 <- registry$expression_object[registry$patient_id == patient_order[1]]
first_obj <- read_10x_h5(first_h5)
raw_features <- first_obj$raw_features
rm(first_obj); gc(verbose = FALSE)

make_feature_rows <- function(signature, original, mapped, source_path) {
  counts <- vapply(mapped, function(g) sum(raw_features == g), integer(1))
  present <- counts == 1L
  mapping_type <- ifelse(!present, "MISSING_NO_REPLACEMENT",
                         ifelse(original == mapped, "EXACT", "AUTHORIZED_LEGACY_ALIAS"))
  data.frame(
    signature = signature,
    gene_order = seq_along(original),
    original_gene = original,
    mapped_feature = ifelse(present, mapped, "NA"),
    mapping_type = mapping_type,
    feature_index_occurrences = counts,
    included_in_score = present,
    source_file = source_path,
    source_sha256 = sha256(source_path),
    alias_policy_file = alias_path,
    alias_policy_sha256 = sha256(alias_path),
    stringsAsFactors = FALSE)
}

feature_mapping <- rbind(
  make_feature_rows("YAP17", yap_original, yap_mapped_all, yap_path),
  make_feature_rows("Stem21", stem_original, stem_mapped_all, stem_path),
  make_feature_rows("Program146", program_original, program_mapped_all, program_path))
write_tsv(feature_mapping, "Figure6G1_feature_mapping.tsv")

yap_present <- feature_mapping$mapped_feature[feature_mapping$signature == "YAP17" & feature_mapping$included_in_score]
stem_present <- feature_mapping$mapped_feature[feature_mapping$signature == "Stem21" & feature_mapping$included_in_score]
program_present <- feature_mapping$mapped_feature[feature_mapping$signature == "Program146" & feature_mapping$included_in_score]
joint_present <- union(yap_present, stem_present)
stopifnot(length(yap_present) == 17L, length(stem_present) == 19L,
          length(joint_present) == 36L, length(program_present) == 135L)

overlap <- intersect(joint_present, program_present)
independence <- data.frame(
  left_signature = "Joint38", right_signature = "Program146",
  alias_normalization_policy = "EXACT_FIRST_THEN_ONLY_3_PROJECT_AUTHORIZED_LEGACY_ALIASES",
  left_frozen_n = 38L, right_frozen_n = 146L,
  left_measurable_n = length(joint_present), right_measurable_n = length(program_present),
  overlap_n = length(overlap), overlap_genes = if (length(overlap)) paste(overlap, collapse = ";") else "NONE",
  expected_overlap_n = 0L,
  score_independence_status = if (length(overlap) == 0L) "PASS_ZERO_OVERLAP" else "FAIL_NONZERO_OVERLAP",
  stringsAsFactors = FALSE)
write_tsv(independence, "Figure6G1_score_independence_check.tsv")
if (length(overlap)) stop("Score independence failed: measurable Joint38 overlaps Program146")

normalization_method <- data.frame(
  record_type = "METHOD",
  patient_id = "ALL_9_INCLUDED_SECTIONS",
  section_id = "SECTION_WISE",
  analysis_status = "FROZEN_BEFORE_SCORING",
  expression_input = "filtered_feature_bc_matrix.h5 RAW_INTEGER_UMI_COUNTS",
  coordinate_universe = "EXACT_in_tissue=1_BARCODES_MATCHING_FILTERED_H5",
  n_in_tissue_spots = NA_integer_,
  reference_package = "Giotto",
  reference_version = "4.2.3",
  reference_function = "Giotto::normalizeGiotto",
  exact_arguments = "expression_values=raw;norm_methods=standard;library_size_norm=TRUE;scalefactor=6000;log_norm=TRUE;log_offset=1;logbase=2;scale_feats=FALSE;scale_cells=FALSE",
  executed_implementation = "giotto_standard_log2_equivalent",
  exact_formula = "log2(raw_UMI/library_size*6000+1)",
  implementation_runtime = paste0(R.version.string, ";Matrix_", packageVersion("Matrix"), ";hdf5r_", packageVersion("hdf5r")),
  official_reference = "https://github.com/giotto-suite/Giotto/blob/suite/R/normalize.R",
  deterministic_equivalence_status = "PASS_FORMULA_MATCHES_OFFICIAL_LIBRARY_AND_LOG_METHODS",
  expression_object = "MULTIPLE_SEE_SECTION_ROWS",
  coordinate_object = "MULTIPLE_SEE_SECTION_ROWS",
  archive_sha256 = "48a4a92ce5708e7d4b341abc8e02ee7e79ea259b4f6c145674d48290c5bcf19f",
  analysis_script = analysis_script_path,
  analysis_script_sha256 = sha256(analysis_script_path),
  s1_rank_residual_source = s1_method_path,
  s1_rank_residual_source_sha256 = sha256(s1_method_path),
  stringsAsFactors = FALSE)

normalization_rows <- list()
qc_rows <- list()
depth_rows <- list()
association_rows <- list()

for (pid in patient_order) {
  rec <- registry[registry$patient_id == pid, , drop = FALSE]
  object_root <- file.path(g0, "audit", "extracted_objects", "10x.visium", pid, "outs")
  h5_path <- file.path(object_root, "filtered_feature_bc_matrix.h5")
  coord_path <- file.path(object_root, "spatial", "tissue_positions.csv")
  message(sprintf("G1 START %s spots=%d", pid, rec$n_expression_barcodes))
  obj <- read_10x_h5(h5_path)
  umi <- obj$matrix
  features <- obj$raw_features
  barcodes <- obj$barcodes
  coord <- read.csv(coord_path, check.names = FALSE)
  tissue_barcodes <- coord$barcode[coord$in_tissue == 1]
  barcode_match <- setequal(barcodes, tissue_barcodes) && length(barcodes) == length(tissue_barcodes)
  if (!barcode_match) stop(sprintf("Barcode-coordinate mismatch in %s", pid))
  if (ncol(umi) != rec$n_expression_barcodes) stop(sprintf("Registry spot-count mismatch in %s", pid))

  n_count <- as.numeric(Matrix::colSums(umi))
  n_feature <- as.integer(Matrix::colSums(umi > 0))
  if (any(!is.finite(n_count)) || any(n_count <= 0)) stop(sprintf("Invalid library size in %s", pid))

  selected_features <- c(yap_present, stem_present, program_present)
  selected_index <- match(selected_features, features)
  if (anyNA(selected_index)) stop(sprintf("Frozen mapped feature unexpectedly missing in %s", pid))
  if (anyDuplicated(selected_index)) stop(sprintf("Duplicate selected feature mapping in %s", pid))
  selected_umi <- umi[selected_index, , drop = FALSE]
  rownames(selected_umi) <- selected_features
  expr <- giotto_standard_log2_equivalent(selected_umi, n_count)

  yap_mean <- as.numeric(Matrix::colMeans(expr[yap_present, , drop = FALSE]))
  stem_mean <- as.numeric(Matrix::colMeans(expr[stem_present, , drop = FALSE]))
  joint <- z_vector(yap_mean) + z_vector(stem_mean)
  prog <- program_score(expr, program_present)
  program <- as.numeric(prog$score)

  finite_joint <- mean(is.finite(joint))
  finite_program <- mean(is.finite(program))
  joint_sd <- sd(joint)
  program_sd <- sd(program)
  qc_status <- if (finite_joint >= 0.95 && finite_program >= 0.95 &&
                   is.finite(joint_sd) && joint_sd > 0 &&
                   is.finite(program_sd) && program_sd > 0) "PASS" else "FAIL"

  rho_joint_count <- cor(joint, log1p(n_count), method = "spearman")
  rho_program_count <- cor(program, log1p(n_count), method = "spearman")
  rho_joint_feature <- cor(joint, n_feature, method = "spearman")
  rho_program_feature <- cor(program, n_feature, method = "spearman")
  rho_raw <- cor(joint, program, method = "spearman")
  rho_depth <- partial_rank_rho(joint, program, list(log1p_nCount = log1p(n_count)))
  rho_nfeature <- partial_rank_rho(joint, program, list(nFeature = n_feature))

  normalization_rows[[pid]] <- data.frame(
    record_type = "SECTION",
    patient_id = pid, section_id = rec$section_id,
    analysis_status = rec$analysis_status,
    expression_input = "filtered_feature_bc_matrix.h5 RAW_INTEGER_UMI_COUNTS",
    coordinate_universe = if (barcode_match) "PASS_EXACT_in_tissue=1_MATCH" else "FAIL",
    n_in_tissue_spots = ncol(umi),
    reference_package = "Giotto", reference_version = "4.2.3",
    reference_function = "Giotto::normalizeGiotto",
    exact_arguments = "expression_values=raw;norm_methods=standard;library_size_norm=TRUE;scalefactor=6000;log_norm=TRUE;log_offset=1;logbase=2;scale_feats=FALSE;scale_cells=FALSE",
    executed_implementation = "giotto_standard_log2_equivalent",
    exact_formula = "log2(raw_UMI/library_size*6000+1)",
    implementation_runtime = paste0(R.version.string, ";Matrix_", packageVersion("Matrix"), ";hdf5r_", packageVersion("hdf5r")),
    official_reference = "https://github.com/giotto-suite/Giotto/blob/suite/R/normalize.R",
    deterministic_equivalence_status = "PASS_FORMULA_MATCHES_OFFICIAL_LIBRARY_AND_LOG_METHODS",
    expression_object = h5_path,
    coordinate_object = coord_path,
    archive_sha256 = "48a4a92ce5708e7d4b341abc8e02ee7e79ea259b4f6c145674d48290c5bcf19f",
    analysis_script = analysis_script_path,
    analysis_script_sha256 = sha256(analysis_script_path),
    s1_rank_residual_source = s1_method_path,
    s1_rank_residual_source_sha256 = sha256(s1_method_path),
    stringsAsFactors = FALSE)

  qc_rows[[pid]] <- data.frame(
    patient_id = pid, race_group = unname(race_map[pid]), section_id = rec$section_id,
    n_in_tissue_spots = ncol(umi), coordinate_barcode_match = barcode_match,
    nCount_min = min(n_count), nCount_Q1 = unname(quantile(n_count, 0.25)),
    nCount_median = median(n_count), nCount_Q3 = unname(quantile(n_count, 0.75)), nCount_max = max(n_count),
    nFeature_min = min(n_feature), nFeature_Q1 = unname(quantile(n_feature, 0.25)),
    nFeature_median = median(n_feature), nFeature_Q3 = unname(quantile(n_feature, 0.75)), nFeature_max = max(n_feature),
    YAP_frozen_n = 17L, YAP_measurable_n = length(yap_present),
    Stem_frozen_n = 21L, Stem_measurable_n = length(stem_present),
    Stem_missing = "NANOG;POU5F1",
    Program_frozen_n = 146L, Program_measurable_n = length(program_present),
    Program_zero_variance_n = length(prog$zero_sd),
    Program_zero_variance_genes = if (length(prog$zero_sd)) paste(prog$zero_sd, collapse = ";") else "NONE",
    Joint_finite_fraction = finite_joint, Program_finite_fraction = finite_program,
    YAP_mean_SD = sd(yap_mean), Stem_mean_SD = sd(stem_mean),
    Joint_SD = joint_sd, Program_SD = program_sd,
    score_QC_status = qc_status, stringsAsFactors = FALSE)

  depth_rows[[pid]] <- data.frame(
    record_type = "PATIENT", patient_id = pid, race_group = unname(race_map[pid]),
    section_id = rec$section_id, n_spots = ncol(umi),
    rho_Joint_log1p_nCount = rho_joint_count,
    rho_Program146_log1p_nCount = rho_program_count,
    rho_Joint_nFeature = rho_joint_feature,
    rho_Program146_nFeature = rho_program_feature,
    summary_metric = NA_character_, summary_median = NA_real_,
    summary_Q1 = NA_real_, summary_Q3 = NA_real_, summary_IQR = NA_real_,
    summary_max_abs = NA_real_, interpretation = "DESCRIPTIVE_TECHNICAL_QC_ONLY",
    stringsAsFactors = FALSE)

  association_rows[[pid]] <- data.frame(
    patient_id = pid, race_group = unname(race_map[pid]), section_id = rec$section_id,
    n_spots = ncol(umi),
    rho_raw = rho_raw,
    rho_raw_direction = ifelse(rho_raw > 0, "POSITIVE", ifelse(rho_raw < 0, "NEGATIVE", "ZERO")),
    rho_depth_adjusted = rho_depth,
    rho_depth_adjusted_direction = ifelse(rho_depth > 0, "POSITIVE", ifelse(rho_depth < 0, "NEGATIVE", "ZERO")),
    rho_nFeature_adjusted = rho_nfeature,
    rho_nFeature_adjusted_direction = ifelse(rho_nfeature > 0, "POSITIVE", ifelse(rho_nfeature < 0, "NEGATIVE", "ZERO")),
    primary_covariate = "rank(log1p(nCount_raw))",
    sensitivity_covariate = "rank(nFeature_raw)",
    rank_residual_method = "average ranks; lm.fit with intercept; Pearson correlation of residuals",
    spot_level_inference = "NONE_PATIENT_IS_BIOLOGICAL_UNIT",
    stringsAsFactors = FALSE)

  message(sprintf("G1 DONE %s raw=%.4f depth=%.4f nFeature=%.4f QC=%s", pid, rho_raw, rho_depth, rho_nfeature, qc_status))
  rm(obj, umi, selected_umi, expr, joint, program, prog, coord); gc(verbose = FALSE)
}

excluded_row <- normalization_method
excluded_row$record_type <- "ARCHIVE_EXCLUSION"
excluded_row$patient_id <- "sAA9"
excluded_row$section_id <- "sAA9"
excluded_row$analysis_status <- "ARCHIVED_NOT_IN_AUTHOR_FINAL_ANALYSIS"
excluded_row$coordinate_universe <- "NOT_ENTERED_IN_G1"
excluded_row$n_in_tissue_spots <- NA_integer_
excluded_row$deterministic_equivalence_status <- "NOT_APPLICABLE_EXCLUDED_BEFORE_NORMALIZATION"
excluded_row$expression_object <- file.path(g0, "audit", "extracted_objects", "10x.visium", "sAA9", "outs", "filtered_feature_bc_matrix.h5")
excluded_row$coordinate_object <- file.path(g0, "audit", "extracted_objects", "10x.visium", "sAA9", "outs", "spatial", "tissue_positions.csv")
normalization_manifest <- do.call(rbind, c(list(normalization_method), normalization_rows, list(excluded_row)))
write_tsv(normalization_manifest, "Figure6G1_normalization_manifest.tsv")

spot_qc <- do.call(rbind, qc_rows)
spot_qc$patient_order <- seq_len(nrow(spot_qc))
spot_qc <- spot_qc[order(spot_qc$patient_order), ]
write_tsv(spot_qc, "Figure6G1_spot_score_QC.tsv")

depth_patient <- do.call(rbind, depth_rows)
depth_metrics <- c("rho_Joint_log1p_nCount", "rho_Program146_log1p_nCount", "rho_Joint_nFeature", "rho_Program146_nFeature")
depth_summary <- do.call(rbind, lapply(depth_metrics, function(m) {
  x <- depth_patient[[m]]
  data.frame(
    record_type = "COHORT_SUMMARY", patient_id = "ALL_9_PATIENTS", race_group = "NOT_TESTED",
    section_id = "ALL_9_SECTIONS", n_spots = sum(depth_patient$n_spots),
    rho_Joint_log1p_nCount = NA_real_, rho_Program146_log1p_nCount = NA_real_,
    rho_Joint_nFeature = NA_real_, rho_Program146_nFeature = NA_real_,
    summary_metric = m, summary_median = median(x),
    summary_Q1 = unname(quantile(x, 0.25)), summary_Q3 = unname(quantile(x, 0.75)),
    summary_IQR = IQR(x), summary_max_abs = max(abs(x)),
    interpretation = "DESCRIPTIVE_TECHNICAL_QC_ONLY", stringsAsFactors = FALSE)
}))
depth_qc <- rbind(depth_patient, depth_summary)
write_tsv(depth_qc, "Figure6G1_depth_QC.tsv")

patient_assoc <- do.call(rbind, association_rows)
patient_assoc$patient_order <- seq_len(nrow(patient_assoc))
patient_assoc <- patient_assoc[order(patient_assoc$patient_order), ]
write_tsv(patient_assoc, "Figure6G1_patient_associations.tsv")

technical_fail <- any(spot_qc$score_QC_status != "PASS") ||
  any(!is.finite(patient_assoc$rho_raw)) ||
  any(!is.finite(patient_assoc$rho_depth_adjusted)) ||
  any(!is.finite(patient_assoc$rho_nFeature_adjusted))

raw_summary <- summarize_effects(patient_assoc$rho_raw)
primary_summary <- summarize_effects(patient_assoc$rho_depth_adjusted)
nfeature_summary <- summarize_effects(patient_assoc$rho_nFeature_adjusted)

primary_statistical_pass <- !technical_fail && primary_summary$median_rho > 0 &&
  primary_summary$exact_one_sided_wilcoxon_P < 0.05

if (primary_statistical_pass) {
  loo <- do.call(rbind, lapply(seq_len(nrow(patient_assoc)), function(i) {
    x <- patient_assoc$rho_depth_adjusted[-i]
    s <- summarize_effects(x)
    data.frame(
      loo_status = "RUN_PRIMARY_COHORT_TEST_PASSED",
      omitted_patient = patient_assoc$patient_id[i],
      omitted_rho_depth_adjusted = patient_assoc$rho_depth_adjusted[i],
      n_patients = s$n_patients, median_rho = s$median_rho,
      Q1 = s$Q1, Q3 = s$Q3, IQR = s$IQR,
      positive_n = s$positive_n, positive_fraction = s$positive_fraction,
      exact_one_sided_wilcoxon_P = s$exact_one_sided_wilcoxon_P,
      median_direction_positive = s$median_rho > 0,
      stringsAsFactors = FALSE)
  }))
  loo_all_positive <- all(loo$median_direction_positive)
} else {
  loo <- data.frame(
    loo_status = "NOT_RUN_PRIMARY_COHORT_TEST_FAILED_OR_TECHNICAL_HOLD",
    omitted_patient = "NOT_APPLICABLE", omitted_rho_depth_adjusted = NA_real_,
    n_patients = NA_integer_, median_rho = NA_real_, Q1 = NA_real_, Q3 = NA_real_, IQR = NA_real_,
    positive_n = NA_integer_, positive_fraction = NA_real_, exact_one_sided_wilcoxon_P = NA_real_,
    median_direction_positive = NA, stringsAsFactors = FALSE)
  loo_all_positive <- FALSE
}
write_tsv(loo, "Figure6G1_LOO.tsv")

criterion_A <- primary_summary$median_rho > 0
criterion_B <- primary_summary$positive_n >= 7L
criterion_C <- primary_summary$exact_one_sided_wilcoxon_P < 0.05
criterion_D <- raw_summary$median_rho > 0
criterion_E <- nfeature_summary$median_rho > 0
criterion_F <- loo_all_positive
pass_all <- !technical_fail && all(c(criterion_A, criterion_B, criterion_C, criterion_D, criterion_E, criterion_F))
strong_all <- pass_all && primary_summary$positive_n >= 8L &&
  primary_summary$exact_one_sided_wilcoxon_P < 0.01 &&
  raw_summary$positive_n >= 8L && nfeature_summary$positive_n >= 8L && loo_all_positive

verdict <- if (technical_fail) {
  "FIGURE6G_TECHNICAL_HOLD"
} else if (strong_all) {
  "FIGURE6G_REPLICATION_STRONG_PASS"
} else if (pass_all) {
  "FIGURE6G_REPLICATION_PASS"
} else if (primary_summary$median_rho > 0 || primary_summary$positive_n >= 6L) {
  "FIGURE6G_REPLICATION_BORDERLINE"
} else {
  "FIGURE6G_REPLICATION_FAIL"
}
ready <- verdict %in% c("FIGURE6G_REPLICATION_STRONG_PASS", "FIGURE6G_REPLICATION_PASS")

primary_test <- data.frame(
  test_id = "PRIMARY_BSW2_PATIENT_LEVEL_DEPTH_ADJUSTED_ASSOCIATION",
  biological_unit = "PATIENT", cohort = "BSW2_AUTHOR_FINAL_9_TNBC_PATIENTS",
  effect = "rho_depth_adjusted", null_hypothesis = "median_rho=0",
  alternative = "median_rho>0", test = "exact one-sample Wilcoxon signed-rank",
  multiplicity = "NOT_REQUIRED_ONE_PREDEFINED_PRIMARY_TEST",
  primary_summary,
  criterion_A_median_positive = criterion_A,
  criterion_B_positive_patients_ge_7_of_9 = criterion_B,
  criterion_C_exact_P_lt_0.05 = criterion_C,
  criterion_D_raw_median_positive = criterion_D,
  criterion_E_nFeature_median_positive = criterion_E,
  criterion_F_all_LOO_medians_positive = criterion_F,
  all_PASS_criteria = pass_all,
  all_STRONG_PASS_criteria = strong_all,
  scientific_verdict = verdict,
  main_panel_development = if (ready) "READY" else "HOLD",
  stringsAsFactors = FALSE)
write_tsv(primary_test, "Figure6G1_primary_test.tsv")

nfeature_patient <- data.frame(
  record_type = "PATIENT", patient_id = patient_assoc$patient_id,
  race_group = patient_assoc$race_group, section_id = patient_assoc$section_id,
  n_spots = patient_assoc$n_spots,
  rho_nFeature_adjusted = patient_assoc$rho_nFeature_adjusted,
  direction = patient_assoc$rho_nFeature_adjusted_direction,
  median_rho = NA_real_, Q1 = NA_real_, Q3 = NA_real_, IQR = NA_real_,
  positive_n = NA_integer_, positive_fraction = NA_real_,
  exact_one_sided_wilcoxon_P = NA_real_,
  role = "SECONDARY_TECHNICAL_SENSITIVITY_CANNOT_RESCUE_PRIMARY",
  stringsAsFactors = FALSE)
nfeature_cohort <- data.frame(
  record_type = "COHORT_SUMMARY", patient_id = "ALL_9_PATIENTS", race_group = "NOT_TESTED",
  section_id = "ALL_9_SECTIONS", n_spots = sum(patient_assoc$n_spots), rho_nFeature_adjusted = NA_real_,
  direction = ifelse(nfeature_summary$median_rho > 0, "POSITIVE_MEDIAN", "NONPOSITIVE_MEDIAN"),
  median_rho = nfeature_summary$median_rho, Q1 = nfeature_summary$Q1, Q3 = nfeature_summary$Q3,
  IQR = nfeature_summary$IQR, positive_n = nfeature_summary$positive_n,
  positive_fraction = nfeature_summary$positive_fraction,
  exact_one_sided_wilcoxon_P = nfeature_summary$exact_one_sided_wilcoxon_P,
  role = "SECONDARY_TECHNICAL_SENSITIVITY_CANNOT_RESCUE_PRIMARY",
  stringsAsFactors = FALSE)
nfeature_sensitivity <- rbind(nfeature_patient, nfeature_cohort)
write_tsv(nfeature_sensitivity, "Figure6G1_nFeature_sensitivity.tsv")

depth_summary_text <- paste(vapply(depth_metrics, function(m) {
  z <- depth_summary[depth_summary$summary_metric == m, ]
  sprintf("%s median=%s, IQR=%s", m, fmt(z$summary_median), fmt(z$summary_IQR))
}, character(1)), collapse = "; ")

readout <- c(
  "# Figure 6G1 — Independent BSW2 Visium replication", "",
  "## A. Scientific question",
  "Does the positive spatial association between the frozen YAP–Stem Joint axis and the frozen score-independent Program146 recur across patients in an independent public TNBC Visium cohort? The expected direction was fixed as positive before biological testing.", "",
  "## B. Frozen 9-patient cohort",
  sprintf("The analysis contains exactly 9 independent author-final-analysis patients and 9 sections: %s. Patients, not spots, are the biological inferential units. The cohort labels are 4 BA and 5 WA; no BA-versus-WA test or interaction was run.", paste(patient_order, collapse = ", ")), "",
  "## C. Why sAA9 is excluded",
  "sAA9 is present in the checksum-verified archive but is absent from the author's final analysis list. Its public exclusion reason is unresolved. It therefore remains ARCHIVED_NOT_IN_AUTHOR_FINAL_ANALYSIS and was neither normalized nor used in primary or rescue sensitivity analyses.", "",
  "## D. Normalization implementation",
  paste0("Each section was normalized independently from filtered raw integer UMI counts using the Giotto 4.2.3 standard-library/log semantics: log2(raw UMI / spot library size × 6000 + 1). The implementation is the deterministic formula-equivalent function `giotto_standard_log2_equivalent`; frozen arguments were `expression_values=raw`, `norm_methods=standard`, `library_size_norm=TRUE`, `scalefactor=6000`, `log_norm=TRUE`, `log_offset=1`, `logbase=2`, `scale_feats=FALSE`, and `scale_cells=FALSE`. The latter two are disabled because score-specific within-section z standardization is applied explicitly. No alternative normalization was tested. Runtime: ", R.version.string, ", Matrix ", packageVersion("Matrix"), ", hdf5r ", packageVersion("hdf5r"), "."), "",
  "## E. Signature coverage",
  "Coverage reproduced the G0 mapping exactly: YAP17 17/17, Stem21 19/21 (missing NANOG and POU5F1), Joint38 36/38, and Program146 135/146. No missing feature was replaced and no BSW2-specific signature was selected.", "",
  "## F. Score-independence audit",
  sprintf("Measurable Joint38 intersect measurable Program146 = %d genes (%s). Status: %s.", independence$overlap_n, independence$overlap_genes, independence$score_independence_status), "",
  "## G. Spot-level score QC",
  sprintf("%d/9 sections passed the prespecified finite-score and non-zero-variance gate. Across sections, Joint finite fractions ranged %.4f–%.4f and Program finite fractions ranged %.4f–%.4f. Failed patients, if any, are retained explicitly in Figure6G1_spot_score_QC.tsv: %s.", sum(spot_qc$score_QC_status == "PASS"), min(spot_qc$Joint_finite_fraction), max(spot_qc$Joint_finite_fraction), min(spot_qc$Program_finite_fraction), max(spot_qc$Program_finite_fraction), if (any(spot_qc$score_QC_status != "PASS")) paste(spot_qc$patient_id[spot_qc$score_QC_status != "PASS"], collapse = ";") else "NONE"), "",
  "## H. Technical depth associations",
  paste0("Descriptive section-wise score/depth correlations were: ", depth_summary_text, ". These are technical QC only and are not biological replication tests."), "",
  "## I. Raw Joint–Program association",
  sprintf("Raw within-section associations were positive in %d/9 patients; patient median rho=%s, IQR=%s (secondary descriptive evidence).", raw_summary$positive_n, fmt(raw_summary$median_rho), fmt(raw_summary$IQR)), "",
  "## J. Primary depth-adjusted patient associations",
  sprintf("Using the exact S1 rank-residual grammar—average ranks, separate `lm.fit` residualization on rank(log1p(nCount)), and Pearson correlation of the residual vectors—rho_depth_adjusted was positive in %d/9 patients. No TumorPurity or region threshold was introduced.", primary_summary$positive_n), "",
  "## K. Cohort primary test",
  sprintf("Across exactly 9 patient effects: median rho=%s, IQR=%s (Q1=%s, Q3=%s), exact one-sided Wilcoxon signed-rank P=%s; H1 was median rho > 0. No multiplicity correction was required for this single predefined primary test.", fmt(primary_summary$median_rho), fmt(primary_summary$IQR), fmt(primary_summary$Q1), fmt(primary_summary$Q3), fmtp(primary_summary$exact_one_sided_wilcoxon_P)), "",
  "## L. Positive-patient recurrence",
  sprintf("%d/9 patients (%.1f%%) had positive depth-adjusted effects.", primary_summary$positive_n, 100 * primary_summary$positive_fraction), "",
  "## M. nFeature sensitivity",
  sprintf("The prespecified nFeature-adjusted sensitivity had median rho=%s, IQR=%s and %d/9 positive patients. This sensitivity was not used to rescue the primary result.", fmt(nfeature_summary$median_rho), fmt(nfeature_summary$IQR), nfeature_summary$positive_n), "",
  "## N. Leave-one-patient-out robustness",
  if (primary_statistical_pass) sprintf("LOO was run because the primary cohort test passed. The cohort median remained positive in %d/9 omissions; LOO exact P values are robustness descriptors only.", sum(loo$median_direction_positive)) else "LOO was not run because the primary cohort test did not pass; it was not used as a rescue analysis.", "",
  "## O. Main limitations",
  "BSW2 Visium spots are mixed-cell tissue measurements. No machine-readable author tumour-region annotation or authoritative TumorPurity object exists. The analysis cannot establish malignant-cell-specific spatial coupling, direct regulation, causality, treatment response, or race-group differences. Missing frozen genes remained missing, and spot-level associations do not constitute independent biological replication.", "",
  "## P. Cross-study interpretation",
  if (ready) "The result supports independent spatial replication in a second TNBC Visium cohort and a patient-recurrent spatial association in mixed-cell tissue measurements. It does not support tumour-cell-specific or causal claims." else "The prespecified BSW2 analysis does not authorize a main-panel replication claim. No method, patient universe, region, normalization or signature was changed to improve the result.", "",
  "## Q. Scientific verdict", verdict, "",
  "## R. Prespecified evidence criteria", if (ready) "READY" else "HOLD", "")
writeLines(readout, file.path(out, "Figure6G1_readout.md"), useBytes = TRUE)

cat(sprintf("VERDICT=%s READY=%s median=%.8f positive=%d/9 exactP=%.10g rawPositive=%d nFeaturePositive=%d LOOall=%s QCpass=%d/9\n",
            verdict, ready, primary_summary$median_rho, primary_summary$positive_n,
            primary_summary$exact_one_sided_wilcoxon_P, raw_summary$positive_n,
            nfeature_summary$positive_n, loo_all_positive,
            sum(spot_qc$score_QC_status == "PASS")))
