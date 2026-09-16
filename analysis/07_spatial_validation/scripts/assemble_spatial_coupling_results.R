# Purpose: Assemble patient-recurrent spatial coupling results
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 6 and Supplementary Fig. S4.
options(stringsAsFactors = FALSE)
.libPaths(c(".runtime/Figure6_S0Q3_Rlib", .Library))
suppressPackageStartupMessages(library(hdf5r))

input_root <- ".runtime/Figure6_S1_inputs"
chunk_root <- ".runtime/Figure6_S1_chunks"
stage_root <- ".runtime/Figure6_S1_final"
dir.create(stage_root, recursive = TRUE, showWarnings = FALSE)

section_registry <- read.delim(file.path(input_root, "section_registry.tsv"), check.names = FALSE)
patient_registry <- read.delim(file.path(input_root, "patient_registry.tsv"), check.names = FALSE)
score_genes <- read.delim(file.path(input_root, "score_gene_list.tsv"), check.names = FALSE)
program146_original <- scan(file.path(input_root, "146_gene_program.txt"), what = character(), quiet = TRUE)
q3 <- read.delim(file.path(input_root, "Q3_ESTIMATE_spot_scores.tsv"), check.names = FALSE)
source_hashes <- read.delim(file.path(input_root, "source_hashes.tsv"), check.names = FALSE)

alias_map <- c(CYR61 = "CCN1", CTGF = "CCN2", `HNRNPU-AS1` = "HNRNPU")
map_alias <- function(x) {
  y <- x
  hit <- x %in% names(alias_map)
  y[hit] <- unname(alias_map[x[hit]])
  y
}
yap17_original <- score_genes$gene[score_genes$in_YAP_score]
stem21_original <- score_genes$gene[score_genes$in_Stemness_score]
joint38_original <- score_genes$gene[score_genes$in_YAP_Stem_construction]
yap17 <- unique(map_alias(yap17_original))
stem21 <- unique(map_alias(stem21_original))
joint38 <- unique(map_alias(joint38_original))
program146 <- unique(map_alias(program146_original))
program142 <- setdiff(program146, "HNRNPU")
stopifnot(length(yap17_original) == 17L, length(stem21_original) == 21L,
          length(joint38_original) == 38L, length(program146_original) == 146L,
          length(intersect(joint38, program146)) == 0L)

section_registry$patient_number <- as.integer(sub("^P", "", section_registry$patient_id))
section_registry <- section_registry[order(section_registry$patient_number, section_registry$section_number), ]
patient_registry$patient_number <- as.integer(sub("^P", "", patient_registry$patient_id))
patient_registry <- patient_registry[order(patient_registry$patient_number), ]

score_files <- file.path(chunk_root, paste0(section_registry$sample_id, "_", section_registry$section_id, "_scores.tsv"))
meta_files <- file.path(chunk_root, paste0(section_registry$sample_id, "_", section_registry$section_id, "_score_meta.tsv"))
error_files <- file.path(chunk_root, paste0(section_registry$sample_id, "_", section_registry$section_id, "_ERROR.txt"))
if (any(file.exists(error_files))) stop("Worker errors exist")
if (!all(file.exists(score_files)) || !all(file.exists(meta_files))) stop("Not all 43 S1 chunks are complete")
scores <- do.call(rbind, lapply(score_files, function(f) read.delim(f, check.names = FALSE)))
meta <- do.call(rbind, lapply(meta_files, function(f) read.delim(f, check.names = FALSE)))
stopifnot(nrow(scores) == 56567L, nrow(meta) == 43L, !anyDuplicated(paste(scores$sample_id, scores$barcode)))

q3_keep <- q3[, c("patient_id", "sample_id", "section_id", "barcode", "nCount_raw", "nFeature_raw", "TumorPurity")]
names(q3_keep)[5:6] <- c("Q3_nCount_raw", "Q3_nFeature_raw")
dat <- merge(scores, q3_keep, by = c("patient_id", "sample_id", "section_id", "barcode"), all = FALSE, sort = FALSE)
stopifnot(nrow(dat) == 56567L, !anyDuplicated(paste(dat$sample_id, dat$barcode)))
stopifnot(all(dat$nCount_raw == dat$Q3_nCount_raw), all(dat$nFeature_raw == dat$Q3_nFeature_raw), all(is.finite(dat$TumorPurity)))

meta$Q3_spot_join_n <- vapply(meta$sample_id, function(s) sum(dat$sample_id == s), integer(1))
meta$Q3_spot_join_status <- ifelse(meta$Q3_spot_join_n == meta$n_spots, "PASS", "FAIL")
meta$score_QC_status <- ifelse(meta$score_QC_status == "PASS" & meta$Q3_spot_join_status == "PASS", "PASS", "FAIL")
failed_sections <- sum(meta$score_QC_status != "PASS")
if (failed_sections > 4L) stop(sprintf("More than 10%% sections failed score QC: %d/43", failed_sections))

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
section_effect <- function(d, program_col) {
  ok <- is.finite(d$Joint_spatial) & is.finite(d[[program_col]]) & is.finite(d$TumorPurity) &
    is.finite(d$nCount_raw) & is.finite(d$nFeature_raw)
  z <- d[ok, ]
  ct <- suppressWarnings(cor.test(z$Joint_spatial, z[[program_col]], method = "spearman", exact = FALSE))
  c(n_spots = nrow(z), rho_raw = unname(ct$estimate), P_raw_section = ct$p.value,
    rho_tumour_adjusted = partial_rank_rho(z$Joint_spatial, z[[program_col]], list(TumorPurity = z$TumorPurity)),
    rho_tumour_depth_adjusted = partial_rank_rho(z$Joint_spatial, z[[program_col]], list(TumorPurity = z$TumorPurity, log1p_nCount = log1p(z$nCount_raw))),
    rho_tumour_nfeature_adjusted = partial_rank_rho(z$Joint_spatial, z[[program_col]], list(TumorPurity = z$TumorPurity, nFeature = z$nFeature_raw)))
}

section_rows <- lapply(seq_len(nrow(section_registry)), function(i) {
  sr <- section_registry[i, ]
  d <- dat[dat$sample_id == sr$sample_id, ]
  e <- section_effect(d, "Program146_spatial")
  data.frame(patient_id = sr$patient_id, sample_id = sr$sample_id, section_id = sr$section_id,
             section_number = sr$section_number, n_spots = as.integer(e["n_spots"]),
             rho_raw = e["rho_raw"], P_raw_section = e["P_raw_section"],
             rho_tumour_adjusted = e["rho_tumour_adjusted"],
             rho_tumour_depth_adjusted = e["rho_tumour_depth_adjusted"],
             rho_tumour_nfeature_adjusted = e["rho_tumour_nfeature_adjusted"],
             section_P_interpretation = "DESCRIPTIVE_ONLY_NOT_BIOLOGICAL_REPLICATION_TEST",
             stringsAsFactors = FALSE)
})
section_assoc <- do.call(rbind, section_rows)

section142_rows <- lapply(seq_len(nrow(section_registry)), function(i) {
  sr <- section_registry[i, ]
  d <- dat[dat$sample_id == sr$sample_id, ]
  e <- section_effect(d, "Program142_spatial")
  data.frame(record_type = "SECTION", patient_id = sr$patient_id, sample_id = sr$sample_id,
             section_id = sr$section_id, section_number = sr$section_number,
             n_spots = as.integer(e["n_spots"]), rho_tumour_depth_adjusted = e["rho_tumour_depth_adjusted"],
             patient_section_count = NA_integer_, patient_weight_sum = NA_real_,
             median_patient_rho = NA_real_, positive_patient_n = NA_integer_, positive_patient_fraction = NA_real_,
             one_sided_wilcoxon_P = NA_real_, status = "SECTION_SENSITIVITY", stringsAsFactors = FALSE)
})
section142 <- do.call(rbind, section142_rows)

aggregate_patient <- function(section_table, rho_col) {
  rows <- lapply(patient_registry$patient_id, function(pid) {
    s <- section_table[section_table$patient_id == pid & is.finite(section_table[[rho_col]]), ]
    if (!nrow(s)) return(data.frame(patient_id = pid, n_sections = 0L, total_spots = 0L, weight_sum = 0, patient_rho = NA_real_))
    rho <- pmax(pmin(s[[rho_col]], 1 - 1e-12), -1 + 1e-12)
    w <- s$n_spots - 3
    pz <- sum(w * atanh(rho)) / sum(w)
    data.frame(patient_id = pid, n_sections = nrow(s), total_spots = sum(s$n_spots), weight_sum = sum(w), patient_rho = tanh(pz))
  })
  do.call(rbind, rows)
}

model_cols <- c("rho_raw", "rho_tumour_adjusted", "rho_tumour_depth_adjusted", "rho_tumour_nfeature_adjusted")
patient_base <- data.frame(patient_id = patient_registry$patient_id, stringsAsFactors = FALSE)
for (m in model_cols) {
  a <- aggregate_patient(section_assoc, m)
  if (!"n_sections" %in% names(patient_base)) {
    patient_base$n_sections <- a$n_sections
    patient_base$total_spots <- a$total_spots
    patient_base$weight_sum <- a$weight_sum
  }
  patient_base[[m]] <- a$patient_rho
  patient_base[[paste0(m, "_direction")]] <- ifelse(a$patient_rho > 0, "POSITIVE", ifelse(a$patient_rho < 0, "NEGATIVE", "ZERO"))
}
patient_assoc <- patient_base

a142 <- aggregate_patient(section142, "rho_tumour_depth_adjusted")
patient142 <- data.frame(record_type = "PATIENT", patient_id = a142$patient_id,
                         sample_id = "PATIENT_AGGREGATE", section_id = "PATIENT_AGGREGATE", section_number = NA_integer_,
                         n_spots = a142$total_spots, rho_tumour_depth_adjusted = a142$patient_rho,
                         patient_section_count = a142$n_sections, patient_weight_sum = a142$weight_sum,
                         median_patient_rho = NA_real_, positive_patient_n = NA_integer_, positive_patient_fraction = NA_real_,
                         one_sided_wilcoxon_P = NA_real_, status = "PATIENT_SENSITIVITY", stringsAsFactors = FALSE)

run_wilcox <- function(x) {
  x <- x[is.finite(x)]
  exact_ok <- all(x != 0) && !anyDuplicated(abs(x))
  wt <- suppressWarnings(wilcox.test(x, mu = 0, alternative = "greater", exact = exact_ok, correct = FALSE))
  list(P = wt$p.value, method = wt$method, exact = exact_ok)
}
summarize_patient <- function(x) {
  x <- x[is.finite(x)]
  w <- run_wilcox(x)
  data.frame(n_patients = length(x), median_rho = median(x), Q1 = unname(quantile(x, 0.25)),
             Q3 = unname(quantile(x, 0.75)), IQR = IQR(x), minimum = min(x), maximum = max(x),
             positive_n = sum(x > 0), positive_fraction = mean(x > 0),
             one_sided_wilcoxon_P = w$P, wilcoxon_mode = if (w$exact) "EXACT" else "ASYMPTOTIC_TIES_OR_ZERO",
             hypothesis = "median patient rho > 0", stringsAsFactors = FALSE)
}

summaries <- lapply(model_cols, function(m) summarize_patient(patient_assoc[[m]]))
names(summaries) <- model_cols
summary142 <- summarize_patient(a142$patient_rho)
primary_summary <- summaries[["rho_tumour_depth_adjusted"]]
primary_basic_pass <- primary_summary$median_rho > 0 && primary_summary$one_sided_wilcoxon_P < 0.05

if (primary_basic_pass) {
  loo <- do.call(rbind, lapply(seq_len(nrow(patient_assoc)), function(i) {
    x <- patient_assoc$rho_tumour_depth_adjusted[-i]
    s <- summarize_patient(x)
    data.frame(status = "RUN_PRIMARY_PASSED", omitted_patient = patient_assoc$patient_id[i],
               n_patients = s$n_patients, median_rho = s$median_rho,
               positive_n = s$positive_n, positive_fraction = s$positive_fraction,
               one_sided_wilcoxon_P = s$one_sided_wilcoxon_P,
               direction_positive = s$median_rho > 0,
               primary_recurrence_criteria_preserved = s$median_rho > 0 && s$positive_fraction >= 0.70 && s$one_sided_wilcoxon_P < 0.05,
               stringsAsFactors = FALSE)
  }))
  loo_direction_all <- all(loo$direction_positive)
} else {
  loo <- data.frame(status = "NOT_RUN_PRIMARY_FAILED", omitted_patient = "NOT_APPLICABLE",
                    n_patients = NA_integer_, median_rho = NA_real_, positive_n = NA_integer_,
                    positive_fraction = NA_real_, one_sided_wilcoxon_P = NA_real_,
                    direction_positive = NA, primary_recurrence_criteria_preserved = NA,
                    stringsAsFactors = FALSE)
  loo_direction_all <- FALSE
}

criteria <- c(
  A = primary_summary$median_rho > 0,
  B = primary_summary$positive_fraction >= 0.70,
  C = primary_summary$one_sided_wilcoxon_P < 0.05,
  D = summaries[["rho_raw"]]$median_rho > 0,
  E = summaries[["rho_tumour_adjusted"]]$median_rho > 0,
  F = summaries[["rho_tumour_nfeature_adjusted"]]$median_rho > 0,
  G = summary142$median_rho > 0,
  H = loo_direction_all)
core_pass <- all(criteria)
strong_pass <- core_pass && primary_summary$positive_fraction >= 0.80 && primary_summary$one_sided_wilcoxon_P < 0.01
verdict <- if (strong_pass) "CORE_SPATIAL_STRONG_PASS" else if (core_pass) "CORE_SPATIAL_PASS" else if (primary_summary$median_rho > 0) "CORE_SPATIAL_BORDERLINE" else "CORE_SPATIAL_FAIL"
ready <- verdict %in% c("CORE_SPATIAL_STRONG_PASS", "CORE_SPATIAL_PASS")
retention_ratio <- primary_summary$median_rho / summaries[["rho_raw"]]$median_rho

primary_test <- do.call(rbind, lapply(c(model_cols, "Program142_rho_tumour_depth_adjusted"), function(m) {
  s <- if (m == "Program142_rho_tumour_depth_adjusted") summary142 else summaries[[m]]
  data.frame(model = m,
             role = if (m == "rho_tumour_depth_adjusted") "PRIMARY" else if (m == "rho_raw" || m == "rho_tumour_adjusted") "SUPPORTIVE" else "SENSITIVITY",
             s,
             primary_multiplicity_correction = if (m == "rho_tumour_depth_adjusted") "NOT_REQUIRED_ONE_PRIMARY_TEST" else "NOT_APPLICABLE_SUPPORTIVE",
             adjusted_to_raw_median_ratio = if (m == "rho_tumour_depth_adjusted") retention_ratio else NA_real_,
             criterion_A_median_positive = criteria["A"], criterion_B_positive_ge_70pct = criteria["B"],
             criterion_C_primary_P_lt_0.05 = criteria["C"], criterion_D_raw_median_positive = criteria["D"],
             criterion_E_tumour_adjusted_median_positive = criteria["E"],
             criterion_F_nfeature_median_positive = criteria["F"],
             criterion_G_program142_direction_positive = criteria["G"],
             criterion_H_no_single_patient_direction_determines = criteria["H"],
             core_spatial_verdict = verdict, stringsAsFactors = FALSE)
}))

cohort142 <- data.frame(record_type = "COHORT", patient_id = "ALL_22_PATIENTS", sample_id = "COHORT_SUMMARY",
                        section_id = "COHORT_SUMMARY", section_number = NA_integer_, n_spots = sum(a142$total_spots),
                        rho_tumour_depth_adjusted = NA_real_, patient_section_count = sum(a142$n_sections),
                        patient_weight_sum = sum(a142$weight_sum), median_patient_rho = summary142$median_rho,
                        positive_patient_n = summary142$positive_n, positive_patient_fraction = summary142$positive_fraction,
                        one_sided_wilcoxon_P = summary142$one_sided_wilcoxon_P,
                        status = if (summary142$median_rho > 0) "SAME_POSITIVE_COHORT_DIRECTION" else "DIRECTION_NOT_PRESERVED",
                        stringsAsFactors = FALSE)
program142_sensitivity <- rbind(section142, patient142, cohort142)

paired_ids <- patient_registry$patient_id[patient_registry$n_sections == 2L]
paired <- do.call(rbind, lapply(paired_ids, function(pid) {
  s <- section_assoc[section_assoc$patient_id == pid, ]
  s <- s[order(s$section_number), ]
  r1 <- s$rho_tumour_depth_adjusted[1]; r2 <- s$rho_tumour_depth_adjusted[2]
  data.frame(patient_id = pid, section_1_sample = s$sample_id[1], section_1_id = s$section_id[1],
             section_1_rho = r1, section_1_direction = ifelse(r1 > 0, "POSITIVE", ifelse(r1 < 0, "NEGATIVE", "ZERO")),
             section_2_sample = s$sample_id[2], section_2_id = s$section_id[2],
             section_2_rho = r2, section_2_direction = ifelse(r2 > 0, "POSITIVE", ifelse(r2 < 0, "NEGATIVE", "ZERO")),
             same_sign = sign(r1) == sign(r2), stringsAsFactors = FALSE)
}))
paired$same_sign_n_of_21 <- sum(paired$same_sign)
paired$same_sign_fraction <- mean(paired$same_sign)

depth_rows <- lapply(seq_len(nrow(section_registry)), function(i) {
  sr <- section_registry[i, ]
  d <- dat[dat$sample_id == sr$sample_id, ]
  data.frame(record_type = "SECTION", patient_id = sr$patient_id, sample_id = sr$sample_id, section_id = sr$section_id,
             n_spots = nrow(d),
             rho_Joint_log1p_nCount = cor(d$Joint_spatial, log1p(d$nCount_raw), method = "spearman"),
             rho_Program146_log1p_nCount = cor(d$Program146_spatial, log1p(d$nCount_raw), method = "spearman"),
             rho_Joint_nFeature = cor(d$Joint_spatial, d$nFeature_raw, method = "spearman"),
             rho_Program146_nFeature = cor(d$Program146_spatial, d$nFeature_raw, method = "spearman"),
             summary_metric = NA_character_, summary_median = NA_real_, summary_Q1 = NA_real_,
             summary_Q3 = NA_real_, summary_IQR = NA_real_, summary_max_abs = NA_real_, stringsAsFactors = FALSE)
})
depth_section <- do.call(rbind, depth_rows)
depth_metrics <- c("rho_Joint_log1p_nCount", "rho_Program146_log1p_nCount", "rho_Joint_nFeature", "rho_Program146_nFeature")
depth_summary <- do.call(rbind, lapply(depth_metrics, function(m) {
  x <- depth_section[[m]]
  data.frame(record_type = "SUMMARY", patient_id = "ALL", sample_id = "ALL_43_SECTIONS", section_id = "ALL",
             n_spots = sum(depth_section$n_spots),
             rho_Joint_log1p_nCount = NA_real_, rho_Program146_log1p_nCount = NA_real_,
             rho_Joint_nFeature = NA_real_, rho_Program146_nFeature = NA_real_,
             summary_metric = m, summary_median = median(x), summary_Q1 = unname(quantile(x, 0.25)),
             summary_Q3 = unname(quantile(x, 0.75)), summary_IQR = IQR(x), summary_max_abs = max(abs(x)),
             stringsAsFactors = FALSE)
}))
depth_qc <- rbind(depth_section, depth_summary)

hash_value <- function(label) {
  z <- source_hashes$sha256[source_hashes$label == label]
  if (length(z)) z[1] else "NOT_AVAILABLE"
}
first_h5_dir <- file.path(".runtime/Figure6_S0Q2_geo", section_registry$sample_id[1])
first_h5 <- list.files(first_h5_dir, pattern = "filtered_feature_bc_matrix[.]h5$", full.names = TRUE)
h <- hdf5r::H5File$new(first_h5, mode = "r")
raw_features <- make.unique(h[["matrix/features/name"]][])
h$close_all()

is_missing_in_meta <- function(gene, strings) {
  vapply(strings, function(s) s != "NONE" && gene %in% strsplit(s, ";", fixed = TRUE)[[1]], logical(1))
}
manifest_object <- function(object, original, normalized, missing_col, role) {
  do.call(rbind, lapply(seq_along(original), function(i) {
    g0 <- original[i]; g <- normalized[i]
    missing_n <- sum(is_missing_in_meta(g, meta[[missing_col]]))
    data.frame(object = object, original_gene = g0, normalized_feature = g,
               alias_action = if (g0 == g) "NONE" else paste0(g0, "->", g),
               frozen_membership = TRUE,
               primary_score_membership = if (object == "Program146") TRUE else NA,
               Program142_sensitivity_membership = if (object == "Program146") g != "HNRNPU" else NA,
               raw_feature_available = g %in% raw_features,
               SCT_measurable_section_n = 43L - missing_n,
               SCT_measurable_section_fraction = (43L - missing_n) / 43,
               score_role = role,
               normalization = "section-wise SCTransform 0.3.2 data-slot log1p(corrected UMI)",
               source_file = if (object == "Program146") "146_gene_program.txt" else "score_gene_list.tsv",
               source_sha256 = if (object == "Program146") hash_value("146_gene_program.txt") else hash_value("score_gene_list.tsv"),
               stringsAsFactors = FALSE)
  }))
}
score_manifest <- rbind(
  manifest_object("YAP17", yap17_original, map_alias(yap17_original), "YAP17_missing", "mean normalized expression then within-section z"),
  manifest_object("Stem21", stem21_original, map_alias(stem21_original), "Stem21_missing", "mean normalized expression then within-section z"),
  manifest_object("Program146", program146_original, map_alias(program146_original), "Program146_missing", "mean of gene-wise within-section z scores"))

shared_lists <- list(
  intersect(yap17, program146),
  intersect(stem21, program146),
  intersect(joint38, program146),
  intersect(joint38, program142)
)
shared_audit <- data.frame(
  comparison = c("YAP17_vs_Program146", "Stem21_vs_Program146", "Joint38_vs_Program146", "Joint38_vs_Program142_sensitivity"),
  identity_policy = "AUTHORIZED_ALIAS_NORMALIZED",
  left_n = c(length(yap17), length(stem21), length(joint38), length(joint38)),
  right_frozen_definition_n = c(length(program146), length(program146), length(program146), length(program142)),
  right_raw_measurable_n = c(sum(program146 %in% raw_features), sum(program146 %in% raw_features),
                             sum(program146 %in% raw_features), sum(program142 %in% raw_features)),
  n_shared = lengths(shared_lists),
  shared_genes = vapply(shared_lists, function(x) if (length(x)) paste(x, collapse = ";") else "NONE", character(1)),
  status = "PASS_ZERO_SHARED_GENES", stringsAsFactors = FALSE)

fmt <- function(x, d = 4) formatC(x, format = "f", digits = d)
readout <- c(
  "# Figure 6 S1 — Patient-recurrent spatial coupling", "",
  "## A. Scientific question", "Does the frozen YAP–Stem Joint axis spatially co-vary with the frozen score-independent Program146 within TNBC tissues, recurrently across patients after continuous tumour-composition and technical-depth adjustment?", "",
  "## B. Cohort and biological replication unit", "GSE210616: 22 TNBC patients, 43 Visium sections and 56,567 in-tissue spots. The biological replication unit is the patient; spots are nested within sections and sections within patients.", "",
  "## C. Frozen score definitions", "Joint_spatial = within-section z(YAP17 mean normalized expression) + within-section z(Stem21 mean normalized expression). Program146_spatial = the mean of measurable gene-wise within-section z scores. No control genes, gene reselection, integration or batch correction was used.", "",
  "## D. Gene coverage", sprintf("Raw-feature coverage was YAP17 17/17, Stem21 21/21 and Program146 143/146 under the authorized alias policy. The conservative sensitivity contained 142/146 raw-measurable features after excluding the HNRNPU mapping. Section-wise SCT-measurable Program146 genes ranged from %d to %d; the conservative sensitivity ranged from %d to %d because the fixed section-wise SCT model does not model features detected in fewer than five spots. Three unavailable frozen programme features remained unavailable; no genes were selected after viewing the association results.", min(meta$Program146_measurable_n), max(meta$Program146_measurable_n), min(meta$Program142_measurable_n), max(meta$Program142_measurable_n)), "",
  "## E. Shared-gene audit", "After authorized alias normalization, Joint38 ∩ Program146 = 0 genes and Joint38 ∩ Program142 sensitivity = 0 genes. No defining-score genes overlap the tested programs.", "",
  "## F. Section-level score QC", sprintf("%d/43 sections passed the prespecified finite-score and non-zero-variance gate; %d failed sections were retained in the QC table and not silently removed.", sum(meta$score_QC_status == "PASS"), failed_sections), "",
  "## G. Raw spatial association", sprintf("Across 22 patient-level Fisher-z aggregates: median rho=%s, IQR=%s, positive=%d/22; one-sided Wilcoxon P=%s (supportive).", fmt(summaries[["rho_raw"]]$median_rho), fmt(summaries[["rho_raw"]]$IQR), summaries[["rho_raw"]]$positive_n, format(summaries[["rho_raw"]]$one_sided_wilcoxon_P, scientific = TRUE, digits = 4)), "",
  "## H. Tumour-composition-adjusted association", sprintf("Median patient rho=%s, IQR=%s, positive=%d/22; one-sided Wilcoxon P=%s (supportive). TumorPurity was continuous and never thresholded.", fmt(summaries[["rho_tumour_adjusted"]]$median_rho), fmt(summaries[["rho_tumour_adjusted"]]$IQR), summaries[["rho_tumour_adjusted"]]$positive_n, format(summaries[["rho_tumour_adjusted"]]$one_sided_wilcoxon_P, scientific = TRUE, digits = 4)), "",
  "## I. Tumour+depth-adjusted association", sprintf("Primary model: median patient rho=%s, IQR=%s, positive=%d/22 (%.1f%%), one-sided Wilcoxon P=%s. Adjusted/raw median ratio=%s; no post-hoc retention cutoff was imposed.", fmt(primary_summary$median_rho), fmt(primary_summary$IQR), primary_summary$positive_n, 100*primary_summary$positive_fraction, format(primary_summary$one_sided_wilcoxon_P, scientific = TRUE, digits = 4), fmt(retention_ratio)), "",
  "## J. nFeature sensitivity", sprintf("Median patient rho=%s, IQR=%s, positive=%d/22; one-sided Wilcoxon P=%s.", fmt(summaries[["rho_tumour_nfeature_adjusted"]]$median_rho), fmt(summaries[["rho_tumour_nfeature_adjusted"]]$IQR), summaries[["rho_tumour_nfeature_adjusted"]]$positive_n, format(summaries[["rho_tumour_nfeature_adjusted"]]$one_sided_wilcoxon_P, scientific = TRUE, digits = 4)), "",
  "## K. Patient-level primary result", sprintf("The primary inference used exactly 22 patient-level rho values and a one-sided one-sample Wilcoxon signed-rank test against zero. Result: median=%s; P=%s.", fmt(primary_summary$median_rho), format(primary_summary$one_sided_wilcoxon_P, scientific = TRUE, digits = 4)), "",
  "## L. Positive-patient recurrence", sprintf("%d/22 patients (%.1f%%) had positive tumour+depth-adjusted rho.", primary_summary$positive_n, 100*primary_summary$positive_fraction), "",
  "## M. Paired-section consistency", sprintf("The two sections had the same primary-effect sign in %d/21 paired patients (%.1f%%). This is descriptive secondary evidence, not the cohort test.", sum(paired$same_sign), 100*mean(paired$same_sign)), "",
  "## N. Program142 alias sensitivity", sprintf("After excluding the HNRNPU-AS1/HNRNPU mapped feature, median patient rho=%s, positive=%d/22 (%.1f%%), one-sided Wilcoxon P=%s; cohort direction status=%s.", fmt(summary142$median_rho), summary142$positive_n, 100*summary142$positive_fraction, format(summary142$one_sided_wilcoxon_P, scientific = TRUE, digits = 4), cohort142$status), "",
  "## O. Leave-one-patient-out robustness", if (primary_basic_pass) sprintf("LOO was run because the primary test passed. Median direction remained positive in %d/22 omissions; full primary recurrence criteria were retained in %d/22 omissions.", sum(loo$direction_positive), sum(loo$primary_recurrence_criteria_preserved)) else "NOT_RUN_PRIMARY_FAILED", "",
  "## P. Technical-depth score QC", paste(vapply(depth_metrics, function(m) { z <- depth_summary[depth_summary$summary_metric == m, ]; sprintf("%s: median rho=%s, IQR=%s, max |rho|=%s", m, fmt(z$summary_median), fmt(z$summary_IQR), fmt(z$summary_max_abs)) }, character(1)), collapse = "; "), "",
  "## Q. Main limitations", "Visium spots are multicellular tissue measurements, not malignant cells. ESTIMATE TumorPurity is expression-derived and not a direct malignant-cell fraction. Some frozen programme genes were absent from the per-section SCT model because they were not sufficiently detected; this was recorded without post-result gene selection. Section-level P values are descriptive only. The analysis supports association, not causality, direct YAP regulation, metastasis or clinical prediction.", "",
  "## R. Core spatial verdict", verdict, "",
  "## S. Prespecified evidence criteria", if (ready) "READY" else "HOLD")

write_tsv <- function(x, name) write.table(x, file.path(stage_root, name), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
write_tsv(score_manifest, "Figure6_S1_score_manifest.tsv")
write_tsv(meta, "Figure6_S1_section_score_QC.tsv")
write_tsv(shared_audit, "Figure6_S1_shared_gene_audit.tsv")
write_tsv(section_assoc, "Figure6_S1_section_associations.tsv")
write_tsv(patient_assoc, "Figure6_S1_patient_associations.tsv")
write_tsv(primary_test, "Figure6_S1_primary_test.tsv")
write_tsv(paired, "Figure6_S1_paired_section_consistency.tsv")
write_tsv(program142_sensitivity, "Figure6_S1_program142_sensitivity.tsv")
write_tsv(depth_qc, "Figure6_S1_depth_score_QC.tsv")
write_tsv(loo, "Figure6_S1_LOO.tsv")
writeLines(readout, file.path(stage_root, "Figure6_S1_readout.md"), useBytes = TRUE)
cat(sprintf("VERDICT=%s READY=%s spots=%d sections=%d patients=%d files=%d primary_median=%.6f primary_P=%.10g positive=%d\n",
            verdict, ready, nrow(dat), nrow(section_assoc), nrow(patient_assoc), length(list.files(stage_root)),
            primary_summary$median_rho, primary_summary$one_sided_wilcoxon_P, primary_summary$positive_n))
