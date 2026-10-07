options(stringsAsFactors = FALSE, scipen = 999)
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(hdf5r)
  library(digest)
})

root <- Sys.getenv("AHIPPO_ANALYSIS_ROOT", unset = ".")
out <- Sys.getenv("AHIPPO_RESULTS_DIR", unset = file.path(root, "outputs/artemis_functional_programs"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) fread(path, sep = "\t", na.strings = c("NA", "NaN", ""))
write_tsv <- function(x, name) {
  fwrite(as.data.table(x), file.path(out, name), sep = "\t", quote = FALSE, na = "NA")
}
bh <- function(p) {
  ans <- rep(NA_real_, length(p)); ok <- is.finite(p)
  ans[ok] <- p.adjust(p[ok], method = "BH"); ans
}
safe_wilcox <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x) || all(abs(x) < 1e-15)) return(1)
  suppressWarnings(wilcox.test(x, mu = 0, exact = FALSE)$p.value)
}
safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(c(rho = NA_real_, P = NA_real_, n = sum(ok)))
  z <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  c(rho = unname(z$estimate), P = z$p.value, n = sum(ok))
}

paths <- c(
  h5ad = file.path(root, "scouter/af8c4fce-4c63-4671-b339-91a383cf36f6.h5ad"),
  cache = file.path(root, "0716_yan2026_validation/inputs/yan_cancer_cell_scores.tsv.gz"),
  module_audit = file.path(root, "Figure4D_continuous_functional_trajectories/Figure4D_module_definition_audit.tsv"),
  wu_module_rho = file.path(root, "Figure4D_continuous_functional_trajectories/Figure4D_patient_program_rho.tsv"),
  program_source = file.path(root, "Figure3_bottom_panels_HI_feasibility_and_plotting/Figure3_I_effect_consistency_source.tsv"),
  yan_state_sensitivity = file.path(root, "0728_sensitivity_and_cnv_closure/06_Yan_sensitivity/Yan_primary146_state_definition_sensitivity.tsv"),
  gse_transfer = file.path(root, "Figure4F_GSE180286_program_transfer/Figure4F_program_gene_transfer.tsv")
)
stopifnot(all(file.exists(paths)))

cat("Reading frozen definitions...\n")
module_audit <- read_tsv(paths[["module_audit"]])
modules <- module_audit[main_panel_decision == "SELECTED_MAIN"]
stopifnot(nrow(modules) == 6L)
modules[, module_order := match(internal_id, c(
  "hallmark_tnfa_nfkb", "hallmark_emt", "hallmark_hypoxia", "hallmark_apoptosis",
  "custom_survival_stress", "custom_adhesion_remodeling"
))]
setorder(modules, module_order)
stopifnot(!anyNA(modules$module_order))

program_source <- read_tsv(paths[["program_source"]])
program <- program_source[program_146 == TRUE]
stopifnot(nrow(program) == 146L, all(program$overall_log2FC_High_vs_Other > 0))
program_genes <- program$gene

module_genes <- setNames(
  lapply(modules$retained_genes_after_exclusion, function(x) unique(strsplit(x, ";", fixed = TRUE)[[1]])),
  modules$internal_id
)
desired_genes <- unique(c(unlist(module_genes, use.names = FALSE), program_genes))

cat("Reading score cache and selecting the frozen Yan eligibility set...\n")
cache <- read_tsv(paths[["cache"]])
eligible_patients <- cache[, .N, by = patient][N >= 50, patient]
eligible_pos <- which(cache$patient %in% eligible_patients)
d <- cache[eligible_pos]
stopifnot(uniqueN(d$patient) == 78L, nrow(d) == 48885L)

cat("Opening H5AD and resolving cell/gene indices...\n")
f <- H5File$new(paths[["h5ad"]], mode = "r")
on.exit(f$close_all(), add = TRUE)
obs_index <- f[["obs/_index"]][]
var_symbol <- f[["var/gene_symbols"]][]
cell_rows <- match(d$cell_id, obs_index)
stopifnot(!anyNA(cell_rows), !anyDuplicated(cell_rows))
symbol_to_col <- match(desired_genes, var_symbol)
coverage <- data.table(
  gene = desired_genes,
  present_in_Yan_H5AD = !is.na(symbol_to_col),
  h5ad_var_index_zero_based = fifelse(!is.na(symbol_to_col), symbol_to_col - 1L, NA_integer_),
  in_primary146 = desired_genes %in% program_genes,
  module_membership = vapply(desired_genes, function(g) {
    hit <- names(module_genes)[vapply(module_genes, function(x) g %in% x, logical(1))]
    if (length(hit)) paste(hit, collapse = ";") else NA_character_
  }, character(1))
)
write_tsv(coverage, "Fig5_Yan_gene_coverage.tsv")

present <- coverage[present_in_Yan_H5AD == TRUE]
desired_cols0 <- present$h5ad_var_index_zero_based
gene_order <- present$gene
stopifnot(!anyDuplicated(desired_cols0))

extract_sparse <- function(f, rows0, cols0) {
  input_order <- order(rows0)
  sorted_rows <- rows0[input_order]
  breaks <- which(diff(sorted_rows) > 1L)
  starts <- c(1L, breaks + 1L)
  ends <- c(breaks, length(sorted_rows))
  indptr <- as.numeric(f[["X/indptr"]][])
  i_parts <- vector("list", length(sorted_rows))
  j_parts <- vector("list", length(sorted_rows))
  x_parts <- vector("list", length(sorted_rows))
  out_cursor <- 0L
  for (run_index in seq_along(starts)) {
    s <- starts[run_index]; e <- ends[run_index]
    h5_start <- sorted_rows[s]
    h5_end <- sorted_rows[e] + 1L
    pointer_start <- indptr[h5_start + 1L]
    pointer_end <- indptr[h5_end + 1L]
    if (pointer_end <= pointer_start) next
    slice <- (pointer_start + 1):pointer_end
    run_data <- f[["X/data"]][slice]
    run_cols <- f[["X/indices"]][slice]
    for (local in seq_len(h5_end - h5_start)) {
      global0 <- h5_start + local - 1L
      p0 <- indptr[global0 + 1L] - pointer_start + 1L
      p1 <- indptr[global0 + 2L] - pointer_start
      row_pos <- s + local - 1L
      if (p1 < p0) next
      idx <- p0:p1
      mapped <- match(run_cols[idx], cols0, nomatch = 0L)
      keep <- mapped > 0L
      if (any(keep)) {
        out_cursor <- out_cursor + 1L
        i_parts[[out_cursor]] <- rep.int(row_pos, sum(keep))
        j_parts[[out_cursor]] <- mapped[keep]
        x_parts[[out_cursor]] <- as.numeric(run_data[idx][keep])
      }
    }
    if (run_index %% 10L == 0L || run_index == length(starts)) {
      cat(sprintf("  extracted run %d/%d\n", run_index, length(starts)))
    }
  }
  i <- unlist(i_parts[seq_len(out_cursor)], use.names = FALSE)
  j <- unlist(j_parts[seq_len(out_cursor)], use.names = FALSE)
  x <- unlist(x_parts[seq_len(out_cursor)], use.names = FALSE)
  mat_sorted <- sparseMatrix(i = i, j = j, x = x, dims = c(length(rows0), length(cols0)))
  inverse <- order(input_order)
  mat_sorted[inverse, , drop = FALSE]
}

cat(sprintf("Extracting %d genes for %d eligible malignant cells...\n", length(gene_order), nrow(d)))
expr <- extract_sparse(f, cell_rows - 1L, desired_cols0)
colnames(expr) <- gene_order
rownames(expr) <- d$cell_id
f$close_all()
cat(sprintf("Extraction complete: %d x %d, nonzero=%d\n", nrow(expr), ncol(expr), length(expr@x)))

residualize_vector <- function(values, frame) {
  outv <- rep(NA_real_, length(values))
  for (pt in unique(frame$patient)) {
    pos <- which(frame$patient == pt)
    X <- cbind(
      1,
      log1p(frame$nCount_RNA[pos]),
      log1p(frame$nFeature_RNA[pos]),
      frame$S_score[pos],
      frame$G2M_score[pos]
    )
    fit <- lm.fit(X, values[pos])
    outv[pos] <- fit$residuals
  }
  outv
}

patient_correlations <- function(feature_name, values, internal_id = NA_character_, display_name = NA_character_) {
  rec <- vector("list", uniqueN(d$patient) * 2L)
  k <- 0L
  for (pt in sort(unique(d$patient))) {
    pos <- which(d$patient == pt)
    x <- d$YAP_Stem_axis[pos]
    y <- values[pos]
    raw <- safe_spearman(x, y)
    X <- cbind(1, log1p(d$nCount_RNA[pos]), log1p(d$nFeature_RNA[pos]), d$S_score[pos], d$G2M_score[pos])
    xr <- lm.fit(X, x)$residuals
    yr <- lm.fit(X, y)$residuals
    adj <- safe_spearman(xr, yr)
    for (model in c("raw", "technical_adjusted")) {
      z <- if (model == "raw") raw else adj
      k <- k + 1L
      rec[[k]] <- data.table(
        record_type = "patient", patient = pt, n_cells = length(pos),
        feature = feature_name, internal_id = internal_id, display_name = display_name,
        model = model, rho = z[["rho"]], descriptive_cell_level_P = z[["P"]],
        n_complete_cells = z[["n"]],
        inference_guardrail = "patient-level effect retained; cell-level P is descriptive and not used as cross-patient biological replication"
      )
    }
  }
  rbindlist(rec[seq_len(k)], fill = TRUE)
}

# -----------------------------------------------------------------------------
# Candidates 6 and 7: current six frozen modules in Yan and cross-cohort mapping.
# -----------------------------------------------------------------------------
cat("Scoring the six current frozen modules...\n")
module_scores <- list()
module_coverage <- vector("list", nrow(modules))
for (i in seq_len(nrow(modules))) {
  id <- modules$internal_id[i]
  requested <- module_genes[[id]]
  available <- requested[requested %in% colnames(expr)]
  if (!length(available)) stop("No Yan genes available for module: ", id)
  module_scores[[id]] <- Matrix::rowMeans(expr[, available, drop = FALSE])
  module_coverage[[i]] <- data.table(
    internal_id = id, display_name = modules$display_name[i], requested_gene_n = length(requested),
    detected_gene_n = length(available), missing_gene_n = length(setdiff(requested, available)),
    detected_genes = paste(available, collapse = ";"),
    missing_genes = if (length(setdiff(requested, available))) paste(setdiff(requested, available), collapse = ";") else NA_character_
  )
}
write_tsv(rbindlist(module_coverage), "Fig5_current_module_Yan_coverage.tsv")

module_patient <- rbindlist(lapply(seq_len(nrow(modules)), function(i) {
  id <- modules$internal_id[i]
  patient_correlations(modules$display_name[i], module_scores[[id]], id, modules$display_name[i])
}))

module_summary <- module_patient[, .(
  n_testable = sum(is.finite(rho)), n_positive = sum(rho > 0, na.rm = TRUE),
  n_negative = sum(rho < 0, na.rm = TRUE), median_rho = median(rho, na.rm = TRUE),
  q1_rho = quantile(rho, .25, na.rm = TRUE), q3_rho = quantile(rho, .75, na.rm = TRUE),
  min_rho = min(rho, na.rm = TRUE), max_rho = max(rho, na.rm = TRUE),
  patient_unit_wilcoxon_P = safe_wilcox(rho)
), by = .(internal_id, display_name, model)]
module_summary[, BH_FDR := bh(patient_unit_wilcoxon_P), by = model]
module_summary[, `:=`(record_type = "module_summary", patient = "ALL_ELIGIBLE")]
write_tsv(rbindlist(list(module_patient, module_summary), fill = TRUE), "Fig5_module_patient_landscape.tsv")

wu <- read_tsv(paths[["wu_module_rho"]])
wu <- wu[internal_id %in% modules$internal_id & valid_status == "VALID"]
wu_summary <- wu[, .(
  Wu_n_patients = uniqueN(patient), Wu_median_rho = median(spearman_rho, na.rm = TRUE),
  Wu_q1_rho = quantile(spearman_rho, .25, na.rm = TRUE), Wu_q3_rho = quantile(spearman_rho, .75, na.rm = TRUE),
  Wu_positive_n = sum(spearman_rho > 0, na.rm = TRUE)
), by = .(internal_id, Wu_display_name = module)]
yan_raw <- module_summary[model == "raw", .(
  internal_id, Yan_raw_n_patients = n_testable, Yan_raw_median_rho = median_rho,
  Yan_raw_q1_rho = q1_rho, Yan_raw_q3_rho = q3_rho, Yan_raw_positive_n = n_positive,
  Yan_raw_BH_FDR = BH_FDR
)]
yan_adj <- module_summary[model == "technical_adjusted", .(
  internal_id, Yan_adjusted_n_patients = n_testable, Yan_adjusted_median_rho = median_rho,
  Yan_adjusted_q1_rho = q1_rho, Yan_adjusted_q3_rho = q3_rho, Yan_adjusted_positive_n = n_positive,
  Yan_adjusted_BH_FDR = BH_FDR
)]
cross <- Reduce(function(x, y) merge(x, y, by = "internal_id", all = TRUE), list(wu_summary, yan_raw, yan_adj))
cross[, `:=`(
  Wu_Yan_raw_same_direction = sign(Wu_median_rho) == sign(Yan_raw_median_rho),
  Wu_Yan_adjusted_same_direction = sign(Wu_median_rho) == sign(Yan_adjusted_median_rho),
  scale_guardrail = "within-cohort Spearman rho; compare direction/rank, not absolute module expression"
)]
raw_conc <- safe_spearman(cross$Wu_median_rho, cross$Yan_raw_median_rho)
adj_conc <- safe_spearman(cross$Wu_median_rho, cross$Yan_adjusted_median_rho)
cross_summary <- data.table(
  internal_id = "ALL_SIX_MODULES", record_type = "cross_module_summary",
  Wu_Yan_raw_rank_spearman = raw_conc[["rho"]], Wu_Yan_raw_rank_P = raw_conc[["P"]],
  Wu_Yan_adjusted_rank_spearman = adj_conc[["rho"]], Wu_Yan_adjusted_rank_P = adj_conc[["P"]],
  Wu_Yan_raw_direction_concordant_n = sum(cross$Wu_Yan_raw_same_direction, na.rm = TRUE),
  Wu_Yan_adjusted_direction_concordant_n = sum(cross$Wu_Yan_adjusted_same_direction, na.rm = TRUE),
  n_modules = nrow(cross),
  inference_guardrail = "six prespecified modules are features; rank correlation is descriptive with n=6"
)
cross[, record_type := "module"]
write_tsv(rbindlist(list(cross, cross_summary), fill = TRUE), "Fig5_module_crosscohort_concordance.tsv")

# -----------------------------------------------------------------------------
# Candidate 8: fixed strict_v2 146-gene program transfer to Yan.
# -----------------------------------------------------------------------------
program_present <- program_genes[program_genes %in% colnames(expr)]
stopifnot(length(program_present) == 141L)
program_score <- Matrix::rowMeans(expr[, program_present, drop = FALSE])
program_cont <- patient_correlations("strict_v2 146-gene program score", program_score,
                                     "strict_v2_146", "Strict-v2 program")
program_cont[, analysis_metric := "continuous_joint_axis_correlation"]
program_summary <- program_cont[, .(
  n_testable = sum(is.finite(rho)), n_positive = sum(rho > 0, na.rm = TRUE),
  n_negative = sum(rho < 0, na.rm = TRUE), median_effect = median(rho, na.rm = TRUE),
  q1_effect = quantile(rho, .25, na.rm = TRUE), q3_effect = quantile(rho, .75, na.rm = TRUE),
  patient_unit_wilcoxon_P = safe_wilcox(rho)
), by = model]
program_summary[, `:=`(
  record_type = "summary", patient = "ALL_ELIGIBLE", feature = "strict_v2 146-gene program score",
  internal_id = "strict_v2_146", display_name = "Strict-v2 program",
  analysis_metric = "continuous_joint_axis_correlation"
)]

state_existing <- read_tsv(paths[["yan_state_sensitivity"]])[
  definition == "Primary pooled tertile" & model %in% c("raw", "technical_adjusted")
]
state_existing[, `:=`(
  analysis_metric = "High_vs_Other_program_score_difference",
  source_note = "copied unchanged from authoritative Yan state-definition sensitivity",
  feature = "strict_v2 146-gene program score", internal_id = "strict_v2_146",
  display_name = "Strict-v2 program"
)]
write_tsv(rbindlist(list(program_cont, program_summary, state_existing), fill = TRUE), "Fig5_strictv2_program_transfer.tsv")

# -----------------------------------------------------------------------------
# Candidate 9: gene-level Yan conservation, with patient-level High-vs-Other
# mean-expression effects. The same frozen state and technical covariates are used.
# -----------------------------------------------------------------------------
cat("Computing patient-level effects for all 141 measurable program genes...\n")
gene_effect_parts <- vector("list", uniqueN(d$patient) * 2L)
k <- 0L
for (pt in sort(unique(d$patient))) {
  pos <- which(d$patient == pt)
  high <- d$YS_state[pos] == "High"
  other <- !high
  Y <- as.matrix(expr[pos, program_present, drop = FALSE])
  raw_delta <- if (any(high) && any(other)) colMeans(Y[high, , drop = FALSE]) - colMeans(Y[other, , drop = FALSE]) else rep(NA_real_, ncol(Y))
  X <- cbind(1, log1p(d$nCount_RNA[pos]), log1p(d$nFeature_RNA[pos]), d$S_score[pos], d$G2M_score[pos])
  Yres <- lm.fit(X, Y)$residuals
  adj_delta <- if (any(high) && any(other)) colMeans(Yres[high, , drop = FALSE]) - colMeans(Yres[other, , drop = FALSE]) else rep(NA_real_, ncol(Y))
  for (model in c("raw", "technical_adjusted")) {
    k <- k + 1L
    gene_effect_parts[[k]] <- data.table(
      record_type = "patient_gene", patient = pt, n_cells = length(pos),
      n_high = sum(high), n_other = sum(other), gene = program_present, model = model,
      patient_High_vs_Other_effect = if (model == "raw") raw_delta else adj_delta,
      effect_direction = fifelse((if (model == "raw") raw_delta else adj_delta) > 0, "positive",
                          fifelse((if (model == "raw") raw_delta else adj_delta) < 0, "negative", "zero")),
      valid = any(high) && any(other),
      patient_effect_definition = if (model == "raw") "mean log-normalized expression in High minus Other within Yan patient" else "mean technical-residual expression in High minus Other within Yan patient"
    )
  }
}
gene_patient <- rbindlist(gene_effect_parts[seq_len(k)])
gene_summary <- gene_patient[, .(
  n_testable = sum(valid & is.finite(patient_High_vs_Other_effect)),
  n_positive = sum(patient_High_vs_Other_effect > 0, na.rm = TRUE),
  n_negative = sum(patient_High_vs_Other_effect < 0, na.rm = TRUE),
  n_zero = sum(patient_High_vs_Other_effect == 0, na.rm = TRUE),
  median_Yan_effect = median(patient_High_vs_Other_effect, na.rm = TRUE),
  q1_Yan_effect = quantile(patient_High_vs_Other_effect, .25, na.rm = TRUE),
  q3_Yan_effect = quantile(patient_High_vs_Other_effect, .75, na.rm = TRUE),
  patient_unit_wilcoxon_P = safe_wilcox(patient_High_vs_Other_effect)
), by = .(gene, model)]
gene_summary[, BH_FDR := bh(patient_unit_wilcoxon_P), by = model]
gene_summary <- merge(gene_summary, program[, .(
  gene, Wu_overall_log2FC = overall_log2FC_High_vs_Other,
  Wu_BH_FDR = BH_FDR, Wu_positive_patient_count = positive_patient_count,
  representative_gene
)], by = "gene", all.x = TRUE)
gene_summary[, `:=`(
  record_type = "gene_summary", patient = "ALL_ELIGIBLE",
  Yan_conservation_class = fifelse(median_Yan_effect > 0, "direction_preserved",
                            fifelse(median_Yan_effect < 0, "direction_reversed", "near_zero"))
)]

missing_program <- data.table(
  record_type = "gene_summary", patient = "ALL_ELIGIBLE",
  gene = setdiff(program_genes, program_present), model = c("raw"),
  n_testable = 0L, n_positive = NA_integer_, n_negative = NA_integer_, n_zero = NA_integer_,
  median_Yan_effect = NA_real_, q1_Yan_effect = NA_real_, q3_Yan_effect = NA_real_,
  patient_unit_wilcoxon_P = NA_real_, BH_FDR = NA_real_, Yan_conservation_class = "unavailable"
)
missing_program <- merge(missing_program, program[, .(
  gene, Wu_overall_log2FC = overall_log2FC_High_vs_Other,
  Wu_BH_FDR = BH_FDR, Wu_positive_patient_count = positive_patient_count,
  representative_gene
)], by = "gene", all.x = TRUE)
write_tsv(rbindlist(list(gene_patient, gene_summary, missing_program), fill = TRUE), "Fig5_strictv2_gene_conservation.tsv")

# -----------------------------------------------------------------------------
# Candidate 10: three-cohort directional conservation, without gene reselection.
# -----------------------------------------------------------------------------
yan_raw_gene <- gene_summary[model == "raw", .(
  gene, Yan_n_testable = n_testable, Yan_n_positive = n_positive,
  Yan_median_effect = median_Yan_effect, Yan_direction = fifelse(median_Yan_effect > 0, "positive",
                                                        fifelse(median_Yan_effect < 0, "negative", "zero"))
)]
gse <- read_tsv(paths[["gse_transfer"]])
three <- merge(program[, .(
  gene, Wu_effect = overall_log2FC_High_vs_Other, Wu_BH_FDR = BH_FDR,
  Wu_positive_patient_count = positive_patient_count,
  representative_gene
)], yan_raw_gene, by = "gene", all.x = TRUE)
three <- merge(three, gse[, .(
  gene, detectable_in_GSE180286, GSE180286_effect = GSE180286_High_vs_Other_effect,
  GSE180286_direction = GSE180286_effect_direction,
  GSE180286_expected_direction_preserved = expected_direction_preserved,
  GSE180286_unavailability_reason = unavailability_reason
)], by = "gene", all.x = TRUE)
three[, `:=`(
  Wu_direction = "positive_by_frozen_program_definition",
  Yan_measurable = is.finite(Yan_median_effect),
  GSE180286_positive = is.finite(GSE180286_effect) & GSE180286_effect > 0,
  Yan_positive = is.finite(Yan_median_effect) & Yan_median_effect > 0,
  three_cohort_measurable = is.finite(Yan_median_effect) & is.finite(GSE180286_effect),
  direction_preserved_in_all_three = is.finite(Yan_median_effect) & is.finite(GSE180286_effect) & Yan_median_effect > 0 & GSE180286_effect > 0
)]
three[, conservation_class := fifelse(!three_cohort_measurable, "unavailable_in_one_or_more_validation_datasets",
                               fifelse(direction_preserved_in_all_three, "positive_in_all_three",
                               fifelse(Yan_positive & !GSE180286_positive, "positive_in_Yan_only",
                               fifelse(!Yan_positive & GSE180286_positive, "positive_in_GSE180286_only", "not_positive_in_either_validation_dataset"))))]
write_tsv(three, "Fig5_threecohort_conservation.tsv")

source_hash <- data.table(
  source_role = names(paths), source_file = unname(paths),
  sha256 = vapply(unname(paths), digest, character(1), algo = "sha256", file = TRUE)
)
write_tsv(source_hash, "Fig5_stage2_h5ad_source_hashes.tsv")

cat("Stage 2 H5AD-derived candidates completed.\n")
cat("Modules:", nrow(modules), "\n")
cat("Program genes detected:", length(program_present), "/146\n")
cat("Eligible patients:", uniqueN(d$patient), "\n")
