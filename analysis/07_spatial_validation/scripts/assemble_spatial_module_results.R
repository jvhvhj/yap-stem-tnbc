# Purpose: Assemble patient-level spatial module results
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 6 and Supplementary Fig. S4.
options(stringsAsFactors = FALSE)

out_root <- "Figure6_FG_biological_enrichment_phase"
chunk_root <- file.path(out_root, "intermediate", "module_score_chunks")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

registry <- read.delim(".runtime/Figure6_S1_inputs/section_registry.tsv", check.names = FALSE)
registry$patient_number <- as.integer(sub("^P", "", registry$patient_id))
registry <- registry[order(registry$patient_number, registry$section_number), ]
module_table <- read.delim("0716_yan2026_validation/tables/module_score_independent_gene_lists.tsv", check.names = FALSE)
identity_audit <- read.delim("0719_frozen_score_definition_audit/tables/module_gene_identity_checks.tsv", check.names = FALSE)
score_genes <- read.delim(".runtime/Figure6_S1_inputs/score_gene_list.tsv", check.names = FALSE)
program146 <- scan(".runtime/Figure6_S1_inputs/146_gene_program.txt", what = character(), quiet = TRUE)
composition <- read.delim("Figure6_S0Q3_GSE210616_ESTIMATE_tumour_composition_reconstruction/Figure6_S0Q3_ESTIMATE_spot_scores.tsv", check.names = FALSE)

module_order <- c("UPR", "TNFA NFKB", "Hypoxia", "Adhesion Remodeling", "Survival Stress")
module_cols <- c("UPR", "TNFA_NFKB", "Hypoxia", "Adhesion_Remodeling", "Survival_Stress")
names(module_cols) <- module_order
patient_order <- unique(registry$patient_id)
stopifnot(length(patient_order) == 22L)

hash_file <- function(path) unname(tools::md5sum(path))
sha256_file <- function(path) {
  cmd <- sprintf('certutil -hashfile "%s" SHA256', normalizePath(path, winslash = "\\", mustWork = TRUE))
  z <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
  z <- gsub(" ", "", z[grepl("^[0-9A-Fa-f ]{64,}$", z)])
  if (length(z)) toupper(z[1]) else "NOT_AVAILABLE"
}

alias_map <- c(CYR61 = "CCN1", CTGF = "CCN2", `HNRNPU-AS1` = "HNRNPU")
map_alias <- function(x) {
  y <- x; hit <- x %in% names(alias_map); y[hit] <- unname(alias_map[x[hit]]); y
}
yap17 <- unique(map_alias(score_genes$gene[score_genes$in_YAP_score]))
stem21 <- unique(map_alias(score_genes$gene[score_genes$in_Stemness_score]))
joint38 <- union(yap17, stem21)
program146_norm <- unique(map_alias(program146))

score_files <- file.path(".runtime/Figure6_S1_chunks",
                         paste0(registry$sample_id, "_", registry$section_id, "_scores.tsv"))
module_files <- file.path(chunk_root, paste0(registry$sample_id, "_", registry$section_id, "_module_scores.tsv"))
meta_files <- file.path(chunk_root, paste0(registry$sample_id, "_", registry$section_id, "_module_meta.tsv"))
error_files <- list.files(chunk_root, pattern = "_ERROR[.]txt$", full.names = TRUE)
if (length(error_files)) stop("Module-score worker errors exist: ", paste(error_files, collapse = "; "))
if (!all(file.exists(score_files)) || !all(file.exists(module_files)) || !all(file.exists(meta_files))) stop("Not all 43 required chunks exist")

module_spots <- do.call(rbind, lapply(module_files, function(f) read.delim(f, check.names = FALSE)))
module_meta <- do.call(rbind, lapply(meta_files, function(f) read.delim(f, check.names = FALSE)))
s1_spots <- do.call(rbind, lapply(score_files, function(f) read.delim(f, check.names = FALSE)))

key <- function(x) paste(x$sample_id, x$section_id, x$barcode, sep = "|")
km <- key(module_spots); ks <- key(s1_spots); kc <- key(composition)
stopifnot(nrow(module_spots) == 56567L, nrow(s1_spots) == 56567L, nrow(composition) == 56567L)
stopifnot(!anyDuplicated(km), !anyDuplicated(ks), !anyDuplicated(kc))
mi_s <- match(km, ks); mi_c <- match(km, kc)
stopifnot(!anyNA(mi_s), !anyNA(mi_c))
dat <- cbind(
  module_spots,
  s1_spots[mi_s, c("Joint_spatial", "Program146_spatial"), drop = FALSE],
  composition[mi_c, c("TumorPurity", "StromalScore", "ImmuneScore"), drop = FALSE])
stopifnot(all(dat$nCount_raw == s1_spots$nCount_raw[mi_s]), all(dat$nFeature_raw == s1_spots$nFeature_raw[mi_s]))

rank_residual <- function(y, covars) {
  yr <- rank(y, ties.method = "average")
  xr <- as.data.frame(lapply(covars, rank, ties.method = "average"))
  residuals(lm(yr ~ ., data = xr))
}
partial_spearman_equiv <- function(x, y, covars) {
  ok <- is.finite(x) & is.finite(y) & Reduce(`&`, lapply(covars, is.finite))
  x <- x[ok]; y <- y[ok]; covars <- lapply(covars, `[`, ok)
  if (length(x) < 4L || sd(x) == 0 || sd(y) == 0) return(NA_real_)
  cor(rank_residual(x, covars), rank_residual(y, covars), method = "pearson")
}
safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 4L || sd(x) == 0 || sd(y) == 0) return(NA_real_)
  cor(x, y, method = "spearman")
}
aggregate_fisher <- function(rho, n) {
  ok <- is.finite(rho) & is.finite(n) & n > 3
  if (!any(ok)) return(NA_real_)
  rr <- pmax(pmin(rho[ok], 1 - 1e-12), -1 + 1e-12)
  tanh(weighted.mean(atanh(rr), w = n[ok] - 3))
}
summ <- function(x) {
  c(median = median(x), Q1 = unname(quantile(x, .25, type = 7)),
    Q3 = unname(quantile(x, .75, type = 7)), IQR = IQR(x),
    minimum = min(x), maximum = max(x), positive_n = sum(x > 0), tested_n = length(x))
}
wilcox_exact <- function(x, alternative = "greater") {
  x <- x[is.finite(x)]
  w <- suppressWarnings(wilcox.test(x, mu = 0, alternative = alternative, exact = TRUE, correct = FALSE))
  c(statistic = unname(w$statistic), P = w$p.value)
}
write_tsv <- function(x, name) write.table(x, file.path(out_root, name), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")

# Reproduce the accepted S1 residualization numerically before extending it.
s1_reference <- read.delim("Figure6_S1_GSE210616_patient_recurrent_spatial_coupling/Figure6_S1_section_associations.tsv", check.names = FALSE)
resid_check <- lapply(seq_len(nrow(registry)), function(i) {
  sr <- registry[i, ]
  d <- dat[dat$sample_id == sr$sample_id & dat$section_id == sr$section_id, ]
  ref <- s1_reference[s1_reference$sample_id == sr$sample_id & s1_reference$section_id == sr$section_id, ]
  data.frame(
    depth_difference = partial_spearman_equiv(d$Joint_spatial, d$Program146_spatial, list(d$TumorPurity, log1p(d$nCount_raw))) - ref$rho_tumour_depth_adjusted,
    nfeature_difference = partial_spearman_equiv(d$Joint_spatial, d$Program146_spatial, list(d$TumorPurity, d$nFeature_raw)) - ref$rho_tumour_nfeature_adjusted)
})
resid_check <- do.call(rbind, resid_check)
residualization_max_abs_difference <- max(abs(unlist(resid_check)))
stopifnot(is.finite(residualization_max_abs_difference), residualization_max_abs_difference < 1e-12)

# -----------------------------------------------------------------------------
# Functional-program associations
# -----------------------------------------------------------------------------
module_src <- "0716_yan2026_validation/tables/module_score_independent_gene_lists.tsv"
module_hash <- sha256_file(module_src)
prov_rows <- list()
for (m in module_order) {
  genes <- unique(map_alias(module_table$gene[module_table$module == m]))
  for (g in genes) {
    prov_rows[[length(prov_rows) + 1L]] <- data.frame(
      module = m, gene = g, frozen_source_version = unique(module_table$source_version[module_table$module == m])[1],
      primary_spatial_score_membership = TRUE, source_file = normalizePath(module_src, winslash = "/"),
      source_sha256 = module_hash, scoring_rule = "within-section z-standardize each measurable gene across in-tissue spots; mean gene z",
      stringsAsFactors = FALSE)
  }
}
f_prov <- do.call(rbind, prov_rows)
f_prov$identity_audit_file <- normalizePath("0719_frozen_score_definition_audit/tables/module_gene_identity_checks.tsv", winslash = "/")
f_prov$identity_audit_status <- identity_audit$status[match(f_prov$module, identity_audit$module)]
f_prov$identity_audit_sha256 <- sha256_file("0719_frozen_score_definition_audit/tables/module_gene_identity_checks.tsv")
stopifnot(all(f_prov$identity_audit_status == "PASS_EXACT"))
write_tsv(f_prov, "Figure6F_module_provenance.tsv")

coverage_rows <- list()
for (m in module_order) {
  genes <- f_prov$gene[f_prov$module == m]
  mm <- module_meta[module_meta$object == m, ]
  coverage_rows[[length(coverage_rows) + 1L]] <- data.frame(
    module = m, frozen_gene_n = length(genes), section_n = nrow(mm),
    n_measurable_GSE210616 = min(mm$measurable_n), coverage_fraction = min(mm$measurable_n) / length(genes),
    measurable_gene_n_min = min(mm$measurable_n), measurable_gene_n_median = median(mm$measurable_n), measurable_gene_n_max = max(mm$measurable_n),
    full_gene_coverage_section_n = sum(mm$measurable_n == length(genes)),
    YAP17_overlap_n = length(intersect(genes, yap17)), YAP17_overlap_genes = paste(intersect(genes, yap17), collapse = ";"),
    Stem21_overlap_n = length(intersect(genes, stem21)), Stem21_overlap_genes = paste(intersect(genes, stem21), collapse = ";"),
    Joint38_overlap_n = length(intersect(genes, joint38)), Joint38_overlap_genes = paste(intersect(genes, joint38), collapse = ";"),
    Program146_overlap_n = length(intersect(genes, program146_norm)), Program146_overlap_genes = paste(intersect(genes, program146_norm), collapse = ";"),
    primary_overlap_gate = if (length(intersect(genes, joint38)) == 0L) "PASS_ZERO_JOINT38_OVERLAP" else "FAIL",
    stringsAsFactors = FALSE)
}
f_coverage <- do.call(rbind, coverage_rows)
for (nm in grep("_overlap_genes$", names(f_coverage), value = TRUE)) f_coverage[[nm]][f_coverage[[nm]] == ""] <- "NONE"
write_tsv(f_coverage, "Figure6F_module_coverage_overlap.tsv")
if (any(f_coverage$primary_overlap_gate != "PASS_ZERO_JOINT38_OVERLAP")) stop("F overlap gate failed")

f_section_rows <- list()
for (i in seq_len(nrow(registry))) {
  sr <- registry[i, ]
  d <- dat[dat$sample_id == sr$sample_id & dat$section_id == sr$section_id, ]
  for (m in module_order) {
    y <- d[[module_cols[[m]]]]
    f_section_rows[[length(f_section_rows) + 1L]] <- data.frame(
      patient_id = sr$patient_id, sample_id = sr$sample_id, section_id = sr$section_id,
      section_number = sr$section_number, n_spots = nrow(d), module = m,
      rho_raw = safe_spearman(d$Joint_spatial, y),
      rho_tumour_depth_adjusted = partial_spearman_equiv(d$Joint_spatial, y, list(d$TumorPurity, log1p(d$nCount_raw))),
      rho_tumour_nfeature_adjusted = partial_spearman_equiv(d$Joint_spatial, y, list(d$TumorPurity, d$nFeature_raw)),
      inference_role = "SECTION_DESCRIPTIVE_ONLY", stringsAsFactors = FALSE)
  }
}
f_section <- do.call(rbind, f_section_rows)
write_tsv(f_section, "Figure6F_section_module_associations.tsv")

f_patient_rows <- list()
for (pid in patient_order) for (m in module_order) {
  z <- f_section[f_section$patient_id == pid & f_section$module == m, ]
  f_patient_rows[[length(f_patient_rows) + 1L]] <- data.frame(
    patient_id = pid, n_sections = nrow(z), total_spots = sum(z$n_spots), weight_sum = sum(z$n_spots - 3), module = m,
    rho_raw = aggregate_fisher(z$rho_raw, z$n_spots),
    rho_tumour_depth_adjusted = aggregate_fisher(z$rho_tumour_depth_adjusted, z$n_spots),
    rho_tumour_nfeature_adjusted = aggregate_fisher(z$rho_tumour_nfeature_adjusted, z$n_spots),
    biological_replication_unit = "PATIENT", stringsAsFactors = FALSE)
}
f_patient <- do.call(rbind, f_patient_rows)
write_tsv(f_patient, "Figure6F_patient_module_associations.tsv")

f_tests <- list()
for (m in module_order) {
  x <- f_patient$rho_tumour_depth_adjusted[f_patient$module == m]
  z <- summ(x); wt <- wilcox_exact(x, "greater")
  f_tests[[length(f_tests) + 1L]] <- data.frame(module = m, n_patients = z["tested_n"], median_rho = z["median"], Q1 = z["Q1"], Q3 = z["Q3"], IQR = z["IQR"], minimum = z["minimum"], maximum = z["maximum"], positive_n = z["positive_n"], positive_fraction = z["positive_n"] / z["tested_n"], wilcoxon_statistic = wt["statistic"], one_sided_exact_P = wt["P"], stringsAsFactors = FALSE)
}
f_tests <- do.call(rbind, f_tests)
f_tests$Holm_P <- p.adjust(f_tests$one_sided_exact_P, method = "holm")
f_tests$decision <- ifelse(f_tests$median_rho > 0 & f_tests$positive_n >= 18 & f_tests$Holm_P < .01, "STRONG_PASS", ifelse(f_tests$median_rho > 0 & f_tests$positive_n >= 16 & f_tests$Holm_P < .05, "PASS", "NOT_PASS"))
write_tsv(f_tests, "Figure6F_module_primary_tests.tsv")

f_sens <- list()
for (model in c("rho_raw", "rho_tumour_nfeature_adjusted")) for (m in module_order) {
  x <- f_patient[f_patient$module == m, model]
  z <- summ(x); wt <- wilcox_exact(x, "greater")
  f_sens[[length(f_sens) + 1L]] <- data.frame(model = model, module = m, n_patients = z["tested_n"], median_rho = z["median"], Q1 = z["Q1"], Q3 = z["Q3"], IQR = z["IQR"], minimum = z["minimum"], maximum = z["maximum"], positive_n = z["positive_n"], positive_fraction = z["positive_n"] / z["tested_n"], one_sided_exact_P = wt["P"], stringsAsFactors = FALSE)
}
f_sens <- do.call(rbind, f_sens)
f_sens$Holm_P_within_model <- ave(f_sens$one_sided_exact_P, f_sens$model, FUN = function(x) p.adjust(x, method = "holm"))
write_tsv(f_sens, "Figure6F_module_sensitivity.tsv")

f_pass_n <- sum(f_tests$decision %in% c("PASS", "STRONG_PASS"))
f_strong_n <- sum(f_tests$decision == "STRONG_PASS")
f_opposite <- sum(f_tests$median_rho < 0)
f_verdict <- if (f_pass_n >= 4 && f_strong_n >= 3 && f_opposite == 0) "FIGURE6F_MAIN_PANEL_STRONG" else if (f_pass_n >= 4) "FIGURE6F_MAIN_PANEL_PASS" else "FIGURE6F_EXTENDED_DATA_ONLY"

# -----------------------------------------------------------------------------
# Remodeling-program associations
# -----------------------------------------------------------------------------
g_prov <- data.frame(
  object = c("Joint_spatial", "Program146_spatial", "TumorPurity", "StromalScore", "ImmuneScore", "Joint_spatial_no_ZEB2", "Program146_spatial_no_DDR2_TNFAIP3"),
  role = c("frozen primary score", "frozen primary score", "primary tissue-composition context", "secondary tissue context", "secondary tissue context", "overlap-clean sensitivity", "overlap-clean sensitivity"),
  definition = c("z_section(YAP17 mean normalized expression)+z_section(Stem21 mean normalized expression)", "mean within-section gene z across frozen Program146", "ESTIMATE-derived continuous TumorPurity", "ESTIMATE StromalScore", "ESTIMATE ImmuneScore", "frozen Joint score with ZEB2 excluded from Stem21 only", "frozen Program146 score with DDR2 and TNFAIP3 excluded only"),
  source_file = c(rep(".runtime/Figure6_S1_chunks/*_scores.tsv", 2), rep(normalizePath("Figure6_S0Q3_GSE210616_ESTIMATE_tumour_composition_reconstruction/Figure6_S0Q3_ESTIMATE_spot_scores.tsv", winslash = "/"), 3), rep(normalizePath(chunk_root, winslash = "/"), 2)),
  adjustment = c("log1p(nCount) primary; nFeature sensitivity", "log1p(nCount) primary; nFeature sensitivity", rep("outcome/context variable; rank-residualized with score against technical covariate", 3), rep("log1p(nCount) primary", 2)),
  stringsAsFactors = FALSE)
g_prov$residualization_validation <- sprintf("Reproduced all 43 accepted S1 section estimates; max absolute difference %.3g", residualization_max_abs_difference)
write_tsv(g_prov, "Figure6G_composition_provenance.tsv")

score_map <- c(Joint = "Joint_spatial", Program146 = "Program146_spatial")
clean_map <- c(Joint = "Joint_spatial_no_ZEB2", Program146 = "Program146_spatial_no_DDR2_TNFAIP3")
g_section_rows <- list()
for (i in seq_len(nrow(registry))) {
  sr <- registry[i, ]; d <- dat[dat$sample_id == sr$sample_id & dat$section_id == sr$section_id, ]
  for (sn in names(score_map)) {
    x <- d[[score_map[[sn]]]]
    for (metric in c("TumorPurity", "StromalScore", "ImmuneScore")) {
      y <- d[[metric]]
      g_section_rows[[length(g_section_rows) + 1L]] <- data.frame(
        patient_id = sr$patient_id, sample_id = sr$sample_id, section_id = sr$section_id, section_number = sr$section_number,
        n_spots = nrow(d), score = sn, context_metric = metric,
        rho_raw = safe_spearman(x, y),
        rho_log1p_nCount_adjusted = partial_spearman_equiv(x, y, list(log1p(d$nCount_raw))),
        rho_nFeature_adjusted = partial_spearman_equiv(x, y, list(d$nFeature_raw)),
        inference_role = "SECTION_DESCRIPTIVE_ONLY", stringsAsFactors = FALSE)
    }
  }
}
g_section <- do.call(rbind, g_section_rows)
write_tsv(g_section, "Figure6G_section_composition_associations.tsv")

g_patient_rows <- list()
for (pid in patient_order) for (sn in names(score_map)) for (metric in c("TumorPurity", "StromalScore", "ImmuneScore")) {
  z <- g_section[g_section$patient_id == pid & g_section$score == sn & g_section$context_metric == metric, ]
  g_patient_rows[[length(g_patient_rows) + 1L]] <- data.frame(
    patient_id = pid, n_sections = nrow(z), total_spots = sum(z$n_spots), weight_sum = sum(z$n_spots - 3), score = sn, context_metric = metric,
    rho_raw = aggregate_fisher(z$rho_raw, z$n_spots),
    rho_log1p_nCount_adjusted = aggregate_fisher(z$rho_log1p_nCount_adjusted, z$n_spots),
    rho_nFeature_adjusted = aggregate_fisher(z$rho_nFeature_adjusted, z$n_spots),
    biological_replication_unit = "PATIENT", stringsAsFactors = FALSE)
}
g_patient <- do.call(rbind, g_patient_rows)
write_tsv(g_patient, "Figure6G_patient_composition_associations.tsv")

g_primary <- list()
for (sn in names(score_map)) {
  z0 <- g_patient[g_patient$score == sn & g_patient$context_metric == "TumorPurity", ]
  x <- z0$rho_log1p_nCount_adjusted; z <- summ(x); wt <- wilcox_exact(x, "greater")
  g_primary[[length(g_primary) + 1L]] <- data.frame(score = sn, context_metric = "TumorPurity", model = "log1p_nCount_adjusted", n_patients = z["tested_n"], median_rho = z["median"], Q1 = z["Q1"], Q3 = z["Q3"], IQR = z["IQR"], minimum = z["minimum"], maximum = z["maximum"], positive_n = z["positive_n"], positive_fraction = z["positive_n"] / z["tested_n"], one_sided_exact_P = wt["P"], stringsAsFactors = FALSE)
}
g_primary <- do.call(rbind, g_primary)
g_primary$Holm_P <- p.adjust(g_primary$one_sided_exact_P, method = "holm")

g_nfeature <- list()
for (sn in names(score_map)) {
  z0 <- g_patient[g_patient$score == sn & g_patient$context_metric == "TumorPurity", ]
  x <- z0$rho_nFeature_adjusted; z <- summ(x); wt <- wilcox_exact(x, "greater")
  g_nfeature[[length(g_nfeature) + 1L]] <- data.frame(
    score = sn, nfeature_median_rho = z["median"], nfeature_Q1 = z["Q1"], nfeature_Q3 = z["Q3"],
    nfeature_IQR = z["IQR"], nfeature_positive_n = z["positive_n"],
    nfeature_one_sided_exact_P = wt["P"], stringsAsFactors = FALSE)
}
g_nfeature <- do.call(rbind, g_nfeature)
g_nfeature$nfeature_Holm_P <- p.adjust(g_nfeature$nfeature_one_sided_exact_P, method = "holm")
g_primary <- merge(g_primary, g_nfeature, by = "score", sort = FALSE)

g_secondary <- list()
for (sn in names(score_map)) for (metric in c("StromalScore", "ImmuneScore")) {
  z0 <- g_patient[g_patient$score == sn & g_patient$context_metric == metric, ]
  x <- z0$rho_log1p_nCount_adjusted; z <- summ(x); wt <- wilcox_exact(x, "two.sided")
  g_secondary[[length(g_secondary) + 1L]] <- data.frame(score = sn, context_metric = metric, model = "log1p_nCount_adjusted", n_patients = z["tested_n"], median_rho = z["median"], Q1 = z["Q1"], Q3 = z["Q3"], IQR = z["IQR"], minimum = z["minimum"], maximum = z["maximum"], positive_n = z["positive_n"], positive_fraction = z["positive_n"] / z["tested_n"], two_sided_exact_P = wt["P"], stringsAsFactors = FALSE)
}
g_secondary <- do.call(rbind, g_secondary)
g_secondary$Holm_P <- p.adjust(g_secondary$two_sided_exact_P, method = "holm")
write_tsv(g_secondary, "Figure6G_secondary_context_tests.tsv")

clean_section <- list()
for (i in seq_len(nrow(registry))) {
  sr <- registry[i, ]; d <- dat[dat$sample_id == sr$sample_id & dat$section_id == sr$section_id, ]
  for (sn in names(clean_map)) {
    clean_section[[length(clean_section) + 1L]] <- data.frame(
      patient_id = sr$patient_id, sample_id = sr$sample_id, section_id = sr$section_id, n_spots = nrow(d), score = sn,
      rho_tumour_log1p_nCount_adjusted = partial_spearman_equiv(d[[clean_map[[sn]]]], d$TumorPurity, list(log1p(d$nCount_raw))),
      stringsAsFactors = FALSE)
  }
}
clean_section <- do.call(rbind, clean_section)
clean_section_out <- transform(clean_section, record_type = "SECTION", n_sections = NA_integer_, total_spots = NA_integer_, median_rho = NA_real_, positive_n = NA_integer_, n_patients = NA_integer_, one_sided_exact_P = NA_real_, Holm_P = NA_real_)
clean_section_out <- clean_section_out[, c("record_type", "patient_id", "sample_id", "section_id", "score", "n_spots", "n_sections", "total_spots", "rho_tumour_log1p_nCount_adjusted", "median_rho", "positive_n", "n_patients", "one_sided_exact_P", "Holm_P")]
clean_patient <- list()
for (pid in patient_order) for (sn in names(clean_map)) {
  z <- clean_section[clean_section$patient_id == pid & clean_section$score == sn, ]
  clean_patient[[length(clean_patient) + 1L]] <- data.frame(record_type = "PATIENT", patient_id = pid, sample_id = NA_character_, section_id = NA_character_, score = sn, n_spots = NA_integer_, n_sections = nrow(z), total_spots = sum(z$n_spots), rho_tumour_log1p_nCount_adjusted = aggregate_fisher(z$rho_tumour_log1p_nCount_adjusted, z$n_spots), median_rho = NA_real_, positive_n = NA_integer_, n_patients = NA_integer_, one_sided_exact_P = NA_real_, Holm_P = NA_real_, stringsAsFactors = FALSE)
}
clean_patient <- do.call(rbind, clean_patient)
clean_summary <- list()
for (sn in names(clean_map)) {
  x <- clean_patient$rho_tumour_log1p_nCount_adjusted[clean_patient$score == sn]
  z <- summ(x); wt <- wilcox_exact(x, "greater")
  clean_summary[[length(clean_summary) + 1L]] <- data.frame(record_type = "SUMMARY", patient_id = "ALL", sample_id = NA_character_, section_id = NA_character_, score = sn, n_spots = NA_integer_, n_sections = NA_integer_, total_spots = NA_integer_, rho_tumour_log1p_nCount_adjusted = NA_real_, median_rho = z["median"], positive_n = z["positive_n"], n_patients = z["tested_n"], one_sided_exact_P = wt["P"], Holm_P = NA_real_, stringsAsFactors = FALSE)
}
clean_summary <- do.call(rbind, clean_summary)
clean_summary$Holm_P <- p.adjust(clean_summary$one_sided_exact_P, method = "holm")
g_clean <- rbind(clean_section_out, clean_patient, clean_summary)
write_tsv(g_clean, "Figure6G_overlap_clean_sensitivity.tsv")

clean_positive <- setNames(clean_summary$median_rho > 0 & clean_summary$positive_n >= 16, clean_summary$score)
g_primary$overlap_clean_sensitivity_positive <- clean_positive[g_primary$score]
g_primary$decision <- ifelse(g_primary$median_rho > 0 & g_primary$positive_n >= 18 & g_primary$Holm_P < .01 & g_primary$overlap_clean_sensitivity_positive, "STRONG_PASS", ifelse(g_primary$median_rho > 0 & g_primary$positive_n >= 16 & g_primary$Holm_P < .05, "PASS", "NOT_PASS"))
write_tsv(g_primary, "Figure6G_primary_tumour_context_tests.tsv")

g_strong_n <- sum(g_primary$decision == "STRONG_PASS")
g_pass_n <- sum(g_primary$decision %in% c("PASS", "STRONG_PASS"))
g_verdict <- if (g_strong_n == 2L) "FIGURE6G_MAIN_PANEL_STRONG" else if (g_pass_n == 2L) "FIGURE6G_MAIN_PANEL_PASS" else if (g_pass_n == 1L) "FIGURE6G_CONTEXT_ONLY" else "FIGURE6G_EXTENDED_DATA_ONLY"

f_lines <- c(
  "# Figure 6F — spatial recurrence of five core functional programmes", "",
  "## Frozen inputs", sprintf("Five score-independent module definitions were recovered from `%s` (SHA-256 `%s`).", module_src, module_hash),
  "All five primary spatial module scores have zero overlap with the frozen Joint38 score-definition universe. No AddModuleScore, control genes, feature selection, smoothing or spot thresholding was used.", "",
  "## Statistical unit and model", sprintf("All 56,567 in-tissue spots from 43 sections were retained. Within each section, every measurable module gene was z-standardized across spots and module score was the mean gene z. Raw Spearman and rank-residual partial-Spearman-equivalent associations were calculated. The primary model residualized both Joint and module ranks on TumorPurity and log1p(nCount). The implementation reproduced all 43 accepted S1 section estimates with maximum absolute difference %.3g. Section effects were aggregated to 22 patients using Fisher-z weights n_spots-3. One-sided exact Wilcoxon signed-rank tests used patient as the biological unit; five primary P values were Holm-adjusted.", residualization_max_abs_difference), "",
  "## Primary results", paste(apply(f_tests, 1, function(r) sprintf("- %s: median rho %.4f (IQR %.4f), positive %d/%d, exact P %.4g, Holm P %.4g, %s.", r[["module"]], as.numeric(r[["median_rho"]]), as.numeric(r[["IQR"]]), as.integer(r[["positive_n"]]), as.integer(r[["n_patients"]]), as.numeric(r[["one_sided_exact_P"]]), as.numeric(r[["Holm_P"]]), r[["decision"]])), collapse = "\n"), "",
  "## Prespecified evidence assessment", f_verdict, "", "Interpretation is limited to the reported patient-level evidence.")
writeLines(f_lines, file.path(out_root, "Figure6F_readout.md"), useBytes = TRUE)

g_lines <- c(
  "# Figure 6G — tissue-composition context of frozen spatial scores", "",
  "## Frozen inputs", "Primary scores were the frozen Joint_spatial and Program146_spatial values from S1. TumorPurity, StromalScore and ImmuneScore were the authoritative continuous ESTIMATE values for all 56,567 spots. No tumour threshold was created.", "",
  "## Statistical unit and model", sprintf("Within-section associations used the exact S1 partial-Spearman-equivalent implementation: Pearson correlation of score-rank and context-rank residuals after rank-residualization on the specified technical covariate. Reproduction against all 43 accepted S1 section estimates had maximum absolute difference %.3g. The primary TumorPurity model adjusted for log1p(nCount); nFeature was sensitivity. Section effects were Fisher-z aggregated to 22 patients with n_spots-3 weights. Primary one-sided tests were Holm-adjusted across two frozen scores; secondary two-sided tests were Holm-adjusted across four score-context pairs.", residualization_max_abs_difference), "",
  "## Primary TumorPurity results", paste(apply(g_primary, 1, function(r) sprintf("- %s: median rho %.4f (IQR %.4f), positive %d/%d, exact P %.4g, Holm P %.4g, overlap-clean positive=%s, %s.", r[["score"]], as.numeric(r[["median_rho"]]), as.numeric(r[["IQR"]]), as.integer(r[["positive_n"]]), as.integer(r[["n_patients"]]), as.numeric(r[["one_sided_exact_P"]]), as.numeric(r[["Holm_P"]]), r[["overlap_clean_sensitivity_positive"]], r[["decision"]])), collapse = "\n"), "",
  "## nFeature sensitivity", paste(apply(g_primary, 1, function(r) sprintf("- %s: median rho %.4f, positive %d/22, exact P %.4g, Holm P %.4g.", r[["score"]], as.numeric(r[["nfeature_median_rho"]]), as.integer(r[["nfeature_positive_n"]]), as.numeric(r[["nfeature_one_sided_exact_P"]]), as.numeric(r[["nfeature_Holm_P"]]))), collapse = "\n"), "",
  "## Secondary stromal/immune context", paste(apply(g_secondary, 1, function(r) sprintf("- %s vs %s: median rho %.4f, positive %d/22, two-sided exact P %.4g, Holm P %.4g.", r[["score"]], r[["context_metric"]], as.numeric(r[["median_rho"]]), as.integer(r[["positive_n"]]), as.numeric(r[["two_sided_exact_P"]]), as.numeric(r[["Holm_P"]]))), collapse = "\n"), "",
  "## Boundary", "TumorPurity is an ESTIMATE-derived expression context, not a malignant-cell fraction or histologic tumour mask. Stromal and immune analyses are secondary contextual associations and do not establish cellular origin.", "",
  "## Prespecified evidence assessment", g_verdict, "", "Interpretation is limited to the reported patient-level evidence.")
writeLines(g_lines, file.path(out_root, "Figure6G_readout.md"), useBytes = TRUE)

cat(sprintf("F=%s G=%s rows=%d spots=%d\n", f_verdict, g_verdict, nrow(dat), nrow(dat)))
