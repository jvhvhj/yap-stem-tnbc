#!/usr/bin/env Rscript
# Purpose: Patient-level scoring and covariate robustness
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 3 and Supplementary Fig. S2.

options(stringsAsFactors = FALSE, future.globals.maxSize = 8 * 1024^3)
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
})
future::plan("sequential")

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else getwd()
default_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
root <- Sys.getenv("AHIPPO_YAP_ROOT", unset = default_root)
out <- file.path(root, "0717_methodological_framework_benchmark")
dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "logs"), recursive = TRUE, showWarnings = FALSE)

log_con <- file(file.path(out, "logs", "02_wu_framework_effects.log"), open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message", append = TRUE)
on.exit({sink(type = "message"); sink(type = "output"); close(log_con)}, add = TRUE)

write_tsv <- function(x, name) {
  path <- file.path(out, "tables", name)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  cat(sprintf("Wrote %s (%d rows)\n", path, nrow(x)))
}

safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 5L || sd(x[keep]) == 0 || sd(y[keep]) == 0) return(NA_real_)
  suppressWarnings(cor(x[keep], y[keep], method = "spearman"))
}

z <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

residualize <- function(y, covariates) {
  covariates <- as.matrix(covariates)
  keep <- apply(covariates, 2, function(v) sd(v, na.rm = TRUE) > 1e-12)
  covariates <- covariates[, keep, drop = FALSE]
  if (ncol(covariates)) covariates <- apply(covariates, 2, z)
  design <- cbind(Intercept = 1, covariates)
  as.numeric(y - qr.Q(qr(design)) %*% crossprod(qr.Q(qr(design)), y))
}

feature_score <- function(genes, expr) {
  present <- intersect(genes, rownames(expr))
  if (!length(present)) return(rep(NA_real_, ncol(expr)))
  Matrix::colMeans(expr[present, , drop = FALSE])
}

cat("Loading frozen Wu compact object\n")
obj <- readRDS(file.path(root, "output_step0", "Wu2021_TNBC_malignant.rds"))
expr <- GetAssayData(obj, assay = "RNA", layer = "data")
meta <- read.csv(file.path(root, "0703_rebuild", "continuous_state_main_figures", "Wu2021_continuous_state_cell_metadata.csv"), check.names = FALSE)
meta <- meta[match(colnames(expr), meta$cell_id), , drop = FALSE]
stopifnot(all(meta$cell_id == colnames(expr)))

manifest <- read.delim(file.path(root, "0716_yan2026_validation", "tables", "frozen_signature_manifest.tsv"), check.names = FALSE)
modules <- c("UPR", "TNFA NFKB", "Hypoxia", "Adhesion Remodeling", "Wound Healing", "Survival Stress", "Anoikis Resistance", "Integrin Adhesion")
module_manifest <- manifest %>% filter(signature_type == "module", module %in% modules)
module_sets <- list()
for (m in modules) {
  d <- module_manifest %>% filter(module == m)
  module_sets[[paste(m, "original", sep = "||")]] <- unique(d$gene)
  module_sets[[paste(m, "score_independent", sep = "||")]] <- unique(d$gene[d$score_independent_gene %in% c(TRUE, "True", "TRUE")])
  module_sets[[paste(m, "pairwise_unique", sep = "||")]] <- unique(d$gene[d$pairwise_unique_gene %in% c(TRUE, "True", "TRUE")])
}
rep_manifest <- manifest %>% filter(signature_type == "representative_gene") %>% distinct(gene, module, in_yap_stem_score)
rep_genes <- intersect(rep_manifest$gene, rownames(expr))

score_rows <- list()
for (key in names(module_sets)) {
  parts <- strsplit(key, "\\|\\|", fixed = FALSE)[[1]]
  score_rows[[key]] <- data.frame(
    feature_type = "module", feature = parts[1], version = parts[2],
    n_genes = length(module_sets[[key]]),
    score = I(list(feature_score(module_sets[[key]], expr)))
  )
}
feature_defs <- bind_rows(lapply(score_rows, function(x) x[, 1:4]))
write_tsv(feature_defs, "wu_framework_feature_gene_counts.tsv")

# Materialize only the small frozen feature matrices.
feature_vectors <- lapply(module_sets, feature_score, expr = expr)
for (g in rep_genes) feature_vectors[[paste(g, "representative_gene", sep = "||")]] <- as.numeric(expr[g, ])

global_group <- ifelse(meta$YAP_Stem_Score >= median(meta$YAP_Stem_Score), "High", "Low")
hh_group <- ifelse(
  meta$YAP_score >= median(meta$YAP_score) & meta$Stemness_score >= median(meta$Stemness_score),
  "HH", "Other"
)

effect_rows <- list(); k <- 1L
add_conventional <- function(feature_type, feature, version, values, group, workflow) {
  high_label <- if (workflow == "conventional_global_median") "High" else "HH"
  low_label <- if (workflow == "conventional_global_median") "Low" else "Other"
  keep <- is.finite(values) & group %in% c(high_label, low_label)
  high <- values[keep & group == high_label]
  low <- values[keep & group == low_label]
  p <- if (length(high) && length(low)) wilcox.test(high, low, exact = FALSE)$p.value else NA_real_
  data.frame(
    record_type = "pooled", cohort = "Wu2021", feature_type = feature_type,
    feature = feature, workflow = workflow, version = version, patient = "POOLED",
    n_cells = sum(keep), n_high = length(high), n_low = length(low),
    effect = median(high, na.rm = TRUE) - median(low, na.rm = TRUE),
    p_value = p, n_patients = NA_integer_, n_positive = NA_integer_, BH_FDR = NA_real_
  )
}

for (key in names(feature_vectors)) {
  parts <- strsplit(key, "\\|\\|", fixed = FALSE)[[1]]
  if (length(parts) == 2 && parts[2] == "representative_gene") {
    feature <- parts[1]; feature_type <- "representative_gene"; version <- "gene"
  } else {
    feature <- parts[1]; feature_type <- "module"; version <- parts[2]
  }
  values <- feature_vectors[[key]]
  if (all(!is.finite(values))) next
  effect_rows[[k]] <- add_conventional(feature_type, feature, version, values, global_group, "conventional_global_median"); k <- k + 1L
  effect_rows[[k]] <- add_conventional(feature_type, feature, version, values, hh_group, "conventional_HH_other"); k <- k + 1L

  # Patient-aware raw and technical-adjusted high-low effects.
  for (workflow in c("patient_raw", "technical_adjusted")) {
    patient_rows <- list()
    for (p in sort(unique(meta$patient))) {
      idx <- which(meta$patient == p)
      state <- meta$YS_tertile[idx]
      val <- values[idx]
      if (workflow == "technical_adjusted") {
        covs <- cbind(
          log_nCount = log1p(meta$nCount_RNA[idx]),
          nFeature = meta$nFeature_RNA[idx],
          Cell_cycle = meta$Cell_cycle[idx]
        )
        val <- residualize(val, covs)
      }
      high <- val[state == "High" & is.finite(val)]
      low <- val[state == "Low" & is.finite(val)]
      patient_rows[[p]] <- data.frame(
        record_type = "patient", cohort = "Wu2021", feature_type = feature_type,
        feature = feature, workflow = workflow, version = version, patient = p,
        n_cells = length(idx), n_high = length(high), n_low = length(low),
        effect = if (length(high) && length(low)) median(high) - median(low) else NA_real_,
        p_value = NA_real_, n_patients = NA_integer_, n_positive = NA_integer_, BH_FDR = NA_real_
      )
    }
    patient_df <- bind_rows(patient_rows)
    vals <- patient_df$effect[is.finite(patient_df$effect)]
    pval <- if (length(vals)) wilcox.test(vals, mu = 0, exact = FALSE)$p.value else NA_real_
    summary <- data.frame(
      record_type = "summary", cohort = "Wu2021", feature_type = feature_type,
      feature = feature, workflow = workflow, version = version, patient = "ALL",
      n_cells = sum(patient_df$n_cells), n_high = sum(patient_df$n_high), n_low = sum(patient_df$n_low),
      effect = if (length(vals)) median(vals) else NA_real_, p_value = pval,
      n_patients = length(vals), n_positive = sum(vals > 0), BH_FDR = NA_real_
    )
    effect_rows[[k]] <- patient_df; k <- k + 1L
    effect_rows[[k]] <- summary; k <- k + 1L
  }
}
effects <- bind_rows(effect_rows)
effects <- effects %>% group_by(workflow, feature_type, version, record_type) %>%
  mutate(BH_FDR = ifelse(record_type == "summary", p.adjust(p_value, method = "BH"), BH_FDR)) %>% ungroup()
write_tsv(effects, "framework_feature_effects_wu.tsv")

# Continuous score-module correlations on patient-aware versions, including
# technical adjustment of YAP, Stemness and the tested feature.
cor_rows <- list(); kc <- 1L
for (method in c("score_independent", "pairwise_unique")) {
  for (m in modules) {
    values <- feature_vectors[[paste(m, method, sep = "||")]]
    if (all(!is.finite(values))) next
    for (p in sort(unique(meta$patient))) {
      idx <- which(meta$patient == p)
      covs <- cbind(log_nCount = log1p(meta$nCount_RNA[idx]), nFeature = meta$nFeature_RNA[idx], Cell_cycle = meta$Cell_cycle[idx])
      for (model in c("raw", "technical_adjusted")) {
        y <- meta$YAP_score[idx]; s <- meta$Stemness_score[idx]; v <- values[idx]
        if (model == "technical_adjusted") {
          y <- residualize(y, covs); s <- residualize(s, covs); v <- residualize(v, covs)
        }
        joint <- z(y) + z(s)
        cor_rows[[kc]] <- data.frame(
          cohort = "Wu2021", patient = p, module = m, version = method, model = model,
          n_cells = length(idx), rho_with_joint = safe_cor(joint, v)
        ); kc <- kc + 1L
      }
    }
  }
}
write_tsv(bind_rows(cor_rows), "framework_module_correlations_wu.tsv")

# Technical-confounding audit for every predefined scoring axis.
technical_rows <- list(); kt <- 1L
score_vars <- list(YAP = meta$YAP_score, Stemness = meta$Stemness_score, Joint = meta$YAP_Stem_Score)
tech_vars <- list(log_nCount = log1p(meta$nCount_RNA), nFeature = meta$nFeature_RNA, Cell_cycle = meta$Cell_cycle)
for (scope in c("pooled", sort(unique(meta$patient)))) {
  idx <- if (scope == "pooled") seq_len(nrow(meta)) else which(meta$patient == scope)
  for (score_name in names(score_vars)) {
    for (covariate in names(tech_vars)) {
      technical_rows[[kt]] <- data.frame(
        cohort = "Wu2021", scope = ifelse(scope == "pooled", "pooled", "patient"),
        patient = ifelse(scope == "pooled", "POOLED", scope), score = score_name,
        covariate = covariate, n_cells = length(idx),
        spearman_rho = safe_cor(score_vars[[score_name]][idx], tech_vars[[covariate]][idx])
      ); kt <- kt + 1L
    }
  }
}
write_tsv(bind_rows(technical_rows), "score_technical_association_wu.tsv")

cat("Wu framework effect calculations completed\n")
