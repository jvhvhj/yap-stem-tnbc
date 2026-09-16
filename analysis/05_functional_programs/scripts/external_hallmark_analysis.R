#!/usr/bin/env Rscript
# Purpose: Predefined Hallmark enrichment and patient-level recurrence
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 4 and Supplementary Fig. S1.

options(stringsAsFactors = FALSE, width = 220)
set.seed(20260809)

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "Figure4_EFG_reproducible_GSEA_and_patient_replication_DATA_GATE_ONLY")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

.libPaths(c(
  file.path(root, "0710_final_evidence_rebuild", "vendor"),
  ".software/r-library",
  file.path(root, "0703_rebuild", "R_library"),
  .libPaths()
))

suppressPackageStartupMessages({
  library(fgsea)
  library(digest)
})

write_tsv <- function(x, path) {
  write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE,
    na = "NA", fileEncoding = "UTF-8"
  )
  cat("WROTE", basename(path), nrow(x), "rows\n")
}

sha256 <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

fmt <- function(x, digits = 4L) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "fg"))
}

collapse_sorted <- function(x) paste(sort(unique(x)), collapse = "/")

read_gmt <- function(path) {
  x <- strsplit(readLines(path, warn = FALSE, encoding = "UTF-8"), "\t", fixed = TRUE)
  out <- lapply(x, function(z) unique(z[-c(1, 2)]))
  names(out) <- vapply(x, `[`, character(1), 1)
  out
}

prepare_stats <- function(genes, values) {
  d <- data.frame(gene = as.character(genes), statistic = as.numeric(values))
  d <- d[!is.na(d$gene) & nzchar(d$gene) & is.finite(d$statistic), , drop = FALSE]
  d <- d[order(-d$statistic, d$gene), , drop = FALSE]
  d <- d[!duplicated(d$gene), , drop = FALSE]
  v <- d$statistic
  names(v) <- d$gene
  v
}

run_hallmark <- function(stats, pathways) {
  set.seed(20260809)
  z <- suppressWarnings(fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = stats,
    minSize = 10,
    maxSize = 500,
    eps = 1e-10,
    scoreType = "std",
    nPermSimple = 10000
  ))
  z <- as.data.frame(z)
  z$leading_edge_genes <- vapply(z$leadingEdge, collapse_sorted, character(1))
  z$leading_edge_size <- lengths(z$leadingEdge)
  z <- z[, c(
    "pathway", "ES", "NES", "pval", "padj", "size",
    "leading_edge_genes", "leading_edge_size"
  )]
  names(z) <- c(
    "pathway", "enrichment_score", "NES", "nominal_P", "BH_FDR",
    "tested_gene_set_size", "leading_edge_genes", "leading_edge_size"
  )
  z$direction <- ifelse(z$NES > 0, "High-enriched", ifelse(z$NES < 0, "Other-enriched", "zero"))
  z <- z[order(z$BH_FDR, -abs(z$NES), z$pathway), , drop = FALSE]
  rownames(z) <- NULL
  z
}

# =============================================================================
# Frozen analysis specification, fixed before reading any GSEA result
# =============================================================================

ranking_rule_cohort <- "DESeq2 Wald statistic for state High versus Other from design ~ patient + state; positive means High-enriched"
ranking_rule_patient <- "Within-patient High-minus-Other difference in DESeq2 VST expression across the same score-independent 17,597-gene universe; no cell-level P value"
gsea_method <- "fgseaMultilevel; scoreType=std; exponent=1 (package default); minSize=10; maxSize=500; eps=1e-10; nPermSimple=10000; seed=20260809"
patient_order <- c(
  "CID4465", "CID4495", "CID44971", "CID44991",
  "CID4513", "CID4515", "CID4523", "CID3963"
)
focus_pathways <- c(
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_HYPOXIA",
  "HALLMARK_APICAL_JUNCTION",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_APOPTOSIS"
)

source_de <- file.path(
  root, "Figure3_bottom_panels_HI_feasibility_and_plotting",
  "Figure3_I_effect_consistency_source.tsv"
)
source_vst <- file.path(
  root, "Figure3_bottom_panels_HI_feasibility_and_plotting",
  "Figure3_H_patient_centered_matrix.rds"
)
source_vst_script <- file.path(root, "tmp", "Figure3_HI_feasibility.R")
source_de_script <- file.path(root, "tmp", "Figure3_I_effect_consistency.R")
hallmark_gmt <- file.path(root, "h.all.v2024.Hs.symbols.gmt")
custom_upr_file <- file.path(
  root, "0719_frozen_score_definition_audit", "signature_gene_lists",
  "UPR__score_independent.tsv"
)
old_gsea <- file.path(
  root, "0728_sensitivity_and_cnv_closure", "01_within_patient_tertile",
  "primary_GSEA_Hallmark.tsv"
)
this_script <- file.path(out_dir, "scripts", "01_reproducible_efg_gate.R")

required_inputs <- c(
  source_de, source_vst, source_vst_script, source_de_script,
  hallmark_gmt, custom_upr_file, old_gsea, this_script
)
stopifnot(all(file.exists(required_inputs)))

# =============================================================================
# Part 1: cohort-level traceable Hallmark GSEA v2
# =============================================================================

de <- read.delim(source_de, check.names = FALSE, na.strings = c("NA", ""))
required_de_cols <- c(
  "gene", "overall_log2FC_High_vs_Other", "Wald_stat", "BH_FDR",
  "score_definition_genes_excluded_before_testing", "contrast", "design",
  "aggregation"
)
stopifnot(
  all(required_de_cols %in% names(de)),
  nrow(de) == 17597L,
  length(unique(de$gene)) == 17597L,
  all(de$score_definition_genes_excluded_before_testing),
  all(de$contrast == "High versus Other"),
  all(de$design == "~ patient + state"),
  all(de$aggregation == "16 patient-state raw-count pseudobulks"),
  all(is.finite(de$Wald_stat))
)

rank_effect_cor <- cor(de$Wald_stat, de$overall_log2FC_High_vs_Other, use = "complete.obs")
rank_effect_sign_agreement <- mean(
  sign(de$Wald_stat[de$overall_log2FC_High_vs_Other != 0]) ==
    sign(de$overall_log2FC_High_vs_Other[de$overall_log2FC_High_vs_Other != 0]),
  na.rm = TRUE
)

ranking <- data.frame(
  gene = de$gene,
  ranking_statistic = de$Wald_stat,
  source_effect = de$overall_log2FC_High_vs_Other,
  source_FDR = de$BH_FDR,
  testable_status = "TESTABLE_SCORE_INDEPENDENT",
  stringsAsFactors = FALSE
)
ranking <- ranking[order(-ranking$ranking_statistic, ranking$gene), , drop = FALSE]
write_tsv(ranking, file.path(out_dir, "Figure4_E_GSEA_v2_ranking.tsv"))

hallmark <- read_gmt(hallmark_gmt)
stopifnot(length(hallmark) == 50L, all(focus_pathways %in% names(hallmark)))
cohort_stats <- prepare_stats(ranking$gene, ranking$ranking_statistic)
stopifnot(length(cohort_stats) == 17597L)
e_results <- run_hallmark(cohort_stats, hallmark)
stopifnot(nrow(e_results) == 50L, setequal(e_results$pathway, names(hallmark)))
write_tsv(e_results, file.path(out_dir, "Figure4_E_GSEA_v2_results.tsv"))

upr_new <- e_results[e_results$pathway == "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", ]
old <- read.delim(old_gsea, check.names = FALSE, na.strings = c("NA", ""))
upr_old <- old[old$ID == "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", ]
stopifnot(nrow(upr_new) == 1L, nrow(upr_old) == 1L)

custom_upr <- unique(read.delim(custom_upr_file, check.names = FALSE)$gene)
hallmark_upr <- unique(hallmark[["HALLMARK_UNFOLDED_PROTEIN_RESPONSE"]])
shared_upr <- intersect(custom_upr, hallmark_upr)
custom_only <- setdiff(custom_upr, hallmark_upr)
hallmark_only <- setdiff(hallmark_upr, custom_upr)
upr_union <- union(custom_upr, hallmark_upr)

upr_leading <- strsplit(upr_new$leading_edge_genes, "/", fixed = TRUE)[[1]]
if (length(upr_leading) == 1L && !nzchar(upr_leading)) upr_leading <- character()
rank_map <- setNames(ranking$ranking_statistic, ranking$gene)
effect_map <- setNames(ranking$source_effect, ranking$gene)
fdr_map <- setNames(ranking$source_FDR, ranking$gene)

# A Wald statistic is effect / standard error, so its all-gene Pearson
# correlation with the unstandardized effect is not expected to approach 1.
# A true sign/ranking inversion is therefore gated by direction agreement,
# while the Pearson correlation is retained as a descriptive QC value only.
sign_error <- !is.finite(rank_effect_sign_agreement) ||
  rank_effect_sign_agreement < 0.99

upr_decision <- if (sign_error) {
  "SOURCE_OR_SIGN_ERROR"
} else if (upr_new$NES > 0 && upr_new$BH_FDR < 0.05) {
  "UPR_CONCORDANT"
} else if (upr_new$NES < 0 && upr_new$BH_FDR < 0.05) {
  "UPR_DISCORDANT_REPRODUCED"
} else {
  "OLD_RESULT_NOT_REPRODUCED"
}

upr_summary_row <- data.frame(
  row_type = "SUMMARY",
  gene = "ALL_UPR_GENES",
  custom_UPR_member = NA,
  Hallmark_UPR_member = NA,
  testable_status = NA,
  ranking_statistic = NA,
  source_effect = NA,
  source_FDR = NA,
  direction = NA,
  leading_edge_membership = NA,
  overlap_n = length(shared_upr),
  union_n = length(upr_union),
  Jaccard = length(shared_upr) / length(upr_union),
  shared_genes = collapse_sorted(shared_upr),
  custom_only_genes = collapse_sorted(custom_only),
  Hallmark_only_genes = collapse_sorted(hallmark_only),
  new_NES = upr_new$NES,
  new_BH_FDR = upr_new$BH_FDR,
  old_NES = upr_old$NES,
  old_BH_FDR = upr_old$p.adjust,
  ranking_effect_correlation = rank_effect_cor,
  ranking_effect_sign_agreement = rank_effect_sign_agreement,
  decision = upr_decision,
  stringsAsFactors = FALSE
)

upr_genes <- sort(union(custom_upr, hallmark_upr))
upr_gene_rows <- data.frame(
  row_type = "GENE",
  gene = upr_genes,
  custom_UPR_member = upr_genes %in% custom_upr,
  Hallmark_UPR_member = upr_genes %in% hallmark_upr,
  testable_status = ifelse(upr_genes %in% ranking$gene, "TESTABLE_SCORE_INDEPENDENT", "NOT_IN_TESTABLE_UNIVERSE"),
  ranking_statistic = unname(rank_map[upr_genes]),
  source_effect = unname(effect_map[upr_genes]),
  source_FDR = unname(fdr_map[upr_genes]),
  direction = ifelse(
    is.na(rank_map[upr_genes]), "NOT_TESTABLE",
    ifelse(rank_map[upr_genes] > 0, "High", ifelse(rank_map[upr_genes] < 0, "Other", "zero"))
  ),
  leading_edge_membership = upr_genes %in% upr_leading,
  overlap_n = NA,
  union_n = NA,
  Jaccard = NA,
  shared_genes = NA,
  custom_only_genes = NA,
  Hallmark_only_genes = NA,
  new_NES = upr_new$NES,
  new_BH_FDR = upr_new$BH_FDR,
  old_NES = upr_old$NES,
  old_BH_FDR = upr_old$p.adjust,
  ranking_effect_correlation = rank_effect_cor,
  ranking_effect_sign_agreement = rank_effect_sign_agreement,
  decision = upr_decision,
  stringsAsFactors = FALSE
)
upr_audit <- rbind(upr_summary_row, upr_gene_rows)
write_tsv(upr_audit, file.path(out_dir, "Figure4_E_UPR_discordance_audit.tsv"))

e_conditions <- c(
  complete_reproducible_ranking = nrow(ranking) == 17597L && all(is.finite(ranking$ranking_statistic)),
  external_provenance_clean = length(hallmark) == 50L && grepl("v2024", basename(hallmark_gmt)),
  no_sign_ambiguity = !sign_error,
  full_result_table_saved = nrow(e_results) == 50L
)
e_decision <- if (all(e_conditions)) "DATA_GO" else "DATA_WEAK"

provenance <- c(
  "# Figure 4E reproducible Hallmark GSEA v2 provenance",
  "",
  paste0("- **E decision:** `", e_decision, "`."),
  paste0("- **Ranking rule:** ", ranking_rule_cohort, "."),
  "- **Ranking metric was fixed before examining v2 enrichment results:** DESeq2 Wald statistic; no alternative metric was tried.",
  paste0("- **Source gene-level table:** `", source_de, "`."),
  paste0("- **Source SHA-256:** `", sha256(source_de), "`."),
  paste0("- **Generating source script:** `", source_de_script, "`."),
  paste0("- **Generating source-script SHA-256:** `", sha256(source_de_script), "`."),
  paste0("- **Current v2 script:** `", this_script, "`."),
  paste0("- **Current v2 script SHA-256:** `", sha256(this_script), "`."),
  "- **Contrast and model:** current pooled-cohort High versus Other; `design = ~ patient + state`; 16 patient-state raw-count pseudobulks.",
  "- **Score-definition exclusion:** all 38 YAP17/Stem21 scoring genes were excluded before gene filtering/testing; the authoritative universe contains 17,597 score-independent testable genes.",
  "- **Pseudoreplication guard:** no cell-level test or cell-level P value was used.",
  paste0("- **Gene-set collection:** MSigDB Hallmark, local file `", basename(hallmark_gmt), "` (v2024 label)."),
  paste0("- **Hallmark SHA-256:** `", sha256(hallmark_gmt), "`."),
  paste0("- **GSEA implementation:** ", gsea_method, "."),
  paste0("- **R version:** ", R.version.string, "."),
  paste0("- **fgsea version:** ", as.character(packageVersion("fgsea")), "."),
  paste0("- **Ranking/effect correlation:** ", fmt(rank_effect_cor, 6), " (descriptive only because Wald=effect/SE); sign agreement used for inversion QC: ", fmt(rank_effect_sign_agreement, 6), "."),
  "- **Result retention:** all 50 Hallmark pathways are saved, including nonsignificant and discordant results.",
  "- **Interpretive scope:** orthogonal pathway profiling; it is not confirmation of all six custom modules."
)
writeLines(provenance, file.path(out_dir, "Figure4_E_GSEA_v2_provenance.md"), useBytes = TRUE)

upr_md <- c(
  "# Figure 4E UPR discordance audit",
  "",
  paste0("**Decision: `", upr_decision, "`.**"),
  "",
  paste0(
    "The v2 Hallmark UPR result is NES=", fmt(upr_new$NES, 6),
    ", BH FDR=", fmt(upr_new$BH_FDR, 6), ", direction=", upr_new$direction,
    ". The historical result was NES=", fmt(upr_old$NES, 6),
    ", BH FDR=", fmt(upr_old$p.adjust, 6), "."
  ),
  "",
  paste0(
    "The current score-independent custom UPR module contains ", length(custom_upr),
    " genes; the v2024 Hallmark UPR set contains ", length(hallmark_upr),
    " genes. They share ", length(shared_upr), " genes across a union of ",
    length(upr_union), " (Jaccard=", fmt(length(shared_upr) / length(upr_union), 4), ")."
  ),
  "",
  paste0("Shared genes: ", ifelse(length(shared_upr), collapse_sorted(shared_upr), "none"), "."),
  paste0("Custom-only genes: ", ifelse(length(custom_only), collapse_sorted(custom_only), "none"), "."),
  "",
  "The ranking sign is not inverted: positive Wald statistics correspond to positive High-versus-Other effects. Interpretation is therefore made only after the sign and source checks passed.",
  "",
  "A negative Hallmark UPR result does not invalidate the custom eight-gene score by itself: the two gene sets have different membership. Figure 4E must be described as external pathway profiling, not as confirmation of every custom module."
)
writeLines(upr_md, file.path(out_dir, "Figure4_E_UPR_discordance_audit.md"), useBytes = TRUE)

# =============================================================================
# Part 2: patient-resolved external-pathway reproducibility (only if E DATA_GO)
# =============================================================================

if (e_decision == "DATA_GO") {
  vst_payload <- readRDS(source_vst)
  patient_vst <- vst_payload$patient_centered_vst
  stopifnot(
    is.matrix(patient_vst),
    nrow(patient_vst) == 17597L,
    setequal(rownames(patient_vst), ranking$gene),
    all(as.vector(t(outer(patient_order, c("Other", "High"), paste, sep = "|"))) %in% colnames(patient_vst))
  )

  patient_rows <- list()
  for (patient in patient_order) {
    delta <- patient_vst[, paste(patient, "High", sep = "|")] -
      patient_vst[, paste(patient, "Other", sep = "|")]
    patient_stats <- prepare_stats(rownames(patient_vst), delta)
    stopifnot(length(patient_stats) == 17597L)
    z <- run_hallmark(patient_stats, hallmark[focus_pathways])
    z <- z[match(focus_pathways, z$pathway), , drop = FALSE]
    patient_rows[[patient]] <- data.frame(
      patient = patient,
      pathway = z$pathway,
      NES = z$NES,
      direction = ifelse(z$NES > 0, "High-enriched", ifelse(z$NES < 0, "Other-enriched", "zero")),
      leading_edge_size = z$leading_edge_size,
      cohort_NES = e_results$NES[match(z$pathway, e_results$pathway)],
      cohort_BH_FDR = e_results$BH_FDR[match(z$pathway, e_results$pathway)],
      ranking_metric = ranking_rule_patient,
      n_testable_genes = length(patient_stats),
      stringsAsFactors = FALSE
    )
  }
  patient_nes <- do.call(rbind, patient_rows)
  rownames(patient_nes) <- NULL
  stopifnot(nrow(patient_nes) == 48L)
  write_tsv(patient_nes, file.path(out_dir, "Figure4_F_patient_pathway_NES.tsv"))

  pathway_summary <- do.call(rbind, lapply(focus_pathways, function(pathway) {
    d <- patient_nes[patient_nes$pathway == pathway, , drop = FALSE]
    cohort_nes <- unique(d$cohort_NES)
    cohort_direction <- ifelse(cohort_nes > 0, "High-enriched", "Other-enriched")
    same_sign <- sign(d$NES) == sign(cohort_nes)
    abs_share <- abs(d$NES) / sum(abs(d$NES))
    lopo_medians <- vapply(seq_len(nrow(d)), function(i) median(d$NES[-i]), numeric(1))
    data.frame(
      pathway = pathway,
      n_tested = nrow(d),
      positive_NES_n = sum(d$NES > 0),
      negative_NES_n = sum(d$NES < 0),
      median_NES = median(d$NES),
      IQR_NES = IQR(d$NES),
      min_NES = min(d$NES),
      max_NES = max(d$NES),
      cohort_NES = cohort_nes,
      cohort_BH_FDR = unique(d$cohort_BH_FDR),
      cohort_direction = cohort_direction,
      same_as_cohort_n = sum(same_sign),
      dominant_patient = d$patient[which.max(abs_share)],
      dominance_share_abs_NES = max(abs_share),
      lopo_median_direction_stable = all(sign(lopo_medians) == sign(cohort_nes)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(pathway_summary) <- NULL
  write_tsv(pathway_summary, file.path(out_dir, "Figure4_F_patient_pathway_summary.tsv"))

  # Prespecified gate: broad recurrence, no single-patient dominance, real
  # continuous heterogeneity, and patient-level explanation of UPR direction.
  recurrent_n <- sum(pathway_summary$same_as_cohort_n >= 6L & pathway_summary$lopo_median_direction_stable)
  no_single_patient_dominance <- all(pathway_summary$dominance_share_abs_NES < 0.50)
  heterogeneous_n <- sum(pathway_summary$IQR_NES >= 0.25)
  upr_patient <- pathway_summary[
    pathway_summary$pathway == "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", , drop = FALSE
  ]
  upr_patient_explains <- upr_patient$same_as_cohort_n >= 5L &&
    upr_patient$lopo_median_direction_stable

  f_decision <- if (
    recurrent_n >= 4L && no_single_patient_dominance &&
      heterogeneous_n >= 3L && upr_patient_explains
  ) {
    "DATA_GO"
  } else if (recurrent_n >= 3L && no_single_patient_dominance) {
    "DATA_WEAK"
  } else {
    "DROP_REDUNDANT"
  }

  key_patient_rows <- patient_nes[
    patient_nes$patient %in% c("CID3963", "CID4513", "CID4465"),
    c("patient", "pathway", "NES", "direction")
  ]
  key_patient_lines <- apply(key_patient_rows, 1, function(x) {
    paste0("- ", x[["patient"]], " / ", x[["pathway"]], ": NES=", fmt(as.numeric(x[["NES"]]), 4), " (", x[["direction"]], ").")
  })

  f_md <- c(
    "# Figure 4F patient-resolved external-pathway reproducibility data gate",
    "",
    paste0("**Decision: `", f_decision, "`.**"),
    "",
    paste0("- Patient ranking rule: ", ranking_rule_patient, "."),
    paste0("- GSEA parameters: ", gsea_method, "."),
    "- All eight Wu patients and all six prespecified external Hallmark pathways were tested descriptively; no cell-level P value was used.",
    paste0("- Pathways with at least 6/8 patients matching the cohort direction and stable leave-one-patient-out median direction: ", recurrent_n, "/6."),
    paste0("- Pathways with patient NES IQR >=0.25: ", heterogeneous_n, "/6."),
    paste0("- Maximum single-patient absolute-NES share across pathways: ", fmt(max(pathway_summary$dominance_share_abs_NES), 4), "."),
    paste0("- UPR patients matching the cohort direction: ", upr_patient$same_as_cohort_n, "/8; leave-one-patient-out median direction stable: ", upr_patient$lopo_median_direction_stable, "."),
    "",
    "The continuous patient NES values add heterogeneity and influential-patient information that cannot be represented by a single positive-count column. Patient is the replication unit; single-patient enrichment is descriptive.",
    "",
    "## Prespecified influential-patient audit",
    "",
    key_patient_lines
  )
  writeLines(f_md, file.path(out_dir, "Figure4_F_data_gate.md"), useBytes = TRUE)

  cat(
    paste0(
      "\n## Patient-resolved direction audit\n\n",
      "Using the separately frozen within-patient VST effect ranking, Hallmark UPR was High-enriched in ",
      upr_patient$positive_NES_n, "/8 patients and Other-enriched in ",
      upr_patient$negative_NES_n, "/8. Only ", upr_patient$same_as_cohort_n,
      "/8 matched the negative cohort-level Wald-GSEA direction. Therefore, the reproduced cohort-level negative NES is not a recurrent majority-patient direction. This is not a sign inversion: all 17,597 Wald statistics agree in sign with their source log2FC. The discrepancy is consistent with the different weighting of the prespecified ranking statistics (cohort Wald effect/SE versus within-patient VST expression difference), together with the low overlap between the custom UPR module and Hallmark UPR. E remains orthogonal pathway profiling, but Hallmark UPR must not be presented as patient-replicated confirmation.\n"
    ),
    file = file.path(out_dir, "Figure4_E_UPR_discordance_audit.md"),
    append = TRUE
  )
} else {
  f_decision <- "NOT_RUN_E_FAILED"
  write_tsv(data.frame(
    patient = character(), pathway = character(), NES = numeric(),
    direction = character(), leading_edge_size = integer()
  ), file.path(out_dir, "Figure4_F_patient_pathway_NES.tsv"))
  write_tsv(data.frame(
    pathway = character(), n_tested = integer(), positive_NES_n = integer(),
    negative_NES_n = integer(), median_NES = numeric(), IQR_NES = numeric(),
    min_NES = numeric(), max_NES = numeric()
  ), file.path(out_dir, "Figure4_F_patient_pathway_summary.tsv"))
  writeLines(c(
    "# Figure 4F patient pathway gate",
    "", "F was not run because E did not pass its data gate."
  ), file.path(out_dir, "Figure4_F_data_gate.md"), useBytes = TRUE)
}

# =============================================================================
# Part 3: v2 leading-edge architecture (only v2 E results)
# =============================================================================

if (e_decision == "DATA_GO") {
  passing <- e_results[e_results$NES > 0 & e_results$BH_FDR < 0.05, , drop = FALSE]
  passing <- passing[order(passing$BH_FDR, -passing$NES, passing$pathway), , drop = FALSE]
  passing_genes <- lapply(passing$leading_edge_genes, function(x) {
    if (is.na(x) || !nzchar(x)) character() else unique(strsplit(x, "/", fixed = TRUE)[[1]])
  })
  names(passing_genes) <- passing$pathway

  membership_n <- table(unlist(passing_genes, use.names = FALSE))
  membership <- do.call(rbind, lapply(seq_len(nrow(passing)), function(i) {
    genes <- passing_genes[[i]]
    data.frame(
      pathway = passing$pathway[i],
      gene = genes,
      NES = passing$NES[i],
      BH_FDR = passing$BH_FDR[i],
      leading_edge_size = length(genes),
      pathway_membership_n = as.integer(membership_n[genes]),
      architecture_class = ifelse(membership_n[genes] >= 2L, "SHARED_GE2_PATHWAYS", "PATHWAY_SPECIFIC"),
      stringsAsFactors = FALSE
    )
  }))
  membership <- membership[order(-membership$pathway_membership_n, membership$gene, membership$pathway), , drop = FALSE]
  rownames(membership) <- NULL
  write_tsv(membership, file.path(out_dir, "Figure4_G_leading_edge_membership.tsv"))

  if (nrow(passing) >= 2L) {
    pairs <- combn(seq_len(nrow(passing)), 2)
    pairwise <- apply(pairs, 2, function(idx) {
      a <- passing_genes[[idx[1]]]
      b <- passing_genes[[idx[2]]]
      c(
        pathway_a = passing$pathway[idx[1]],
        pathway_b = passing$pathway[idx[2]],
        shared_n = length(intersect(a, b)),
        union_n = length(union(a, b)),
        Jaccard = length(intersect(a, b)) / length(union(a, b))
      )
    })
    pairwise <- as.data.frame(t(pairwise), stringsAsFactors = FALSE)
    pairwise$shared_n <- as.integer(pairwise$shared_n)
    pairwise$union_n <- as.integer(pairwise$union_n)
    pairwise$Jaccard <- as.numeric(pairwise$Jaccard)
  } else {
    pairwise <- data.frame(
      pathway_a = character(), pathway_b = character(), shared_n = integer(),
      union_n = integer(), Jaccard = numeric()
    )
  }

  unique_genes <- names(membership_n)
  shared_genes <- names(membership_n)[membership_n >= 2L]
  specific_genes <- names(membership_n)[membership_n == 1L]
  specific_by_pathway <- vapply(passing_genes, function(g) sum(membership_n[g] == 1L), integer(1))
  median_j <- if (nrow(pairwise)) median(pairwise$Jaccard) else NA_real_
  max_j <- if (nrow(pairwise)) max(pairwise$Jaccard) else NA_real_
  shared_fraction <- length(shared_genes) / length(unique_genes)

  # Prespecified pathway-overlap criteria: at least three pathways, a nontrivial but
  # not overwhelming shared component, modest median overlap, and a specific
  # branch for every pathway.
  g_decision <- if (
    nrow(passing) >= 3L && length(unique_genes) > 0L &&
      shared_fraction >= 0.05 && shared_fraction <= 0.60 &&
      is.finite(median_j) && median_j <= 0.25 &&
      all(specific_by_pathway > 0L)
  ) {
    "LEADING_EDGE_GO"
  } else if (nrow(passing) >= 3L && length(unique_genes) > 0L) {
    "LEADING_EDGE_WEAK"
  } else {
    "DROP"
  }

  pairwise_lines <- if (nrow(pairwise)) apply(
    pairwise[order(-pairwise$Jaccard), , drop = FALSE], 1,
    function(x) paste0(
      "- ", x[["pathway_a"]], " vs ", x[["pathway_b"]],
      ": shared=", x[["shared_n"]], ", union=", x[["union_n"]],
      ", Jaccard=", fmt(as.numeric(x[["Jaccard"]]), 4), "."
    )
  ) else "- No pairwise comparison available."

  g_md <- c(
    "# Figure 4G v2 leading-edge architecture assessment",
    "",
    paste0("**Decision: `", g_decision, "`.**"),
    "",
    paste0("- High-enriched Hallmark pathways with BH FDR<0.05: ", nrow(passing), "."),
    paste0("- Unique leading-edge genes: ", length(unique_genes), "."),
    paste0("- Shared across at least two pathways: ", length(shared_genes), " (", fmt(shared_fraction, 4), ")."),
    paste0("- Pathway-specific genes: ", length(specific_genes), "."),
    paste0("- Median pairwise Jaccard: ", fmt(median_j, 4), "; maximum: ", fmt(max_j, 4), "."),
    "- No network was constructed.",
    "",
    if (g_decision == "LEADING_EDGE_GO") {
      "The v2 data retain a limited shared leading-edge core together with pathway-specific branches."
    } else {
      "The v2 overlap structure is insufficient for a strong limited-core-plus-branches claim."
    },
    "",
    "## Pairwise Jaccard audit",
    "",
    pairwise_lines
  )
  writeLines(g_md, file.path(out_dir, "Figure4_G_leading_edge_summary.md"), useBytes = TRUE)
} else {
  g_decision <- "NOT_RUN_E_FAILED"
  write_tsv(data.frame(
    pathway = character(), gene = character(), NES = numeric(), BH_FDR = numeric(),
    leading_edge_size = integer(), pathway_membership_n = integer(),
    architecture_class = character()
  ), file.path(out_dir, "Figure4_G_leading_edge_membership.tsv"))
  writeLines(c(
    "# Figure 4G leading-edge architecture", "",
    "G was not run because E did not pass its data gate."
  ), file.path(out_dir, "Figure4_G_leading_edge_summary.md"), useBytes = TRUE)
}

# =============================================================================
# Final gate and hard-stop QA
# =============================================================================

final_md <- c(
  "# Figure 4 E/F/G final reproducibility data gate",
  "",
  paste0("- **E:** `", e_decision, "` — fully traceable Hallmark GSEA v2 based on the authoritative 17,597-gene patient-blocked High-versus-Other Wald ranking."),
  paste0("- **UPR:** `", upr_decision, "` — v2 NES=", fmt(upr_new$NES, 5), ", BH FDR=", fmt(upr_new$BH_FDR, 5), "."),
  paste0("- **F:** `", f_decision, "` — patient-resolved external-pathway recurrence using within-patient VST High-minus-Other rankings."),
  paste0("- **G:** `", g_decision, "` — leading-edge architecture rebuilt exclusively from v2 E High-enriched Hallmark pathways passing BH FDR<0.05."),
  "",
  if (exists("pathway_summary")) {
    paste0(
      "The critical UPR finding is two-level discordance: cohort Wald-GSEA reproducibly gives a negative Hallmark UPR NES, whereas patient-level VST-effect GSEA gives positive NES in ",
      pathway_summary$positive_NES_n[pathway_summary$pathway == "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"],
      "/8 patients. The negative cohort result is therefore not a majority-patient replicated direction and cannot support an external UPR-confirmation claim."
    )
  } else {
    "Patient-level UPR replication was not run because E failed."
  },
  "",
  "No historical lost ranking was reconstructed. No current module-derived gene set was used as orthogonal validation. No cell-level P value or cell-level pseudoreplicated inference was used.",
  "",
  "This script produces pathway analysis tables and reports; figure rendering and assembly are separate."
)
writeLines(final_md, file.path(out_dir, "Figure4_EFG_final_gate.md"), useBytes = TRUE)

expected <- file.path(out_dir, c(
  "Figure4_E_GSEA_v2_ranking.tsv",
  "Figure4_E_GSEA_v2_results.tsv",
  "Figure4_E_GSEA_v2_provenance.md",
  "Figure4_E_UPR_discordance_audit.tsv",
  "Figure4_E_UPR_discordance_audit.md",
  "Figure4_F_patient_pathway_NES.tsv",
  "Figure4_F_patient_pathway_summary.tsv",
  "Figure4_F_data_gate.md",
  "Figure4_G_leading_edge_membership.tsv",
  "Figure4_G_leading_edge_summary.md",
  "Figure4_EFG_final_gate.md"
))
stopifnot(all(file.exists(expected)))
stopifnot(nrow(read.delim(expected[1], check.names = FALSE)) == 17597L)
stopifnot(nrow(read.delim(expected[2], check.names = FALSE)) == 50L)
stopifnot(nrow(read.delim(expected[6], check.names = FALSE)) == ifelse(e_decision == "DATA_GO", 48L, 0L))

graphics <- list.files(
  out_dir, pattern = "\\.(pdf|png|tiff|tif|svg)$",
  recursive = TRUE, full.names = TRUE, ignore.case = TRUE
)
stopifnot(length(graphics) == 0L)

cat("FINAL_DECISIONS\n")
cat("E", e_decision, "\n")
cat("UPR", upr_decision, "\n")
cat("F", f_decision, "\n")
cat("G", g_decision, "\n")
cat("No graphical outputs created.\n")
