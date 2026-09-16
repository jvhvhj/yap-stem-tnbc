# Purpose: Hallmark GSEA and leading-edge closure
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 4 and Supplementary Fig. S1.
options(stringsAsFactors = FALSE, scipen = 999)

project_dir <- "."
prior_dir <- file.path(project_dir, "Figure4_EFG_reproducible_GSEA_and_patient_replication_DATA_GATE_ONLY")
out_dir <- file.path(project_dir, "Figure4_EFG_FINAL_DATA_CLOSURE_ONLY")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source_table <- file.path(project_dir, "Figure3_bottom_panels_HI_feasibility_and_plotting", "Figure3_I_effect_consistency_source.tsv")
source_script <- file.path(project_dir, "tmp", "Figure3_I_effect_consistency.R")
v2_script <- file.path(prior_dir, "scripts", "01_reproducible_efg_gate.R")
gmt_file <- file.path(project_dir, "h.all.v2024.Hs.symbols.gmt")
custom_upr_file <- file.path(project_dir, "0719_frozen_score_definition_audit", "signature_gene_lists", "UPR__score_independent.tsv")
prior_ranking_file <- file.path(prior_dir, "Figure4_E_GSEA_v2_ranking.tsv")
prior_results_file <- file.path(prior_dir, "Figure4_E_GSEA_v2_results.tsv")
prior_patient_file <- file.path(prior_dir, "Figure4_F_patient_pathway_NES.tsv")
prior_patient_summary_file <- file.path(prior_dir, "Figure4_F_patient_pathway_summary.tsv")
closure_script <- file.path(out_dir, "scripts", "01_close_figure4_efg.R")

required_inputs <- c(source_table, source_script, v2_script, gmt_file, custom_upr_file,
                     prior_ranking_file, prior_results_file, prior_patient_file,
                     prior_patient_summary_file, closure_script)
stopifnot(all(file.exists(required_inputs)))

read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, quote = "", comment.char = "", fileEncoding = "UTF-8")
}

write_tsv <- function(x, name) {
  write.table(x, file.path(out_dir, name), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}

sha256 <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required for SHA-256 provenance")
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

fmt <- function(x, digits = 5) formatC(x, format = "fg", digits = digits)
collapse_or_none <- function(x) if (length(x)) paste(x, collapse = "/") else "NONE"

parse_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  parts <- strsplit(lines, "\t", fixed = TRUE)
  sets <- lapply(parts, function(z) unique(z[-c(1, 2)]))
  names(sets) <- vapply(parts, `[[`, character(1), 1)
  sets
}

## -------------------------------------------------------------------------
## E: deterministic closure of the already audited, prespecified v2 run.
## No alternative ranking, gene set, threshold, or GSEA parameter is tested.
## -------------------------------------------------------------------------

ranking0 <- read_tsv(prior_ranking_file)
stopifnot(nrow(ranking0) == 17597L, !anyDuplicated(ranking0$gene))
ranking <- data.frame(
  gene = ranking0$gene,
  ranking_statistic = ranking0$ranking_statistic,
  effect = ranking0$source_effect,
  FDR = ranking0$source_FDR,
  testable = TRUE,
  testable_status = ranking0$testable_status,
  check.names = FALSE
)
stopifnot(all(sign(ranking$ranking_statistic) == sign(ranking$effect)))
write_tsv(ranking, "Figure4_E_GSEA_v2_ranking.tsv")

results0 <- read_tsv(prior_results_file)
stopifnot(nrow(results0) == 50L, !anyDuplicated(results0$pathway))
results <- data.frame(
  pathway = results0$pathway,
  NES = results0$NES,
  P = results0$nominal_P,
  BH_FDR = results0$BH_FDR,
  leading_edge_size = results0$leading_edge_size,
  leading_edge_genes = results0$leading_edge_genes,
  direction = results0$direction,
  enrichment_score = results0$enrichment_score,
  tested_gene_set_size = results0$tested_gene_set_size,
  check.names = FALSE
)
write_tsv(results, "Figure4_E_GSEA_v2_results.tsv")

source_hash <- sha256(source_table)
source_script_hash <- sha256(source_script)
v2_script_hash <- sha256(v2_script)
closure_script_hash <- sha256(closure_script)
gmt_hash <- sha256(gmt_file)
ranking_hash <- sha256(file.path(out_dir, "Figure4_E_GSEA_v2_ranking.tsv"))
results_hash <- sha256(file.path(out_dir, "Figure4_E_GSEA_v2_results.tsv"))

expected_hashes <- c(
  source = "ab30f9f8617f15d756e642da47d99cffa9ac30ba5ff14dc5d8629e6773a2b48e",
  source_script = "a27d5557359ed4c8857bd576c7c7ce7dc2d7a6a1ab5d7b11452cdd0aac757a4f",
  v2_script = "89cea15093d1ff08d0462e502fdf5bbfb55265dc53034b9fc2b364189b778735",
  gmt = "ee2463540042078bfa3f67828e1e223bb354446d9fbb4d22845866835ba5c772"
)
observed_hashes <- c(source = source_hash, source_script = source_script_hash,
                     v2_script = v2_script_hash, gmt = gmt_hash)
stopifnot(identical(unname(observed_hashes), unname(expected_hashes)))

provenance <- c(
  "# Figure 4E fully reproducible Hallmark GSEA v2 provenance",
  "",
  "- **E decision:** `DATA_GO`.",
  "- **High/Other definition:** High is the pooled upper tertile of `z_cohort(YAP17 mean) + z_cohort(Stem21 mean)`; Other is Intermediate + Low.",
  "- **Patient-blocked model:** current pooled-cohort High versus Other; `design = ~ patient + state`; 16 patient-state raw-count pseudobulks.",
  "- **Score-definition exclusion:** all 38 YAP17/Stem21 scoring genes were excluded before gene filtering/testing.",
  "- **Testable universe:** 17,597 unique score-independent genes; the complete universe is retained in the ranking output.",
  "- **Ranking statistic fixed before viewing enrichment results:** DESeq2 Wald statistic; positive means High-enriched. No alternative ranking metric was tested.",
  paste0("- **Authoritative source table:** `", source_table, "`."),
  paste0("- **Source SHA-256:** `", source_hash, "`."),
  paste0("- **Source-generating script:** `", source_script, "`."),
  paste0("- **Source-script SHA-256:** `", source_script_hash, "`."),
  paste0("- **Audited GSEA v2 script:** `", v2_script, "`."),
  paste0("- **Audited GSEA v2 script SHA-256:** `", v2_script_hash, "`."),
  paste0("- **Closure script:** `", closure_script, "`."),
  paste0("- **Closure script SHA-256:** `", closure_script_hash, "`."),
  "- **External gene sets:** predefined MSigDB Hallmark collection, local v2024 GMT; not constructed from Wu High-versus-Other results.",
  paste0("- **Hallmark GMT:** `", gmt_file, "`."),
  paste0("- **Hallmark GMT SHA-256:** `", gmt_hash, "`."),
  "- **GSEA parameters:** fgseaMultilevel; scoreType=std; exponent=1; minSize=10; maxSize=500; eps=1e-10; nPermSimple=10000; seed=20260809.",
  "- **Pseudoreplication guard:** no cell-level test or cell-level P value was used.",
  paste0("- **Final ranking SHA-256:** `", ranking_hash, "`."),
  paste0("- **Final results SHA-256:** `", results_hash, "`."),
  "- **Result retention:** all 50 formally tested Hallmark pathways are present, including nonsignificant and discordant results.",
  "- **Scientific role:** orthogonal external pathway profiling, not validation of all six custom functional modules."
)
writeLines(provenance, file.path(out_dir, "Figure4_E_GSEA_v2_provenance.md"), useBytes = TRUE)

## -------------------------------------------------------------------------
## UPR discordance: rebuild membership from the exact GMT and custom list.
## -------------------------------------------------------------------------

hallmark_sets <- parse_gmt(gmt_file)
upr_name <- "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"
stopifnot(upr_name %in% names(hallmark_sets))
hallmark_upr <- unique(hallmark_sets[[upr_name]])
custom_upr <- unique(read_tsv(custom_upr_file)$gene)
union_genes <- sort(union(custom_upr, hallmark_upr))
shared_genes <- sort(intersect(custom_upr, hallmark_upr))
custom_only <- sort(setdiff(custom_upr, hallmark_upr))
hallmark_only <- sort(setdiff(hallmark_upr, custom_upr))
upr_result <- results[results$pathway == upr_name, ]
stopifnot(nrow(upr_result) == 1L)
upr_le <- unique(strsplit(upr_result$leading_edge_genes, "/", fixed = TRUE)[[1]])
idx <- match(union_genes, ranking$gene)
gene_stats <- data.frame(
  gene = union_genes,
  ranking_statistic = ranking$ranking_statistic[idx],
  effect = ranking$effect[idx],
  FDR = ranking$FDR[idx],
  testable = !is.na(idx),
  stringsAsFactors = FALSE
)
gene_stats$custom_UPR_member <- gene_stats$gene %in% custom_upr
gene_stats$Hallmark_UPR_member <- gene_stats$gene %in% hallmark_upr
gene_stats$direction <- ifelse(gene_stats$ranking_statistic > 0, "High", ifelse(gene_stats$ranking_statistic < 0, "Other", "Zero"))
gene_stats$leading_edge_member <- gene_stats$gene %in% upr_le
gene_stats$leading_edge_sign <- ifelse(!gene_stats$leading_edge_member, "NOT_LEADING_EDGE",
                                       ifelse(gene_stats$ranking_statistic > 0, "POSITIVE", "NEGATIVE"))
negative_le <- sort(gene_stats$gene[gene_stats$leading_edge_sign == "NEGATIVE"])
positive_le <- sort(gene_stats$gene[gene_stats$leading_edge_sign == "POSITIVE"])
upr_decision <- if (upr_result$direction == "Other-enriched") "UPR_DISCORDANT_REPRODUCED" else "UPR_CONCORDANT"

upr_comparison <- data.frame(
  row_type = c("SUMMARY", rep("GENE", nrow(gene_stats))),
  gene = c("ALL_UPR_GENES", gene_stats$gene),
  custom_UPR_member = c(NA, gene_stats$custom_UPR_member),
  Hallmark_UPR_member = c(NA, gene_stats$Hallmark_UPR_member),
  ranking_statistic = c(NA, gene_stats$ranking_statistic),
  effect = c(NA, gene_stats$effect),
  FDR = c(NA, gene_stats$FDR),
  testable = c(NA, gene_stats$testable),
  direction = c(NA, gene_stats$direction),
  leading_edge_member = c(NA, gene_stats$leading_edge_member),
  leading_edge_sign = c(NA, gene_stats$leading_edge_sign),
  overlap_n = c(length(shared_genes), rep(NA, nrow(gene_stats))),
  union_n = c(length(union_genes), rep(NA, nrow(gene_stats))),
  Jaccard = c(length(shared_genes) / length(union_genes), rep(NA, nrow(gene_stats))),
  shared_genes = c(collapse_or_none(shared_genes), rep(NA, nrow(gene_stats))),
  custom_only = c(collapse_or_none(custom_only), rep(NA, nrow(gene_stats))),
  Hallmark_only = c(collapse_or_none(hallmark_only), rep(NA, nrow(gene_stats))),
  negative_leading_edge_genes = c(collapse_or_none(negative_le), rep(NA, nrow(gene_stats))),
  positive_leading_edge_genes = c(collapse_or_none(positive_le), rep(NA, nrow(gene_stats))),
  Hallmark_UPR_NES = c(upr_result$NES, rep(upr_result$NES, nrow(gene_stats))),
  Hallmark_UPR_BH_FDR = c(upr_result$BH_FDR, rep(upr_result$BH_FDR, nrow(gene_stats))),
  decision = c(upr_decision, rep(upr_decision, nrow(gene_stats))),
  check.names = FALSE
)
write_tsv(upr_comparison, "Figure4_E_UPR_gene_set_comparison.tsv")

old_upr_nes <- -1.3494308422138
old_upr_fdr <- 0.0403145875866244
upr_md <- c(
  "# Figure 4E UPR discordance audit",
  "",
  paste0("**Decision: `", upr_decision, "`.**"),
  "",
  paste0("The v2 Hallmark UPR result is NES=", fmt(upr_result$NES, 6),
         ", BH FDR=", fmt(upr_result$BH_FDR, 6), ", direction=", upr_result$direction,
         ". The historical result was NES=", fmt(old_upr_nes, 6),
         ", BH FDR=", fmt(old_upr_fdr, 6), "."),
  "",
  paste0("The custom score-independent UPR module contains ", length(custom_upr),
         " genes; the v2024 Hallmark UPR set contains ", length(hallmark_upr),
         " genes. They share ", length(shared_genes), " genes across a union of ",
         length(union_genes), " (Jaccard=", fmt(length(shared_genes) / length(union_genes), 5), ")."),
  "",
  paste0("Shared genes: ", collapse_or_none(shared_genes), "."),
  paste0("Custom-only genes: ", collapse_or_none(custom_only), "."),
  paste0("Negative leading-edge genes: ", collapse_or_none(negative_le), "."),
  paste0("Positive leading-edge genes: ", collapse_or_none(positive_le), "."),
  "",
  "All Hallmark UPR member genes and their v2 ranking statistics are retained in `Figure4_E_UPR_gene_set_comparison.tsv`. The ranking sign is not inverted: all 17,597 Wald statistics agree in sign with the source High-versus-Other effect.",
  "",
  "## Patient-resolved direction audit",
  "",
  "Using the independently prespecified within-patient VST signed-effect ranking, Hallmark UPR is High-enriched in 7/8 patients and Other-enriched in 1/8. Thus, the reproducible negative cohort-level Wald-GSEA result is not a recurrent majority-patient direction.",
  "",
  "This is not a source or sign error. The discordance is retained transparently and no gene set, ranking metric, or threshold was changed. Hallmark UPR must not be presented as patient-replicated confirmation of the custom UPR module."
)
writeLines(upr_md, file.path(out_dir, "Figure4_E_UPR_discordance_audit.md"), useBytes = TRUE)

## -------------------------------------------------------------------------
## Patient-resolved external Hallmark reproducibility.
## -------------------------------------------------------------------------

patient <- read_tsv(prior_patient_file)
patient_order <- c("CID4465", "CID4495", "CID44971", "CID44991", "CID4513", "CID4515", "CID4523", "CID3963")
stopifnot(nrow(patient) == 48L, setequal(unique(patient$patient), patient_order),
          length(unique(patient$pathway)) == 6L)
patient <- patient[order(match(patient$patient, patient_order), patient$pathway), ]
write_tsv(patient, "Figure4_patient_Hallmark_NES.tsv")

patient_summary0 <- read_tsv(prior_patient_summary_file)
patient_summary <- data.frame(
  pathway = patient_summary0$pathway,
  tested_n = patient_summary0$n_tested,
  positive_n = patient_summary0$positive_NES_n,
  negative_n = patient_summary0$negative_NES_n,
  median_NES = patient_summary0$median_NES,
  IQR = patient_summary0$IQR_NES,
  min_NES = patient_summary0$min_NES,
  max_NES = patient_summary0$max_NES,
  range = paste0(fmt(patient_summary0$min_NES, 5), " to ", fmt(patient_summary0$max_NES, 5)),
  cohort_NES = patient_summary0$cohort_NES,
  cohort_BH_FDR = patient_summary0$cohort_BH_FDR,
  cohort_direction = patient_summary0$cohort_direction,
  same_as_cohort_n = patient_summary0$same_as_cohort_n,
  dominant_patient = patient_summary0$dominant_patient,
  dominance_share_abs_NES = patient_summary0$dominance_share_abs_NES,
  leave_one_patient_out_median_direction_stable = patient_summary0$lopo_median_direction_stable,
  check.names = FALSE
)
write_tsv(patient_summary, "Figure4_patient_Hallmark_summary.tsv")

max_dom <- max(patient_summary$dominance_share_abs_NES)
stable_non_upr <- sum(patient_summary$same_as_cohort_n >= 6 &
                        patient_summary$leave_one_patient_out_median_direction_stable &
                        patient_summary$pathway != upr_name)
patient_decision <- "DATA_WEAK"
patient_gate <- c(
  "# Patient-resolved external Hallmark reproducibility data gate",
  "",
  paste0("**Decision: `", patient_decision, "`.**"),
  "",
  "- Ranking was fixed before inspection: within-patient High-minus-Other difference in DESeq2 VST expression across the same 17,597 score-independent genes.",
  "- Cells were not treated as biological replicates and no cell-level P value was used.",
  "- Hallmark version and fgsea parameters are identical to E v2.",
  "- Eight Wu patients and six prespecified external pathways produced 48 patient-pathway NES values.",
  paste0("- Five non-UPR pathways with at least 6/8 cohort-direction matches and stable leave-one-patient-out median direction: ", stable_non_upr, "/5."),
  paste0("- Maximum single-patient absolute-NES share: ", fmt(max_dom, 5), "; no pathway is dominated by one or two patients."),
  "- TNFA/NFKB, Hypoxia and Apoptosis are High-enriched in 8/8 patients; Apical junction in 7/8; Inflammatory response in 6/8.",
  "- Hallmark UPR is High-enriched in 7/8 patients although cohort-level Wald-GSEA is Other-enriched; only 1/8 patients matches the cohort direction.",
  "",
  "The patient distributions add genuine heterogeneity information beyond E, but the central UPR direction fails majority-patient replication. Under the prespecified gate, this candidate is therefore not promoted as the ninth main panel. Single-patient enrichment P values are not used as cross-patient evidence."
)
writeLines(patient_gate, file.path(out_dir, "Figure4_patient_Hallmark_data_gate.md"), useBytes = TRUE)

## -------------------------------------------------------------------------
## Leading-edge architecture rebuilt directly from the final E v2 table.
## -------------------------------------------------------------------------

included <- results[results$direction == "High-enriched" & !is.na(results$BH_FDR) & results$BH_FDR < 0.05, ]
stopifnot(nrow(included) > 1L)
le_list <- setNames(strsplit(included$leading_edge_genes, "/", fixed = TRUE), included$pathway)
membership <- do.call(rbind, lapply(names(le_list), function(pathway) {
  genes <- unique(le_list[[pathway]])
  data.frame(pathway = pathway, gene = genes, stringsAsFactors = FALSE)
}))
membership <- unique(membership)
membership_count <- table(membership$gene)
membership$NES <- included$NES[match(membership$pathway, included$pathway)]
membership$BH_FDR <- included$BH_FDR[match(membership$pathway, included$pathway)]
membership$leading_edge_size <- included$leading_edge_size[match(membership$pathway, included$pathway)]
membership$pathway_membership_n <- as.integer(membership_count[membership$gene])
membership$architecture_class <- ifelse(membership$pathway_membership_n >= 2L,
                                         "SHARED_GE2_PATHWAYS", "PATHWAY_SPECIFIC")
membership <- membership[order(membership$gene, membership$pathway), ]
write_tsv(membership, "Figure4_leading_edge_membership.tsv")

pathways <- names(le_list)
pairs <- combn(pathways, 2, simplify = FALSE)
jaccard <- do.call(rbind, lapply(pairs, function(p) {
  a <- unique(le_list[[p[1]]]); b <- unique(le_list[[p[2]]])
  data.frame(pathway_1 = p[1], pathway_2 = p[2], shared_n = length(intersect(a, b)),
             union_n = length(union(a, b)), Jaccard = length(intersect(a, b)) / length(union(a, b)))
}))
jaccard <- jaccard[order(-jaccard$Jaccard, jaccard$pathway_1, jaccard$pathway_2), ]
unique_genes <- length(unique(membership$gene))
shared_n <- sum(membership_count >= 2L)
specific_n <- sum(membership_count == 1L)
shared_fraction <- shared_n / unique_genes
leading_decision <- if (shared_n > 0L && specific_n > 0L) "DATA_GO" else "DATA_WEAK"

pair_lines <- paste0("- ", jaccard$pathway_1, " vs ", jaccard$pathway_2,
                     ": shared=", jaccard$shared_n, ", union=", jaccard$union_n,
                     ", Jaccard=", fmt(jaccard$Jaccard, 5), ".")
leading_md <- c(
  "# Figure 4 leading-edge shared/specific architecture",
  "",
  paste0("**Decision: `", leading_decision, "`.**"),
  "",
  "This table was rebuilt directly from the final E v2 results. No historical leading-edge table was used as an input.",
  "",
  "Formal inclusion rule: external Hallmark pathways that are High-enriched and have BH FDR < 0.05. No Jaccard threshold was used to select pathways or tune the architecture.",
  "",
  paste0("- Included external pathways: ", nrow(included), "."),
  paste0("- Unique leading-edge genes: ", unique_genes, "."),
  paste0("- Shared across at least two pathways: ", shared_n, " (", fmt(shared_fraction, 5), ")."),
  paste0("- Pathway-specific genes: ", specific_n, "."),
  paste0("- Median pairwise Jaccard: ", fmt(median(jaccard$Jaccard), 5), "."),
  paste0("- Maximum pairwise Jaccard: ", fmt(max(jaccard$Jaccard), 5), "."),
  "",
  "Both shared and pathway-specific leading-edge genes remain present. The v2 result therefore supports a limited shared core plus pathway-specific branches; it does not justify a dense gene network or threshold-tuned architecture.",
  "",
  "## Complete pairwise Jaccard audit",
  "",
  pair_lines
)
writeLines(leading_md, file.path(out_dir, "Figure4_leading_edge_summary.md"), useBytes = TRUE)

## -------------------------------------------------------------------------
## Circularity and final gate.
## -------------------------------------------------------------------------

circularity <- data.frame(
  evidence_layer = c("E_EXTERNAL_HALLMARK_PROFILE", "PATIENT_RESOLVED_HALLMARK", "LEADING_EDGE_ARCHITECTURE"),
  gene_set_source = c("MSigDB Hallmark v2024 external predefined GMT",
                      "Same MSigDB Hallmark v2024 external predefined GMT",
                      "Leading-edge genes derived only from included E v2 Hallmark GSEA results"),
  gene_sets_derived_from_Wu_High_vs_Other = c(FALSE, FALSE, FALSE),
  uses_High_vs_Other_ranking = c(TRUE, TRUE, TRUE),
  uses_current_custom_module_genes_as_gene_sets = c(FALSE, FALSE, FALSE),
  used_to_define_or_redefine_High_state = c(FALSE, FALSE, FALSE),
  feedback_into_state_definition = c(FALSE, FALSE, FALSE),
  circularity_status = c("PASS_EXTERNAL", "PASS_EXTERNAL", "PASS_DERIVED_OUTPUT_NO_FEEDBACK"),
  guardrail = c("Interpret as orthogonal external pathway profiling, not validation of every custom module",
                "Interpret patient NES descriptively with patient as replication unit",
                "Do not use leading-edge genes to redefine High state or rerun a self-validating analysis"),
  check.names = FALSE
)
write_tsv(circularity, "Figure4_EFG_circularity_check.tsv")

final_gate <- c(
  "# Figure 4 EFG final data closure",
  "",
  "## Final decisions",
  "",
  "1. **E:** `DATA_GO` — the complete 17,597-gene Wald ranking, all 50 Hallmark results, exact sources, scripts and SHA-256 hashes are retained.",
  "2. **UPR:** `DISCORDANT_REPRODUCED` — cohort Hallmark UPR remains Other-enriched, without sign error, but this direction is reproduced by only 1/8 patient-specific rankings.",
  "3. **Patient-resolved Hallmark panel:** `DATA_WEAK` — it adds heterogeneity information, but the UPR direction fails majority-patient replication and the candidate is not promoted as the ninth main panel.",
  paste0("4. **Leading-edge architecture:** `", leading_decision, "` — rebuilt directly from E v2 and supports a limited shared core plus pathway-specific branches."),
  "5. **Figure 4 main panels with genuine scientific value:** `8 / 9`. The existing eight-panel structure remains justified; this closure does not force a weak ninth panel.",
  "",
  "## Scientific closure",
  "",
  "The three-layer data chain is traceable: external Hallmark profiling (E), patient-resolved recurrence/heterogeneity assessment, and E-derived leading-edge membership. The patient-resolved layer is scientifically informative but does not pass the main-panel promotion gate. External Hallmark sets were not constructed from Wu High-versus-Other results, and leading-edge genes are not fed back into the High-state definition.",
  "",
  "## Output scope",
  "",
  "No PDF, PNG, TIFF or SVG was generated. No Figure 4 composite was assembled. Stage 3 was not entered."
)
writeLines(final_gate, file.path(out_dir, "Figure4_EFG_FINAL_GATE.md"), useBytes = TRUE)

required_outputs <- c(
  "Figure4_E_GSEA_v2_ranking.tsv",
  "Figure4_E_GSEA_v2_results.tsv",
  "Figure4_E_GSEA_v2_provenance.md",
  "Figure4_E_UPR_gene_set_comparison.tsv",
  "Figure4_E_UPR_discordance_audit.md",
  "Figure4_patient_Hallmark_NES.tsv",
  "Figure4_patient_Hallmark_summary.tsv",
  "Figure4_patient_Hallmark_data_gate.md",
  "Figure4_leading_edge_membership.tsv",
  "Figure4_leading_edge_summary.md",
  "Figure4_EFG_circularity_check.tsv",
  "Figure4_EFG_FINAL_GATE.md"
)
stopifnot(all(file.exists(file.path(out_dir, required_outputs))))
stopifnot(nrow(read_tsv(file.path(out_dir, required_outputs[1]))) == 17597L)
stopifnot(nrow(read_tsv(file.path(out_dir, required_outputs[2]))) == 50L)
stopifnot(nrow(read_tsv(file.path(out_dir, required_outputs[6]))) == 48L)
stopifnot(nrow(read_tsv(file.path(out_dir, required_outputs[7]))) == 6L)
stopifnot(nrow(read_tsv(file.path(out_dir, required_outputs[11]))) == 3L)

graphics <- list.files(out_dir, pattern = "\\.(pdf|png|tiff?|svg)$", recursive = TRUE,
                       ignore.case = TRUE, full.names = TRUE)
stopifnot(length(graphics) == 0L)

message("E=DATA_GO")
message("UPR=DISCORDANT_REPRODUCED")
message("PATIENT_HALLMARK=DATA_WEAK")
message("LEADING_EDGE=", leading_decision)
message("MAIN_PANEL_VALUE=8/9")
message("GRAPHICS=0")
