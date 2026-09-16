#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 200)
set.seed(20260808)

root <- "F:/桌面/AHIPPO-YAP"
out_dir <- file.path(root, "Figure3_bottom_panels_HI_feasibility_and_plotting")
.libPaths(c(
  file.path(root, "0710_final_evidence_rebuild/vendor"),
  "C:/Users/86156/R/bioc-3.23",
  file.path(root, "0703_rebuild/R_library"),
  .libPaths()
))
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(DESeq2)
})

patient_order <- c("CID4465", "CID4495", "CID44971", "CID44991",
                   "CID4513", "CID4515", "CID4523", "CID3963")
state_order <- c("Other", "High")

meta_path <- file.path(
  root, "0703_rebuild/continuous_state_main_figures/Wu2021_continuous_state_cell_metadata.csv"
)
object_path <- file.path(root, "output_step0/Wu2021_TNBC_malignant.rds")
yap_path <- file.path(
  root, "0719_frozen_score_definition_audit/signature_gene_lists/YAP_clean__clean_score.tsv"
)
stem_path <- file.path(
  root, "0719_frozen_score_definition_audit/signature_gene_lists/Stemness_clean__clean_score.tsv"
)
tier_path <- file.path(
  root, "0728_sensitivity_and_cnv_closure/02_direction_gate_sensitivity/strict_gene_tiers.tsv"
)

meta <- read.csv(meta_path, check.names = FALSE)
obj <- readRDS(object_path)
DefaultAssay(obj) <- "RNA"
stopifnot(nrow(meta) == 10836L, all(meta$cell_id %in% colnames(obj)))

score_genes <- unique(c(
  read.delim(yap_path, check.names = FALSE)$gene,
  read.delim(stem_path, check.names = FALSE)$gene
))
stopifnot(length(score_genes) == 38L)

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")[, meta$cell_id, drop = FALSE]
state <- factor(ifelse(meta$YS_tertile == "High", "High", "Other"), levels = state_order)
sample_id <- paste(meta$patient, state, sep = "|")
sample_levels <- as.vector(t(outer(patient_order, state_order, paste, sep = "|")))
mm <- sparse.model.matrix(~ 0 + factor(sample_id, levels = sample_levels))
colnames(mm) <- sample_levels
pb <- counts %*% mm
colnames(pb) <- sample_levels

split_id <- strsplit(sample_levels, "\\|")
sm <- data.frame(
  sample = sample_levels,
  patient = factor(vapply(split_id, `[`, character(1), 1), levels = patient_order),
  state = factor(vapply(split_id, `[`, character(1), 2), levels = state_order)
)
rownames(sm) <- sm$sample
design <- model.matrix(~ patient + state, data = sm)
pb0 <- pb[!rownames(pb) %in% score_genes, , drop = FALSE]
keep <- filterByExpr(DGEList(pb0), design = design, min.count = 5)
pb1 <- pb0[keep, , drop = FALSE]
stopifnot(nrow(pb1) == 17597L)

dense <- round(as.matrix(pb1))
storage.mode(dense) <- "integer"
dds <- DESeqDataSetFromMatrix(dense, sm, design = ~ patient + state)
dds <- DESeq(dds, quiet = TRUE)
rr <- as.data.frame(results(
  dds, contrast = c("state", "High", "Other"),
  alpha = 0.05, independentFiltering = FALSE
))

norm <- counts(dds, normalized = TRUE)
delta <- sapply(patient_order, function(p) {
  norm[, paste(p, "High", sep = "|")] - norm[, paste(p, "Other", sep = "|")]
})
colnames(delta) <- patient_order

tier <- read.delim(tier_path, check.names = FALSE)
program_genes <- tier$gene[as.logical(tier$in_primary146)]
stopifnot(length(program_genes) == 146L)
representative_genes <- c(
  "MCL1", "PPP1R15A", "ATF3", "BIRC3", "TNFAIP3", "TNF", "CXCL2",
  "ICAM1", "JUN", "FOS", "TACSTD2", "LAMB3", "LAMC2"
)

out <- data.frame(
  gene = rownames(rr),
  overall_log2FC_High_vs_Other = rr$log2FoldChange,
  lfcSE = rr$lfcSE,
  Wald_stat = rr$stat,
  p_value = rr$pvalue,
  BH_FDR = rr$padj,
  baseMean = rr$baseMean,
  positive_patient_count = rowSums(delta > 0),
  negative_patient_count = rowSums(delta < 0),
  zero_patient_count = rowSums(delta == 0),
  median_patient_normalized_count_delta = apply(delta, 1, median),
  patient_direction_pattern = apply(delta, 1, function(v) {
    paste(ifelse(v > 0, "+", ifelse(v < 0, "-", "0")), collapse = "")
  }),
  program_146 = rownames(rr) %in% program_genes,
  representative_gene = rownames(rr) %in% representative_genes,
  score_definition_genes_excluded_before_testing = TRUE,
  contrast = "High versus Other",
  design = "~ patient + state",
  aggregation = "16 patient-state raw-count pseudobulks",
  stringsAsFactors = FALSE
)

criteria <- !is.na(out$BH_FDR) & out$BH_FDR < 0.05 &
  out$overall_log2FC_High_vs_Other >= 0.5 & out$positive_patient_count >= 6
out$meets_full_numeric_146_gate <- criteria
if (!setequal(out$gene[criteria], program_genes)) {
  missing_from_recomputed <- setdiff(program_genes, out$gene[criteria])
  extra_in_recomputed <- setdiff(out$gene[criteria], program_genes)
  stop(sprintf(
    "Recomputed formal gate differs from authoritative 146: missing=%d extra=%d",
    length(missing_from_recomputed), length(extra_in_recomputed)
  ))
}
stopifnot(all(representative_genes %in% out$gene))

out <- out[order(out$BH_FDR, -out$overall_log2FC_High_vs_Other, na.last = TRUE), ]
write.table(
  out, file.path(out_dir, "Figure3_I_effect_consistency_source.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

cat("N_GENES", nrow(out), "\n")
cat("PROGRAM146", sum(out$program_146), "\n")
cat("REPRESENTATIVE", sum(out$representative_gene), "\n")
cat("POSITIVE_PATIENT_COUNT_DISTRIBUTION\n")
print(table(out$positive_patient_count))
