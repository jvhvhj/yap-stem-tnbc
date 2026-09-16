#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220)
set.seed(20260727)

root <- normalizePath("F:/桌面/AHIPPO-YAP", winslash = "/", mustWork = TRUE)
out <- file.path(root, "0727_analysis_closure")
dirs <- file.path(out, c(
  "00_manifest", "01_definition_audit", "02_pseudobulk_DE",
  "03_strict142_reconciliation", "04_pathway_refresh",
  "05_CNV_crosscaller", "06_Yan_locked_revalidation",
  "07_QC_report", "logs", "tmp"
))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
Sys.setenv(TMPDIR = file.path(out, "tmp"), TEMP = file.path(out, "tmp"), TMP = file.path(out, "tmp"))
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
  library(ggplot2)
  library(ggrepel)
})

log_path <- file.path(out, "logs", "01_wu_pseudobulk_and_reconciliation.log")
logcon <- file(log_path, "wt")
sink(logcon, type = "output", split = TRUE)
sink(logcon, type = "message", append = TRUE)
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(logcon)
}, add = TRUE)

cat("Started:", format(Sys.time()), "\n")
cat("R:", R.version.string, "\n")

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  cat("Wrote", path, nrow(x), "rows\n")
}

read_gene_file <- function(path) {
  x <- read.delim(path, check.names = FALSE)
  col <- intersect(c("gene", "Gene", "symbol", "SYMBOL"), names(x))[1]
  if (is.na(col)) stop("No gene column: ", path)
  unique(trimws(as.character(x[[col]])))
}

sha256_file <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) return(NA_character_)
  digest::digest(file = path, algo = "sha256")
}

safe_pdf <- function(path, width, height, code) {
  grDevices::cairo_pdf(path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
}

meta_path <- file.path(root, "0703_rebuild/continuous_state_main_figures/Wu2021_continuous_state_cell_metadata.csv")
obj_path <- file.path(root, "output_step0/Wu2021_TNBC_malignant.rds")
full_obj_path <- file.path(root, "Wu2021_seurat.rds")
audit_path <- file.path(root, "0719_frozen_score_definition_audit/exact_signature_provenance.tsv")
old_all_path <- file.path(root, "0703_rebuild/continuous_state_main_figures/Figure4_DEG_pathway/tables/High_vs_Low_DE_all.csv")
old_score_ind_path <- file.path(root, "0713_score_independent_rebuild/tables/Fig4_DEG_score_independent.tsv")
old_script_path <- file.path(root, "0703_rebuild/scripts/12_build_continuous_state_figures.R")

yap_path <- file.path(root, "0719_frozen_score_definition_audit/signature_gene_lists/YAP_clean__clean_score.tsv")
stem_path <- file.path(root, "0719_frozen_score_definition_audit/signature_gene_lists/Stemness_clean__clean_score.tsv")
module_names <- c("UPR", "TNFA_NFKB", "Hypoxia", "Adhesion_Remodeling", "Wound_Healing", "Survival_Stress")
module_paths <- setNames(
  file.path(root, "0719_frozen_score_definition_audit/signature_gene_lists",
            paste0(module_names, "__score_independent.tsv")),
  module_names
)
inputs <- c(meta_path, obj_path, full_obj_path, audit_path, old_all_path, old_score_ind_path,
            old_script_path, yap_path, stem_path, unname(module_paths))
input_manifest <- data.frame(
  role = c(
    "frozen_cell_metadata", "frozen_malignant_object", "full_Wu_atlas",
    "definition_provenance", "historical_High_vs_Low_all",
    "historical_score_independent", "historical_generation_script",
    "YAP17", "Stem21", paste0("module_", module_names)
  ),
  path = normalizePath(inputs, winslash = "/", mustWork = TRUE),
  bytes = file.info(inputs)$size,
  modified = format(file.info(inputs)$mtime, "%Y-%m-%d %H:%M:%S"),
  sha256 = vapply(inputs, sha256_file, character(1)),
  stringsAsFactors = FALSE
)
write_tsv(input_manifest, file.path(out, "00_manifest", "input_file_manifest.tsv"))

meta <- read.csv(meta_path, check.names = FALSE)
obj <- readRDS(obj_path)
DefaultAssay(obj) <- "RNA"
if (nrow(meta) != 10836L) stop("Expected exactly 10,836 frozen malignant epithelial cells")
if (!all(meta$cell_id %in% colnames(obj))) stop("Frozen metadata barcodes do not match malignant object")
if (length(unique(meta$patient)) != 8L) stop("Expected exactly 8 evaluable Wu patients")
if (!all(c("YS_tertile", "YAP_score", "Stemness_score", "YAP_Stem_Score") %in% names(meta))) {
  stop("Frozen score/state columns are missing")
}
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")[, meta$cell_id, drop = FALSE]
if (any(counts@x < 0) || any(abs(counts@x - round(counts@x)) > 1e-8)) {
  stop("RNA counts layer is not non-negative integer raw counts")
}

yap <- read_gene_file(yap_path)
stem <- read_gene_file(stem_path)
score_genes <- unique(c(yap, stem))
if (length(yap) != 17L || length(stem) != 21L || length(score_genes) != 38L) {
  stop("Frozen YAP17/Stem21/union38 invariant failed")
}
modules <- lapply(module_paths, read_gene_file)

gene_rows <- list(
  data.frame(set = "YAP clean", version = "clean_score", gene = yap),
  data.frame(set = "Stemness clean", version = "clean_score", gene = stem)
)
for (nm in names(modules)) {
  gene_rows[[length(gene_rows) + 1L]] <- data.frame(
    set = gsub("_", " ", nm), version = "score_independent", gene = modules[[nm]]
  )
}
frozen_gene_sets <- do.call(rbind, gene_rows)
frozen_gene_sets$in_RNA_counts <- frozen_gene_sets$gene %in% rownames(counts)
frozen_gene_sets$set_gene_count <- ave(frozen_gene_sets$gene, frozen_gene_sets$set, FUN = length)
write_tsv(frozen_gene_sets, file.path(out, "01_definition_audit", "frozen_gene_sets.tsv"))

barcodes <- data.frame(
  barcode = meta$cell_id,
  patient = meta$patient,
  YS_tertile = meta$YS_tertile,
  stringsAsFactors = FALSE
)
con <- gzfile(file.path(out, "01_definition_audit", "core_cell_barcodes.tsv.gz"), "wt")
write.table(barcodes, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)

patient_order <- unique(as.character(meta$patient))
conflicts <- data.frame(
  item = c(
    "YAP gene set", "Stemness gene set", "Joint-axis standardization",
    "continuous trend bins", "state contrast", "historical DEG contrast",
    "historical score-gene removal", "current score-gene removal"
  ),
  historical_variant = c(
    "4-gene variants exist in exploratory scripts", "5-gene variants exist in exploratory scripts",
    "benchmark sensitivity used within-patient z in selected subanalyses",
    "top-decile/top-vs-bottom variants exist", "High vs Low used by historical Figure 4",
    "High vs Low", "post-test removal in 0713 rebuild", "pre-test removal"
  ),
  frozen_current = c(
    "YAP clean 17 genes", "Stemness clean 21 genes",
    "cohort-pooled z(YAP mean)+z(Stemness mean)",
    "patient-resolved pooled-axis deciles; descriptive only",
    "High vs Other; High=upper tertile, Other=middle+lower tertiles",
    "High vs Other", "not accepted as circularity-safe", "all 38 genes excluded before filter/model"
  ),
  decision = c(
    rep("USE_FROZEN_CURRENT", 5),
    "HISTORICAL_ONLY_NOT_DIRECTLY_REPRODUCIBLE_UNDER_NEW_CONTRAST",
    "AUDIT_AS_HISTORICAL_LIMITATION", "USE_CURRENT"
  ),
  source = c(
    rep("0719_frozen_score_definition_audit", 5),
    "0703_rebuild/scripts/12_build_continuous_state_figures.R",
    "0713_score_independent_rebuild/scripts/00_score_independent_reanalysis.R",
    "0727 task prespecification"
  ),
  stringsAsFactors = FALSE
)
write_tsv(conflicts, file.path(out, "01_definition_audit", "definition_conflicts.tsv"))

yaml_lines <- c(
  "analysis_id: 0727_analysis_closure",
  "random_seed: 20260727",
  "cohort: Wu2021_TNBC",
  "patients:",
  paste0("  - ", patient_order),
  "n_patients: 8",
  "n_malignant_epithelial_cells: 10836",
  paste0("cell_barcode_sha256: ", sha256_file(file.path(out, "01_definition_audit", "core_cell_barcodes.tsv.gz"))),
  "scores:",
  "  expression_layer: RNA/data_log_normalized_for_score_assignment",
  "  YAP: arithmetic_mean_of_17_frozen_genes",
  "  Stemness: arithmetic_mean_of_21_frozen_genes",
  "  standardization: cohort_pooled_z_across_10836_cells",
  "  joint_axis: z(YAP_mean)+z(Stemness_mean)",
  "states:",
  "  Low: pooled_lower_tertile",
  "  Intermediate: pooled_middle_tertile",
  "  High: pooled_upper_tertile",
  "  DE_contrast: High_vs_Other",
  "  Other: Intermediate_plus_Low",
  "  deciles: pooled_joint_axis_deciles_used_only_for_descriptive_continuous_trends",
  "differential_expression:",
  "  input: raw_RNA_counts",
  "  aggregation: patient_x_state_sum",
  "  design: '~ patient + state'",
  "  score_genes_excluded_before_testing: 38",
  "Yan:",
  "  main_patient_threshold: '>=50 author-labelled cancer cells'",
  "  sensitivity_thresholds: [20, 100]",
  "representative_genes: [MCL1, PPP1R15A, ATF3]"
)
writeLines(yaml_lines, file.path(out, "01_definition_audit", "frozen_definitions.yaml"), useBytes = TRUE)

audit_report <- c(
  "# Frozen-definition audit",
  "",
  "**Status: PASS for the prespecified 0727 High-versus-Other analysis.**",
  "",
  paste0("- Wu analysis population: 10,836 malignant epithelial cells from ", length(patient_order), " evaluable patients."),
  paste0("- Frozen patient order: ", paste(patient_order, collapse = ", "), "."),
  "- Primary scores: arithmetic means of the frozen YAP17 and Stem21 genes on log-normalized RNA expression.",
  "- Joint axis: cohort-pooled z(YAP mean) + cohort-pooled z(Stemness mean).",
  "- Low/Intermediate/High: pooled tertiles. High-versus-Other uses High against Intermediate+Low.",
  "- Deciles are retained only for descriptive continuous trends.",
  "- Differential expression uses raw RNA counts, patient×state sums and design ~ patient + state.",
  "- All 38 score-defining genes are removed before expression filtering and model construction.",
  "",
  "Historical definitions are not silently substituted. The conflict table records older 4/5-gene, top-decile and High-versus-Low variants."
)
writeLines(audit_report, file.path(out, "01_definition_audit", "definition_audit_report.md"), useBytes = TRUE)

# ---------------------------------------------------------------------
# Phase 1: frozen High vs Other patient × state raw-count pseudobulk.
# ---------------------------------------------------------------------
meta$state_binary <- factor(ifelse(meta$YS_tertile == "High", "High", "Other"),
                            levels = c("Other", "High"))
sample_id <- paste(meta$patient, meta$state_binary, sep = "|")
sid_levels <- as.vector(t(outer(patient_order, c("Other", "High"), paste, sep = "|")))
mm <- sparse.model.matrix(~ 0 + factor(sample_id, levels = sid_levels))
pb <- counts %*% mm
if (ncol(pb) != 16L) stop("Expected 16 patient × state pseudobulks")
colnames(pb) <- sid_levels

split_id <- strsplit(colnames(pb), "\\|")
pb_meta <- data.frame(
  sample = colnames(pb),
  patient = vapply(split_id, `[`, character(1), 1),
  state = vapply(split_id, `[`, character(1), 2),
  stringsAsFactors = FALSE
)
pb_meta$patient <- factor(pb_meta$patient, levels = patient_order)
pb_meta$state <- factor(pb_meta$state, levels = c("Other", "High"))
rownames(pb_meta) <- pb_meta$sample
cell_tab <- table(factor(sample_id, levels = colnames(pb)))
pb_meta$cell_count <- as.integer(cell_tab[pb_meta$sample])
pb_meta$total_UMI <- as.numeric(Matrix::colSums(pb))
pb_meta$library_size <- pb_meta$total_UMI
pb_meta$detected_genes <- as.numeric(Matrix::colSums(pb > 0))
if (any(pb_meta$cell_count <= 0) || any(table(pb_meta$patient) != 2L)) {
  stop("Every patient must contribute non-empty High and Other pseudobulks")
}
write_tsv(pb_meta, file.path(out, "02_pseudobulk_DE", "pseudobulk_metadata.tsv"))

con <- gzfile(file.path(out, "02_pseudobulk_DE", "pseudobulk_raw_counts.tsv.gz"), "wt")
write.table(cbind(gene = rownames(pb), as.matrix(pb)), con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)

design <- model.matrix(~ patient + state, data = pb_meta)
if (qr(design)$rank != ncol(design)) stop("Design matrix is not full rank")
write_tsv(cbind(sample = rownames(design), as.data.frame(design, check.names = FALSE)),
          file.path(out, "02_pseudobulk_DE", "design_matrix.tsv"))

safe_pdf(file.path(out, "02_pseudobulk_DE", "pseudobulk_QC.pdf"), 8.2, 4.2, {
  par(mfrow = c(1, 3), mar = c(7, 4, 2, 1))
  cols <- ifelse(pb_meta$state == "High", "#B23A48", "#9AA3A8")
  barplot(pb_meta$cell_count, names.arg = pb_meta$sample, las = 2, col = cols,
          ylab = "Cells", main = "Patient × state cell counts", cex.names = 0.65)
  barplot(pb_meta$total_UMI / 1e6, names.arg = pb_meta$sample, las = 2, col = cols,
          ylab = "Library size (million UMI)", main = "Raw-count library size", cex.names = 0.65)
  barplot(pb_meta$detected_genes, names.arg = pb_meta$sample, las = 2, col = cols,
          ylab = "Detected genes", main = "Pseudobulk gene detection", cex.names = 0.65)
})
writeLines(c(
  paste("Started", format(Sys.time())),
  "Input: raw integer RNA counts from output_step0/Wu2021_TNBC_malignant.rds.",
  "Population: exact 10,836 frozen malignant epithelial barcodes.",
  "State: High=upper tertile; Other=Intermediate+Low.",
  "Aggregation: raw-count sums by patient × state.",
  "Design: ~ patient + state.",
  paste("Design rank:", qr(design)$rank, "of", ncol(design)),
  paste("Completed", format(Sys.time()))
), file.path(out, "02_pseudobulk_DE", "pseudobulk_build_log.txt"))

# ---------------------------------------------------------------------
# Phase 2: score-gene decoupled DESeq2 using filterByExpr fallback.
# ---------------------------------------------------------------------
score_present <- intersect(score_genes, rownames(pb))
pb_decoupled <- pb[!rownames(pb) %in% score_genes, , drop = FALSE]
y_filter <- DGEList(pb_decoupled)
keep <- filterByExpr(y_filter, design = design, min.count = 5)
pb_test <- pb_decoupled[keep, , drop = FALSE]
pb_dense <- round(as.matrix(pb_test))
storage.mode(pb_dense) <- "integer"

exclusion_audit <- data.frame(
  stage = c(
    "raw_RNA_gene_space", "score_genes_present", "after_score_gene_exclusion",
    "after_edgeR_filterByExpr", "score_genes_in_tested_space"
  ),
  n_genes = c(nrow(pb), length(score_present), nrow(pb_decoupled), nrow(pb_test),
              sum(rownames(pb_test) %in% score_genes)),
  rule = c(
    "all genes in RNA counts", "YAP17 union Stem21 intersect RNA counts",
    "remove all 38 frozen score genes before filtering",
    "edgeR::filterByExpr(DGEList, design=~patient+state, min.count=5)",
    "must be zero"
  ),
  status = c("PASS", "PASS", "PASS", "PASS",
             ifelse(any(rownames(pb_test) %in% score_genes), "FAIL", "PASS"))
)
write_tsv(exclusion_audit, file.path(out, "02_pseudobulk_DE", "score_gene_exclusion_audit.tsv"))

dds <- DESeqDataSetFromMatrix(pb_dense, pb_meta, design = ~ patient + state)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("state", "High", "Other"), alpha = 0.05,
               independentFiltering = FALSE)
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
names(res_df) <- c("gene", "baseMean", "log2FC_unshrunk", "lfcSE", "Wald_stat", "pvalue", "padj")

coef_name <- grep("^state_High_vs_Other$", resultsNames(dds), value = TRUE)
if (length(coef_name) == 1L) {
  shr <- tryCatch(
    lfcShrink(dds, coef = coef_name, type = "normal"),
    error = function(e) NULL
  )
} else {
  shr <- NULL
}
if (is.null(shr)) {
  res_df$log2FC_shrunken <- res_df$log2FC_unshrunk
  shrink_method <- "unavailable; unshrunk copied for display"
} else {
  res_df$log2FC_shrunken <- as.data.frame(shr)[res_df$gene, "log2FoldChange"]
  shrink_method <- "DESeq2 lfcShrink type=normal"
}
res_df$contrast <- "High_vs_Other"
res_df$design <- "~ patient + state"
res_df$filter_rule <- "strict_v2 edgeR::filterByExpr"
res_df$score_genes_excluded_before_testing <- TRUE
res_df <- res_df[order(res_df$padj, -res_df$log2FC_unshrunk, na.last = TRUE), ]
write_tsv(res_df, file.path(out, "02_pseudobulk_DE", "deseq2_all_score_decoupled.tsv"))
write_tsv(res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ],
          file.path(out, "02_pseudobulk_DE", "deseq2_significant_FDR05.tsv"))

safe_pdf(file.path(out, "02_pseudobulk_DE", "deseq2_model_diagnostics.pdf"), 8.2, 4.4, {
  par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
  plotMA(res, alpha = 0.05, ylim = c(-5, 5), main = "DESeq2 MA")
  plotDispEsts(dds, main = "Dispersion estimates")
  hist(res_df$pvalue, breaks = 40, col = "#AAB3B7", border = "white",
       main = "Wald P-value distribution", xlab = "P value")
})

# ---------------------------------------------------------------------
# Phase 3: patient directions and LOPO.
# ---------------------------------------------------------------------
norm_counts <- counts(dds, normalized = TRUE)
log_norm <- log2(norm_counts + 1)
delta <- matrix(
  NA_real_, nrow = nrow(log_norm), ncol = length(patient_order),
  dimnames = list(rownames(log_norm), patient_order)
)
for (p in patient_order) {
  delta[, p] <- log_norm[, paste(p, "High", sep = "|")] -
    log_norm[, paste(p, "Other", sep = "|")]
}
direction_wide <- data.frame(gene = rownames(delta), delta, check.names = FALSE)
con <- gzfile(file.path(out, "03_strict142_reconciliation", "patient_direction_matrix.tsv.gz"), "wt")
write.table(direction_wide, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)

direction_summary <- data.frame(
  gene = rownames(delta),
  n_positive_patients = rowSums(delta > 0, na.rm = TRUE),
  n_negative_patients = rowSums(delta < 0, na.rm = TRUE),
  n_zero_patients = rowSums(delta == 0, na.rm = TRUE),
  median_patient_delta = apply(delta, 1, median, na.rm = TRUE),
  min_patient_delta = apply(delta, 1, min, na.rm = TRUE),
  max_patient_delta = apply(delta, 1, max, na.rm = TRUE),
  stringsAsFactors = FALSE
)
write_tsv(direction_summary, file.path(out, "03_strict142_reconciliation", "patient_direction_summary.tsv"))

lopo <- do.call(rbind, lapply(patient_order, function(left_out) {
  d <- delta[, setdiff(patient_order, left_out), drop = FALSE]
  data.frame(
    gene = rownames(d),
    left_out_patient = left_out,
    n_positive_of_7 = rowSums(d > 0, na.rm = TRUE),
    median_delta_of_7 = apply(d, 1, median, na.rm = TRUE),
    positive_direction_retained = apply(d, 1, median, na.rm = TRUE) > 0,
    stringsAsFactors = FALSE
  )
}))
write_tsv(lopo, file.path(out, "03_strict142_reconciliation", "LOPO_gene_stability.tsv"))

# ---------------------------------------------------------------------
# Phase 4: exact historical 142 provenance and strict_v2.
# ---------------------------------------------------------------------
old_ind <- read.delim(old_score_ind_path, check.names = FALSE)
old142 <- unique(old_ind$gene[!is.na(old_ind$adj.P.Val) & old_ind$adj.P.Val < 0.05])
if (length(old142) != 142L) stop("Historical score-independent FDR<0.05 set is not 142")
writeLines(old142, file.path(out, "03_strict142_reconciliation", "historical_old142_genes.txt"))

all_genes <- union(union(old_ind$gene, rownames(pb)), res_df$gene)
rec <- data.frame(gene = all_genes, stringsAsFactors = FALSE)
rec$in_old_142 <- rec$gene %in% old142
rec$tested_in_new_DESeq2 <- rec$gene %in% res_df$gene
rec$passed_expression_filter <- rec$gene %in% rownames(pb_test)
rec$scoring_gene <- rec$gene %in% score_genes
rec <- merge(rec, res_df[, c("gene", "baseMean", "log2FC_unshrunk", "log2FC_shrunken", "pvalue", "padj")],
             by = "gene", all.x = TRUE, sort = FALSE)
rec <- merge(rec, direction_summary, by = "gene", all.x = TRUE, sort = FALSE)
rec$passed_FDR <- !is.na(rec$padj) & rec$padj < 0.05
rec$passed_effect_gate <- !is.na(rec$log2FC_unshrunk) & rec$log2FC_unshrunk >= 0.5
rec$passed_direction_gate <- !is.na(rec$n_positive_patients) & rec$n_positive_patients >= 6
rec$passed_other_historical_gate <- NA
rec$final_strict <- rec$tested_in_new_DESeq2 & !rec$scoring_gene & rec$passed_FDR &
  rec$passed_effect_gate & rec$passed_direction_gate
rec$exclusion_reason <- ifelse(
  rec$final_strict, "included_strict_v2",
  ifelse(rec$scoring_gene, "score_defining_gene_preexcluded",
         ifelse(!rec$passed_expression_filter, "failed_expression_filter_or_absent",
                ifelse(!rec$passed_FDR, "failed_BH_FDR_0.05",
                       ifelse(!rec$passed_effect_gate, "failed_log2FC_ge_0.5",
                              ifelse(!rec$passed_direction_gate, "failed_positive_patients_ge_6",
                                     "other_or_not_tested")))))
)
rec <- rec[order(!rec$final_strict, rec$padj, rec$gene, na.last = TRUE), ]
write_tsv(rec, file.path(out, "03_strict142_reconciliation", "strict_gene_reconciliation.tsv"))

strict_up <- rec$gene[rec$final_strict]
strict_down <- rec$gene[
  rec$tested_in_new_DESeq2 & !rec$scoring_gene & rec$passed_FDR &
    !is.na(rec$log2FC_unshrunk) & rec$log2FC_unshrunk <= -0.5 &
    !is.na(rec$n_negative_patients) & rec$n_negative_patients >= 6
]
writeLines(strict_up, file.path(out, "03_strict142_reconciliation", "final_strict_score_independent_genes.txt"))
writeLines(strict_up, file.path(out, "03_strict142_reconciliation", "final_strict_score_independent_genes_strict_v2.txt"))
writeLines(strict_down, file.path(out, "03_strict142_reconciliation", "final_strict_downregulated_genes_strict_v2.txt"))

cascade <- data.frame(
  order = 1:7,
  stage = c(
    "Raw RNA genes", "After expression filter", "After excluding score genes",
    "BH FDR < 0.05", "High-vs-Other log2FC >= 0.5",
    "Positive in >=6/8 patients", "Final strict_v2 up"
  ),
  n_genes = c(
    nrow(pb), sum(filterByExpr(DGEList(pb), design = design, min.count = 5)),
    nrow(pb_test), sum(rec$tested_in_new_DESeq2 & rec$passed_FDR, na.rm = TRUE),
    sum(rec$tested_in_new_DESeq2 & rec$passed_FDR & rec$passed_effect_gate, na.rm = TRUE),
    sum(rec$tested_in_new_DESeq2 & rec$passed_FDR & rec$passed_effect_gate &
          rec$passed_direction_gate, na.rm = TRUE),
    length(strict_up)
  ),
  stringsAsFactors = FALSE
)
write_tsv(cascade, file.path(out, "03_strict142_reconciliation", "threshold_cascade.tsv"))

p_cascade <- ggplot(cascade, aes(n_genes, reorder(stage, order))) +
  geom_col(fill = "#496B7A", width = 0.70) +
  geom_text(aes(label = n_genes), hjust = -0.15, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Prespecified High-versus-Other strict_v2 cascade",
       subtitle = "Raw-count patient-blocked pseudobulk; all 38 scoring genes excluded",
       x = "Genes retained", y = NULL) +
  theme_classic(base_family = "Arial", base_size = 9)
ggsave(file.path(out, "03_strict142_reconciliation", "threshold_cascade.pdf"),
       p_cascade, width = 180, height = 95, units = "mm", device = cairo_pdf)

overlap <- length(intersect(old142, strict_up))
upset_data <- data.frame(
  set = c("Historical old142", "strict_v2", "Intersection", "Old142 only", "strict_v2 only"),
  n = c(length(old142), length(strict_up), overlap,
        length(setdiff(old142, strict_up)), length(setdiff(strict_up, old142)))
)
p_overlap <- ggplot(upset_data, aes(n, reorder(set, n))) +
  geom_col(aes(fill = set == "Intersection"), width = 0.68, show.legend = FALSE) +
  geom_text(aes(label = n), hjust = -0.15, size = 3) +
  scale_fill_manual(values = c(`TRUE` = "#B23A48", `FALSE` = "#788A93")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(title = "Historical 142 versus current strict_v2",
       subtitle = "Different contrast and different circularity control; not expected to be identical",
       x = "Genes", y = NULL) +
  theme_classic(base_family = "Arial", base_size = 9)
ggsave(file.path(out, "03_strict142_reconciliation", "old142_vs_new_upset.pdf"),
       p_overlap, width = 150, height = 90, units = "mm", device = cairo_pdf)

old_effect <- old_ind[, c("gene", "logFC")]
names(old_effect)[2] <- "historical_limma_logFC_High_vs_Low"
effect_compare <- merge(old_effect, res_df[, c("gene", "log2FC_unshrunk")], by = "gene")
effect_compare$old142 <- effect_compare$gene %in% old142
effect_compare$strict_v2 <- effect_compare$gene %in% strict_up
write_tsv(effect_compare, file.path(out, "03_strict142_reconciliation", "effect_size_comparison.tsv"))
p_effect <- ggplot(effect_compare, aes(historical_limma_logFC_High_vs_Low, log2FC_unshrunk)) +
  geom_hline(yintercept = 0, colour = "#BBBBBB") +
  geom_vline(xintercept = 0, colour = "#BBBBBB") +
  geom_point(colour = "#A8B0B4", size = 0.7, alpha = 0.45) +
  geom_point(data = effect_compare[effect_compare$old142 | effect_compare$strict_v2, ],
             aes(colour = strict_v2), size = 1.2, alpha = 0.8) +
  scale_colour_manual(values = c(`FALSE` = "#63839A", `TRUE` = "#B23A48"), name = "strict_v2") +
  labs(title = "Effect-size comparison across analytical definitions",
       subtitle = "Historical limma High-vs-Low versus DESeq2 High-vs-Other",
       x = "Historical limma log2FC", y = "Current DESeq2 log2FC") +
  theme_classic(base_family = "Arial", base_size = 9)
ggsave(file.path(out, "03_strict142_reconciliation", "effect_size_comparison.pdf"),
       p_effect, width = 145, height = 110, units = "mm", device = cairo_pdf)

focus <- rec[match(c("MCL1", "PPP1R15A", "ATF3"), rec$gene), ]
focus_status <- focus[, c(
  "gene", "baseMean", "log2FC_unshrunk", "log2FC_shrunken", "pvalue", "padj",
  "n_positive_patients", "median_patient_delta", "passed_FDR",
  "passed_effect_gate", "passed_direction_gate", "final_strict", "exclusion_reason"
)]
write_tsv(focus_status, file.path(out, "03_strict142_reconciliation", "MCL1_PPP1R15A_ATF3_gate_status.tsv"))

report <- c(
  "# Historical 142 versus High-versus-Other strict_v2",
  "",
  "## How the old 142 was produced",
  "",
  "The exact generation script was recovered at `0703_rebuild/scripts/12_build_continuous_state_figures.R`. It selected frozen Low and High cells only, summed raw RNA counts by patient×state, used `edgeR::filterByExpr(group, min.count=5)`, TMM normalization, voom and limma with `~patient+group`, and called 142 score-independent genes at BH FDR < 0.05 after downstream score-gene removal.",
  "",
  "## Why the current result is not forced to 142",
  "",
  "The prespecified current analysis changes the contrast to High versus Other (Intermediate+Low), excludes all 38 scoring genes before filtering/model construction, and uses DESeq2. Therefore the historical 142 is not an exact target and threshold tuning to recover 142 would be invalid.",
  "",
  paste0("- Historical set: ", length(old142), " genes."),
  paste0("- Current strict_v2 upregulated set: ", length(strict_up), " genes."),
  paste0("- Overlap: ", overlap, " genes."),
  paste0("- Historical-only: ", length(setdiff(old142, strict_up)), " genes."),
  paste0("- strict_v2-only: ", length(setdiff(strict_up, old142)), " genes."),
  paste0("- Score genes in strict_v2: ", sum(strict_up %in% score_genes), " (must be zero)."),
  "",
  "## Recommended naming",
  "",
  "Use **High-versus-Other score-decoupled strict_v2** for the current set. Retain **historical High-versus-Low 142-gene set** only as provenance/sensitivity. The two sets are not interchangeable.",
  "",
  "## Representative genes",
  "",
  paste(apply(focus_status, 1, function(x) paste0(
    "- ", x[["gene"]], ": log2FC=", signif(as.numeric(x[["log2FC_unshrunk"]]), 3),
    ", BH FDR=", signif(as.numeric(x[["padj"]]), 3),
    ", positive patients=", x[["n_positive_patients"]], "/8, final_strict=", x[["final_strict"]],
    "."
  )), collapse = "\n")
)
writeLines(report, file.path(out, "03_strict142_reconciliation", "reconciliation_report.md"), useBytes = TRUE)

# ---------------------------------------------------------------------
# Phase 5: analysis-only figures.
# ---------------------------------------------------------------------
strict_plot <- head(strict_up, 40)
if (length(strict_plot) < 10L) {
  strict_plot <- head(rec$gene[rec$tested_in_new_DESeq2 & rec$passed_FDR &
                                rec$log2FC_unshrunk > 0], 40)
}
heat_genes <- unique(c("MCL1", "PPP1R15A", "ATF3", strict_plot))
heat_genes <- intersect(heat_genes, rownames(log_norm))
z_by_gene <- t(scale(t(log_norm[heat_genes, , drop = FALSE])))
z_by_gene[!is.finite(z_by_gene)] <- 0

safe_pdf(file.path(out, "03_strict142_reconciliation", "Fig_state_heatmap_analysis.pdf"), 9.2, 7.2, {
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    ann <- pb_meta[, c("patient", "state"), drop = FALSE]
    pheatmap::pheatmap(z_by_gene, cluster_cols = FALSE, cluster_rows = TRUE,
                      annotation_col = ann, show_colnames = TRUE, fontsize = 7,
                      main = "High vs Other score-decoupled genes\nn=8 patients; patient-blocked raw-count pseudobulk")
  } else {
    heatmap(z_by_gene, Colv = NA, scale = "none")
  }
})

delta_long <- data.frame(
  gene = rep(heat_genes, each = length(patient_order)),
  patient = rep(patient_order, times = length(heat_genes)),
  delta = as.vector(t(delta[heat_genes, patient_order, drop = FALSE])),
  stringsAsFactors = FALSE
)
top_effect <- head(heat_genes, 18)
p_pair <- ggplot(delta_long[delta_long$gene %in% top_effect, ],
                 aes(delta, reorder(gene, delta, median), colour = patient)) +
  geom_vline(xintercept = 0, colour = "#BDBDBD") +
  geom_point(size = 1.6, alpha = 0.78) +
  stat_summary(aes(group = gene), fun = median, geom = "point", shape = 23,
               fill = "white", colour = "black", size = 2.2) +
  labs(title = "Patient-resolved effects for top strict genes",
       subtitle = "High−Other log2 normalized pseudobulk counts; diamonds are patient medians",
       x = "Patient-specific High−Other", y = NULL, colour = "Patient") +
  theme_classic(base_family = "Arial", base_size = 8) +
  theme(legend.position = "bottom")
ggsave(file.path(out, "03_strict142_reconciliation", "Fig_patient_effects_analysis.pdf"),
       p_pair, width = 180, height = 125, units = "mm", device = cairo_pdf)

focus_long <- delta_long[delta_long$gene %in% c("MCL1", "PPP1R15A", "ATF3"), ]
p_focus <- ggplot(focus_long, aes(gene, delta, group = patient, colour = patient)) +
  geom_hline(yintercept = 0, colour = "#BBBBBB") +
  geom_line(alpha = 0.35) +
  geom_point(size = 2) +
  stat_summary(aes(group = gene), fun = median, geom = "point", shape = 23,
               fill = "white", colour = "black", size = 3) +
  labs(title = "Predefined representative genes",
       subtitle = "n=8 patients; High vs Other; score genes excluded; patient-blocked pseudobulk",
       x = NULL, y = "High−Other log2 normalized count", colour = "Patient") +
  theme_classic(base_family = "Arial", base_size = 9) +
  theme(legend.position = "bottom")
ggsave(file.path(out, "03_strict142_reconciliation", "Fig_MCL1_PPP1R15A_ATF3_analysis.pdf"),
       p_focus, width = 150, height = 105, units = "mm", device = cairo_pdf)

vol <- res_df
vol$category <- ifelse(vol$gene %in% strict_up, "strict_v2 up",
                       ifelse(!is.na(vol$padj) & vol$padj < 0.05, "FDR<0.05 other", "not significant"))
vol$label <- ifelse(vol$gene %in% c("MCL1", "PPP1R15A", "ATF3") |
                      vol$gene %in% head(strict_up, 8), vol$gene, "")
p_vol <- ggplot(vol, aes(log2FC_unshrunk, -log10(padj))) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = 2, colour = "#AAAAAA") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "#AAAAAA") +
  geom_point(aes(colour = category), size = 0.8, alpha = 0.62) +
  geom_text_repel(aes(label = label), seed = 20260727, size = 2.7,
                  max.overlaps = Inf, min.segment.length = 0) +
  scale_colour_manual(values = c("strict_v2 up" = "#B23A48", "FDR<0.05 other" = "#527C93",
                                 "not significant" = "#C5CACC")) +
  labs(title = "High-versus-Other score-decoupled DESeq2",
       subtitle = "n=8 patients; ~patient+state; raw-count pseudobulk; 38 scoring genes pre-excluded",
       x = "Unshrunk log2FC (High vs Other)", y = expression(-log[10]("BH FDR")), colour = NULL) +
  theme_classic(base_family = "Arial", base_size = 9) +
  theme(legend.position = "bottom")
ggsave(file.path(out, "03_strict142_reconciliation", "Fig_score_decoupled_volcano_analysis.pdf"),
       p_vol, width = 165, height = 115, units = "mm", device = cairo_pdf)

# ---------------------------------------------------------------------
# Phase 6: pathway refresh using complete signed Wald ranking.
# ---------------------------------------------------------------------
rank_vec <- res_df$Wald_stat
names(rank_vec) <- res_df$gene
rank_vec <- sort(rank_vec[is.finite(rank_vec) & !duplicated(names(rank_vec))], decreasing = TRUE)

empty_enrichment <- function(reason) {
  data.frame(ID = NA, Description = NA, setSize = NA, enrichmentScore = NA,
             NES = NA, pvalue = NA, p.adjust = NA, qvalue = NA,
             status = "NOT_RUN", reason = reason, stringsAsFactors = FALSE)
}

run_gsea <- function(t2g, label) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    return(empty_enrichment("clusterProfiler unavailable"))
  }
  ans <- tryCatch(
    clusterProfiler::GSEA(rank_vec, TERM2GENE = t2g, pvalueCutoff = 1,
                          minGSSize = 10, maxGSSize = 500, seed = TRUE,
                          verbose = FALSE),
    error = function(e) e
  )
  if (inherits(ans, "error")) return(empty_enrichment(conditionMessage(ans)))
  d <- as.data.frame(ans)
  if (!nrow(d)) return(empty_enrichment(paste(label, "returned zero terms")))
  d$status <- "PASS"
  d$reason <- ""
  d
}

get_msig <- function(collection, subcollection = NULL) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) return(NULL)
  f <- msigdbr::msigdbr
  x <- tryCatch({
    if (is.null(subcollection)) {
      f(species = "Homo sapiens", collection = collection)
    } else {
      f(species = "Homo sapiens", collection = collection, subcollection = subcollection)
    }
  }, error = function(e) NULL)
  if (is.null(x)) {
    x <- tryCatch({
      if (is.null(subcollection)) {
        f(species = "Homo sapiens", category = collection)
      } else {
        f(species = "Homo sapiens", category = collection, subcategory = subcollection)
      }
    }, error = function(e) NULL)
  }
  if (is.null(x)) return(NULL)
  term_col <- intersect(c("gs_name", "gs_id"), names(x))[1]
  gene_col <- intersect(c("gene_symbol", "human_gene_symbol"), names(x))[1]
  unique(x[, c(term_col, gene_col)])
}

hall <- get_msig("H")
gobp <- get_msig("C5", "GO:BP")
if (is.null(gobp)) gobp <- get_msig("C5", "BP")
react <- get_msig("C2", "CP:REACTOME")
if (is.null(react)) react <- get_msig("C2", "REACTOME")

gsea_h <- if (is.null(hall)) empty_enrichment("Hallmark gene sets unavailable") else run_gsea(hall, "Hallmark")
gsea_g <- if (is.null(gobp)) empty_enrichment("GO BP gene sets unavailable") else run_gsea(gobp, "GO BP")
gsea_r <- if (is.null(react)) empty_enrichment("Reactome gene sets unavailable") else run_gsea(react, "Reactome")
write_tsv(gsea_h, file.path(out, "04_pathway_refresh", "GSEA_Hallmark.tsv"))
write_tsv(gsea_g, file.path(out, "04_pathway_refresh", "GSEA_GO_BP.tsv"))
write_tsv(gsea_r, file.path(out, "04_pathway_refresh", "GSEA_Reactome.tsv"))

universe <- rownames(pb_test)
ora <- empty_enrichment("No strict genes or gene-set package unavailable")
if (length(strict_up) && !is.null(gobp) && requireNamespace("clusterProfiler", quietly = TRUE)) {
  ora_ans <- tryCatch(
    clusterProfiler::enricher(strict_up, universe = universe, TERM2GENE = gobp,
                              pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5),
    error = function(e) e
  )
  if (!inherits(ora_ans, "error") && nrow(as.data.frame(ora_ans))) {
    ora <- as.data.frame(ora_ans)
    ora$status <- "PASS"
    ora$reason <- ""
  } else if (inherits(ora_ans, "error")) {
    ora <- empty_enrichment(conditionMessage(ora_ans))
  }
}
write_tsv(ora, file.path(out, "04_pathway_refresh", "ORA_strict_genes.tsv"))

module_eval <- do.call(rbind, lapply(names(modules), function(nm) {
  genes <- modules[[nm]]
  tested <- intersect(genes, res_df$gene)
  r <- res_df[match(tested, res_df$gene), ]
  data.frame(
    module = gsub("_", " ", nm),
    n_frozen = length(genes),
    n_tested = length(tested),
    median_log2FC = if (nrow(r)) median(r$log2FC_unshrunk, na.rm = TRUE) else NA,
    n_FDR05_up = if (nrow(r)) sum(r$padj < 0.05 & r$log2FC_unshrunk > 0, na.rm = TRUE) else 0,
    n_in_strict_v2 = sum(genes %in% strict_up),
    stringsAsFactors = FALSE
  )
}))
write_tsv(module_eval, file.path(out, "04_pathway_refresh", "pathway_comparison_old_vs_new.tsv"))
path_report <- c(
  "# Pathway refresh",
  "",
  "The current High-versus-Other DESeq2 ranking differs by design from the historical High-versus-Low limma ranking, so pathway analysis was refreshed without using Yan to choose genes or pathways.",
  "",
  "- GSEA ranking: complete signed DESeq2 Wald statistic after expression filtering and pre-test score-gene exclusion.",
  "- Collections: Hallmark, GO Biological Process and Reactome.",
  "- ORA query: final strict_v2 upregulated genes.",
  "- ORA universe: all expression-filtered, score-decoupled genes tested by DESeq2.",
  "- Six frozen modules are evaluated separately in `pathway_comparison_old_vs_new.tsv`.",
  "",
  paste0("Hallmark status: ", paste(unique(gsea_h$status), collapse = ","), "."),
  paste0("GO BP status: ", paste(unique(gsea_g$status), collapse = ","), "."),
  paste0("Reactome status: ", paste(unique(gsea_r$status), collapse = ","), "."),
  paste0("ORA status: ", paste(unique(ora$status), collapse = ","), ".")
)
writeLines(path_report, file.path(out, "04_pathway_refresh", "pathway_refresh_report.md"), useBytes = TRUE)

# Save compact computational state for downstream CNV/Yan/QC steps.
saveRDS(list(
  patient_order = patient_order,
  score_genes = score_genes,
  strict_up = strict_up,
  strict_down = strict_down,
  old142 = old142,
  overlap = overlap,
  res = res_df,
  focus = focus_status,
  module_eval = module_eval,
  shrink_method = shrink_method,
  counts_dim = dim(counts),
  pseudobulk_dim = dim(pb),
  tested_genes = rownames(pb_test)
), file.path(out, "00_manifest", "wu_closure_state.rds"))

cat("Completed:", format(Sys.time()), "\n")
