#!/usr/bin/env Rscript
# Program146 derivation: portable copy of the corrected July 27 original.
# Scientific calculations are retained from the provenance archive. See README.md.
# --counts is the ORIGINAL Seurat RDS with the raw RNA/counts layer, not a
# normalized matrix. No state assignment, signature scoring or CNV is run here.

options(stringsAsFactors = FALSE, width = 220)
set.seed(20260727)

usage <- function() cat(paste(
  "Usage: Rscript analysis/04_program146/scripts/01_derive_program146.R",
  "  --counts <original_Seurat_RDS_with_RNA_counts>",
  "  --cell-metadata <frozen_cell_state_metadata.csv>",
  "  --historical-score-independent <historical_High_Low_score_independent.tsv>",
  "  --output-dir <NEW_isolated_output_directory>",
  "  [--yap-signature config/signatures/YAP17.tsv]",
  "  [--stem-signature config/signatures/Stem21.tsv]",
  "  [--dry-run]",
  "Run from the repository root. Dry-run validates small inputs only; it does not",
  "load the Seurat object, load analysis packages, calculate statistics, or write outputs.",
  sep = "\n"), "\n")

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) { usage(); quit(status = 0L) }
flags <- list()
allowed <- c("counts", "cell-metadata", "historical-score-independent",
             "output-dir", "yap-signature", "stem-signature", "dry-run")
i <- 1L
while (i <= length(args)) {
  key <- sub("^--", "", args[i])
  if (!startsWith(args[i], "--") || !key %in% allowed || key %in% names(flags))
    stop("Unknown or duplicate option: ", args[i])
  if (key == "dry-run") { flags[[key]] <- TRUE; i <- i + 1L; next }
  if (i == length(args) || startsWith(args[i + 1L], "--"))
    stop("Missing value for --", key)
  flags[[key]] <- args[i + 1L]
  i <- i + 2L
}
required <- c("counts", "cell-metadata", "historical-score-independent", "output-dir")
if (!all(required %in% names(flags))) {
  usage(); stop("Required options missing: ", paste(setdiff(required, names(flags)), collapse = ", "))
}
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
signature_dir <- file.path(root, "config", "signatures")
if (!file.exists(file.path(signature_dir, "Program146.tsv")))
  stop("Run this command from the publication repository root.")
obj_path <- normalizePath(flags[["counts"]], winslash = "/", mustWork = TRUE)
meta_path <- normalizePath(flags[["cell-metadata"]], winslash = "/", mustWork = TRUE)
old_score_ind_path <- normalizePath(flags[["historical-score-independent"]], winslash = "/", mustWork = TRUE)
yap_path <- if (is.null(flags[["yap-signature"]])) file.path(signature_dir, "YAP17.tsv") else flags[["yap-signature"]]
stem_path <- if (is.null(flags[["stem-signature"]])) file.path(signature_dir, "Stem21.tsv") else flags[["stem-signature"]]
yap_path <- normalizePath(yap_path, winslash = "/", mustWork = TRUE)
stem_path <- normalizePath(stem_path, winslash = "/", mustWork = TRUE)
input_paths <- c(obj_path, meta_path, old_score_ind_path, yap_path, stem_path)
if (any(file.info(input_paths)$isdir)) stop("Every input must be a file.")

# Reject reuse even of an empty directory. The parent must already exist.
# Never write under the source repository or into a historical output directory.
requested_out <- path.expand(flags[["output-dir"]])
if (file.exists(requested_out) || dir.exists(requested_out))
  stop("Output path already exists; choose a NEW isolated directory. Nothing was overwritten.")
parent_out <- normalizePath(dirname(requested_out), winslash = "/", mustWork = TRUE)
out <- file.path(parent_out, basename(requested_out))
if (out == root || startsWith(paste0(out, "/"), paste0(root, "/")))
  stop("Choose an output directory outside the source repository.")

read_gene_file <- function(path) {
  x <- read.delim(path, check.names = FALSE)
  col <- intersect(c("gene", "Gene", "symbol", "SYMBOL"), names(x))[1]
  if (is.na(col)) stop("No gene column: ", path)
  unique(trimws(as.character(x[[col]])))
}
write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  cat("Wrote", path, nrow(x), "rows\n")
}

# Added contract checks do not change eligible cells, gene definitions or order.
input_meta <- read.csv(meta_path, check.names = FALSE)
metadata_fields <- c("cell_id", "patient", "YS_tertile", "YAP_score", "Stemness_score", "YAP_Stem_Score")
if (!all(metadata_fields %in% names(input_meta))) stop("Required frozen metadata fields missing.")
expected_patients <- c("CID3963", "CID4465", "CID4495", "CID44971", "CID44991", "CID4513", "CID4515", "CID4523")
if (nrow(input_meta) != 10836L || anyNA(input_meta$cell_id) || any(input_meta$cell_id == "") ||
    anyDuplicated(input_meta$cell_id) || anyNA(input_meta$patient) ||
    !setequal(as.character(input_meta$patient), expected_patients))
  stop("Frozen cell/patient identity contract failed.")
if (anyNA(input_meta$YS_tertile) || any(input_meta$YS_tertile == ""))
  stop("State labels must already exist and be nonmissing; no assignment is inferred.")
for (nm in c("YAP17", "Stem21")) {
  supplied_path <- if (nm == "YAP17") yap_path else stem_path
  supplied <- read.delim(supplied_path, check.names = FALSE)
  g <- read_gene_file(supplied_path)
  ref <- read_gene_file(file.path(signature_dir, paste0(nm, ".tsv")))
  if (nrow(supplied) != length(g) || anyNA(g) || any(g == "") || !identical(g, ref))
    stop("Signature must exactly match the fixed repository file, including gene order: ", nm)
}
if (length(read_gene_file(yap_path)) != 17L || length(read_gene_file(stem_path)) != 21L ||
    length(union(read_gene_file(yap_path), read_gene_file(stem_path))) != 38L)
  stop("YAP17/Stem21/union38 contract failed.")
historical_check <- read.delim(old_score_ind_path, check.names = FALSE)
if (!all(c("gene", "adj.P.Val") %in% names(historical_check)))
  stop("Historical reconciliation input must contain gene and adj.P.Val.")
if (anyNA(historical_check$gene) || any(historical_check$gene == ""))
  stop("Historical reconciliation gene identifiers must be nonmissing.")
if (length(unique(historical_check$gene[!is.na(historical_check$adj.P.Val) & historical_check$adj.P.Val < 0.05])) != 142L)
  stop("Historical reconciliation input does not contain the original 142-gene comparison set.")
reference_program <- read_gene_file(file.path(signature_dir, "Program146.tsv"))
if (length(reference_program) != 146L || any(reference_program %in% union(read_gene_file(yap_path), read_gene_file(stem_path))))
  stop("Repository Program146 definition failed integrity checks.")
if (isTRUE(flags[["dry-run"]])) {
  cat("INPUT_PREFLIGHT_ONLY: 10836 unique cells; 8 expected patients; fixed 17/21/146 signatures.\n")
  cat("Seurat expression object NOT LOADED. No aggregation, DE, scoring, statistics or files generated.\n")
  quit(status = 0L)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(DESeq2)
})
if (!dir.create(out, recursive = FALSE, showWarnings = FALSE))
  stop("Could not exclusively create the fresh output directory.")
input_manifest <- data.frame(
  role = c("RNA_counts_Seurat_RDS", "frozen_cell_metadata", "historical_score_independent_comparison", "YAP17", "Stem21"),
  path = input_paths, bytes = file.info(input_paths)$size,
  md5 = unname(tools::md5sum(input_paths)), stringsAsFactors = FALSE
)
write_tsv(input_manifest, file.path(out, "input_manifest.tsv"))
cat("Started:", format(Sys.time()), "\n")
cat("R:", R.version.string, "\n")

# BEGIN ORIGINAL DERIVATION CALCULATIONS

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

patient_order <- unique(as.character(meta$patient))

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
write_tsv(pb_meta, file.path(out, "pseudobulk_samples.tsv"))

con <- gzfile(file.path(out, "pseudobulk_counts.tsv.gz"), "wt")
write.table(cbind(gene = rownames(pb), as.matrix(pb)), con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)

design <- model.matrix(~ patient + state, data = pb_meta)
if (qr(design)$rank != ncol(design)) stop("Design matrix is not full rank")
write_tsv(cbind(sample = rownames(design), as.data.frame(design, check.names = FALSE)),
          file.path(out, "pseudobulk_design.tsv"))

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
write_tsv(exclusion_audit, file.path(out, "score_gene_exclusion_summary.tsv"))

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
write_tsv(res_df, file.path(out, "Program146_DE_all.tsv"))
write_tsv(res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ],
          file.path(out, "Program146_DE_FDR05.tsv"))

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
con <- gzfile(file.path(out, "patient_log2_normalized_difference.tsv.gz"), "wt")
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
write_tsv(direction_summary, file.path(out, "patient_direction_summary.tsv"))

old_ind <- read.delim(old_score_ind_path, check.names = FALSE)
old142 <- unique(old_ind$gene[!is.na(old_ind$adj.P.Val) & old_ind$adj.P.Val < 0.05])
if (length(old142) != 142L) stop("Historical score-independent FDR<0.05 set is not 142")
writeLines(old142, file.path(out, "historical142_comparison_genes.txt"))

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
write_tsv(rec, file.path(out, "gene_selection_audit.tsv"))

strict_up <- rec$gene[rec$final_strict]
strict_down <- rec$gene[
  rec$tested_in_new_DESeq2 & !rec$scoring_gene & rec$passed_FDR &
    !is.na(rec$log2FC_unshrunk) & rec$log2FC_unshrunk <= -0.5 &
    !is.na(rec$n_negative_patients) & rec$n_negative_patients >= 6
]
writeLines(strict_up, file.path(out, "Program146_genes.txt"))
writeLines(strict_up, file.path(out, "Program146_genes_historical_alias.txt"))
writeLines(strict_down, file.path(out, "recurrent_downregulated_genes.txt"))

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
write_tsv(cascade, file.path(out, "threshold_cascade.tsv"))

focus <- rec[match(c("MCL1", "PPP1R15A", "ATF3"), rec$gene), ]
focus_status <- focus[, c(
  "gene", "baseMean", "log2FC_unshrunk", "log2FC_shrunken", "pvalue", "padj",
  "n_positive_patients", "median_patient_delta", "passed_FDR",
  "passed_effect_gate", "passed_direction_gate", "final_strict", "exclusion_reason"
)]
write_tsv(focus_status, file.path(out, "representative_gene_results.tsv"))

# END ORIGINAL DERIVATION CALCULATIONS

# Added reporting/validation only: the reference is NEVER used to select genes.
score_gene_audit <- data.frame(
  gene = score_genes,
  present_in_raw_counts = score_genes %in% rownames(pb),
  retained_after_score_exclusion = score_genes %in% rownames(pb_decoupled),
  retained_in_tested_universe = score_genes %in% rownames(pb_test),
  stringsAsFactors = FALSE
)
write_tsv(score_gene_audit, file.path(out, "score_gene_exclusion.tsv"))
observed_checks <- c(
  nrow(meta), length(patient_order), ncol(pb), nrow(pb_test),
  sum(!is.na(res_df$padj) & res_df$padj < 0.05), length(strict_up),
  sum(rownames(pb_test) %in% score_genes), as.integer(setequal(strict_up, reference_program))
)
expected_checks <- c(10836, 8, 16, 17597, 362, 146, 0, 1)
checks <- data.frame(
  check = c("cells", "patients", "pseudobulk_samples", "tested_genes", "FDR_significant_genes",
            "Program146_genes", "score_genes_tested", "exact_reference_membership"),
  observed = observed_checks, expected = expected_checks,
  matches_reference = observed_checks == expected_checks
)
write_tsv(checks, file.path(out, "derivation_reference_checks.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
writeLines(paste("Historical display shrinkage:", shrink_method), file.path(out, "shrinkage_note.txt"))
if (!all(checks$matches_reference))
  stop("Derived results differ from the historical reference; preserve outputs for investigation. Do not tune thresholds.")
cat("Completed derivation. Independent existing state assignment and all original primary selection rules retained.\n")
