#!/usr/bin/env Rscript
# Purpose: Effect-size threshold and LFC-shrinkage sensitivity
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Supplementary Fig. S1–S4.

options(stringsAsFactors = FALSE, scipen = 999, width = 240, encoding = "UTF-8")
set.seed(20260810)

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "Figure4_effect_size_threshold_sensitivity_audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tmp"), recursive = TRUE, showWarnings = FALSE)

.libPaths(c(
  file.path(root, "0703_rebuild/R_library"),
  file.path(root, "0710_final_evidence_rebuild/vendor"),
  ".software/r-library",
  .libPaths()
))

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(DESeq2)
  library(apeglm)
  library(ggplot2)
  library(patchwork)
  library(digest)
})

paths <- list(
  authoritative = file.path(root, "Figure3_bottom_panels_HI_feasibility_and_plotting/Figure3_I_effect_consistency_source.tsv"),
  primary146 = file.path(root, "0728_sensitivity_and_cnv_closure/02_direction_gate_sensitivity/gate_6of8_primary_genes.txt"),
  direction_distribution = file.path(root, "0728_sensitivity_and_cnv_closure/02_direction_gate_sensitivity/direction_distribution.tsv"),
  meta = file.path(root, "0703_rebuild/continuous_state_main_figures/Wu2021_continuous_state_cell_metadata.csv"),
  object = file.path(root, "output_step0/Wu2021_TNBC_malignant.rds"),
  yap = file.path(root, "0719_frozen_score_definition_audit/signature_gene_lists/YAP_clean__clean_score.tsv"),
  stem = file.path(root, "0719_frozen_score_definition_audit/signature_gene_lists/Stemness_clean__clean_score.tsv"),
  hallmark_gmt = file.path(root, "h.all.v2024.Hs.symbols.gmt"),
  principal_pathways = file.path(root, "Figure4F_cross_dataset_pathway_concordance/Figure4F_pathway_metrics.tsv"),
  palette = file.path(root, "project_semantic_palette.tsv"),
  style_mapping = file.path(out_dir, "panel_style_mapping.tsv"),
  script = file.path(out_dir, "scripts/01_run_effect_size_and_shrinkage_sensitivity.R")
)
missing <- names(paths)[!vapply(paths, file.exists, logical(1))]
if (length(missing)) stop("Missing required input(s): ", paste(missing, collapse = ", "))

sha256 <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
read_tsv <- function(path) read.delim(path, check.names = FALSE, quote = "", comment.char = "", fileEncoding = "UTF-8")
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
collapse_or_none <- function(x) if (length(x)) paste(x, collapse = "/") else "NONE"

input_hash_before <- vapply(paths[names(paths) != "script"], sha256, character(1))
cat("Reading frozen authoritative tables...\n")

auth <- read_tsv(paths$authoritative)
primary146 <- unique(trimws(readLines(paths$primary146, warn = FALSE, encoding = "UTF-8")))
primary146 <- primary146[nzchar(primary146)]
direction_ref <- read_tsv(paths$direction_distribution)
principal <- read_tsv(paths$principal_pathways)

required_auth <- c(
  "gene", "overall_log2FC_High_vs_Other", "lfcSE", "Wald_stat", "p_value", "BH_FDR", "baseMean",
  "positive_patient_count", "negative_patient_count", "zero_patient_count", "program_146"
)
stopifnot(
  nrow(auth) == 17597L,
  !anyDuplicated(auth$gene),
  all(required_auth %in% names(auth)),
  length(primary146) == 146L,
  nrow(direction_ref) == 17597L,
  setequal(auth$gene, direction_ref$gene),
  nrow(principal) == 10L,
  !anyDuplicated(principal$pathway)
)

primary_from_column <- auth$gene[auth$program_146 %in% TRUE]
primary_from_gate <- auth$gene[
  !is.na(auth$BH_FDR) & auth$BH_FDR < 0.05 &
    auth$overall_log2FC_High_vs_Other >= 0.5 &
    auth$positive_patient_count >= 6
]
stopifnot(
  length(primary_from_column) == 146L,
  length(primary_from_gate) == 146L,
  setequal(primary146, primary_from_column),
  setequal(primary146, primary_from_gate)
)

## ---------------------------------------------------------------------
## Rebuild the exact patient-blocked DESeq2 model for apeglm sensitivity.
## This does not overwrite or redefine the frozen unshrunk results.
## ---------------------------------------------------------------------

patient_order <- c("CID4465", "CID4495", "CID44971", "CID44991", "CID4513", "CID4515", "CID4523", "CID3963")
state_order <- c("Other", "High")
meta <- read.csv(paths$meta, check.names = FALSE)
obj <- readRDS(paths$object)
DefaultAssay(obj) <- "RNA"
stopifnot(nrow(meta) == 10836L, all(meta$cell_id %in% colnames(obj)))

score_genes <- unique(c(read_tsv(paths$yap)$gene, read_tsv(paths$stem)$gene))
stopifnot(length(score_genes) == 38L)

counts_sc <- GetAssayData(obj, assay = "RNA", layer = "counts")[, meta$cell_id, drop = FALSE]
state <- factor(ifelse(meta$YS_tertile == "High", "High", "Other"), levels = state_order)
stopifnot(sum(state == "High") == 3612L, sum(state == "Other") == 7224L)

sample_id <- paste(meta$patient, state, sep = "|")
sample_levels <- as.vector(t(outer(patient_order, state_order, paste, sep = "|")))
sample_factor <- factor(sample_id, levels = sample_levels)
stopifnot(!anyNA(sample_factor), all(table(sample_factor) > 0))
mm <- sparse.model.matrix(~ 0 + sample_factor)
colnames(mm) <- sample_levels
pb <- counts_sc %*% mm
colnames(pb) <- sample_levels

bits <- strsplit(sample_levels, "\\|")
sample_info <- data.frame(
  sample = sample_levels,
  patient = factor(vapply(bits, `[`, character(1), 1), levels = patient_order),
  state = factor(vapply(bits, `[`, character(1), 2), levels = state_order),
  stringsAsFactors = FALSE
)
rownames(sample_info) <- sample_info$sample
design_matrix <- model.matrix(~ patient + state, data = sample_info)

pb_score_independent <- pb[!rownames(pb) %in% score_genes, , drop = FALSE]
keep <- filterByExpr(DGEList(pb_score_independent), design = design_matrix, min.count = 5)
pb_testable <- pb_score_independent[keep, , drop = FALSE]
stopifnot(nrow(pb_testable) == 17597L, setequal(rownames(pb_testable), auth$gene))

dense <- round(as.matrix(pb_testable))
storage.mode(dense) <- "integer"
dds <- DESeqDataSetFromMatrix(countData = dense, colData = sample_info, design = ~ patient + state)
cat("Fitting the frozen DESeq2 design...\n")
dds <- DESeq(dds, quiet = TRUE)
coef_name <- grep("^state_High_vs_Other$", resultsNames(dds), value = TRUE)
stopifnot(length(coef_name) == 1L)

unshrunk <- as.data.frame(results(dds, name = coef_name))
unshrunk$gene <- rownames(unshrunk)
idx_u <- match(auth$gene, unshrunk$gene)
stopifnot(!anyNA(idx_u))
unshrunk <- unshrunk[idx_u, , drop = FALSE]

# The authoritative strict_v2 result stores BH adjustment across the complete
# fixed 17,597-gene testable universe. DESeq2's default independent-filtered
# padj is therefore not the quantity used by the frozen gate.
rebuilt_BH_FDR <- p.adjust(unshrunk$pvalue, method = "BH")

max_abs_difference <- c(
  log2FC = max(abs(unshrunk$log2FoldChange - auth$overall_log2FC_High_vs_Other), na.rm = TRUE),
  lfcSE = max(abs(unshrunk$lfcSE - auth$lfcSE), na.rm = TRUE),
  Wald_stat = max(abs(unshrunk$stat - auth$Wald_stat), na.rm = TRUE),
  p_value = max(abs(unshrunk$pvalue - auth$p_value), na.rm = TRUE),
  BH_FDR = max(abs(rebuilt_BH_FDR - auth$BH_FDR), na.rm = TRUE)
)
if (max_abs_difference[["log2FC"]] > 1e-8 || max_abs_difference[["Wald_stat"]] > 1e-8) {
  stop("Rebuilt DESeq2 model does not reproduce the authoritative unshrunk effects/statistics: ",
       paste(names(max_abs_difference), signif(max_abs_difference, 6), collapse = "; "))
}

cat("Running apeglm shrinkage with coef=", coef_name, "...\n", sep = "")
shrunk <- as.data.frame(lfcShrink(dds, coef = coef_name, type = "apeglm"))
shrunk$gene <- rownames(shrunk)
idx_s <- match(auth$gene, shrunk$gene)
stopifnot(!anyNA(idx_s))
shrunk <- shrunk[idx_s, , drop = FALSE]
stopifnot(all(is.finite(shrunk$log2FoldChange)))

## ---------------------------------------------------------------------
## Independent nested unshrunk-effect sensitivity subsets.
## ---------------------------------------------------------------------

thresholds <- c(0.50, 0.75, 1.00, 1.50)
threshold_labels <- c("0.5 (frozen primary)", "0.75", "1.0", "1.5 (large-effect core)")
names(threshold_labels) <- format(thresholds, nsmall = 2)
key_representatives <- c("MCL1", "PPP1R15A", "ATF3", "BIRC3", "TNFAIP3", "TNF", "CXCL2", "ICAM1", "JUN", "FOS", "TACSTD2", "LAMB3", "LAMC2")
stopifnot(all(key_representatives %in% primary146))

subsets <- setNames(lapply(thresholds, function(thr) {
  auth$gene[
    !is.na(auth$BH_FDR) & auth$BH_FDR < 0.05 &
      auth$positive_patient_count >= 6 &
      auth$overall_log2FC_High_vs_Other >= thr
  ]
}), format(thresholds, nsmall = 2))
stopifnot(setequal(subsets[["0.50"]], primary146))
for (i in 2:length(subsets)) stopifnot(all(subsets[[i]] %in% subsets[[i - 1L]]))

threshold_summary <- do.call(rbind, lapply(seq_along(thresholds), function(i) {
  thr <- thresholds[i]
  nm <- format(thr, nsmall = 2)
  genes <- subsets[[nm]]
  reps <- key_representatives[key_representatives %in% genes]
  data.frame(
    threshold = thr,
    threshold_label = unname(threshold_labels[nm]),
    role = ifelse(thr == 0.5, "FROZEN_PRIMARY_DO_NOT_REDEFINE", ifelse(thr == 1.5, "SENSITIVITY_LARGE_EFFECT_CORE", "SENSITIVITY_ONLY")),
    retained_gene_count = length(genes),
    percentage_retained_relative_to_primary146 = 100 * length(genes) / length(primary146),
    overlap_with_primary146_n = length(intersect(genes, primary146)),
    union_with_primary146_n = length(union(genes, primary146)),
    Jaccard_with_primary146 = length(intersect(genes, primary146)) / length(union(genes, primary146)),
    all_retained_are_primary146 = all(genes %in% primary146),
    key_representative_retained_n = length(reps),
    key_representative_total_n = length(key_representatives),
    key_representative_genes_retained = collapse_or_none(reps),
    retained_genes = collapse_or_none(sort(genes)),
    stringsAsFactors = FALSE
  )
}))

## ---------------------------------------------------------------------
## Hallmark competitive ORA on the same fixed 17,597-gene universe.
## Direction is the sign of a Haldane-corrected log2 odds ratio.
## ---------------------------------------------------------------------

parse_gmt <- function(path) {
  z <- strsplit(readLines(path, warn = FALSE, encoding = "UTF-8"), "\t", fixed = TRUE)
  sets <- lapply(z, function(x) unique(x[-c(1, 2)]))
  names(sets) <- vapply(z, `[[`, character(1), 1)
  sets
}
hallmark <- parse_gmt(paths$hallmark_gmt)
stopifnot(length(hallmark) == 50L)
universe <- auth$gene

ora_rows <- list()
for (i in seq_along(thresholds)) {
  thr <- thresholds[i]
  nm <- format(thr, nsmall = 2)
  selected <- subsets[[nm]]
  z <- do.call(rbind, lapply(names(hallmark), function(pathway) {
    geneset <- intersect(unique(hallmark[[pathway]]), universe)
    a <- length(intersect(selected, geneset))
    b <- length(setdiff(selected, geneset))
    c <- length(setdiff(geneset, selected))
    d <- length(setdiff(universe, union(selected, geneset)))
    tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    ft_two <- fisher.test(tab, alternative = "two.sided")
    ft_greater <- fisher.test(tab, alternative = "greater")
    log2_or <- log2(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)))
    data.frame(
      threshold = thr,
      threshold_label = unname(threshold_labels[nm]),
      retained_gene_count = length(selected),
      pathway = pathway,
      pathway_gene_count_in_universe = length(geneset),
      overlap_gene_count = a,
      overlap_genes = collapse_or_none(sort(intersect(selected, geneset))),
      odds_ratio_exact = unname(ft_two$estimate),
      log2_odds_ratio_Haldane = log2_or,
      two_sided_P = ft_two$p.value,
      overrepresentation_P = ft_greater$p.value,
      stringsAsFactors = FALSE
    )
  }))
  z$two_sided_BH_FDR <- p.adjust(z$two_sided_P, method = "BH")
  z$overrepresentation_BH_FDR <- p.adjust(z$overrepresentation_P, method = "BH")
  z$pathway_direction <- ifelse(z$log2_odds_ratio_Haldane > 0, "High-associated", ifelse(z$log2_odds_ratio_Haldane < 0, "Other-associated", "Neutral"))
  ora_rows[[nm]] <- z
}
hallmark_sensitivity <- do.call(rbind, ora_rows)

principal_idx <- match(hallmark_sensitivity$pathway, principal$pathway)
hallmark_sensitivity$principal_Figure4_pathway <- !is.na(principal_idx)
hallmark_sensitivity$principal_Figure4_order <- principal$frozen_order[principal_idx]
hallmark_sensitivity$principal_Figure4_direction <- principal$biological_group[principal_idx]
hallmark_sensitivity$principal_display_pathway <- principal$display_pathway[principal_idx]
primary_dir <- hallmark_sensitivity$pathway_direction[match(
  hallmark_sensitivity$pathway,
  hallmark_sensitivity$pathway[hallmark_sensitivity$threshold == 0.5]
)]
hallmark_sensitivity$primary_0.5_pathway_direction <- primary_dir
hallmark_sensitivity$directional_concordance_with_primary_0.5 <- hallmark_sensitivity$pathway_direction == primary_dir
hallmark_sensitivity$directional_concordance_with_Figure4B <- ifelse(
  hallmark_sensitivity$principal_Figure4_pathway,
  hallmark_sensitivity$pathway_direction == hallmark_sensitivity$principal_Figure4_direction,
  NA
)
hallmark_sensitivity$significant_two_sided_BH_FDR_0.05 <- hallmark_sensitivity$two_sided_BH_FDR < 0.05
hallmark_sensitivity <- hallmark_sensitivity[order(hallmark_sensitivity$threshold, hallmark_sensitivity$pathway), ]

principal_hallmark <- hallmark_sensitivity[hallmark_sensitivity$principal_Figure4_pathway, ]
principal_counts <- do.call(rbind, lapply(thresholds, function(thr) {
  z <- principal_hallmark[principal_hallmark$threshold == thr, ]
  data.frame(
    threshold = thr,
    principal_pathways_direction_concordant_with_primary_n = sum(z$directional_concordance_with_primary_0.5),
    principal_pathways_total_n = nrow(z),
    principal_pathways_direction_concordant_with_Figure4B_n = sum(z$directional_concordance_with_Figure4B),
    principal_pathways_BH_FDR_0.05_n = sum(z$significant_two_sided_BH_FDR_0.05),
    stringsAsFactors = FALSE
  )
}))
threshold_summary <- merge(threshold_summary, principal_counts, by = "threshold", all.x = TRUE, sort = FALSE)
threshold_summary <- threshold_summary[match(thresholds, threshold_summary$threshold), ]

## ---------------------------------------------------------------------
## All-gene and frozen-program apeglm comparison.
## ---------------------------------------------------------------------

gene_shrink <- data.frame(
  row_type = "GENE",
  scope = "GENE_LEVEL",
  gene = auth$gene,
  program_146 = auth$gene %in% primary146,
  BH_FDR = auth$BH_FDR,
  positive_patient_count = auth$positive_patient_count,
  unshrunk_log2FC = auth$overall_log2FC_High_vs_Other,
  apeglm_shrunken_log2FC = shrunk$log2FoldChange,
  apeglm_posterior_SE = shrunk$lfcSE,
  shrinkage_delta = shrunk$log2FoldChange - auth$overall_log2FC_High_vs_Other,
  unshrunk_direction = ifelse(auth$overall_log2FC_High_vs_Other > 0, "Positive", ifelse(auth$overall_log2FC_High_vs_Other < 0, "Negative", "Zero")),
  shrunken_direction = ifelse(shrunk$log2FoldChange > 0, "Positive", ifelse(shrunk$log2FoldChange < 0, "Negative", "Zero")),
  direction_concordant = sign(auth$overall_log2FC_High_vs_Other) == sign(shrunk$log2FoldChange),
  large_effect_core_unshrunk = auth$gene %in% subsets[["1.50"]],
  remains_positive_after_shrinkage = auth$gene %in% subsets[["1.50"]] & shrunk$log2FoldChange > 0,
  remains_ge_1.5_after_shrinkage = auth$gene %in% subsets[["1.50"]] & shrunk$log2FoldChange >= 1.5,
  stringsAsFactors = FALSE
)

summarize_scope <- function(scope_name, use) {
  z <- gene_shrink[use, ]
  large <- z[z$large_effect_core_unshrunk, ]
  data.frame(
    row_type = "SUMMARY",
    scope = scope_name,
    gene = NA_character_,
    program_146 = ifelse(scope_name == "PROGRAM_146", TRUE, NA),
    BH_FDR = NA_real_,
    positive_patient_count = NA_real_,
    unshrunk_log2FC = NA_real_,
    apeglm_shrunken_log2FC = NA_real_,
    apeglm_posterior_SE = NA_real_,
    shrinkage_delta = NA_real_,
    unshrunk_direction = NA_character_,
    shrunken_direction = NA_character_,
    direction_concordant = NA,
    large_effect_core_unshrunk = NA,
    remains_positive_after_shrinkage = NA,
    remains_ge_1.5_after_shrinkage = NA,
    n_genes = nrow(z),
    Spearman_rho_unshrunk_vs_apeglm = cor(z$unshrunk_log2FC, z$apeglm_shrunken_log2FC, method = "spearman", use = "complete.obs"),
    direction_concordant_n = sum(z$direction_concordant, na.rm = TRUE),
    direction_evaluable_n = sum(is.finite(z$unshrunk_log2FC) & is.finite(z$apeglm_shrunken_log2FC)),
    direction_concordance_fraction = mean(z$direction_concordant, na.rm = TRUE),
    large_effect_unshrunk_n = nrow(large),
    large_effect_remains_positive_n = sum(large$apeglm_shrunken_log2FC > 0),
    large_effect_remains_ge_1.5_n = sum(large$apeglm_shrunken_log2FC >= 1.5),
    stringsAsFactors = FALSE
  )
}

shrink_summary <- rbind(
  summarize_scope("ALL_17597", rep(TRUE, nrow(gene_shrink))),
  summarize_scope("PROGRAM_146", gene_shrink$program_146)
)
for (nm in setdiff(names(shrink_summary), names(gene_shrink))) gene_shrink[[nm]] <- NA
for (nm in setdiff(names(gene_shrink), names(shrink_summary))) shrink_summary[[nm]] <- NA
shrink_output <- rbind(shrink_summary[, names(gene_shrink)], gene_shrink)

## ---------------------------------------------------------------------
## Write requested tables before plotting.
## ---------------------------------------------------------------------

threshold_file <- file.path(out_dir, "log2FC_threshold_sensitivity.tsv")
hallmark_file <- file.path(out_dir, "Hallmark_threshold_sensitivity.tsv")
shrink_file <- file.path(out_dir, "LFC_shrinkage_sensitivity.tsv")
write_tsv(threshold_summary, threshold_file)
write_tsv(hallmark_sensitivity, hallmark_file)
write_tsv(shrink_output, shrink_file)

stopifnot(
  nrow(read_tsv(threshold_file)) == 4L,
  nrow(read_tsv(hallmark_file)) == 200L,
  nrow(read_tsv(shrink_file)) == 17599L,
  read_tsv(threshold_file)$retained_gene_count[1] == 146L
)

## ---------------------------------------------------------------------
## Compact Extended Data QC figure.
## ---------------------------------------------------------------------

palette <- read_tsv(paths$palette)
get_col <- function(token, fallback) {
  z <- palette$hex[match(token, palette$token)]
  ifelse(length(z) == 1L && !is.na(z), z, fallback)
}
col_pos <- get_col("Positive", "#D7605C")
col_neg <- get_col("Negative", "#2F6DB3")
col_neutral <- get_col("Neutral", "#F7F7F7")
col_mid <- get_col("Intermediate", "#8B9298")
col_text <- "#3D4650"
col_grid <- "#E2E5E8"
base_family <- "Arial"

theme_pub <- theme_classic(base_size = 7.2, base_family = base_family) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(size = 8.5, face = "bold", hjust = 0, colour = col_text, margin = margin(b = 4)),
    axis.title = element_text(size = 7.1, colour = col_text),
    axis.text = element_text(size = 6.5, colour = col_text),
    axis.line = element_line(linewidth = 0.35, colour = col_text),
    axis.ticks = element_line(linewidth = 0.3, colour = col_text),
    legend.title = element_text(size = 6.5, face = "bold", colour = col_text),
    legend.text = element_text(size = 6.3, colour = col_text),
    plot.margin = margin(3, 4, 3, 4)
  )

p_a_data <- threshold_summary
p_a_data$threshold_factor <- factor(p_a_data$threshold_label, levels = threshold_labels)
p_a <- ggplot(p_a_data, aes(x = threshold_factor, y = retained_gene_count)) +
  geom_col(width = 0.68, fill = col_pos, alpha = 0.88) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", retained_gene_count, percentage_retained_relative_to_primary146)),
            vjust = -0.18, size = 2.15, lineheight = 0.9, family = base_family, colour = col_text) +
  scale_y_continuous(limits = c(0, max(p_a_data$retained_gene_count) * 1.30), expand = c(0, 0)) +
  labs(title = "A  Nested effect-size subsets", x = "Unshrunk log2FC threshold", y = "Retained genes") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 28, hjust = 1), axis.ticks.x = element_blank())

p_b_data <- principal_hallmark
p_b_data$threshold_factor <- factor(p_b_data$threshold_label, levels = threshold_labels)
p_b_data$pathway_factor <- factor(
  p_b_data$principal_display_pathway,
  levels = rev(principal$display_pathway[order(principal$frozen_order)])
)
p_b_data$star <- ifelse(p_b_data$two_sided_BH_FDR < 0.001, "***",
                        ifelse(p_b_data$two_sided_BH_FDR < 0.01, "**",
                               ifelse(p_b_data$two_sided_BH_FDR < 0.05, "*", "")))
or_lim <- max(abs(p_b_data$log2_odds_ratio_Haldane))
p_b <- ggplot(p_b_data, aes(x = threshold_factor, y = pathway_factor, fill = log2_odds_ratio_Haldane)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  geom_text(aes(label = star), size = 2.25, family = base_family, fontface = "bold", colour = col_text) +
  scale_fill_gradient2(low = col_neg, mid = col_neutral, high = col_pos, midpoint = 0,
                       limits = c(-or_lim, or_lim), name = expression(log[2]~"odds ratio")) +
  labs(title = "B  Principal Hallmark directions", x = "Unshrunk log2FC threshold", y = NULL) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 28, hjust = 1),
    axis.text.y = element_text(size = 6.2),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.key.height = grid::unit(8, "mm")
  )

all_sum <- shrink_summary[shrink_summary$scope == "ALL_17597", ]
prog_sum <- shrink_summary[shrink_summary$scope == "PROGRAM_146", ]
plot_scatter <- gene_shrink
plot_scatter$class <- ifelse(plot_scatter$large_effect_core_unshrunk, "Large-effect core",
                             ifelse(plot_scatter$program_146, "Frozen 146", "All testable genes"))
plot_scatter$class <- factor(plot_scatter$class, levels = c("All testable genes", "Frozen 146", "Large-effect core"))
lfc_lim <- ceiling(max(abs(c(plot_scatter$unshrunk_log2FC, plot_scatter$apeglm_shrunken_log2FC))) * 2) / 2
ann <- sprintf("All genes: rho = %.3f, direction = %.1f%%\nFrozen 146: rho = %.3f, direction = %.1f%%",
               all_sum$Spearman_rho_unshrunk_vs_apeglm, 100 * all_sum$direction_concordance_fraction,
               prog_sum$Spearman_rho_unshrunk_vs_apeglm, 100 * prog_sum$direction_concordance_fraction)
p_c <- ggplot() +
  geom_abline(slope = 1, intercept = 0, colour = col_mid, linewidth = 0.4, linetype = "dashed") +
  geom_point(data = plot_scatter[plot_scatter$class == "All testable genes", ],
             aes(unshrunk_log2FC, apeglm_shrunken_log2FC), colour = "#AEB5BA", alpha = 0.22, size = 0.45) +
  geom_point(data = plot_scatter[plot_scatter$class == "Frozen 146", ],
             aes(unshrunk_log2FC, apeglm_shrunken_log2FC), colour = col_pos, alpha = 0.8, size = 0.9) +
  geom_point(data = plot_scatter[plot_scatter$class == "Large-effect core", ],
             aes(unshrunk_log2FC, apeglm_shrunken_log2FC), shape = 21, fill = col_pos, colour = "#8F2F35", stroke = 0.35, size = 1.35) +
  annotate("text", x = -0.96 * lfc_lim, y = 0.96 * lfc_lim, label = ann,
           hjust = 0, vjust = 1, size = 2.15, lineheight = 0.95, family = base_family, colour = col_text) +
  scale_x_continuous(limits = c(-lfc_lim, lfc_lim), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-lfc_lim, lfc_lim), expand = c(0, 0)) +
  coord_equal(clip = "off") +
  labs(title = "C  apeglm shrinkage sensitivity", x = "Unshrunk log2FC", y = "apeglm-shrunken log2FC") +
  theme_pub

left_block <- p_a / p_b + plot_layout(heights = c(0.72, 1.28))
combined <- (left_block | p_c) +
  plot_layout(widths = c(1.15, 1)) +
  plot_annotation(
    title = "Effect-size threshold and LFC-shrinkage sensitivity",
    caption = "The 146-gene program remains the frozen primary definition. Hallmark tiles show competitive ORA log2 odds ratios; * BH FDR < 0.05, ** < 0.01, *** < 0.001.",
    theme = theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(family = base_family, size = 10, face = "bold", hjust = 0, colour = col_text, margin = margin(b = 4)),
      plot.caption = element_text(family = base_family, size = 6.3, hjust = 0, colour = col_text, margin = margin(t = 4)),
      plot.margin = margin(4, 5, 4, 5, unit = "mm")
    )
  )

pdf_file <- file.path(out_dir, "effect_size_threshold_sensitivity_QC.pdf")
png_file <- file.path(out_dir, "effect_size_threshold_sensitivity_QC.png")
ggsave(pdf_file, combined, width = 180, height = 125, units = "mm", device = grDevices::cairo_pdf, bg = "white")
ggsave(png_file, combined, width = 180, height = 125, units = "mm", dpi = 600, bg = "white")
stopifnot(file.exists(pdf_file), file.info(pdf_file)$size > 10000, file.exists(png_file), file.info(png_file)$size > 10000)

## ---------------------------------------------------------------------
## Provenance and hard non-modification checks.
## ---------------------------------------------------------------------

input_hash_after <- vapply(paths[names(paths) != "script"], sha256, character(1))
stopifnot(identical(input_hash_before, input_hash_after))

source_hashes <- data.frame(
  file_role = c(names(input_hash_before), "analysis_script", "threshold_output", "hallmark_output", "shrinkage_output", "qc_pdf", "qc_png"),
  file = c(unname(unlist(paths[names(paths) != "script"])), paths$script, threshold_file, hallmark_file, shrink_file, pdf_file, png_file),
  sha256 = c(unname(input_hash_before), sha256(paths$script), sha256(threshold_file), sha256(hallmark_file), sha256(shrink_file), sha256(pdf_file), sha256(png_file)),
  stringsAsFactors = FALSE
)
write_tsv(source_hashes, file.path(out_dir, "source_hashes.tsv"))

large_summary <- shrink_summary[shrink_summary$scope == "PROGRAM_146", ]
qc_lines <- c(
  "Figure 4 effect-size threshold and apeglm sensitivity QC",
  "=========================================================",
  "DATA_GATE: PASS",
  "Frozen 146-gene program modified: NO",
  "Figure 4B classification modified: NO",
  "Authoritative testable genes: 17597",
  "Formal model: 16 patient-state raw-count pseudobulks; design = ~ patient + state",
  "High/Other: pooled upper-tertile High versus Intermediate + Low",
  "Score-definition genes excluded before filtering/testing: 38",
  paste0("Primary >=0.5 gate exact membership reproduced: ", setequal(subsets[["0.50"]], primary146)),
  paste0("Primary gene count: ", length(primary146)),
  paste0("Threshold retained counts: ", paste(format(thresholds), vapply(subsets, length, integer(1)), sep = "=", collapse = "; ")),
  paste0("DESeq2 coefficient: ", coef_name),
  paste0("DESeq2 version: ", as.character(packageVersion("DESeq2"))),
  paste0("apeglm version: ", as.character(packageVersion("apeglm"))),
  "BH FDR reproduction: direct BH adjustment across all 17,597 authoritative testable genes",
  "DESeq2 default independent-filtered padj used for the frozen gate: NO",
  paste0("Maximum rebuilt-authoritative unshrunk differences: ", paste(names(max_abs_difference), signif(max_abs_difference, 6), sep = "=", collapse = "; ")),
  paste0("All-gene Spearman unshrunk vs shrunken: ", signif(all_sum$Spearman_rho_unshrunk_vs_apeglm, 7)),
  paste0("All-gene direction concordance: ", all_sum$direction_concordant_n, "/", all_sum$direction_evaluable_n, " (", signif(100 * all_sum$direction_concordance_fraction, 6), "%)"),
  paste0("Frozen-146 Spearman unshrunk vs shrunken: ", signif(prog_sum$Spearman_rho_unshrunk_vs_apeglm, 7)),
  paste0("Frozen-146 direction concordance: ", prog_sum$direction_concordant_n, "/", prog_sum$direction_evaluable_n, " (", signif(100 * prog_sum$direction_concordance_fraction, 6), "%)"),
  paste0("Unshrunk >=1.5 large-effect core: ", large_summary$large_effect_unshrunk_n),
  paste0("Large-effect genes remaining positive after shrinkage: ", large_summary$large_effect_remains_positive_n),
  paste0("Large-effect genes remaining >=1.5 after shrinkage: ", large_summary$large_effect_remains_ge_1.5_n),
  "Hallmark method: two-sided competitive Fisher ORA on the fixed 17,597-gene universe; direction is sign of Haldane-corrected log2 odds ratio",
  "Hallmark version: MSigDB Hallmark v2024 local GMT",
  paste0("Hallmark GMT SHA-256: ", input_hash_before[["hallmark_gmt"]]),
  "Plot classification: Supplementary/Extended Data QC candidate",
  "PDF: vector, 180 x 125 mm",
  "PNG: 600 dpi, white background",
  "Clustering: NO",
  "Threshold tuning after results: NO",
  "Primary program redefinition: NO",
  "PLOT_STATUS: GENERATED_PENDING_VISUAL_QC"
)
writeLines(qc_lines, file.path(out_dir, "effect_size_sensitivity_QC.txt"), useBytes = TRUE)

cat("Completed effect-size threshold and apeglm sensitivity audit.\n")
cat("Retained counts: ", paste(vapply(subsets, length, integer(1)), collapse = ", "), "\n", sep = "")
cat("All-gene rho: ", all_sum$Spearman_rho_unshrunk_vs_apeglm, "\n", sep = "")
cat("Program-146 rho: ", prog_sum$Spearman_rho_unshrunk_vs_apeglm, "\n", sep = "")
cat("Large-effect remains >=1.5: ", large_summary$large_effect_remains_ge_1.5_n, "\n", sep = "")
