#!/usr/bin/env Rscript
# Purpose: Patient-blocked High-versus-Other effect verification
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 4.

options(stringsAsFactors = FALSE, width = 220, scipen = 999)
set.seed(20260810)

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "Figure4D_patient_effect_consistency")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(digest)
})

patient_order <- c(
  "CID4465", "CID4495", "CID44971", "CID44991",
  "CID4513", "CID4515", "CID4523", "CID3963"
)
state_order <- c("Other", "High")

paths <- list(
  fig4c_audit = file.path(root, "Figure-cursor/figure4c/Figure4C_audit_table.tsv"),
  strict_overall = file.path(root, "Figure3_bottom_panels_HI_feasibility_and_plotting/Figure3_I_effect_consistency_source.tsv"),
  direction_distribution = file.path(root, "0728_sensitivity_and_cnv_closure/02_direction_gate_sensitivity/direction_distribution.tsv"),
  strict_tiers = file.path(root, "0728_sensitivity_and_cnv_closure/02_direction_gate_sensitivity/strict_gene_tiers.tsv"),
  gate_summary = file.path(root, "0728_sensitivity_and_cnv_closure/02_direction_gate_sensitivity/gate_sensitivity_summary.tsv"),
  gate_5 = file.path(root, "0728_sensitivity_and_cnv_closure/02_direction_gate_sensitivity/gate_5of8_genes.txt"),
  gate_6 = file.path(root, "0728_sensitivity_and_cnv_closure/02_direction_gate_sensitivity/gate_6of8_primary_genes.txt"),
  gate_7 = file.path(root, "0728_sensitivity_and_cnv_closure/02_direction_gate_sensitivity/gate_7of8_genes.txt"),
  gate_8 = file.path(root, "0728_sensitivity_and_cnv_closure/02_direction_gate_sensitivity/gate_8of8_genes.txt"),
  strict_script = file.path(root, "tmp/Figure3_I_effect_consistency.R"),
  direction_script = file.path(root, "Figure3_I_direction_recompute.R"),
  meta = file.path(root, "0703_rebuild/continuous_state_main_figures/Wu2021_continuous_state_cell_metadata.csv"),
  object = file.path(root, "output_step0/Wu2021_TNBC_malignant.rds"),
  yap = file.path(root, "0719_frozen_score_definition_audit/signature_gene_lists/YAP_clean__clean_score.tsv"),
  stem = file.path(root, "0719_frozen_score_definition_audit/signature_gene_lists/Stemness_clean__clean_score.tsv"),
  palette = file.path(root, "project_semantic_palette.tsv")
)
missing_paths <- names(paths)[!vapply(paths, file.exists, logical(1))]
if (length(missing_paths)) {
  stop("Missing authoritative inputs: ", paste(missing_paths, collapse = ", "))
}

sha256 <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
source_hash <- vapply(paths, sha256, character(1))

fig4c <- read.delim(paths$fig4c_audit, check.names = FALSE)
required_fig4c_columns <- c(
  "gene", "strict_v2_log2FC", "strict_v2_FDR", "recurrence_count",
  "patient_positive_n", "patient_total_n"
)
stopifnot(all(required_fig4c_columns %in% names(fig4c)))

fig4c_genes <- unique(fig4c$gene)
expected_text_genes <- c(
  "IL6", "SERPINE1", "JUN", "PLAUR", "CDKN1A", "KLF6", "IRF1",
  "TNFAIP3", "ATF3", "ICAM1", "PPP1R15A", "GADD45A", "TRIB1", "CXCL2",
  "COX5A", "VDAC1", "HPRT1", "PHB2", "PRDX3", "NCBP2", "YWHAQ", "NME4"
)
stopifnot(length(fig4c_genes) == 22L, identical(fig4c_genes, expected_text_genes))

fig4c_gene_audit <- do.call(rbind, lapply(fig4c_genes, function(g) {
  z <- fig4c[fig4c$gene == g, , drop = FALSE]
  fields <- c("strict_v2_log2FC", "strict_v2_FDR", "recurrence_count", "patient_positive_n", "patient_total_n")
  for (f in fields) {
    if (length(unique(z[[f]])) != 1L) stop("Figure 4C audit is inconsistent for ", g, " / ", f)
  }
  data.frame(
    gene = g,
    strict_v2_log2FC = z$strict_v2_log2FC[1],
    strict_v2_FDR = z$strict_v2_FDR[1],
    recurrence_count = z$recurrence_count[1],
    Figure4C_order = match(g, fig4c_genes),
    stringsAsFactors = FALSE
  )
}))

overall <- read.delim(paths$strict_overall, check.names = FALSE)
direction_ref <- read.delim(paths$direction_distribution, check.names = FALSE)
tier_ref <- read.delim(paths$strict_tiers, check.names = FALSE)
gate_summary <- read.delim(paths$gate_summary, check.names = FALSE)

stopifnot(nrow(overall) == 17597L, nrow(direction_ref) == 17597L)
stopifnot(length(unique(overall$gene)) == 17597L, length(unique(direction_ref$gene)) == 17597L)
stopifnot(setequal(overall$gene, direction_ref$gene))
stopifnot(all(fig4c_genes %in% overall$gene))

meta <- read.csv(paths$meta, check.names = FALSE)
obj <- readRDS(paths$object)
DefaultAssay(obj) <- "RNA"
stopifnot(
  nrow(meta) == 10836L,
  all(meta$cell_id %in% colnames(obj)),
  identical(unique(as.character(meta$patient[match(patient_order, meta$patient)])), patient_order)
)

score_genes <- unique(c(
  read.delim(paths$yap, check.names = FALSE)$gene,
  read.delim(paths$stem, check.names = FALSE)$gene
))
stopifnot(length(score_genes) == 38L)

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")[, meta$cell_id, drop = FALSE]
state <- factor(ifelse(meta$YS_tertile == "High", "High", "Other"), levels = state_order)
stopifnot(sum(state == "High") == 3612L, sum(state == "Other") == 7224L)

sample_id <- paste(meta$patient, state, sep = "|")
sample_levels <- as.vector(t(outer(patient_order, state_order, paste, sep = "|")))
sample_factor <- factor(sample_id, levels = sample_levels)
stopifnot(!anyNA(sample_factor), all(table(sample_factor) > 0))
mm <- sparse.model.matrix(~ 0 + sample_factor)
colnames(mm) <- sample_levels
pb <- counts %*% mm
colnames(pb) <- sample_levels

bits <- strsplit(sample_levels, "\\|")
sample_info <- data.frame(
  sample = sample_levels,
  patient = factor(vapply(bits, `[`, character(1), 1), levels = patient_order),
  state = factor(vapply(bits, `[`, character(1), 2), levels = state_order),
  stringsAsFactors = FALSE
)
rownames(sample_info) <- sample_info$sample
design <- model.matrix(~ patient + state, data = sample_info)

pb0 <- pb[!rownames(pb) %in% score_genes, , drop = FALSE]
keep <- filterByExpr(DGEList(pb0), design = design, min.count = 5)
pb1 <- pb0[keep, , drop = FALSE]
stopifnot(nrow(pb1) == 17597L, setequal(rownames(pb1), overall$gene))

dense <- round(as.matrix(pb1))
storage.mode(dense) <- "integer"
dds <- DESeqDataSetFromMatrix(dense, sample_info, design = ~ patient + state)
dds <- estimateSizeFactors(dds)
norm <- counts(dds, normalized = TRUE)

patient_effect <- sapply(patient_order, function(p) {
  log2(norm[, paste(p, "High", sep = "|")] + 1) -
    log2(norm[, paste(p, "Other", sep = "|")] + 1)
})
colnames(patient_effect) <- patient_order

recomputed <- data.frame(
  gene = rownames(patient_effect),
  n_positive_patients = rowSums(patient_effect > 0, na.rm = TRUE),
  n_negative_patients = rowSums(patient_effect < 0, na.rm = TRUE),
  n_zero_patients = rowSums(patient_effect == 0, na.rm = TRUE),
  median_patient_delta = apply(patient_effect, 1, median, na.rm = TRUE),
  min_patient_delta = apply(patient_effect, 1, min, na.rm = TRUE),
  max_patient_delta = apply(patient_effect, 1, max, na.rm = TRUE),
  stringsAsFactors = FALSE
)

idx_ref <- match(recomputed$gene, direction_ref$gene)
idx_overall <- match(recomputed$gene, overall$gene)
stopifnot(!anyNA(idx_ref), !anyNA(idx_overall))

count_match <-
  all(recomputed$n_positive_patients == as.integer(direction_ref$n_positive_patients[idx_ref])) &&
  all(recomputed$n_negative_patients == as.integer(direction_ref$n_negative_patients[idx_ref])) &&
  all(recomputed$n_zero_patients == as.integer(direction_ref$n_zero_patients[idx_ref])) &&
  all(recomputed$n_positive_patients == as.integer(overall$positive_patient_count[idx_overall])) &&
  all(recomputed$n_negative_patients == as.integer(overall$negative_patient_count[idx_overall])) &&
  all(recomputed$n_zero_patients == as.integer(overall$zero_patient_count[idx_overall]))

numeric_diff <- c(
  median = max(abs(recomputed$median_patient_delta - direction_ref$median_patient_delta[idx_ref]), na.rm = TRUE),
  minimum = max(abs(recomputed$min_patient_delta - direction_ref$min_patient_delta[idx_ref]), na.rm = TRUE),
  maximum = max(abs(recomputed$max_patient_delta - direction_ref$max_patient_delta[idx_ref]), na.rm = TRUE)
)
numeric_match <- all(numeric_diff < 1e-10)

gate_files <- list(`5of8` = paths$gate_5, `6of8` = paths$gate_6, `7of8` = paths$gate_7, `8of8` = paths$gate_8)
gate_expected_n <- c(`5of8` = 147L, `6of8` = 146L, `7of8` = 130L, `8of8` = 62L)
gate_checks <- do.call(rbind, lapply(5:8, function(k) {
  gate_name <- paste0(k, "of8")
  criteria <- !is.na(overall$BH_FDR) & overall$BH_FDR < 0.05 &
    overall$overall_log2FC_High_vs_Other >= 0.5 &
    recomputed$n_positive_patients[match(overall$gene, recomputed$gene)] >= k
  calculated <- overall$gene[criteria]
  authoritative <- trimws(readLines(gate_files[[gate_name]], warn = FALSE))
  authoritative <- authoritative[nzchar(authoritative)]
  data.frame(
    gate = gate_name,
    calculated_n = length(calculated),
    expected_n = unname(gate_expected_n[gate_name]),
    authoritative_file_n = length(authoritative),
    exact_membership_match = setequal(calculated, authoritative),
    missing_from_calculated = length(setdiff(authoritative, calculated)),
    extra_in_calculated = length(setdiff(calculated, authoritative)),
    stringsAsFactors = FALSE
  )
}))

gate_summary_match <- all(gate_checks$expected_n == gate_summary$n_genes[match(gate_checks$gate, gate_summary$gate)])
strict_verified <- count_match && numeric_match && gate_summary_match &&
  all(gate_checks$calculated_n == gate_checks$expected_n) &&
  all(gate_checks$authoritative_file_n == gate_checks$expected_n) &&
  all(gate_checks$exact_membership_match)

if (!strict_verified) {
  stop(
    "STRICT_V2 REPRODUCTION FAILED: count_match=", count_match,
    "; numeric_match=", numeric_match,
    "; gate_summary_match=", gate_summary_match,
    "; gate_checks=", paste(gate_checks$exact_membership_match, collapse = ",")
  )
}

effect_22 <- patient_effect[fig4c_genes, patient_order, drop = FALSE]
stopifnot(all(is.finite(effect_22)), nrow(effect_22) == 22L, ncol(effect_22) == 8L)

overall_idx_22 <- match(fig4c_genes, overall$gene)
summary_22 <- merge(
  fig4c_gene_audit,
  data.frame(
    gene = fig4c_genes,
    overall_log2FC_High_vs_Other = overall$overall_log2FC_High_vs_Other[overall_idx_22],
    BH_FDR = overall$BH_FDR[overall_idx_22],
    positive_patient_count = overall$positive_patient_count[overall_idx_22],
    negative_patient_count = overall$negative_patient_count[overall_idx_22],
    zero_patient_count = overall$zero_patient_count[overall_idx_22],
    stringsAsFactors = FALSE
  ),
  by = "gene", sort = FALSE
)
summary_22 <- summary_22[match(fig4c_genes, summary_22$gene), ]

stopifnot(
  max(abs(summary_22$strict_v2_log2FC - summary_22$overall_log2FC_High_vs_Other)) < 1e-12,
  max(abs(summary_22$strict_v2_FDR - summary_22$BH_FDR)) < 1e-12
)

summary_22$association_class <- ifelse(
  summary_22$overall_log2FC_High_vs_Other > 0,
  "High-associated", "Other-associated"
)
summary_22$valid_patient_count <- rowSums(is.finite(effect_22))
summary_22$concordant_patient_count <- ifelse(
  summary_22$association_class == "High-associated",
  rowSums(effect_22 > 0), rowSums(effect_22 < 0)
)
summary_22$concordance_fraction <- paste0(
  summary_22$concordant_patient_count, "/", summary_22$valid_patient_count
)
summary_22$concordance_proportion <- summary_22$concordant_patient_count / summary_22$valid_patient_count

summary_22$provisional_selected <- FALSE
summary_22$provisional_selection_rank_within_class <- NA_integer_
summary_22$provisional_selection_reason <- "Not selected for provisional main-figure subset"

select_within_class <- function(dat, cls, n) {
  z <- dat[dat$association_class == cls, , drop = FALSE]
  z$significant <- !is.na(z$BH_FDR) & z$BH_FDR < 0.05
  z <- z[order(
    -as.integer(z$significant),
    -z$recurrence_count,
    -z$concordant_patient_count,
    z$BH_FDR,
    -abs(z$overall_log2FC_High_vs_Other),
    z$Figure4C_order,
    na.last = TRUE
  ), , drop = FALSE]
  head(z$gene, n)
}

selected_high <- select_within_class(summary_22, "High-associated", 10L)
selected_other <- select_within_class(summary_22, "Other-associated", 6L)
selected_genes <- c(selected_high, selected_other)
stopifnot(length(selected_genes) == 16L, length(unique(selected_genes)) == 16L)

for (cls in c("High-associated", "Other-associated")) {
  genes_cls <- if (cls == "High-associated") selected_high else selected_other
  ix <- match(genes_cls, summary_22$gene)
  summary_22$provisional_selected[ix] <- TRUE
  summary_22$provisional_selection_rank_within_class[ix] <- seq_along(ix)
  summary_22$provisional_selection_reason[ix] <- paste0(
    "Provisional ", cls,
    " ranking: BH FDR<0.05 priority, then recurrence, concordance, BH FDR, |overall log2FC|"
  )
}

effect_source <- do.call(rbind, lapply(seq_along(fig4c_genes), function(i) {
  g <- fig4c_genes[i]
  s <- summary_22[summary_22$gene == g, , drop = FALSE]
  data.frame(
    gene = g,
    Figure4C_order = i,
    patient = patient_order,
    patient_order = seq_along(patient_order),
    patient_level_High_vs_Other_effect = as.numeric(effect_22[g, patient_order]),
    effect_metric = "log2(DESeq2 size-factor-normalized High pseudobulk + 1) - log2(DESeq2 size-factor-normalized Other pseudobulk + 1)",
    effect_direction = ifelse(
      effect_22[g, patient_order] > 0, "High-associated",
      ifelse(effect_22[g, patient_order] < 0, "Other-associated", "Zero")
    ),
    valid_status = ifelse(is.finite(effect_22[g, patient_order]), "VALID", "MISSING"),
    concordant_with_overall_direction = ifelse(
      s$association_class == "High-associated",
      effect_22[g, patient_order] > 0,
      effect_22[g, patient_order] < 0
    ),
    overall_association = s$association_class,
    strict_v2_overall_log2FC = s$overall_log2FC_High_vs_Other,
    strict_v2_BH_FDR = s$BH_FDR,
    Figure4C_leading_edge_recurrence = s$recurrence_count,
    patient_effect_source_file = paths$object,
    patient_effect_source_sha256 = source_hash["object"],
    strict_logic_source_file = paths$direction_script,
    strict_logic_source_sha256 = source_hash["direction_script"],
    stringsAsFactors = FALSE
  )
}))

write.table(
  effect_source,
  file.path(out_dir, "Figure4D_patient_effect_source.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
write.table(
  summary_22,
  file.path(out_dir, "Figure4D_gene_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

palette <- read.delim(paths$palette, check.names = FALSE)
pal <- setNames(palette$hex, palette$token)
pick <- function(token, fallback) if (token %in% names(pal)) unname(pal[token]) else fallback
col_positive <- pick("Positive", "#D7605C")
col_negative <- pick("Negative", "#2F6DB3")
col_neutral <- pick("Neutral", "#F7F7F7")

full_limit <- max(abs(c(effect_22, summary_22$overall_log2FC_High_vs_Other)))
full_limit <- ceiling(full_limit * 10) / 10
effect_col_fun <- colorRamp2(
  c(-full_limit, 0, full_limit),
  c(col_negative, col_neutral, col_positive)
)

render_matrix <- function(display_genes, output_path, height_mm) {
  zsum <- summary_22[match(display_genes, summary_22$gene), , drop = FALSE]
  zmat <- effect_22[display_genes, patient_order, drop = FALSE]
  zsplit <- factor(
    zsum$association_class,
    levels = c("High-associated", "Other-associated")
  )
  zcons <- matrix(0, nrow = length(display_genes), ncol = 1,
                  dimnames = list(display_genes, "Consistency"))
  zoverall <- matrix(
    zsum$overall_log2FC_High_vs_Other,
    nrow = length(display_genes), ncol = 1,
    dimnames = list(display_genes, "Overall log2FC")
  )

  ht_main <- Heatmap(
    zmat,
    name = "High - Other effect",
    col = effect_col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_order = seq_along(display_genes),
    column_order = seq_along(patient_order),
    row_split = zsplit,
    row_gap = unit(1.5, "mm"),
    row_title_rot = 0,
    row_title_side = "left",
    row_title_gp = gpar(fontfamily = "Arial", fontsize = 6.6, fontface = "bold", col = "#343A40"),
    row_names_side = "left",
    row_names_gp = gpar(fontfamily = "Arial", fontsize = 7.0, fontface = "italic", col = "#232629"),
    column_names_gp = gpar(fontfamily = "Arial", fontsize = 6.8, fontface = "bold", col = "#232629"),
    column_names_rot = 45,
    border = FALSE,
    rect_gp = gpar(col = "white", lwd = 0.15),
    width = unit(103, "mm"),
    heatmap_legend_param = list(
      title = "High - Other effect",
      at = c(-full_limit, 0, full_limit),
      labels = format(c(-full_limit, 0, full_limit), trim = TRUE),
      title_gp = gpar(fontfamily = "Arial", fontsize = 6.8, fontface = "bold"),
      labels_gp = gpar(fontfamily = "Arial", fontsize = 6.5),
      legend_width = unit(32, "mm"),
      direction = "horizontal"
    )
  )

  ht_cons <- Heatmap(
    zcons,
    name = "Consistency",
    col = c("0" = "#FFFFFF"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = zsplit,
    row_gap = unit(1.5, "mm"),
    show_row_names = FALSE,
    row_title = NULL,
    show_heatmap_legend = FALSE,
    column_names_gp = gpar(fontfamily = "Arial", fontsize = 6.6, fontface = "bold"),
    column_names_rot = 0,
    rect_gp = gpar(col = "#D6DADF", lwd = 0.35),
    width = unit(17, "mm"),
    cell_fun = function(j, i, x, y, w, h, fill) {
      grid.text(
        zsum$concordance_fraction[i], x, y,
        gp = gpar(fontfamily = "Arial", fontsize = 6.6, fontface = "bold", col = "#343A40")
      )
    }
  )

  ht_overall <- Heatmap(
    zoverall,
    name = "High - Other effect",
    col = effect_col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = zsplit,
    row_gap = unit(1.5, "mm"),
    show_row_names = FALSE,
    row_title = NULL,
    show_heatmap_legend = FALSE,
    column_names_gp = gpar(fontfamily = "Arial", fontsize = 6.6, fontface = "bold"),
    column_names_rot = 0,
    rect_gp = gpar(col = "white", lwd = 0.15),
    width = unit(22, "mm"),
    cell_fun = function(j, i, x, y, w, h, fill) {
      val <- zsum$overall_log2FC_High_vs_Other[i]
      text_col <- if (abs(val) / full_limit > 0.62) "white" else "#222222"
      grid.text(
        sprintf("%.2f", val), x, y,
        gp = gpar(fontfamily = "Arial", fontsize = 6.4, fontface = "bold", col = text_col)
      )
    }
  )

  grDevices::cairo_pdf(
    filename = output_path,
    width = 180 / 25.4,
    height = height_mm / 25.4,
    family = "Arial",
    bg = "white",
    onefile = FALSE
  )
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = 2, ncol = 1,
    heights = unit.c(unit(10, "mm"), unit(1, "null"))
  )))
  pushViewport(viewport(layout.pos.row = 1))
  grid.text(
    "D", x = unit(1.0, "mm"), y = unit(5.5, "mm"), just = c("left", "center"),
    gp = gpar(fontfamily = "Arial", fontsize = 11, fontface = "bold", col = "#1E2429")
  )
  grid.text(
    "Patient-resolved gene-effect consistency",
    x = unit(8.5, "mm"), y = unit(5.5, "mm"), just = c("left", "center"),
    gp = gpar(fontfamily = "Arial", fontsize = 9.2, fontface = "bold", col = "#1E2429")
  )
  popViewport()
  pushViewport(viewport(layout.pos.row = 2))
  draw(
    ht_main + ht_cons + ht_overall,
    newpage = FALSE,
    heatmap_legend_side = "bottom",
    merge_legends = TRUE,
    padding = unit(c(2, 2, 2, 2), "mm")
  )
  popViewport(2)
  dev.off()
}

all_display <- summary_22$gene[order(
  factor(summary_22$association_class, levels = c("High-associated", "Other-associated")),
  summary_22$Figure4C_order
)]
provisional_display <- summary_22$gene[
  summary_22$provisional_selected
][order(
  factor(summary_22$association_class[summary_22$provisional_selected],
         levels = c("High-associated", "Other-associated")),
  summary_22$Figure4C_order[summary_22$provisional_selected]
)]

render_matrix(
  all_display,
  file.path(out_dir, "Figure4D_all_Fig4C_genes_candidate.pdf"),
  height_mm = 127
)
render_matrix(
  provisional_display,
  file.path(out_dir, "Figure4D_provisional_main.pdf"),
  height_mm = 105
)

logic_lines <- c(
  "# Figure 4D strict_v2 patient-direction definition",
  "",
  "- Contrast: pooled upper-tertile High versus Other (Intermediate + Low).",
  "- Biological unit for the formal overall effect: 16 patient-state raw-count pseudobulks; DESeq2 design `~ patient + state`.",
  "- Score-definition genes: the 38 frozen YAP17/Stem21 genes were excluded before expression filtering and testing.",
  "- Testable universe: 17,597 genes after the existing edgeR `filterByExpr(..., min.count = 5)` filter.",
  "- Patient-level direction quantity: `log2(DESeq2 size-factor-normalized High pseudobulk + 1) - log2(DESeq2 size-factor-normalized Other pseudobulk + 1)` for the same patient.",
  "- Direction: positive = High-associated; negative = Other-associated; exact zero = valid but non-concordant; non-finite value = missing/invalid.",
  "- The strict_v2 gate uses the sign of this within-patient quantity. The untransformed normalized-count difference has the identical sign because log2(x + 1) is monotonic.",
  "- No patient-level P values are calculated or displayed.",
  "- Strict primary gate: overall patient-blocked BH FDR < 0.05, unshrunk overall log2FC >= 0.5, and positive direction in at least 6/8 patients.",
  "",
  "## Reproduction result",
  paste0("- Direction counts across all 17,597 genes: ", if (count_match) "EXACT" else "FAILED", "."),
  paste0("- Patient-effect median/minimum/maximum: maximum absolute differences = ", paste(names(numeric_diff), format(numeric_diff, scientific = TRUE), collapse = "; "), "."),
  paste0("- Direction gates: ", paste0(gate_checks$gate, "=", gate_checks$calculated_n, collapse = ", "), "."),
  paste0("- Exact membership reproduction for every gate: ", all(gate_checks$exact_membership_match), "."),
  "",
  "## Authoritative sources",
  paste0("- `", paths$strict_script, "` (SHA-256: `", source_hash["strict_script"], "`)"),
  paste0("- `", paths$direction_script, "` (SHA-256: `", source_hash["direction_script"], "`)"),
  paste0("- `", paths$direction_distribution, "` (SHA-256: `", source_hash["direction_distribution"], "`)"),
  paste0("- `", paths$strict_overall, "` (SHA-256: `", source_hash["strict_overall"], "`)"),
  paste0("- `", paths$fig4c_audit, "` (SHA-256: `", source_hash["fig4c_audit"], "`)"),
  "",
  "Version B is provisional and can be regenerated by subsetting `Figure4D_patient_effect_source.tsv`; no pseudobulk or DE analysis needs to be rerun unless the underlying gene definitions change."
)
writeLines(logic_lines, file.path(out_dir, "Figure4D_strict_v2_direction_logic.md"), useBytes = TRUE)

output_paths <- c(
  file.path(out_dir, "Figure4D_all_Fig4C_genes_candidate.pdf"),
  file.path(out_dir, "Figure4D_provisional_main.pdf"),
  file.path(out_dir, "Figure4D_patient_effect_source.tsv"),
  file.path(out_dir, "Figure4D_gene_summary.tsv"),
  file.path(out_dir, "Figure4D_strict_v2_direction_logic.md")
)
output_hashes <- vapply(output_paths, sha256, character(1))

qc_lines <- c(
  "Figure 4D patient-effect consistency QC",
  "========================================",
  paste0("Status: ", if (strict_verified) "PASS" else "FAIL"),
  "",
  "AUTHORITATIVE DEFINITION",
  "Contrast: High versus Other (Intermediate + Low)",
  "Formal model: 16 patient-state raw-count pseudobulks; design = ~ patient + state",
  "Patient effect: log2(DESeq2 size-factor-normalized High + 1) - log2(DESeq2 size-factor-normalized Other + 1)",
  "Direction: positive / negative / exact zero; exact zero is valid but not concordant; NA/non-finite is missing",
  "Patient-level P values: NOT CALCULATED / NOT DISPLAYED",
  "",
  "STRICT_V2 DIRECTION-GATE REPRODUCTION",
  paste0("17,597-gene direction counts exact: ", count_match),
  paste0("17,597-gene median/min/max exact within 1e-10: ", numeric_match),
  paste0("Maximum absolute numerical differences: ", paste(names(numeric_diff), format(numeric_diff, scientific = TRUE), collapse = "; ")),
  paste0("Gate 5/8: ", gate_checks$calculated_n[gate_checks$gate == "5of8"], " genes; exact membership = ", gate_checks$exact_membership_match[gate_checks$gate == "5of8"]),
  paste0("Gate 6/8: ", gate_checks$calculated_n[gate_checks$gate == "6of8"], " genes; exact membership = ", gate_checks$exact_membership_match[gate_checks$gate == "6of8"]),
  paste0("Gate 7/8: ", gate_checks$calculated_n[gate_checks$gate == "7of8"], " genes; exact membership = ", gate_checks$exact_membership_match[gate_checks$gate == "7of8"]),
  paste0("Gate 8/8: ", gate_checks$calculated_n[gate_checks$gate == "8of8"], " genes; exact membership = ", gate_checks$exact_membership_match[gate_checks$gate == "8of8"]),
  paste0("EXACTLY REPRODUCES EXISTING STRICT_V2 >=6/8 DIRECTION-GATE RESULTS: ", strict_verified),
  "",
  "FIGURE 4C INPUT COVERAGE",
  paste0("Actual Figure 4C genes: ", length(fig4c_genes)),
  paste0("Patient-effect records: ", nrow(effect_source), " (expected 22 x 8 = 176)"),
  paste0("Valid patient-effect records: ", sum(effect_source$valid_status == "VALID")),
  paste0("Missing patient-effect records: ", sum(effect_source$valid_status == "MISSING")),
  paste0("Fixed patient order: ", paste(patient_order, collapse = ", ")),
  "Patient clustering/reordering: NO",
  "Gene clustering: NO",
  "Row/patient scaling: NO",
  "Winsorization/clipping: NO",
  paste0("Shared symmetric effect scale: ", -full_limit, " to ", full_limit),
  "",
  "VERSION B (PROVISIONAL)",
  paste0("Selected genes (", length(provisional_display), "): ", paste(provisional_display, collapse = ", ")),
  "Selection rule: within association class, prioritize BH FDR < 0.05, then leading-edge recurrence, concordance, BH FDR, |overall log2FC|, and Figure 4C order",
  "Status: PROVISIONAL PENDING FINAL FIGURE 4C GENE-LIST FREEZE",
  "If Figure 4C typography/layout changes only: no Figure 4D recalculation required",
  "If Figure 4C displayed gene list changes: regenerate Version B by subsetting the completed full table",
  "",
  "OUTPUTS",
  paste0(basename(output_paths), "\tSHA-256=", output_hashes),
  "",
  paste0("R version: ", R.version.string),
  paste0("DESeq2 version: ", as.character(packageVersion("DESeq2"))),
  paste0("ComplexHeatmap version: ", as.character(packageVersion("ComplexHeatmap"))),
  "Enrichment analysis is not performed by this script."
)
writeLines(qc_lines, file.path(out_dir, "Figure4D_QC.txt"), useBytes = TRUE)

cat("STRICT_VERIFIED\t", strict_verified, "\n", sep = "")
print(gate_checks)
cat("FIG4C_GENES\t", length(fig4c_genes), "\n", sep = "")
cat("PATIENT_EFFECT_RECORDS\t", nrow(effect_source), "\n", sep = "")
cat("PROVISIONAL_GENES\t", paste(provisional_display, collapse = ","), "\n", sep = "")
cat("OUTPUT_DIR\t", out_dir, "\n", sep = "")
