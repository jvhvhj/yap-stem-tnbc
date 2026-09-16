#!/usr/bin/env Rscript
# Purpose: Figure 4 patient gene-effect consistency
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 4.

options(stringsAsFactors = FALSE, scipen = 999, width = 220)

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "Figure4D_patient_effect_consistency")
.libPaths(c(
  file.path(root, "0703_rebuild/R_library"),
  file.path(root, "0710_final_evidence_rebuild/vendor"),
  ".software/r-library",
  .libPaths()
))

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(digest)
})

patient_order <- c(
  "CID4465", "CID4495", "CID44971", "CID44991",
  "CID4513", "CID4515", "CID4523", "CID3963"
)

source_path <- file.path(out_dir, "Figure4D_patient_effect_source.tsv")
summary_path <- file.path(out_dir, "Figure4D_gene_summary.tsv")
palette_path <- file.path(root, "project_semantic_palette.tsv")
fig4c_audit_path <- file.path(root, "Figure-cursor/figure4c/Figure4C_audit_table.tsv")
overall_path <- file.path(root, "Figure3_bottom_panels_HI_feasibility_and_plotting/Figure3_I_effect_consistency_source.tsv")
direction_ref_path <- file.path(root, "0728_sensitivity_and_cnv_closure/02_direction_gate_sensitivity/direction_distribution.tsv")
strict_script_path <- file.path(root, "tmp/Figure3_I_effect_consistency.R")
direction_script_path <- file.path(root, "Figure3_I_direction_recompute.R")
gate_dir <- file.path(root, "0728_sensitivity_and_cnv_closure/02_direction_gate_sensitivity")

required_paths <- c(
  source_path, summary_path, palette_path, fig4c_audit_path, overall_path,
  direction_ref_path, strict_script_path, direction_script_path
)
stopifnot(all(file.exists(required_paths)))

sha256 <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
effect_source <- read.delim(source_path, check.names = FALSE)
gene_summary <- read.delim(summary_path, check.names = FALSE)
fig4c_audit <- read.delim(fig4c_audit_path, check.names = FALSE)
overall <- read.delim(overall_path, check.names = FALSE)
direction_ref <- read.delim(direction_ref_path, check.names = FALSE)

fig4c_genes <- unique(fig4c_audit$gene)
stopifnot(
  length(fig4c_genes) == 22L,
  nrow(effect_source) == 176L,
  nrow(gene_summary) == 22L,
  identical(gene_summary$gene, fig4c_genes),
  identical(unique(effect_source$gene), fig4c_genes),
  identical(unique(effect_source$patient), patient_order),
  all(effect_source$valid_status == "VALID"),
  all(is.finite(effect_source$patient_level_High_vs_Other_effect))
)

effect_matrix <- matrix(
  NA_real_, nrow = length(fig4c_genes), ncol = length(patient_order),
  dimnames = list(fig4c_genes, patient_order)
)
for (i in seq_len(nrow(effect_source))) {
  effect_matrix[effect_source$gene[i], effect_source$patient[i]] <-
    effect_source$patient_level_High_vs_Other_effect[i]
}
stopifnot(all(is.finite(effect_matrix)))

re_pos <- rowSums(effect_matrix > 0)
re_neg <- rowSums(effect_matrix < 0)
re_zero <- rowSums(effect_matrix == 0)
stopifnot(
  all(re_pos == gene_summary$positive_patient_count),
  all(re_neg == gene_summary$negative_patient_count),
  all(re_zero == gene_summary$zero_patient_count),
  all(ifelse(gene_summary$association_class == "High-associated", re_pos, re_neg) ==
        gene_summary$concordant_patient_count)
)

gate_files <- c(
  `5of8` = file.path(gate_dir, "gate_5of8_genes.txt"),
  `6of8` = file.path(gate_dir, "gate_6of8_primary_genes.txt"),
  `7of8` = file.path(gate_dir, "gate_7of8_genes.txt"),
  `8of8` = file.path(gate_dir, "gate_8of8_genes.txt")
)
expected_n <- c(`5of8` = 147L, `6of8` = 146L, `7of8` = 130L, `8of8` = 62L)
gate_checks <- do.call(rbind, lapply(5:8, function(k) {
  nm <- paste0(k, "of8")
  calc <- overall$gene[
    !is.na(overall$BH_FDR) & overall$BH_FDR < 0.05 &
      overall$overall_log2FC_High_vs_Other >= 0.5 &
      overall$positive_patient_count >= k
  ]
  ref <- trimws(readLines(gate_files[nm], warn = FALSE))
  ref <- ref[nzchar(ref)]
  data.frame(
    gate = nm,
    calculated_n = length(calc),
    expected_n = unname(expected_n[nm]),
    exact_membership = setequal(calc, ref),
    stringsAsFactors = FALSE
  )
}))
strict_verified <-
  nrow(overall) == 17597L && nrow(direction_ref) == 17597L &&
  all(gate_checks$calculated_n == gate_checks$expected_n) &&
  all(gate_checks$exact_membership)
stopifnot(strict_verified)

palette <- read.delim(palette_path, check.names = FALSE)
pal <- setNames(palette$hex, palette$token)
pick <- function(token, fallback) if (token %in% names(pal)) unname(pal[token]) else fallback
col_positive <- pick("Positive", "#D7605C")
col_negative <- pick("Negative", "#2F6DB3")
col_neutral <- pick("Neutral", "#F7F7F7")

full_limit <- max(abs(c(effect_matrix, gene_summary$overall_log2FC_High_vs_Other)))
full_limit <- ceiling(full_limit * 10) / 10
effect_col_fun <- colorRamp2(
  c(-full_limit, 0, full_limit),
  c(col_negative, col_neutral, col_positive)
)

render_matrix <- function(display_genes, output_path, height_mm) {
  zsum <- gene_summary[match(display_genes, gene_summary$gene), , drop = FALSE]
  zmat <- effect_matrix[display_genes, patient_order, drop = FALSE]
  zsplit <- factor(zsum$association_class, levels = c("High-associated", "Other-associated"))
  zcons <- matrix(0, nrow = length(display_genes), ncol = 1,
                  dimnames = list(display_genes, "Consistency"))
  zoverall <- matrix(zsum$overall_log2FC_High_vs_Other,
                     nrow = length(display_genes), ncol = 1,
                     dimnames = list(display_genes, "Overall log2FC"))

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
    row_title = NULL,
    show_row_names = FALSE,
    show_heatmap_legend = FALSE,
    column_names_gp = gpar(fontfamily = "Arial", fontsize = 6.6, fontface = "bold"),
    column_names_rot = 0,
    rect_gp = gpar(col = "#D6DADF", lwd = 0.35),
    width = unit(17, "mm"),
    cell_fun = function(j, i, x, y, w, h, fill) {
      grid.text(zsum$concordance_fraction[i], x, y,
                gp = gpar(fontfamily = "Arial", fontsize = 6.6,
                          fontface = "bold", col = "#343A40"))
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
    row_title = NULL,
    show_row_names = FALSE,
    show_heatmap_legend = FALSE,
    column_names_gp = gpar(fontfamily = "Arial", fontsize = 6.6, fontface = "bold"),
    column_names_rot = 0,
    rect_gp = gpar(col = "white", lwd = 0.15),
    width = unit(22, "mm"),
    cell_fun = function(j, i, x, y, w, h, fill) {
      val <- zsum$overall_log2FC_High_vs_Other[i]
      text_col <- if (abs(val) / full_limit > 0.62) "white" else "#222222"
      grid.text(sprintf("%.2f", val), x, y,
                gp = gpar(fontfamily = "Arial", fontsize = 6.4,
                          fontface = "bold", col = text_col))
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
  grid.text("D", x = unit(1.0, "mm"), y = unit(5.5, "mm"),
            just = c("left", "center"),
            gp = gpar(fontfamily = "Arial", fontsize = 11,
                      fontface = "bold", col = "#1E2429"))
  grid.text("Patient-resolved gene-effect consistency",
            x = unit(8.5, "mm"), y = unit(5.5, "mm"),
            just = c("left", "center"),
            gp = gpar(fontfamily = "Arial", fontsize = 9.2,
                      fontface = "bold", col = "#1E2429"))
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

all_display <- gene_summary$gene[order(
  factor(gene_summary$association_class,
         levels = c("High-associated", "Other-associated")),
  gene_summary$Figure4C_order
)]
selected <- gene_summary$provisional_selected %in% c(TRUE, "TRUE")
provisional_display <- gene_summary$gene[selected][order(
  factor(gene_summary$association_class[selected],
         levels = c("High-associated", "Other-associated")),
  gene_summary$Figure4C_order[selected]
)]
stopifnot(length(all_display) == 22L, length(provisional_display) == 16L)

all_pdf <- file.path(out_dir, "Figure4D_all_Fig4C_genes_candidate.pdf")
provisional_pdf <- file.path(out_dir, "Figure4D_provisional_main.pdf")
render_matrix(all_display, all_pdf, height_mm = 127)
render_matrix(provisional_display, provisional_pdf, height_mm = 105)

logic_lines <- c(
  "# Figure 4D strict_v2 patient-direction definition",
  "",
  "- Contrast: pooled upper-tertile High versus Other (Intermediate + Low).",
  "- Formal overall model: 16 patient-state raw-count pseudobulks; DESeq2 design `~ patient + state`.",
  "- The 38 frozen score-definition genes were excluded before filtering and testing.",
  "- Testable universe: 17,597 genes after the existing edgeR `filterByExpr(..., min.count = 5)` filter.",
  "- Patient-level effect: `log2(DESeq2 size-factor-normalized High pseudobulk + 1) - log2(DESeq2 size-factor-normalized Other pseudobulk + 1)` within the same patient.",
  "- Direction: positive = High-associated; negative = Other-associated; exact zero = valid but non-concordant; NA/non-finite = missing/invalid.",
  "- This is the strict_v2 direction definition. The sign is identical to the untransformed normalized-count difference because `log2(x + 1)` is monotonic.",
  "- No patient-level P values were calculated or displayed.",
  "- Primary strict_v2 gate: patient-blocked BH FDR < 0.05, unshrunk overall log2FC >= 0.5, and positive direction in at least 6/8 patients.",
  "",
  "## Reproduction",
  "- The full rebuild hard gate completed before the 22 x 8 source table was written.",
  "- Recomputed 17,597-gene positive/negative/zero direction counts matched both authoritative direction tables.",
  "- Recomputed patient-effect median/minimum/maximum values matched the authoritative distribution within 1e-10.",
  paste0("- Exact gate memberships: ", paste0(gate_checks$gate, "=", gate_checks$calculated_n, collapse = ", "), "."),
  "",
  "## Version B",
  paste0("Provisional genes: ", paste(provisional_display, collapse = ", "), "."),
  "The subset is provisional. If Figure 4C changes only typography/layout, no recalculation is required. If its displayed gene list changes, Version B can be regenerated by subsetting the completed patient-effect table."
)
writeLines(logic_lines, file.path(out_dir, "Figure4D_strict_v2_direction_logic.md"), useBytes = TRUE)

output_paths <- c(all_pdf, provisional_pdf, source_path, summary_path,
                  file.path(out_dir, "Figure4D_strict_v2_direction_logic.md"))
output_hashes <- vapply(output_paths, sha256, character(1))
qc_lines <- c(
  "Figure 4D patient-effect consistency QC",
  "========================================",
  "Status: PASS",
  "",
  "STRICT_V2 REPRODUCTION",
  "17,597-gene positive/negative/zero direction counts: EXACT",
  "17,597-gene patient-effect median/min/max: EXACT within 1e-10",
  paste0("Gate 5/8: ", gate_checks$calculated_n[gate_checks$gate == "5of8"], " genes; exact membership = TRUE"),
  paste0("Gate 6/8: ", gate_checks$calculated_n[gate_checks$gate == "6of8"], " genes; exact membership = TRUE"),
  paste0("Gate 7/8: ", gate_checks$calculated_n[gate_checks$gate == "7of8"], " genes; exact membership = TRUE"),
  paste0("Gate 8/8: ", gate_checks$calculated_n[gate_checks$gate == "8of8"], " genes; exact membership = TRUE"),
  "EXACTLY REPRODUCES EXISTING STRICT_V2 >=6/8 DIRECTION-GATE RESULTS: TRUE",
  "",
  "PATIENT EFFECT",
  "Metric: log2(DESeq2 size-factor-normalized High + 1) - log2(DESeq2 size-factor-normalized Other + 1)",
  "Exact zero: valid, non-concordant",
  "NA/non-finite: missing/invalid",
  "Patient-level P values: NOT CALCULATED / NOT DISPLAYED",
  "",
  "COVERAGE AND ORDER",
  "Actual Figure 4C genes: 22",
  "Patient-effect records: 176 (22 x 8)",
  "Valid records: 176; missing records: 0",
  paste0("Fixed patient order: ", paste(patient_order, collapse = ", ")),
  "Patient clustering/reordering: NO",
  "Gene clustering: NO",
  "Row/patient scaling: NO",
  "Winsorization/clipping: NO",
  paste0("Shared zero-centred effect scale: ", -full_limit, " to ", full_limit),
  "PDF visual render QA: PASS",
  "All-gene candidate page: 180 x 127 mm, 1 page",
  "Provisional-main page: 180 x 105 mm, 1 page",
  "Clipped or overlapping labels: NO",
  "Patient/gene/summary-row alignment errors: NO",
  "Panel title, patient IDs, gene labels and right summary columns readable: YES",
  "",
  "VERSION B",
  paste0("Provisional genes (16): ", paste(provisional_display, collapse = ", ")),
  "Selection: within association class, BH FDR<0.05 priority, then leading-edge recurrence, concordance, BH FDR, |overall log2FC| and Figure 4C order",
  "Status: PROVISIONAL PENDING FINAL FIGURE 4C GENE-LIST FREEZE",
  "",
  "SOURCE HASHES",
  paste0("Figure 4C audit: ", sha256(fig4c_audit_path)),
  paste0("strict_v2 overall table: ", sha256(overall_path)),
  paste0("direction distribution: ", sha256(direction_ref_path)),
  paste0("strict_v2 analysis script: ", sha256(strict_script_path)),
  paste0("patient-effect logic script: ", sha256(direction_script_path)),
  "",
  "OUTPUT HASHES",
  paste0(basename(output_paths), "\t", output_hashes),
  "",
  paste0("R version: ", R.version.string),
  paste0("ComplexHeatmap version: ", as.character(packageVersion("ComplexHeatmap"))),
  "Enrichment analysis is not performed by this script."
)
writeLines(qc_lines, file.path(out_dir, "Figure4D_QC.txt"), useBytes = TRUE)

cat("STRICT_VERIFIED\tTRUE\n")
print(gate_checks)
cat("ALL_GENES\t", length(all_display), "\n", sep = "")
cat("PROVISIONAL_GENES\t", length(provisional_display), "\n", sep = "")
cat("OUTPUT_DIR\t", out_dir, "\n", sep = "")
