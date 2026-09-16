#!/usr/bin/env Rscript
# Purpose: Figure 4 leading-edge architecture
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 4 and Supplementary Fig. S1.
# =============================================================================
# Figure 4C — Leading-edge gene–pathway architecture of the score-independent
#             YAP–Stem-high program
#
# Directly follows Figure 4B (directional Hallmark landscape): shows WHICH
# concrete score-independent genes drive the Hallmark enrichments shown there.
#
# Layout: structured pathway x gene circle-matrix (no network, no generic heatmap).
#   rows    = 10 selected Hallmark pathways (7 High-enriched, divider, 3 Other-enriched)
#   columns = ~22 selected score-independent leading-edge genes
#   circle  = present only where the gene is a TRUE leadingEdge/core_enrichment member
#             of that pathway; FIXED circle size
#   circle colour = strict_v2 pseudobulk log2FC (High vs Other): warm = positive,
#             cool = negative (project palette #D7605C / #2F6DB3)
#   row-side annotation : NES and BH FDR per pathway
#   column-top annotation: recurrence count (number of selected pathways containing
#             the gene)
#
# Authoritative inputs (all frozen):
#   Figure4_EFG_FINAL_DATA_CLOSURE_ONLY/Figure4_E_GSEA_v2_results.tsv   (leading-edge genes)
#   Figure3_bottom_panels_HI_feasibility_and_plotting/Figure3_I_effect_consistency_source.tsv (strict_v2 DE)
#   0713_score_independent_rebuild/tables/score_gene_list.tsv           (38 score genes, excluded)
# =============================================================================

suppressMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(dplyr)
  library(Cairo)
  library(svglite)
})

DATA_DIR <- "."
OUTDIR   <- file.path(DATA_DIR, "Figure-cursor", "figure4c")
GSEA_T  <- file.path(DATA_DIR, "Figure4_EFG_FINAL_DATA_CLOSURE_ONLY", "Figure4_E_GSEA_v2_results.tsv")
DE_T    <- file.path(DATA_DIR, "Figure3_bottom_panels_HI_feasibility_and_plotting",
                     "Figure3_I_effect_consistency_source.tsv")
SCORE_T <- file.path(DATA_DIR, "0713_score_independent_rebuild/tables", "score_gene_list.tsv")

PDF_OUT   <- file.path(OUTDIR, "Figure4C_leading_edge_gene_pathway_architecture_refined.pdf")
SVG_OUT   <- file.path(OUTDIR, "Figure4C_leading_edge_gene_pathway_architecture_refined.svg")
PNG_OUT   <- file.path(OUTDIR, "Figure4C_leading_edge_gene_pathway_architecture_refined.png")
AUDIT_OUT <- file.path(OUTDIR, "Figure4C_audit_table.tsv")
FULL_OUT  <- file.path(OUTDIR, "Figure4C_full_pathway_leading_edge_table.tsv")

CANVAS_W_MM <- 180
CANVAS_H_MM <- 120
FONT <- "Helvetica"
HIGH_C <- "#D7605C"; OTHER_C <- "#2F6DB3"; NEU_C <- "#F7F7F7"

# -----------------------------------------------------------------------------
# 1. Frozen inputs
# -----------------------------------------------------------------------------
gsea <- read.delim(GSEA_T, stringsAsFactors = FALSE, check.names = FALSE)
de   <- read.delim(DE_T, stringsAsFactors = FALSE, check.names = FALSE)
score<- read.delim(SCORE_T, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(gsea) == 50, nrow(score) == 38)

# selected pathways: 7 High-enriched then 3 Other-enriched (Figure 4B direction)
high_pw <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_INFLAMMATORY_RESPONSE",
             "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_TGF_BETA_SIGNALING",
             "HALLMARK_APOPTOSIS","HALLMARK_HYPOXIA","HALLMARK_APICAL_JUNCTION")
other_pw <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_MYC_TARGETS_V1","HALLMARK_DNA_REPAIR")
sel_pw <- c(high_pw, other_pw)
stopifnot(all(sel_pw %in% gsea$pathway))

# -----------------------------------------------------------------------------
# 2. Candidate genes = union of TRUE leading-edge genes, score genes excluded
# -----------------------------------------------------------------------------
pw_of_gene <- list()
for (i in seq_len(nrow(gsea))) if (gsea$pathway[i] %in% sel_pw)
  for (ge in strsplit(gsea$leading_edge_genes[i], "/")[[1]])
    pw_of_gene[[ge]] <- c(pw_of_gene[[ge]], gsea$pathway[i])
score_genes <- score$gene
cand <- setdiff(names(pw_of_gene), score_genes)
cat("[QA] leading-edge union:", length(names(pw_of_gene)), "-> score-independent:", length(cand), "\n")
de2 <- de[match(cand, de$gene), ]
dat <- data.frame(gene = cand,
  recurrence = sapply(pw_of_gene[cand], length),
  recurrence_high = sapply(pw_of_gene[cand], function(p) sum(p %in% high_pw)),
  recurrence_other = sapply(pw_of_gene[cand], function(p) sum(p %in% other_pw)),
  log2FC = de2$overall_log2FC_High_vs_Other, BH_FDR = de2$BH_FDR,
  pos_pat = de2$positive_patient_count, neg_pat = de2$negative_patient_count,
  stringsAsFactors = FALSE)
dat$total_pat <- dat$pos_pat + dat$neg_pat
dat$patient_consistency <- ifelse(dat$total_pat == 0, 0, dat$pos_pat / dat$total_pat)

# -----------------------------------------------------------------------------
# 3. Gene selection: balanced, reproducible (~22 genes)
#    12 High-shared + 2 High-specific + 6 Other-shared + 2 Other-specific
# -----------------------------------------------------------------------------
rank_by <- function(d) d[order(-d$recurrence, d$BH_FDR, -abs(d$log2FC), -d$patient_consistency), ]
hs  <- rank_by(dat[dat$recurrence_high >= 2, ])
os  <- rank_by(dat[dat$recurrence_other >= 2, ])
hsp <- rank_by(dat[dat$recurrence_high == 1, ])
osp <- rank_by(dat[dat$recurrence_other == 1, ])
sel_g <- c(head(hs$gene, 12), head(hsp$gene, 2), head(os$gene, 6), head(osp$gene, 2))
sel_dat <- dat[dat$gene %in% sel_g, ]
cat("[QA] selected genes:", nrow(sel_dat),
    " (High", sum(sel_dat$recurrence_high > 0), " / Other", sum(sel_dat$recurrence_other > 0), ")\n")

# column order: High-shared, High-specific, Other-shared, Other-specific (each by recurrence desc)
sel_dat$blk <- ifelse(sel_dat$recurrence_high > 0, "High",
               ifelse(sel_dat$recurrence_other > 0, "Other", "High"))
gene_order <- sel_dat[order(ifelse(sel_dat$blk == "High", 0, 1),
                            -sel_dat$recurrence, sel_dat$BH_FDR), "gene"]
# gene block (High-associated genes then Other-associated genes) for a subtle
# vertical separator between the block ending at CXCL2 and the block starting at COX5A
gene_block <- sel_dat$blk[match(gene_order, sel_dat$gene)]

# -----------------------------------------------------------------------------
# 4. Matrix: rows = pathways, cols = genes; value = log2FC where true leading-edge
# -----------------------------------------------------------------------------
PW_DISP <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB"="TNFα–NF-κB",
             "HALLMARK_INFLAMMATORY_RESPONSE"="Inflammatory response",
             "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"="EMT",
             "HALLMARK_TGF_BETA_SIGNALING"="TGF-β",
             "HALLMARK_APOPTOSIS"="Apoptosis",
             "HALLMARK_HYPOXIA"="Hypoxia",
             "HALLMARK_APICAL_JUNCTION"="Apical junction",
             "HALLMARK_OXIDATIVE_PHOSPHORYLATION"="Oxidative phosphorylation",
             "HALLMARK_MYC_TARGETS_V1"="MYC targets V1",
             "HALLMARK_DNA_REPAIR"="DNA repair")
mat <- matrix(NA_real_, nrow = length(sel_pw), ncol = length(gene_order),
              dimnames = list(unname(PW_DISP[sel_pw]), gene_order))
for (i in seq_len(nrow(gsea))) if (gsea$pathway[i] %in% sel_pw) {
  gs <- strsplit(gsea$leading_edge_genes[i], "/")[[1]]
  for (ge in intersect(gs, gene_order)) {
    mat[PW_DISP[gsea$pathway[i]], ge] <- dat$log2FC[dat$gene == ge]
  }
}
cat("[QA] matrix", nrow(mat), "x", ncol(mat), "; filled cells:", sum(!is.na(mat)), "\n")

# -----------------------------------------------------------------------------
# 5. Annotations
# -----------------------------------------------------------------------------
g_sel <- gsea[match(sel_pw, gsea$pathway), ]
nes_vals <- sprintf("%.2f", g_sel$NES)
fdr_vals <- ifelse(g_sel$BH_FDR < 0.001, "<0.001", sprintf("%.1e", g_sel$BH_FDR))
rec_vals <- dat$recurrence[match(gene_order, dat$gene)]

ra <- rowAnnotation(
  "NES" = anno_text(nes_vals, gp = gpar(fontfamily = FONT, fontsize = 6.5, col = "grey15"),
                    location = 0.5, just = "center"),
  "BH FDR" = anno_text(fdr_vals, gp = gpar(fontfamily = FONT, fontsize = 6.5, col = "grey40"),
                       location = 0.5, just = "center"),
  show_annotation_name = TRUE,
  annotation_name_side = "top",    # explicit headers "NES" / "BH FDR" above the stats
  annotation_name_rot = 0,
  annotation_name_gp = gpar(fontfamily = FONT, fontsize = 6.2, fontface = "bold"),
  annotation_name_offset = unit(3, "mm"),
  annotation_width = unit(c(9, 11), "mm"))

ca <- columnAnnotation(
  "recurrence" = anno_text(rec_vals, gp = gpar(fontfamily = FONT, fontsize = 6.5, col = "grey15"),
                           location = 0.5, just = "center"),
  show_annotation_name = FALSE,
  annotation_height = unit(5, "mm"))

col_fun <- colorRamp2(c(-0.5, 0, 1.5), c(OTHER_C, NEU_C, HIGH_C))

# -----------------------------------------------------------------------------
# 6. Circle-matrix heatmap
# -----------------------------------------------------------------------------
ht <- Heatmap(mat,
  col = col_fun, name = "log2FC\n(High vs Other)",
  na_col = "white",
  cell_fun = function(j, i, x, y, w, h, fill) {
    v <- mat[i, j]
    if (!is.na(v)) grid.circle(x, y, r = unit(1.45, "mm"),
                               gp = gpar(fill = col_fun(v), col = "grey55", lwd = 0.5))
  },
  row_split = factor(c(rep("High-enriched", 7), rep("Other-enriched", 3)),
                     levels = c("High-enriched", "Other-enriched")),
  row_gap = unit(4, "mm"),
  column_split = factor(gene_block, levels = c("High", "Other")),  # subtle separator High/Other gene blocks
  column_gap = unit(3, "mm"),
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_dend = FALSE, show_column_dend = FALSE,
  row_title = c("High-enriched", "Other-enriched"),
  row_title_side = "left", row_title_rot = 0,
  row_title_gp = gpar(fontfamily = FONT, fontsize = 6.5, fontface = "bold", col = "grey25"),
  row_names_side = "left",
  row_names_gp = gpar(fontfamily = FONT, fontsize = 6.8, fontface = "bold"),
  row_names_max_width = unit(34, "mm"),
  right_annotation = ra,
  column_names_side = "top", column_names_rot = 90,
  column_names_gp = gpar(fontfamily = FONT, fontsize = 6.5),
  column_names_max_height = unit(11, "mm"),
  top_annotation = ca,
  width = unit(96, "mm"), height = unit(72, "mm"),
  border = FALSE,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title = "log2FC\n(High vs Other)",
    title_gp = gpar(fontfamily = FONT, fontsize = 6.8, fontface = "bold"),
    labels_gp = gpar(fontfamily = FONT, fontsize = 6.3),
    at = c(-0.5, 0, 0.5, 1.0, 1.5),
    grid_width = unit(3.5, "mm"), legend_height = unit(34, "mm"),
    border = NA))

# -----------------------------------------------------------------------------
# 7. Assembly (panel letter + title + panel + caption note)
# -----------------------------------------------------------------------------
LAYOUT <- list(panel = c(y = 14, h = 100), ti = c(y = 111, h = 9))

upr_note <- "Genes are strict_v2 score-independent leading-edge members; every YAP17/Stem21 score gene is excluded. Circle colour = pseudobulk log2FC (High vs Other)."

draw_f4c <- function() {
  grid.newpage()
  pushViewport(viewport(x = unit(0, "mm"), y = unit(LAYOUT$ti["y"], "mm"),
                        width = unit(CANVAS_W_MM, "mm"), height = unit(LAYOUT$ti["h"], "mm"),
                        just = c("left", "bottom")))
  grid.text("C", x = unit(3, "mm"), y = unit(0.5, "npc"), just = "left",
            gp = gpar(fontfamily = FONT, fontsize = 11.5, fontface = "bold"))
  grid.text("Leading-edge gene–pathway architecture of the YAP–Stem-high program",
            x = unit(11, "mm"), y = unit(0.5, "npc"), just = "left",
            gp = gpar(fontfamily = FONT, fontsize = 9.5, fontface = "bold"))
  popViewport()
  pushViewport(viewport(x = unit(0, "mm"), y = unit(LAYOUT$panel["y"], "mm"),
                        width = unit(CANVAS_W_MM, "mm"), height = unit(LAYOUT$panel["h"], "mm"),
                        just = c("left", "bottom")))
  draw(ht, newpage = FALSE, heatmap_legend_side = "right",
       padding = unit(c(2, 2, 2, 2), "mm"))
  # Explicit headers above the right-side statistics (NES / BH FDR).
  # Fixed deterministic geometry: body left ~46.2 mm, body width 96 mm;
  # right annotation = NES column (9 mm) then BH FDR column (11 mm); body top at 95 mm.
  body_left_mm <- 46.2
  nes_cx <- body_left_mm + 96 + 4.5
  fdr_cx <- body_left_mm + 96 + 9 + 5.5
  hdr_y  <- 97
  grid.text("NES", x = unit(nes_cx, "mm"), y = unit(hdr_y, "mm"), just = "center",
            gp = gpar(fontfamily = FONT, fontsize = 6.2, fontface = "bold", col = "grey15"))
  grid.text("BH FDR", x = unit(fdr_cx, "mm"), y = unit(hdr_y, "mm"), just = "center",
            gp = gpar(fontfamily = FONT, fontsize = 6.2, fontface = "bold", col = "grey40"))
  popViewport()
}

# -----------------------------------------------------------------------------
# 8. Render PDF + SVG + PNG
# -----------------------------------------------------------------------------
cairo_pdf(PDF_OUT, width = CANVAS_W_MM / 25.4, height = CANVAS_H_MM / 25.4, family = FONT)
draw_f4c()
dev.off()
cat("[PDF] wrote", PDF_OUT, "\n")

svglite::svglite(SVG_OUT, width = CANVAS_W_MM / 25.4, height = CANVAS_H_MM / 25.4)
draw_f4c()
dev.off()
cat("[SVG] wrote", SVG_OUT, "\n")

CairoPNG(PNG_OUT, width = round(CANVAS_W_MM / 25.4 * 600), height = round(CANVAS_H_MM / 25.4 * 600),
         res = 600, bg = "white")
draw_f4c()
dev.off()
cat("[PNG] wrote", PNG_OUT, " (600 dpi)\n")

# -----------------------------------------------------------------------------
# 9. Audit table (DISPLAYED genes x 10 pathways, membership TRUE/FALSE)
#    + full unfiltered pathway–leading-edge table (supplementary source data)
# -----------------------------------------------------------------------------
audit <- do.call(rbind, lapply(seq_len(nrow(gsea)), function(i) {
  if (!gsea$pathway[i] %in% sel_pw) return(NULL)
  gs <- strsplit(gsea$leading_edge_genes[i], "/")[[1]]
  data.frame(pathway = unname(PW_DISP[gsea$pathway[i]]),
             NES = gsea$NES[i], BH_FDR = gsea$BH_FDR[i],
             gene = sel_g,
             core_enrichment_membership = sel_g %in% gs,
             strict_v2_log2FC = de$overall_log2FC_High_vs_Other[match(sel_g, de$gene)],
             strict_v2_FDR = de$BH_FDR[match(sel_g, de$gene)],
             recurrence_count = sapply(sel_g, function(ge) length(pw_of_gene[[ge]])),
             patient_positive_n = de$positive_patient_count[match(sel_g, de$gene)],
             patient_total_n = de$positive_patient_count[match(sel_g, de$gene)] + de$negative_patient_count[match(sel_g, de$gene)],
             stringsAsFactors = FALSE)
}))
write.table(audit, AUDIT_OUT, sep = "\t", row.names = FALSE, quote = FALSE)
cat("[TSV] audit table wrote", basename(AUDIT_OUT), "(", nrow(audit), "rows = 22 displayed genes x 10 pathways )\n")

full <- do.call(rbind, lapply(seq_len(nrow(gsea)), function(i) {
  if (!gsea$pathway[i] %in% sel_pw) return(NULL)
  gs <- strsplit(gsea$leading_edge_genes[i], "/")[[1]]
  data.frame(pathway = unname(PW_DISP[gsea$pathway[i]]), gene = gs,
             stringsAsFactors = FALSE)
}))
write.table(full, FULL_OUT, sep = "\t", row.names = FALSE, quote = FALSE)
cat("[TSV] full pathway–leading-edge table wrote", basename(FULL_OUT), "(", nrow(full), "rows )\n")

stopifnot(file.exists(PDF_OUT), file.exists(SVG_OUT), file.exists(PNG_OUT))
cat("[DONE] Figure 4C rendered.\n")
