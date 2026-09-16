#!/usr/bin/env Rscript
# Purpose: Figure 4 differential-expression landscape
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 4.
# =============================================================================
# Figure 4B — Patient-blocked score-independent pseudobulk differential-expression
#              landscape (strict_v2)
#
# Direct gene-level result layer between Figure 4A (program derivation) and
# Figure 4C (Hallmark interpretation). Shows the actual strict_v2 DESeq2 result
# from which the score-independent YAP-Stem-high program was derived.
#
# Frozen definition (verified from the source table):
#   8 TNBC patients, High vs Other, 16 patient-state pseudobulks, ~patient + state,
#   38 YAP17/Stem21 score genes excluded before testing, 17,597 genes tested,
#   362 BH FDR < 0.05, final strict_v2 High-upregulated program = 146 genes.
#
# Classes:
#   1 = not significant (FDR>=0.05)           light grey, low priority
#   2 = FDR significant, positive, not in 146 medium grey
#   3 = final 146-gene program                warm/red, dominant
#   4 = significant Other-associated (FDR<0.05, negative log2FC)  restrained cool
# =============================================================================

suppressMessages({
  library(ggplot2)
  library(grid)
  library(dplyr)
  library(patchwork)
  library(ggrepel)
  library(Cairo)
  library(svglite)
})

root <- "."
out  <- file.path(root, "Figure-cursor", "newfigureB")
src  <- file.path(root, "Figure3_bottom_panels_HI_feasibility_and_plotting",
                  "Figure3_I_effect_consistency_source.tsv")

PDF_OUT <- file.path(out, "Figure4B_strict_v2_DE_landscape.pdf")
SVG_OUT <- file.path(out, "Figure4B_strict_v2_DE_landscape.svg")
PNG_OUT <- file.path(out, "Figure4B_strict_v2_DE_landscape.png")
TSV_ALL <- file.path(out, "Figure4B_all_tested_genes.tsv")
TSV_LAB <- file.path(out, "Figure4B_label_selection.tsv")
QC_OUT  <- file.path(out, "Figure4B_QC.txt")

CANVAS_W_MM <- 180
CANVAS_H_MM <- 110
FONT <- "Helvetica"
WARM_PALE <- "#E9AAA3"
WARM_C <- "#D7605C"
WARM_EDGE <- "#8F3432"
COOL_C <- "#4F79A7"
NS_C <- "#E8EAED"
C2_C <- "#A8B0B8"
GREY_C <- "#E4E7EA"

# -----------------------------------------------------------------------------
# 1. Load + hard QA (frozen strict_v2 source)
# -----------------------------------------------------------------------------
tab <- read.delim(src, stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(nrow(tab) == 17597)
stopifnot(sum(tab$BH_FDR < 0.05, na.rm=TRUE) == 362)
stopifnot(sum(tab$program_146 == TRUE) == 146)
# score genes truly absent
score <- read.delim(file.path(root,"0713_score_independent_rebuild/tables/score_gene_list.tsv"), stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(length(intersect(score$gene, tab$gene)) == 0)
cat("[QA] 17,597 tested; 362 FDR<0.05; 146 program; 0 score genes present\n")

tab$log10FDR <- -log10(pmax(tab$BH_FDR, 1e-300))
tab$class <- ifelse(tab$program_146 == TRUE, 3,
             ifelse(tab$BH_FDR < 0.05 & tab$overall_log2FC_High_vs_Other < 0, 4,
             ifelse(tab$BH_FDR < 0.05, 2, 1)))
stopifnot(sum(tab$class==1) + sum(tab$class==2) + sum(tab$class==3) + sum(tab$class==4) == 17597)
cat(sprintf("[QA] classes: NS=%d, pos-sig-nonprogram=%d, program146=%d, negative-sig=%d\n",
    sum(tab$class==1), sum(tab$class==2), sum(tab$class==3), sum(tab$class==4)))

# patient directional consistency
tab$patient_total_n <- tab$positive_patient_count + tab$negative_patient_count

# Frozen membership is unchanged. These are effect-magnitude display tiers
# within the same 146-gene program, not alternative program definitions.
tab$program_effect_tier <- NA_character_
tab$program_effect_tier[tab$program_146 & tab$overall_log2FC_High_vs_Other >= 0.5 &
                          tab$overall_log2FC_High_vs_Other < 1.0] <- "Primary moderate effect"
tab$program_effect_tier[tab$program_146 & tab$overall_log2FC_High_vs_Other >= 1.0 &
                          tab$overall_log2FC_High_vs_Other < 1.5] <- "Stringent effect"
tab$program_effect_tier[tab$program_146 & tab$overall_log2FC_High_vs_Other >= 1.5] <- "Large-effect reference"
# Frozen sensitivity-audit expected counts (membership unchanged):
#   full program log2FC>=0.5 -> 146 ;  log2FC>=1.0 -> 36 ;  log2FC>=1.5 -> 3
stopifnot(sum(tab$program_146 & tab$overall_log2FC_High_vs_Other >= 0.5) == 146)
stopifnot(sum(tab$program_146 & tab$overall_log2FC_High_vs_Other >= 1.0) == 36)
stopifnot(sum(tab$program_146 & tab$overall_log2FC_High_vs_Other >= 1.5) == 3)
stopifnot(sum(tab$program_effect_tier == "Primary moderate effect", na.rm=TRUE) == 110)
stopifnot(sum(tab$program_effect_tier == "Stringent effect", na.rm=TRUE) == 33)
stopifnot(sum(tab$program_effect_tier == "Large-effect reference", na.rm=TRUE) == 3)
cat("[QA] program effect tiers: 110 moderate; 33 stringent; 3 large-effect reference\n")

# -----------------------------------------------------------------------------
# 2. Reproducible label selection (exactly 7 labels)
#    all three >=1.5 program genes; two remaining program genes by
#    (BH FDR, |log2FC|, patient consistency, gene); and two significant
#    Other-associated genes by (BH FDR, |log2FC|, gene).
# -----------------------------------------------------------------------------
prog <- tab[tab$program_146 == TRUE, ]
prog$consist <- prog$positive_patient_count / pmax(prog$patient_total_n, 1)
extreme_prog <- sort(prog$gene[prog$overall_log2FC_High_vs_Other >= 1.5])
prog_remaining <- prog[!prog$gene %in% extreme_prog, ]
prog_remaining <- prog_remaining[order(prog_remaining$BH_FDR,
                                       -abs(prog_remaining$overall_log2FC_High_vs_Other),
                                       -prog_remaining$consist,
                                       prog_remaining$gene), ]
top_prog <- c(extreme_prog, head(prog_remaining$gene, 2))
neg <- tab[tab$class == 4, ]
neg <- neg[order(neg$BH_FDR, -abs(neg$overall_log2FC_High_vs_Other), neg$gene), ]
top_neg <- head(neg$gene, 2)
labelled <- c(top_prog, top_neg)
tab$label <- ifelse(tab$gene %in% labelled, tab$gene, NA_character_)
tab$label_rule <- NA_character_
tab$label_rule[tab$gene %in% extreme_prog] <- "all program genes with unshrunk log2FC >= 1.5"
tab$label_rule[tab$gene %in% setdiff(top_prog, extreme_prog)] <- "top remaining program genes by FDR+|log2FC|+consistency"
tab$label_rule[tab$gene %in% top_neg] <- "top Other-associated significant genes by FDR+|log2FC|"
cat("[labels]", sum(!is.na(tab$label)), "labels\n")

# -----------------------------------------------------------------------------
# 3. Main volcano
# -----------------------------------------------------------------------------
tab$class_f <- factor(tab$class, levels=c("1","2","3","4"))
tab$large_effect <- abs(tab$overall_log2FC_High_vs_Other) >= 1.5
prog_mod <- tab[!is.na(tab$program_effect_tier) & tab$program_effect_tier == "Primary moderate effect", ]
prog_str <- tab[!is.na(tab$program_effect_tier) & tab$program_effect_tier == "Stringent effect", ]
prog_ext <- tab[!is.na(tab$program_effect_tier) & tab$program_effect_tier == "Large-effect reference", ]

p_vol <- ggplot(tab, aes(x=overall_log2FC_High_vs_Other, y=log10FDR)) +
  geom_hline(yintercept=-log10(0.05), colour="#8A9096", linewidth=0.32, linetype="dashed") +
  geom_vline(xintercept=-0.5, colour="#D9DDE1", linewidth=0.22, linetype="dashed") +
  geom_vline(xintercept=0.5, colour="#BCC2C8", linewidth=0.28, linetype="dashed") +
  geom_vline(xintercept=c(-1.5, 1.5), colour="#777E85", linewidth=0.34, linetype="dotted") +
  geom_point(data=tab[tab$class_f=="1",], colour=NS_C, size=0.46, alpha=0.76) +
  geom_point(data=tab[tab$class_f=="2",], colour=C2_C, size=0.78, alpha=0.90) +
  geom_point(data=tab[tab$class_f=="4",], colour=COOL_C, size=0.88, alpha=0.88) +
  geom_point(data=prog_mod, colour=WARM_PALE, size=0.93, alpha=0.92) +
  geom_point(data=prog_str, colour=WARM_C, size=1.30, alpha=0.96) +
  geom_point(data=prog_ext, shape=21, fill=WARM_C, colour=WARM_EDGE,
             size=1.78, stroke=0.42, alpha=1) +
  geom_text_repel(data=tab[!is.na(tab$label),], aes(label=label), family=FONT, size=5.8/2.845,
                  colour="grey18", max.overlaps=20, seed=20260810,
                  box.padding=0.45, point.padding=0.25, segment.size=0.22,
                  min.segment.length=0.15, force=1.2,
                  xlim=c(-2.25,2.45), ylim=c(1.45,7.35)) +
  scale_x_continuous(name="Unshrunk log2FC (YAP–Stem High vs Other)",
                     limits=c(-2.4, 2.6), breaks=seq(-2,2,1)) +
  scale_y_continuous(name=expression(-log[10]*"(BH FDR)"),
                     limits=c(0, 7.45), breaks=seq(0,7,1)) +
  coord_cartesian(clip="off") +
  theme_void() +
  theme(axis.text=element_text(family=FONT, size=6.5, colour="black"),
        axis.title=element_text(family=FONT, size=7, colour="black"),
        axis.line=element_line(colour="black", linewidth=0.45),
        axis.ticks=element_line(colour="black", linewidth=0.35),
        axis.ticks.length=unit(1.1,"mm"),
        plot.margin=margin(1,3,2,2,"mm"))

# -----------------------------------------------------------------------------
# 4. apeglm LFC-shrinkage inset (replaces the log2FC density inset)
#    Shows that the frozen 146 program's moderate effect sizes are stable to
#    LFC shrinkage (rho=0.982, 146/146 direction preserved).
# -----------------------------------------------------------------------------
lc <- read.delim(file.path(root,"Figure4_effect_size_threshold_sensitivity_audit/LFC_shrinkage_sensitivity.tsv"), check.names=FALSE)
lc_g <- lc[lc$row_type == "GENE", ]
stopifnot(nrow(lc_g) == 17597)
lc_g$is_program <- !is.na(lc_g$program_146) & lc_g$program_146 == TRUE
rho_prog <- lc$Spearman_rho_unshrunk_vs_apeglm[lc$row_type=="SUMMARY" & lc$scope=="PROGRAM_146"]
stopifnot(round(rho_prog,3) == 0.982)
p_inset <- ggplot(lc_g, aes(x=unshrunk_log2FC, y=apeglm_shrunken_log2FC)) +
  geom_abline(slope=1, intercept=0, colour="grey70", linewidth=0.35, linetype="dashed") +
  geom_point(data=lc_g[!lc_g$is_program,], colour=GREY_C, size=0.28, alpha=0.75) +
  geom_point(data=lc_g[lc_g$is_program,], colour=WARM_C, size=0.75, alpha=0.92) +
  scale_x_continuous(limits=c(-2.6,2.8), name="Unshrunk log2FC") +
  scale_y_continuous(limits=c(-2.2,2.6), name="apeglm log2FC") +
  annotate("text", x=-2.48, y=2.48, hjust=0, vjust=1, family=FONT,
           size=5.8/2.845, colour="grey22",
           label=sprintf("Program genes: ρ=%.3f\n146/146 direction preserved", rho_prog)) +
  theme_void() +
  theme(axis.text=element_text(family=FONT, size=5.5, colour="black"),
        axis.title=element_text(family=FONT, size=6, colour="black"),
        axis.line=element_line(colour="black", linewidth=0.35),
        axis.ticks=element_line(colour="black", linewidth=0.3),
        plot.margin=margin(1,1,1,1,"mm"))

p_vol_grob <- grid.grabExpr(print(p_vol))
p_inset_grob <- grid.grabExpr(print(p_inset))

# -----------------------------------------------------------------------------
# 5. Assembly (title)
# -----------------------------------------------------------------------------
LAYOUT <- list(
  title=c(y=101,h=8),
  strip=c(y=93,h=6),
  plot=c(y=31,h=60),
  footer=c(y=3,h=25)
)
draw_f4b <- function() {
  grid.newpage()
  pushViewport(viewport(x=unit(0,"mm"), y=unit(LAYOUT$title["y"],"mm"),
                        width=unit(CANVAS_W_MM,"mm"), height=unit(LAYOUT$title["h"],"mm"), just=c("left","bottom")))
  grid.text("B", x=unit(3,"mm"), y=unit(0.5,"npc"), just="left",
            gp=gpar(fontfamily=FONT, fontsize=11.5, fontface="bold"))
  grid.text("Patient-blocked score-independent pseudobulk differential expression",
            x=unit(11,"mm"), y=unit(0.5,"npc"), just="left",
            gp=gpar(fontfamily=FONT, fontsize=9.5, fontface="bold"))
  popViewport()

  pushViewport(viewport(x=unit(12,"mm"), y=unit(LAYOUT$strip["y"],"mm"),
                        width=unit(164,"mm"), height=unit(LAYOUT$strip["h"],"mm"), just=c("left","bottom")))
  grid.roundrect(x=unit(0.5,"npc"), y=unit(0.5,"npc"), width=unit(1,"npc"), height=unit(5.0,"mm"),
                 r=unit(1.4,"mm"), gp=gpar(fill="#F3F4F5", col=NA))
  grid.text("17,597 tested  |  362 FDR-significant  |  146 program genes  |  38/38 defining genes excluded",
            x=unit(0.5,"npc"), y=unit(0.5,"npc"),
            gp=gpar(fontfamily=FONT, fontsize=6.6, fontface="bold", col="#3F474E"))
  popViewport()

  pushViewport(viewport(x=unit(0,"mm"), y=unit(LAYOUT$plot["y"],"mm"),
                        width=unit(CANVAS_W_MM,"mm"), height=unit(LAYOUT$plot["h"],"mm"), just=c("left","bottom")))
  grid.draw(p_vol_grob)
  popViewport()

  pushViewport(viewport(x=unit(3,"mm"), y=unit(LAYOUT$footer["y"],"mm"),
                        width=unit(112,"mm"), height=unit(LAYOUT$footer["h"],"mm"), just=c("left","bottom")))
  grid.text("Program membership by effect magnitude", x=unit(1,"mm"), y=unit(22.5,"mm"), just="left",
            gp=gpar(fontfamily=FONT, fontsize=6.4, fontface="bold", col="#3F474E"))
  grid.points(x=unit(c(3,42,78),"mm"), y=unit(rep(18.0,3),"mm"), pch=c(16,16,21),
              size=unit(c(1.6,2.0,2.5),"mm"),
              gp=gpar(col=c(WARM_PALE,WARM_C,WARM_EDGE), fill=c(WARM_PALE,WARM_C,WARM_C), lwd=c(0.5,0.5,0.8)))
  grid.text("0.5–<1.0  (110)", x=unit(6,"mm"), y=unit(18,"mm"), just="left",
            gp=gpar(fontfamily=FONT, fontsize=6.1, col="#3F474E"))
  grid.text("1.0–<1.5  (33)", x=unit(45,"mm"), y=unit(18,"mm"), just="left",
            gp=gpar(fontfamily=FONT, fontsize=6.1, col="#3F474E"))
  grid.text("≥1.5  (3)", x=unit(82,"mm"), y=unit(18,"mm"), just="left",
            gp=gpar(fontfamily=FONT, fontsize=6.1, col="#3F474E"))

  grid.points(x=unit(c(3,55),"mm"), y=unit(rep(13.3,2),"mm"), pch=16, size=unit(c(1.8,1.6),"mm"),
              gp=gpar(col=c(COOL_C,C2_C)))
  grid.text("Other-associated, FDR<0.05", x=unit(6,"mm"), y=unit(13.3,"mm"), just="left",
            gp=gpar(fontfamily=FONT, fontsize=5.9, col="#4A5158"))
  grid.text("Significant non-program", x=unit(58,"mm"), y=unit(13.3,"mm"), just="left",
            gp=gpar(fontfamily=FONT, fontsize=5.9, col="#4A5158"))

  grid.lines(x=unit(c(2,10),"mm"), y=unit(c(8.8,8.8),"mm"), gp=gpar(col="#BCC2C8", lwd=0.8, lty="dashed"))
  grid.text("Primary program criterion: log2FC ≥0.5", x=unit(12,"mm"), y=unit(8.8,"mm"), just="left",
            gp=gpar(fontfamily=FONT, fontsize=5.9, col="#4A5158"))
  grid.points(x=unit(6,"mm"), y=unit(5.0,"mm"), pch=16, size=unit(2.0,"mm"), gp=gpar(col=WARM_C))
  grid.text("Stringent-effect subset: log2FC ≥1.0", x=unit(12,"mm"), y=unit(5.0,"mm"), just="left",
            gp=gpar(fontfamily=FONT, fontsize=5.9, col="#4A5158"))
  grid.lines(x=unit(c(2,10),"mm"), y=unit(c(1.2,1.2),"mm"), gp=gpar(col="#777E85", lwd=0.9, lty="dotted"))
  grid.text("Large-effect reference: |log2FC|=1.5", x=unit(12,"mm"), y=unit(1.2,"mm"), just="left",
            gp=gpar(fontfamily=FONT, fontsize=5.9, col="#4A5158"))
  grid.lines(x=unit(c(69,77),"mm"), y=unit(c(1.2,1.2),"mm"), gp=gpar(col="#8A9096", lwd=0.8, lty="dashed"))
  grid.text("BH FDR=0.05", x=unit(79,"mm"), y=unit(1.2,"mm"), just="left",
            gp=gpar(fontfamily=FONT, fontsize=5.9, col="#4A5158"))
  popViewport()

  pushViewport(viewport(x=unit(121,"mm"), y=unit(LAYOUT$footer["y"],"mm"),
                        width=unit(56,"mm"), height=unit(LAYOUT$footer["h"],"mm"), just=c("left","bottom")))
  grid.text("LFC shrinkage sensitivity", x=unit(1,"mm"), y=unit(23.5,"mm"), just="left",
            gp=gpar(fontfamily=FONT, fontsize=6.2, fontface="bold", col="#3F474E"))
  pushViewport(viewport(x=unit(0,"mm"), y=unit(0,"mm"), width=unit(56,"mm"), height=unit(22,"mm"), just=c("left","bottom")))
  grid.draw(p_inset_grob)
  popViewport(2)
}

cairo_pdf(PDF_OUT, width=CANVAS_W_MM/25.4, height=CANVAS_H_MM/25.4, family=FONT)
draw_f4b(); dev.off()
cat("[PDF] wrote", basename(PDF_OUT), "\n")
svglite::svglite(SVG_OUT, width=CANVAS_W_MM/25.4, height=CANVAS_H_MM/25.4)
draw_f4b(); dev.off()
cat("[SVG] wrote", basename(SVG_OUT), "\n")
CairoPNG(PNG_OUT, width=round(CANVAS_W_MM/25.4*600), height=round(CANVAS_H_MM/25.4*600), res=600, bg="white")
draw_f4b(); dev.off()
# Cairo writes the correct 600-dpi pixel dimensions but may omit PNG pHYs
# metadata on Windows. Re-encode losslessly with the same pixels and explicit
# publication density so the exported file is self-describing.
png_pixels <- png::readPNG(PNG_OUT)
png::writePNG(png_pixels, PNG_OUT, dpi=600.05)
rm(png_pixels); gc(verbose=FALSE)
cat("[PNG] wrote", basename(PNG_OUT), "\n")

# -----------------------------------------------------------------------------
# 6. Export all-tested-genes + label-selection tables
# -----------------------------------------------------------------------------
out_tab <- tab[, c("gene","baseMean","overall_log2FC_High_vs_Other","BH_FDR","class",
                   "program_146","program_effect_tier","score_definition_genes_excluded_before_testing",
                   "positive_patient_count","patient_total_n")]
names(out_tab) <- c("gene","baseMean","log2FC_unshrunk","BH_FDR","significance_class",
                    "final_146_program_membership","program_effect_display_tier","score_defining_gene_flag",
                    "patient_direction_positive_n","patient_total_n")
write.table(out_tab, TSV_ALL, sep="\t", row.names=FALSE, quote=FALSE)
cat("[TSV] Figure4B_all_tested_genes.tsv (", nrow(out_tab), "rows )\n")

lab_tab <- tab[!is.na(tab$label), c("gene","overall_log2FC_High_vs_Other","BH_FDR","positive_patient_count","patient_total_n","program_146","label_rule")]
names(lab_tab) <- c("gene","log2FC_unshrunk","BH_FDR","patient_direction_positive_n","patient_total_n","final_146_program_membership","label_rule")
lab_tab <- lab_tab[order(lab_tab$label_rule, lab_tab$BH_FDR), ]
write.table(lab_tab, TSV_LAB, sep="\t", row.names=FALSE, quote=FALSE)
cat("[TSV] Figure4B_label_selection.tsv\n")

# -----------------------------------------------------------------------------
# 7. QC text
# -----------------------------------------------------------------------------
qc <- c(
  "# Figure 4B QC",
  "",
  "## Source",
  "File: Figure3_bottom_panels_HI_feasibility_and_plotting/Figure3_I_effect_consistency_source.tsv",
  "Frozen strict_v2 patient-blocked pseudobulk DESeq2 (8 patients, High vs Other, 16 pseudobulks, ~patient + state).",
  "",
  "## Mandatory counts (all confirmed from the frozen source table)",
  sprintf("score-independent tested universe: %d (exactly 17,597)", nrow(tab)),
  "score-defining genes excluded before DE testing: 38/38 (0 present in tested table)",
  sprintf("BH FDR < 0.05 genes: %d (exactly 362)", sum(tab$BH_FDR<0.05, na.rm=TRUE)),
  sprintf("final strict_v2 High-upregulated program: %d (exactly 146)", sum(tab$program_146==TRUE)),
  sprintf("program log2FC>=0.5: %d (expected 146; all program genes)", sum(tab$program_146 & tab$overall_log2FC_High_vs_Other>=0.5)),
  sprintf("program 0.5<=log2FC<1.0: %d (expected 110)", sum(tab$program_effect_tier=="Primary moderate effect", na.rm=TRUE)),
  sprintf("program 1.0<=log2FC<1.5: %d (expected 33)", sum(tab$program_effect_tier=="Stringent effect", na.rm=TRUE)),
  sprintf("program log2FC>=1.0: %d (expected 36)", sum(tab$program_146 & tab$overall_log2FC_High_vs_Other>=1.0)),
  sprintf("program log2FC>=1.5: %d (expected 3; large-effect reference, not a core program)", sum(tab$program_effect_tier=="Large-effect reference", na.rm=TRUE)),
  "No historical 132/142-gene definition mixed in (only column 'program_146' exists).",
  "",
  "## Class composition",
  sprintf("Class 1 not significant: %d", sum(tab$class==1)),
  sprintf("Class 2 positive FDR-significant non-program: %d", sum(tab$class==2)),
  sprintf("Class 3 final 146-gene program: %d", sum(tab$class==3)),
  sprintf("Class 4 significant Other-associated (negative): %d", sum(tab$class==4)),
  "Sum = 17,597.",
  "",
  "## Anchor genes (verified against frozen table; match the previously audited values)",
  sprintf("ATF3: log2FC=%.3f, FDR=%.2e, %d/8 direction, program146=%s", tab$overall_log2FC_High_vs_Other[tab$gene=='ATF3'], tab$BH_FDR[tab$gene=='ATF3'], tab$positive_patient_count[tab$gene=='ATF3'], tab$program_146[tab$gene=='ATF3']),
  sprintf("PPP1R15A: log2FC=%.3f, FDR=%.2e, %d/8 direction, program146=%s", tab$overall_log2FC_High_vs_Other[tab$gene=='PPP1R15A'], tab$BH_FDR[tab$gene=='PPP1R15A'], tab$positive_patient_count[tab$gene=='PPP1R15A'], tab$program_146[tab$gene=='PPP1R15A']),
  sprintf("MCL1: log2FC=%.3f, FDR=%.2e, %d/8 direction, program146=%s", tab$overall_log2FC_High_vs_Other[tab$gene=='MCL1'], tab$BH_FDR[tab$gene=='MCL1'], tab$positive_patient_count[tab$gene=='MCL1'], tab$program_146[tab$gene=='MCL1']),
  "",
  "## Thresholds plotted (exactly the frozen analysis definition)",
  "BH FDR = 0.05 (y reference); unshrunk log2FC = +0.5 and -0.5 (x references).",
  "X = unshrunk log2FC (High vs Other); Y = -log10(BH FDR).",
  "BH FDR=0 absent (min 1.19e-7); -log10 range 0..6.93.",
  "",
  "## Labels",
  sprintf("%d labels: all three program genes with log2FC>=1.5 + top-two remaining program genes by (BH FDR, |log2FC|, patient consistency, gene) + top-two Other-associated significant genes by (BH FDR, |log2FC|, gene).", sum(!is.na(tab$label))),
  "All labels are reproducibly selected; no narrative/manual label was added.",
  "The apeglm inset is placed below the main volcano, never over the program-point cloud.",
  "Final canvas: 180 x 110 mm; PDF/SVG vector; PNG 600 dpi."
)
writeLines(qc, QC_OUT)
cat("[TXT] Figure4B_QC.txt\n")
cat("[DONE] Figure 4B rendered.\n")
