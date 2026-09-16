# Purpose: Figure 2 volcano threshold panel
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 1–2.
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(data.table)
})

root <- "."
out <- file.path(root, "Figure2_pseudobulk_DE_threshold_patch")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

source_file <- file.path(
  root,
  "0723_Figure1_and_ED1_targeted_scientific_visual_closure_v5",
  "plotting_data",
  "ExtendedData1_panel_b_full_malignant_normal_volcano.tsv"
)
palette_file <- file.path(root, "project_semantic_palette.tsv")

d <- fread(source_file, data.table = FALSE)
pal <- fread(palette_file, data.table = FALSE)
col <- setNames(pal$hex, pal$token)

fc_cutoff <- 1.5
fdr_cutoff <- 0.05

d$display_status <- "Other genes"
d$display_status[
  d$gene_class == "Normal-high epithelial/luminal" &
    d$logFC <= -fc_cutoff & d$adj.P.Val < fdr_cutoff
] <- "Normal-high epithelial/luminal"
d$display_status[
  d$gene_class == "Malignant-high proliferation" &
    d$logFC >= fc_cutoff & d$adj.P.Val < fdr_cutoff
] <- "Malignant-high proliferation"
d$display_status <- factor(
  d$display_status,
  levels = c(
    "Other genes",
    "Normal-high epithelial/luminal",
    "Malignant-high proliferation"
  )
)
d$minus_log10_FDR_plot <- -log10(pmax(d$adj.P.Val, .Machine$double.xmin))

representative <- c("PIGR", "AGR3", "TFF1", "PGR", "BIRC5", "UBE2C", "CDK1", "TOP2A")
lab <- d[d$gene %in% representative & d$display_status != "Other genes", , drop = FALSE]
stopifnot(nrow(lab) == length(representative))

audit <- d[d$gene %in% representative, c("gene", "logFC", "adj.P.Val", "display_status")]
names(audit) <- c("gene", "log2FC", "FDR", "highlight_status")
audit$passes_display_rule <- with(
  audit,
  (highlight_status == "Normal-high epithelial/luminal" & log2FC <= -fc_cutoff & FDR < fdr_cutoff) |
    (highlight_status == "Malignant-high proliferation" & log2FC >= fc_cutoff & FDR < fdr_cutoff)
)
audit <- audit[match(representative, audit$gene), ]
write.table(
  audit,
  file.path(out, "Figure2_panelH_labeled_gene_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

blue <- unname(col["Normal"])
coral <- unname(col["Tumour"])
light_gray <- "#D9DDE1"
dark <- "#37434B"

p <- ggplot(d, aes(logFC, minus_log10_FDR_plot)) +
  geom_point(
    data = d[d$display_status == "Other genes", , drop = FALSE],
    colour = light_gray, alpha = 0.42, size = 0.48, stroke = 0
  ) +
  geom_point(
    data = d[d$display_status != "Other genes", , drop = FALSE],
    aes(colour = display_status), alpha = 0.94, size = 1.25, stroke = 0
  ) +
  geom_vline(
    xintercept = c(-fc_cutoff, fc_cutoff),
    linetype = "22", linewidth = 0.36, colour = "#7F878D"
  ) +
  geom_hline(
    yintercept = -log10(fdr_cutoff),
    linetype = "22", linewidth = 0.36, colour = "#7F878D"
  ) +
  geom_text_repel(
    data = lab,
    aes(label = gene, colour = display_status),
    family = "Arial", fontface = "italic", size = 2.35,
    max.overlaps = Inf, min.segment.length = 0,
    segment.size = 0.24, segment.colour = "#7F878D",
    box.padding = 0.24, point.padding = 0.12, show.legend = FALSE,
    seed = 20260808
  ) +
  annotate(
    "text", x = -fc_cutoff, y = 0.035, label = "-1.5",
    family = "Arial", size = 2.0, vjust = 1.1, colour = "#596168"
  ) +
  annotate(
    "text", x = fc_cutoff, y = 0.035, label = "+1.5",
    family = "Arial", size = 2.0, vjust = 1.1, colour = "#596168"
  ) +
  annotate(
    "text", x = max(d$logFC, na.rm = TRUE), y = -log10(fdr_cutoff),
    label = "BH FDR = 0.05", hjust = 1.02, vjust = -0.45,
    family = "Arial", size = 2.0, colour = "#596168"
  ) +
  scale_colour_manual(
    values = c(
      "Normal-high epithelial/luminal" = blue,
      "Malignant-high proliferation" = coral
    ),
    breaks = c("Normal-high epithelial/luminal", "Malignant-high proliferation"),
    labels = c("Normal-high", "Malignant-high"),
    name = "Blue and red points indicate genes passing both\n|log2FC| and FDR thresholds."
  ) +
  labs(
    title = "H  Pseudobulk DE audit",
    x = expression(log[2]~fold~change~("malignant - normal")),
    y = expression(-log[10]~("BH FDR")),
    caption = "Display cutoff: |log2FC| >= 1.5 and BH FDR < 0.05."
  ) +
  guides(colour = guide_legend(
    title.position = "top", title.hjust = 0, nrow = 1, byrow = TRUE,
    override.aes = list(size = 2.6, alpha = 1)
  )) +
  theme_classic(base_family = "Arial", base_size = 7.1) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(size = 8.4, face = "bold", colour = dark, hjust = 0),
    axis.title = element_text(size = 7.2, colour = dark),
    axis.text = element_text(size = 6.6, colour = dark),
    axis.line = element_line(linewidth = 0.42, colour = dark),
    axis.ticks = element_line(linewidth = 0.32, colour = dark),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.title = element_text(size = 6.5, face = "plain", colour = dark),
    legend.text = element_text(size = 6.5, colour = dark),
    legend.key.width = grid::unit(4.2, "mm"),
    legend.spacing.x = grid::unit(1.0, "mm"),
    plot.caption = element_text(size = 6.5, colour = "#596168", hjust = 0),
    plot.margin = margin(3.5, 4.0, 2.5, 3.5)
  )

ggsave(
  file.path(out, "Figure2_panelH_pseudobulk_DE_audit_thresholded.pdf"),
  p, width = 120, height = 70, units = "mm",
  device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE
)
ggsave(
  file.path(out, "Figure2_panelH_pseudobulk_DE_audit_thresholded.png"),
  p, width = 120, height = 70, units = "mm", dpi = 601,
  device = "png", type = "cairo", bg = "white", limitsize = FALSE
)

cat("Rows:", nrow(d), "\n")
cat("Blue:", sum(d$display_status == "Normal-high epithelial/luminal"), "\n")
cat("Red:", sum(d$display_status == "Malignant-high proliferation"), "\n")
cat("FDR-significant genes:", sum(d$adj.P.Val < fdr_cutoff, na.rm = TRUE), "\n")
cat("Minimum |log2FC| among FDR-significant genes:",
    min(abs(d$logFC[d$adj.P.Val < fdr_cutoff]), na.rm = TRUE), "\n")
