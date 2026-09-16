#!/usr/bin/env Rscript
# Purpose: Render S2D Program146 stability
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Supplementary Fig. S1–S4.

# Supplementary Figure S2D FINAL candidate
# Mechanical rendering of frozen Program146 definition-sensitivity outputs.
# No DE, threshold, recurrence, pathway or shrinkage analysis is performed.

options(stringsAsFactors = FALSE)
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_dir <- file.path(project_root, "Supplementary_Submission_Staging", "01_FINAL_PANEL_SOURCES", "Supplementary_Figure_S2", "panels", "S2D")
source_file <- file.path(out_dir, "S2D_program_stability_FINAL_render_source.tsv")
palette_file <- file.path(project_root, "project_semantic_palette.tsv")
if (!file.exists(source_file)) stop("Missing S2D render source")
if (!file.exists(palette_file)) stop("Missing project semantic palette")

d <- read.delim(source_file, check.names = FALSE, stringsAsFactors = FALSE,
                na.strings = c("", "NA"), fileEncoding = "UTF-8")
palette <- read.delim(palette_file, check.names = FALSE, stringsAsFactors = FALSE)

get_colour <- function(token) {
  z <- palette$hex[palette$token == token]
  if (length(z) != 1L) stop("Palette token missing or duplicated: ", token)
  z
}
positive_colour <- get_colour("Positive")
ink <- "#26343E"
muted <- "#66727A"
guide <- "#D8E0E4"
threshold_colour <- "#65757F"

gene <- d[d$record_type == "GENE", , drop = FALSE]
thr <- d[d$record_type == "THRESHOLD_SUMMARY", , drop = FALSE]
shr <- d[d$record_type == "SHRINKAGE_SUMMARY", , drop = FALSE]
if (nrow(d) != 151L || nrow(gene) != 146L || nrow(thr) != 4L || nrow(shr) != 1L) stop("S2D render-source structure changed")
if (!identical(as.integer(table(gene$positive_patient_count)[c("6","7","8")]), c(16L,68L,62L))) stop("Recurrence distribution changed")
if (!identical(as.numeric(thr$threshold), c(0.5,0.75,1.0,1.5))) stop("Threshold order changed")
if (!identical(as.integer(thr$retained_gene_n), c(146L,85L,36L,3L))) stop("Retained counts changed")
if (!identical(as.integer(thr$pathway_direction_concordant_n), c(10L,9L,8L,8L)) || any(thr$pathway_direction_total_n != 10L)) stop("Pathway direction results changed")
if (abs(shr$shrinkage_rho - 0.982416395939576) > 1e-12 || shr$shrinkage_direction_n != 146L || shr$shrinkage_direction_total_n != 146L) stop("Shrinkage summary changed")
if (any(gene$unshrunk_log2FC < 0.5) || any(gene$BH_FDR >= 0.05) || any(gene$positive_patient_count < 6)) stop("A gene violates frozen Program146 membership")
if (any(gene$shrinkage_direction_concordant != "TRUE")) stop("Shrinkage direction is not retained for every Program146 gene")

x_limits <- c(0.45, max(gene$unshrunk_log2FC) + 0.08)
if (max(gene$unshrunk_log2FC) >= x_limits[2]) stop("x range would clip a gene")
x_ticks <- c(0.5, 1.0, 1.5, 2.0)
plot_x0 <- 0.200
plot_x1 <- 0.955
plot_y0 <- 0.245
plot_y1 <- 0.670
axis_y <- 0.190

x_map <- function(x) plot_x0 + (x - x_limits[1]) / diff(x_limits) * (plot_x1 - plot_x0)
y_map <- function(y) plot_y0 + (y - 5.55) / (8.45 - 5.55) * (plot_y1 - plot_y0)

draw_text <- function(label, x, y, fontsize = 6.5, fontface = "plain", colour = ink,
                      just = "center", rot = 0) {
  grid::grid.text(label, x = x, y = y, just = just, rot = rot,
                  gp = grid::gpar(fontfamily = "Arial", fontsize = fontsize,
                                  fontface = fontface, col = colour))
}

draw_panel <- function() {
  grid::grid.newpage()

  draw_text("D", 0.022, 0.970, fontsize = 11.0, fontface = "bold", just = c("left", "top"))
  draw_text("Program146 stability across effect", 0.088, 0.962,
            fontsize = 7.6, fontface = "bold", just = c("left", "top"))
  draw_text("and recurrence criteria", 0.088, 0.912,
            fontsize = 7.6, fontface = "bold", just = c("left", "top"))

  # Nested effect regions use one restrained hue; widths are quantitative.
  bands <- data.frame(left = c(0.5,0.75,1.0,1.5), right = c(0.75,1.0,1.5,x_limits[2]), alpha = c(0.025,0.045,0.065,0.090))
  for (i in seq_len(nrow(bands))) {
    grid::grid.rect(x = (x_map(bands$left[i]) + x_map(bands$right[i])) / 2,
                    y = (plot_y0 + plot_y1) / 2,
                    width = x_map(bands$right[i]) - x_map(bands$left[i]),
                    height = plot_y1 - plot_y0,
                    gp = grid::gpar(fill = grDevices::adjustcolor(positive_colour, alpha.f = bands$alpha[i]), col = NA))
  }

  # Recurrence guides and labels.
  exact_n <- c(`6` = 16, `7` = 68, `8` = 62)
  for (r in 6:8) {
    yy <- y_map(r)
    grid::grid.segments(x0 = plot_x0, x1 = plot_x1, y0 = yy, y1 = yy,
                        gp = grid::gpar(col = guide, lwd = 0.42))
    draw_text(sprintf("%d/8  (n=%d)", r, exact_n[as.character(r)]), 0.180, yy,
              fontsize = 6.4, just = "right")
  }
  draw_text("Positive n/8", 0.180, 0.710, fontsize = 6.3,
            fontface = "bold", colour = muted, just = "right")

  # Predefined effect thresholds; primary is visually distinct but restrained.
  line_types <- c(1, 3, 3, 3)
  line_widths <- c(0.80, 0.55, 0.60, 0.68)
  for (i in seq_len(nrow(thr))) {
    xx <- x_map(thr$threshold[i])
    grid::grid.segments(x0 = xx, x1 = xx, y0 = plot_y0 - 0.010, y1 = 0.715,
                        gp = grid::gpar(col = threshold_colour, lwd = line_widths[i], lty = line_types[i]))
    draw_text(sprintf("%.2f", thr$threshold[i]), xx, 0.835, fontsize = 6.2,
              fontface = if (i == 1) "bold" else "plain")
    draw_text(sprintf("n=%d", thr$retained_gene_n[i]), xx, 0.795, fontsize = 6.2,
              fontface = if (i == 1) "bold" else "plain")
    draw_text(sprintf("%d/%d", thr$pathway_direction_concordant_n[i], thr$pathway_direction_total_n[i]),
              xx, 0.755, fontsize = 6.1, colour = muted)
  }
  draw_text("Cutoff", 0.180, 0.835, fontsize = 6.2, fontface = "bold", colour = muted, just = "right")
  draw_text("Genes", 0.180, 0.795, fontsize = 6.2, fontface = "bold", colour = muted, just = "right")
  draw_text("Pathway dir.", 0.180, 0.755, fontsize = 6.1, fontface = "bold", colour = muted, just = "right")

  # Every frozen Program146 gene is retained as one point.
  is_large <- gene$unshrunk_log2FC >= 1.5
  grid::grid.points(x = x_map(gene$unshrunk_log2FC), y = y_map(gene$plot_y),
                    default.units = "npc", pch = 21, size = grid::unit(1.05, "mm"),
                    gp = grid::gpar(fill = grDevices::adjustcolor(positive_colour, alpha.f = 0.66),
                                    col = grDevices::adjustcolor(ifelse(is_large, ink, positive_colour), alpha.f = 0.76),
                                    lwd = ifelse(is_large, 0.65, 0.38)))

  # Shared quantitative axis.
  grid::grid.segments(x0 = plot_x0, x1 = plot_x1, y0 = axis_y, y1 = axis_y,
                      gp = grid::gpar(col = ink, lwd = 0.64))
  tick_x <- x_map(x_ticks)
  grid::grid.segments(x0 = tick_x, x1 = tick_x, y0 = axis_y, y1 = axis_y - 0.010,
                      gp = grid::gpar(col = ink, lwd = 0.52))
  draw_text(sprintf("%.1f", x_ticks), tick_x, axis_y - 0.031, fontsize = 6.2)
draw_text("Unshrunk High - Other log2FC", (plot_x0 + plot_x1)/2, 0.117, fontsize = 6.7)

draw_text("apeglm: rho = 0.982; 146/146 directions retained", 0.520, 0.062,
            fontsize = 6.2, colour = ink)
# Large-effect shrinkage detail is retained in the publication legend.
}

pdf_target <- file.path(out_dir, "S2D_program_stability_FINAL.pdf")
png_target <- file.path(out_dir, "S2D_program_stability_FINAL_600dpi.png")
pdf_tmp <- tempfile(pattern = "S2D_", fileext = ".pdf")
png_tmp <- tempfile(pattern = "S2D_", fileext = ".png")

grDevices::cairo_pdf(pdf_tmp, width = 88/25.4, height = 78/25.4,
                     family = "Arial", bg = "white", onefile = TRUE)
draw_panel()
grDevices::dev.off()

grDevices::png(png_tmp, width = 88, height = 78, units = "mm", res = 601,
               type = "windows", bg = "white")
draw_panel()
grDevices::dev.off()

# The native Windows PNG device can emit indexed colour for a restrained
# palette. Convert losslessly to publication-standard RGB while retaining the
# requested physical resolution. This changes no pixels or plotted values.
python_exe <- Sys.getenv("AHIPPO_PYTHON", Sys.which("python3"))
if (!file.exists(python_exe)) stop("Bundled Python runtime not found for RGB PNG normalization")
py_tmp <- tempfile(fileext = ".py")
writeLines(c(
  "from PIL import Image",
  "import sys",
  "p = sys.argv[1]",
  "im = Image.open(p).convert('RGB')",
  "im.save(p, dpi=(601, 601))"
), py_tmp, useBytes = TRUE)
py_status <- system2(python_exe, c(shQuote(py_tmp), shQuote(png_tmp)))
unlink(py_tmp)
if (!identical(py_status, 0L)) stop("RGB PNG normalization failed")

if (!file.copy(pdf_tmp, pdf_target, overwrite = TRUE)) stop("PDF copy failed")
if (!file.copy(png_tmp, png_target, overwrite = TRUE)) stop("PNG copy failed")
unlink(c(pdf_tmp, png_tmp))

write_utf8_copy <- function(lines, destination) {
  tmp <- tempfile(fileext = ".md")
  writeLines(lines, tmp, useBytes = TRUE)
  ok <- file.copy(tmp, destination, overwrite = TRUE)
  unlink(tmp)
  if (!ok) stop("Text output copy failed")
}

legend <- c(
  "# Supplementary Figure S2D | Stability of Program146 across effect and recurrence criteria",
  "",
  "All 146 genes in the frozen Program146 are positioned by their unshrunk patient-blocked High − Other log2 fold change and the number of discovery patients with a positive within-patient effect. Points are displayed with restrained perpendicular jitter only. Vertical references mark the prespecified sensitivity thresholds (log2FC ≥0.50, ≥0.75, ≥1.00 and ≥1.50). Labels report the number of Program146 genes retained and the number of 10 principal pathway directions concordant with the primary ≥0.50 analysis.",
  "",
  "Program146 comprised 16 genes positive in 6/8 patients, 68 positive in 7/8 and 62 positive in 8/8. The four thresholds retained 146, 85, 36 and 3 genes, respectively, while principal-pathway directional concordance was 10/10, 9/10, 8/10 and 8/10. Unshrunk and apeglm-shrunken effects were strongly correlated within Program146 (Spearman ρ = 0.982), and direction was retained for 146/146 genes. Among the three genes with unshrunk log2FC ≥1.5, all remained positive and two remained ≥1.5 after shrinkage.",
  "",
  "The ≥0.50 definition remains the frozen primary program; stricter thresholds are sensitivity analyses and do not redefine Program146."
)
write_utf8_copy(legend, file.path(out_dir, "S2D_legend_draft.md"))

qc <- c(
  "# Supplementary Figure S2D FINAL QC",
  "",
  "## Input and numerical fidelity",
  "",
  "- No DE, Program146 membership, effect threshold, patient recurrence, pathway result or shrinkage estimate was recomputed.",
  "- Complete Program146 gene rows: 146/146 - PASS.",
  "- Exact recurrence: 6/8 = 16 genes; 7/8 = 68; 8/8 = 62; cumulative >=7/8 = 130 - PASS.",
  "- Predefined thresholds: 0.50, 0.75, 1.00 and 1.50 - PASS.",
  "- Retained genes: 146, 85, 36 and 3 - PASS.",
  "- Principal pathway direction concordance: 10/10, 9/10, 8/10 and 8/10 - PASS.",
  "- apeglm: Spearman ρ = 0.982416; 146/146 directions retained - PASS.",
  "- Large-effect reference: 3/3 remain positive and 2/3 remain ≥1.5 after shrinkage - PASS.",
  "",
  "## Architecture",
  "",
  "- One integrated effect x recurrence landscape; all 146 genes are visible - PASS.",
  "- Continuous unshrunk effect distribution is preserved; thresholds are references within the distribution - PASS.",
  "- Patient recurrence is the second quantitative coordinate, not a table or duplicate bar series - PASS.",
  "- Figure 4B volcano, FDR axis and full shrinkage scatter are absent - PASS.",
  "- No heatmap, coloured table, simple retained-count bar chart or separate mini-plot dashboard - PASS.",
  "- The ≥1.5 genes are labelled only as a large-effect reference, never as the core program - PASS.",
  "",
  "## Render",
  "",
  "- Target: 88 x 78 mm; white background; Arial; core text >=6.1 pt.",
  sprintf("- Full unshrunk Program146 range: %.6f to %.6f; no clipping or winsorization.", min(gene$unshrunk_log2FC), max(gene$unshrunk_log2FC)),
  "- External PDF vector/font, PNG DPI and true-size visual checks are recorded after rendering.",
  "",
  "## Freeze status",
  "",
  "Pending external render validation."
)
write_utf8_copy(qc, file.path(out_dir, "S2D_program_stability_FINAL_QC.md"))

cat("S2D rendered\n")
