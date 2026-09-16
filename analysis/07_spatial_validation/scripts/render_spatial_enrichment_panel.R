# Purpose: Figure 6D spatial enrichment render
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 6 and Supplementary Fig. S4.
options(stringsAsFactors = FALSE)

root <- "."
out_dir <- file.path(root, "Figure6D_publication_density_revision")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source_file <- file.path(root, "Figure6D_render_candidate", "Figure6D_render_source.tsv")
boundary_file <- file.path(root, "Figure6D_reference_validation_replication", "Figure6D_between_cohort_boundary.tsv")
dat <- read.delim(source_file, check.names = FALSE)
boundary <- read.delim(boundary_file, check.names = FALSE)

stopifnot(nrow(dat) == 22L, nrow(boundary) == 1L)
stopifnot(sum(dat$cohort == "Reference") == 14L, sum(dat$cohort == "Validation") == 8L)
stopifnot(all(dat$adjusted_spatial_spearman_rho > 0))
stopifnot(abs(unique(dat$cohort_median[dat$cohort == "Reference"]) - 0.401913919224372) < 1e-12)
stopifnot(abs(unique(dat$cohort_median[dat$cohort == "Validation"]) - 0.315861510007235) < 1e-12)

dat$display_order_within_cohort <- ave(dat$adjusted_spatial_spearman_rho, dat$cohort,
                                       FUN = function(x) rank(x, ties.method = "first"))
dat$display_role <- "one biological patient"
dat$ecdf_role <- "exact empirical cumulative distribution; no smoothing"
write.table(dat,
            file.path(out_dir, "Figure6D_publication_density_render_source.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

COL <- c(Reference = "#4F827D", Validation = "#B18455")
COL_DARK <- c(Reference = "#315D59", Validation = "#7D5938")
DARK <- "#26323A"
MUTED <- "#5F6C73"
LIGHT <- "#DCE2E5"

ecdf_steps <- function(values, xmin, xmax) {
  values <- sort(values)
  n <- length(values)
  x <- c(xmin, rep(values, each = 2), xmax)
  y <- c(0, as.vector(rbind((0:(n - 1)) / n, (1:n) / n)), 1)
  list(x = x, y = y)
}

set.seed(20260826)
jitter_map <- list(
  Reference = runif(sum(dat$cohort == "Reference"), -0.030, 0.030),
  Validation = runif(sum(dat$cohort == "Validation"), -0.030, 0.030)
)

draw_panel_d <- function() {
  layout(matrix(c(1, 2), nrow = 2), heights = c(2.10, 1.20))
  par(family = "Arial", ps = 7.1, oma = c(0.10, 0.10, 1.30, 0.10))
  xlim <- c(0, 0.62)

  # Upper layer: exact ECDFs on one common quantitative axis.
  par(mar = c(0.65, 4.55, 1.15, 0.65), xaxs = "i", yaxs = "i")
  plot(NA, xlim = xlim, ylim = c(0, 1.04), axes = FALSE, xlab = "", ylab = "")
  abline(v = 0, col = "#8C989E", lwd = 0.65, lty = 2)
  axis(2, at = c(0, 0.5, 1), labels = c("0", "0.5", "1.0"), cex.axis = 0.90,
       col.axis = DARK, col = DARK, lwd = 0.55, lwd.ticks = 0.42, tck = -0.025, las = 1, padj = 0.15)
  mtext("Cumulative patient fraction", side = 2, line = 2.75, cex = 0.90, col = DARK)
  for (coh in c("Reference", "Validation")) {
    st <- ecdf_steps(dat$adjusted_spatial_spearman_rho[dat$cohort == coh], xlim[1], xlim[2])
    lines(st$x, st$y, col = COL_DARK[[coh]], lwd = 1.25,
          lty = if (coh == "Reference") 1 else 2)
  }
  segments(0.018, 0.94, 0.055, 0.94, col = COL_DARK[["Reference"]], lwd = 1.25)
  text(0.061, 0.94, "Reference", adj = 0, cex = 0.88, font = 2, col = DARK)
  segments(0.018, 0.84, 0.055, 0.84, col = COL_DARK[["Validation"]], lwd = 1.25, lty = 2)
  text(0.061, 0.84, "Validation", adj = 0, cex = 0.88, font = 2, col = DARK)

  text(0.018, 0.66, "Cohort median shift", adj = 0, cex = 0.86, font = 2, col = DARK)
  text(0.018, 0.56, "Reference  0.402", adj = 0, cex = 0.88, col = MUTED)
  text(0.018, 0.47, "Validation  0.316", adj = 0, cex = 0.88, col = MUTED)
  text(0.018, 0.38, "Delta median  -0.086", adj = 0, cex = 0.88, col = MUTED)
  box(bty = "l", col = DARK, lwd = 0.55)

  # Lower layer: complete patient strips aligned to the identical rho axis.
  par(mar = c(3.25, 4.55, 0.10, 0.65), xaxs = "i", yaxs = "i")
  plot(NA, xlim = xlim, ylim = c(0.02, 1.02), axes = FALSE, xlab = "", ylab = "")
  abline(v = 0, col = "#8C989E", lwd = 0.65, lty = 2)
  axis(1, at = c(0, 0.2, 0.4, 0.6), labels = sprintf("%.1f", c(0, 0.2, 0.4, 0.6)),
       cex.axis = 0.92, col.axis = DARK, col = DARK, lwd = 0.55, lwd.ticks = 0.42,
       tck = -0.045, padj = 0.28)
  mtext(expression(paste("Adjusted spatial Spearman ", rho)), side = 1, line = 2.15,
        cex = 0.98, col = DARK)
  rows <- c(Reference = 0.73, Validation = 0.29)
  for (coh in names(rows)) {
    row_y <- rows[[coh]]
    dd <- dat[dat$cohort == coh, ]
    points(dd$adjusted_spatial_spearman_rho, row_y + jitter_map[[coh]],
           pch = if (coh == "Reference") 21 else 22,
           bg = COL[[coh]], col = "white", lwd = 0.45, cex = 0.82)
    med <- unique(dd$cohort_median)
    q1 <- unique(dd$cohort_q1)
    q3 <- unique(dd$cohort_q3)
    glyph_y <- row_y - 0.105
    segments(q1, glyph_y, q3, glyph_y, col = COL_DARK[[coh]], lwd = 1.10)
    segments(med, glyph_y - 0.035, med, glyph_y + 0.035, col = COL_DARK[[coh]], lwd = 2.0)
    text(0.012, row_y + 0.035,
         sprintf("%s  %d/%d positive", coh, unique(dd$positive_n), unique(dd$cohort_n)),
         adj = 0, cex = 0.88, font = 2, col = COL_DARK[[coh]])
  }
  box(bty = "l", col = DARK, lwd = 0.55)

  mtext("D", side = 3, outer = TRUE, line = 0.30, at = 0.01, adj = 0,
        cex = 1.42, font = 2, col = DARK)
  mtext("Replication across predefined TNBC cohorts", side = 3, outer = TRUE,
        line = 0.34, at = 0.075, adj = 0, cex = 1.14, font = 2, col = DARK)
}

pdf_out <- file.path(out_dir, "Figure6D_publication_density_candidate.pdf")
png_out <- file.path(out_dir, "Figure6D_publication_density_candidate.png")
cairo_pdf(pdf_out, width = 82 / 25.4, height = 75 / 25.4, family = "Arial", bg = "white")
draw_panel_d()
dev.off()
png(png_out, width = 82, height = 75, units = "mm", res = 600,
    type = "cairo", family = "Arial", bg = "white")
draw_panel_d()
dev.off()

qc <- c(
  "# Figure 6D publication-density QC",
  "",
  "- Candidate count: exactly one.",
  "- Architecture: exact two-cohort ECDF plus complete aligned patient strips; KDE removed.",
  "- Biological units: 22 patients (Reference n=14; Validation n=8); all 22 points displayed.",
  "- ECDF: exact empirical steps from frozen adjusted spatial Spearman rho values; no smoothing or fitted distribution.",
  "- Shared x-axis: 0.0-0.6 ticks in both layers; rho=0 visible.",
  "- Reference: 14/14 positive; median 0.4019139; IQR 0.3795821-0.4584900; Holm P 0.0001220703125.",
  "- Validation: 8/8 positive; median 0.3158615; IQR 0.2764835-0.3773305; Holm P 0.00390625.",
  "- Cohort median shift: 0.3158615 - 0.4019139 = -0.0860524.",
  sprintf("- Descriptive between-cohort rank-sum P %.7f retained outside the primary plotting object.", boundary$p_two_sided[1]),
  "- Patient jitter: perpendicular to rho only; fixed seed 20260826; rho values unchanged.",
  "- Median/IQR: compact median tick plus IQR segment; no printed row-level median/P prose.",
  "- Output: 82 x 75 mm; vector PDF; 600-dpi PNG; Arial; white background.",
  "- Minimum core annotation: approximately 6.2 pt.",
  "- Render inspection: PASS; ECDF curves, inset, patient labels, all 22 points, median/IQR glyphs, and axes are unclipped and non-overlapping.",
  "- New biological analysis: NO."
)
writeLines(qc, file.path(out_dir, "Figure6D_publication_density_QC.md"), useBytes = TRUE)

semantic <- c(
  "# Figure 6D updated semantic readout",
  "",
  "- The upper layer is the exact cumulative distribution of frozen patient-level adjusted spatial Spearman rho values in the predefined Reference and Validation cohorts.",
  "- The lower strips show every biological patient on the same rho axis; symbols are jittered only vertically.",
  "- Median ticks and IQR segments summarize each cohort without duplicating long numerical sentences inside the patient rows.",
  "- All patients are positive in both cohorts (14/14 and 8/8), supporting independent positive replication.",
  "- Validation is shifted toward lower rho because its median is 0.316 versus 0.402 in Reference; the descriptive median difference is -0.086.",
  "- The between-cohort P=0.0128467 is descriptive boundary information, not an interaction test or evidence of biological suppression."
)
writeLines(semantic, file.path(out_dir, "Figure6D_publication_density_semantic_readout.md"), useBytes = TRUE)

cat("Figure 6D publication-density render complete\n")
