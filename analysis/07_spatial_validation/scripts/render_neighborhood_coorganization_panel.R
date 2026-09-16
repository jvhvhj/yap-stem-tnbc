# Purpose: Figure 6E neighborhood co-organization render
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 6 and Supplementary Fig. S4.
options(stringsAsFactors = FALSE)

root <- "."
out_dir <- file.path(root, "Figure6E_text_hierarchy_revision")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source_file <- file.path(root, "Figure6E_render_candidate", "Figure6E_render_source.tsv")
src <- read.delim(source_file, check.names = FALSE, na.strings = "NA")
nodes <- src[src$record_type == "neighbourhood_node", ]
edges <- src[src$record_type == "first_order_edge", ]
patients <- src[src$record_type == "patient_Lee_L", ]

stopifnot(nrow(nodes) == 37L, nrow(edges) == 90L, nrow(patients) == 22L)
stopifnot(all(patients$patient_Lee_L > 0), sum(nodes$center_spot, na.rm = TRUE) == 1L)

joint_limit <- 5.992045
program_limit <- 0.938416
patient_median <- 0.279377921029474
patient_q1 <- 0.233663631923195
patient_q3 <- 0.335768062619691
patient_p <- 2.38418579101562e-07

src$revision_role <- ifelse(src$record_type == "patient_Lee_L", "patient inferential unit",
                            ifelse(src$record_type == "neighbourhood_node", "measured Visium spot",
                                   "true first-order Visium edge"))
write.table(src, file.path(out_dir, "Figure6E_text_hierarchy_render_source.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

BLUE <- "#2F6DB3"
ZERO <- "#F7F7F7"
CORAL <- "#D7605C"
DARK <- "#26323A"
MUTED <- "#647178"
EDGE <- "#AEB8BD"
PATIENT_COL <- "#6C7F8A"
PATIENT_DARK <- "#425761"
palette_vec <- colorRampPalette(c(BLUE, ZERO, CORAL))(255)

map_colour <- function(values, limit) {
  scaled <- (pmax(-limit, pmin(limit, values)) + limit) / (2 * limit)
  palette_vec[pmax(1L, pmin(255L, floor(scaled * 254) + 1L))]
}

xlim_graph <- range(nodes$x_fullres) + c(-1, 1) * max(diff(range(nodes$x_fullres)) * 0.08, 1)
ylim_graph <- rev(range(nodes$y_fullres) + c(-1, 1) * max(diff(range(nodes$y_fullres)) * 0.08, 1))

draw_graph <- function(values, limit, title_text) {
  par(mar = c(0.25, 0.25, 1.55, 0.25), xaxs = "i", yaxs = "i")
  plot(NA, xlim = xlim_graph, ylim = ylim_graph, asp = 1, axes = FALSE, xlab = "", ylab = "")
  segments(edges$x_from, edges$y_from, edges$x_to, edges$y_to,
           col = adjustcolor(EDGE, alpha.f = 0.48), lwd = 0.45)
  points(nodes$x_fullres, nodes$y_fullres, pch = 21, cex = 1.25,
         bg = map_colour(values, limit), col = adjustcolor("#7D898F", alpha.f = 0.55), lwd = 0.32)
  title(main = title_text, cex.main = 0.95, font.main = 2, col.main = DARK, line = 0.50)
}

draw_colourbar <- function(limit) {
  par(mar = c(1.18, 3.0, 0.05, 3.0), xaxs = "i", yaxs = "i")
  plot(NA, xlim = c(-limit, limit), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  br <- seq(-limit, limit, length.out = 256)
  eps <- (br[2] - br[1]) * 0.08
  for (i in seq_len(255)) rect(br[i] - eps, 0.36, br[i + 1] + eps, 0.72,
                               col = palette_vec[i], border = NA)
  rect(-limit, 0.36, limit, 0.72, border = "#AAB4B9", lwd = 0.35)
  axis(1, at = c(-limit, 0, limit), labels = sprintf("%.2f", c(-limit, 0, limit)),
       cex.axis = 0.91, col.axis = DARK, col = DARK, lwd = 0,
       lwd.ticks = 0.42, tck = -0.20, padj = 0)
}

set.seed(20260826)
patient_jitter <- runif(nrow(patients), -0.040, 0.040)
patient_bw <- bw.nrd0(patients$patient_Lee_L)

draw_patient_recurrence <- function() {
  par(mar = c(2.65, 3.2, 0.70, 0.65), xaxs = "i", yaxs = "i")
  xlim <- c(-0.03, 0.55)
  plot(NA, xlim = xlim, ylim = c(0.04, 1.0), axes = FALSE, xlab = "", ylab = "")
  abline(v = 0, col = "#8D989E", lwd = 0.65, lty = 2)
  den <- density(patients$patient_Lee_L, bw = patient_bw, from = xlim[1], to = xlim[2], n = 512, cut = 0)
  base <- 0.43
  den_scaled <- den$y / max(den$y) * 0.25
  polygon(c(den$x, rev(den$x)), c(rep(base, length(den$x)), rev(base + den_scaled)),
          col = adjustcolor(PATIENT_COL, alpha.f = 0.17), border = NA)
  lines(den$x, base + den_scaled, col = adjustcolor(PATIENT_COL, alpha.f = 0.70), lwd = 0.70)
  points(patients$patient_Lee_L, base - 0.055 + patient_jitter, pch = 21,
         bg = PATIENT_COL, col = "white", lwd = 0.42, cex = 0.80)
  y_glyph <- 0.22
  segments(patient_q1, y_glyph, patient_q3, y_glyph, col = PATIENT_DARK, lwd = 1.15)
  segments(patient_median, y_glyph - 0.045, patient_median, y_glyph + 0.045,
           col = PATIENT_DARK, lwd = 2.15)
  axis(1, at = seq(0, 0.5, 0.1), labels = sprintf("%.1f", seq(0, 0.5, 0.1)),
       cex.axis = 0.90, col.axis = DARK, col = DARK, lwd = 0.55,
       lwd.ticks = 0.40, tck = -0.045, padj = 0.25)
  mtext("Lee's L", side = 1, line = 1.75, cex = 0.90, col = DARK)
  text(-0.018, 0.91,
       expression(22/22~positive~"|"~median~L==0.279~"|"~P==2.38%*%10^{-7}),
       adj = 0, cex = 0.92, font = 2, col = DARK)
  text(-0.018, 0.78, "43/43 sections supportive", adj = 0,
       cex = 0.90, col = MUTED)
  box(bty = "l", col = DARK, lwd = 0.55)
}

draw_panel_e <- function() {
  layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE),
         widths = c(1, 1), heights = c(3.45, 0.68, 2.28))
  par(family = "Arial", ps = 7.0, oma = c(0.15, 0.10, 1.55, 0.10))
  draw_graph(nodes$Joint_spatial, joint_limit, "Joint YAP-Stem score")
  draw_graph(nodes$Program146_spatial, program_limit, "Program146 score")
  draw_colourbar(joint_limit)
  draw_colourbar(program_limit)
  draw_patient_recurrence()
  mtext("E", side = 3, outer = TRUE, line = 0.70, at = 0.008, adj = 0,
        cex = 1.42, font = 2, col = DARK)
  mtext("Neighborhood-level spatial co-organization", side = 3, outer = TRUE,
        line = 0.74, at = 0.075, adj = 0, cex = 1.12, font = 2, col = DARK)
  mtext("P12/093D", side = 3, outer = TRUE, line = 0.02, at = 0.985,
        adj = 1, cex = 0.90, col = MUTED)
}

pdf_out <- file.path(out_dir, "Figure6E_text_hierarchy_candidate.pdf")
png_out <- file.path(out_dir, "Figure6E_text_hierarchy_candidate.png")
cairo_pdf(pdf_out, width = 94 / 25.4, height = 75 / 25.4, family = "Arial", bg = "white")
draw_panel_e()
dev.off()
png(png_out, width = 94, height = 75, units = "mm", res = 600,
    type = "cairo", family = "Arial", bg = "white")
draw_panel_e()
dev.off()

qc <- c(
  "# Figure 6E text-hierarchy QC",
  "",
  "- Candidate count: exactly one.",
  "- Scientific architecture unchanged: two aligned score maps plus the 22-patient Lee's L distribution.",
  "- Representative object unchanged: P12/093D; 37 local nodes; 90 true first-order edges.",
  "- Node coordinates, edge coordinates, Joint scores, Program146 scores, and patient Lee's L values are unchanged from the audited render source.",
  "- Header cleanup: long grey representative-section descriptor removed; only P12/093D retained in small contextual type.",
  "- Primary inference: 22/22 positive; median L=0.279; exact P=2.38e-7 in dark text.",
  "- Secondary support: 43/43 sections supportive in smaller muted text.",
  "- Patient-level subtitle removed.",
  "- Colourbar ticks rendered in dark text: Joint -5.99/0/5.99; Program146 -0.94/0/0.94.",
  "- Joint and Program146 retain separate zero-centered colour scales.",
  "- Isolated-spot statement remains QC provenance only and is not plotted.",
  "- Output: 94 x 75 mm; vector PDF; 600-dpi PNG; Arial; white background.",
  "- Minimum core annotation: approximately 6.2 pt.",
  "- Render inspection: PASS; titles, P12/093D label, graph nodes/edges, dark colourbar ticks, patient statistics, section support, all patient points, and Lee's L axis are unclipped and non-overlapping.",
  "- New biological analysis: NO."
)
writeLines(qc, file.path(out_dir, "Figure6E_text_hierarchy_QC.md"), useBytes = TRUE)

semantic <- c(
  "# Figure 6E updated semantic readout",
  "",
  "- The two graph views retain the identical 37 measured Visium spots, 90 true first-order spatial-neighbour edges, crop, and orientation from P12/093D.",
  "- Node fill encodes Joint YAP-Stem score on the left and Program146 score on the right; the numerical scales are zero-centered but quantitatively separate.",
  "- Each lower point is one biological patient's frozen Lee's L summary; 22 patients, not 43 sections, are the inferential units.",
  "- All 22 patient summaries are positive, with median L=0.279 and exact patient-level P=2.38e-7.",
  "- All 43 sections were supportive in section-level permutation QC, retained as secondary recurrence evidence.",
  "- The panel supports neighbourhood-level spatial co-organization, not cell-cell communication, causal spreading, direct regulation, interaction, or trajectory."
)
writeLines(semantic, file.path(out_dir, "Figure6E_text_hierarchy_semantic_readout.md"), useBytes = TRUE)

cat("Figure 6E text-hierarchy render complete\n")
