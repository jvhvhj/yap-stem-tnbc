# Purpose: Render S3B patient-level distribution
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Supplementary Fig. S1–S4.
# Integrated distribution display from existing results; no scientific reanalysis.
# R/grid owns the complete shared category geometry and all scientific graphics.
options(stringsAsFactors = FALSE, digits = 17)
suppressWarnings(try(Sys.setlocale("LC_CTYPE", "Chinese_China.utf8"), silent = TRUE))
library(grid)
args <- commandArgs(trailingOnly = TRUE)
root_arg <- args[!grepl("^--", args)]
root <- normalizePath(if (length(root_arg)) root_arg[1] else getwd(), winslash = "/", mustWork = TRUE)
out <- file.path(root, "Supplementary_S3B_multitrack_TRUE_SIZE")
prefix <- "Supplementary_Figure_S3B_multitrack"
src <- file.path(root, "Supplementary_S3ABC_redesign_gate/S3B_data_topology_audit.tsv")
ref <- file.path(root, "Supplementary_S3B_multitrack_architecture_gate/S3B_multitrack_reference_audit.md")
stopifnot(file.exists(src), file.exists(ref), dir.exists(out))
SOURCE_SHA256 <- "6176580d98cbd8df74ef842689ac1bb42dc977babe72d76fbbc32ed9f9e1b9ad"
stopifnot(tolower(unname(tools::md5sum(src))) == "cb38705a772af7c940fc659fbff9a9f4")
paths <- setNames(file.path(out, paste0(prefix, c(".pdf", "_600dpi.png", "_plotdata.tsv", "_QC.md", "_readout.md"))),
                  c("pdf", "png", "data", "qc", "readout"))
if (!"--audit-only" %in% args && any(file.exists(paths)))
  stop("Output directory already exists. Choose a new output directory to preserve existing files.")
protected_dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
protected_dirs <- protected_dirs[grepl("^Supplementary_S3", basename(protected_dirs)) & basename(protected_dirs) != basename(out)]
protected_dirs <- c(protected_dirs, file.path(root, "Supplementary_Submission_Staging/01_FINAL_PANEL_SOURCES/Supplementary_Figure_S3"))
protected <- unlist(lapply(protected_dirs, list.files, recursive = TRUE, full.names = TRUE))
protected <- protected[!file.info(protected)$isdir]
before_md5 <- tools::md5sum(protected)

# Read the sole scientific source; group summaries are extracted, not re-estimated.
d <- read.delim(src, check.names = FALSE, fileEncoding = "UTF-8")
states <- c("G2.M", "S.G1", "Stress", "Hypoxia", "StressER", "Interferon", "HLA", "Basal", "EMT", "LumSec", "Cholesterol")
models <- c("Raw", "Technical-adjusted")
summary_fields <- c("tested_n", "positive_n", "negative_n", "zero_n", "min", "q25", "median", "q75", "max")
required <- c("patient_id", "state", "state_display", "state_order", "method", "rho", "positive_flag", "negative_flag", "source_column", summary_fields)
stopifnot(all(required %in% names(d)), nrow(d) == 1470L, all(is.finite(d$rho)),
  identical(sort(unique(d$state_order)), 1:11), setequal(unique(d$method), models),
  all(d$state_display == states[d$state_order]),
  sum(d$method == "Raw") == 735L, sum(d$method == "Technical-adjusted") == 735L,
  length(unique(paste(d$patient_id, d$state))) == 735L,
  !anyDuplicated(d[c("patient_id", "state", "method")]),
  all(d$source_column[d$method == "Raw"] == "rho_raw_within_state"),
  all(d$source_column[d$method == "Technical-adjusted"] == "rho_technical_frozen_patient_fit"),
  all(as.logical(d$positive_flag) == (d$rho > 0)),
  all(as.logical(d$negative_flag) == (d$rho < 0)))
groups <- list()
checks <- list()
for (i in seq_along(states)) {
  raw <- d[d$state_order == i & d$method == "Raw", ]
  adj <- d[d$state_order == i & d$method == "Technical-adjusted", ]
  stopifnot(identical(sort(raw$patient_id), sort(adj$patient_id)))
  for (model in models) {
    g <- d[d$state_order == i & d$method == model, ]
    formal <- vapply(summary_fields, function(nm) {
      u <- unique(g[[nm]]); stopifnot(length(u) == 1L, is.finite(u)); u
    }, numeric(1))
    # Read-only count/range reconciliation is not a new hypothesis test or summary fit.
    stopifnot(formal["tested_n"] == nrow(g), formal["positive_n"] == sum(g$rho > 0),
      formal["negative_n"] == sum(g$rho < 0), formal["zero_n"] == sum(g$rho == 0),
      formal["min"] == min(g$rho), formal["max"] == max(g$rho), formal["negative_n"] > 0,
      formal["min"] <= formal["q25"], formal["q25"] <= formal["median"],
      formal["median"] <= formal["q75"], formal["q75"] <= formal["max"])
    checks[[length(checks) + 1L]] <- data.frame(state = states[i], model = model,
      as.list(formal), bw_nrd0 = bw.nrd0(g$rho), check.names = FALSE)
    groups[[paste(i, model, sep = ":")]] <- list(values = g$rho, summary = formal)
  }
}
audit <- do.call(rbind, checks)
# The only median computed here is the explicitly approved DISPLAY bandwidth.
COMMON_BW <- median(audit$bw_nrd0)
stopifnot(abs(COMMON_BW - 0.073118533228091) < 1e-14, nrow(audit) == 22L, all(audit$median > 0),
          sum(d$rho[d$method == "Raw"] < 0) == 94L,
          sum(d$rho[d$method == "Technical-adjusted"] < 0) == 145L)
cat("NUMERIC_GATE_PASS; 1470 rho; 735 paired units; 22 source medians; 22 source fractions\n")
cat(sprintf("COMMON_BW=%.15f; OBSERVED_RANGE=[%.15f, %.15f]\n", COMMON_BW, min(d$rho), max(d$rho)))
if ("--audit-only" %in% args) { print(audit); quit(save = "no", status = 0) }

# All geometry and display parameters fixed BEFORE opening a graphics device.
WIDTH_MM <- 140
HEIGHT_MM <- 104
WMAX_MM <- 3.1
FONT <- "Arial"
INK <- "#29343B"
COLORS <- c(Raw = "#6EB6E4", `Technical-adjusted` = "#E7837D")
rho_range <- c(-.72, .88)
median_range <- c(0, .40)
positive_range <- c(0, 1)
axis_left <- 18
axis_right <- 137
# This single vector is used in EVERY track and the only state-label lane.
centers <- axis_left + (seq_along(states) - .5) * (axis_right - axis_left) / length(states)
model_offset <- c(Raw = -1.1, `Technical-adjusted` = 1.1)
middle_y <- c(28, 80)
top_y <- c(82, 97)
bottom_y <- c(14, 26)
map_value <- function(v, domain, span) span[1] + diff(span) * (v - domain[1]) / diff(domain)
rho_y <- function(v) map_value(v, rho_range, middle_y)
median_y <- function(v) map_value(v, median_range, top_y)
positive_y <- function(v) map_value(v, positive_range, bottom_y)
stopifnot(min(d$rho) > rho_range[1], max(d$rho) < rho_range[2],
  all(audit$median >= median_range[1] & audit$median <= median_range[2]),
  all(diff(centers) > 2 * WMAX_MM), min(centers) - WMAX_MM > axis_left,
  max(centers) + WMAX_MM < axis_right)
rho_grid <- seq(rho_range[1], rho_range[2], length.out = 4097L)
for (key in names(groups)) {
  g <- groups[[key]]
  # Equal-patient Gaussian display density; no filtered points or state-specific bandwidth.
  den <- rowMeans(dnorm(outer(rho_grid, g$values, "-"), sd = COMMON_BW))
  width <- WMAX_MM * den / max(den)
  q <- g$summary[c("q25", "q75")]
  in_iqr <- approx(rho_grid, width, xout = seq(q[1], q[2], length.out = 101L))$y
  g$half_width <- width
  g$iqr_offset <- min(.75, .4 * min(in_iqr))
  stopifnot(g$iqr_offset > .2, abs(max(width) - WMAX_MM) < 1e-12)
  groups[[key]] <- g
}
fmt <- function(x) formatC(x, format = "f", digits = 15)
summary_lines <- c("| Cancer state | Model | Positive/tested | Negative | Q25 | Median | Q75 |",
  "|---|---|---:|---:|---:|---:|---:|",
  sprintf("| %s | %s | %d/%d | %d | %.12f | %.12f | %.12f |", audit$state, audit$model,
    audit$positive_n, audit$tested_n, audit$negative_n, audit$q25, audit$median, audit$q75))
counts_lines <- c("| Cancer state | Raw positive/tested | Technical-adjusted positive/tested |",
  "|---|---|---|", vapply(seq_along(states), function(i) {
    r <- groups[[paste(i, "Raw", sep = ":")]]$summary
    a <- groups[[paste(i, "Technical-adjusted", sep = ":")]]$summary
    sprintf("| %s | %d/%d | %d/%d |", states[i], r["positive_n"], r["tested_n"], a["positive_n"], a["tested_n"])
  }, character(1)))
bw_lines <- c("| Cancer state | Model | bw.nrd0 |", "|---|---|---:|",
  sprintf("| %s | %s | %s |", audit$state, audit$model, fmt(audit$bw_nrd0)))
geometry_lines <- c("| State | Shared center (mm) | Raw median/bar offset (mm) | Technical median/bar offset (mm) |",
  "|---|---:|---:|---:|", sprintf("| %s | %.12f | -1.1 | +1.1 |", states, centers))
qc <- c("# Supplementary S3B - integrated multi-track true-size QC", "",
  "Numerical and export checks are separate from visual inspection at the final display size.", "",
  "## Input and nonmodification gate", "",
  "- Sole scientific input: Supplementary_S3ABC_redesign_gate/S3B_data_topology_audit.tsv.",
  paste0("- Source SHA-256: ", SOURCE_SHA256, "."),
  "- 11 states; 735 patient-state units; 735 Raw + 735 Technical-adjusted = 1470 finite rho values; 74 distinct patients.",
  "- 22 source median fields and 22 source positive_n/tested_n fractions. Repeated summary fields are unique within each group.",
  "- Zero missing rho, zero duplicate patient-state-model keys, identical Raw/Technical patient sets in all 11 states.",
  "- Technical source field is rho_technical_frozen_patient_fit; no within-state refit substituted.",
  "- Raw negatives: 94; Technical-adjusted negatives: 145. Every state/model group includes negative values, all used in its density.",
  sprintf("- Exact observed range: [%s, %s].", fmt(min(d$rho)), fmt(max(d$rho))),
  "- Source median/Q25/Q75 are directly extracted, NOT recomputed. Count and extrema reconciliation only; no new scientific estimates.",
  "- Plotdata is an exact byte copy, including existing source provenance. No score, rho, state, threshold, eligibility, P, CI or hypothesis test recalculation.", "",
  "## Display lock", "",
  sprintf("- One common Gaussian display bandwidth: %s = median of the 22 bw.nrd0 values. Locked before rendering, no tuning.", fmt(COMMON_BW)),
  "- Rho viewport is exactly [-0.72, +0.88], not the rejected two-bandwidth viewport. Every observed value is inside; no observation clipping or winsorization.",
  "- The mathematical Gaussian tails beyond the finite display viewport are not drawn; they are not additional observations. Patient values and density inputs are unchanged.",
  "- Maximum half-width = 3.1 mm in all 22 groups; width encodes within-group density shape only, not N/cells/absolute cross-state density.",
  "- Top: 22 median symbols at source medians; per-state Raw/Technical connector only, no across-state connector.",
  "- Middle: all observations enter split densities; source Q25-Q75 thin dark-grey segments; duplicate internal median ticks omitted.",
  "- Bottom: exactly 22 narrow bars at source positive_n/tested_n; common true 0-100% axis and zero baseline, no per-bar number labels.",
  "- Common colors across all tracks: Raw #6EB6E4, Technical-adjusted #E7837D. Density alpha 0.82, thin same-color edge; no shaded positive region.",
  "- No patient cloud, boxplot, whisker, stars, P values, CI, table lane, decorative boxes or grey explanatory paragraph.", "",
  "## Geometry and typography lock", "",
  "- Direct 140 x 104 mm canvas, white background; no render-larger-then-shrink.",
  "- One shared state center vector; top symbols and bottom bars use identical +/-1.1 mm model offsets around it; middle halves meet at that same center.",
  "- Top y=82..97 mm (15 mm); middle y=28..80 mm (52 mm); bottom y=14..26 mm (12 mm); two 2 mm gaps.",
  "- Top scale 0..0.40; middle -0.72..+0.88; positive fraction 0..100%. No independent state scales.",
  "- Arial: letter 10 pt bold; title 8.7 pt bold; all axis titles, ticks, model legend and state labels 7 pt. Minimum 7 pt.",
  "- Full official state names appear once, at bottom, rotated 45 degrees. One shared model legend at the top.", "",
  geometry_lines, "", "## Source summaries (outside the plotting area)", "", summary_lines, "",
  "## Deterministic bandwidth audit", "", bw_lines, "",
  "## Design references", "",
  "Approved A-prime is a CLOSE combined adaptation, not an exact published template. Huang et al., Nature Communications 2023 Fig.4e-f supplies category-aligned proportions and paired distribution grammar; Xu et al., Nature Genetics 2025 Fig.4a supplies flush multi-track alignment hierarchy only, NOT its heatmap.",
  "Design references include the local Nature archetypes and Wanjie visual-language and plot-recipes resources. The display uses shared geometry, restrained colors and a consistent hierarchy. Heatmap, automatic-significance and sample-weight encodings are not used. Scientific graphics are rendered with R/grid and exported as vector PDF.",
  "Reference audit: Supplementary_S3B_multitrack_architecture_gate/S3B_multitrack_reference_audit.md.", "",
  "## Post-export inspection", "", "Awaiting file-level inspection and visual assessment of this single render.", "",
  "## True-size visual inspection", "",
  "1. Are all 22 medians immediately identifiable?", "2. Does the middle show meaningful distribution heterogeneity, not generic violins?",
  "3. Are negative tails visible without large blank space?", "4. Do proportions add useful directional-support information?",
  "5. Do bars remain secondary?", "6. Are all three tracks perfectly aligned by state?", "7. Does it read as one scientific object?",
  "8. Is it materially richer than violin-only?", "9. Is it elegant rather than crowded?", "10. Is it worthy of a high-impact Supplementary Figure?", "",
  "Inspect the exported panel at its final display size. Bandwidth, axes, tracks and fonts are fixed in the rendering code.")
readout <- c("# Supplementary S3B - within-state YAP-Stem coupling", "",
  "Standalone within-state patient-distribution display from the recorded source results.", "",
  "## Scientific reading and boundaries", "",
  "All 22 source state/model medians are positive, but negative patient observations remain in every distribution. The top shows central location, the middle patient heterogeneity/shape/negative tails, and the bottom directional patient support. These are complementary descriptive views of the SAME data, not three independent validations.",
  "G2.M has higher source medians than Stress under both models. Technical-adjusted medians are lower than Raw medians in these states; no new significance or equivalence claim is made.",
  "The within-state top connector links two group medians. It is not a patient trajectory, a median of paired differences or a confidence interval. Positive fraction is not a significance rate, FDR or cell fraction. No claim of every-patient positivity, state independence or causal cancer-state effects is made.", "",
  "## Figure legend draft", "",
  "Within-state YAP\u2013Stem coupling was evaluated across 11 independently annotated cancer states. The upper track shows Raw and Technical-adjusted state-level median Spearman correlations; lines connect the two descriptive medians within each state. Split violins in the middle show the complete corresponding patient-level distributions, with interquartile ranges indicated internally. The lower track shows the proportion of tested patients with positive correlations under each model. All testable observations, including negative correlations, were retained. Positive/tested patient counts are provided in Supplementary Table S3.",
  "Each patient has equal weight within a state/model distribution. One common Gaussian bandwidth and one maximum half-width are used; width expresses within-group shape only. Technical adjustment was performed within each cohort using the prespecified cohort-specific technical covariates.", "",
  "## Positive/tested counts for the legend and table", "", counts_lines, "",
  "Counts are recorded here without modifying Supplementary Table S3; its cross-reference can be checked during later assembly.", "",
  "## Display and provenance", "",
  sprintf("Direct size: 140 x 104 mm. Top median scale 0..0.40; middle rho scale -0.72..+0.88; bottom fraction scale 0..100%%. Common bandwidth %s; maximum half-width 3.1 mm.", fmt(COMMON_BW)),
  "Fixed state order: G2.M, S.G1, Stress, Hypoxia, StressER, Interferon, HLA, Basal, EMT, LumSec, Cholesterol.",
  "Sole source: Supplementary_S3ABC_redesign_gate/S3B_data_topology_audit.tsv. The output plotdata retains its exact bytes. Technical field: rho_technical_frozen_patient_fit.",
  "Architecture reference: [Huang et al., Nature Communications 2023, Fig.4e-f](https://www.nature.com/articles/s41467-023-36310-9/figures/4), plus the alignment-only principle in [Xu et al., Nature Genetics 2025, Fig.4a](https://www.nature.com/articles/s41588-025-02128-y/figures/4). CLOSE adaptation, not an exact template; no biological interpretation, test or smoothing parameter copied.", "",
  "## Output scope", "", "Only the standalone S3B panel is rendered; figure-page assembly is separate.")
stopifnot(file.copy(src, paths["data"], overwrite = FALSE),
  identical(unname(tools::md5sum(src)), unname(tools::md5sum(paths["data"])) ))
writeLines(qc, paths["qc"], useBytes = TRUE)
writeLines(readout, paths["readout"], useBytes = TRUE)
cat("PRE_RENDER_LOCK_WRITTEN; CANVAS=140x104; VIEWPORT=-0.72,0.88; WMAX=3.1\n")

mm <- function(x) unit(x, "mm")
txt <- function(label, x, y, size = 7, just = "centre", face = "plain", rot = 0) {
  grid.text(label, x = mm(x), y = mm(y), just = just, rot = rot,
    gp = gpar(fontfamily = FONT, fontsize = size, fontface = face, col = INK, lineheight = .95))
}
seg <- function(x0, y0, x1, y1, col = INK, lwd = .5) {
  grid.segments(mm(x0), mm(y0), mm(x1), mm(y1), gp = gpar(col = col, lwd = lwd, lineend = "butt"))
}
symbol <- function(cx, cy, model) {
  if (model == "Raw") grid.circle(mm(cx), mm(cy), r = mm(.47), gp = gpar(fill = COLORS[model], col = NA))
  else grid.polygon(mm(cx + c(0, .60, 0, -.60)), mm(cy + c(.60, 0, -.60, 0)), gp = gpar(fill = COLORS[model], col = NA))
}
axis <- function(span, values, map, labels) {
  seg(axis_left, span[1], axis_left, span[2], INK, .45)
  for (j in seq_along(values)) {
    seg(axis_left - .7, map(values[j]), axis_left, map(values[j]), INK, .45)
    txt(labels[j], axis_left - 1.5, map(values[j]), 7, "right")
  }
}
draw_panel <- function() {
  grid.newpage()
  pushViewport(viewport(width = mm(WIDTH_MM), height = mm(HEIGHT_MM), x = 0, y = 0, just = c("left", "bottom"), clip = "off"))
  grid.rect(gp = gpar(fill = "white", col = NA))
  txt("B", 1.8, 100.5, 10, "left", "bold")
  txt("Within-state YAP\u2013Stem coupling", 6.4, 100.5, 8.7, "left", "bold")
  symbol(86, 100.5, "Raw"); txt("Raw", 87.8, 100.5, 7, "left")
  symbol(99, 100.5, "Technical-adjusted"); txt("Technical-adjusted", 100.9, 100.5, 7, "left")
  seg(axis_left, rho_y(0), axis_right, rho_y(0), "#C0C4C7", .40)
  for (i in seq_along(states)) {
    cx <- centers[i]
    r <- groups[[paste(i, "Raw", sep = ":")]]$summary
    a <- groups[[paste(i, "Technical-adjusted", sep = ":")]]$summary
    seg(cx + model_offset["Raw"], median_y(r["median"]),
        cx + model_offset["Technical-adjusted"], median_y(a["median"]), "#B6BABD", .4)
    for (model in models) {
      g <- groups[[paste(i, model, sep = ":")]]
      side <- if (model == "Raw") -1 else 1
      symbol(cx + model_offset[model], median_y(g$summary["median"]), model)
      ox <- cx + side * g$half_width
      grid.polygon(mm(c(cx, ox, cx)), mm(c(rho_y(rho_grid[1]), rho_y(rho_grid), rho_y(tail(rho_grid, 1)))),
        gp = gpar(fill = adjustcolor(COLORS[model], alpha.f = .82), col = NA))
      grid.lines(mm(ox), mm(rho_y(rho_grid)), gp = gpar(col = COLORS[model], lwd = .32, lineend = "round"))
      sx <- cx + side * g$iqr_offset
      seg(sx, rho_y(g$summary["q25"]), sx, rho_y(g$summary["q75"]), adjustcolor(INK, alpha.f = .8), .55)
      fraction <- g$summary["positive_n"] / g$summary["tested_n"]
      grid.rect(x = mm(cx + model_offset[model]), y = mm(bottom_y[1]), width = mm(.70),
        height = mm(positive_y(fraction) - bottom_y[1]), just = c("centre", "bottom"), gp = gpar(fill = COLORS[model], col = NA))
    }
    txt(states[i], cx, 12.1, 7, c("right", "centre"), rot = 45)
  }
  axis(top_y, c(0, .2, .4), median_y, c("0", "0.2", "0.4"))
  axis(middle_y, c(-.6, -.3, 0, .3, .6), rho_y, c("\u22120.6", "\u22120.3", "0", "0.3", "0.6"))
  axis(bottom_y, c(0, .5, 1), positive_y, c("0", "50", "100"))
  seg(axis_left, bottom_y[1], axis_right, bottom_y[1], INK, .45)
  txt("Median\nSpearman \u03c1", 4.4, mean(top_y), 7, rot = 90)
  txt("Within-patient Spearman \u03c1", 4.4, mean(middle_y), 7, rot = 90)
  txt("Positive\npatients (%)", 4.4, mean(bottom_y), 7, rot = 90)
  popViewport()
}

# Windows Cairo can round PDF page dimensions. Only normalize page boxes;
# retain the complete vector drawing content verbatim.
normalize_pdf_page <- function(path) {
  py <- Sys.getenv("S3B_PYTHON", unset = Sys.getenv("AHIPPO_PYTHON", Sys.which("python3")))
  stopifnot(file.exists(py))
  code <- paste("import sys,io", "from pathlib import Path", "from pypdf import PdfReader,PdfWriter",
    "from pypdf.generic import RectangleObject", "p=Path(sys.argv[1])",
    "r=PdfReader(io.BytesIO(p.read_bytes()))", "old=r.pages[0].get_contents().get_data()",
    "w=PdfWriter();w.clone_document_from_reader(r);w.pdf_header='%PDF-1.7'",
    "box=RectangleObject([0,0,140/25.4*72,104/25.4*72])", "w.pages[0].mediabox=box;w.pages[0].cropbox=box",
    "b=io.BytesIO();w.write(b)", "assert PdfReader(io.BytesIO(b.getvalue())).pages[0].get_contents().get_data()==old",
    "p.write_bytes(b.getvalue())", "print('140x104_MM_PAGE_BOX_VECTOR_STREAM_UNCHANGED')", sep = ";")
  result <- system2(py, c("-c", shQuote(code), shQuote(normalizePath(path, winslash = "/"))), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(result, "status")) && attr(result, "status") != 0) stop(paste(result, collapse = "\n"))
  cat(paste(result, collapse = "\n"), "\n")
}
stopifnot(capabilities("cairo"))
cairo_pdf(paths["pdf"], width = WIDTH_MM/25.4, height = HEIGHT_MM/25.4, family = FONT, bg = "white")
draw_panel()
invisible(dev.off())
normalize_pdf_page(paths["pdf"])
png(paths["png"], width = WIDTH_MM, height = HEIGHT_MM, units = "mm", res = 600, type = "cairo", bg = "white")
draw_panel()
invisible(dev.off())
stopifnot(identical(before_md5, tools::md5sum(protected)))
cat(sprintf("Panel rendered; %d protected files unchanged\n", length(protected)))
