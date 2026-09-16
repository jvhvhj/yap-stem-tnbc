# Purpose: Render S3C state sensitivity
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Supplementary Fig. S1–S4.
# Single authorized S3C compact dual-sensitivity candidate.
# All medians/IQRs are READ from the audited source; no effect/quantile recalculation.
# R/grid is used for the complete scientific graphic. No density, CI or new tests.
options(stringsAsFactors = FALSE, digits = 17)
suppressWarnings(try(Sys.setlocale("LC_CTYPE", "Chinese_China.utf8"), silent = TRUE))
library(grid)
args <- commandArgs(trailingOnly = TRUE)
root_args <- args[!grepl("^--", args)]
root <- normalizePath(if (length(root_args)) root_args[1] else getwd(), winslash = "/", mustWork = TRUE)
out <- file.path(root, "Supplementary_S3C_compact_dual_sensitivity_TRUE_SIZE")
prefix <- "Supplementary_Figure_S3C_compact_dual_sensitivity"
src <- file.path(root, "Supplementary_S3ABC_redesign_gate/S3C_data_topology_audit.tsv")
source_sha256 <- "f478041e51c994f724cc10331ec7ab1d0d981d8073324e39697b2699848a88f9"
stopifnot(file.exists(src), dir.exists(out),
          tolower(unname(tools::md5sum(src))) == "9e70378b63a771aefd0099a01628be34")
paths <- setNames(file.path(out, paste0(prefix, c(".pdf", "_600dpi.png", "_plotdata.tsv", "_QC.md", "_readout.md"))),
                  c("pdf", "png", "data", "qc", "readout"))
if (!"--audit-only" %in% args && any(file.exists(paths)))
  stop("Existing candidate output: do not overwrite or generate another automatic version.")
protected_dirs <- list.dirs(root, full.names = TRUE, recursive = FALSE)
protected_dirs <- protected_dirs[grepl("^Supplementary_S3", basename(protected_dirs)) & basename(protected_dirs) != basename(out)]
protected_dirs <- c(protected_dirs, file.path(root, "Supplementary_Submission_Staging/01_FINAL_PANEL_SOURCES/Supplementary_Figure_S3"))
protected <- unlist(lapply(protected_dirs, list.files, full.names = TRUE, recursive = TRUE))
protected <- protected[!file.info(protected)$isdir]
before_md5 <- tools::md5sum(protected)

states <- c("G2.M", "S.G1", "Stress", "Hypoxia", "StressER", "Interferon", "HLA", "Basal", "EMT", "LumSec", "Cholesterol")
variants <- c("Raw Rule C", "Technical Rule C", "Raw Rule D")
summary_fields <- c("tested_n", "q25", "median", "q75")
d <- read.delim(src, check.names = FALSE, fileEncoding = "UTF-8", na.strings = "NA")
d_text <- read.delim(src, check.names = FALSE, fileEncoding = "UTF-8", colClasses = "character", na.strings = character())
stopifnot(all(c(summary_fields, "patient_id", "state", "state_display", "state_order", "variant", "testable", "effect_high_minus_other") %in% names(d)),
          nrow(d) == 2574L, length(unique(d$patient_id)) == 78L,
          !anyDuplicated(d[c("patient_id", "state", "variant")]),
          all(d$state_display == states[d$state_order]), setequal(unique(d$variant), variants),
          all(d$testable %in% 0:1), all(d$program_genes == 146), all(d$measurable_program_genes == 141),
          all(d$source_column == "effect_high_minus_other"),
          all(is.finite(d$effect_high_minus_other[d$testable == 1])),
          all(is.na(d$effect_high_minus_other[d$testable == 0])),
          all(is.na(d$positive_flag[d$testable == 0])), all(is.na(d$negative_flag[d$testable == 0])))
groups <- list()
summary_text <- list()
retention <- numeric(11)
subset_n <- 0L
for (i in seq_along(states)) {
  rc <- d[d$state_order == i & d$variant == "Raw Rule C", ]
  tc <- d[d$state_order == i & d$variant == "Technical Rule C", ]
  rd <- d[d$state_order == i & d$variant == "Raw Rule D", ]
  stopifnot(identical(sort(rc$patient_id[rc$testable == 1]), sort(tc$patient_id[tc$testable == 1])),
            all(rd$patient_id[rd$testable == 1] %in% rc$patient_id[rc$testable == 1]),
            sum(rd$testable) < sum(rc$testable))
  cr <- d_text[d_text$state_order == as.character(i) & d_text$variant == "Raw Rule C", ]
  dr <- d_text[d_text$state_order == as.character(i) & d_text$variant == "Raw Rule D" & d_text$testable == "1", ]
  stopifnot(identical(dr$effect_high_minus_other, cr$effect_high_minus_other[match(dr$patient_id, cr$patient_id)]))
  subset_n <- subset_n + nrow(dr)
  for (variant in variants) {
    z <- d[d$state_order == i & d$variant == variant, ]
    zt <- d_text[d_text$state_order == as.character(i) & d_text$variant == variant, ]
    stopifnot(nrow(z) == 78L, all(z$summary_reproduced == "PASS within source precision"))
    formal <- vapply(summary_fields, function(field) {
      v <- unique(z[[field]]); stopifnot(length(v) == 1L, is.finite(v)); v
    }, numeric(1))
    exact <- vapply(summary_fields, function(field) {
      v <- unique(zt[[field]]); stopifnot(length(v) == 1L); v
    }, character(1))
    stopifnot(formal["tested_n"] == sum(z$testable), formal["q25"] <= formal["median"], formal["median"] <= formal["q75"])
    groups[[paste(i, variant, sep = ":")]] <- formal
    summary_text[[paste(i, variant, sep = ":")]] <- exact
  }
  retention[i] <- groups[[paste(i, "Raw Rule D", sep = ":")]]["tested_n"] /
                  groups[[paste(i, "Raw Rule C", sep = ":")]]["tested_n"]
}
totals <- as.integer(vapply(variants, function(v) sum(d$testable[d$variant == v]), numeric(1)))
negative_n <- as.integer(vapply(variants, function(v) sum(d$negative_flag[d$variant == v], na.rm = TRUE), numeric(1)))
stopifnot(identical(totals, c(490L, 490L, 325L)), subset_n == 325L,
          identical(negative_n, c(5L, 18L, 2L)), sum(is.na(d$effect_high_minus_other)) == 1269L,
          all(retention > 0 & retention < 1))
all_summary_values <- unlist(lapply(groups, function(g) g[c("q25", "median", "q75")]))
effect_range <- c(-.03, .25)
stopifnot(min(all_summary_values) > effect_range[1], max(all_summary_values) < effect_range[2])
cat(sprintf("INPUT_PASS; 33 summary groups; totals=%s; identical Rule D/C effects=%d; NT=%d\n", paste(totals,collapse="/"), subset_n, sum(is.na(d$effect_high_minus_other))))
cat(sprintf("SUMMARY_RANGE=[%.17g,%.17g]; AXES=[-.03,.25]\n", min(all_summary_values), max(all_summary_values)))
if ("--audit-only" %in% args) quit(save = "no", status = 0)

# Prespecified single geometry, fixed before rendering.
WIDTH_MM <- 130; HEIGHT_MM <- 86
FONT <- "Arial"; INK <- "#27323A"; BLUE <- "#6EB6E4"; CORAL <- "#E7837D"; CHARCOAL <- "#4E5961"
left_x <- c(22, 64); right_x <- c(71, 113); retention_x <- c(116, 125.5)
row_y <- 64.3 - (seq_along(states)-1) * 4.65
axis_y <- 12
dodge <- .72
xmap <- function(v, domain) domain[1] + diff(domain) * (v-effect_range[1])/diff(effect_range)
stopifnot(diff(left_x) == diff(right_x), diff(retention_x) <= 10,
          tail(row_y,1) - dodge - .7 > axis_y, head(row_y,1)+dodge+.7 < 67)

mm <- function(x) unit(x, "mm")
txt <- function(label,x,y,size=7,just="centre",face="plain") {
  grid.text(label,x=mm(x),y=mm(y),just=just,
            gp=gpar(fontfamily=FONT,fontsize=size,fontface=face,col=INK,lineheight=.95))
}
seg <- function(x0,y0,x1,y1,col=INK,lwd=.5,lty=1) {
  grid.segments(mm(x0),mm(y0),mm(x1),mm(y1),gp=gpar(col=col,lwd=lwd,lty=lty,lineend="butt"))
}
circle <- function(x,y,r=.49,col=BLUE) grid.circle(mm(x),mm(y),r=mm(r),gp=gpar(fill=col,col=NA))
diamond <- function(x,y,r=.57) grid.polygon(mm(x+c(0,r,0,-r)),mm(y+c(r,0,-r,0)),gp=gpar(fill=CORAL,col=NA))
square <- function(x,y,s=.93) grid.rect(mm(x),mm(y),width=mm(s),height=mm(s),gp=gpar(fill="white",col=CHARCOAL,lwd=.6))
draw_panel <- function() {
  grid.newpage()
  pushViewport(viewport(x=0,y=0,width=mm(WIDTH_MM),height=mm(HEIGHT_MM),just=c("left","bottom"),clip="off"))
  grid.rect(gp=gpar(fill="white",col=NA))
  txt("C",1.8,82.3,10,"left","bold")
  txt("State-stratified Program146 robustness",6.3,82.3,8.7,"left","bold")
  circle(24.5,75.7)
  txt("Raw Rule C",26.3,75.7,7,"left")
  diamond(57.3,75.7)
  txt("Technical Rule C",59.2,75.7,7,"left")
  square(98.5,75.7)
  txt("Raw Rule D",100.4,75.7,7,"left")
  txt("Technical adjustment",mean(left_x),69.5,7.8,"centre","bold")
  txt("Stricter support",mean(right_x),69.5,7.8,"centre","bold")
  txt("Retained",mean(retention_x),69.5,6.6)
  for (domain in list(left_x,right_x)) {
    seg(xmap(0,domain),axis_y,xmap(0,domain),66.5,"#B8BEC3",.4,"22")
  }
  for (i in seq_along(states)) {
    y <- row_y[i]
    c <- groups[[paste(i,"Raw Rule C",sep=":")]]
    t <- groups[[paste(i,"Technical Rule C",sep=":")]]
    r <- groups[[paste(i,"Raw Rule D",sep=":")]]
    txt(states[i],19.5,y,7,"right")
    # Descriptive median displacement connectors, behind both median/IQR pairs.
    seg(xmap(c["median"],left_x),y+dodge,xmap(t["median"],left_x),y-dodge,"#D4D8DA",.32)
    seg(xmap(c["median"],right_x),y+dodge,xmap(r["median"],right_x),y-dodge,"#D4D8DA",.32)
    seg(xmap(c["q25"],left_x),y+dodge,xmap(c["q75"],left_x),y+dodge,BLUE,.72)
    circle(xmap(c["median"],left_x),y+dodge)
    seg(xmap(t["q25"],left_x),y-dodge,xmap(t["q75"],left_x),y-dodge,CORAL,.72)
    diamond(xmap(t["median"],left_x),y-dodge)
    light_blue <- adjustcolor(BLUE,alpha.f=.58)
    seg(xmap(c["q25"],right_x),y+dodge,xmap(c["q75"],right_x),y+dodge,light_blue,.60)
    circle(xmap(c["median"],right_x),y+dodge,.42,light_blue)
    seg(xmap(r["q25"],right_x),y-dodge,xmap(r["q75"],right_x),y-dodge,CHARCOAL,.62)
    square(xmap(r["median"],right_x),y-dodge)
    # Full-length hairline is 100%; filled length is existing Rule-D/Rule-C tested N.
    seg(retention_x[1],y-dodge,retention_x[2],y-dodge,"#D9DDE1",.35)
    grid.rect(x=mm(retention_x[1]),y=mm(y-dodge),width=mm(diff(retention_x)*retention[i]),height=mm(.75),
              just=c("left","centre"),gp=gpar(fill=adjustcolor(CHARCOAL,alpha.f=.58),col=NA))
  }
  ticks <- seq(0,.25,.05)
  labels <- c("0","0.05","0.10","0.15","0.20","0.25")
  for (domain in list(left_x,right_x)) {
    seg(domain[1],axis_y,domain[2],axis_y,INK,.45)
    for (j in seq_along(ticks)) {
      xx <- xmap(ticks[j],domain)
      seg(xx,axis_y,xx,axis_y-.65,INK,.45)
      txt(labels[j],xx,axis_y-2.2,7)
    }
  }
  txt("Program146 High \u2212 Other effect",mean(c(left_x[1],right_x[2])),4.1,7)
  popViewport()
}

# Retain the complete audit as byte-identical Source Data; do not author a reduced table.
stopifnot(file.copy(src,paths["data"],overwrite=FALSE),
          identical(unname(tools::md5sum(src)),unname(tools::md5sum(paths["data"]))))
summary_lines <- c("| Cancer state | Variant | Tested N | Q25 | Median | Q75 |", "|---|---|---:|---:|---:|---:|")
counts_lines <- c("| Cancer state | Rule C | Rule D | D/C retention |", "|---|---:|---:|---:|")
for (i in seq_along(states)) {
  for (variant in variants) {
    s <- summary_text[[paste(i,variant,sep=":")]]
    summary_lines <- c(summary_lines,sprintf("| %s | %s | %s | %s | %s | %s |",states[i],variant,s[1],s[2],s[3],s[4]))
  }
  counts_lines <- c(counts_lines,sprintf("| %s | %d | %d | %.12f |",states[i],
      groups[[paste(i,"Raw Rule C",sep=":")]]["tested_n"],groups[[paste(i,"Raw Rule D",sep=":")]]["tested_n"],retention[i]))
}
qc <- c("# Supplementary Figure S3C - compact dual-sensitivity QC", "",
  "Numerical and export checks are separate from visual inspection at the final display size.", "",
  "## Data and numerical checks", "",
  paste0("Sole scientific input: Supplementary_S3ABC_redesign_gate/S3C_data_topology_audit.tsv. SHA-256: `",source_sha256,"`."),
  "- 2574 unique patient-state-variant records; 78 patients x 11 states x 3 variants. 0 duplicate keys.",
  "- Existing testable flags: Raw Rule C 490; Technical Rule C 490; Raw Rule D 325. Raw C and Technical C retain identical patient sets within each state.",
  "- All 325 Rule D raw effects exactly match the corresponding Rule C source strings. Rule D remains a strict subset in every state; no eligibility flags are generated or changed.",
  "- All 33 sets of tested_n/q25/median/q75 are internally constant within their source group and copied directly into the plot objects. No quantile/effect/score/rho recalculation.",
  "- 1269 NT effects remain NA in the byte-identical plotting source. Existing negative patient counts remain 5/18/2 for Raw C/Technical C/Raw D. No patient effects were deleted from Source Data.",
  "- Exactly 44 median positions and 44 uncapped Q25-Q75 segments: Raw C appears once in each of two comparisons as their common reference, Technical C once, Raw D once per state.",
  "- 11 retention lengths = 9.5 mm x (existing Rule D tested_n / existing Rule C tested_n). No percent rounding enters geometry; no 11-row printed counts/percentages.",
  sprintf("- All displayed summary values lie in [%.17g, %.17g]. Both x axes use [-0.03, +0.25] and exactly 42 mm. No summary clipping or winsorization.",min(all_summary_values),max(all_summary_values)),
  "- Patient-level extrema do not determine a summary-only axis. The figure does not claim to display individual tails; the legend explicitly directs readers to complete patient-level Source Data.",
  "- No new hypothesis tests, CIs, whisker caps, KDE, density, violin, half-eye, re-tertiling, gene selection, rescaling or patient filtering.", "",
  "## Exact displayed source summaries", "",summary_lines,"","## Exact support retention", "",counts_lines,"",
  "## Display layout", "",
  "- 130 x 86 mm, white, Arial. State order: G2.M, S.G1, Stress, Hypoxia, StressER, Interferon, HLA, Basal, EMT, LumSec, Cholesterol. One state-label column.",
  "- Two adjacent equal-width domains: Technical adjustment (22-64 mm), Stricter support (71-113 mm). Shared effect range and ticks; one zero reference per domain.",
  "- Filled blue Raw C circles (#6EB6E4), coral Technical C diamonds (#E7837D), open charcoal Rule D squares (#4E5961). Right Raw C reference uses the identical blue at reduced alpha and slightly smaller size. No new biological color mapping.",
  "- Slight vertical dodge +/-0.72 mm. Neutral 0.32-lwd median connectors are descriptive, not patient trajectories or tests. Horizontal segments are IQRs, not CIs.",
  "- Retained is the sole retention heading; 9.5-mm support marks start 3 mm beyond the right effect domain, aligned with the Rule D row. No separate retention panel/axis.",
  "- Fixed row pitch 4.65 mm. State labels 7 pt, axis/ticks 7 pt, legend 7 pt, domain headings 7.8 pt, retention heading 6.6 pt, main title 8.7 pt bold, panel letter 10 pt bold.",
  "- No density tails, point clouds, boxes, large explanatory legends, gray title bands, decorative cards, significance text or large blank center region.","",
  "## Design references", "",
  "Forest-family design reference: nature_figure_skill-master/archetypes/clinical-forest/plot.R and out/ref.png were inspected. Only aligned estimates, uncapped horizontal intervals and a quantitative reference axis are borrowed. The clinical template's HR=1, logarithmic axis, 95% CI and significance colors are NOT applicable and are not run. The dual-sensitivity display preserves separate IQR and retention encodings.",
  "Wanjie visual-language/plot-recipes/review-checklist/palette.json supplied the existing blue/coral/neutral semantics, consistent typography and non-decorative shared alignment. User constraints override auto-stars, gray facet bands and other irrelevant recipes. All primitives are R/grid, not manual image editing.",
  "The plot is implemented in R/grid with explicit color and typography settings; it does not depend on the legacy nature_theme.R or save_nature interface. The source TSV is copied without numerical modification.","",
  "## Export and protection", "",
  sprintf("Before rendering, %d pre-existing S3-related files are protected by content hashes; final comparisons follow after export.",length(protected)),
  "Only this new candidate directory is written. No S3A/S3B/S3D revision, A4 assembly, old-file deletion or automatic alternative.",
  "Post-export vector, resolution, text bounds, glyph-coordinate and actual-size review: pending.","",
  "## True-size visual inspection", "",
  "Inspect the exported 11-row display for readable estimates, intervals and retention annotations. The technical and support sensitivity comparisons use the same Raw C reference. Visual inspection does not change the underlying estimates or testability criteria.")
readout <- c("# Supplementary Figure S3C - state-stratified Program146 robustness", "",
  "Standalone dual-sensitivity display from the recorded source results.", "",
  "## Figure legend draft", "",
  "**C, State-stratified Program146 robustness.** State-stratified Program146 effects are summarized as patient-level medians and interquartile ranges. The left domain compares the primary Raw Rule C analysis with the prespecified technical adjustment. The right domain compares Raw Rule C with the stricter-support Rule D subset. Rule D requires at least 10 High and 10 Other cells per patient-state, compared with at least 5 and 5 for Rule C; adjacent retention marks indicate the fraction of Rule C-testable patients retained under Rule D. Complete patient-level distributions, including negative effects and non-testable entries, are reported in Supplementary Table S3 and Source Data.",
  "Filled blue circles indicate Raw Rule C; coral diamonds indicate Technical Rule C; open charcoal squares indicate Raw Rule D. The lighter blue right-hand reference is the same Raw Rule C summary, not another analysis. Thin connecting segments link medians descriptively. Horizontal intervals are Q25-Q75, not confidence intervals. Zero marks no High-Other score difference. The full-length neutral retention guide represents 100% of Rule C-testable patients.",
  "Rule D is a stricter-support subset using unchanged raw patient effects, not a third independent model. Technical adjustment uses the prespecified cohort-specific technical covariates. Non-testable entries are NA, never zero. No significance testing or equality claim is made by this display.", "",
  "## Numerical and interpretive boundary", "",
  "All 33 audited summary medians are positive; this does not mean every patient effect is positive. Existing negative counts are 5 Raw Rule C, 18 Technical Rule C and 2 Raw Rule D. No individual negative effects or NT entries are removed from the accompanying byte-identical audit copy.",
  "Technical adjustment produces lower state-level medians in all 11 comparisons. The Rule D comparison reflects changing support/eligible subset rather than refitting another model. The summaries do not establish state independence or statistical equivalence.",
  "The 490/490/325 testable entries are patient-state observations, not counts of distinct patients across the 11 states. Exact support counts are listed below for legend/Table S3 linkage.","",counts_lines,"",
  "## Display and provenance", "",
  "130 x 86 mm; both summary axes [-0.03, +0.25]; 600 dpi PNG plus vector PDF. Source medians/IQRs are not recomputed. Full patient extrema are intentionally not depicted in this summary-only graphic; all data remain in Source Data.",
  paste0("Sole source: Supplementary_S3ABC_redesign_gate/S3C_data_topology_audit.tsv; SHA-256 ",source_sha256,"."),
  "The local forest archetype informs row/estimate alignment only; its clinical statistics are not imported. Wanjie informs the established palette and restrained hierarchy. No exact published dual-sensitivity template is claimed.","",
  "Stop after the single candidate. S3A, S3B, S3D, prior S3C outputs and composites remain unchanged.")
writeLines(qc,paths["qc"],useBytes=TRUE)
writeLines(readout,paths["readout"],useBytes=TRUE)

normalize_pdf_page <- function(path) {
  py <- Sys.getenv("S3C_PYTHON",unset=Sys.getenv("AHIPPO_PYTHON", Sys.which("python3")))
  code <- paste("import sys,io","from pathlib import Path","from pypdf import PdfReader,PdfWriter",
    "from pypdf.generic import RectangleObject","p=Path(sys.argv[1])","r=PdfReader(io.BytesIO(p.read_bytes()))",
    "old=r.pages[0].get_contents().get_data()","w=PdfWriter();w.clone_document_from_reader(r);w.pdf_header='%PDF-1.7'",
    "box=RectangleObject([0,0,130/25.4*72,86/25.4*72])","w.pages[0].mediabox=box;w.pages[0].cropbox=box",
    "b=io.BytesIO();w.write(b)","assert PdfReader(io.BytesIO(b.getvalue())).pages[0].get_contents().get_data()==old",
    "p.write_bytes(b.getvalue())","print('130x86_MM_VECTOR_STREAM_UNCHANGED')",sep=";")
  result <- system2(py,c("-c",shQuote(code),shQuote(normalizePath(path,winslash="/"))),stdout=TRUE,stderr=TRUE)
  if (!is.null(attr(result,"status")) && attr(result,"status")!=0) stop(paste(result,collapse="\n"))
  cat(paste(result,collapse="\n"),"\n")
}
stopifnot(capabilities("cairo"))
png(paths["png"],width=WIDTH_MM,height=HEIGHT_MM,units="mm",res=600,type="cairo",bg="white")
draw_panel(); invisible(dev.off())
cairo_pdf(paths["pdf"],width=WIDTH_MM/25.4,height=HEIGHT_MM/25.4,family=FONT,bg="white")
draw_panel(); invisible(dev.off())
normalize_pdf_page(paths["pdf"])
stopifnot(identical(before_md5,tools::md5sum(protected)))
cat(sprintf("Panel rendered; %d protected files unchanged\n",length(protected)))
