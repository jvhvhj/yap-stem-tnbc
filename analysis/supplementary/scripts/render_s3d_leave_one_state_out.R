# Purpose: Render S3D leave-one-cancer-state-out robustness
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Supplementary Fig. S1–S4.
# Supplementary Figure S3D: one approved 140 x 86 mm split-distribution render.
# Run from project root:
# Rscript --vanilla Supplementary_S3D_split_distribution_TRUE_SIZE/Supplementary_Figure_S3D_split_distribution_plot.R
# Only descriptive plotting summaries/KDE are computed. No patient effects or tests.
suppressPackageStartupMessages(library(grid))
options(stringsAsFactors = FALSE, digits = 17)
if (.Platform$OS.type == "windows") {
  suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese_China.utf8"))
}
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".", winslash = "/", mustWork = TRUE)
outdir <- file.path(root, "Supplementary_S3D_split_distribution_TRUE_SIZE")
prefix <- "Supplementary_Figure_S3D_split_distribution"
file_out <- function(suffix) file.path(outdir, paste0(prefix, suffix))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
# Prevent an unreviewed second image or silent overwrite.
stopifnot(!file.exists(file_out(".pdf")), !file.exists(file_out("_600dpi.png")),
          capabilities("cairo"))
source_dir <- file.path(root, "Supplementary_S3D_distributional_architecture_gate")
source_file <- file.path(source_dir, "S3D_patient_paired_displacement_audit.tsv")
reference_file <- file.path(source_dir, "S3D_distributional_architecture_candidates.md")
stopifnot(file.exists(source_file), file.exists(reference_file))
read_tsv <- function(f) read.delim(f, colClasses = "character", check.names = FALSE,
                                  na.strings = character(), fileEncoding = "UTF-8")
src <- read_tsv(source_file)
protected_dirs <- file.path(root, c(
  "Supplementary_Submission_Staging/01_FINAL_PANEL_SOURCES/Supplementary_Figure_S3",
  "Supplementary_S3D_matched_influence_FINAL",
  "Supplementary_S3D_matched_reference_gate",
  "Supplementary_S3D_distributional_architecture_gate"))
protected <- unique(c(source_file, reference_file,
  unlist(lapply(protected_dirs, list.files, recursive = TRUE, full.names = TRUE))))
protected <- protected[file.exists(protected) & !dir.exists(protected)]
protected_md5 <- tools::md5sum(protected)

states <- c("G2.M", "S.G1", "Stress", "Hypoxia", "StressER", "Interferon",
            "HLA", "Basal", "EMT", "LumSec", "Cholesterol")
metrics <- c("rho_raw", "rho_technical", "program_raw_high_minus_other",
             "program_technical_high_minus_other")
keys <- paste(src$omitted_state, src$domain, src$patient_id, sep = "|")
stopifnot(nrow(src) == 3296L, !anyDuplicated(keys),
          setequal(unique(src$cancer_state), states),
          setequal(unique(src$domain), metrics),
          all(src$patient_pairing == "PASS"),
          all(is.finite(as.numeric(src$paired_delta))),
          all(is.finite(as.numeric(src$full_effect))),
          all(is.finite(as.numeric(src$omission_effect))))
stopifnot(identical(as.integer(table(factor(src$domain, levels = metrics))),
                    c(825L, 825L, 823L, 823L)))
for (st in states) {
  stopifnot(length(unique(src$omitted_state[src$cancer_state == st])) == 1L)
  for (pair in list(metrics[1:2], metrics[3:4])) {
    # Within each domain only; never intersect the two biological domains.
    a <- src$patient_id[src$cancer_state == st & src$domain == pair[1]]
    b <- src$patient_id[src$cancer_state == st & src$domain == pair[2]]
    stopifnot(setequal(a, b))
  }
}
group_id <- paste(src$cancer_state, src$domain, sep = "|")
idx <- unlist(lapply(states, function(st) paste(st, metrics, sep = "|")))
groups <- lapply(idx, function(k) which(group_id == k))
names(groups) <- idx
stopifnot(length(groups) == 44L, all(lengths(groups) > 0L))
num <- function(x) as.numeric(x)
fmt <- function(x) format(x, digits = 17, scientific = FALSE, trim = TRUE)
summaries <- lapply(groups, function(ii) {
  x <- num(src$paired_delta[ii])
  qq <- quantile(x, c(0, .10, .25, .50, .75, .90, 1), type = 7, names = FALSE)
  expected <- vapply(c("minimum", "q10", "q25", "median_paired_delta",
                       "q75", "q90", "maximum"), function(col) {
    v <- unique(src[[col]][ii]); stopifnot(length(v) == 1L); num(v)
  }, numeric(1))
  stopifnot(max(abs(qq - expected)) < 2e-15,
            all(num(src$paired_n[ii]) == length(ii)),
            all(num(src$omission_absolute_tested_n[ii]) == length(ii)),
            sum(num(src$omission_effect[ii]) > 0) == num(src$omission_absolute_positive_n[ii[1]]),
            all(src$omission_absolute_positive[ii] ==
                ifelse(num(src$omission_effect[ii]) > 0, "TRUE", "FALSE")))
  data.frame(cancer_state = src$cancer_state[ii[1]], domain = src$domain[ii[1]],
    N = length(ii), minimum = qq[1], Q10 = qq[2], Q25 = qq[3], median = qq[4],
    Q75 = qq[5], Q90 = qq[6], maximum = qq[7],
    group_bw_nrd0 = stats::bw.nrd0(x),
    positive_absolute = num(src$omission_absolute_positive_n[ii[1]]))
})
summary_df <- do.call(rbind, summaries)
rownames(summary_df) <- NULL
summary_df$biological_domain <- ifelse(grepl("^rho_", summary_df$domain),
                                      "YAP-Stem coupling", "Program146 effect")
summary_df$model <- ifelse(summary_df$domain %in% metrics[c(1,3)], "Raw", "Technical-adjusted")
summary_df$state_order <- match(summary_df$cancer_state, states)
# Prespecified bandwidth rule. These two bandwidths are immutable after this calculation.
bandwidths <- tapply(summary_df$group_bw_nrd0, summary_df$biological_domain, median)
stopifnot(length(bandwidths) == 2L, all(is.finite(bandwidths)), all(bandwidths > 0))
summary_df$bandwidth_used <- unname(bandwidths[summary_df$biological_domain])
p57 <- src[src$cancer_state == "Stress" & src$domain == metrics[4] & src$patient_id == "P57", ]
stopifnot(nrow(p57) == 1L,
  p57$omission_effect == "-0.013021691225597543",
  p57$paired_delta == "-0.0398367950697579",
  p57$omission_absolute_positive_n == "73",
  p57$omission_absolute_tested_n == "74")
coupling_n <- c(75,75,74,74,75,76,76,75,75,75,75)
program_n <- c(75,74,74,74,75,76,75,75,75,75,75)
stopifnot(identical(as.integer(summary_df$N[summary_df$domain == metrics[1]]), as.integer(coupling_n)),
          identical(as.integer(summary_df$N[summary_df$domain == metrics[4]]), as.integer(program_n)))
stopifnot(sum(!as.logical(summary_df$median == 0) &
  sign(summary_df$median) != sign(vapply(groups, function(ii)
    num(src$prior_difference_of_matched_medians[ii[1]]), numeric(1)))) == 12L)

# Finite numerical Gaussian support covers ALL observed values plus four bandwidths
# on both sides. Each numerical density integrates to one. No observed tail omitted.
trapezoid <- function(x,y) sum(diff(x) * (head(y,-1L) + tail(y,-1L))/2)
domains <- c("YAP-Stem coupling", "Program146 effect")
domain_specs <- lapply(domains, function(dom) {
  ii <- which(summary_df$biological_domain == dom)
  rr <- range(c(summary_df$minimum[ii], summary_df$maximum[ii]))
  bw <- unname(bandwidths[dom])
  padding <- max(4*bw, .02*diff(rr))
  list(observed = rr, limits = rr + c(-padding, padding), bandwidth = bw)
})
names(domain_specs) <- domains
densities <- vector("list", nrow(summary_df))
for (i in seq_len(nrow(summary_df))) {
  dom <- summary_df$biological_domain[i]
  spec <- domain_specs[[dom]]
  x <- num(src$paired_delta[groups[[i]]])
  xx <- sort(unique(c(seq(spec$limits[1], spec$limits[2], length.out = 4096),
                      x, summary_df$Q25[i], summary_df$median[i], summary_df$Q75[i])))
  # Evaluate the Gaussian KDE directly; avoids FFT approximations and preserves area.
  yy <- rowMeans(dnorm(outer(xx, x, "-") / spec$bandwidth)) / spec$bandwidth
  pre_area <- trapezoid(xx, yy)
  yy <- yy/pre_area
  stopifnot(abs(trapezoid(xx,yy)-1) < 1e-12,
            min(x) > spec$limits[1], max(x) < spec$limits[2])
  # Close each polygon only beyond min/max + 4 bandwidths.
  # Integration/storage still uses the full common domain.
  support <- xx >= min(x)-4*spec$bandwidth & xx <= max(x)+4*spec$bandwidth
  densities[[i]] <- list(x=xx, y=yy, support=support, area=trapezoid(xx,yy),
                         unnormalized_area=pre_area)
}
height_max_mm <- 2.02
for (dom in domains) {
  ii <- which(summary_df$biological_domain == dom)
  peak <- max(vapply(densities[ii], function(z) max(z$y), numeric(1)))
  domain_specs[[dom]]$height_mm_per_density <- height_max_mm / peak
}
summary_df$density_area <- vapply(densities, function(z) z$area, numeric(1))
summary_df$gaussian_area_before_normalization <- vapply(densities, function(z) z$unnormalized_area, numeric(1))
summary_df$height_mm_per_density <- vapply(summary_df$biological_domain,
  function(dom) domain_specs[[dom]]$height_mm_per_density, numeric(1))

# Preserve every original decimal string; append descriptive/display metadata only.
plotdata <- src
mi <- match(group_id, paste(summary_df$cancer_state,summary_df$domain,sep="|"))
extra <- summary_df[mi,c("N","minimum","Q25","median","Q75","maximum",
  "group_bw_nrd0","bandwidth_used","biological_domain","model","state_order",
  "density_area","gaussian_area_before_normalization","height_mm_per_density")]
names(extra) <- paste0("render_", names(extra))
for (col in names(extra)) {
  plotdata[[col]] <- if (is.numeric(extra[[col]])) fmt(extra[[col]]) else extra[[col]]
}
plotdata$render_source_sha256 <- "F86773B0D6E970DA28775B58C60F31BB3877CC350299BACD7F3E1390EE1797F2"
plotdata$render_summary_semantics <- "Patient paired delta distribution; not difference of cohort medians"
write.table(plotdata, file_out("_plotdata.tsv"), sep="\t", row.names=FALSE,
            quote=FALSE, na="", fileEncoding="UTF-8")
reread <- read_tsv(file_out("_plotdata.tsv"))
stopifnot(identical(reread[names(src)], src))
bandwidth_lock <- c(
  "# S3D split-distribution true-size render QC", "",
  "Numerical source checks completed; true-size inspection is separate.",
  "Bandwidths locked before the first figure device was opened.",
  paste0("- Coupling bandwidth: ", fmt(bandwidths["YAP-Stem coupling"])),
  paste0("- Program146 bandwidth: ", fmt(bandwidths["Program146 effect"])),
  "- Rule: median of 22 group-specific stats::bw.nrd0 values per biological domain.",
  "- No alternate bandwidth versions permitted.", "")
writeLines(bandwidth_lock,file_out("_QC.md"), useBytes=TRUE)
writeLines(c("# S3D scientific readout", "", "Pre-render bandwidth lock",
  bandwidth_lock[5:8]),file_out("_readout.md"),useBytes=TRUE)
cat("PRE-RENDER LOCK\n")
print(bandwidths)

# All geometry is in millimetres, authored at the requested size.
W <- 140; H <- 86
palette <- c("Raw"="#6EB6E4", "Technical-adjusted"="#E7837D")
ink <- "#202B33"; quantile_ink <- "#37434B"
row_y <- 22 + (0:10)*4.7
plot_x <- list(c(24,69), c(79,123.5))
names(plot_x) <- domains
n_x <- 73
recurrence_x <- 131.5
text_records <- list()
text_mm <- function(label,x,y,size=7,bold=FALSE,just="centre",record=TRUE) {
  grob <- textGrob(label,x=unit(x,"mm"),y=unit(H-y,"mm"),just=just,
    gp=gpar(fontfamily="Arial",fontsize=size,fontface=if(bold)"bold" else "plain",
            col=ink,lineheight=1.08))
  width <- convertWidth(grobWidth(grob),"mm",valueOnly=TRUE)
  height <- convertHeight(grobHeight(grob),"mm",valueOnly=TRUE)
  lo <- if(just=="left") x else if(just=="right") x-width else x-width/2
  if(record) text_records[[length(text_records)+1L]] <<-
    data.frame(label=gsub("\n"," / ",label),x0=lo,x1=lo+width,
               y0=y-height/2,y1=y+height/2,size=size)
  if(lo < .3 || lo+width > W-.3 || y-height/2 < .3 || y+height/2 > H-.3)
    stop("Text outside canvas: ",label)
  grid.draw(grob)
}
line_mm <- function(x1,y1,x2,y2,col=ink,pt=.4,lty=1) {
  grid.segments(unit(x1,"mm"),unit(H-y1,"mm"),unit(x2,"mm"),unit(H-y2,"mm"),
                gp=gpar(col=col,lwd=pt/.75,lty=lty,lineend="butt"))
}
# grid lwd uses R's 1/96 inch, so pt / .75 gives the requested physical point width.
x_position <- function(x,dom) {
  xx <- plot_x[[dom]]; lim <- domain_specs[[dom]]$limits
  xx[1] + (x-lim[1])/diff(lim)*diff(xx)
}
legend_half <- function(x,y,upper,col) {
  tt <- seq(0,pi,length.out=64)
  xx <- x+1.5*cos(tt)
  yy <- y+if(upper) -0.8*sin(tt) else 0.8*sin(tt)
  grid.polygon(unit(c(xx,xx[1]),"mm"),unit(H-c(yy,yy[1]),"mm"),
               gp=gpar(fill=col,col=NA))
}
quantile_geometry <- vector("list",44)
for (i in seq_len(nrow(summary_df))) {
  s <- summary_df[i,]; dom <- s$biological_domain; de <- densities[[i]]
  hfun <- approxfun(de$x,de$y * domain_specs[[dom]]$height_mm_per_density)
  qx <- seq(s$Q25,s$Q75,length.out=128)
  min_iqr_height <- min(hfun(qx))
  offset <- min(.34, min_iqr_height*.40)
  median_height <- hfun(s$median)
  quantile_geometry[[i]] <- c(offset=offset,tick_low=offset*.18,
        tick_high=min(median_height*.86, offset+.42))
  stopifnot(offset>0, quantile_geometry[[i]]["tick_high"]>quantile_geometry[[i]]["tick_low"])
}
draw_panel <- function() {
  grid.newpage()
  pushViewport(viewport(width=unit(W,"mm"),height=unit(H,"mm")))
  grid.rect(gp=gpar(fill="white",col=NA))
  text_records <<- list()
  text_mm("D",2.4,3.8,9.7,TRUE,"left")
  text_mm("Patient-level influence of single cancer-state omission",8.1,3.8,8.4,TRUE,"left")
  legend_half(25.6,9.6,TRUE,palette["Raw"])
  text_mm("Raw (upper)",28,9.4,6.8,FALSE,"left")
  legend_half(48,9.2,FALSE,palette["Technical-adjusted"])
  text_mm("Technical-adjusted (lower)",50.5,9.4,6.8,FALSE,"left")
  text_mm("Omitted\ncancer state",2.6,15.4,6.8,FALSE,"left")
  text_mm("YAP\u2013Stem coupling",46.5,15.1,7.3,TRUE)
  text_mm("Paired\nN",n_x,15.1,6.5)
  text_mm("Program146 effect",101.25,15.1,7.3,TRUE)
  text_mm("Adjusted\nabsolute\npositive n/N",recurrence_x,14.6,6.5)
  for(i in seq_along(states)) text_mm(states[i],2.6,row_y[i],7,FALSE,"left")
  for (i in seq_len(nrow(summary_df))) {
    s <- summary_df[i,]; dom <- s$biological_domain; de <- densities[[i]]
    side <- if(s$model=="Raw") -1 else 1
    y0 <- row_y[s$state_order]
    use <- de$support
    xp <- x_position(de$x[use],dom)
    yp <- y0 + side*de$y[use]*domain_specs[[dom]]$height_mm_per_density
    col <- palette[s$model]
    grid.polygon(unit(c(xp[1],xp,tail(xp,1)),"mm"),
      unit(H-c(y0,yp,y0),"mm"),gp=gpar(fill=col,col=NA))
    # A very thin same-colour contour makes real low-density tails visible.
    grid.lines(unit(xp,"mm"),unit(H-yp,"mm"),
      gp=gpar(col=col,lwd=.25/.75,linejoin="round"))
  }
  for(dom in domains) {
    line_mm(x_position(0,dom),19.7,x_position(0,dom),71.1,
            col="#6C737970",pt=.32)
  }
  for(i in seq_len(nrow(summary_df))) {
    s <- summary_df[i,]; dom <- s$biological_domain
    side <- if(s$model=="Raw") -1 else 1
    y0 <- row_y[s$state_order]; g <- quantile_geometry[[i]]
    line_mm(x_position(s$Q25,dom),y0+side*g["offset"],
            x_position(s$Q75,dom),y0+side*g["offset"],quantile_ink,.45)
    line_mm(x_position(s$median,dom),y0+side*g["tick_low"],
            x_position(s$median,dom),y0+side*g["tick_high"],quantile_ink,.55)
  }
  for(i in seq_along(states)) {
    text_mm(as.character(coupling_n[i]),n_x,row_y[i],6.8)
    s <- summary_df[summary_df$cancer_state==states[i] & summary_df$domain==metrics[4],]
    lab <- paste0(s$positive_absolute,"/",s$N,if(states[i]=="Stress")"\u2020" else "")
    text_mm(lab,recurrence_x,row_y[i],6.8)
  }
  ticks <- list(c(-.10,-.05,0,.05,.10),c(-.10,-.05,0,.05))
  for(j in seq_along(domains)) {
    dom <- domains[j]; px <- plot_x[[dom]]
    line_mm(px[1],72.6,px[2],72.6,ink,.4)
    for(tk in ticks[[j]]) {
      xpos <- x_position(tk,dom)
      line_mm(xpos,72.6,xpos,73.25,ink,.4)
      lab <- if(tk==0) "0" else paste0(if(tk<0)"\u2212" else "+",sprintf("%.2f",abs(tk)))
      text_mm(lab,xpos,74.7,6.8)
    }
  }
  text_mm("Patient-level \u0394 Spearman \u03c1\nafter state omission",46.5,79.0,6.8)
  text_mm("Patient-level \u0394 Program146 effect\nafter state omission",101.25,79.0,6.8)
  text_mm("\u2020 After Stress omission, P57: adjusted absolute Program146 effect \u22120.013; paired \u0394 = \u22120.040.",
          2.6,84.0,6.5,FALSE,"left")
  popViewport()
  invisible(do.call(rbind,text_records))
}
grDevices::cairo_pdf(file_out(".pdf"),width=W/25.4,height=H/25.4,
                    family="Arial",pointsize=7,bg="white",onefile=TRUE)
pdf_text <- draw_panel()
dev.off()
# Cairo rounds its MediaBox down to whole points. Correct only the PDF page box;
# keep every drawing/text content byte unchanged. This is not another render.
python_bin <- Sys.getenv("S3D_QC_PYTHON", Sys.getenv("AHIPPO_PYTHON", Sys.which("python3")))
stopifnot(file.exists(python_bin))
pagebox_code <- 'import sys,os; from pypdf import PdfReader,PdfWriter; from pypdf.generic import RectangleObject; p=sys.argv[1]; r=PdfReader(p); before=r.pages[0].get_contents().get_data(); w=PdfWriter(); w.clone_document_from_reader(r); page=w.pages[0]; page.mediabox=RectangleObject([0,0,140/25.4*72,86/25.4*72]); page.cropbox=RectangleObject([0,0,140/25.4*72,86/25.4*72]); tmp=p+".pagebox.tmp"; f=open(tmp,"wb"); w.write(f); f.close(); os.replace(tmp,p); after=PdfReader(p); assert after.pages[0].get_contents().get_data()==before; print("PAGEBOX_140x86_MM_CONTENT_UNCHANGED")'
stopifnot(system2(python_bin,c("-c",shQuote(pagebox_code),shQuote(file_out(".pdf")))) == 0L)
grDevices::png(file_out("_600dpi.png"),width=W,height=H,units="mm",res=600,
               type="cairo",bg="white",pointsize=7)
png_text <- draw_panel()
dev.off()
stopifnot(identical(protected_md5,tools::md5sum(protected)))
stopifnot(min(pdf_text$size)>=6.5,min(png_text$size)>=6.5)
# Text overlap is also independently checked against final PDF glyph boxes later.
summary_table <- c(
 "| State | Coupling N | Program N | Adjusted absolute Program positive |",
 "|---|---:|---:|---:|",
 vapply(seq_along(states),function(i) {
   s <- summary_df[summary_df$cancer_state==states[i]&summary_df$domain==metrics[4],]
   paste0("| ",states[i]," | ",coupling_n[i]," | ",s$N," | ",s$positive_absolute,"/",s$N,
          if(states[i]=="Stress")"\u2020" else ""," |")
 },character(1)))
density_table <- c("| Domain | Observed min | Observed max | Linear display limits | Bandwidth | mm per density unit |",
 "|---|---:|---:|---|---:|---:|",
 vapply(domains,function(dom) {
   z <- domain_specs[[dom]]
   paste0("| ",dom," | ",fmt(z$observed[1])," | ",fmt(z$observed[2])," | ",
     fmt(z$limits[1])," to ",fmt(z$limits[2])," | ",fmt(z$bandwidth)," | ",
     fmt(z$height_mm_per_density)," |")
 },character(1)))
qc <- c(bandwidth_lock,
 "## Numeric integrity", "",
 "- 3,296/3,296 patient x omission x domain records; 44/44 distributions.",
 "- Coupling: 825 Raw + 825 Technical-adjusted. Program146: 823 + 823.",
 "- All original input columns are copied as unchanged character strings, including all decimal strings.",
 "- Read-back comparison of every original field: exact equality.",
 "- Every group N, min, Q10, Q25, median, Q75, Q90, max reproduces the audited fields (tolerance 2e-15).",
 "- Separate audited domain-specific eligibility; no common four-domain intersection, no imputation, no zero insertion.",
 "- All observed extremes enter the density calculation; no clipping/winsorization/broken/nonlinear axis.",
 "- 12/44 signs differ between median(patient deltas) and difference of matched medians; expected distinct summaries, not a QC failure.",
 paste0("- P57 after Stress: absolute adjusted Program146 = ",p57$omission_effect,
        "; paired delta = ",p57$paired_delta,"; adjusted absolute positive = 73/74."),
 "- No new rho, score, state, signature, threshold, hypothesis test, P, CI or significance annotation.",
 "", "## Locked display calculation", "",
 density_table, "",
 "- Gaussian KDE evaluated directly on >=4,096 coordinates. A four-bandwidth numerical support margin is used.",
 "- Each numerical density integrates to 1, after correcting negligible numerical-quadrature/finite-support deviations.",
 paste0("- Before normalization, integrated Gaussian mass range: ",
        fmt(min(summary_df$gaussian_area_before_normalization))," to ",
        fmt(max(summary_df$gaussian_area_before_normalization)),"."),
 "- One domain-wide density-to-mm conversion; individual peaks are NOT normalized.",
 "- Maximum half-height 2.02 mm; row pitch 4.7 mm; minimum inter-row gap 0.66 mm.",
 "- Thin same-colour density contours aid low-density tail visibility; no patient points.",
 "- Q25-Q75 dark segment and median vertical tick sit within each half; they are descriptive, not CI.",
 "- Q10/Q90 retained in data only; not additional plot lines.",
 "", "## Geometry / typography / exports", "",
 "- Exactly one layout; 140 x 86 mm, authored at true size.",
"- Vector Cairo PDF; PNG 600 dpi, white background.",
"- Cairo page-box integer-point rounding corrected to exact 140 x 86 mm using pypdf; drawing/text content bytes unchanged (no scaling or second render).",
 "- Arial throughout. Panel letter 9.7 pt; title 8.4 pt; domain titles 7.3 pt; state labels 7 pt; other core labels 6.8 pt; dense headers/footnote 6.5 pt.",
 "- No text bounding box extends past page boundaries in either R device.",
 "- Eleven state rows follow the exact requested fixed order.",
 "- One shared model legend: Raw #6EB6E4 upper, Technical-adjusted #E7837D lower. Colour is not delta direction.",
 "- Coupling paired N and Program absolute positive n/N are separate. S.G1 and HLA retain their distinct denominators.",
 "", "## Scope / provenance", "",
 paste0("- Input: ",source_file),
 "- Input SHA-256: F86773B0D6E970DA28775B58C60F31BB3877CC350299BACD7F3E1390EE1797F2",
 paste0("- ",length(protected)," protected input/archived S3D/staged S3/reference files: before/after MD5 identical."),
 "- No S3A-S3C modification, no old S3D overwrite, no A4 assembly, no S4.",
 "- Design references: Wanjie typography/alignment guidance and the Nature ridge archetype, adapted to split distributions and aligned rows rather than copied as an exact template.",
 "- User-authorized R implementation; no workbook or image-generation backend.",
 "", "## Visual review", "",
 "Numerical integrity checks do not replace inspection of the exported panel at its final display size.",
 "Final PDF/PNG inspection notes will be appended after independent rendered-output checks.",
 "")
writeLines(qc,file_out("_QC.md"),useBytes=TRUE)
readout <- c("# Supplementary Figure S3D", "",
 "## Panel title", "",
 "D  Patient-level influence of single cancer-state omission", "",
 "## Caption draft", "",
 "Distributions show the patient-level change after omission of each independently annotated cancer state relative to the same patient's full-analysis effect. Raw estimates are shown in upper blue half-violins and technical-adjusted estimates in lower coral half-violins. Left: change in within-patient YAP-Stem Spearman correlation. Right: change in the Program146 High-Other effect. The zero line indicates no change relative to the same patient's full-analysis effect. Short dark segments denote the patient-distribution IQR and vertical ticks its median, not confidence intervals. All testable patient-level changes are included; the two biological domains retain their own audited patient sets.",
 "",
 "The coupling lane reports paired N. The Program146 lane reports the number with a positive absolute technical-adjusted Program146 effect after omission, not the number with a positive change.",
 "",
 "\u2020 After Stress-state omission, P57 retained a negative adjusted Program146 effect (-0.013; paired delta = -0.040); 73/74 patients had a positive absolute effect.",
 "",
 "Patient-level responses to single-state omission were heterogeneous, while no single independently annotated cancer state accounted for the cohort-level recurrence.",
 "",
 "## Interpretation boundaries", "",
 "- This is descriptive sensitivity visualization, not equivalence, absence of an effect, causal contribution or state independence.",
 "- Median(omission - same-patient FULL) is not median(omission) - median(FULL). Twelve of 44 summary signs differ between these two descriptions; no old cohort-median displacement is plotted.",
 "- Negative delta denotes a reduction from that patient's FULL estimate; it does not itself mean a negative absolute Program146 effect.",
 "- The same patients recur in multiple omissions; 3,296 records are not 3,296 independent patients.",
 "", "## Audited denominators and absolute recurrence", "", summary_table,
 "", "## Bandwidth and scale record (locked before render)", "", density_table,
 "", "- Prespecified median of 22 group-specific bw.nrd0 values per domain; no tuning after render.",
 "- All groups have unit numerical density area; common density-height conversion within each domain.",
 "- Numerical Gaussian support is finite; all real patient values are included, with four bandwidths of additional margin.",
 "", "## Visual provenance", "",
 "This display adapts the split-distribution and aligned-ridgeline principles recorded in the design references. It is not an exact reproduction of a published layout.",
 "Design references: Forde et al., Nature Medicine 2021, Fig. 5b,c (doi:10.1038/s41591-021-01541-0); Proctor et al., Nature Medicine 2021, Extended Data Fig. 3 (doi:10.1038/s41591-021-01383-w). Their biological group meanings and inference were not transferred.",
 "",
 "## Export", "", "True-size display; numerical checks do not replace visual inspection.",
 "Standalone 140 x 86 mm panel; figure-page assembly is separate.")
writeLines(readout,file_out("_readout.md"),useBytes=TRUE)
cat("Numerical checks completed; standalone panel rendered.\n")
print(summary_df[,c("cancer_state","domain","N","minimum","Q25","median","Q75","maximum","bandwidth_used")])
