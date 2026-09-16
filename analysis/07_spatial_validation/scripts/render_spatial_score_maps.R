#!/usr/bin/env Rscript
# Purpose: Figure 6A left spatial map display reset
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 6 and Supplementary Fig. S4.
# Display-only replacement of the two P1 / 117B spatial-score blocks.
# No score estimation, correlation, patient statistic, or composite rendering.
suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
args <- commandArgs(FALSE)
script <- normalizePath(sub("^--file=", "", grep("^--file=", args, value=TRUE)[1]), winslash="/")
out <- dirname(script)
root <- dirname(out)
.libPaths(c(file.path(root, "0703_rebuild/R_library"), .libPaths()))
stopifnot(all(vapply(c("ragg", "png", "digest"), requireNamespace, logical(1), quietly=TRUE)))
library(grid)
options(digits=17, scipen=999)

input <- file.path(root, "Figure6A_final_reset/Figure6A_final_reset_source.tsv")
upstream <- file.path(root, "Figure6A_R1_premium_candidate/Figure6A_R1_spatial_plot_source.tsv")
old_script <- file.path(root, ".tmp/Figure6A_final_reset/render_final_reset.py")
other_layout <- file.path(root, "Figure6_ABC_lowfidelity/Figure6_ABC_layout_source.tsv")
read_chars <- function(p) read.delim(p, colClasses="character", check.names=FALSE,
                                    quote="", na.strings=character(), fileEncoding="UTF-8")
hash <- function(p) digest::digest(file=p, algo="sha256", serialize=FALSE)

# Snapshot the unchanged complete Figure 6A and its source/split-layout chain.
protected_dirs <- file.path(root,c("Figure6A_final_reset", "Figure6A_R1_premium_candidate",
                                   "Figure6_ABC_lowfidelity"))
protected <- unique(c(unlist(lapply(protected_dirs, list.files, recursive=TRUE, full.names=TRUE)), old_script))
protected <- protected[file.exists(protected) & !dir.exists(protected)]
protected <- protected[grepl("\\.(pdf|png|svg|tsv|md|R|py)$", protected, ignore.case=TRUE) &
                         file.info(protected)$size < 20*1024^2]
before_hash <- vapply(protected, hash, character(1))
cat("Verified source-chain snapshot:",length(protected),"files\n")

src <- read_chars(input)
s <- src[src$record_type == "REPRESENTATIVE_SPOT", , drop=FALSE]
u <- read_chars(upstream)
stopifnot(nrow(s)==1109L, !anyDuplicated(s$barcode),
          identical(unique(s$patient_id), "P1"), identical(unique(s$section_id), "117B"),
          identical(unique(s$sample_id), "GSM6433597"))
u <- u[match(s$barcode, u$barcode), , drop=FALSE]
pairs <- c(hires_x="hires_x", hires_y="hires_y", joint_spatial="Joint_spatial",
           program146_spatial="Program146_spatial")
for (nm in names(pairs)) {
  stopifnot(identical(as.numeric(s[[nm]]), as.numeric(u[[pairs[[nm]]]])))
}
if (file.exists(other_layout)) {
  abc <- read_chars(other_layout)
  abc <- abc[abc$record_type == "REPRESENTATIVE_SPOT", , drop=FALSE]
  abc <- abc[match(s$barcode, abc$barcode), , drop=FALSE]
  stopifnot(nrow(abc)==1109L, identical(s$section_id, abc$section_id))
  for (nm in names(pairs)) stopifnot(identical(as.numeric(s[[nm]]), as.numeric(abc[[nm]])))
}
stopifnot(hash(upstream) == src$source_sha256[src$record_type=="SOURCE_REGISTRY" &
            basename(src$source_file)==basename(upstream)])

# Display reference: spatial-feature/plot.R (ADAPT); fixed physical coordinates,
# circular spot marks and raster spot layer. No synthetic data or template scales.
# Typography and endpoints come from the current Figure 6A render script.
# Its unavailable external nature_theme.R dependency is not used or reconstructed.
old_code <- readLines(old_script, warn=FALSE, encoding="UTF-8")
get_token <- function(nm) {
  ln <- grep(paste0("^", nm, " = "), old_code, value=TRUE)
  stopifnot(length(ln)==1L)
  sub('.*"(#[A-Fa-f0-9]{6})".*', '\\1', ln)
}
blue <- get_token("BLUE")
red <- get_token("CORAL")
dark <- get_token("DARK")
neutral <- "#FFFFFF" # Explicit white midpoint requested for these replacements.
palette <- grDevices::colorRamp(c(blue, neutral, red), space="rgb")
colour_at <- function(v, L) {
  clipped <- pmax(-L, pmin(L, v))
  rgb <- palette((clipped/L+1)/2)
  grDevices::rgb(rgb[,1], rgb[,2], rgb[,3], maxColorValue=255)
}

probs <- c(0, .01, .025, .05, .5, .95, .975, .99, 1)
q_names <- c("min", "q01", "q025", "q05", "median", "q95", "q975", "q99", "max")
columns <- c("joint_spatial", "program146_spatial")
titles <- c("Joint YAP\u2013Stem score", "Program146 score")
stems <- c("Figure6A_leftmap_JointYAPStem_RESET", "Figure6A_leftmap_Program146_RESET")
audit <- do.call(rbind, lapply(seq_along(columns), function(i) {
  v <- as.numeric(s[[columns[i]]])
  stopifnot(length(v)==1109L, all(is.finite(v)))
  q <- unname(quantile(v, probs=probs, type=7))
  L <- max(abs(q[c(3,7)]))
  stopifnot(is.finite(L), L>0)
  cbind(data.frame(score=columns[i], patient_id="P1", section_id="117B", n_spots=length(v)),
        as.data.frame(as.list(setNames(q,q_names))),
        data.frame(display_min=-L, display_mid=0, display_max=L,
                   below_display_min_n=sum(v< -L), above_display_max_n=sum(v>L),
                   clipped_display_n=sum(abs(v)>L), clipped_display_percent=100*mean(abs(v)>L),
                   quantile_method="R type 7", raw_values_changed=FALSE))
}))

# Preserve original numerical strings byte-for-byte in source columns.
render_source <- s[,c("patient_id","sample_id","section_id","barcode",
                      "hires_x","hires_y","joint_spatial","program146_spatial")]
render_source$source_file <- basename(input)
render_source$source_sha256 <- hash(input)
for (i in seq_along(columns)) {
  nm <- columns[i]
  v <- as.numeric(s[[nm]])
  L <- audit$display_max[i]
  render_source[[paste0(nm,"_display")]] <- sprintf("%.17g", pmax(-L,pmin(L,v)))
  render_source[[paste0(nm,"_display_clipped")]] <- abs(v)>L
  render_source[[paste0(nm,"_display_limit")]] <- sprintf("%.17g",L)
}
write.table(render_source, file.path(out,"Figure6A_left_spatial_maps_render_source.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE, fileEncoding="UTF-8")
write.table(audit, file.path(out,"Figure6A_left_spatial_maps_display_audit.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE, fileEncoding="UTF-8")
check <- read_chars(file.path(out,"Figure6A_left_spatial_maps_render_source.tsv"))
for(nm in names(pairs)) stopifnot(identical(check[[nm]], s[[nm]]))

x <- as.numeric(s$hires_x)
y <- as.numeric(s$hires_y)
pad_x <- max(diff(range(x))*.08,18)
pad_y <- max(diff(range(y))*.08,18)
xlim <- c(min(x)-pad_x,max(x)+pad_x)
ylim <- c(min(y)-pad_y,max(y)+pad_y)
width_mm <- 60
height_mm <- 65
map_width_mm <- 54
map_height_mm <- map_width_mm * diff(ylim)/diff(xlim)
map_top_mm <- 58
map_bottom_mm <- map_top_mm-map_height_mm

# Matplotlib s=5.8 means squared-point diameter parameter. Preserve semantics,
# increase marker area by 20%, and use the enlarged independent block to maintain gaps.
marker_line <- grep("x, y, c=values, s=", old_code, value=TRUE)
stopifnot(length(marker_line)==1L)
old_area_pt2 <- as.numeric(sub(".*s=([0-9.]+),.*", "\\1", marker_line))
new_area_pt2 <- old_area_pt2*1.20
radius_mm <- sqrt(new_area_pt2)/2*25.4/72
pixel_dpi <- 600

make_spot_raster <- function(v,L) {
  tmp <- tempfile(fileext=".png")
  on.exit(unlink(tmp),add=TRUE)
  # Integer pixel size ceil ensures embedded resolution is at least 600 dpi.
  ragg::agg_png(tmp, width=ceiling(map_width_mm/25.4*pixel_dpi),
                height=ceiling(map_height_mm/25.4*pixel_dpi), units="px", res=pixel_dpi,
                background="white", scaling=1)
  grid.newpage()
  grid.circle(x=unit((x-xlim[1])/diff(xlim),"npc"),
              y=unit((ylim[2]-y)/diff(ylim),"npc"), r=unit(radius_mm,"mm"),
              gp=gpar(fill=colour_at(v,L),col=NA))
  dev.off()
  png::readPNG(tmp,native=TRUE)
}

draw_block <- function(i, raster) {
  grid.newpage()
  grid.rect(gp=gpar(fill="white",col=NA))
  grid.text(titles[i],x=unit(width_mm/2,"mm"),y=unit(62,"mm"),
            gp=gpar(fontfamily="Arial",fontsize=9,fontface="bold",col=dark))
  grid.raster(raster,x=unit(width_mm/2,"mm"),
              y=unit(map_bottom_mm+map_height_mm/2,"mm"),
              width=unit(map_width_mm,"mm"),height=unit(map_height_mm,"mm"),interpolate=FALSE)
  L <- audit$display_max[i]
  bar_x <- 12
  bar_w <- 36
  n_rect <- 512L
  grid.rect(x=unit(bar_x+(seq_len(n_rect)-.5)*bar_w/n_rect,"mm"),y=unit(4.9,"mm"),
            width=unit(bar_w/n_rect+.001,"mm"),height=unit(1.65,"mm"),
            gp=gpar(fill=colour_at(seq(-L,L,length.out=n_rect),L),col=NA))
  tick_x <- bar_x+c(0,.5,1)*bar_w
  grid.segments(x0=unit(tick_x,"mm"),x1=unit(tick_x,"mm"),
                y0=unit(4.075,"mm"),y1=unit(3.5,"mm"),gp=gpar(col=dark,lwd=.5))
  grid.text(c(sprintf("%.3f",-L),"0",sprintf("%.3f",L)),
            x=unit(tick_x,"mm"),y=unit(2,"mm"),
            gp=gpar(fontfamily="Arial",fontsize=7,col=dark))
}

for (i in seq_along(columns)) {
  v <- as.numeric(s[[columns[i]]])
  raster <- make_spot_raster(v,audit$display_max[i])
  cairo_pdf(file.path(out,paste0(stems[i],".pdf")),width=width_mm/25.4,
            height=height_mm/25.4,family="Arial",bg="white",fallback_resolution=600)
  draw_block(i,raster)
  dev.off()
  ragg::agg_png(file.path(out,paste0(stems[i],"_600dpi.png")),width=width_mm,
                height=height_mm,units="mm",res=600,background="white")
  draw_block(i,raster)
  dev.off()
}

after_hash <- vapply(protected,hash,character(1))
stopifnot(identical(before_hash,after_hash))
qc <- c("# Figure 6A spatial-map QC", "", "## Numerical and source checks",
        "- Representative: P1 / 117B / GSM6433597, unchanged.",
        "- 1,109/1,109 spots per map; zero missing, zero duplicated barcodes.",
        "- Both maps use identical coordinates, orientation and 8% coordinate padding.",
        "- Scores and coordinates match the current Figure6A source, upstream spatial table and later A-C layout table exactly.",
        "- Raw source numeric strings survive TSV export unchanged; display-clipped values are separate columns.",
        "- No scores, correlations, P values, patient statistics or biological classifications were computed or modified.",
        sprintf("- SHA-256 unchanged for %d pre-existing Figure 6 figure/source/script/report files checked (each <20 MiB, including the original render script).",length(protected)),
        "", "## Rendering specification",
        sprintf("- Two independent blocks only: %g x %g mm each; no composite created.",width_mm,height_mm),
        "- Arial; titles 9 pt and colour-bar ticks 7 pt; white background.",
        sprintf("- Current point area %.3f pt^2 -> %.3f pt^2 (+20%%); diameter %.4f mm.",old_area_pt2,new_area_pt2,2*radius_mm),
        "- Circular spot glyphs; no borders, opacity reduction, H&E, interpolation, smoothing or contouring.",
        "- Spot layer rasterized at >=600 dpi in PDF; titles and colour-bar ticks remain vector text.",
        "- PNG export 600 dpi. Separate symmetric colour limits; explicit zero midpoint.",
        "- Outliers saturate only in the colour encoding; no raw-score winsorization.",
        "", "## Reference adaptation",
        "- ADAPT: nature-figure-archetypes spatial-feature continuous-feature layer; fixed aspect and rasterized spatial circles retained.",
        "- Template synthetic data, sequential magma palette, categorical domain panel and captions are not used.",
        "- Wanjie review principles used for typography, white background and no overlap; current Figure 6A endpoint colours retained.",
        "- Missing external nature_theme.R not reconstructed; R grid/ragg uses the already approved Figure 6A parameters.",
        "", "## File identity",
        paste0("- Input: ",basename(input),"; SHA-256 ",hash(input)),
        paste0("- Upstream: ",basename(upstream),"; SHA-256 ",hash(upstream)),
        paste0("- Prior render: ",basename(old_script),"; SHA-256 ",hash(old_script)),
        "", "## Visual verification", "- Pending inspection of final PDF rasterizations and 600 dpi PNG files.")
writeLines(qc,file.path(out,"Figure6A_left_spatial_maps_QC.md"),useBytes=TRUE)
readout <- c("# Figure 6A spatial-score maps", "",
             "The two maps show the same 1,109 Visium spots from patient P1, section 117B (GSE210616; GSM6433597).",
             "Each circle represents one mixed-cell tissue measurement at its recorded spatial coordinate.","",
             "Blue indicates negative scores, white indicates zero, and red indicates positive scores. Colour is score magnitude, not statistical significance.",
             "The scores use separate units and colour limits; equal colours must not be interpreted as equal scores between maps.","",
             "## Display limits",
             "Each limit uses L = max(abs(Q2.5), abs(Q97.5)), with R type-7 quantiles over the displayed section. The scale is [-L, 0, +L].",
             sprintf("- Joint YAP-Stem score: [%.9f, 0, +%.9f]; %d spots below and %d above the display limits (%d/%d; %.3f%% total saturated).",
                     audit$display_min[1],audit$display_max[1],audit$below_display_min_n[1],audit$above_display_max_n[1],
                     audit$clipped_display_n[1],audit$n_spots[1],audit$clipped_display_percent[1]),
             sprintf("- Program146 score: [%.9f, 0, +%.9f]; %d spots below and %d above the display limits (%d/%d; %.3f%% total saturated).",
                     audit$display_min[2],audit$display_max[2],audit$below_display_min_n[2],audit$above_display_max_n[2],
                     audit$clipped_display_n[2],audit$n_spots[2],audit$clipped_display_percent[2]),
             "Extreme colours saturate at these limits for display only. All original values, including values outside these limits, remain in the accompanying source table.","",
             "## Suggested caption wording",
             "Computational Joint YAP-Stem and Program146 spot-level scores are shown at the corresponding Visium coordinates for section 117B. Separate zero-centred colour limits are based on the 2.5th and 97.5th percentiles; extreme values are colour-saturated for display only. Visium spots are mixed-cell tissue measurements.","",
             "Only the two map blocks are replaced. No surrounding figure content or statistical result has changed.")
writeLines(readout,file.path(out,"Figure6A_left_spatial_maps_readout.md"),useBytes=TRUE)
print(audit,row.names=FALSE)
cat("PROTECTED_UNCHANGED",length(protected),"\n")
cat("SPOT_DIAMETER_MM",2*radius_mm,"MAP_SIZE_MM",map_width_mm,map_height_mm,"\n")
