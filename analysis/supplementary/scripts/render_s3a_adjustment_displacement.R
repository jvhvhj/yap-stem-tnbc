# Purpose: Render S3A adjustment displacement
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Supplementary Fig. S1–S4.
# Standalone S3A display from existing results; no model fitting or new tests.
# Run with Rscript from any working directory; PDF is native vector graphics.
options(stringsAsFactors=FALSE, warn=1)
if (.Platform$OS.type=='windows') invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
suppressPackageStartupMessages(library(grid))
args <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '', args[grepl('^--file=',args)][1])
out <- dirname(normalizePath(script, winslash='/', mustWork=TRUE))
root <- dirname(out)
src <- file.path(root,'Supplementary_S3ABC_redesign_gate','S3A_data_topology_audit.tsv')
plot_src <- file.path(out,'Supplementary_Figure_S3A_integrated_plotdata.tsv')
d <- read.delim(plot_src,check.names=FALSE,na.strings='NA')
a <- read.delim(src,check.names=FALSE)
stopifnot(nrow(d)==78L,!anyDuplicated(d$patient_id),identical(d$patient_id,a$patient_id))
for(nm in c('technical_rho','extended_rho','delta_rho','delta_q25','delta_median','delta_q75')) {
  stopifnot(identical(d[[nm]],a[[nm]]),all(is.finite(d[[nm]])))
}
stopifnot(sum(d$display_category=='Positive in both')==68L,
          sum(d$display_category=='Negative in both')==6L,
          sum(d$display_category=='Positive-to-negative')==4L,
          setequal(d$patient_id[d$display_category=='Positive-to-negative'],c('P15','P23','P67','P69')),
          sum(d$sign_retained)==74L,
          max(abs(d$technical_rho+d$delta_rho-d$extended_rho))<1e-14)
W <- 130; H <- 92
main <- c(left=18,right=105,top=16,bottom=78)
marg <- c(left=106.5,right=125,top=16,bottom=78)
xlim <- c(-.28,.54); ylim <- c(-.29,.05)
stopifnot(all(d$technical_rho>xlim[1]&d$technical_rho<xlim[2]),
          all(d$delta_rho>ylim[1]&d$delta_rho<ylim[2]))
mx <- function(x) main['left']+(x-xlim[1])/diff(xlim)*diff(main[c('left','right')])
my <- function(y) main['bottom']-(y-ylim[1])/diff(ylim)*diff(main[c('top','bottom')])
# Frozen descriptive summaries, not re-estimated.
q1 <- unique(d$delta_q25); med <- unique(d$delta_median); q3 <- unique(d$delta_q75)
stopifnot(length(q1)==1,length(med)==1,length(q3)==1)
# The sole authorized smoothing is the marginal density; exactly one fixed rule.
den <- density(d$delta_rho,bw='nrd0')
# Interpolate the density *rendering grid* to the shared display extent only;
# this does not interpolate patient values. All 78 original deltas enter density().
dy <- sort(unique(c(ylim[1],den$x[den$x>ylim[1]&den$x<ylim[2]],ylim[2])))
dv <- approx(den$x,den$y,xout=dy,rule=2)$y
ink <- '#37434B'; neutral <- '#8A8F94'; light <- '#D9DDE1'
blue <- '#6EB6E4'; coral <- '#E7837D'
mm <- function(x) unit(x,'mm')
gx <- function(x) mm(x)
gy <- function(y) mm(H-y)
line <- function(x0,y0,x1,y1,col=ink,lwd=.55,lty='solid',name=NULL) {
  grid.segments(gx(x0),gy(y0),gx(x1),gy(y1),gp=gpar(col=col,lwd=lwd,lty=lty,lineend='butt'),name=name)
}
txt <- function(s,x,y,size=7,face='plain',col=ink,just='centre',rot=0,name=NULL) {
  grid.text(s,x=gx(x),y=gy(y),just=just,rot=rot,gp=gpar(fontfamily='Arial',fontsize=size,fontface=face,col=col),name=name)
}
circle <- function(x,y,r=.44,col=neutral,fill=col,lwd=.35,name=NULL) {
  grid.circle(gx(x),gy(y),r=mm(r),gp=gpar(col=col,fill=fill,lwd=lwd),name=name)
}
diamond <- function(x,y,r=.63,col=coral,name=NULL) {
  grid.polygon(gx(x+c(-r,0,r,0)),gy(y+c(0,-r,0,r)),
               gp=gpar(fill=adjustcolor(col,alpha.f=.6),col=col,lwd=.65),name=name)
}

draw_panel <- function() {
  grid.newpage()
  txt('A',4.5,5.5,11,'bold',just='left',name='panel_letter')
  txt('Additional-adjustment sensitivity',10,5.5,9,'bold',just='left',name='panel_title')
  # A single compact legend; minority classes only carry chromatic emphasis.
  circle(18.7,11.5,.44,neutral,adjustcolor(neutral,alpha.f=.80),name='legend_positive')
  txt('Positive in both',20.4,11.5,6.8,just='left')
  circle(53,11.5,.47,blue,'white',.9,name='legend_negative')
  txt('Negative in both',54.7,11.5,6.8,just='left')
  diamond(90,11.5,.60,coral,name='legend_switch')
  txt('Positive-to-negative',91.7,11.5,6.8,just='left')

  # Mathematical references drawn behind all patient observations.
  line(main['left'],my(0),marg['right'],my(0),col=light,lwd=.6,lty='22',name='delta_zero')
  line(mx(0),main['top'],mx(0),main['bottom'],col=light,lwd=.6,lty='22',name='technical_zero')
  line(mx(-.05),my(.05),mx(.29),my(-.29),col=neutral,lwd=.65,lty='43',name='extended_zero')
  # The label sits beside a sparsely occupied lower part of the exact boundary.
  angle <- -atan((62/.34)/(87/.82))*180/pi
  txt('Extended \u03c1 = 0',mx(.231)-1.8,my(-.231)-.8,6.7,rot=angle,name='extended_boundary_label')

  # Retain coordinates exactly. Transparent fills allow natural near-overlaps.
  for(i in which(d$display_category=='Positive in both')) {
    circle(mx(d$technical_rho[i]),my(d$delta_rho[i]),r=.44,
           col=adjustcolor(neutral,alpha.f=.9),fill=adjustcolor(neutral,alpha.f=.65),lwd=.3,
           name=paste0('patient_',d$patient_id[i]))
  }
  for(i in which(d$display_category=='Negative in both')) {
    circle(mx(d$technical_rho[i]),my(d$delta_rho[i]),r=.48,col=blue,fill='white',lwd=.9,
           name=paste0('patient_',d$patient_id[i]))
  }
  for(i in which(d$display_category=='Positive-to-negative')) {
    diamond(mx(d$technical_rho[i]),my(d$delta_rho[i]),r=.60,col=coral,
            name=paste0('patient_',d$patient_id[i]))
  }
  # Label displacement only, never patient-coordinate displacement.
  label_pos <- data.frame(patient_id=c('P67','P69','P15','P23'),
                          x=c(46,45,50,52),y=c(31.5,40.5,47,55))
  for(i in seq_len(nrow(label_pos))) {
    j <- match(label_pos$patient_id[i],d$patient_id)
    px <- mx(d$technical_rho[j]); py <- my(d$delta_rho[j])
    ex <- label_pos$x[i]+2.7; ey <- label_pos$y[i]
    dr <- sqrt((ex-px)^2+(ey-py)^2)
    line(px+(ex-px)*.8/dr,py+(ey-py)*.8/dr,ex,ey,col=neutral,lwd=.42)
    txt(label_pos$patient_id[i],label_pos$x[i],label_pos$y[i],7,col=ink,
        name=paste0('label_',label_pos$patient_id[i]))
  }
  txt('74/78 retained sign',main['right']-1,main['top']+3,7,'bold',just='right',name='headline')

  # Attached one-sided marginal, identical physical/data y mapping.
  base <- marg['left']; density_width <- 17.6
  dx <- base + dv/max(dv)*density_width
  grid.polygon(gx(c(base,dx,base)),gy(c(my(dy[1]),my(dy),my(tail(dy,1)))),
               gp=gpar(fill=adjustcolor(light,alpha.f=.8),col=NA),name='delta_density_fill')
  grid.lines(gx(dx),gy(my(dy)),gp=gpar(col=neutral,lwd=.65),name='delta_density_outline')
  line(base,my(ylim[1]),base,my(ylim[2]),col=light,lwd=.35,name='density_baseline')
  line(base+.5,my(q1),base+.5,my(q3),col=ink,lwd=2.1,name='delta_iqr')
  line(base-.4,my(med),base+1.6,my(med),col=ink,lwd=1.1,name='delta_median')

  # Only one y axis and one x axis; no enclosing rectangle or plot grid.
  line(main['left'],main['top'],main['left'],main['bottom'],lwd=.55,name='main_y_axis')
  line(main['left'],main['bottom'],main['right'],main['bottom'],lwd=.55,name='main_x_axis')
  for(x in seq(-.2,.5,.1)) {
    line(mx(x),main['bottom'],mx(x),main['bottom']+1,lwd=.5)
    txt(if(abs(x)<1e-10)'0' else sub('^-','\u2212',sprintf('%.1f',x)),mx(x),main['bottom']+3.2,6.7)
  }
  for(y in c(-.25,-.20,-.15,-.10,-.05,0,.05)) {
    line(main['left']-1,my(y),main['left'],my(y),lwd=.5)
    txt(if(y==0)'0' else sub('^-','\u2212',sprintf('%.2f',y)),main['left']-1.7,my(y),6.7,just='right')
  }
  txt('Technical-adjusted Spearman \u03c1',(main['left']+main['right'])/2,87,7.4,name='x_title')
  txt('\u0394\u03c1 (Extended \u2212 Technical)',5,(main['top']+main['bottom'])/2,7.4,rot=90,name='y_title')
}

pdf_file <- file.path(out,'Supplementary_Figure_S3A_integrated.pdf')
stopifnot(!file.exists(pdf_file))
cairo_pdf(pdf_file,width=W/25.4,height=H/25.4,family='Arial',pointsize=7,bg='white',fallback_resolution=600)
draw_panel()
all_names <- grid.ls(print=FALSE)$name
stopifnot(sum(grepl('^patient_',all_names))==78L,
          setequal(sub('^label_','',all_names[grepl('^label_',all_names)]),c('P15','P23','P67','P69')))
dev.off()
# Cairo rounds page bounds to integer PostScript points. Restore the requested
# physical canvas only, without scaling/moving a single vector/data coordinate.
python <- Sys.getenv("AHIPPO_PYTHON", Sys.which("python3"))
poppler <- Sys.getenv("AHIPPO_PDFTOPPM", Sys.which("pdftoppm"))
stopifnot(file.exists(python),file.exists(poppler))
normalise_page <- paste(
  'import io,sys; from pathlib import Path; from pypdf import PdfReader,PdfWriter; from pypdf.generic import RectangleObject;',
  'p=Path(sys.argv[1]); w=PdfWriter(clone_from=PdfReader(io.BytesIO(p.read_bytes())));',
  'w.pages[0].mediabox=RectangleObject([0,0,130*72/25.4,92*72/25.4]);',
  'w.pages[0].cropbox=RectangleObject([0,0,130*72/25.4,92*72/25.4]);',
  'w.write(str(p))'
)
stopifnot(system2(python,c('-c',shQuote(normalise_page),shQuote(pdf_file)))==0L)
stopifnot(system2(poppler,c('-png','-singlefile','-r','600',shQuote(pdf_file),
                  shQuote(file.path(out,'Supplementary_Figure_S3A_integrated_600dpi'))))==0L)
cat('patients=78; numerical_changes=0; categories=68/6/4; sign_retained=74/78\n')
cat('delta_density_bw=',format(den$bw,digits=17),'; method=nrd0; n=',den$n,'\n',sep='')
cat('density_kernel_extent=',paste(format(range(den$x),digits=17),collapse=','),'\n',sep='')
cat('IQR=',paste(format(c(q1,q3),digits=17),collapse=','),'; median=',format(med,digits=17),'\n',sep='')
cat('PDF=',pdf_file,'\n',sep='')
cat('Standalone adjustment-displacement panel rendered.\n')
