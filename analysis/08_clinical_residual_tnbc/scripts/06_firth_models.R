# Portable reproduction of documented clinical models; not executed during release editing.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2L) stop('Usage: Rscript 06_firth_models.R merged.tsv output_dir')
if (getRversion()!=numeric_version('4.3.3')) stop('Documented model environment requires R 4.3.3')
suppressPackageStartupMessages(library(logistf))
if (packageVersion('logistf')!=numeric_version('1.26.1')) stop('Requires logistf 1.26.1')
d <- read.delim(args[1], check.names=FALSE)
stopifnot(nrow(d)==48L, !anyDuplicated(d$Public_ID), sum(d$Outcome_group=='D')==34L)
d$D_event <- as.integer(d$Outcome_group=='D')
d$RCB_class <- factor(d$RCB_class, levels=c(2,3))
keep <- d$Operative_subtype=='HR-HER2-';stopifnot(sum(keep)==45L)
out <- args[2];dir.create(out,recursive=TRUE,showWarnings=FALSE)
export <- function(fit,label,n) {
 data.frame(Model=label,Patients_n=n,Term=names(coef(fit)),Beta=unname(coef(fit)),
 OR=exp(unname(coef(fit))),CI95_lower=exp(unname(fit$ci.lower)),CI95_upper=exp(unname(fit$ci.upper)),
 Penalized_likelihood_ratio_P=unname(fit$prob),Interval_method='95% profile penalized-likelihood')
}
primary <- logistf(D_event ~ RCB_class + Program146_1SD,data=d,pl=TRUE)
unadjusted <- logistf(D_event ~ Program146_1SD,data=d,pl=TRUE)
sensitivity <- logistf(D_event ~ RCB_class + Program146_1SD,data=d[keep,],pl=TRUE)
res <- rbind(export(primary,'RCB-adjusted primary',48),export(unadjusted,'Unadjusted',48),export(sensitivity,'Operative HR-negative HER2-negative sensitivity',45))
write.table(res,file.path(out,'clinical_model_results.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
saveRDS(list(data=d,primary=primary,unadjusted=unadjusted,sensitivity=sensitivity),file.path(out,'clinical_models.rds'))
