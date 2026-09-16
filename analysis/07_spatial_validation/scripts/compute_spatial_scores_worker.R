# Purpose: Compute section-level YAP-Stem and Program146 spot scores
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 6 and Supplementary Fig. S4.
options(stringsAsFactors = FALSE)
.libPaths(c(Sys.getenv("R_LIBS_USER"), .Library))

suppressPackageStartupMessages({
  library(Matrix)
  library(MASS)
  library(future)
  library(future.apply)
  library(ggplot2)
  library(reshape2)
  library(gridExtra)
  library(matrixStats)
  library(Rcpp)
  library(RcppArmadillo)
  library(hdf5r)
})

sct_root <- ".runtime/Figure6_S0Q3_sctransform_src"
dyn.load(file.path(sct_root, "src", "sctransform.dll"))
for (f in sort(list.files(file.path(sct_root, "R"), full.names = TRUE))) {
  source(f, local = .GlobalEnv)
}

worker_id <- as.integer(commandArgs(trailingOnly = TRUE)[1])
n_workers <- as.integer(commandArgs(trailingOnly = TRUE)[2])
stopifnot(is.finite(worker_id), is.finite(n_workers), worker_id >= 1L, worker_id <= n_workers)

input_root <- ".runtime/Figure6_S1_inputs"
output_root <- ".runtime/Figure6_S1_chunks"
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

registry <- read.delim(file.path(input_root, "section_registry.tsv"), check.names = FALSE)
score_genes <- read.delim(file.path(input_root, "score_gene_list.tsv"), check.names = FALSE)
program146_original <- scan(file.path(input_root, "146_gene_program.txt"), what = character(), quiet = TRUE)
alias_map <- c(CYR61 = "CCN1", CTGF = "CCN2", `HNRNPU-AS1` = "HNRNPU")
map_alias <- function(x) {
  y <- x
  hit <- x %in% names(alias_map)
  y[hit] <- unname(alias_map[x[hit]])
  y
}
yap17_original <- score_genes$gene[score_genes$in_YAP_score]
stem21_original <- score_genes$gene[score_genes$in_Stemness_score]
yap17 <- unique(map_alias(yap17_original))
stem21 <- unique(map_alias(stem21_original))
program146 <- unique(map_alias(program146_original))
program142 <- setdiff(program146, "HNRNPU")

stopifnot(length(yap17_original) == 17L, length(stem21_original) == 21L,
          length(program146_original) == 146L)
stopifnot(length(intersect(union(yap17, stem21), program146)) == 0L)

read_10x_h5 <- function(path) {
  infile <- hdf5r::H5File$new(filename = path, mode = "r")
  on.exit(infile$close_all(), add = TRUE)
  g <- "matrix"
  features <- make.unique(infile[[paste0(g, "/features/name")]][])
  barcodes <- infile[[paste0(g, "/barcodes")]][]
  mat <- Matrix::sparseMatrix(
    i = infile[[paste0(g, "/indices")]][] + 1L,
    p = infile[[paste0(g, "/indptr")]][],
    x = as.numeric(infile[[paste0(g, "/data")]][]),
    dims = infile[[paste0(g, "/shape")]][], repr = "T")
  rownames(mat) <- features
  colnames(mat) <- barcodes
  as(mat, "CsparseMatrix")
}

z_vector <- function(x) {
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x)) / s
}

program_score <- function(expr, requested) {
  present <- intersect(requested, rownames(expr))
  missing <- setdiff(requested, present)
  if (!length(present)) stop("No requested program genes in SCT expression")
  x <- as.matrix(expr[present, , drop = FALSE])
  sds <- matrixStats::rowSds(x)
  nonzero <- is.finite(sds) & sds > 0
  z <- matrix(0, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
  if (any(nonzero)) {
    z[nonzero, ] <- sweep(sweep(x[nonzero, , drop = FALSE], 1, rowMeans(x[nonzero, , drop = FALSE]), "-"), 1, sds[nonzero], "/")
  }
  list(score = colMeans(z), present = present, missing = missing,
       nonzero = present[nonzero], zero_sd = present[!nonzero])
}

selected <- which(((seq_len(nrow(registry)) - 1L) %% n_workers) == (worker_id - 1L))
for (ri in selected) {
  rec <- registry[ri, , drop = FALSE]
  sample_id <- rec$sample_id
  section_id <- rec$section_id
  chunk_path <- file.path(output_root, paste0(sample_id, "_", section_id, "_scores.tsv"))
  meta_path <- file.path(output_root, paste0(sample_id, "_", section_id, "_score_meta.tsv"))
  err_path <- file.path(output_root, paste0(sample_id, "_", section_id, "_ERROR.txt"))
  if (file.exists(chunk_path) && file.exists(meta_path)) {
    message(sprintf("WORKER %d SKIP %s %s", worker_id, sample_id, section_id))
    next
  }
  message(sprintf("WORKER %d START %d/%d %s %s", worker_id, ri, nrow(registry), sample_id, section_id))
  started <- Sys.time()
  tryCatch({
    h5_dir <- file.path(".runtime/Figure6_S0Q2_geo", sample_id)
    h5_files <- list.files(h5_dir, pattern = "filtered_feature_bc_matrix[.]h5$", full.names = TRUE)
    if (length(h5_files) != 1L) stop(sprintf("Expected one H5, found %d", length(h5_files)))
    umi <- read_10x_h5(h5_files)
    if (ncol(umi) != rec$n_expression_spots) stop("H5 spot count disagrees with registry")
    n_count <- Matrix::colSums(umi)
    n_feature <- Matrix::colSums(umi > 0)

    set.seed(1448145)
    fit <- vst(
      umi = umi,
      cell_attr = data.frame(row.names = colnames(umi)),
      n_cells = min(5000L, ncol(umi)),
      residual_type = "none",
      return_cell_attr = TRUE,
      return_gene_attr = TRUE,
      return_corrected_umi = FALSE,
      verbosity = 0)
    corrected <- as(correct_counts(fit, umi, verbosity = 0), "CsparseMatrix")
    expr <- log1p(corrected)

    yap_present <- intersect(yap17, rownames(expr))
    stem_present <- intersect(stem21, rownames(expr))
    if (!length(yap_present) || !length(stem_present)) stop("No measurable YAP17 or Stem21 genes")
    yap_mean <- Matrix::colMeans(expr[yap_present, , drop = FALSE])
    stem_mean <- Matrix::colMeans(expr[stem_present, , drop = FALSE])
    joint <- z_vector(yap_mean) + z_vector(stem_mean)

    p146 <- program_score(expr, program146)
    p142 <- program_score(expr, program142)
    tab <- data.frame(
      patient_id = rec$patient_id,
      sample_id = sample_id,
      section_id = section_id,
      barcode = colnames(umi),
      nCount_raw = as.numeric(n_count),
      nFeature_raw = as.integer(n_feature),
      YAP17_mean_normalized = as.numeric(yap_mean),
      Stem21_mean_normalized = as.numeric(stem_mean),
      Joint_spatial = as.numeric(joint),
      Program146_spatial = as.numeric(p146$score),
      Program142_spatial = as.numeric(p142$score),
      stringsAsFactors = FALSE)
    meta <- data.frame(
      patient_id = rec$patient_id, sample_id = sample_id, section_id = section_id,
      n_spots = nrow(tab),
      YAP17_requested_n = length(yap17), YAP17_measurable_n = length(yap_present),
      YAP17_missing = if (length(setdiff(yap17, yap_present))) paste(setdiff(yap17, yap_present), collapse = ";") else "NONE",
      Stem21_requested_n = length(stem21), Stem21_measurable_n = length(stem_present),
      Stem21_missing = if (length(setdiff(stem21, stem_present))) paste(setdiff(stem21, stem_present), collapse = ";") else "NONE",
      Program146_requested_feature_n = length(program146), Program146_measurable_n = length(p146$present),
      Program146_nonzero_variance_n = length(p146$nonzero),
      Program146_missing = if (length(p146$missing)) paste(p146$missing, collapse = ";") else "NONE",
      Program146_zero_sd = if (length(p146$zero_sd)) paste(p146$zero_sd, collapse = ";") else "NONE",
      Program142_requested_feature_n = length(program142), Program142_measurable_n = length(p142$present),
      Program142_nonzero_variance_n = length(p142$nonzero),
      Program142_missing = if (length(p142$missing)) paste(p142$missing, collapse = ";") else "NONE",
      Program142_zero_sd = if (length(p142$zero_sd)) paste(p142$zero_sd, collapse = ";") else "NONE",
      finite_Joint_fraction = mean(is.finite(joint)),
      finite_Program146_fraction = mean(is.finite(p146$score)),
      finite_Program142_fraction = mean(is.finite(p142$score)),
      Joint_mean = mean(joint), Joint_SD = sd(joint),
      Program146_mean = mean(p146$score), Program146_SD = sd(p146$score),
      Program142_mean = mean(p142$score), Program142_SD = sd(p142$score),
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
      score_QC_status = if (mean(is.finite(joint)) >= 0.95 && mean(is.finite(p146$score)) >= 0.95 && sd(joint) > 0 && sd(p146$score) > 0) "PASS" else "FAIL",
      stringsAsFactors = FALSE)
    write.table(tab, chunk_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
    write.table(meta, meta_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
    if (file.exists(err_path)) unlink(err_path)
    message(sprintf("WORKER %d DONE %s %s spots=%d p146=%d sec=%.1f", worker_id, sample_id, section_id, nrow(tab), length(p146$present), meta$elapsed_seconds))
    rm(umi, fit, corrected, expr, tab, meta, p146, p142)
    gc(verbose = FALSE)
  }, error = function(e) {
    writeLines(c(conditionMessage(e), capture.output(traceback())), err_path)
    message(sprintf("WORKER %d ERROR %s %s: %s", worker_id, sample_id, section_id, conditionMessage(e)))
  })
}
message(sprintf("WORKER %d COMPLETE", worker_id))
