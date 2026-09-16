# Purpose: Score predefined spatial functional modules
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 6 and Supplementary Fig. S4.
options(stringsAsFactors = FALSE)
.libPaths(c(".runtime/Figure6_S0Q3_Rlib", Sys.getenv("R_LIBS_USER"), .Library))

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

args <- commandArgs(trailingOnly = TRUE)
worker_id <- as.integer(args[1])
n_workers <- as.integer(args[2])
stopifnot(is.finite(worker_id), is.finite(n_workers), worker_id >= 1L, worker_id <= n_workers)

input_root <- ".runtime/Figure6_S1_inputs"
output_root <- "Figure6_FG_biological_enrichment_phase/intermediate/module_score_chunks"
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

registry <- read.delim(file.path(input_root, "section_registry.tsv"), check.names = FALSE)
registry$patient_number <- as.integer(sub("^P", "", registry$patient_id))
registry <- registry[order(registry$patient_number, registry$section_number), ]
module_table <- read.delim("0716_yan2026_validation/tables/module_score_independent_gene_lists.tsv", check.names = FALSE)
score_genes <- read.delim(file.path(input_root, "score_gene_list.tsv"), check.names = FALSE)
program146_original <- scan(file.path(input_root, "146_gene_program.txt"), what = character(), quiet = TRUE)

module_order <- c("UPR", "TNFA NFKB", "Hypoxia", "Adhesion Remodeling", "Survival Stress")
stopifnot(all(module_order %in% unique(module_table$module)))
module_table <- module_table[module_table$module %in% module_order, ]

alias_map <- c(CYR61 = "CCN1", CTGF = "CCN2", `HNRNPU-AS1` = "HNRNPU")
map_alias <- function(x) {
  y <- x
  hit <- x %in% names(alias_map)
  y[hit] <- unname(alias_map[x[hit]])
  y
}

modules <- setNames(lapply(module_order, function(m) unique(map_alias(module_table$gene[module_table$module == m]))), module_order)
yap17 <- unique(map_alias(score_genes$gene[score_genes$in_YAP_score]))
stem21 <- unique(map_alias(score_genes$gene[score_genes$in_Stemness_score]))
stem20_clean <- setdiff(stem21, "ZEB2")
program146 <- unique(map_alias(program146_original))
program144_clean <- setdiff(program146, c("DDR2", "TNFAIP3"))

stopifnot(length(yap17) == 17L, length(stem21) == 21L, length(stem20_clean) == 20L)
stopifnot(length(program146) == 146L, length(program144_clean) == 144L)
stopifnot(all(vapply(modules, function(x) length(intersect(x, union(yap17, stem21))) == 0L, logical(1))))

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

mean_gene_z_score <- function(expr, requested) {
  present <- intersect(requested, rownames(expr))
  missing <- setdiff(requested, present)
  if (!length(present)) stop("No requested genes measurable")
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
  prefix <- paste0(rec$sample_id, "_", rec$section_id)
  chunk_path <- file.path(output_root, paste0(prefix, "_module_scores.tsv"))
  meta_path <- file.path(output_root, paste0(prefix, "_module_meta.tsv"))
  err_path <- file.path(output_root, paste0(prefix, "_ERROR.txt"))
  if (file.exists(chunk_path) && file.exists(meta_path)) next
  message(sprintf("WORKER %d START %d/%d %s", worker_id, ri, nrow(registry), prefix))
  started <- Sys.time()
  tryCatch({
    h5_dir <- file.path(".runtime/Figure6_S0Q2_geo", rec$sample_id)
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
      n_cells = min(5000L, ncol(umi)), residual_type = "none",
      return_cell_attr = TRUE, return_gene_attr = TRUE,
      return_corrected_umi = FALSE, verbosity = 0)
    corrected <- as(correct_counts(fit, umi, verbosity = 0), "CsparseMatrix")
    expr <- log1p(corrected)

    module_results <- lapply(modules, function(g) mean_gene_z_score(expr, g))
    yap_present <- intersect(yap17, rownames(expr))
    stem_clean_present <- intersect(stem20_clean, rownames(expr))
    if (!length(yap_present) || !length(stem_clean_present)) stop("No measurable clean Joint-axis genes")
    yap_mean <- Matrix::colMeans(expr[yap_present, , drop = FALSE])
    stem_clean_mean <- Matrix::colMeans(expr[stem_clean_present, , drop = FALSE])
    joint_clean <- z_vector(yap_mean) + z_vector(stem_clean_mean)
    program_clean <- mean_gene_z_score(expr, program144_clean)

    tab <- data.frame(
      patient_id = rec$patient_id, sample_id = rec$sample_id, section_id = rec$section_id,
      barcode = colnames(umi), nCount_raw = as.numeric(n_count), nFeature_raw = as.integer(n_feature),
      UPR = as.numeric(module_results[["UPR"]]$score),
      TNFA_NFKB = as.numeric(module_results[["TNFA NFKB"]]$score),
      Hypoxia = as.numeric(module_results[["Hypoxia"]]$score),
      Adhesion_Remodeling = as.numeric(module_results[["Adhesion Remodeling"]]$score),
      Survival_Stress = as.numeric(module_results[["Survival Stress"]]$score),
      Joint_spatial_no_ZEB2 = as.numeric(joint_clean),
      Program146_spatial_no_DDR2_TNFAIP3 = as.numeric(program_clean$score),
      stringsAsFactors = FALSE)

    meta_rows <- lapply(names(module_results), function(m) {
      z <- module_results[[m]]
      data.frame(patient_id = rec$patient_id, sample_id = rec$sample_id, section_id = rec$section_id,
                 object = m, requested_n = length(modules[[m]]), measurable_n = length(z$present),
                 nonzero_variance_n = length(z$nonzero), measurable_genes = paste(z$present, collapse = ";"),
                 missing_genes = if (length(z$missing)) paste(z$missing, collapse = ";") else "NONE",
                 zero_sd_genes = if (length(z$zero_sd)) paste(z$zero_sd, collapse = ";") else "NONE",
                 stringsAsFactors = FALSE)
    })
    meta_rows[[length(meta_rows) + 1L]] <- data.frame(
      patient_id = rec$patient_id, sample_id = rec$sample_id, section_id = rec$section_id,
      object = "Joint_spatial_no_ZEB2", requested_n = length(yap17) + length(stem20_clean),
      measurable_n = length(yap_present) + length(stem_clean_present), nonzero_variance_n = NA_integer_,
      measurable_genes = paste(c(yap_present, stem_clean_present), collapse = ";"),
      missing_genes = paste(c(setdiff(yap17, yap_present), setdiff(stem20_clean, stem_clean_present)), collapse = ";"),
      zero_sd_genes = "NOT_APPLICABLE", stringsAsFactors = FALSE)
    meta_rows[[length(meta_rows) + 1L]] <- data.frame(
      patient_id = rec$patient_id, sample_id = rec$sample_id, section_id = rec$section_id,
      object = "Program146_spatial_no_DDR2_TNFAIP3", requested_n = length(program144_clean),
      measurable_n = length(program_clean$present), nonzero_variance_n = length(program_clean$nonzero),
      measurable_genes = paste(program_clean$present, collapse = ";"),
      missing_genes = if (length(program_clean$missing)) paste(program_clean$missing, collapse = ";") else "NONE",
      zero_sd_genes = if (length(program_clean$zero_sd)) paste(program_clean$zero_sd, collapse = ";") else "NONE",
      stringsAsFactors = FALSE)
    meta <- do.call(rbind, meta_rows)
    meta$n_spots <- nrow(tab)
    meta$elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    meta$status <- ifelse(meta$measurable_n > 0, "PASS", "FAIL")

    write.table(tab, chunk_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
    write.table(meta, meta_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
    if (file.exists(err_path)) unlink(err_path)
    message(sprintf("WORKER %d DONE %s spots=%d sec=%.1f", worker_id, prefix, nrow(tab), max(meta$elapsed_seconds)))
    rm(umi, fit, corrected, expr, tab, meta, module_results, program_clean)
    gc(verbose = FALSE)
  }, error = function(e) {
    writeLines(c(conditionMessage(e), capture.output(traceback())), err_path)
    message(sprintf("WORKER %d ERROR %s: %s", worker_id, prefix, conditionMessage(e)))
  })
}
message(sprintf("WORKER %d COMPLETE", worker_id))
