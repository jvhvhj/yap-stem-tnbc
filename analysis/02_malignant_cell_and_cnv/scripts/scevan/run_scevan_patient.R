#!/usr/bin/env Rscript
# Purpose: Anchored patient-wise SCEVAN inference
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 2 and Supplementary Fig. S1.

options(stringsAsFactors = FALSE, width = 220)
set.seed(20260728)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: Rscript 07_run_scevan_patient.R PATIENT [anchored|unanchored]")
patient <- args[1]
mode <- if (length(args) >= 2L) args[2] else "anchored"
patients <- c(
  "CID3963", "CID4465", "CID4495", "CID44971",
  "CID44991", "CID4513", "CID4515", "CID4523"
)
stopifnot(patient %in% patients, mode %in% c("anchored", "unanchored"))

default_root <- if (.Platform$OS.type == "windows") "." else "."
root <- normalizePath(
  Sys.getenv("AHIPPO_ROOT", default_root),
  winslash = "/",
  mustWork = TRUE
)
project <- file.path(root, "0728_sensitivity_and_cnv_closure")
input_dir <- file.path(project, "00_manifest", "per_patient_cnv_inputs", patient)
patient_dir <- file.path(project, "03_SCEVAN", patient)
if (mode == "unanchored") patient_dir <- file.path(patient_dir, "unanchored_sensitivity")
native_dir <- file.path(patient_dir, "native_output")
dir.create(native_dir, recursive = TRUE, showWarnings = FALSE)
tmp <- file.path(project, "tmp", paste0("SCEVAN_", patient, "_", mode))
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = tmp, TEMP = tmp, TMP = tmp)
.libPaths(c(
  file.path(project, "R_library"),
  file.path(root, "0703_rebuild", "R_library"),
  ".software/r-library",
  .libPaths()
))

log_path <- file.path(patient_dir, "run_log.txt")
logcon <- file(log_path, "wt")
sink(logcon, split = TRUE)
sink(logcon, type = "message", append = TRUE)
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(logcon), silent = TRUE)
}, add = TRUE)

write_status <- function(status, reason, result_class = NA_character_, n_calls = NA_integer_) {
  row <- data.frame(
    patient = patient,
    caller = "SCEVAN",
    mode = mode,
    status = status,
    reason = reason,
    result_class = result_class,
    n_cell_calls = n_calls,
    seed = 20260728,
    cores = 12,
    stringsAsFactors = FALSE
  )
  write.table(
    row,
    file.path(patient_dir, "patient_status.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}

cat("Started", format(Sys.time()), "\n")
cat("Patient", patient, "mode", mode, "seed 20260728\n")
if (!requireNamespace("SCEVAN", quietly = TRUE)) {
  write_status(
    "NOT_TESTABLE",
    "SCEVAN package unavailable in the F-drive project library; run installation on the 128-GB workstation"
  )
  writeLines(capture.output(sessionInfo()), file.path(patient_dir, "sessionInfo.txt"))
  cat("NOT_TESTABLE: SCEVAN package unavailable\n")
  quit(save = "no", status = 0)
}

counts_file <- file.path(input_dir, paste0(patient, "_raw_counts_sparse.rds"))
annotation_file <- file.path(input_dir, paste0(patient, "_cell_annotations.tsv"))
counts <- readRDS(counts_file)
annotation <- read.delim(annotation_file, check.names = FALSE)
stopifnot(identical(colnames(counts), annotation$cell_id))
stopifnot(all(counts@x >= 0), all(abs(counts@x - round(counts@x)) < 1e-8))
observation <- annotation$cell_id[annotation$caller_role == "frozen_malignant"]
reference <- annotation$cell_id[annotation$caller_role == "same_patient_reference"]
if (!length(observation) || length(reference) < 100L) {
  write_status("NOT_TESTABLE", "Insufficient frozen malignant or same-patient reference cells")
  quit(save = "no", status = 0)
}

max_reference <- suppressWarnings(as.integer(Sys.getenv("SCEVAN_MAX_REFERENCE", "NA")))
if (!is.na(max_reference) && length(reference) > max_reference) {
  groups <- split(reference, annotation$reference_group[match(reference, annotation$cell_id)])
  allocation <- pmax(1L, floor(max_reference * lengths(groups) / length(reference)))
  while (sum(allocation) > max_reference) {
    index <- which.max(allocation)
    allocation[index] <- allocation[index] - 1L
  }
  while (sum(allocation) < max_reference) {
    capacity <- lengths(groups) - allocation
    index <- which.max(capacity)
    allocation[index] <- allocation[index] + 1L
  }
  reference <- unlist(Map(function(x, n) sample(x, n), groups, allocation), use.names = FALSE)
  cat("Reference downsampling enabled:", length(reference), "cells; stratified by reference_group\n")
} else {
  cat("All same-patient reference cells retained:", length(reference), "\n")
}
writeLines(observation, file.path(patient_dir, "actual_frozen_malignant_barcodes.txt"))
writeLines(reference, file.path(patient_dir, "actual_reference_barcodes.txt"))

native <- tryCatch(
  SCEVAN::pipelineCNA(
    count_mtx = counts,
    sample = patient,
    par_cores = 12,
    norm_cell = if (mode == "anchored") reference else NULL,
    SUBCLONES = FALSE,
    beta_vega = 0.5,
    ClonalCN = TRUE,
    plotTree = FALSE,
    organism = "human",
    FIXED_NORMAL_CELLS = FALSE,
    output_dir = native_dir
  ),
  error = function(error) error
)
saveRDS(native, file.path(patient_dir, "native_result.rds"), compress = FALSE)
writeLines(
  capture.output(str(native, max.level = 4)),
  file.path(patient_dir, "native_result_structure.txt")
)

if (inherits(native, "error")) {
  write_status("NOT_TESTABLE", conditionMessage(native), class(native)[1])
  writeLines(capture.output(sessionInfo()), file.path(patient_dir, "sessionInfo.txt"))
  quit(save = "no", status = 0)
}

calls_testable <- is.data.frame(native) &&
  !is.null(rownames(native)) &&
  "class" %in% colnames(native) &&
  all(unique(na.omit(native$class)) %in% c("tumor", "normal", "filtered")) &&
  sum(rownames(native) %in% colnames(counts)) > 0
if (calls_testable) {
  calls <- data.frame(
    cell_id = rownames(native),
    patient = patient,
    scevan_native_class = as.character(native$class),
    scevan_tumor_positive = as.character(native$class) == "tumor",
    caller_role = annotation$caller_role[match(rownames(native), annotation$cell_id)],
    returned_by_native_result = TRUE,
    stringsAsFactors = FALSE
  )
  write.table(
    calls,
    file.path(patient_dir, "cell_calls.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
} else {
  calls <- data.frame()
  write.table(
    data.frame(
      status = "NOT_TESTABLE",
      reason = "Native object did not expose an auditable cell-indexed class column with documented tumor/normal/filtered semantics"
    ),
    file.path(patient_dir, "cell_calls.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

# pipelineCNA 1.0.3 returns a classification data.frame and does not return its
# cell-level CNA matrix. Preserve an explicit status object instead of guessing
# from unrelated generated files.
cna_status <- list(
  status = "NOT_TESTABLE_FROM_NATIVE_RETURN",
  reason = "SCEVAN pipelineCNA 1.0.3 native return contains cell classes but no cell-level CNA matrix",
  source_version = as.character(utils::packageVersion("SCEVAN"))
)
saveRDS(cna_status, file.path(patient_dir, "CNA_matrix.rds"))
writeLines(capture.output(sessionInfo()), file.path(patient_dir, "sessionInfo.txt"))
write_status(
  if (calls_testable) "PASS_CLASSIFICATION" else "NOT_TESTABLE",
  if (calls_testable) "Audited native tumor/normal/filtered class table" else "Unresolved native cell-call semantics",
  paste(class(native), collapse = ";"),
  if (calls_testable) nrow(calls) else 0L
)
cat("Completed", format(Sys.time()), "\n")
