#!/usr/bin/env Rscript
# Purpose: CNA matrix and identity preparation
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 1–2.

options(stringsAsFactors = FALSE, width = 220)
set.seed(20260728)
root <- Sys.getenv("AHIPPO_ROOT", ".")
project <- file.path(root, "0728_sensitivity_and_cnv_closure")
input_root <- file.path(project, "00_manifest", "per_patient_cnv_inputs")
dir.create(input_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project, "logs"), recursive = TRUE, showWarnings = FALSE)
tmp <- file.path(project, "tmp", "prepare_inputs_linux")
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = tmp, TEMP = tmp, TMP = tmp)

# Honour an explicitly supplied AHIPPO_R_LIBRARY before R_LIBS_USER and the
# project-local library.  Some workstation Rscript wrappers prepend a shared
# /opt library to R_LIBS_USER; the dedicated override preserves the user's
# requested first-priority library in that situation.
preferred_library <- Sys.getenv("AHIPPO_R_LIBRARY", unset = "")
preferred_libraries <- if (nzchar(preferred_library)) {
  strsplit(preferred_library, .Platform$path.sep, fixed = TRUE)[[1L]]
} else {
  character()
}
missing_preferred_libraries <- preferred_libraries[
  nzchar(preferred_libraries) & !dir.exists(preferred_libraries)
]
if (length(missing_preferred_libraries)) {
  stop(
    "AHIPPO_R_LIBRARY does not exist: ",
    paste(missing_preferred_libraries, collapse = .Platform$path.sep)
  )
}
preferred_libraries <- preferred_libraries[
  nzchar(preferred_libraries) & dir.exists(preferred_libraries)
]
requested_user_library <- Sys.getenv("R_LIBS_USER", unset = "")
requested_user_libraries <- if (nzchar(requested_user_library)) {
  strsplit(requested_user_library, .Platform$path.sep, fixed = TRUE)[[1L]]
} else {
  character()
}
requested_user_libraries <- requested_user_libraries[
  nzchar(requested_user_libraries) & dir.exists(requested_user_libraries)
]
.libPaths(unique(c(
  preferred_libraries,
  requested_user_libraries,
  file.path(project, "R_library"),
  .libPaths()
)))

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

log_path <- file.path(project, "logs", "linux_prepare_per_patient_inputs.log")
logcon <- file(log_path, "wt")
sink(logcon, split = TRUE)
sink(logcon, type = "message", append = TRUE)
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(logcon), silent = TRUE)
}, add = TRUE)

cat("Started", format(Sys.time()), "\n")
cat(
  "AHIPPO_R_LIBRARY:",
  Sys.getenv("AHIPPO_R_LIBRARY", unset = "<unset>"),
  "\n"
)
cat("R_LIBS_USER:", Sys.getenv("R_LIBS_USER", unset = "<unset>"), "\n")
cat("Resolved .libPaths():\n")
cat(paste0("  - ", .libPaths(), collapse = "\n"), "\n")
cat("Seurat loaded from:", find.package("Seurat"), "\n")
cat("Matrix loaded from:", find.package("Matrix"), "\n")
flush.console()

patients <- c(
  "CID3963", "CID4465", "CID4495", "CID44971",
  "CID44991", "CID4513", "CID4515", "CID4523"
)
object_path <- file.path(root, "Wu2021_seurat.rds")
frozen_path <- file.path(root, "Wu2021_continuous_state_cell_metadata.csv")
full_path <- file.path(
  project, "00_manifest", "Wu2021_full_atlas_umap_metadata.csv"
)
stopifnot(file.exists(object_path), file.exists(frozen_path), file.exists(full_path))

full_atlas <- read.csv(full_path, check.names = FALSE)
frozen <- read.csv(frozen_path, check.names = FALSE)

# The transferred metadata file is the complete Wu atlas (100,064 cells),
# not an already subsetted TNBC table.  Establish the audited hierarchy
# explicitly before constructing the eight evaluable-patient CNA inputs:
# complete atlas -> all TNBC cells (10 donors) -> 8 evaluable donors.
stopifnot(nrow(full_atlas) == 100064L)
stopifnot(!anyDuplicated(full_atlas$cell_id))
stopifnot(all(c(
  "cell_id", "patient", "subtype", "celltype_major", "celltype_minor"
) %in% colnames(full_atlas)))
full <- full_atlas[full_atlas$subtype == "TNBC", , drop = FALSE]
stopifnot(nrow(full) == 42512L)
stopifnot(length(unique(full$patient)) == 10L)
stopifnot(nrow(frozen) == 10836L)
stopifnot(setequal(unique(frozen$patient), patients))
stopifnot(all(frozen$cell_id %in% full$cell_id))

evaluable_full <- full[full$patient %in% patients, , drop = FALSE]
stopifnot(nrow(evaluable_full) == 39607L)
stopifnot(setequal(unique(evaluable_full$patient), patients))
cat(
  "Metadata hierarchy:",
  nrow(full_atlas), "complete-atlas cells ->",
  nrow(full), "TNBC cells across", length(unique(full$patient)), "donors ->",
  nrow(evaluable_full), "cells across 8 evaluable donors\n"
)
flush.console()

cat("Reading Seurat object", format(Sys.time()), "\n")
flush.console()
object <- readRDS(object_path)
counts <- tryCatch(
  SeuratObject::LayerData(object, assay = "RNA", layer = "counts"),
  error = function(error) Seurat::GetAssayData(object, assay = "RNA", slot = "counts")
)
stopifnot(inherits(counts, "sparseMatrix"))
stopifnot(all(counts@x >= 0))
stopifnot(all(abs(counts@x - round(counts@x)) < 1e-8))
stopifnot(all(evaluable_full$cell_id %in% colnames(counts)))
stopifnot(all(frozen$cell_id %in% colnames(counts)))

# Only the eight evaluable TNBC donors are used by the CNA callers.  Drop the
# remaining atlas cells and the full Seurat object before gene-symbol collapse
# and per-patient splitting to keep workstation memory bounded.
counts <- counts[, evaluable_full$cell_id, drop = FALSE]
rm(object, full_atlas)
invisible(gc())
cat(
  "Retained sparse raw-count matrix:",
  nrow(counts), "genes x", ncol(counts), "cells\n"
)
flush.console()

collapse_duplicate_symbols <- function(matrix) {
  symbols <- rownames(matrix)
  if (!anyDuplicated(symbols)) {
    return(list(matrix = matrix, n_duplicate_rows = 0L))
  }
  unique_symbols <- unique(symbols)
  triplet <- summary(matrix)
  new_row <- match(symbols, unique_symbols)
  collapsed <- sparseMatrix(
    i = new_row[triplet$i],
    j = triplet$j,
    x = triplet$x,
    dims = c(length(unique_symbols), ncol(matrix)),
    dimnames = list(unique_symbols, colnames(matrix))
  )
  list(
    matrix = collapsed,
    n_duplicate_rows = length(symbols) - length(unique_symbols)
  )
}
collapsed <- collapse_duplicate_symbols(counts)
counts <- collapsed$matrix

manifest <- list()
all_barcodes <- list()
partial_manifest_path <- file.path(
  project, "00_manifest", "CNV_per_patient_input_manifest.partial.tsv"
)
for (patient in patients) {
  cat("Preparing", patient, format(Sys.time()), "\n")
  flush.console()
  patient_dir <- file.path(input_root, patient)
  dir.create(patient_dir, recursive = TRUE, showWarnings = FALSE)
  patient_meta <- full[full$patient == patient, ]
  patient_meta$caller_role <- ifelse(
    patient_meta$cell_id %in% frozen$cell_id,
    "frozen_malignant",
    "same_patient_reference"
  )
  patient_meta$reference_group <- ifelse(
    patient_meta$caller_role == "frozen_malignant",
    "observation",
    as.character(patient_meta$celltype_major)
  )
  keep <- patient_meta$cell_id
  patient_counts <- counts[, keep, drop = FALSE]
  stopifnot(identical(colnames(patient_counts), keep))
  stopifnot(
    sum(patient_meta$caller_role == "frozen_malignant") ==
      sum(frozen$patient == patient)
  )

  counts_file <- file.path(patient_dir, paste0(patient, "_raw_counts_sparse.rds"))
  annotation_file <- file.path(patient_dir, paste0(patient, "_cell_annotations.tsv"))
  observation_file <- file.path(patient_dir, paste0(patient, "_frozen_malignant_barcodes.txt"))
  reference_file <- file.path(patient_dir, paste0(patient, "_same_patient_reference_barcodes.txt"))
  saveRDS(patient_counts, counts_file, compress = FALSE)
  write.table(
    patient_meta[, c(
      "cell_id", "patient", "subtype", "celltype_major", "celltype_minor",
      "caller_role", "reference_group"
    )],
    annotation_file,
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  writeLines(
    patient_meta$cell_id[patient_meta$caller_role == "frozen_malignant"],
    observation_file
  )
  writeLines(
    patient_meta$cell_id[patient_meta$caller_role == "same_patient_reference"],
    reference_file
  )
  manifest[[patient]] <- data.frame(
    patient = patient,
    genes = nrow(patient_counts),
    cells = ncol(patient_counts),
    frozen_malignant = sum(patient_meta$caller_role == "frozen_malignant"),
    same_patient_reference = sum(patient_meta$caller_role == "same_patient_reference"),
    duplicate_gene_rows_collapsed = collapsed$n_duplicate_rows,
    counts_file = counts_file,
    annotation_file = annotation_file,
    stringsAsFactors = FALSE
  )
  all_barcodes[[patient]] <- patient_meta[, c(
    "cell_id", "patient", "caller_role", "reference_group"
  )]
  write.table(
    do.call(rbind, manifest),
    partial_manifest_path,
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  rm(patient_counts)
  invisible(gc())
  cat(
    "Completed patient", patient,
    "(", nrow(manifest[[patient]]), "manifest row written )\n"
  )
  flush.console()
}

manifest <- do.call(rbind, manifest)
stopifnot(nrow(manifest) == 8L)
stopifnot(sum(manifest$frozen_malignant) == 10836L)
stopifnot(all(manifest$same_patient_reference >= 100L))
write.table(
  manifest,
  file.path(project, "00_manifest", "CNV_per_patient_input_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  do.call(rbind, all_barcodes),
  file.path(project, "00_manifest", "CNV_all_input_barcodes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(project, "00_manifest", "linux_input_preparation_sessionInfo.txt")
)
if (file.exists(partial_manifest_path)) {
  file.remove(partial_manifest_path)
}
cat("Completed", format(Sys.time()), "\n")
flush.console()
