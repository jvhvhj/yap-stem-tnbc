#!/usr/bin/env Rscript
# Purpose: Reference-subsampled CopyKAT inference
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 2 and Supplementary Fig. S1.

options(stringsAsFactors = FALSE, width = 220)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 05c_copykat_one_patient_subsampled_robust.R <PATIENT_ID> <PROJECT_ROOT> [N_CORES] [REF_CAP] [SEED]")
}

patient_id <- args[1]
project_root <- normalizePath(args[2], winslash = "/", mustWork = TRUE)
n_cores <- if (length(args) >= 3) as.integer(args[3]) else 1L
ref_cap <- if (length(args) >= 4) as.integer(args[4]) else 1000L
seed <- if (length(args) >= 5) as.integer(args[5]) else 123L
if (is.na(n_cores) || n_cores < 1L) stop("N_CORES must be a positive integer.")
if (is.na(ref_cap) || ref_cap < 50L) stop("REF_CAP must be an integer >= 50.")
if (is.na(seed)) stop("SEED must be an integer.")

candidate_libs <- unique(c(
  Sys.getenv("R_LIBS_USER", unset = NA_character_),
  ".software/r-library",
  file.path(project_root, "0703_rebuild", "R_library")
))
candidate_libs <- candidate_libs[!is.na(candidate_libs) & dir.exists(candidate_libs)]
if (length(candidate_libs)) .libPaths(c(candidate_libs, .libPaths()))

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(copykat)
})

seurat_path <- file.path(project_root, "Wu2021_seurat.rds")
out_root <- file.path(project_root, "0703_rebuild", "workstation_results", "copykat_full")
out_dir <- file.path(out_root, paste0(patient_id, "_native_output"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

prediction_path <- file.path(out_root, paste0(patient_id, "_copykat_prediction_subsampled.csv"))
input_meta_path <- file.path(out_root, paste0(patient_id, "_input_cells_subsampled.csv"))
parameter_path <- file.path(out_root, paste0(patient_id, "_copykat_parameters_subsampled.txt"))
session_path <- file.path(out_root, paste0(patient_id, "_copykat_sessionInfo_subsampled.txt"))

message("Loading: ", seurat_path)
obj <- readRDS(seurat_path)
meta <- obj@meta.data

patient_candidates <- intersect(c("patient", "orig.ident"), colnames(meta))
patient_matches <- patient_candidates[vapply(
  patient_candidates,
  function(x) patient_id %in% as.character(meta[[x]]),
  logical(1)
)]
if (!length(patient_matches)) {
  stop("Patient not found in either patient or orig.ident metadata: ", patient_id)
}
patient_col <- patient_matches[1]
message("Patient column: ", patient_col)

if (!"celltype_major" %in% colnames(meta)) stop("Missing metadata column: celltype_major")
patient_values <- as.character(meta[[patient_col]])
celltype_values <- as.character(meta$celltype_major)
patient_cells <- rownames(meta)[!is.na(patient_values) & patient_values == patient_id]
patient_meta <- meta[patient_cells, , drop = FALSE]
patient_celltype <- as.character(patient_meta$celltype_major)

cancer_cells <- rownames(patient_meta)[!is.na(patient_celltype) & patient_celltype == "Cancer Epithelial"]
ref_cells <- rownames(patient_meta)[!is.na(patient_celltype) & patient_celltype != "Cancer Epithelial"]
if (!length(cancer_cells)) stop("No Cancer Epithelial cells for ", patient_id)
if (length(ref_cells) < 50L) stop("Fewer than 50 same-patient reference cells for ", patient_id)

set.seed(seed)
ref_cells_for_copykat <- if (length(ref_cells) > ref_cap) sample(ref_cells, ref_cap) else ref_cells
cells_for_copykat <- c(cancer_cells, ref_cells_for_copykat)
if (anyDuplicated(cells_for_copykat)) stop("Duplicated cells in CopyKAT input set.")

raw_counts <- tryCatch(
  GetAssayData(obj, assay = "RNA", layer = "counts"),
  error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
)
missing_raw <- setdiff(cells_for_copykat, colnames(raw_counts))
if (length(missing_raw)) {
  stop(length(missing_raw), " selected cells are absent from RNA counts. Example: ", missing_raw[1])
}
raw_counts <- raw_counts[, cells_for_copykat, drop = FALSE]

input_meta <- meta[cells_for_copykat, , drop = FALSE]
input_meta$cell_id <- rownames(input_meta)
input_meta$copykat_input_role <- ifelse(
  input_meta$cell_id %in% cancer_cells,
  "Cancer_Epithelial",
  "Reference_sampled"
)
write.csv(input_meta, input_meta_path, row.names = FALSE, quote = TRUE)

writeLines(c(
  paste("patient_id", patient_id),
  paste("project_root", project_root),
  paste("patient_column", patient_col),
  paste("n_cores", n_cores),
  paste("ref_cap", ref_cap),
  paste("seed", seed),
  paste("cancer_epithelial_cells", length(cancer_cells)),
  paste("reference_cells_available", length(ref_cells)),
  paste("reference_cells_sampled", length(ref_cells_for_copykat)),
  paste("copykat_version", as.character(packageVersion("copykat"))),
  paste("start_time", Sys.time())
), parameter_path)

message("CopyKAT input: ", nrow(raw_counts), " genes x ", ncol(raw_counts), " cells")
copykat_mat <- as.matrix(raw_counts)
storage.mode(copykat_mat) <- "numeric"
rm(raw_counts, obj)
gc()

old_wd <- getwd()
setwd(out_dir)
on.exit(setwd(old_wd), add = TRUE)

ck <- copykat::copykat(
  rawmat = copykat_mat,
  id.type = "S",
  cell.line = "no",
  ngene.chr = 5,
  win.size = 25,
  KS.cut = 0.1,
  sam.name = patient_id,
  distance = "euclidean",
  norm.cell.names = ref_cells_for_copykat,
  n.cores = n_cores,
  genome = "hg20"
)

if (is.null(ck$prediction) || !nrow(as.data.frame(ck$prediction))) {
  stop("CopyKAT returned no prediction rows for ", patient_id)
}
pred <- as.data.frame(ck$prediction)
write.csv(pred, prediction_path, row.names = FALSE, quote = TRUE)
message("Saved prediction before attempting large-object serialization: ", prediction_path)

object_path <- file.path(out_dir, paste0(patient_id, "_copykat_object.rds"))
tryCatch(
  saveRDS(ck, object_path, compress = FALSE),
  error = function(e) warning("Prediction is safe, but CopyKAT object could not be saved: ", conditionMessage(e))
)
writeLines(capture.output(sessionInfo()), session_path)
message("Completed: ", patient_id)
