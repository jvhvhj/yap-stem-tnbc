#!/usr/bin/env Rscript
# Purpose: Continuous inferCNV burden and native HMM status
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 2 and Supplementary Fig. S1.

options(stringsAsFactors = FALSE, width = 220)
set.seed(20260728)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: Rscript 08_run_infercnv_patient.R PATIENT [HG38_GENE_ORDER]")
patient <- args[1]
patients <- c(
  "CID3963", "CID4465", "CID4495", "CID44971",
  "CID44991", "CID4513", "CID4515", "CID4523"
)
stopifnot(patient %in% patients)

default_root <- if (.Platform$OS.type == "windows") "." else "."
root <- normalizePath(
  Sys.getenv("AHIPPO_ROOT", default_root),
  winslash = "/",
  mustWork = TRUE
)
project <- file.path(root, "0728_sensitivity_and_cnv_closure")
gene_order_file <- if (length(args) >= 2L) {
  normalizePath(args[2], winslash = "/", mustWork = TRUE)
} else {
  file.path(project, "00_manifest", "hg38_gene_order_UCSC_refGene.tsv")
}
input_dir <- file.path(project, "00_manifest", "per_patient_cnv_inputs", patient)
patient_dir <- file.path(project, "04_inferCNV", patient)
continuous_dir <- file.path(patient_dir, "continuous_native")
hmm_dir <- file.path(patient_dir, "hmm_native")
dir.create(continuous_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(hmm_dir, recursive = TRUE, showWarnings = FALSE)
tmp <- file.path(project, "tmp", paste0("inferCNV_", patient))
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

write_status <- function(continuous_status, hmm_status, reason, mapped_genes = NA_integer_) {
  row <- data.frame(
    patient = patient,
    caller = "inferCNV",
    continuous_status = continuous_status,
    hmm_status = hmm_status,
    reason = reason,
    mapped_genes = mapped_genes,
    neutral_baseline = 1,
    burden_formula = "mean(abs(signal - 1))",
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
cat("Patient", patient, "seed 20260728\n")
cat("Gene order", gene_order_file, "\n")
if (!requireNamespace("infercnv", quietly = TRUE)) {
  write_status(
    "NOT_TESTABLE",
    "NOT_TESTABLE",
    "infercnv unavailable in F-drive project library; install JAGS 4.x, rjags and infercnv on the 128-GB workstation"
  )
  writeLines(capture.output(sessionInfo()), file.path(patient_dir, "sessionInfo.txt"))
  cat("NOT_TESTABLE: infercnv package unavailable\n")
  quit(save = "no", status = 0)
}

counts <- readRDS(file.path(input_dir, paste0(patient, "_raw_counts_sparse.rds")))
annotation <- read.delim(
  file.path(input_dir, paste0(patient, "_cell_annotations.tsv")),
  check.names = FALSE
)
stopifnot(identical(colnames(counts), annotation$cell_id))
stopifnot(all(counts@x >= 0), all(abs(counts@x - round(counts@x)) < 1e-8))
gene_order <- read.delim(
  gene_order_file,
  header = FALSE,
  col.names = c("gene", "chromosome", "start", "stop"),
  check.names = FALSE
)
gene_order <- gene_order[
  !duplicated(gene_order$gene) & gene_order$gene %in% rownames(counts),
]
if (nrow(gene_order) < 5000L) {
  write_status("NOT_TESTABLE", "NOT_TESTABLE", "Fewer than 5,000 raw-count genes map to audited hg38 gene order", nrow(gene_order))
  quit(save = "no", status = 0)
}
counts <- counts[gene_order$gene, , drop = FALSE]
observation <- annotation$cell_id[annotation$caller_role == "frozen_malignant"]
reference <- annotation$cell_id[annotation$caller_role == "same_patient_reference"]
actual_group <- ifelse(
  annotation$caller_role == "frozen_malignant",
  "observation",
  paste0("reference_", gsub("[^A-Za-z0-9]+", "_", annotation$reference_group))
)
infer_annotation <- data.frame(cell = annotation$cell_id, group = actual_group)
annotation_path <- file.path(patient_dir, "inferCNV_annotations.tsv")
gene_order_path <- file.path(patient_dir, "inferCNV_gene_order_hg38.tsv")
write.table(
  infer_annotation,
  annotation_path,
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)
write.table(
  gene_order,
  gene_order_path,
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)
reference_groups <- sort(unique(actual_group[annotation$caller_role == "same_patient_reference"]))

base_object <- tryCatch(
  infercnv::CreateInfercnvObject(
    raw_counts_matrix = counts,
    annotations_file = annotation_path,
    delim = "\t",
    gene_order_file = gene_order_path,
    ref_group_names = reference_groups
  ),
  error = function(error) error
)
if (inherits(base_object, "error")) {
  write_status("NOT_TESTABLE", "NOT_TESTABLE", conditionMessage(base_object), nrow(gene_order))
  quit(save = "no", status = 0)
}

continuous <- tryCatch(
  infercnv::run(
    base_object,
    cutoff = 0.1,
    out_dir = continuous_dir,
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = FALSE,
    num_threads = 12,
    no_plot = FALSE
  ),
  error = function(error) error
)
if (inherits(continuous, "error")) {
  saveRDS(continuous, file.path(patient_dir, "continuous_object.rds"))
  write_status("NOT_TESTABLE", "NOT_TESTABLE", conditionMessage(continuous), nrow(gene_order))
  quit(save = "no", status = 0)
}
saveRDS(continuous, file.path(patient_dir, "continuous_object.rds"), compress = FALSE)
continuous_matrix <- continuous@expr.data
saveRDS(continuous_matrix, file.path(patient_dir, "continuous_matrix.rds"), compress = FALSE)
burden <- colMeans(abs(continuous_matrix - 1), na.rm = TRUE)
reference_present <- intersect(reference, names(burden))
threshold <- unname(quantile(burden[reference_present], 0.95, na.rm = TRUE, names = FALSE))
burden_table <- data.frame(
  cell_id = names(burden),
  patient = patient,
  caller_role = ifelse(names(burden) %in% observation, "frozen_malignant", "same_patient_reference"),
  inferCNV_continuous_burden = as.numeric(burden),
  neutral_baseline = 1,
  reference_95pct_threshold = threshold,
  inferCNV_descriptive_CNV_high = as.numeric(burden) > threshold,
  threshold_interpretation = "descriptive within-patient reference 95th percentile; not a validated malignancy threshold",
  stringsAsFactors = FALSE
)
write.table(
  burden_table,
  file.path(patient_dir, "cell_cnv_burden.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  burden_table[burden_table$caller_role == "same_patient_reference", ],
  file.path(patient_dir, "reference_burden_distribution.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

hmm_hours <- suppressWarnings(as.numeric(Sys.getenv("INFERCNV_HMM_TIME_LIMIT_HOURS", "24")))
setTimeLimit(elapsed = hmm_hours * 3600, transient = TRUE)
hmm <- tryCatch(
  infercnv::run(
    base_object,
    cutoff = 0.1,
    out_dir = hmm_dir,
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = TRUE,
    HMM_type = "i6",
    num_threads = 12,
    no_plot = FALSE
  ),
  error = function(error) error
)
setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
if (inherits(hmm, "error")) {
  hmm_table <- data.frame(
    cell_id = colnames(counts),
    patient = patient,
    hmm_testable = FALSE,
    inferCNV_HMM_positive = NA,
    reason = conditionMessage(hmm),
    stringsAsFactors = FALSE
  )
  hmm_status <- "NOT_TESTABLE"
  hmm_reason <- conditionMessage(hmm)
} else {
  saveRDS(hmm, file.path(patient_dir, "hmm_object.rds"), compress = FALSE)
  # inferCNV HMM states are gene/segment-level. A binary malignant-cell call
  # requires a separately validated threshold, which is not introduced here.
  hmm_table <- data.frame(
    cell_id = colnames(hmm@expr.data),
    patient = patient,
    hmm_testable = FALSE,
    inferCNV_HMM_positive = NA,
    reason = "HMM object succeeded, but inferCNV provides state-level CNA predictions rather than a validated binary malignant-cell call",
    stringsAsFactors = FALSE
  )
  hmm_status <- "PASS_NATIVE_OBJECT_BINARY_CALL_NOT_TESTABLE"
  hmm_reason <- unique(hmm_table$reason)
}
write.table(
  hmm_table,
  file.path(patient_dir, "hmm_cell_calls.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
writeLines(capture.output(sessionInfo()), file.path(patient_dir, "sessionInfo.txt"))
write_status("PASS", hmm_status, hmm_reason, nrow(gene_order))
cat("Completed", format(Sys.time()), "\n")
