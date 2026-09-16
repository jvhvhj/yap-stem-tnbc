#!/usr/bin/env Rscript
# YAP17 and Stem21 component scoring.
#
# Definition (fixed; do not reselect genes):
#   score = unweighted arithmetic mean of the measurable signature genes in the
#           log-normalised RNA expression layer of the dataset being scored.
#
# For the Wu discovery cohort this is the RNA assay "data" layer (log-normalised).
# Raw counts are NOT used. The signature files are consumed exactly as shipped in
# config/signatures/; no alias substitution, fuzzy matching or subtype-specific
# optimisation is applied, and genes absent from the object remain absent from the
# mean rather than being replaced.
#
# This script produces the two component scores only. The Joint axis and the
# Low/Intermediate/High state assignment are derived separately in
# derive_joint_axis_and_tertiles.R.

options(stringsAsFactors = FALSE, width = 220)

usage <- function() cat(paste(
  "Usage: Rscript analysis/03_yap_stem_continuum/scripts/score_yap17_stem21.R",
  "  --counts <Seurat_RDS>",
  "  --output-dir <NEW_output_directory>",
  "  [--yap-signature config/signatures/YAP17.tsv]",
  "  [--stem-signature config/signatures/Stem21.tsv]",
  "  [--assay RNA]",
  "  [--layer data]",
  "  [--dry-run]",
  "",
  "Writes component_scores.tsv with one row per cell:",
  "  cell_id, patient, n_yap_measurable, n_stem_measurable, YAP_score, Stemness_score",
  "Run from the repository root. Dry-run validates the signature files only.",
  sep = "\n"), "\n")

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) { usage(); quit(status = 0L) }
if (length(args) == 0L) { usage(); quit(status = 2L) }

flags <- list()
allowed <- c("counts", "output-dir", "yap-signature", "stem-signature", "assay", "layer", "dry-run")
i <- 1L
while (i <= length(args)) {
  key <- sub("^--", "", args[i])
  if (!startsWith(args[i], "--") || !key %in% allowed || key %in% names(flags))
    stop("Unknown or duplicate option: ", args[i])
  if (key == "dry-run") { flags[[key]] <- TRUE; i <- i + 1L; next }
  if (i == length(args) || startsWith(args[i + 1L], "--"))
    stop("Missing value for --", key)
  flags[[key]] <- args[i + 1L]
  i <- i + 2L
}
required <- c("counts", "output-dir")
if (!all(required %in% names(flags)))
  stop("Required options missing: ", paste(setdiff(required, names(flags)), collapse = ", "))

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "config", "signatures", "YAP17.tsv")))
  stop("Run this command from the publication repository root.")

read_sig <- function(path) {
  if (!file.exists(path)) stop("Signature file not found: ", path)
  genes <- read.delim(path, check.names = FALSE)$gene
  genes <- unique(as.character(genes[nzchar(genes)]))
  if (!length(genes)) stop("Signature file is empty: ", path)
  genes
}

yap_path <- flags[["yap-signature"]]
if (is.null(yap_path)) yap_path <- file.path(root, "config", "signatures", "YAP17.tsv")

stem_path <- flags[["stem-signature"]]
if (is.null(stem_path)) stem_path <- file.path(root, "config", "signatures", "Stem21.tsv")

yap_genes  <- read_sig(normalizePath(yap_path,  winslash = "/", mustWork = TRUE))
stem_genes <- read_sig(normalizePath(stem_path, winslash = "/", mustWork = TRUE))

overlap <- intersect(yap_genes, stem_genes)
if (length(overlap)) stop("Signature files overlap: ", paste(overlap, collapse = ", "))

cat("YAP signature genes :", length(yap_genes), "\n")
cat("Stem signature genes:", length(stem_genes), "\n")

if (isTRUE(flags[["dry-run"]])) {
  cat("\n[dry-run] signature validation only; no outputs written.\n")
  quit(status = 0L)
}

suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(Matrix) })

obj_path <- normalizePath(flags[["counts"]], winslash = "/", mustWork = TRUE)
assay <- if (!is.null(flags[["assay"]])) flags[["assay"]] else "RNA"
layer <- if (!is.null(flags[["layer"]])) flags[["layer"]] else "data"

obj <- readRDS(obj_path)
if (!assay %in% names(obj@assays)) stop("Assay not present in object: ", assay)
mat <- SeuratObject::LayerData(obj, assay = assay, layer = layer)
available <- rownames(mat)

yap_use  <- intersect(yap_genes,  available)
stem_use <- intersect(stem_genes, available)
cat("measurable YAP genes :", length(yap_use), "/", length(yap_genes), "\n")
cat("measurable Stem genes:", length(stem_use), "/", length(stem_genes), "\n")
if (!length(yap_use) || !length(stem_use))
  stop("No measurable genes for at least one signature in this dataset.")

# Unweighted arithmetic mean over the measurable members of each signature.
yap_score  <- Matrix::colMeans(mat[yap_use,  , drop = FALSE])
stem_score <- Matrix::colMeans(mat[stem_use, , drop = FALSE])

out_dir <- flags[["output-dir"]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)

# Patient column: prefer an explicit 'patient' field, else the object's
# orig.ident, which is where the Wu per-donor identifier is stored.
patient_vec <- NULL
for (cand in c("patient", "orig.ident")) {
  if (cand %in% names(obj@meta.data)) { patient_vec <- as.character(obj@meta.data[[cand]]); break }
}
if (is.null(patient_vec)) {
  warning("No 'patient' or 'orig.ident' column found; patient column left blank.")
  patient_vec <- rep(NA_character_, ncol(obj))
}

df <- data.frame(
  cell_id = colnames(obj),
  patient = patient_vec,
  n_yap_measurable = length(yap_use),
  n_stem_measurable = length(stem_use),
  YAP_score = as.numeric(yap_score),
  Stemness_score = as.numeric(stem_score),
  stringsAsFactors = FALSE)

write.table(df, file.path(out_dir, "component_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(data.frame(
  signature = c("YAP17", "Stem21"),
  declared_genes = c(length(yap_genes), length(stem_genes)),
  measurable_genes = c(length(yap_use), length(stem_use)),
  missing_genes = c(paste(setdiff(yap_genes, available), collapse = ","),
                    paste(setdiff(stem_genes, available), collapse = ",")),
  assay = assay, layer = layer, stringsAsFactors = FALSE),
  file.path(out_dir, "signature_coverage.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)

cat("\ncells scored:", nrow(df), "\n")
cat("outputs written to:", out_dir, "\n")
