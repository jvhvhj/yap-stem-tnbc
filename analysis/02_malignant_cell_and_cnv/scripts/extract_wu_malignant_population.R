#!/usr/bin/env Rscript
# Wu malignant epithelial population extraction (Fig. 2).
#
# This is POPULATION EXTRACTION, not a biological derivation algorithm. The
# malignant identity is pre-existing: it comes from the source object's own
# author annotation. Two annotation fields select the population and nothing
# else is applied.
#
#   subtype        == "TNBC"                 (source-object TNBC annotation)
#   celltype_major == "Cancer Epithelial"    (source-object epithelial annotation)
#
# The eight evaluable discovery donors are not filtered separately. Donors whose
# TNBC subset contains no Cancer Epithelial cells drop out on their own, which is
# why the 10-donor TNBC subset yields an 8-donor malignant population.
#
# RNA-CNA caller results (CopyKAT / SCEVAN / inferCNV) do NOT define inclusion.
# They are supporting sensitivity evidence and are computed downstream in
# analysis/02_malignant_cell_and_cnv.
#
# No score is calculated, no cell is relabelled, and no donor is reclassified here.

options(stringsAsFactors = FALSE, width = 220)

usage <- function() cat(paste(
  "Usage: Rscript analysis/02_malignant_cell_and_cnv/scripts/extract_wu_malignant_population.R",
  "  --atlas-metadata <Wu2021_full_atlas_umap_metadata.csv>",
  "  --output-dir <NEW_output_directory>",
  "  [--atlas-object <Wu2021_full_atlas.rds>]",
  "  [--dry-run]",
  "",
  "The atlas metadata table must contain at least the columns:",
  "  cell_id, patient, subtype, celltype_major",
  "",
  "Writes cell_ids.tsv, per_patient_counts.tsv and extraction_summary.tsv into",
  "--output-dir. With --atlas-object, also writes the subset Seurat object.",
  "Dry-run validates the input table only; it writes no outputs.",
  sep = "\n"), "\n")

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) { usage(); quit(status = 0L) }
if (length(args) == 0L) { usage(); quit(status = 2L) }

flags <- list()
allowed <- c("atlas-metadata", "output-dir", "atlas-object", "dry-run")
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
required <- c("atlas-metadata", "output-dir")
if (!all(required %in% names(flags)))
  stop("Required options missing: ", paste(setdiff(required, names(flags)), collapse = ", "))

atlas_path <- normalizePath(flags[["atlas-metadata"]], winslash = "/", mustWork = TRUE)
out_dir <- flags[["output-dir"]]

atlas <- read.csv(atlas_path, stringsAsFactors = FALSE, check.names = FALSE)
needed <- c("cell_id", "patient", "subtype", "celltype_major")
missing_cols <- setdiff(needed, names(atlas))
if (length(missing_cols))
  stop("Atlas metadata is missing required column(s): ", paste(missing_cols, collapse = ", "))

cat("full atlas cells :", nrow(atlas), "\n")
cat("full atlas donors:", length(unique(atlas$patient)), "\n")

# --- Step 1: source-object TNBC annotation --------------------------------
tnbc <- atlas[atlas$subtype == "TNBC", , drop = FALSE]
cat("TNBC subset      :", nrow(tnbc), "cells /", length(unique(tnbc$patient)), "donors\n")

# --- Step 2: pre-existing Cancer Epithelial annotation ---------------------
malig <- tnbc[tnbc$celltype_major == "Cancer Epithelial", , drop = FALSE]
cat("malignant cells  :", nrow(malig), "cells /", length(unique(malig$patient)), "donors\n")

if (nrow(malig) == 0L) stop("Extraction produced zero cells; check the annotation values.")

# Donors in the TNBC subset that contribute no malignant cells (informational).
dropped <- setdiff(unique(tnbc$patient), unique(malig$patient))

if (isTRUE(flags[["dry-run"]])) {
  cat("\n[dry-run] validation only; no outputs written.\n")
  quit(status = 0L)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)

write.table(data.frame(cell_id = malig$cell_id), file.path(out_dir, "cell_ids.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

counts <- as.data.frame(table(patient = malig$patient), stringsAsFactors = FALSE)
colnames(counts) <- c("patient", "malignant_cells")
write.table(counts, file.path(out_dir, "per_patient_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

summary_tbl <- data.frame(
  step = c("full_atlas", "tnbc_subset", "malignant_epithelial"),
  selection = c("all cells",
                "subtype == 'TNBC'",
                "subtype == 'TNBC' & celltype_major == 'Cancer Epithelial'"),
  cells = c(nrow(atlas), nrow(tnbc), nrow(malig)),
  donors = c(length(unique(atlas$patient)), length(unique(tnbc$patient)),
             length(unique(malig$patient))),
  stringsAsFactors = FALSE)
write.table(summary_tbl, file.path(out_dir, "extraction_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

if (length(dropped)) {
  cat("\nTNBC donors contributing no Cancer Epithelial cells (drop out automatically):\n")
  cat(" ", paste(sort(dropped), collapse = ", "), "\n")
}

if (!is.null(flags[["atlas-object"]])) {
  suppressPackageStartupMessages({ library(Seurat); library(SeuratObject) })
  obj_path <- normalizePath(flags[["atlas-object"]], winslash = "/", mustWork = TRUE)
  obj <- readRDS(obj_path)
  keep <- intersect(colnames(obj), malig$cell_id)
  if (!length(keep)) stop("No extracted cell_ids matched the atlas object.")
  sub <- subset(obj, cells = keep)
  saveRDS(sub, file.path(out_dir, "Wu2021_TNBC_malignant.rds"))
  cat("subset object written:", ncol(sub), "cells\n")
}

cat("\noutputs written to:", out_dir, "\n")
