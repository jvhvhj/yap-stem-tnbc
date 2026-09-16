#!/usr/bin/env Rscript
# Purpose: Export cell-cycle reference genes
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 5; orthogonal transfer in Fig. 4.

options(stringsAsFactors = FALSE)

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else getwd()
default_root <- normalizePath(file.path(dirname(script_path), "..", ".."),
                              winslash = "/", mustWork = TRUE)
root <- Sys.getenv("AHIPPO_YAP_ROOT", unset = default_root)
out <- file.path(root, "0716_yan2026_validation", "inputs",
                 "seurat_cell_cycle_genes_updated_2019.tsv")
lib <- file.path(root, "0703_rebuild", "R_library")
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages(library(Seurat))
data("cc.genes.updated.2019", package = "Seurat")

tab <- rbind(
  data.frame(phase = "S", gene = cc.genes.updated.2019$s.genes),
  data.frame(phase = "G2M", gene = cc.genes.updated.2019$g2m.genes)
)
write.table(tab, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote", out, "with", nrow(tab), "rows\n")
