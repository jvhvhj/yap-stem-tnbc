#!/usr/bin/env Rscript
# Purpose: Read existing normalized expression by exact cell ID
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 2 and Supplementary Fig. S1.

# Read-only recovery of the already normalized MCL1 RNA/data layer.
# This script does not run NormalizeData, rescale expression, or alter the RDS.

suppressPackageStartupMessages(library(Matrix))

root <- Sys.getenv("AHIPPO_ROOT", unset = ".")
source_rds <- file.path(root, "output_step0", "Wu2021_TNBC_malignant.rds")
output_tsv <- file.path(root, "0911_FigureS1_preplot_gate", "tables", "MCL1_recovered_normalized_expression.tsv.gz")

stopifnot(file.exists(source_rds))
dir.create(dirname(output_tsv), recursive = TRUE, showWarnings = FALSE)

object <- readRDS(source_rds)
object_attributes <- attributes(object)
stopifnot("assays" %in% names(object_attributes), "meta.data" %in% names(object_attributes))

rna <- object_attributes$assays[["RNA"]]
rna_attributes <- attributes(rna)
stopifnot("layers" %in% names(rna_attributes), "features" %in% names(rna_attributes), "cells" %in% names(rna_attributes))
stopifnot("data" %in% names(rna_attributes$layers))

data_layer <- rna_attributes$layers[["data"]]
feature_map <- rna_attributes$features
cell_map <- rna_attributes$cells
feature_names <- attr(feature_map, "dimnames")[[1]]
cell_names <- attr(cell_map, "dimnames")[[1]]

stopifnot(length(feature_names) == nrow(data_layer))
stopifnot(length(cell_names) == ncol(data_layer))
stopifnot("MCL1" %in% feature_names)
stopifnot(length(cell_names) == 10836L)
stopifnot(!anyDuplicated(cell_names))

mcl1_index <- match("MCL1", feature_names)
mcl1_value <- as.numeric(data_layer[mcl1_index, ])
stopifnot(length(mcl1_value) == 10836L, !anyNA(mcl1_value), length(unique(mcl1_value)) > 1L)

result <- data.frame(
  cell_id = cell_names,
  patient = sub("_.*$", "", cell_names),
  MCL1 = mcl1_value,
  source_assay = "RNA",
  source_layer = "data",
  extraction = "existing_normalized_value_no_recomputation",
  stringsAsFactors = FALSE
)

connection <- gzfile(output_tsv, open = "wt")
write.table(result, file = connection, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
close(connection)

cat("Recovered existing MCL1 RNA/data values:", nrow(result), "cells\n")
cat("Non-zero values:", sum(result$MCL1 != 0), "\n")
cat("Source:", source_rds, "\n")
cat("Output:", output_tsv, "\n")
