#!/usr/bin/env Rscript
# Parse scripts only. Does not source, evaluate or run any analysis.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript parse_R_scripts.R REPOSITORY_ROOT")
root <- args[1]
if (!dir.exists(root)) stop("Repository directory does not exist")
files <- list.files(file.path(root, "analysis"), pattern = "\\.[Rr]$", full.names = TRUE, recursive = TRUE)
if (!length(files)) stop("No R analysis scripts found; directory/encoding check required")
for (p in files) parse(file = p, keep.source = FALSE, encoding = "UTF-8")
cat("R syntax:", length(files), "files. No script evaluated.\n")
