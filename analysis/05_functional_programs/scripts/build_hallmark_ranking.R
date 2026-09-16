#!/usr/bin/env Rscript
# Build the Figure 4E/F Hallmark GSEA ranking from the authoritative gene-level
# DESeq2 table.
#
# This closes the previously unresolved upstream link in the Figure 4E/F chain.
# It is a deterministic PROJECTION of an existing result table, not a
# re-derivation: no model is fitted, no gene is filtered, and no statistic is
# recomputed here.
#
# Chain:
#   Figure3_I_effect_consistency_source.tsv        (gene-level DESeq2 Wald table,
#                                                   17,597 score-excluded genes;
#                                                   produced by the programme
#                                                   derivation module)
#     -> this script -> Figure4_E_GSEA_v2_ranking.tsv      (ranking input)
#     -> external_hallmark_analysis.R /
#        hallmark_gsea_and_leading_edge.R                  (fgseaMultilevel)
#     -> Figure4_E_GSEA_v2_results.tsv                     (50 Hallmark pathways)
#     -> Figure 4E / 4F
#
# Ranking rule (fixed before the enrichment result was examined):
#   DESeq2 Wald statistic for state High versus Other, design ~ patient + state.
#   Positive values mean High-enriched.
#
# Score-definition exclusion: all 38 YAP17/Stem21 genes were removed before gene
# filtering and testing, so the testable universe is 17,597 score-independent
# genes. The source table carries `score_definition_genes_excluded_before_testing`
# and this script asserts it is TRUE for every retained row rather than trusting it.

options(stringsAsFactors = FALSE, width = 220)

usage <- function() cat(paste(
  "Usage: Rscript analysis/05_functional_programs/scripts/build_hallmark_ranking.R",
  "  --source-table <Figure3_I_effect_consistency_source.tsv>",
  "  --output-dir <NEW_output_directory>",
  "  [--expected-genes 17597]",
  "  [--dry-run]",
  "",
  "Writes Figure4_E_GSEA_v2_ranking.tsv with columns:",
  "  gene, ranking_statistic, source_effect, source_FDR, testable_status",
  sep = "\n"), "\n")

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) { usage(); quit(status = 0L) }
if (length(args) == 0L) { usage(); quit(status = 2L) }

flags <- list()
allowed <- c("source-table", "output-dir", "expected-genes", "dry-run")
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
required <- c("source-table", "output-dir")
if (!all(required %in% names(flags)))
  stop("Required options missing: ", paste(setdiff(required, names(flags)), collapse = ", "))

src_path <- normalizePath(flags[["source-table"]], winslash = "/", mustWork = TRUE)
src <- read.delim(src_path, check.names = FALSE, stringsAsFactors = FALSE)

need <- c("gene", "Wald_stat", "overall_log2FC_High_vs_Other", "BH_FDR",
          "score_definition_genes_excluded_before_testing")
if (length(setdiff(need, names(src))))
  stop("Source table missing column(s): ", paste(setdiff(need, names(src)), collapse = ", "))

cat("source genes:", nrow(src), "\n")

# The score-definition exclusion must already have been applied upstream. Assert
# rather than assume: a FALSE or missing value here would silently change the
# tested universe.
if (anyNA(src$score_definition_genes_excluded_before_testing) ||
    !all(as.logical(src$score_definition_genes_excluded_before_testing)))
  stop("Source table contains genes that were not score-definition excluded.")

if (anyNA(src$Wald_stat)) stop("Source table contains missing Wald statistics.")
if (anyDuplicated(src$gene)) stop("Source table contains duplicated gene symbols.")

expected <- if (!is.null(flags[["expected-genes"]])) as.integer(flags[["expected-genes"]]) else 17597L
if (nrow(src) != expected)
  stop("Expected ", expected, " score-independent testable genes, found ", nrow(src), ".")

if (isTRUE(flags[["dry-run"]])) {
  cat("[dry-run] validation only; no outputs written.\n")
  quit(status = 0L)
}

ranking <- data.frame(
  gene = src$gene,
  ranking_statistic = src$Wald_stat,
  source_effect = src$overall_log2FC_High_vs_Other,
  source_FDR = src$BH_FDR,
  testable_status = "TESTABLE_SCORE_INDEPENDENT",
  stringsAsFactors = FALSE)

# Descending Wald statistic defines the GSEA ranking order.
ranking <- ranking[order(-ranking$ranking_statistic), , drop = FALSE]

out_dir <- flags[["output-dir"]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)

write.table(ranking, file.path(out_dir, "Figure4_E_GSEA_v2_ranking.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("ranked genes written:", nrow(ranking), "\n")
cat("top gene:", ranking$gene[1], "(", ranking$ranking_statistic[1], ")\n")
cat("outputs written to:", out_dir, "\n")
