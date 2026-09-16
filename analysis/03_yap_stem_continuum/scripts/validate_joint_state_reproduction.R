#!/usr/bin/env Rscript
# Validate the public Joint-axis / state derivation against the frozen columns.
#
# This script exists so that the released derivation can be checked rather than
# trusted. It reports, against the frozen authoritative metadata:
#
#   n_cells
#   max_abs_joint_difference          (derived YAP_Stem_Score vs frozen)
#   number_of_state_mismatches        (derived YS_tertile vs frozen)
#   High_n, Middle_n, Low_n
#
# Release criterion: number_of_state_mismatches == 0.
#
# It additionally reports every equal-thirds rule tested, so that the residual
# ambiguity in the original implementation is explicit rather than hidden. Several
# standard rules agree exactly on this population because n = 10,836 is exactly
# divisible by 3 and no tie spans either boundary; that is a property of the data,
# not evidence about which code path the original authors ran.

options(stringsAsFactors = FALSE, width = 220)

usage <- function() cat(paste(
  "Usage: Rscript analysis/03_yap_stem_continuum/scripts/validate_joint_state_reproduction.R",
  "  --derived <joint_state.tsv from derive_joint_axis_and_tertiles.R>",
  "  --frozen  <frozen_cell_state_metadata.csv>",
  "  --output-dir <NEW_output_directory>",
  "",
  "The frozen metadata must contain: cell_id, YAP_Stem_Score, YS_tertile",
  "Exit status is 0 when state mismatches are zero, 1 otherwise.",
  sep = "\n"), "\n")

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) { usage(); quit(status = 0L) }
if (length(args) == 0L) { usage(); quit(status = 2L) }

flags <- list()
allowed <- c("derived", "frozen", "output-dir")
i <- 1L
while (i <= length(args)) {
  key <- sub("^--", "", args[i])
  if (!startsWith(args[i], "--") || !key %in% allowed || key %in% names(flags))
    stop("Unknown or duplicate option: ", args[i])
  if (i == length(args) || startsWith(args[i + 1L], "--"))
    stop("Missing value for --", key)
  flags[[key]] <- args[i + 1L]
  i <- i + 2L
}
required <- c("derived", "frozen", "output-dir")
if (!all(required %in% names(flags)))
  stop("Required options missing: ", paste(setdiff(required, names(flags)), collapse = ", "))

derived <- read.delim(normalizePath(flags[["derived"]], winslash = "/", mustWork = TRUE),
                      check.names = FALSE, stringsAsFactors = FALSE)
frozen  <- read.csv(normalizePath(flags[["frozen"]], winslash = "/", mustWork = TRUE),
                    check.names = FALSE, stringsAsFactors = FALSE)

for (col in c("cell_id", "YAP_Stem_Score", "YS_tertile"))
  if (!col %in% names(frozen)) stop("Frozen table missing column: ", col)
for (col in c("cell_id", "YAP_Stem_Score", "YS_tertile"))
  if (!col %in% names(derived)) stop("Derived table missing column: ", col)

# Align on cell_id; order must not be assumed.
m <- merge(derived[, c("cell_id", "YAP_Stem_Score", "YS_tertile")],
           frozen[, c("cell_id", "YAP_Stem_Score", "YS_tertile")],
           by = "cell_id", suffixes = c("_derived", "_frozen"))
if (nrow(m) != nrow(frozen))
  stop("Cell-id join lost rows: derived ", nrow(m), " vs frozen ", nrow(frozen))

fz <- as.character(m$YS_tertile_frozen)
dz <- as.character(m$YS_tertile_derived)

out_dir <- flags[["output-dir"]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)

primary <- data.frame(
  n_cells = nrow(m),
  max_abs_joint_difference = max(abs(m$YAP_Stem_Score_derived - m$YAP_Stem_Score_frozen)),
  number_of_state_mismatches = sum(dz != fz),
  High_n = sum(dz == "High"),
  Middle_n = sum(dz == "Intermediate"),
  Low_n = sum(dz == "Low"),
  stringsAsFactors = FALSE)
write.table(primary, file.path(out_dir, "state_reproduction_validation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- tertile-rule equivalence evidence -------------------------------------
x <- m$YAP_Stem_Score_derived
n <- length(x); s <- sort(x)
rules <- list(
  "cut(quantile type=1)"        = as.character(cut(x, c(-Inf, quantile(x, c(1/3,2/3), type=1), Inf), labels=c("Low","Intermediate","High"))),
  "cut(quantile type=7)"        = as.character(cut(x, c(-Inf, quantile(x, c(1/3,2/3), type=7), Inf), labels=c("Low","Intermediate","High"))),
  "cut(sorted index n/3, 2n/3)" = as.character(cut(x, c(-Inf, s[n/3], s[2*n/3], Inf), labels=c("Low","Intermediate","High"))),
  "rank(ties=first) equal thirds" = rep(c("Low","Intermediate","High"), each = n/3)[rank(x, ties.method = "first")],
  "order() equal thirds"        = { o <- order(x); g <- character(n); g[o[1:(n/3)]] <- "Low"; g[o[(n/3+1):(2*n/3)]] <- "Intermediate"; g[o[(2*n/3+1):n]] <- "High"; g }
)
eq <- do.call(rbind, lapply(names(rules), function(nm) data.frame(
  rule = nm, mismatches_vs_frozen = sum(rules[[nm]] != fz), stringsAsFactors = FALSE)))
write.table(eq, file.path(out_dir, "tertile_rule_equivalence.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("n_cells                    :", primary$n_cells, "\n")
cat("max_abs_joint_difference   :", sprintf("%.3e", primary$max_abs_joint_difference), "\n")
cat("number_of_state_mismatches :", primary$number_of_state_mismatches, "\n")
cat("High / Middle / Low        :", primary$High_n, "/", primary$Middle_n, "/", primary$Low_n, "\n")
cat("\nequal-thirds rules tested:\n"); print(eq)

if (primary$number_of_state_mismatches != 0) {
  cat("\nFAIL: state assignment does not reproduce the frozen column.\n"); quit(status = 1L)
}
cat("\nPASS: state assignment reproduces the frozen column exactly.\n")
