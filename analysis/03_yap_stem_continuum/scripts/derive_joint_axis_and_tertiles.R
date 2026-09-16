#!/usr/bin/env Rscript
# Joint YAP-Stem axis and Low/Intermediate/High state assignment (Wu discovery).
#
# Recovered and verification-checked construction:
#
#   1. Standardise each component score across the frozen malignant population
#      using the SAMPLE standard deviation (R scale(), n - 1 denominator).
#   2. Joint = z(YAP17 mean) + z(Stem21 mean).
#   3. Split the Joint value into equal thirds of the frozen population.
#
# The cells used for standardisation are exactly the frozen malignant epithelial
# population from extract_wu_malignant_population.R (10,836 cells / 8 donors).
# Standardisation is cohort-level, not per-patient.
#
# On the tertile rule
# -------------------
# The cut points are the sorted Joint values at the one-third and two-thirds
# positions, i.e. quantile(x, c(1/3, 2/3), type = 1). For the Wu population
# n = 10,836 is exactly divisible by 3 and no tie spans either boundary, so every
# equal-thirds rule that places the breaks between the same two adjacent sorted
# values yields the identical assignment. This was verified against the frozen
# YS_tertile column for cut(quantile type=1), cut(quantile type=7), explicit
# sorted-index breaks, rank-based splitting and order-based splitting: all five
# produce zero mismatches on all 10,836 cells. The rule is therefore reproducible
# even though the original implementation's exact code path is not uniquely
# determined by the data. See validate_joint_state_reproduction.R.
#
# Expected frozen outcome: High = 3,612, Other = 7,224.
#
# All 38 YAP17/Stem21 genes are excluded from Program146 by construction, so the
# Joint axis is score-defining and must not be used as an outcome of itself.

options(stringsAsFactors = FALSE, width = 220)

usage <- function() cat(paste(
  "Usage: Rscript analysis/03_yap_stem_continuum/scripts/derive_joint_axis_and_tertiles.R",
  "  --component-scores <component_scores.tsv from score_yap17_stem21.R>",
  "  --output-dir <NEW_output_directory>",
  "  [--dry-run]",
  "",
  "Writes joint_state.tsv with one row per cell:",
  "  cell_id, patient, YAP_score, Stemness_score, z_YAP, z_Stem,",
  "  YAP_Stem_Score, YS_tertile",
  "Dry-run validates the input table only.",
  sep = "\n"), "\n")

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) { usage(); quit(status = 0L) }
if (length(args) == 0L) { usage(); quit(status = 2L) }

flags <- list()
allowed <- c("component-scores", "output-dir", "dry-run")
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
required <- c("component-scores", "output-dir")
if (!all(required %in% names(flags)))
  stop("Required options missing: ", paste(setdiff(required, names(flags)), collapse = ", "))

cs_path <- normalizePath(flags[["component-scores"]], winslash = "/", mustWork = TRUE)
cs <- read.delim(cs_path, check.names = FALSE, stringsAsFactors = FALSE)
need <- c("cell_id", "YAP_score", "Stemness_score")
if (length(setdiff(need, names(cs))))
  stop("component_scores table is missing: ", paste(setdiff(need, names(cs)), collapse = ", "))
if (anyNA(cs$YAP_score) || anyNA(cs$Stemness_score))
  stop("Component scores contain missing values.")

cat("cells:", nrow(cs), "\n")

if (isTRUE(flags[["dry-run"]])) {
  cat("[dry-run] validation only; no outputs written.\n")
  quit(status = 0L)
}

# --- 1. cohort-level z standardisation (sample SD) -------------------------
z_yap  <- as.numeric(scale(cs$YAP_score))
z_stem <- as.numeric(scale(cs$Stemness_score))

# --- 2. Joint axis ---------------------------------------------------------
joint <- z_yap + z_stem

# --- 3. equal-thirds state assignment --------------------------------------
breaks <- quantile(joint, c(1/3, 2/3), type = 1)
state <- as.character(cut(joint, breaks = c(-Inf, breaks, Inf),
                          labels = c("Low", "Intermediate", "High")))
if (anyNA(state)) stop("State assignment produced missing values.")

out_dir <- flags[["output-dir"]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)

res <- data.frame(
  cell_id = cs$cell_id,
  patient = if ("patient" %in% names(cs)) cs$patient else NA_character_,
  YAP_score = cs$YAP_score,
  Stemness_score = cs$Stemness_score,
  z_YAP = z_yap,
  z_Stem = z_stem,
  YAP_Stem_Score = joint,
  YS_tertile = state,
  stringsAsFactors = FALSE)

write.table(res, file.path(out_dir, "joint_state.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(data.frame(
  cutpoint_low_intermediate = breaks[[1]],
  cutpoint_intermediate_high = breaks[[2]],
  n_cells = nrow(res),
  Low = sum(state == "Low"),
  Intermediate = sum(state == "Intermediate"),
  High = sum(state == "High"),
  Other = sum(state != "High"),
  stringsAsFactors = FALSE),
  file.path(out_dir, "state_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)

cat("Low / Intermediate / High:", sum(state == "Low"), "/",
    sum(state == "Intermediate"), "/", sum(state == "High"), "\n")
cat("High =", sum(state == "High"), " Other =", sum(state != "High"), "\n")
cat("outputs written to:", out_dir, "\n")
