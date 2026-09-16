# Purpose: Compute cohort recurrence and Lee L neighborhood statistics
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 6 and Supplementary Fig. S4.
options(stringsAsFactors = FALSE, warn = 1)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
tmp_out <- file.path(root, ".tmp", "Figure6_enrichment_gate", "analysis_out")
dir.create(tmp_out, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, quote = "", comment.char = "", na.strings = c("NA", ""))
}

write_tsv <- function(x, name) {
  write.table(x, file.path(tmp_out, name), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE, na = "NA")
}

num_id <- function(x) as.integer(sub("^P", "", x))

fmt <- function(x, digits = 6) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "fg", flag = "#"))
}

safe_wilcox <- function(x, y = NULL, alternative = "two.sided", paired = FALSE) {
  exact_ok <- if (is.null(y)) {
    all(is.finite(x)) && !any(x == 0) && !anyDuplicated(abs(x))
  } else {
    all(is.finite(c(x, y))) && !anyDuplicated(c(x, y))
  }
  obj <- suppressWarnings(wilcox.test(x, y, alternative = alternative,
                                      paired = paired, exact = exact_ok,
                                      conf.int = FALSE, correct = !exact_ok))
  list(statistic = unname(obj$statistic), p = obj$p.value,
       exact = exact_ok, method = obj$method)
}

# -----------------------------------------------------------------------------
# Workstream D: author-defined Reference / Validation recurrence
# -----------------------------------------------------------------------------
d_source <- file.path(root, "Figure6_S1_GSE210616_patient_recurrent_spatial_coupling",
                      "Figure6_S1_patient_associations.tsv")
registry_source <- file.path(root, "Figure6_S0Q2_GSE210616_primary_spatial_candidate_integrity_audit",
                             "Figure6_S0Q2_GSE210616_section_registry.tsv")
d <- read_tsv(d_source)
reg <- read_tsv(registry_source)

map <- unique(reg[, c("patient_id", "cohort")])
map$n_sections <- as.integer(table(reg$patient_id)[map$patient_id])
map <- map[order(num_id(map$patient_id)), ]
stopifnot(nrow(map) == 22L,
          length(unique(map$patient_id[map$cohort == "Reference"])) == 14L,
          length(unique(map$patient_id[map$cohort == "Validation"])) == 8L,
          sum(map$n_sections[map$cohort == "Reference"]) == 28L,
          sum(map$n_sections[map$cohort == "Validation"]) == 15L)

map$mapping_source <- registry_source
map$mapping_status <- "AUTHOR_DEFINED_MAPPING_CONFIRMED"
write_tsv(map, "Figure6D_cohort_mapping.tsv")

d2 <- merge(d, map[, c("patient_id", "cohort")], by = "patient_id", all.x = TRUE)
stopifnot(nrow(d2) == 22L, !anyNA(d2$cohort), !anyNA(d2$rho_tumour_depth_adjusted))

cohort_rows <- lapply(c("Reference", "Validation"), function(cc) {
  x <- d2$rho_tumour_depth_adjusted[d2$cohort == cc]
  wt <- safe_wilcox(x, alternative = "greater")
  data.frame(
    cohort = cc,
    n_patients = length(x),
    n_sections = sum(map$n_sections[map$cohort == cc]),
    median_rho = median(x),
    q1_rho = unname(quantile(x, 0.25, type = 7)),
    q3_rho = unname(quantile(x, 0.75, type = 7)),
    positive_n = sum(x > 0),
    positive_fraction = mean(x > 0),
    wilcoxon_statistic = wt$statistic,
    raw_p_one_sided = wt$p,
    exact_p_used = wt$exact,
    test_method = wt$method,
    stringsAsFactors = FALSE
  )
})
cohort_results <- do.call(rbind, cohort_rows)
cohort_results$holm_p_two_cohorts <- p.adjust(cohort_results$raw_p_one_sided, method = "holm")
cohort_results$source_statistic <- "rho_tumour_depth_adjusted"
cohort_results$source_file <- d_source
cohort_results <- cohort_results[, c("cohort", "n_patients", "n_sections", "median_rho",
                                     "q1_rho", "q3_rho", "positive_n", "positive_fraction",
                                     "wilcoxon_statistic", "raw_p_one_sided", "holm_p_two_cohorts",
                                     "exact_p_used", "test_method", "source_statistic", "source_file")]
write_tsv(cohort_results, "Figure6D_reference_validation_results.tsv")

ref <- d2$rho_tumour_depth_adjusted[d2$cohort == "Reference"]
val <- d2$rho_tumour_depth_adjusted[d2$cohort == "Validation"]
between <- safe_wilcox(ref, val, alternative = "two.sided")
boundary <- data.frame(
  comparison = "Reference_vs_Validation",
  test = "two-sided Wilcoxon rank-sum (descriptive)",
  reference_n = length(ref),
  validation_n = length(val),
  reference_median = median(ref),
  validation_median = median(val),
  median_difference_validation_minus_reference = median(val) - median(ref),
  wilcoxon_statistic = between$statistic,
  p_two_sided = between$p,
  exact_p_used = between$exact,
  inferential_role = "DESCRIPTIVE_BOUNDARY_ONLY",
  interpretation_boundary = "Not a formal interaction test; does not establish equivalence or difference in replication strength.",
  stringsAsFactors = FALSE
)
write_tsv(boundary, "Figure6D_between_cohort_boundary.tsv")

d_strong <- all(cohort_results$median_rho > 0 & cohort_results$positive_fraction >= 0.80 & cohort_results$holm_p_two_cohorts < 0.01)
d_pass <- all(cohort_results$median_rho > 0 & cohort_results$positive_fraction >= 0.70 & cohort_results$holm_p_two_cohorts < 0.05)
d_border <- all(cohort_results$median_rho > 0) && any(cohort_results$positive_fraction >= 0.70 | cohort_results$holm_p_two_cohorts < 0.05)
d_verdict <- if (d_strong) "FIGURE6D_REFERENCE_VALIDATION_STRONG_PASS" else if (d_pass) "FIGURE6D_REFERENCE_VALIDATION_PASS" else if (d_border) "FIGURE6D_REFERENCE_VALIDATION_BORDERLINE" else "FIGURE6D_REFERENCE_VALIDATION_FAIL"

writeLines(c(
  "# Figure 6D reference/validation replication readout",
  "",
  paste0("Verdict: **", d_verdict, "**"),
  "",
  "The frozen patient-level `rho_tumour_depth_adjusted` values were mapped to the author-defined Reference and Validation cohorts without recomputing spot scores or patient effects.",
  "",
  paste0("- Reference: n=", cohort_results$n_patients[1], ", median rho=", fmt(cohort_results$median_rho[1], 4),
         ", positive=", cohort_results$positive_n[1], "/", cohort_results$n_patients[1],
         ", one-sided Wilcoxon P=", fmt(cohort_results$raw_p_one_sided[1], 4),
         ", Holm P=", fmt(cohort_results$holm_p_two_cohorts[1], 4), "."),
  paste0("- Validation: n=", cohort_results$n_patients[2], ", median rho=", fmt(cohort_results$median_rho[2], 4),
         ", positive=", cohort_results$positive_n[2], "/", cohort_results$n_patients[2],
         ", one-sided Wilcoxon P=", fmt(cohort_results$raw_p_one_sided[2], 4),
         ", Holm P=", fmt(cohort_results$holm_p_two_cohorts[2], 4), "."),
  paste0("- Descriptive between-cohort rank-sum P=", fmt(boundary$p_two_sided, 4), "."),
  "",
  "The between-cohort test is descriptive only and is not used to rescue either within-cohort replication result."
), file.path(tmp_out, "Figure6D_readout.md"))

# -----------------------------------------------------------------------------
# Workstream E: first-order Visium adjacency and symmetric Lee's L
# -----------------------------------------------------------------------------
score_dir <- ".runtime/Figure6_S1_chunks"
geo_dir <- ".runtime/Figure6_S0Q2_geo"

lee_l <- function(x, y, W, n_nonisolated) {
  xc <- x - mean(x)
  yc <- y - mean(y)
  lx <- as.numeric(W %*% xc)
  ly <- as.numeric(W %*% yc)
  denom <- sqrt(sum(xc^2) * sum(yc^2))
  if (!is.finite(denom) || denom <= 0 || n_nonisolated <= 0) return(NA_real_)
  length(x) / n_nonisolated * sum(lx * ly) / denom
}

graph_rows <- vector("list", nrow(reg))
section_rows <- vector("list", nrow(reg))
perm_rows <- vector("list", nrow(reg))

for (ii in seq_len(nrow(reg))) {
  rr <- reg[ii, ]
  score_path <- file.path(score_dir, paste0(rr$sample_id, "_", rr$section_id, "_scores.tsv"))
  pos_path <- file.path(geo_dir, rr$sample_id,
                        paste0(rr$sample_id, "_", rr$section_id, "_tissue_positions_list.csv.gz"))
  if (!file.exists(score_path) || !file.exists(pos_path)) stop("Missing section source: ", rr$sample_id, "/", rr$section_id)
  sc <- read_tsv(score_path)
  pp <- read.csv(gzfile(pos_path), header = FALSE,
                 col.names = c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"))
  pp <- pp[pp$in_tissue == 1, ]
  dat <- merge(sc[, c("barcode", "Joint_spatial", "Program146_spatial")],
               pp[, c("barcode", "array_row", "array_col")], by = "barcode", all.x = TRUE, sort = FALSE)
  stopifnot(nrow(dat) == nrow(sc), !anyNA(dat$array_row), !anyNA(dat$array_col),
            !anyNA(dat$Joint_spatial), !anyNA(dat$Program146_spatial))
  key <- paste(dat$array_row, dat$array_col, sep = "_")
  stopifnot(!anyDuplicated(key))
  lookup <- setNames(seq_len(nrow(dat)), key)
  # Three forward first-order offsets generate each undirected hex-grid edge once.
  offsets <- matrix(c(0, 2, 1, 1, 1, -1), byrow = TRUE, ncol = 2)
  edge_i <- integer(0); edge_j <- integer(0)
  for (kk in seq_len(nrow(offsets))) {
    target <- paste(dat$array_row + offsets[kk, 1], dat$array_col + offsets[kk, 2], sep = "_")
    jj <- unname(lookup[target])
    ok <- !is.na(jj)
    edge_i <- c(edge_i, which(ok))
    edge_j <- c(edge_j, jj[ok])
  }
  n <- nrow(dat)
  deg <- tabulate(c(edge_i, edge_j), nbins = n)
  n_edges <- length(edge_i)
  stopifnot(n_edges > 0, max(deg) <= 6)
  # Binary symmetric adjacency, row-standardized for Lee's L.
  adj_i <- c(edge_i, edge_j)
  adj_j <- c(edge_j, edge_i)
  w <- 1 / deg[adj_i]
  W <- Matrix::sparseMatrix(i = adj_i, j = adj_j, x = w, dims = c(n, n))
  n_nonisolated <- sum(deg > 0)
  row_sum_error <- if (n_nonisolated > 0) max(abs(Matrix::rowSums(W)[deg > 0] - 1)) else NA_real_
  obs <- lee_l(dat$Joint_spatial, dat$Program146_spatial, W, n_nonisolated)
  symmetry_check <- lee_l(dat$Program146_spatial, dat$Joint_spatial, W, n_nonisolated)
  if (!isTRUE(all.equal(obs, symmetry_check, tolerance = 1e-12))) stop("Lee L symmetry failure: ", rr$section_id)
  seed <- 20260825L + ii
  set.seed(seed)
  null <- numeric(1000)
  for (bb in seq_len(1000)) {
    null[bb] <- lee_l(dat$Joint_spatial, sample(dat$Program146_spatial, replace = FALSE), W, n_nonisolated)
  }
  emp_p <- (1 + sum(null >= obs)) / 1001
  graph_rows[[ii]] <- data.frame(
    patient_id = rr$patient_id, cohort = rr$cohort, sample_id = rr$sample_id,
    section_id = rr$section_id, n_spots = n, n_edges = n_edges,
    median_degree = median(deg), min_degree = min(deg), max_degree = max(deg),
    isolated_spots = sum(deg == 0), row_standardization_max_error = row_sum_error,
    coordinate_join_complete = TRUE, hex_geometry_verified = max(deg) <= 6,
    adjacency_definition = "first-order Visium array neighbors: (0,+/-2),(+/-1,+/-1)",
    score_source = score_path, coordinate_source = pos_path,
    stringsAsFactors = FALSE
  )
  section_rows[[ii]] <- data.frame(
    patient_id = rr$patient_id, cohort = rr$cohort, sample_id = rr$sample_id,
    section_id = rr$section_id, n_spots = n, n_edges = n_edges,
    observed_Lee_L = obs, null_median_Lee_L = median(null),
    empirical_p_one_sided = emp_p, permutations = 1000L,
    permutation_seed = seed, symmetry_check_absolute_difference = abs(obs - symmetry_check),
    section_supportive_p_lt_0_05 = emp_p < 0.05,
    stringsAsFactors = FALSE
  )
  perm_rows[[ii]] <- data.frame(
    patient_id = rr$patient_id, cohort = rr$cohort, sample_id = rr$sample_id,
    section_id = rr$section_id, permutation_id = seq_len(1000),
    Lee_L_null = null, permutation_seed = seed,
    stringsAsFactors = FALSE
  )
  message("E section ", ii, "/", nrow(reg), ": ", rr$patient_id, " ", rr$section_id,
          " Lee L=", signif(obs, 4), " P=", signif(emp_p, 4))
}

graph_qc <- do.call(rbind, graph_rows)
section_lee <- do.call(rbind, section_rows)
perms <- do.call(rbind, perm_rows)
graph_qc <- graph_qc[order(num_id(graph_qc$patient_id), graph_qc$section_id), ]
section_lee <- section_lee[order(num_id(section_lee$patient_id), section_lee$section_id), ]
perms <- perms[order(num_id(perms$patient_id), perms$section_id, perms$permutation_id), ]
write_tsv(graph_qc, "Figure6E_spatial_graph_QC.tsv")
write_tsv(section_lee, "Figure6E_section_LeeL.tsv")
write_tsv(perms, "Figure6E_permutation_results.tsv")

patient_rows <- lapply(split(section_lee, section_lee$patient_id), function(z) {
  data.frame(
    patient_id = z$patient_id[1],
    cohort = z$cohort[1],
    n_sections = nrow(z),
    total_edges = sum(z$n_edges),
    patient_Lee_L = weighted.mean(z$observed_Lee_L, w = z$n_edges),
    positive = weighted.mean(z$observed_Lee_L, w = z$n_edges) > 0,
    aggregation = ifelse(nrow(z) == 1, "single section", "edge-count-weighted section mean"),
    stringsAsFactors = FALSE
  )
})
patient_lee <- do.call(rbind, patient_rows)
patient_lee <- patient_lee[order(num_id(patient_lee$patient_id)), ]
write_tsv(patient_lee, "Figure6E_patient_LeeL.tsv")

e_test <- safe_wilcox(patient_lee$patient_Lee_L, alternative = "greater")
e_positive_fraction <- mean(patient_lee$patient_Lee_L > 0)
e_strong <- median(patient_lee$patient_Lee_L) > 0 && e_positive_fraction >= 0.80 && e_test$p < 0.01
e_pass <- median(patient_lee$patient_Lee_L) > 0 && e_positive_fraction >= 0.70 && e_test$p < 0.05
e_border <- median(patient_lee$patient_Lee_L) > 0 && (e_positive_fraction >= 0.70 || e_test$p < 0.05)
e_verdict <- if (e_strong) "FIGURE6E_NEIGHBORHOOD_SPATIAL_STRONG_PASS" else if (e_pass) "FIGURE6E_NEIGHBORHOOD_SPATIAL_PASS" else if (e_border) "FIGURE6E_NEIGHBORHOOD_SPATIAL_BORDERLINE" else "FIGURE6E_NEIGHBORHOOD_SPATIAL_FAIL"
primary <- data.frame(
  analysis_unit = "patient",
  n_patients = nrow(patient_lee),
  n_sections = nrow(section_lee),
  median_Lee_L = median(patient_lee$patient_Lee_L),
  q1_Lee_L = unname(quantile(patient_lee$patient_Lee_L, 0.25, type = 7)),
  q3_Lee_L = unname(quantile(patient_lee$patient_Lee_L, 0.75, type = 7)),
  positive_n = sum(patient_lee$patient_Lee_L > 0),
  positive_fraction = e_positive_fraction,
  wilcoxon_statistic = e_test$statistic,
  p_one_sided = e_test$p,
  exact_p_used = e_test$exact,
  test_method = e_test$method,
  verdict = e_verdict,
  stringsAsFactors = FALSE
)
write_tsv(primary, "Figure6E_primary_test.tsv")

writeLines(c(
  "# Figure 6E neighbourhood spatial co-organization readout",
  "",
  paste0("Verdict: **", e_verdict, "**"),
  "",
  "A direct first-order Visium hex-grid graph was constructed independently for each of 43 sections using official array coordinates. No expression-neighbour graph, adaptive radius, smoothing, score threshold or spot-effect recomputation was used.",
  "",
  paste0("- Graph QC: ", sum(graph_qc$coordinate_join_complete & graph_qc$hex_geometry_verified), "/", nrow(graph_qc), " sections passed coordinate and geometry checks; total undirected edges=", sum(graph_qc$n_edges), "."),
  paste0("- Section permutation support (one-sided empirical P<0.05): ", sum(section_lee$section_supportive_p_lt_0_05), "/", nrow(section_lee), "."),
  paste0("- Patient-level Lee's L: median=", fmt(primary$median_Lee_L, 4), " (IQR ", fmt(primary$q1_Lee_L, 4), " to ", fmt(primary$q3_Lee_L, 4), "); positive=", primary$positive_n, "/", primary$n_patients, "."),
  paste0("- One-sided patient-level Wilcoxon signed-rank P=", fmt(primary$p_one_sided, 4), "."),
  "",
  "Lee's L was computed with row-standardized first-order adjacency using the symmetric Lee (2001) cross-product formula. Patient values are edge-count-weighted section means when paired sections exist."
), file.path(tmp_out, "Figure6E_readout.md"))

writeLines(c(d_verdict, e_verdict), file.path(tmp_out, "DE_status.txt"))
