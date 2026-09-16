# Purpose: Figure 5 final data closure
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 5; orthogonal transfer in Fig. 4.
options(stringsAsFactors = FALSE, scipen = 999)
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

root <- "."
interim <- file.path(root, "Figure5_final_execution_protocol")
working <- file.path(root, "Figure5_working")
final <- file.path(root, "Figure5_final_closure")
source_out <- file.path(final, "source_data")
report_out <- file.path(final, "reports")
dir.create(working, recursive = TRUE, showWarnings = FALSE)
dir.create(source_out, recursive = TRUE, showWarnings = FALSE)
dir.create(report_out, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) fread(path, sep = "\t", na.strings = c("NA", "NaN", ""))
write_tsv_both <- function(x, working_name, final_name = working_name) {
  fwrite(as.data.table(x), file.path(working, working_name), sep = "\t", quote = FALSE, na = "NA")
  fwrite(as.data.table(x), file.path(source_out, final_name), sep = "\t", quote = FALSE, na = "NA")
}

paths <- c(
  provenance = file.path(root, "0719_frozen_score_definition_audit/exact_signature_provenance.tsv"),
  patient_rho = file.path(root, "0716_yan2026_validation/tables/yan_patient_yap_stem_rho.tsv"),
  patient_meta = file.path(root, "0716_yan2026_validation/tables/yan_yap_stem_meta_analysis.tsv"),
  module_patient = file.path(root, "0716_yan2026_validation/tables/yan_module_effects_by_patient.tsv"),
  module_robustness = file.path(root, "0716_yan2026_validation/tables/yan_module_robustness.tsv"),
  module_sd = file.path(root, "0720_module_scale_harmonization_and_affected_panel_rebuild/tables/yan_module_cohort_sd.tsv"),
  module_harmonized_summary = file.path(root, "0720_module_scale_harmonization_and_affected_panel_rebuild/tables/plotting_data/Figure5D_harmonized_module_replication.tsv"),
  cache = file.path(root, "0716_yan2026_validation/inputs/yan_cancer_cell_scores.tsv.gz"),
  metaprogram = file.path(root, "0716_yan2026_validation/tables/yan_metaprogram_alignment.tsv"),
  strict_source = file.path(root, "Figure3_bottom_panels_HI_feasibility_and_plotting/Figure3_I_effect_consistency_source.tsv"),
  direction_gate = file.path(root, "0728_sensitivity_and_cnv_closure/02_direction_gate_sensitivity/gate_gene_Wu_Yan_direction.tsv"),
  direction_report = file.path(root, "0728_sensitivity_and_cnv_closure/02_direction_gate_sensitivity/DIRECTION_GATE_REPORT.md"),
  program_patient = file.path(root, "0728_sensitivity_and_cnv_closure/06_Yan_sensitivity/Yan_primary146_state_definition_sensitivity.tsv"),
  gse_transfer = file.path(root, "Figure4F_GSE180286_program_transfer/Figure4F_program_gene_transfer.tsv"),
  gse_provenance = file.path(root, "Figure4F_GSE180286_program_transfer/GSE180286_sample_patient_provenance_audit.tsv"),
  interim_gene = file.path(interim, "Fig5_strictv2_gene_conservation.tsv"),
  null = file.path(root, "0717_methodological_framework_benchmark/tables/random_signature_null.tsv"),
  response = file.path(root, "0716_yan2026_validation/tables/yan_response_association.tsv")
)
stopifnot(all(file.exists(paths)))

# -----------------------------------------------------------------------------
# Progressive association: exact stored patient effects + exact stored DL meta.
# -----------------------------------------------------------------------------
rho <- read_tsv(paths[["patient_rho"]])[included_n50 == TRUE]
meta <- read_tsv(paths[["patient_meta"]])[threshold == 50]
model_order <- c(raw = 1L, technical_primary = 2L, technical_extended = 3L)
model_label <- c(
  raw = "Raw",
  technical_primary = "Technical-adjusted",
  technical_extended = "Technical + Hypoxia/UPR"
)
covariate_exact <- c(
  raw = "none",
  technical_primary = "log1p(nCount_RNA) + nFeature_RNA + S_score + G2M_score",
  technical_extended = "log1p(nCount_RNA) + nFeature_RNA + S_score + G2M_score + Hypoxia + UPR"
)
rho[, `:=`(
  record_type = "patient",
  model_order = unname(model_order[model]),
  model_label = unname(model_label[model]),
  exact_covariates = unname(covariate_exact[model]),
  effect_direction = fifelse(rho > 0, "positive", fifelse(rho < 0, "negative", "zero")),
  biological_replicate = "patient",
  primary_threshold = ">=50 malignant cells"
)]
meta[, `:=`(
  model_order = unname(model_order[model]),
  model_label = unname(model_label[model]),
  exact_covariates = unname(covariate_exact[model]),
  meta_method = "Fisher-z DerSimonian-Laird random effects (exact method in authoritative 0716 script)",
  biological_replicate = "patient",
  primary_threshold = ">=50 malignant cells"
)]
progressive <- rbindlist(list(rho, meta), fill = TRUE)
setorder(progressive, record_type, model_order, patient, omitted_patient)
write_tsv_both(progressive, "Fig5_progressive_adjustment.tsv")

# -----------------------------------------------------------------------------
# Eight frozen-module patient landscape on the 0720 harmonized cohort-z scale.
# Contrast remains High-Low, exactly as the frozen module analysis.
# -----------------------------------------------------------------------------
modules <- c("UPR", "TNFA NFKB", "Hypoxia", "Adhesion Remodeling", "Wound Healing",
             "Survival Stress", "Anoikis Resistance", "Integrin Adhesion")
module_order <- setNames(seq_along(modules), modules)
mp <- read_tsv(paths[["module_patient"]])[
  version == "score_independent" & module %in% modules & testable == TRUE
]
sdtab <- read_tsv(paths[["module_sd"]])[
  version == "score_independent" & module %in% modules,
  .(module, cohort_raw_module_score_sd)
]
mp <- merge(mp, sdtab, by = "module", all.x = TRUE, validate = "many-to-one")
stopifnot(!anyNA(mp$cohort_raw_module_score_sd))
mp[, `:=`(
  module_order = unname(module_order[module]),
  effect_z = high_minus_low / cohort_raw_module_score_sd,
  module_class = fifelse(module == "Integrin Adhesion", "supportive_boundary", "core"),
  state_contrast = "High minus Low",
  score_scale = "cohort-wide z(module raw score) across all 49,275 Yan cancer cells",
  biological_replicate = "patient"
)]

cache <- read_tsv(paths[["cache"]])
eligible <- cache[, .N, by = patient][N >= 50, patient]
cache50 <- cache[patient %in% eligible]
cache50[, metaprogram := sub("__.*$", "", cell_state)]
cache50[!grepl("^M[0-9]{2}$", metaprogram), metaprogram := NA_character_]
patient_annotation <- cache50[, {
  tab <- sort(table(metaprogram, useNA = "no"), decreasing = TRUE)
  dominant <- if (length(tab)) sort(names(tab)[tab == max(tab)])[1] else NA_character_
  .(
    n_cells_verified = .N,
    response_verified = sort(names(sort(table(response), decreasing = TRUE)))[1],
    high_state_fraction = mean(YS_state == "High"),
    dominant_author_metaprogram = dominant,
    dominant_metaprogram_fraction = if (length(tab)) as.numeric(max(tab)) / .N else NA_real_
  )
}, by = patient]
mp <- merge(mp, patient_annotation, by = "patient", all.x = TRUE, validate = "many-to-one")
stopifnot(all(mp$n_cells == mp$n_cells_verified), uniqueN(mp$patient) == 77L)

rob <- read_tsv(paths[["module_robustness"]])[module %in% modules]
rob_wide <- dcast(
  rob,
  module ~ version,
  value.var = c("n_patients_testable", "n_positive", "median_high_minus_low", "BH_FDR_within_version", "testability_note")
)
mp <- merge(mp, rob_wide, by = "module", all.x = TRUE, validate = "many-to-one")
setorder(mp, module_order, patient)

harm <- read_tsv(paths[["module_harmonized_summary"]])
summary_check <- mp[, .(
  derived_median_effect_z = median(effect_z), derived_positive_n = sum(effect_z > 0), derived_n = .N
), by = module]
harm_check <- harm[version == "score_independent", .(
  module, authoritative_median_effect_z = median_delta_z,
  authoritative_positive_n = n_positive, authoritative_n = n_patients_testable
)]
summary_check <- merge(summary_check, harm_check, by = "module", all.x = TRUE)
summary_check[, `:=`(
  median_match = abs(derived_median_effect_z - authoritative_median_effect_z) < 5e-7,
  count_match = derived_positive_n == authoritative_positive_n & derived_n == authoritative_n
)]
print(summary_check)
stopifnot(
  all(summary_check[is.finite(authoritative_median_effect_z), median_match]),
  all(summary_check[is.finite(authoritative_median_effect_z), count_match])
)
write_tsv_both(mp, "Fig5_module_patient_landscape.tsv")
write_tsv_both(summary_check, "Fig5_module_patient_landscape_QC.tsv")

# -----------------------------------------------------------------------------
# Orthogonal Yan author-metaprogram alignment. M02/M03 are non-testable and are
# excluded from the display source but retained as an explicit availability note.
# -----------------------------------------------------------------------------
align <- read_tsv(paths[["metaprogram"]])
display_mps <- c("M01", sprintf("M%02d", 4:13))
feature_order <- c("YAP", "Stemness", "Joint axis", "UPR", "TNFA NFKB", "Hypoxia",
                   "Adhesion Remodeling", "Wound Healing", "Survival Stress", "Anoikis Resistance")
align_use <- align[metaprogram %in% display_mps & feature %in% feature_order]
align_use[, `:=`(
  metaprogram_order = match(metaprogram, display_mps),
  feature_order = match(feature, feature_order),
  orthogonal_annotation_source = "Yan author-defined dominant cell_state labels",
  M02_M03_status = "not observed/testable among Yan malignant cells; not plotted"
)]
setorder(align_use, record_type, feature_order, metaprogram_order, patient)
stopifnot(nrow(align_use[record_type == "summary"]) == length(display_mps) * length(feature_order))
write_tsv_both(align_use, "Fig5_metaprogram_alignment.tsv")

# -----------------------------------------------------------------------------
# Missing analysis 1: complete 146-gene Wu-Yan effect table.
# The interim H5AD-backed rerun is checked against the exact 0728 gate table.
# -----------------------------------------------------------------------------
strict <- read_tsv(paths[["strict_source"]])[program_146 == TRUE]
stopifnot(nrow(strict) == 146L)
gate6 <- read_tsv(paths[["direction_gate"]])[gate == "6of8"]
interim_gene <- read_tsv(paths[["interim_gene"]])[
  record_type == "gene_summary" & model == "raw",
  .(gene, Yan_patient_support_n = n_positive, Yan_testable_patients = n_testable,
    Yan_effect_rerun = median_Yan_effect, Yan_patient_unit_BH_FDR = BH_FDR,
    Yan_conservation_class_rerun = Yan_conservation_class)
]
complete <- merge(strict, gate6[, .(
  gene, Yan_direction_gate_n_patients = n_patients,
  Yan_effect_authoritative = median_Yan_mean_high_minus_other,
  Yan_direction_positive_authoritative = Yan_direction_positive
)], by = "gene", all.x = TRUE)
complete <- merge(complete, interim_gene, by = "gene", all.x = TRUE)
complete[, `:=`(
  detectable_in_Yan = is.finite(Yan_effect_authoritative),
  Yan_effect = Yan_effect_authoritative,
  Yan_effect_direction = fifelse(is.na(Yan_effect_authoritative), "unavailable",
                          fifelse(Yan_effect_authoritative > 0, "positive",
                          fifelse(Yan_effect_authoritative < 0, "negative", "near_zero"))),
  direction_gate_tier = fifelse(positive_patient_count == 8, "8/8",
                         fifelse(positive_patient_count == 7, "7/8", "6/8")),
  effect_tier = fifelse(overall_log2FC_High_vs_Other >= 1.5, ">=1.5 large-effect reference",
                fifelse(overall_log2FC_High_vs_Other >= 1.0, "1.0-<1.5 stringent-effect subset", "0.5-<1.0 primary-program member")),
  availability_status = fifelse(is.finite(Yan_effect_authoritative), "measurable", "not measurable in Yan H5AD feature space"),
  validation_rule = "Yan median of within-patient mean log-normalized expression High minus Other; no Yan gene reselection"
)]
complete[, rerun_difference := Yan_effect_rerun - Yan_effect_authoritative]
gene_qc <- data.table(
  metric = c("program_genes", "detectable_in_Yan", "positive_direction", "near_zero_or_reversed",
             "max_abs_rerun_vs_gate_difference", "exact_6of8", "exact_7of8", "exact_8of8",
             "effect_0.5_to_1.0", "effect_1.0_to_1.5", "effect_ge_1.5"),
  value = c(
    nrow(complete), sum(complete$detectable_in_Yan), sum(complete$Yan_effect_direction == "positive"),
    sum(complete$Yan_effect_direction %in% c("near_zero", "negative")),
    max(abs(complete$rerun_difference), na.rm = TRUE),
    sum(complete$direction_gate_tier == "6/8"), sum(complete$direction_gate_tier == "7/8"),
    sum(complete$direction_gate_tier == "8/8"),
    sum(complete$effect_tier == "0.5-<1.0 primary-program member"),
    sum(complete$effect_tier == "1.0-<1.5 stringent-effect subset"),
    sum(complete$effect_tier == ">=1.5 large-effect reference")
  )
)
stopifnot(
  nrow(complete) == 146L,
  sum(complete$detectable_in_Yan) == 141L,
  sum(complete$Yan_effect_direction == "positive") == 139L,
  max(abs(complete$rerun_difference), na.rm = TRUE) < 1e-6,
  sum(complete$direction_gate_tier == "6/8") == 16L,
  sum(complete$direction_gate_tier == "7/8") == 68L,
  sum(complete$direction_gate_tier == "8/8") == 62L,
  sum(complete$effect_tier == "0.5-<1.0 primary-program member") == 110L,
  sum(complete$effect_tier == "1.0-<1.5 stringent-effect subset") == 33L,
  sum(complete$effect_tier == ">=1.5 large-effect reference") == 3L
)
write_tsv_both(complete, "Fig5_146gene_Wu_Yan_complete.tsv")
write_tsv_both(gene_qc, "Fig5_146gene_Wu_Yan_QC.tsv")

# -----------------------------------------------------------------------------
# Missing analysis 2: unified three-dataset conservation object.
# -----------------------------------------------------------------------------
gse <- read_tsv(paths[["gse_transfer"]])
three <- merge(complete, gse[, .(
  gene, GSE_detectable = detectable_in_GSE180286,
  GSE_effect = GSE180286_High_vs_Other_effect,
  GSE_effect_direction = GSE180286_effect_direction,
  GSE_expected_direction_preserved = expected_direction_preserved,
  GSE_unavailability_reason = unavailability_reason
)], by = "gene", all.x = TRUE)
three[, `:=`(
  Wu_effect = overall_log2FC_High_vs_Other,
  Wu_direction = "positive_by_frozen_program_definition",
  Yan_available = detectable_in_Yan,
  Yan_positive = detectable_in_Yan & Yan_effect > 0,
  GSE_available = GSE_detectable == TRUE & is.finite(GSE_effect),
  GSE_positive = GSE_detectable == TRUE & is.finite(GSE_effect) & GSE_effect > 0
)]
three[, three_dataset_status := fifelse(!(Yan_available & GSE_available), "not_testable_all_three",
                                 fifelse(Yan_positive & GSE_positive, "same_positive_direction_all_available_datasets", "discordant"))]
three[, GSE_caveat := "GSE180286 is an orthogonal cross-dataset transfer dominated by one P5 Primary/LN+ source; not an equally powered patient cohort"]
three_qc <- data.table(
  metric = c("frozen_program_genes", "Yan_available", "GSE_available", "available_all_three",
             "same_positive_direction_all_three", "discordant_among_all_three_testable", "not_testable_all_three"),
  value = c(nrow(three), sum(three$Yan_available), sum(three$GSE_available),
            sum(three$Yan_available & three$GSE_available),
            sum(three$three_dataset_status == "same_positive_direction_all_available_datasets"),
            sum(three$three_dataset_status == "discordant"),
            sum(three$three_dataset_status == "not_testable_all_three"))
)
stopifnot(nrow(three) == 146L, sum(three$GSE_available) == 141L)
write_tsv_both(three, "Fig5_threecohort_gene_conservation.tsv")
write_tsv_both(three_qc, "Fig5_threecohort_gene_conservation_QC.tsv")

# Program-level transfer: keep exact existing patient rows and summaries.
program_transfer <- read_tsv(paths[["program_patient"]])[
  definition == "Primary pooled tertile" & model %in% c("raw", "technical_adjusted")
]
program_transfer[, `:=`(
  strict_v2_size = 146L,
  detectable_genes = 141L,
  state_contrast = "High versus Other",
  technical_covariates = fifelse(model == "technical_adjusted",
    "log1p(nCount_RNA) + nFeature_RNA + S_score + G2M_score", "none"),
  validation_role = "predefined program-level external transfer; no Yan gene reselection"
)]
write_tsv_both(program_transfer, "Fig5_program_transfer.tsv")

# Compact robustness source: exact pre-existing thresholds 20/50/100 + Yan null.
robust_threshold <- read_tsv(paths[["patient_meta"]])[
  record_type == "meta_summary" & threshold %in% c(20, 50, 100)
]
robust_threshold[, robustness_type := "eligibility_threshold"]
null <- read_tsv(paths[["null"]])[
  cohort == "Yan2026" & feature == "YAP-Stem axis" & record_type == "axis_null_summary" &
    model %in% c("raw", "technical_adjusted")
]
null[, robustness_type := "matched_random_signature_null"]
robust <- rbindlist(list(robust_threshold, null), fill = TRUE)
write_tsv_both(robust, "Fig5_robustness_summary.tsv")

# Honest clinical null is supplementary only.
response <- read_tsv(paths[["response"]])
response[, placement := "SUPPLEMENTARY_NULL"]
fwrite(response, file.path(source_out, "Fig5_pCR_RD_supplementary_null.tsv"), sep = "\t", quote = FALSE, na = "NA")

# Provenance hashes.
hashes <- data.table(
  source_role = names(paths), source_file = unname(paths),
  sha256 = vapply(unname(paths), digest, character(1), algo = "sha256", file = TRUE)
)
fwrite(hashes, file.path(source_out, "Figure5_final_source_hashes.tsv"), sep = "\t", quote = FALSE, na = "NA")

# QC reports for the two missing analyses.
qc1 <- c(
  "# Fig5 146-gene Wu-Yan complete table QC",
  "",
  "- Frozen program: strict_v2, 146 genes; High vs Other; 16 Wu patient-state pseudobulks; 38 score genes excluded before testing.",
  sprintf("- Yan measurable: %d/146.", sum(complete$detectable_in_Yan)),
  sprintf("- Yan positive median High-vs-Other direction: %d/%d measurable (%.1f%%).", sum(complete$Yan_effect_direction == "positive"), sum(complete$detectable_in_Yan), 100 * mean(complete$Yan_effect_direction[complete$detectable_in_Yan] == "positive")),
  sprintf("- H5AD-backed rerun vs authoritative 0728 direction-gate effect: maximum absolute difference %.3g.", max(abs(complete$rerun_difference), na.rm = TRUE)),
  "- Direction tiers reproduced: 16 exact 6/8, 68 exact 7/8, 62 exact 8/8.",
  "- Effect tiers reproduced: 110 at 0.5-<1.0, 33 at 1.0-<1.5, 3 at >=1.5.",
  "- No Yan-based gene reselection was performed.",
  "",
  "Decision: PASS."
)
writeLines(qc1, file.path(working, "Fig5_146gene_Wu_Yan_QC.md"), useBytes = TRUE)

qc2 <- c(
  "# Fig5 unified three-dataset conservation QC",
  "",
  sprintf("- Frozen program genes: %d.", nrow(three)),
  sprintf("- Yan available: %d; GSE180286 available: %d.", sum(three$Yan_available), sum(three$GSE_available)),
  sprintf("- Available in both validation datasets: %d.", sum(three$Yan_available & three$GSE_available)),
  sprintf("- Same positive direction across Wu, Yan and GSE among fully testable genes: %d.", sum(three$three_dataset_status == "same_positive_direction_all_available_datasets")),
  sprintf("- Discordant among fully testable genes: %d; not testable across all three: %d.", sum(three$three_dataset_status == "discordant"), sum(three$three_dataset_status == "not_testable_all_three")),
  "- GSE180286 is explicitly treated as cross-dataset transfer dominated by one P5 Primary/LN+ source, not as an equally powered patient cohort.",
  "- The three-dataset table is descriptive and does not redefine the original 146-gene program.",
  "",
  "Decision: PASS WITH GSE SAMPLE-STRUCTURE CAVEAT."
)
writeLines(qc2, file.path(working, "Fig5_threecohort_conservation_QC.md"), useBytes = TRUE)

cat("Final Figure 5 data closure completed.\n")
cat("Progressive patient records:", nrow(rho), "\n")
cat("Module patient records:", nrow(mp), "\n")
cat("Metaprogram summary cells:", nrow(align_use[record_type == "summary"]), "\n")
cat("Yan detectable strict_v2 genes:", sum(complete$detectable_in_Yan), "/146\n")
cat("Three-dataset fully testable genes:", sum(three$Yan_available & three$GSE_available), "\n")
