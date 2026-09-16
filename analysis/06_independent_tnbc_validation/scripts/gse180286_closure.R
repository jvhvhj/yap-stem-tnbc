#!/usr/bin/env Rscript
# Purpose: GSE180286 closure analysis
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 5; orthogonal transfer in Fig. 4.

options(stringsAsFactors = FALSE, width = 220)
set.seed(20260718)

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
work <- file.path(root, "0718_GSE180286_primary_LN_singlecell_closure")
tab <- file.path(work, "tables")
fig <- file.path(work, "figures")
repdir <- file.path(work, "reports")
logdir <- file.path(work, "logs")
plotdir <- file.path(tab, "Supplementary_GSE180286_plotting_data")
for (d in c(tab, fig, repdir, logdir, plotdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

.libPaths(c(file.path(root, "0703_rebuild/R_library"), .libPaths()))
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

log_con <- file(file.path(logdir, "01_GSE180286_closure.log"), "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message", append = TRUE)

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  cat("Wrote", basename(path), "rows=", nrow(x), "\n")
}
zscore <- function(x) {
  y <- as.numeric(scale(x))
  y[!is.finite(y)] <- 0
  y
}
mean_score <- function(mat, genes) {
  present <- intersect(genes, rownames(mat))
  if (!length(present)) return(rep(NA_real_, ncol(mat)))
  Matrix::colMeans(mat[present, , drop = FALSE])
}
safe_rho <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10L || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

cat("GSE180286 closure started", format(Sys.time()), "\n")

# Official 15-sample design from GEO and Xu et al. 2021. For P1-P3 the paper
# reports one LN+ and one LN- but does not map LN1/LN2 in the GEO sample title.
sample_manifest <- tibble::tribble(
  ~gsm, ~local_sample, ~patient, ~subtype, ~site, ~node_index, ~site_status, ~site_status_basis,
  "GSM5457199", "GSM5457199_A2019-1", "P1", "Luminal B", "Primary", NA_integer_, "Primary", "GEO sample title",
  "GSM5457200", "GSM5457200_A2019-2", "P1", "Luminal B", "LN", 1L, "LN status unresolved", "paper states one LN+ and one LN-; per-node mapping absent from GEO title",
  "GSM5457201", "GSM5457201_A2019-3", "P1", "Luminal B", "LN", 2L, "LN status unresolved", "paper states one LN+ and one LN-; per-node mapping absent from GEO title",
  "GSM5457202", "GSM5457202_B2019-1", "P2", "HER2-enriched", "Primary", NA_integer_, "Primary", "GEO sample title",
  "GSM5457203", "GSM5457203_B2019-2", "P2", "HER2-enriched", "LN", 1L, "LN status unresolved", "paper states one LN+ and one LN-; per-node mapping absent from GEO title",
  "GSM5457204", "GSM5457204_B2019-3", "P2", "HER2-enriched", "LN", 2L, "LN status unresolved", "paper states one LN+ and one LN-; per-node mapping absent from GEO title",
  "GSM5457205", "GSM5457205_C2020-1", "P3", "HER2-enriched", "Primary", NA_integer_, "Primary", "GEO sample title",
  "GSM5457206", "GSM5457206_C2020-2", "P3", "HER2-enriched", "LN", 1L, "LN status unresolved", "paper states one LN+ and one LN-; per-node mapping absent from GEO title",
  "GSM5457207", "GSM5457207_C2020-3", "P3", "HER2-enriched", "LN", 2L, "LN status unresolved", "paper states one LN+ and one LN-; per-node mapping absent from GEO title",
  "GSM5457208", "GSM5457208_D2020-1", "P4", "TNBC", "Primary", NA_integer_, "Primary", "GEO sample title",
  "GSM5457209", "GSM5457209_D2020-2", "P4", "TNBC", "LN", 1L, "LN-", "Xu et al.: both P4 nodes are LN-",
  "GSM5457210", "GSM5457210_D2020-3", "P4", "TNBC", "LN", 2L, "LN-", "Xu et al.: both P4 nodes are LN-",
  "GSM5457211", "GSM5457211_E2020-1", "P5", "TNBC", "Primary", NA_integer_, "Primary", "GEO sample title",
  "GSM5457212", "GSM5457212_E2020-2", "P5", "TNBC", "LN", 1L, "LN+", "Xu et al.: both P5 nodes are LN+",
  "GSM5457213", "GSM5457213_E2020-3", "P5", "TNBC", "LN", 2L, "LN+", "Xu et al.: both P5 nodes are LN+"
)

rds_path <- file.path(root, "output_step0/GSE180286_validated.rds")
if (!file.exists(rds_path)) stop("Missing malignant-cell RDS: ", rds_path)
obj <- readRDS(rds_path)
DefaultAssay(obj) <- "RNA"
meta <- obj@meta.data %>% tibble::rownames_to_column("cell_id")
stopifnot(nrow(meta) == 5469L, all(meta$copykat.pred == "aneuploid"))

meta <- meta %>%
  left_join(sample_manifest, by = c("orig.ident" = "local_sample")) %>%
  mutate(
    malignant_label = "CopyKAT aneuploid",
    in_old_5469_object = TRUE,
    tnbc = subtype == "TNBC",
    analysis_site = case_when(site == "Primary" ~ "Primary", site_status == "LN+" ~ "LN+", site_status == "LN-" ~ "LN-", TRUE ~ "LN unresolved")
  )
if (anyNA(meta$patient)) stop("At least one local sample could not be mapped to the official sample manifest")

counts <- meta %>% count(patient, orig.ident, subtype, analysis_site, site, site_status, malignant_label, name = "n_malignant_cells")
manifest_audit <- sample_manifest %>%
  left_join(counts %>% select(orig.ident, n_malignant_cells) %>% rename(local_sample = orig.ident), by = "local_sample") %>%
  mutate(n_malignant_cells = coalesce(n_malignant_cells, 0L), present_in_old_5469_object = n_malignant_cells > 0)
write_tsv(manifest_audit, file.path(tab, "GSE180286_metadata_audit.tsv"))
write_tsv(counts %>% filter(subtype == "TNBC"), file.path(tab, "GSE180286_TNBC_malignant_cell_counts.tsv"))

old_5469_audit <- counts %>%
  group_by(subtype, analysis_site) %>%
  summarise(n_cells = sum(n_malignant_cells), n_samples = n(), n_patients = n_distinct(patient), .groups = "drop") %>%
  mutate(fraction_of_5469 = n_cells / sum(n_cells))
write_tsv(old_5469_audit, file.path(tab, "GSE180286_old_5469_source_audit.tsv"))

# Frozen signatures: score axes and eight score-independent modules come from
# the locked Yan manifest; boundary modules come from the v3.1 overlap audit.
frozen_manifest <- read.delim(file.path(root, "0716_yan2026_validation/tables/frozen_signature_manifest.tsv"), check.names = FALSE)
yap_genes <- unique(frozen_manifest$gene[frozen_manifest$module == "YAP clean" & frozen_manifest$gene != "NA"])
stem_genes <- unique(frozen_manifest$gene[frozen_manifest$module == "Stemness clean" & frozen_manifest$gene != "NA"])
module_lists <- read.delim(file.path(root, "0716_yan2026_validation/tables/module_score_independent_gene_lists.tsv"), check.names = FALSE) %>%
  group_by(module) %>% summarise(genes = list(unique(gene)), .groups = "drop")
overlap <- read.delim(file.path(root, "0713_score_independent_rebuild/tables/module_gene_overlap.tsv"), check.names = FALSE)
parse_genes <- function(x) {
  z <- trimws(strsplit(ifelse(is.na(x), "", x), ",", fixed = TRUE)[[1]])
  z[nzchar(z)]
}
for (m in c("Cell Cycle", "EMT", "Migration Invasion")) {
  g <- parse_genes(overlap$retained_genes[overlap$module == m][1])
  module_lists <- bind_rows(module_lists, data.frame(module = m, genes = I(list(g))))
}
prolif_genes <- c("MKI67", "UBE2C", "CENPF", "BIRC5", "STMN1", "TK1", "PTTG1", "NUSAP1", "TOP2A", "CDK1")
module_lists <- bind_rows(module_lists, data.frame(module = "Proliferation", genes = I(list(prolif_genes))))
representative_genes <- c("MCL1", "PPP1R15A", "ATF3", "LAMB3", "LAMC2", "TNF", "CXCL2")

signature_manifest <- bind_rows(
  data.frame(signature_type = "score", signature = "YAP", gene = yap_genes),
  data.frame(signature_type = "score", signature = "Stemness", gene = stem_genes),
  bind_rows(lapply(seq_len(nrow(module_lists)), function(i) data.frame(signature_type = "module", signature = module_lists$module[i], gene = module_lists$genes[[i]]))),
  data.frame(signature_type = "representative_gene", signature = representative_genes, gene = representative_genes)
) %>% distinct()

all_rna_genes <- rownames(obj[["RNA"]])
signature_manifest <- signature_manifest %>%
  mutate(present = gene %in% all_rna_genes, source = case_when(
    signature %in% c("YAP", "Stemness") ~ "0716 frozen_signature_manifest.tsv",
    signature == "Proliferation" ~ "0713 frozen boundary definition",
    signature %in% c("Cell Cycle", "EMT", "Migration Invasion") ~ "0713 module_gene_overlap retained genes",
    signature_type == "representative_gene" ~ "0718 task locked representative genes",
    TRUE ~ "0716 module_score_independent_gene_lists.tsv"
  ))
write_tsv(signature_manifest, file.path(tab, "GSE180286_frozen_signature_projection_manifest.tsv"))

# Score the CopyKAT-aneuploid cells layer by layer, avoiding any join of the
# full Seurat object. This keeps peak memory bounded by one sample layer.
score_names <- c("YAP", "Stemness", module_lists$module, paste0("gene__", representative_genes))
for (nm in score_names) meta[[nm]] <- NA_real_
rna <- obj[["RNA"]]
layers <- Layers(rna)
data_layers <- layers[grepl("^data\\.", layers)]
cat("RNA layers used:", paste(data_layers, collapse = ", "), "\n")

for (layer in data_layers) {
  mat <- LayerData(rna, layer = layer)
  cells <- intersect(colnames(mat), meta$cell_id)
  if (!length(cells)) next
  mat <- mat[, cells, drop = FALSE]
  idx <- match(cells, meta$cell_id)
  meta$YAP[idx] <- mean_score(mat, yap_genes)
  meta$Stemness[idx] <- mean_score(mat, stem_genes)
  for (i in seq_len(nrow(module_lists))) meta[[module_lists$module[i]]][idx] <- mean_score(mat, module_lists$genes[[i]])
  for (g in representative_genes) {
    meta[[paste0("gene__", g)]][idx] <- if (g %in% rownames(mat)) as.numeric(mat[g, ]) else NA_real_
  }
  rm(mat); gc(verbose = FALSE)
}

if (anyNA(meta$YAP) || anyNA(meta$Stemness)) stop("Some cells were not covered by normalized RNA data layers")
tnbc <- meta %>% filter(tnbc)
tnbc$YAP_z <- zscore(tnbc$YAP)
tnbc$Stemness_z <- zscore(tnbc$Stemness)
tnbc$Joint_axis <- tnbc$YAP_z + tnbc$Stemness_z
joint_cut <- quantile(tnbc$Joint_axis, c(1 / 3, 2 / 3), na.rm = TRUE)
tnbc$Joint_state <- cut(tnbc$Joint_axis, breaks = c(-Inf, joint_cut, Inf), labels = c("Low", "Intermediate", "High"), include.lowest = TRUE)
write_tsv(tnbc %>% select(cell_id, orig.ident, patient, subtype, analysis_site, malignant_label, nCount_RNA, nFeature_RNA,
                          YAP, Stemness, YAP_z, Stemness_z, Joint_axis, Joint_state, all_of(module_lists$module), starts_with("gene__")),
          file.path(tab, "GSE180286_TNBC_cell_scores.tsv"))

gene_vars <- paste0("gene__", representative_genes)
score_vars <- c("YAP", "Stemness", "Joint_axis", module_lists$module, gene_vars)
sample_scores <- tnbc %>%
  group_by(patient, orig.ident, subtype, analysis_site, site_status) %>%
  summarise(n_cells = n(), high_state_fraction = mean(Joint_state == "High"),
            across(all_of(setdiff(score_vars, gene_vars)), ~ median(.x, na.rm = TRUE)),
            across(all_of(gene_vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
write_tsv(sample_scores, file.path(tab, "GSE180286_sample_level_scores.tsv"))

pseudobulk <- tnbc %>%
  group_by(patient, orig.ident, analysis_site) %>%
  summarise(n_cells = n(), across(all_of(c(module_lists$module, paste0("gene__", representative_genes))), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
write_tsv(pseudobulk, file.path(tab, "GSE180286_patient_sample_pseudobulk.tsv"))

site_rho <- tnbc %>%
  group_by(patient, orig.ident, analysis_site) %>%
  summarise(n_cells = n(), rho_YAP_Stem = safe_rho(YAP, Stemness), .groups = "drop")
write_tsv(site_rho, file.path(tab, "GSE180286_site_yap_stem_rho.tsv"))

# With the available malignant-cell object, only P5 has a TNBC Primary/LN+
# pair. P4 contributes LN- cells only, so Primary/LN- is not estimable.
metric_class <- tibble(metric = c("YAP", "Stemness", "Joint_axis", "high_state_fraction", module_lists$module, paste0("gene__", representative_genes)),
                       family = c(rep("axis", 4), rep("module", nrow(module_lists)), rep("gene", length(representative_genes))))
site_pairs <- list(c("Primary", "LN+"), c("Primary", "LN-"), c("LN-", "LN+"))
effect_rows <- list()
for (pat in unique(tnbc$patient)) {
  d <- sample_scores %>% filter(patient == pat)
  for (pair in site_pairs) {
    a <- d %>% filter(analysis_site == pair[1])
    b <- d %>% filter(analysis_site == pair[2])
    for (metric in metric_class$metric) {
      value_a <- if (nrow(a)) median(a[[metric]], na.rm = TRUE) else NA_real_
      value_b <- if (nrow(b)) median(b[[metric]], na.rm = TRUE) else NA_real_
      effect_rows[[length(effect_rows) + 1L]] <- data.frame(
        patient = pat, comparison = paste(pair, collapse = " vs "), reference_site = pair[1], contrast_site = pair[2],
        metric = metric, family = metric_class$family[match(metric, metric_class$metric)],
        n_reference_samples = nrow(a), n_contrast_samples = nrow(b), reference_value = value_a, contrast_value = value_b,
        contrast_minus_reference = ifelse(is.finite(value_a) && is.finite(value_b), value_b - value_a, NA_real_)
      )
    }
  }
}
effects <- bind_rows(effect_rows) %>%
  mutate(testability = ifelse(is.finite(contrast_minus_reference), "descriptively testable; n=1 patient", "not testable"))
write_tsv(effects, file.path(tab, "GSE180286_patient_specific_site_effects.tsv"))

module_groups <- list(
  "stress-adaptive" = c("UPR", "TNFA NFKB", "Hypoxia", "Survival Stress"),
  "adhesion-remodeling" = c("Adhesion Remodeling", "Wound Healing"),
  "associated/supportive" = c("Anoikis Resistance", "Integrin Adhesion"),
  "boundary" = c("EMT", "Migration Invasion", "Cell Cycle", "Proliferation")
)
group_lookup <- bind_rows(lapply(names(module_groups), function(g) data.frame(metric = module_groups[[g]], branch = g)))

p5_effect <- effects %>% filter(patient == "P5", comparison == "Primary vs LN+")
conclusion <- p5_effect %>%
  mutate(display_metric = sub("^gene__", "", metric),
         predefined_class = case_when(
           display_metric %in% c("UPR", "TNFA NFKB", "Hypoxia", "Adhesion Remodeling", "Wound Healing", "Survival Stress", "MCL1", "PPP1R15A", "ATF3") ~ "core",
           display_metric == "Anoikis Resistance" ~ "associated but non-independent",
           display_metric %in% c("Integrin Adhesion", "EMT", "Migration Invasion", "Cell Cycle", "Proliferation") ~ "supportive/boundary",
           display_metric %in% c("LAMB3", "LAMC2", "TNF", "CXCL2") ~ "context-dependent",
           TRUE ~ "axis"),
         conclusion_class = case_when(
           !is.finite(contrast_minus_reference) ~ "not testable",
           predefined_class == "context-dependent" ~ "context-dependent",
           contrast_minus_reference > 0 ~ "site-enriched trend",
           contrast_minus_reference < 0 ~ "inconsistent",
           TRUE ~ "inconsistent"),
         inference_limit = "single testable TNBC patient; descriptive direction only") %>%
  select(metric = display_metric, family, predefined_class, contrast_minus_reference, conclusion_class, inference_limit)
not_testable <- data.frame(
  metric = c("Primary vs LN-", "LN+ vs LN-", "patient-level site inference"), family = "design",
  predefined_class = "design", contrast_minus_reference = NA_real_, conclusion_class = "not testable",
  inference_limit = c("P4 Primary has no CopyKAT-aneuploid cells in the old object", "no TNBC patient has both LN+ and LN- malignant-cell samples", "only P5 has a Primary/LN+ malignant-cell pair")
)
conclusion <- bind_rows(conclusion, not_testable)
write_tsv(conclusion, file.path(tab, "GSE180286_conclusion_classification.tsv"))

# Plotting data
pA <- manifest_audit %>% mutate(patient = factor(patient, levels = paste0("P", 5:1)),
                                site_display = ifelse(site == "Primary", "Primary", site_status))
pB <- sample_scores %>% select(patient, orig.ident, analysis_site, n_cells, Joint_axis, high_state_fraction)
pC <- p5_effect %>% filter(metric %in% c("UPR", "TNFA NFKB", "Hypoxia", "Adhesion Remodeling", "Wound Healing", "Survival Stress")) %>%
  left_join(group_lookup, by = "metric")
pD <- p5_effect %>% filter(metric %in% paste0("gene__", c("MCL1", "PPP1R15A", "ATF3"))) %>%
  mutate(metric = sub("gene__", "", metric))
pE <- p5_effect %>% filter(metric %in% c("EMT", "Migration Invasion", "Cell Cycle", "Proliferation"))
write_tsv(pA, file.path(plotdir, "panelA_cohort_site_structure.tsv"))
write_tsv(pB, file.path(plotdir, "panelB_joint_axis.tsv"))
write_tsv(pC, file.path(plotdir, "panelC_six_core_modules.tsv"))
write_tsv(pD, file.path(plotdir, "panelD_core_genes.tsv"))
write_tsv(pE, file.path(plotdir, "panelE_boundary_audit.tsv"))

theme_pub <- theme_minimal(base_family = "Arial", base_size = 8.5) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold", size = 9.5, hjust = 0),
        plot.subtitle = element_text(size = 7.3, color = "#50565B"),
        axis.title = element_text(size = 7.5), axis.text = element_text(size = 7),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6.8),
        strip.text = element_text(face = "bold", size = 7.5), plot.margin = margin(4, 5, 4, 5))
site_cols <- c("Primary" = "#3B6A8F", "LN+" = "#B33B47", "LN-" = "#6B8E62", "LN status unresolved" = "#B5B8BB")

gA <- ggplot(pA, aes(x = site_display, y = patient)) +
  geom_segment(aes(x = "Primary", xend = site_display, yend = patient), color = "#D7D9DB", linewidth = 0.45) +
  geom_point(aes(size = pmax(n_malignant_cells, 1), fill = site_display, alpha = present_in_old_5469_object), shape = 21, color = "white", stroke = 0.35) +
  scale_fill_manual(values = site_cols, drop = FALSE) + scale_alpha_manual(values = c(`TRUE` = 0.95, `FALSE` = 0.18), guide = "none") +
  scale_size_continuous(range = c(1.7, 7), trans = "sqrt", breaks = c(1, 10, 100, 1000), name = "Aneuploid cells") +
  labs(title = "A  Official cohort and recoverable malignant-cell subset", subtitle = "5 patients · 15 samples; TNBC patients P4/P5 highlighted by subtype label", x = NULL, y = NULL, fill = "Site") +
  facet_grid(subtype ~ ., scales = "free_y", space = "free_y") + theme_pub +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))

gB <- ggplot(pB, aes(x = analysis_site, y = Joint_axis, group = patient)) +
  geom_line(data = pB %>% filter(patient == "P5"), color = "#6E7479", linewidth = 0.7) +
  geom_point(aes(fill = analysis_site, size = n_cells), shape = 21, color = "white", stroke = 0.5) +
  geom_text(aes(label = paste0(patient, "\nn=", n_cells)), nudge_x = 0.08, hjust = 0, size = 2.25, lineheight = 0.92) +
  scale_fill_manual(values = site_cols) + scale_size_continuous(range = c(3, 7), guide = "none") +
  labs(title = "B  Joint-axis projection is only pairwise-testable in P5", subtitle = "Points are sample medians; no pooled-cell significance test", x = NULL, y = "Median joint axis", fill = NULL) +
  theme_pub + theme(legend.position = "bottom")

core_order <- rev(c("UPR", "TNFA NFKB", "Hypoxia", "Adhesion Remodeling", "Wound Healing", "Survival Stress"))
gC <- ggplot(pC, aes(x = contrast_minus_reference, y = factor(metric, levels = core_order), color = branch)) +
  geom_vline(xintercept = 0, color = "#9A9EA2", linewidth = 0.45) + geom_segment(aes(x = 0, xend = contrast_minus_reference, yend = factor(metric, levels = core_order)), linewidth = 1.1) +
  geom_point(size = 2.6) +
  scale_color_manual(values = c("stress-adaptive" = "#A13C48", "adhesion-remodeling" = "#237C78")) +
  labs(title = "C  Six core modules: P5 LN+ minus Primary", subtitle = "One patient only: direction is descriptive, not cohort inference", x = "Median-score difference", y = NULL, color = NULL) + theme_pub + theme(legend.position = "bottom")

gD <- ggplot(pD, aes(x = contrast_minus_reference, y = factor(metric, levels = rev(c("MCL1", "PPP1R15A", "ATF3"))))) +
  geom_vline(xintercept = 0, color = "#9A9EA2", linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = contrast_minus_reference, yend = factor(metric, levels = rev(c("MCL1", "PPP1R15A", "ATF3")))), color = "#69457A", linewidth = 1.1) +
  geom_point(color = "#69457A", size = 2.7) +
  labs(title = "D  Core replicated genes in the site comparison", subtitle = "P5 LN+ minus Primary; patient × sample pseudobulk", x = "Mean-expression difference", y = NULL) + theme_pub

gE <- ggplot(pE, aes(x = contrast_minus_reference, y = factor(metric, levels = rev(c("EMT", "Migration Invasion", "Cell Cycle", "Proliferation"))))) +
  geom_vline(xintercept = 0, color = "#9A9EA2", linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = contrast_minus_reference, yend = factor(metric, levels = rev(c("EMT", "Migration Invasion", "Cell Cycle", "Proliferation")))), color = "#667078", linewidth = 1.0) +
  geom_point(aes(shape = contrast_minus_reference > 0), color = "#4D555B", fill = "#D6D9DB", size = 2.6) +
  scale_shape_manual(values = c(`TRUE` = 24, `FALSE` = 25), guide = "none") +
  labs(title = "E  EMT / migration / proliferation boundary audit", subtitle = "Direction only; n=1 paired TNBC patient", x = "P5 LN+ minus Primary", y = NULL) + theme_pub

final_plot <- (gA | gB) / (gC | gD | gE) +
  plot_layout(heights = c(1.18, 1), widths = c(1.08, 0.92)) +
  plot_annotation(
    title = "Supplementary Figure. GSE180286 TNBC primary–lymph-node closure",
    subtitle = "Frozen scores projected to CopyKAT-aneuploid cells; site effects are supportive and patient-aware",
    caption = "Only P5 provides a TNBC Primary/LN+ malignant-cell pair; P4 contributes LN− cells without a matched malignant Primary sample.",
    theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = 14),
                  plot.subtitle = element_text(family = "Arial", size = 9, color = "#454B50"),
                  plot.caption = element_text(family = "Arial", size = 7, color = "#555A5E", hjust = 0))
  )

base <- file.path(fig, "Supplementary_GSE180286_primary_LN_validation_v1")
ggsave(paste0(base, ".pdf"), final_plot, width = 13.2, height = 8.2, units = "in", device = cairo_pdf)
ggsave(paste0(base, ".png"), final_plot, width = 13.2, height = 8.2, units = "in", dpi = 320, bg = "white")
ggsave(paste0(base, ".svg"), final_plot, width = 13.2, height = 8.2, units = "in", device = svglite::svglite)

panel_qc <- tibble::tribble(
  ~panel, ~status, ~scientific_claim, ~statistical_unit, ~limitation,
  "A", "KEEP", "The official cohort has five patients and fifteen samples; the old 5,469-cell object is a seven-sample mixed-subtype aneuploid subset.", "sample", "not all official samples contain cells in the old malignant subset",
  "B", "KEEP", "Joint-axis site comparison is only pairwise-testable in P5.", "patient × sample", "one paired TNBC patient",
  "C", "KEEP", "Six core modules can be assessed directionally in the P5 Primary/LN+ pair.", "patient × sample", "descriptive direction only",
  "D", "KEEP", "MCL1, PPP1R15A and ATF3 can be assessed directionally in the P5 pair.", "patient × sample", "descriptive direction only",
  "E", "KEEP", "EMT, migration, cell-cycle and proliferation remain a boundary audit.", "patient × sample", "descriptive direction only"
)
write_tsv(panel_qc, file.path(tab, "Supplementary_GSE180286_panel_QC.tsv"))

core_delta <- setNames(conclusion$contrast_minus_reference, conclusion$metric)
fmt_delta <- function(x) ifelse(is.finite(x), sprintf("%.3f", x), "not testable")
module_positive <- sum(pC$contrast_minus_reference > 0, na.rm = TRUE)
gene_positive <- sum(pD$contrast_minus_reference > 0, na.rm = TRUE)
boundary_positive <- setNames(pE$contrast_minus_reference > 0, pE$metric)

report <- c(
  "# GSE180286 primary–lymph-node single-cell closure",
  "",
  "## Identity and analysis gate",
  "",
  "- Official design: five patients, fifteen samples (five Primary, ten lymph nodes). P4 and P5 are TNBC; P4 has two LN− nodes and P5 has two LN+ nodes.",
  sprintf("- The old validated object contains %s CopyKAT-aneuploid cells from %s of the 15 samples, not a TNBC-only validation set.", format(nrow(meta), big.mark = ","), n_distinct(meta$orig.ident)),
  sprintf("- TNBC malignant-cell subset: %s cells from %s samples and %s patients.", format(nrow(tnbc), big.mark = ","), n_distinct(tnbc$orig.ident), n_distinct(tnbc$patient)),
  "- Formal site inference is blocked by sample structure: only P5 has both Primary and LN+ malignant cells; P4 has LN− malignant cells but no malignant Primary cells in this object. Therefore, no cohort-level paired P value is reported.",
  "",
  "## Frozen projection results",
  "",
  sprintf("- P5 LN+ minus Primary joint-axis median difference: %s.", fmt_delta(core_delta[["Joint_axis"]])),
  sprintf("- Six-core modules positive in the P5 direction: %s/6.", module_positive),
  sprintf("- Core genes MCL1/PPP1R15A/ATF3 positive in the P5 direction: %s/3.", gene_positive),
  sprintf("- Boundary directions (P5 LN+ minus Primary): EMT=%s; Migration/Invasion=%s; Cell Cycle=%s; Proliferation=%s.",
          fmt_delta(pE$contrast_minus_reference[pE$metric == "EMT"]),
          fmt_delta(pE$contrast_minus_reference[pE$metric == "Migration Invasion"]),
          fmt_delta(pE$contrast_minus_reference[pE$metric == "Cell Cycle"]),
          fmt_delta(pE$contrast_minus_reference[pE$metric == "Proliferation"])),
  "",
  "## Answers to the seven requested questions",
  "",
  "1. **Are TNBC lymph-node samples present?** Yes. P4 contributes LN− malignant cells and P5 contributes LN+ malignant cells in the old CopyKAT-aneuploid object.",
  sprintf("2. **Is the YAP–Stem axis enhanced in LN?** %s in the only testable P5 Primary/LN+ pair (delta %s); not generalizable beyond one patient.", ifelse(core_delta[["Joint_axis"]] > 0, "Directionally supportive", "Inconsistent"), fmt_delta(core_delta[["Joint_axis"]])),
  sprintf("3. **Do the six core modules differ by site?** %s of six point toward LN enrichment in P5; this is a site-enriched trend, not a cohort result.", module_positive),
  sprintf("4. **Is classic EMT/migration enhancement supported?** EMT is %s and Migration/Invasion is %s in P5; retain as a boundary/context result.", ifelse(boundary_positive[["EMT"]], "higher", "not higher"), ifelse(boundary_positive[["Migration Invasion"]], "higher", "not higher")),
  sprintf("5. **Does the data support 'high proliferation, low metastasis'?** Proliferation is %s while Migration/Invasion is %s in the only pair; this small design cannot establish a general high-proliferation/low-migration state.", ifelse(boundary_positive[["Proliferation"]], "higher", "not higher"), ifelse(boundary_positive[["Migration Invasion"]], "higher", "not higher")),
  "6. **Should this enter a main figure?** No. One paired TNBC patient is insufficient for main-line inference.",
  "7. **Where should it go?** Keep as a dedicated supplementary site-validation figure and cite it as a sample-structure-aware closure of the metastasis question.",
  "",
  "## Sources",
  "",
  "- GEO GSE180286 sample titles and 15-sample design.",
  "- Xu et al., Oncogenesis 2021 (doi:10.1038/s41389-021-00355-6), patient subtype and LN+/LN− structure.",
  "- Local malignant-cell object: output_step0/GSE180286_validated.rds (CopyKAT aneuploid subset)."
)
writeLines(report, file.path(repdir, "GSE180286_data_identity_report.md"), useBytes = TRUE)

cat("GSE180286 closure completed", format(Sys.time()), "\n")
sink(type = "message")
sink(type = "output")
close(log_con)
quit(save = "no", status = 0)
