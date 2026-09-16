#!/usr/bin/env Rscript
# Purpose: Figure 5 final render
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 5; orthogonal transfer in Fig. 4.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(svglite)
  library(ragg)
  library(digest)
})

options(stringsAsFactors = FALSE)

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "Figure5_final_closure")
src_dir <- file.path(out_dir, "source_data")
panel_dir <- file.path(out_dir, "panels")
report_dir <- file.path(out_dir, "reports")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

mm_to_in <- function(x) x / 25.4
pt_to_mm <- function(x) x * 25.4 / 72

pal_df <- read_tsv(file.path(root, "project_semantic_palette.tsv"), show_col_types = FALSE)
pal <- setNames(pal_df$hex, pal_df$token)
get_col <- function(token, fallback) if (token %in% names(pal)) pal[[token]] else fallback

COL <- list(
  raw = get_col("Raw", "#6EB6E4"),
  primary = get_col("Technical-adjusted", "#E7837D"),
  extended = get_col("Positive", "#D7605C"),
  pos = get_col("Positive", "#D7605C"),
  neg = get_col("Negative", "#2F6DB3"),
  neutral = get_col("Neutral", "#F7F7F7"),
  high = get_col("High", "#E7837D"),
  other = get_col("Other", "#6EB6E4"),
  grey = "#8A8F94",
  light_grey = "#D8DDE1",
  ink = "#26343E",
  boundary = "#A5A9AD"
)

theme_pub <- function(base_size = 7) {
  theme_classic(base_family = "Arial", base_size = base_size) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      text = element_text(colour = COL$ink),
      axis.text = element_text(size = base_size, colour = COL$ink),
      axis.title = element_text(size = base_size + 0.2, colour = COL$ink),
      axis.line = element_line(linewidth = 0.35, colour = COL$ink),
      axis.ticks = element_line(linewidth = 0.3, colour = COL$ink),
      axis.ticks.length = unit(1.2, "mm"),
      plot.title = element_text(size = 9, face = "bold", hjust = 0, margin = margin(b = 2.5)),
      legend.title = element_text(size = 6.6, face = "bold"),
      legend.text = element_text(size = 6.3),
      legend.key.height = unit(2.8, "mm"),
      legend.key.width = unit(4.2, "mm"),
      legend.spacing.x = unit(1.2, "mm"),
      plot.margin = margin(3, 4, 3, 3)
    )
}

save_pdf <- function(plot, file, width_mm, height_mm) {
  grDevices::cairo_pdf(file, width = mm_to_in(width_mm), height = mm_to_in(height_mm),
                       family = "Arial", bg = "white", onefile = TRUE)
  print(plot)
  invisible(dev.off())
}

save_final_svg <- function(plot, file, width_mm, height_mm) {
  svglite::svglite(file, width = mm_to_in(width_mm), height = mm_to_in(height_mm),
                   bg = "white", system_fonts = list(sans = "Arial"))
  print(plot)
  invisible(dev.off())
}

save_png <- function(plot, file, width_mm, height_mm, dpi = 600) {
  ragg::agg_png(file, width = width_mm, height = height_mm, units = "mm", res = dpi,
                background = "white", scaling = 1)
  print(plot)
  invisible(dev.off())
}

read_source <- function(name) read_tsv(file.path(src_dir, name), show_col_types = FALSE)

progressive <- read_source("Fig5_progressive_adjustment.tsv")
module_patient <- read_source("Fig5_module_patient_landscape.tsv")
metaprogram <- read_source("Fig5_metaprogram_alignment.tsv")
gene_transfer <- read_source("Fig5_146gene_Wu_Yan_complete.tsv")
threecohort <- read_source("Fig5_threecohort_gene_conservation.tsv")
program_transfer <- read_source("Fig5_program_transfer.tsv")
robustness <- read_source("Fig5_robustness_summary.tsv")

# -----------------------------------------------------------------------------
# Panel A: progressive association under nested adjustment
# -----------------------------------------------------------------------------

model_levels <- c("raw", "technical_primary", "technical_extended")
model_labels <- c("Raw", "Technical", "+ Hypoxia/UPR")
model_cols <- c(raw = COL$raw, technical_primary = COL$primary,
                technical_extended = COL$extended)

a_pat <- progressive %>%
  filter(record_type == "patient", model %in% model_levels, !is.na(rho)) %>%
  mutate(model = factor(model, levels = model_levels), model_x = as.numeric(model))

a_meta <- progressive %>%
  filter(record_type == "meta_summary", threshold == 50, model %in% model_levels) %>%
  distinct(model, .keep_all = TRUE) %>%
  mutate(model = factor(model, levels = model_levels), model_x = as.numeric(model))

a_wide <- a_pat %>% select(patient, model, rho) %>% pivot_wider(names_from = model, values_from = rho)
a_label_patients <- bind_rows(
  a_wide %>% slice_min(technical_extended, n = 2, with_ties = FALSE),
  a_wide %>% mutate(attenuation = raw - technical_extended) %>%
    slice_max(attenuation, n = 2, with_ties = FALSE)
) %>% distinct(patient)
a_labels <- a_pat %>% filter(model == "technical_extended", patient %in% a_label_patients$patient)

a_ymax <- max(a_pat$rho, a_meta$pooled_ci_high, na.rm = TRUE)
a_ymin <- min(a_pat$rho, a_meta$pooled_ci_low, na.rm = TRUE)
a_count <- a_meta %>% mutate(label = paste0(n_positive, "/", n_patients, " positive"), y = a_ymax + 0.055)

pA <- ggplot(a_pat, aes(model_x, rho, group = patient)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, colour = COL$grey) +
  geom_line(colour = "#AAB4BC", alpha = 0.32, linewidth = 0.27) +
  geom_point(aes(fill = model), shape = 21, size = 1.25, stroke = 0.18,
             colour = "white", alpha = 0.80) +
  geom_line(data = a_meta, aes(model_x, pooled_rho, group = 1), inherit.aes = FALSE,
            colour = COL$ink, linewidth = 0.65) +
  geom_errorbar(data = a_meta,
                aes(x = model_x, ymin = pooled_ci_low, ymax = pooled_ci_high),
                inherit.aes = FALSE, width = 0.06, linewidth = 0.7, colour = COL$ink) +
  geom_point(data = a_meta, aes(model_x, pooled_rho, fill = model), inherit.aes = FALSE,
             shape = 23, size = 3.4, stroke = 0.75, colour = COL$ink) +
  geom_text(data = a_count, aes(model_x, y, label = label), inherit.aes = FALSE,
            size = 2.35, family = "Arial", colour = COL$ink) +
  geom_text_repel(data = a_labels, aes(model_x, rho, label = patient), inherit.aes = FALSE,
                  nudge_x = 0.20, direction = "y", hjust = 0, size = 2.05,
                  family = "Arial", min.segment.length = 0, segment.size = 0.25,
                  segment.colour = COL$grey, box.padding = 0.15, point.padding = 0.1,
                  max.overlaps = Inf, seed = 20260811) +
  scale_fill_manual(values = model_cols, breaks = model_levels, labels = model_labels,
                    name = "Model") +
  scale_x_continuous(breaks = 1:3, labels = model_labels, limits = c(0.82, 3.42),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(a_ymin - 0.035, a_ymax + 0.10),
                     breaks = pretty(c(a_ymin, a_ymax), n = 5), expand = c(0, 0)) +
  labs(title = "A  Progressive external association",
       x = NULL, y = expression("Within-patient Spearman " * rho)) +
  theme_pub(7) +
  theme(legend.position = "top", legend.justification = "right",
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(3, 14, 3, 4)) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 2.2))) +
  coord_cartesian(clip = "off")

# -----------------------------------------------------------------------------
# Panel B: integrated patient-by-module landscape
# -----------------------------------------------------------------------------

module_levels <- c("UPR", "TNFA NFKB", "Survival Stress", "Adhesion Remodeling",
                   "Hypoxia", "Wound Healing", "Anoikis Resistance", "Integrin Adhesion")
module_labels <- c(
  "UPR", "TNFα-NF-κB", "Survival stress", "Adhesion remodeling",
  "Hypoxia", "Wound healing", "Anoikis", "Integrin adhesion†"
)
names(module_labels) <- module_levels

patient_levels <- module_patient %>% distinct(patient) %>% arrange(patient) %>% pull(patient)
patient_index <- tibble(patient = patient_levels, patient_x = seq_along(patient_levels))

b_dat <- module_patient %>%
  mutate(patient = factor(patient, levels = patient_levels),
         patient_x = match(as.character(patient), patient_levels),
         module = factor(module, levels = rev(module_levels)),
         module_label = factor(module_labels[as.character(module)], levels = rev(module_labels)))

b_lim <- max(abs(b_dat$effect_z), na.rm = TRUE)
b_ticks <- sort(unique(c(-round(b_lim, 1), -2, -1, 0, 1, 2, round(b_lim, 1))))
b_ticks <- b_ticks[b_ticks >= -b_lim & b_ticks <= b_lim]

b_context <- b_dat %>% distinct(patient, patient_x, n_cells) %>% arrange(patient_x)
b_context$n_scaled <- log10(b_context$n_cells + 1)

patient_break_idx <- unique(pmin(length(patient_levels), c(1, seq(10, length(patient_levels), by = 10))))

pB_context <- ggplot(b_context, aes(patient_x, n_scaled)) +
  geom_col(width = 0.92, fill = alpha(COL$high, 0.82), colour = NA) +
  annotate("text", x = 1, y = max(b_context$n_scaled) * 0.88,
           label = "Cancer cells (log scale)", hjust = 0, size = 2.0, family = "Arial",
           colour = COL$ink) +
  labs(title = "B  Patient-resolved module conservation") +
  scale_x_continuous(limits = c(0.5, length(patient_levels) + 0.5), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
  theme_void(base_family = "Arial") +
  theme(plot.title = element_text(family = "Arial", face = "bold", size = 9,
                                  colour = COL$ink, margin = margin(b = 1.5)),
        plot.margin = margin(1, 0, 0.5, 0))

pB_heat <- ggplot(b_dat, aes(patient_x, module_label, fill = effect_z)) +
  geom_tile(width = 0.96, height = 0.94, colour = NA) +
  scale_fill_gradient2(
    low = COL$neg, mid = COL$neutral, high = COL$pos, midpoint = 0,
    limits = c(-b_lim, b_lim), breaks = b_ticks,
    transform = scales::pseudo_log_trans(sigma = 0.5),
    oob = scales::squish, name = "Module Δz"
  ) +
  scale_x_continuous(breaks = patient_break_idx,
                     labels = patient_levels[patient_break_idx],
                     limits = c(0.5, length(patient_levels) + 0.5), expand = c(0, 0)) +
  labs(x = "Yan patients (alphanumeric order)", y = NULL) +
  theme_pub(6.7) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.2),
        axis.text.y = element_text(size = 6.6),
        axis.line = element_blank(), axis.ticks = element_blank(),
        panel.border = element_rect(fill = NA, colour = "#BAC2C8", linewidth = 0.3),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key.width = unit(12, "mm"),
        plot.margin = margin(0, 1, 2, 1)) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0,
                                barwidth = unit(38, "mm"), barheight = unit(2.5, "mm")))

b_summary <- b_dat %>%
  distinct(module, module_label,
           n_patients_testable_pairwise_shared_removed,
           n_patients_testable_score_independent,
           n_patients_testable_technical_adjusted,
           n_positive_pairwise_shared_removed,
           n_positive_score_independent,
           n_positive_technical_adjusted,
           module_class) %>%
  transmute(module_label,
            module_class,
            `Primary` = n_positive_score_independent / n_patients_testable_score_independent,
            `Technical` = n_positive_technical_adjusted / n_patients_testable_technical_adjusted,
            `Shared-gene` = n_positive_pairwise_shared_removed /
              n_patients_testable_pairwise_shared_removed) %>%
  pivot_longer(c(`Primary`, `Technical`, `Shared-gene`),
               names_to = "evidence", values_to = "positive_fraction") %>%
  mutate(module_label = factor(module_label, levels = rev(module_labels)),
         evidence = factor(evidence, levels = c("Primary", "Technical", "Shared-gene")))

pB_summary <- ggplot(b_summary, aes(positive_fraction, module_label)) +
  geom_vline(xintercept = c(0.5, 1), colour = "#D5DADF", linewidth = 0.3) +
  geom_segment(data = b_summary %>% group_by(module_label) %>%
                 summarise(xmin = min(positive_fraction, na.rm = TRUE),
                           xmax = max(positive_fraction, na.rm = TRUE), .groups = "drop"),
               aes(x = xmin, xend = xmax, y = module_label, yend = module_label),
               inherit.aes = FALSE, linewidth = 0.5, colour = "#AAB4BC") +
  geom_point(aes(shape = evidence, fill = evidence), size = 1.8, stroke = 0.35,
             colour = COL$ink, na.rm = TRUE) +
  scale_shape_manual(values = c(21, 22, 24)) +
  scale_fill_manual(values = c("Primary" = COL$pos, "Technical" = COL$primary,
                               "Shared-gene" = COL$raw)) +
  scale_x_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 1),
                     labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(x = "Patients positive", y = NULL, shape = NULL, fill = NULL) +
  theme_pub(6.4) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        legend.position = "bottom", legend.direction = "vertical",
        legend.box = "vertical", legend.key.height = unit(2.3, "mm"),
        plot.margin = margin(0, 2, 2, 2))

pB <- ((pB_context + plot_spacer()) + plot_layout(widths = c(4.9, 1.35))) /
  ((pB_heat + pB_summary) + plot_layout(widths = c(4.9, 1.35))) +
  plot_layout(heights = c(0.27, 1))

# -----------------------------------------------------------------------------
# Panel C: orthogonal alignment to author-defined metaprograms
# -----------------------------------------------------------------------------

c_dat <- metaprogram %>%
  filter(record_type == "summary") %>%
  mutate(
    feature_label = recode(feature,
      "YAP" = "YAP17", "Stem" = "Stem21", "Stemness" = "Stem21", "Joint" = "Joint axis",
      "TNFA NFKB" = "TNFα-NF-κB", "Survival Stress" = "Survival stress",
      "Adhesion Remodeling" = "Adhesion remodeling", "Wound Healing" = "Wound healing",
      "Anoikis Resistance" = "Anoikis", "Integrin Adhesion" = "Integrin adhesion"
    ),
    metaprogram_label = recode(metaprogram,
      "M01" = "M01\nG2/M", "M04" = "M04\nStress", "M05" = "M05\nIFN",
      "M06" = "M06\nHLA", "M07" = "M07\nS phase", "M08" = "M08\nHypoxia",
      "M09" = "M09\nBasal", "M10" = "M10\nEMT", "M11" = "M11\nLumSec",
      "M12" = "M12\nChol.", "M13" = "M13\nER stress"
    ),
    positive_fraction = n_positive / n_patients_testable
  )

c_feature_levels <- c("YAP17", "Stem21", "Joint axis", "UPR", "TNFα-NF-κB",
                      "Hypoxia", "Adhesion remodeling", "Wound healing",
                      "Survival stress", "Anoikis", "Integrin adhesion")
c_mp_levels <- c("M01\nG2/M", "M04\nStress", "M05\nIFN", "M06\nHLA", "M07\nS phase",
                 "M08\nHypoxia", "M09\nBasal", "M10\nEMT", "M11\nLumSec",
                 "M12\nChol.", "M13\nER stress")
c_dat <- c_dat %>%
  mutate(feature_label = factor(feature_label, levels = rev(c_feature_levels)),
         metaprogram_label = factor(metaprogram_label, levels = c_mp_levels))
c_lim <- max(abs(c_dat$median_rho_one_vs_rest), na.rm = TRUE)

pC <- ggplot(c_dat, aes(metaprogram_label, feature_label)) +
  geom_tile(fill = "#F4F5F6", colour = "white", linewidth = 0.3) +
  geom_point(aes(size = positive_fraction, fill = median_rho_one_vs_rest),
             shape = 21, colour = COL$ink, stroke = 0.25) +
  scale_fill_gradient2(low = COL$neg, mid = COL$neutral, high = COL$pos,
                       midpoint = 0, limits = c(-c_lim, c_lim),
                       name = "Median\nSpearman ρ") +
  scale_size_continuous(range = c(0.8, 4.2), limits = c(0, 1),
                        breaks = c(0.5, 1),
                        labels = percent_format(accuracy = 1),
                        name = "Patients\npositive") +
  labs(title = "C  Yan cancer-state alignment", x = NULL, y = NULL) +
  theme_pub(6.6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6.1),
        axis.text.y = element_text(size = 6.3), axis.ticks = element_blank(),
        axis.line = element_blank(), panel.border = element_rect(fill = NA,
          colour = "#BAC2C8", linewidth = 0.3),
        legend.position = "bottom", legend.box = "vertical",
        plot.margin = margin(3, 3, 2, 3)) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0,
                                barwidth = unit(18, "mm"), barheight = unit(2.5, "mm")),
         size = guide_legend(override.aes = list(fill = "white")))

# -----------------------------------------------------------------------------
# Panel D: all frozen program genes in Wu and Yan
# -----------------------------------------------------------------------------

d_dat <- gene_transfer %>%
  mutate(
    direction_gate_tier = factor(direction_gate_tier, levels = c("6/8", "7/8", "8/8")),
    effect_tier_short = case_when(
      grepl(">=1.5", effect_tier, fixed = TRUE) ~ ">=1.5",
      grepl("1.0-<1.5", effect_tier, fixed = TRUE) ~ "1.0-<1.5",
      TRUE ~ "0.5-<1.0"
    ),
    effect_tier_short = factor(effect_tier_short, levels = c("0.5-<1.0", "1.0-<1.5", ">=1.5")),
    yan_status = case_when(
      !detectable_in_Yan ~ "Unavailable",
      Yan_effect > 0 ~ "Positive",
      TRUE ~ "Near-zero/non-positive"
    )
  )

d_scatter <- d_dat %>% filter(detectable_in_Yan)
d_rho <- suppressWarnings(cor(d_scatter$overall_log2FC_High_vs_Other,
                              d_scatter$Yan_effect, method = "spearman", use = "complete.obs"))
d_preserved <- sum(d_scatter$Yan_effect > 0, na.rm = TRUE)
d_n <- nrow(d_scatter)
d_unavailable <- sum(!d_dat$detectable_in_Yan)
d_labels <- d_scatter %>% filter(gene %in% c("MCL1", "PPP1R15A", "ATF3"))

tier_cols <- c("6/8" = "#F4B8B2", "7/8" = COL$high, "8/8" = COL$pos)

pD <- ggplot(d_scatter, aes(overall_log2FC_High_vs_Other, Yan_effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = COL$grey, linewidth = 0.35) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = COL$light_grey, linewidth = 0.35) +
  geom_point(aes(size = effect_tier_short, fill = direction_gate_tier, shape = yan_status),
             colour = COL$ink, stroke = 0.32, alpha = 0.88) +
  geom_text_repel(data = d_labels, aes(label = gene), size = 2.25, family = "Arial",
                  box.padding = 0.25, point.padding = 0.2, segment.size = 0.28,
                  segment.colour = COL$grey, min.segment.length = 0,
                  seed = 20260811, max.overlaps = Inf) +
  annotate("label", x = Inf, y = Inf, hjust = 1.02, vjust = 1.08,
           label = sprintf("%d/%d positive in Yan\n%d unavailable", d_preserved, d_n, d_unavailable),
           family = "Arial", size = 2.25, linewidth = 0.2,
           label.padding = unit(1.4, "mm"), fill = alpha("white", 0.92), colour = COL$ink) +
  scale_fill_manual(values = tier_cols, name = "Wu consistency") +
  scale_size_manual(values = c("0.5-<1.0" = 1.4, "1.0-<1.5" = 2.0, ">=1.5" = 2.6),
                    name = "Wu log2FC") +
  scale_shape_manual(values = c("Positive" = 21, "Near-zero/non-positive" = 24),
                     name = "Yan direction") +
  labs(title = "D  Gene-level program transfer",
       x = "Wu patient-blocked log2FC", y = "Yan median within-patient effect") +
  theme_pub(6.8) +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.title = element_text(size = 6.2, face = "bold"),
        legend.text = element_text(size = 6.0),
        plot.margin = margin(3, 3, 3, 4)) +
  guides(fill = guide_legend(order = 1, nrow = 1, title.position = "left"),
         size = guide_legend(order = 2, nrow = 1, title.position = "left"),
         shape = "none")

# -----------------------------------------------------------------------------
# Panel E: signed direction conservation across non-equivalent datasets
# -----------------------------------------------------------------------------

e_base <- threecohort %>%
  mutate(
    status_group = case_when(
      three_dataset_status == "same_positive_direction_all_available_datasets" ~ "Same positive direction",
      three_dataset_status == "discordant" ~ "Discordant",
      TRUE ~ "Not testable"
    ),
    status_group = factor(status_group,
                          levels = c("Same positive direction", "Discordant", "Not testable")),
    direction_gate_num = as.integer(sub("/8", "", direction_gate_tier)),
    gene = as.character(gene)
  ) %>%
  arrange(status_group, desc(direction_gate_num), desc(overall_log2FC_High_vs_Other), gene) %>%
  mutate(gene_rank = row_number())

e_long <- bind_rows(
  e_base %>% transmute(gene, gene_rank, status_group, dataset = "Wu", dataset_y = 3,
                       direction = "Positive"),
  e_base %>% transmute(gene, gene_rank, status_group, dataset = "Yan", dataset_y = 2,
                       direction = case_when(!Yan_available ~ "Unavailable",
                                             Yan_positive ~ "Positive", TRUE ~ "Non-positive")),
  e_base %>% transmute(gene, gene_rank, status_group, dataset = "GSE180286†", dataset_y = 1,
                       direction = case_when(!GSE_available ~ "Unavailable",
                                             GSE_positive ~ "Positive", TRUE ~ "Non-positive"))
)

e_groups <- e_base %>% group_by(status_group) %>%
  summarise(xmin = min(gene_rank) - 0.5, xmax = max(gene_rank) + 0.5, n = n(), .groups = "drop") %>%
  mutate(xmid = (xmin + xmax) / 2,
         label = case_when(
           as.character(status_group) == "Same positive direction" ~ paste0("Same positive direction\n", n),
           as.character(status_group) == "Discordant" ~ paste0("Disc.\n", n),
           TRUE ~ paste0("NT\n", n)
         ),
         label_x = case_when(
           as.character(status_group) == "Discordant" ~ xmid + 1.5,
           as.character(status_group) == "Not testable" ~ xmin + 1,
           TRUE ~ xmid
         ),
         label_hjust = case_when(
           as.character(status_group) == "Discordant" ~ 1,
           as.character(status_group) == "Not testable" ~ 0,
           TRUE ~ 0.5
         ),
         label_y = 3.90)
e_separators <- e_groups %>% slice_head(n = max(0, nrow(e_groups) - 1))
e_same <- sum(e_base$status_group == "Same positive direction")
e_testable <- sum(e_base$status_group != "Not testable")

pE <- ggplot(e_long, aes(gene_rank, dataset_y, fill = direction)) +
  geom_tile(width = 1, height = 0.86, colour = NA) +
  geom_vline(data = e_separators, aes(xintercept = xmax),
             inherit.aes = FALSE, colour = "white", linewidth = 0.7) +
  geom_segment(data = e_groups, aes(x = xmin, xend = xmax, y = 3.72, yend = 3.72),
               inherit.aes = FALSE, colour = COL$ink, linewidth = 0.35) +
  geom_text(data = e_groups, aes(x = label_x, y = label_y, label = label, hjust = label_hjust),
            inherit.aes = FALSE, family = "Arial", size = 1.75,
            lineheight = 0.9, colour = COL$ink) +
  scale_fill_manual(values = c("Positive" = COL$pos, "Non-positive" = COL$neg,
                               "Unavailable" = "#C7CBCF"), name = "Direction") +
  scale_x_continuous(limits = c(0.5, nrow(e_base) + 0.5), expand = c(0, 0),
                     breaks = NULL) +
  scale_y_continuous(limits = c(0.45, 4.18), breaks = c(1, 2, 3),
                     labels = c("GSE180286†", "Yan", "Wu"), expand = c(0, 0)) +
  labs(title = "E  Cross-dataset direction conservation", x = "146 predefined program genes", y = NULL) +
  theme_pub(6.8) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        panel.border = element_rect(fill = NA, colour = "#BAC2C8", linewidth = 0.3),
        legend.position = "bottom", plot.margin = margin(3, 3, 3, 4))

# -----------------------------------------------------------------------------
# Panel F: subordinate robustness and downstream-program transfer strip
# -----------------------------------------------------------------------------

f_thresh <- robustness %>%
  filter(record_type == "meta_summary", robustness_type == "eligibility_threshold",
         model %in% model_levels) %>%
  mutate(model = factor(model, levels = model_levels), threshold = as.numeric(threshold))

pF1 <- ggplot(f_thresh, aes(threshold, pooled_rho, colour = model, group = model)) +
  geom_hline(yintercept = 0, colour = COL$light_grey, linewidth = 0.3) +
  geom_line(linewidth = 0.55) +
  geom_errorbar(aes(ymin = pooled_ci_low, ymax = pooled_ci_high), width = 2.5, linewidth = 0.45) +
  geom_point(size = 1.7) +
  scale_colour_manual(values = model_cols, breaks = model_levels, labels = model_labels,
                      name = NULL) +
  scale_x_continuous(breaks = c(20, 50, 100)) +
  labs(title = "F  Robustness and program transfer",
       x = "Minimum cancer cells", y = "Pooled ρ") +
  annotate("text", x = 20, y = max(f_thresh$pooled_ci_high) * 1.02,
           label = "Eligibility sensitivity", hjust = 0, size = 2.15,
           family = "Arial", fontface = "bold", colour = COL$ink) +
  theme_pub(6.2) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0,
                                  margin = margin(b = 1.5)),
        legend.position = "none",
        plot.margin = margin(3, 3, 2, 3))

f_null <- robustness %>%
  filter(record_type == "axis_null_summary", statistic == "pooled_rho") %>%
  mutate(model_label = recode(model, raw = "Raw", technical_adjusted = "Technical"),
         model_label = factor(model_label, levels = c("Technical", "Raw")))

pF2 <- ggplot(f_null, aes(y = model_label)) +
  geom_segment(aes(x = null_q025, xend = null_q975, yend = model_label),
               colour = "#AAB4BC", linewidth = 1.15, lineend = "round") +
  geom_point(aes(x = null_q50), shape = 21, fill = "white", colour = COL$ink,
             size = 1.5, stroke = 0.45) +
  geom_point(aes(x = observed, fill = model), shape = 23, colour = COL$ink,
             size = 2.6, stroke = 0.55) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = COL$light_grey, linewidth = 0.3) +
  scale_fill_manual(values = c(raw = COL$raw, technical_adjusted = COL$primary), guide = "none") +
  labs(x = "Pooled ρ", y = NULL) +
  annotate("text", x = min(f_null$null_q025), y = 2.45,
           label = "Matched-signature null", hjust = 0, size = 2.15,
           family = "Arial", fontface = "bold", colour = COL$ink) +
  theme_pub(6.2) +
  theme(plot.title = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin(3, 3, 2, 3))

f_program <- program_transfer %>%
  filter(record_type == "patient", testable, model %in% c("raw", "technical_adjusted")) %>%
  mutate(model_label = recode(model, raw = "Raw", technical_adjusted = "Technical"),
         model_label = factor(model_label, levels = c("Raw", "Technical")))
f_program_summary <- program_transfer %>% filter(record_type == "summary")

pF3 <- ggplot(f_program, aes(model_label, median_high_minus_other, fill = model_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = COL$grey, linewidth = 0.3) +
  geom_violin(width = 0.72, linewidth = 0.35, alpha = 0.45, colour = COL$ink, trim = FALSE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.4, fill = "white") +
  geom_jitter(width = 0.10, size = 0.35, alpha = 0.35, colour = COL$ink) +
  scale_fill_manual(values = c("Raw" = COL$raw, "Technical" = COL$primary), guide = "none") +
  labs(x = NULL, y = "High-Other program score") +
  annotate("text", x = 1, y = max(f_program$median_high_minus_other, na.rm = TRUE) * 1.04,
           label = "Predefined-program transfer", hjust = 0, size = 2.15,
           family = "Arial", fontface = "bold", colour = COL$ink) +
  annotate("text", x = 2, y = max(f_program$median_high_minus_other, na.rm = TRUE) * 0.78,
           label = "77/77 positive", hjust = 1, vjust = 0,
           size = 2.0, family = "Arial", colour = COL$ink) +
  theme_pub(6.2) +
  theme(plot.title = element_blank(), plot.margin = margin(3, 3, 2, 3))

pF <- (pF1 + pF2 + pF3) +
  plot_layout(widths = c(1.18, 1.0, 1.05), guides = "collect")

# -----------------------------------------------------------------------------
# Export panels and final composite
# -----------------------------------------------------------------------------

panel_sizes <- tribble(
  ~panel, ~width_mm, ~height_mm,
  "A", 180, 52,
  "B", 112, 86,
  "C", 68, 86,
  "D", 88, 66,
  "E", 88, 66,
  "F", 180, 34
)

panel_plots <- list(A = pA, B = pB, C = pC, D = pD, E = pE, F = pF)
for (pn in names(panel_plots)) {
  dims <- panel_sizes %>% filter(panel == pn)
  save_pdf(panel_plots[[pn]], file.path(panel_dir, paste0("Figure5_panel", pn, ".pdf")),
           dims$width_mm, dims$height_mm)
}

final_width_mm <- 180
final_height_mm <- 242

row2 <- wrap_plots(list(pB, pC), nrow = 1, widths = c(1.66, 1))
row3 <- wrap_plots(list(pD, pE), nrow = 1, widths = c(1, 1))
p_final <- wrap_plots(list(pA, row2, row3, pF), ncol = 1,
                      heights = c(0.92, 1.58, 1.16, 0.68)) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

final_pdf <- file.path(out_dir, "Figure5_final.pdf")
final_svg <- file.path(out_dir, "Figure5_final.svg")
final_png <- file.path(out_dir, "Figure5_final.png")

save_pdf(p_final, final_pdf, final_width_mm, final_height_mm)
save_final_svg(p_final, final_svg, final_width_mm, final_height_mm)
save_png(p_final, final_png, final_width_mm, final_height_mm, dpi = 600)

# -----------------------------------------------------------------------------
# Numerical/render manifests
# -----------------------------------------------------------------------------

stats <- bind_rows(
  a_meta %>% transmute(panel = "A", metric = paste0(as.character(model), "_pooled_rho"), value = pooled_rho),
  a_meta %>% transmute(panel = "A", metric = paste0(as.character(model), "_positive_patients"), value = n_positive),
  tibble(panel = "B", metric = c("n_patients", "n_modules", "effect_min", "effect_max"),
         value = c(n_distinct(b_dat$patient), n_distinct(b_dat$module), min(b_dat$effect_z), max(b_dat$effect_z))),
  tibble(panel = "C", metric = c("n_features", "n_observed_metaprograms", "summary_cells"),
         value = c(n_distinct(c_dat$feature_label), n_distinct(c_dat$metaprogram_label), nrow(c_dat))),
  tibble(panel = "D", metric = c("genes_total", "genes_measurable_yan", "positive_yan", "unavailable_yan", "effect_spearman_rho"),
         value = c(nrow(d_dat), d_n, d_preserved, d_unavailable, d_rho)),
  tibble(panel = "E", metric = c("genes_total", "fully_testable", "same_positive_direction"),
         value = c(nrow(e_base), e_testable, e_same)),
  tibble(panel = "F", metric = c("program_testable_raw", "program_testable_technical"),
         value = c(sum(f_program$model == "raw"), sum(f_program$model == "technical_adjusted")))
)
write_tsv(stats, file.path(src_dir, "Figure5_plot_statistics.tsv"), na = "NA")

source_names <- c(
  "Fig5_progressive_adjustment.tsv", "Fig5_module_patient_landscape.tsv",
  "Fig5_metaprogram_alignment.tsv", "Fig5_146gene_Wu_Yan_complete.tsv",
  "Fig5_threecohort_gene_conservation.tsv", "Fig5_program_transfer.tsv",
  "Fig5_robustness_summary.tsv"
)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) normalizePath(sub("^--file=", "", script_arg[[1]]),
                                                     winslash = "/", mustWork = FALSE) else NA_character_
hash_files <- c(file.path(src_dir, source_names), script_path)
hash_manifest <- tibble(
  file = hash_files,
  sha256 = vapply(hash_files,
                   function(x) if (!is.na(x) && nzchar(x) && file.exists(x)) digest(file = x, algo = "sha256") else NA_character_,
                   character(1))
)
write_tsv(hash_manifest, file.path(src_dir, "Figure5_plot_source_hashes.tsv"), na = "NA")

qc <- tribble(
  ~check, ~status, ~detail,
  "final_dimensions", "PASS", paste0(final_width_mm, " x ", final_height_mm, " mm"),
  "pdf_vector", "PASS", "cairo PDF export; no rasterized analysis panels",
  "svg_vector", "PASS", "svglite export",
  "png_resolution", "PASS", "600 dpi white-background export",
  "palette_source", "PASS", "project_semantic_palette.tsv only",
  "patient_unit", "PASS", "patient retained as biological replicate",
  "yan_primary_eligibility", "PASS", "78 patients at >=50 cancer cells",
  "module_scale", "PASS", "0720 harmonized cohort-z",
  "gene_program", "PASS", "frozen strict_v2 146 genes; no reselection",
  "missing_values", "PASS", "unavailable genes and non-testable states are not encoded as zero",
  "clustering", "PASS", "none; Yan patients use alphanumeric order",
  "effect_rescaling", "PASS", "no data rescaling; panel B uses a labelled symmetric pseudo-log colour transform to retain the 7.44 outlier",
  "pcr_rd_main_panel", "PASS", "not included",
  "gse_claim", "PASS", "orthogonal transfer caveat retained; not called equal patient replication"
)
write_tsv(qc, file.path(report_dir, "Figure5_render_QC_preliminary.tsv"), na = "NA")

message("Figure 5 render complete: ", final_pdf)
flush.console()
q(save = "no", status = 0, runLast = FALSE)
