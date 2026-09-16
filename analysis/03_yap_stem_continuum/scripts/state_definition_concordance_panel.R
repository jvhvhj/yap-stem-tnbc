#!/usr/bin/env Rscript
# Purpose: Figure 3E concordance label correction
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 3 and Supplementary Fig. S2.

options(stringsAsFactors = FALSE, width = 220, warn = 1)

root <- "."
out_dir <- file.path(root, "Figure3E_kappa_concordance_patch")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

.libPaths(c(file.path(root, "0703_rebuild/R_library"), .libPaths()))
suppressWarnings(suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(grid)
}))

source_path <- file.path(root, "Figure2_panelE_source.tsv")
palette_path <- file.path(root, "project_semantic_palette.tsv")
stopifnot(file.exists(source_path), file.exists(palette_path))

source_data <- read.delim(source_path, check.names = FALSE, fileEncoding = "UTF-8")
palette_tbl <- read.delim(palette_path, check.names = FALSE, fileEncoding = "UTF-8")
pal <- setNames(palette_tbl$hex, palette_tbl$token)

patient_order <- c(
  "CID4465", "CID4495", "CID44971", "CID44991",
  "CID4513", "CID4515", "CID4523", "CID3963"
)
column_names <- c(patient_order, "Overall")
state_levels <- c("Low", "Intermediate", "High")
definition_levels <- c(
  "Pooled tertile", "Within-patient tertile", "Patient-standardized tertile"
)
comparison_levels <- c(
  "Pooled vs Within-patient",
  "Pooled vs Patient-standardized",
  "Within-patient vs Patient-standardized"
)

composition_df <- source_data[source_data$record_type == "state_composition", ]
concordance_df <- source_data[source_data$record_type == "pairwise_concordance", ]
stopifnot(
  nrow(composition_df) == 81L,
  nrow(concordance_df) == 27L,
  setequal(unique(composition_df$patient), column_names),
  setequal(unique(concordance_df$patient), column_names),
  setequal(unique(composition_df$definition), definition_levels),
  setequal(unique(concordance_df$comparison), comparison_levels),
  all(is.finite(concordance_df$cohen_kappa)),
  all(is.finite(concordance_df$high_jaccard))
)

# This script deliberately does not recalculate scores, state assignments, kappa,
# Jaccard, or any biological result. It only re-renders the authoritative source.
column_x <- setNames(seq_along(column_names), column_names)
definition_y <- c(
  "Pooled tertile" = 6.0,
  "Within-patient tertile" = 5.0,
  "Patient-standardized tertile" = 4.0
)
comparison_y <- c(
  "Pooled vs Within-patient" = 2.45,
  "Pooled vs Patient-standardized" = 1.45,
  "Within-patient vs Patient-standardized" = 0.45
)

composition_df$patient_x <- unname(column_x[composition_df$patient])
composition_df$definition_y <- unname(definition_y[composition_df$definition])
composition_df$definition <- factor(composition_df$definition, levels = definition_levels)
composition_df$patient <- factor(composition_df$patient, levels = column_names)
composition_df$state <- factor(composition_df$state, levels = state_levels)
composition_df <- composition_df[order(composition_df$definition, composition_df$patient,
                                       composition_df$state), ]

composition_rows <- split(composition_df, interaction(
  composition_df$definition, composition_df$patient, drop = TRUE
))
composition_df <- do.call(rbind, lapply(composition_rows, function(d) {
  d <- d[match(state_levels, as.character(d$state)), ]
  cum <- c(0, cumsum(d$state_fraction))
  d$xmin <- d$patient_x - 0.43 + 0.86 * cum[seq_along(state_levels)]
  d$xmax <- d$patient_x - 0.43 + 0.86 * cum[seq_along(state_levels) + 1L]
  d$ymin <- d$definition_y - 0.33
  d$ymax <- d$definition_y + 0.33
  d
}))

concordance_df$patient_x <- unname(column_x[concordance_df$patient])
concordance_df$comparison_y <- unname(comparison_y[concordance_df$comparison])
concordance_df$comparison <- factor(concordance_df$comparison, levels = comparison_levels)
concordance_df$patient <- factor(concordance_df$patient, levels = column_names)
concordance_df <- concordance_df[order(concordance_df$comparison, concordance_df$patient), ]
concordance_df$xmin <- concordance_df$patient_x - 0.43
concordance_df$xmax <- concordance_df$patient_x + 0.43
concordance_df$ymin <- concordance_df$comparison_y - 0.36
concordance_df$ymax <- concordance_df$comparison_y + 0.36
concordance_df$kappa_label <- sprintf("κ=%.2f", concordance_df$cohen_kappa)
concordance_df$text_class <- ifelse(concordance_df$cohen_kappa >= 0.60,
                                    "dark_tile", "light_tile")

# Freeze all displayed numbers against the previously rendered labels.
expected_kappa_labels <- c(
  "κ=0.66", "κ=0.27", "κ=0.83", "κ=0.89", "κ=0.47", "κ=0.34",
  "κ=0.25", "κ=0.27", "κ=0.58",
  "κ=0.66", "κ=0.27", "κ=0.83", "κ=0.89", "κ=0.47", "κ=0.34",
  "κ=0.25", "κ=0.27", "κ=0.58",
  "κ=0.95", "κ=0.94", "κ=0.97", "κ=0.96", "κ=0.98", "κ=0.96",
  "κ=0.99", "κ=0.97", "κ=0.97"
)
stopifnot(identical(concordance_df$kappa_label, expected_kappa_labels))

col_low <- unname(pal[["Low"]])
col_mid <- unname(pal[["Intermediate"]])
col_high <- unname(pal[["High"]])
col_neutral <- unname(pal[["Neutral"]])
col_raw <- unname(pal[["Raw"]])
stopifnot(all(nzchar(c(col_low, col_mid, col_high, col_neutral, col_raw))))

font_family <- "sans"
text_col <- "#242424"
axis_col <- "#333333"
struct_dark <- "#74787B"
state_cols <- c(Low = col_low, Intermediate = col_neutral, High = col_high)

theme_submission <- function(base_size = 7.0) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      plot.title = element_text(face = "bold", size = 8.5, hjust = 0,
                                colour = text_col, margin = margin(b = 5.0)),
      axis.title = element_text(size = 7.0, colour = text_col),
      axis.text = element_text(size = 6.6, colour = text_col),
      axis.line = element_line(linewidth = 0.32, colour = axis_col),
      axis.ticks = element_line(linewidth = 0.28, colour = axis_col),
      legend.title = element_text(size = 6.7, colour = text_col),
      legend.text = element_text(size = 6.5, colour = text_col),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(5.0, 5.0, 5.0, 6.0, unit = "mm")
    )
}

column_bg <- data.frame(
  xmin = seq_along(column_names) - 0.47,
  xmax = seq_along(column_names) + 0.47,
  fill = ifelse(seq_along(column_names) %% 2 == 0, "#FAFAFA", "white")
)
n_track <- unique(composition_df[, c("patient", "patient_x", "malignant_cells")])
n_track <- n_track[match(column_names, as.character(n_track$patient)), ]
n_track$label <- format(n_track$malignant_cells, big.mark = ",",
                        scientific = FALSE, trim = TRUE)
legend_dummy <- data.frame(
  state = factor(state_levels, levels = state_levels), x = 1, y = 1
)

p_e <- ggplot() +
  geom_rect(data = column_bg,
            aes(xmin = xmin, xmax = xmax, ymin = 0.02, ymax = 7.12),
            fill = column_bg$fill, colour = NA) +
  geom_rect(data = composition_df[composition_df$state == "Low", ],
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = col_low, colour = "white", linewidth = 0.18) +
  geom_rect(data = composition_df[composition_df$state == "Intermediate", ],
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = col_neutral, colour = "white", linewidth = 0.18) +
  geom_rect(data = composition_df[composition_df$state == "High", ],
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = col_high, colour = "white", linewidth = 0.18) +
  geom_segment(aes(x = 0.53, xend = 9.47, y = 3.23, yend = 3.23),
               linewidth = 0.55, colour = struct_dark) +
  geom_rect(data = concordance_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = cohen_kappa),
            colour = "white", linewidth = 0.45) +
  geom_text(data = concordance_df[concordance_df$text_class == "light_tile", ],
            aes(x = patient_x, y = (ymin + ymax) / 2, label = kappa_label),
            colour = "#1F2933", size = 2.43,
            family = font_family, fontface = "bold") +
  geom_text(data = concordance_df[concordance_df$text_class == "dark_tile", ],
            aes(x = patient_x, y = (ymin + ymax) / 2, label = kappa_label),
            colour = "#FFFFFF", size = 2.43,
            family = font_family, fontface = "bold") +
  geom_text(data = n_track, aes(x = patient_x, y = 6.88, label = label),
            size = 2.27, family = font_family, colour = text_col) +
  annotate("text", x = 0.49, y = 6.88, label = "n", hjust = 1,
           size = 2.27, family = font_family, fontface = "bold", colour = text_col) +
  geom_point(data = legend_dummy, aes(x = x, y = y, colour = state), alpha = 0) +
  scale_colour_manual(
    values = state_cols,
    breaks = state_levels, name = NULL,
    guide = guide_legend(
      direction = "horizontal", nrow = 1,
      override.aes = list(alpha = 1, shape = 15, size = 4.5)
    )
  ) +
  scale_fill_gradientn(
    colours = c(col_neutral, col_raw, col_low),
    values = c(0, 0.58, 1), limits = c(0, 1),
    breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"),
    oob = squish, name = "State-definition agreement (Cohen's κ)",
    guide = guide_colourbar(
      direction = "horizontal", barheight = unit(2.8, "mm"),
      barwidth = unit(24, "mm"), title.position = "top", title.hjust = 0.5,
      ticks.colour = "white", frame.colour = NA
    )
  ) +
  scale_x_continuous(
    breaks = seq_along(column_names), labels = column_names,
    limits = c(0.48, 9.52), position = "top", expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    breaks = c(6.0, 5.0, 4.0, 2.45, 1.45, 0.45),
    labels = c(definition_levels, comparison_levels),
    limits = c(0.02, 7.14), expand = expansion(mult = 0)
  ) +
  labs(title = "State-definition concordance", x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_submission(6.8) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 6.5, angle = 38, hjust = 0, vjust = 0.5,
                               margin = margin(b = 2.3, unit = "mm")),
    axis.text.y = element_text(size = 6.5, hjust = 1,
                               margin = margin(r = 2.5, unit = "mm")),
    legend.position = "bottom", legend.box = "horizontal",
    legend.spacing.x = unit(4.0, "mm"),
    legend.margin = margin(2.0, 0, 0, 0, unit = "mm")
  )

pdf_path <- file.path(out_dir, "Figure3E_contrast_updated.pdf")
svg_path <- file.path(out_dir, "Figure3E_contrast_updated.svg")

cairo_pdf(pdf_path, width = 180 / 25.4, height = 116 / 25.4,
          family = font_family, bg = "white", onefile = TRUE)
print(p_e)
dev.off()

svglite::svglite(svg_path, width = 180 / 25.4, height = 116 / 25.4,
                 bg = "white", system_fonts = list(sans = "Arial"))
print(p_e)
dev.off()

message("Rendered Figure3E with high-contrast Cohen's kappa labels; composite Figure3 was not touched.")
