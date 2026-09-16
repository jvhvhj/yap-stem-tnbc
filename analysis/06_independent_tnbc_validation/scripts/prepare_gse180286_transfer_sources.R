#!/usr/bin/env Rscript
# Purpose: Prepare GSE180286 transfer source tables
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 5; orthogonal transfer in Fig. 4.

options(stringsAsFactors = FALSE, scipen = 999, width = 240, encoding = "UTF-8")

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
source_dir <- file.path(root, "Figure4G_GSE180286_program_transfer")
out_dir <- file.path(root, "Figure4F_GSE180286_program_transfer")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

.libPaths(c(file.path(root, "0703_rebuild/R_library"), .libPaths()))
suppressPackageStartupMessages(library(digest))

read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
write_tsv <- function(x, path) write.table(
  x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA", fileEncoding = "UTF-8"
)
sha256 <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)

source_paths <- list(
  transfer = file.path(source_dir, "Figure4G_program_gene_transfer.tsv"),
  waterfall = file.path(source_dir, "Figure4G_gene_effect_waterfall.tsv"),
  enrichment = file.path(source_dir, "Figure4G_signature_enrichment.tsv"),
  transfer_audit = file.path(source_dir, "Figure4G_transfer_audit.tsv"),
  metadata_audit = file.path(root, "0718_GSE180286_primary_LN_singlecell_closure/tables/GSE180286_metadata_audit.tsv"),
  metadata_script = file.path(root, "0718_GSE180286_primary_LN_singlecell_closure/scripts/01_GSE180286_closure.R"),
  tnbc_counts = file.path(root, "0718_GSE180286_primary_LN_singlecell_closure/tables/GSE180286_TNBC_malignant_cell_counts.tsv")
)
missing <- names(source_paths)[!vapply(source_paths, file.exists, logical(1))]
if (length(missing)) stop("Missing required source(s): ", paste(missing, collapse = ", "))

transfer <- read_tsv(source_paths$transfer)
waterfall <- read_tsv(source_paths$waterfall)
enrichment <- read_tsv(source_paths$enrichment)
transfer_audit <- read_tsv(source_paths$transfer_audit)
metadata <- read_tsv(source_paths$metadata_audit)
tnbc_counts <- read_tsv(source_paths$tnbc_counts)

stopifnot(
  nrow(transfer) == 146L,
  sum(transfer$detectable_in_GSE180286) == 141L,
  nrow(waterfall) == 141L,
  sum(waterfall$GSE180286_High_vs_Other_effect > 0) == 132L,
  sum(waterfall$GSE180286_High_vs_Other_effect == 0) == 1L,
  sum(waterfall$GSE180286_High_vs_Other_effect < 0) == 8L,
  nrow(enrichment) == 33209L,
  identical(unique(enrichment$enrichment_NES), 1.62757269669218),
  identical(unique(enrichment$enrichment_BH_FDR), 0.00000000331676831473242)
)

## Terminology-only update. No value, order, signature, or ranking changes.
transfer$transfer_class[transfer$transfer_class == "Discordant or zero"] <- "Not preserved / near-zero"
waterfall$transfer_class[waterfall$transfer_class == "Discordant or zero"] <- "Not preserved / near-zero"

write_tsv(transfer, file.path(out_dir, "Figure4F_program_gene_transfer.tsv"))
write_tsv(waterfall, file.path(out_dir, "Figure4F_gene_effect_waterfall.tsv"))
write_tsv(enrichment, file.path(out_dir, "Figure4F_signature_enrichment.tsv"))

## Resolve the complete official TNBC sample manifest, including zero-cell samples.
prov <- metadata[metadata$subtype == "TNBC", , drop = FALSE]
prov$tissue_origin <- ifelse(prov$site == "Primary", "Primary", prov$site_status)
prov$included_in_transfer_effect <- prov$local_sample %in% c(
  "GSM5457211_E2020-1", "GSM5457212_E2020-2"
)
prov$transfer_role <- ifelse(
  prov$included_in_transfer_effect,
  "Included in pooled High-minus-Other transfer effect",
  ifelse(
    prov$local_sample == "GSM5457210_D2020-3",
    "Excluded: only one TNBC malignant cell",
    "Excluded: no TNBC malignant cells in analysis object"
  )
)
prov$patient_mapping_status <- "RESOLVED"
prov$tissue_mapping_status <- "RESOLVED"
prov$mapping_conclusion <- ifelse(
  grepl("E2020", prov$local_sample),
  "E2020 samples map to patient P5; E2020-1 is Primary and E2020-2/E2020-3 are LN+",
  "D2020 samples map to patient P4; D2020-1 is Primary and D2020-2/D2020-3 are LN-"
)
prov$earlier_note_reconciliation <- ifelse(
  prov$patient == "P5",
  "The 5,297-cell Primary and 162-cell LN+ samples are two samples from the same patient P5",
  "Earlier two-patient notes count the single P4 cell in the old TNBC object; P4 is not used in the transfer effect"
)
prov$authoritative_metadata_source <- source_paths$metadata_audit
prov$authoritative_metadata_sha256 <- sha256(source_paths$metadata_audit)
prov$authoritative_mapping_script <- source_paths$metadata_script
prov$authoritative_mapping_script_sha256 <- sha256(source_paths$metadata_script)

prov <- prov[, c(
  "gsm", "local_sample", "patient", "subtype", "tissue_origin", "site", "node_index",
  "site_status", "site_status_basis", "n_malignant_cells", "present_in_old_5469_object",
  "included_in_transfer_effect", "transfer_role", "patient_mapping_status",
  "tissue_mapping_status", "mapping_conclusion", "earlier_note_reconciliation",
  "authoritative_metadata_source", "authoritative_metadata_sha256",
  "authoritative_mapping_script", "authoritative_mapping_script_sha256"
)]
write_tsv(prov, file.path(out_dir, "GSE180286_sample_patient_provenance_audit.tsv"))

inc <- prov[prov$included_in_transfer_effect, , drop = FALSE]
stopifnot(
  nrow(inc) == 2L,
  identical(inc$local_sample, c("GSM5457211_E2020-1", "GSM5457212_E2020-2")),
  all(inc$patient == "P5"),
  identical(inc$tissue_origin, c("Primary", "LN+")),
  identical(as.integer(inc$n_malignant_cells), c(5297L, 162L))
)

extra_audit <- data.frame(
  audit_field = c(
    "final_panel_number",
    "transfer_sample_primary",
    "transfer_sample_ln_positive",
    "included_transfer_samples_same_patient",
    "included_transfer_patient",
    "legacy_two_patient_note_reconciliation",
    "positive_effect_genes",
    "zero_effect_genes",
    "negative_effect_genes",
    "median_transferred_effect",
    "terminology_update"
  ),
  value = c(
    "Main Figure 4F",
    "GSM5457211_E2020-1 | P5 | Primary | 5297 TNBC malignant cells",
    "GSM5457212_E2020-2 | P5 | LN+ | 162 TNBC malignant cells",
    "TRUE",
    "P5",
    "Old object contained P4 (1 cell) and P5; exact transfer effect excluded P4 and used the two P5 samples only",
    "132",
    "1",
    "8",
    format(median(waterfall$GSE180286_High_vs_Other_effect), digits = 12),
    "Discordant or zero renamed Not preserved / near-zero; no numerical change"
  ),
  status = "PASS",
  source_file = c(
    source_paths$transfer_audit,
    source_paths$metadata_audit,
    source_paths$metadata_audit,
    source_paths$metadata_audit,
    source_paths$metadata_audit,
    source_paths$metadata_audit,
    source_paths$waterfall,
    source_paths$waterfall,
    source_paths$waterfall,
    source_paths$waterfall,
    source_paths$waterfall
  ),
  source_sha256 = c(
    sha256(source_paths$transfer_audit),
    rep(sha256(source_paths$metadata_audit), 5),
    rep(sha256(source_paths$waterfall), 5)
  ),
  stringsAsFactors = FALSE
)
transfer_audit <- rbind(transfer_audit, extra_audit)
write_tsv(transfer_audit, file.path(out_dir, "Figure4F_transfer_audit.tsv"))

cat("Figure 4F sources prepared without rerunning transfer analysis.\n")
cat("Included samples: ", paste(inc$local_sample, inc$patient, inc$tissue_origin, inc$n_malignant_cells, sep = " | ", collapse = "; "), "\n", sep = "")
