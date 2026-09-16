#!/usr/bin/env Rscript
# Purpose: Finalize score-independent companion tables
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 4 and Supplementary Fig. S1.

options(stringsAsFactors = FALSE, width = 220, timeout = 900)
set.seed(20260713)
root <- normalizePath(".", winslash = "/", mustWork = TRUE)
src <- file.path(root, "0710_final_evidence_rebuild")
out <- file.path(root, "0713_score_independent_rebuild")
tab <- file.path(out, "tables")
log_dir <- file.path(out, "logs")
dir.create(tab, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(src, "vendor"), ".software/r-library",
            file.path(root, "0703_rebuild/R_library"), .libPaths()))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(AnnotationDbi); library(DOSE); library(hexbin)
})
log_con <- file(file.path(log_dir, "01_finalize_companion_tables.log"), "wt")
sink(log_con, type = "output", split = TRUE); sink(log_con, type = "message", append = TRUE)

write_tsv <- function(x, name, allow_empty = FALSE) {
  if (is.null(x) || (!allow_empty && nrow(x) == 0L)) stop("Empty table: ", name)
  write.table(x, file.path(tab, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  cat("Wrote", name, "(", nrow(x), "rows)\n")
}
read_old <- function(name) read.delim(file.path(src, "tables_final", name), check.names = FALSE)
read_new <- function(name) read.delim(file.path(tab, name), check.names = FALSE)

# Figure 1 unified plotting data and CopyKAT supplement ----------------------
atlas <- read.csv(file.path(root, "0703_rebuild/workstation_results/Wu2021_full_atlas_umap_metadata.csv"), check.names = FALSE)
atlas$identity <- "Other cells"
atlas$identity[atlas$celltype_major == "Cancer Epithelial"] <- "Other cancer epithelial"
atlas$identity[atlas$celltype_major == "Normal Epithelial"] <- "Normal epithelial"
atlas$identity[atlas$subtype == "TNBC" & atlas$celltype_major == "Cancer Epithelial"] <- "TNBC malignant epithelial"
atlas_plot <- atlas %>% transmute(panel = "A_atlas", record_type = "cell", cell_id,
  patient = as.character(patient), subtype = as.character(subtype), celltype_major = as.character(celltype_major),
  identity, UMAP_1, UMAP_2)
contribution <- read_old("Fig1_patient_contribution.tsv") %>%
  transmute(panel = "A_patient_contribution", record_type = "patient", patient, n_cells, fraction)
copy_sum <- read_old("Fig1_copykat_summary.tsv") %>%
  mutate(panel = "B_copykat_summary", record_type = "copykat")
copy_pt <- read_old("Fig1_copykat_patient.tsv")
write_tsv(copy_pt, "Fig1_copykat_supp_table.tsv")
deg1 <- read_old("Fig1_malignant_vs_normal_DEG.tsv")
prolif <- c("MKI67","UBE2C","CENPF","CDK1","BIRC5","TK1","PTTG1","NUSAP1","TOP2A","STMN1")
normal_markers <- c("PIGR","AGR3","TFF1","PGR","PIP","MUCL1","LTF","SLPI","KRT14","KRT17","EPCAM","KRT19")
deg1$gene_class <- ifelse(deg1$gene %in% prolif, "Malignant-high proliferation",
                          ifelse(deg1$gene %in% normal_markers, "Normal-high epithelial/luminal", "Other"))
deg1p <- deg1 %>% mutate(panel = "C_volcano", record_type = "gene")
hm1 <- read_old("Fig1_malignant_vs_normal_marker_heatmap_matrix.tsv") %>%
  mutate(panel = "D_marker_heatmap", record_type = "pseudobulk")
prog1 <- read_old("Fig1_malignant_vs_normal_program_summary.tsv") %>%
  mutate(panel = "E_program_summary", record_type = "program")
fig1_all <- bind_rows(atlas_plot, contribution, copy_sum, deg1p, hm1, prog1)
write_tsv(fig1_all, "Fig1_plot_data_all.tsv")

# Figure 2 exact requested companion tables ---------------------------------
meta <- read.csv(file.path(root, "0703_rebuild/continuous_state_main_figures/Wu2021_continuous_state_cell_metadata.csv"), check.names = FALSE)
raw <- read_old("Fig2_patient_rho.tsv") %>% mutate(estimate_type = "Raw")
adj <- read_old("Fig2_adjusted_residual_rho.tsv") %>% mutate(estimate_type = "Adjusted residual")
rho <- bind_rows(raw, adj) %>% dplyr::select(patient, n_cells, estimate_type, rho, ci_low, ci_high, n_boot)
write_tsv(rho, "Fig2_raw_adjusted_rho.tsv")
dec <- read_old("Fig2_decile_trend.tsv")
pooled <- read_old("Fig2_decile_trend_pooled.tsv") %>% mutate(patient = "Across-patient summary", n_cells = NA_integer_, median_YAP = NA_real_, median_Stemness = patient_median, q25_Stemness = patient_q25, q75_Stemness = patient_q75)
write_tsv(bind_rows(dec, pooled[, names(dec)]), "Fig2_decile_trend.tsv")

hex <- bind_rows(lapply(sort(unique(meta$patient)), function(pt) {
  d <- meta[meta$patient == pt, ]
  hb <- hexbin(d$YAP_score, d$Stemness_score, xbins = 38, IDs = TRUE)
  xy <- hcell2xy(hb)
  data.frame(patient = pt, hex_id = seq_along(hb@count), x = xy$x, y = xy$y,
             count = hb@count, log10_count = log10(hb@count + 1), n_cells = nrow(d))
}))
write_tsv(hex, "Fig2_hexbin_density_data.tsv")

# Restore six strictly trend-stable score-independent representative genes.
score_genes <- read_new("score_gene_list.tsv")$gene
old_core <- read_old("Fig5_core_gene_matrix.tsv")
core_ind <- old_core %>% filter(!gene %in% score_genes,
  gene %in% c("TACSTD2","LAMB3","LAMC2","CXCL2","MCL1","BIRC3"))
write_tsv(core_ind, "Fig3_core_gene_dotplot.tsv")

# Add current HDO enrichment to both completed enrichment tables. ------------
orgdb <- AnnotationDbi::loadDb(file.path(src, "vendor/org.Hs.eg.db/extdata/org.Hs.eg.sqlite"))
# Register the frozen official HDO ontology database before calling enrichDO().
# This mirrors the main reanalysis and avoids machine-specific Windows caches.
hdo_path <- file.path(out, "references", "HDO.sqlite")
if (!file.exists(hdo_path)) stop("Frozen HDO database is missing: ", hdo_path)
hdo_db <- AnnotationDbi::loadDb(hdo_path)
yulab.utils::update_cache_item(".GOSemSimEnv", setNames(list(hdo_db), ".onto_HDO"))
deg_all <- read_new("Fig4_DEG_all.tsv")
universe_symbols <- unique(deg_all$gene)
map <- AnnotationDbi::select(orgdb, keys = universe_symbols, keytype = "SYMBOL", columns = "ENTREZID") %>%
  filter(!is.na(ENTREZID)) %>% distinct(SYMBOL, ENTREZID)
universe_entrez <- unique(map$ENTREZID)
entrez_to_symbol <- setNames(map$SYMBOL, map$ENTREZID)
run_hdo <- function(query_symbols) {
  q <- unique(map$ENTREZID[map$SYMBOL %in% query_symbols])
  x <- DOSE::enrichDO(q, universe = universe_entrez, ont = "HDO", pAdjustMethod = "BH",
                      pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
  d <- as.data.frame(x)
  if (!nrow(d)) return(data.frame())
  if ("geneID" %in% names(d)) d$geneID <- vapply(strsplit(d$geneID, "/", fixed = TRUE),
    function(ids) paste(unique(na.omit(entrez_to_symbol[ids])), collapse = "/"), character(1))
  d$database <- "Disease Ontology"
  keep <- intersect(c("database","ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count"), names(d))
  d[, keep, drop = FALSE]
}
for (nm in c("Fig4_enrichment_all_genes.tsv", "Fig4_enrichment_score_independent.tsv")) {
  e <- read_new(nm) %>% filter(database != "Disease Ontology")
  query <- if (grepl("score_independent", nm))
    deg_all$gene[deg_all$significant_high & !deg_all$score_gene] else deg_all$gene[deg_all$significant_high]
  hdo <- tryCatch(run_hdo(query), error = function(err) {
    cat("HDO error for", nm, ":", conditionMessage(err), "\n"); data.frame()
  })
  write_tsv(bind_rows(e, hdo) %>% arrange(p.adjust), nm)
}

cat("Companion-table finalization finished:", format(Sys.time()), "\n")
sink(type = "message"); sink(type = "output"); close(log_con)
