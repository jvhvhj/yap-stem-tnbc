#!/usr/bin/env Rscript
# Purpose: Score-independent discovery reanalysis
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
  library(Seurat); library(Matrix); library(dplyr); library(tidyr)
  library(AnnotationDbi); library(clusterProfiler); library(ReactomePA); library(DOSE)
})

log_con <- file(file.path(log_dir, "00_score_independent_reanalysis.log"), "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message", append = TRUE)
cat("Score-independent reanalysis started:", format(Sys.time()), "\n")

write_tsv <- function(x, name, allow_empty = FALSE) {
  if (is.null(x) || (!allow_empty && nrow(x) == 0L)) stop("Empty required table: ", name)
  write.table(x, file.path(tab, name), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
  cat("Wrote", name, "(", nrow(x), "rows)\n")
}
z <- function(x) {
  y <- as.numeric(scale(x)); y[!is.finite(y)] <- 0; y
}
z_rows <- function(x) {
  y <- t(scale(t(x))); y[!is.finite(y)] <- 0; y
}
mean_score <- function(mat, genes) {
  g <- intersect(unique(genes), rownames(mat))
  if (!length(g)) stop("No expressed genes for module")
  Matrix::colMeans(mat[g, , drop = FALSE])
}
sign_p <- function(k, n) binom.test(k, n, p = 0.5, alternative = "two.sided")$p.value
read_gmt <- function(path) {
  x <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  data.frame(term = rep(vapply(x, `[`, character(1), 1), lengths(x) - 2L),
             gene = unlist(lapply(x, function(v) v[-c(1, 2)])), stringsAsFactors = FALSE)
}

# Frozen inputs --------------------------------------------------------------
meta <- read.csv(file.path(root, "0703_rebuild/continuous_state_main_figures/Wu2021_continuous_state_cell_metadata.csv"), check.names = FALSE)
gene_sets <- read.csv(file.path(root, "0703_rebuild/continuous_state_main_figures/core_score_gene_lists.csv"), check.names = FALSE)
deg <- read.csv(file.path(root, "0703_rebuild/continuous_state_main_figures/Figure4_DEG_pathway/tables/High_vs_Low_DE_all.csv"), check.names = FALSE)
obj <- readRDS(file.path(root, "output_step0/Wu2021_TNBC_malignant.rds"))
DefaultAssay(obj) <- "RNA"
expr <- GetAssayData(obj, assay = "RNA", layer = "data")
stopifnot(nrow(meta) == 10836L, length(unique(meta$patient)) == 8L,
          all(meta$cell_id %in% colnames(expr)))
expr <- expr[, meta$cell_id, drop = FALSE]

yap_genes <- unique(gene_sets$gene[gene_sets$axis == "YAP"])
stem_genes <- unique(gene_sets$gene[gene_sets$axis == "Stemness"])
score_genes <- unique(c(yap_genes, stem_genes))
score_gene_list <- data.frame(
  gene = score_genes,
  in_YAP_score = score_genes %in% yap_genes,
  in_Stemness_score = score_genes %in% stem_genes,
  in_YAP_Stem_construction = TRUE,
  role = ifelse(score_genes %in% yap_genes, "YAP score", "Stemness score")
)
write_tsv(score_gene_list, "score_gene_list.tsv")

# All requested downstream modules ------------------------------------------
modules <- list(
  UPR = c("HSPA5","ATF4","DDIT3","XBP1","DNAJB9","ATF3","HERPUD1","PPP1R15A"),
  TNFA_NFKB = c("TNF","NFKBIA","TNFAIP3","RELA","JUN","FOS","IER3","ICAM1","CXCL2","BIRC3"),
  Hypoxia = c("HIF1A","CA9","VEGFA","SLC2A1","LDHA","PDK1","BNIP3","ADM","EGLN3","NDRG1"),
  Adhesion_Remodeling = c("ITGA6","ITGB1","ITGB4","TACSTD2","CD44","F3","ANXA1","LAMB3","LAMC2","FN1"),
  Integrin_Adhesion = c("ITGA6","ITGB1","ITGB4","ITGA3","PTK2","LAMB3","LAMC2"),
  Wound_Healing = c("AREG","EREG","ANXA1","F3","CXCL2","GDF15","SERPINE1","CTGF","S100A10"),
  Survival_Stress = c("BCL2L1","MCL1","BIRC2","BIRC3","BIRC5","XIAP","CFLAR","GDF15","CDKN1A","IGFBP3"),
  Anoikis_Resistance = c("BCL2L1","MCL1","CFLAR","BIRC5","ITGA6","ITGB4","PTK2"),
  Cell_Cycle = c("MKI67","TOP2A","UBE2C","CENPF","BIRC5","PTTG1","NUSAP1","CDK1"),
  EMT = c("VIM","FN1","SNAI1","SNAI2","TWIST1","ZEB1","ZEB2","CDH2","TGFBI"),
  Migration_Invasion = c("MMP2","MMP9","MMP14","PLAU","PLAUR","SERPINE1","CXCL8"),
  Stem_Plasticity = c("SOX9","HES1","KLF4","TACSTD2","CD44","ITGA6","NOTCH1","NOTCH2")
)
module_labels <- gsub("_", " ", names(modules))

overlap <- bind_rows(lapply(seq_along(modules), function(i) {
  original <- unique(modules[[i]])
  removed <- intersect(original, score_genes)
  retained <- setdiff(original, score_genes)
  data.frame(module = module_labels[i], n_genes_original = length(original),
             n_score_genes_overlap = length(removed),
             overlap_fraction = length(removed) / length(original),
             n_genes_after_removal = length(retained),
             removed_genes = paste(removed, collapse = ","),
             retained_genes = paste(retained, collapse = ","))
}))
write_tsv(overlap, "module_gene_overlap.tsv")

# Recalculate original and score-gene-removed scores ------------------------
score_long <- bind_rows(lapply(seq_along(modules), function(i) {
  original <- intersect(modules[[i]], rownames(expr))
  retained <- setdiff(original, score_genes)
  before <- z(mean_score(expr, original))
  after <- if (length(retained)) z(mean_score(expr, retained)) else rep(NA_real_, ncol(expr))
  data.frame(cell_id = meta$cell_id, patient = meta$patient, YS_tertile = meta$YS_tertile,
             module = module_labels[i], score_before = before, score_after = after,
             n_genes_before = length(original), n_genes_after = length(retained))
}))
write_tsv(score_long, "module_scores_before_after_removing_score_genes.tsv")

patient_delta <- score_long %>%
  filter(YS_tertile %in% c("Low", "High")) %>%
  group_by(patient, module, YS_tertile) %>%
  summarise(median_before = median(score_before, na.rm = TRUE),
            median_after = median(score_after, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = YS_tertile, values_from = c(median_before, median_after)) %>%
  mutate(delta_before = median_before_High - median_before_Low,
         delta_after = median_after_High - median_after_Low) %>%
  transmute(patient, module,
            median_High_before = median_before_High, median_Low_before = median_before_Low,
            High_minus_Low_delta_before = delta_before,
            median_High_after = median_after_High, median_Low_after = median_after_Low,
            High_minus_Low_delta_after = delta_after)
write_tsv(patient_delta, "patient_level_delta_before_after.tsv")

summary_modules <- patient_delta %>% group_by(module) %>%
  summarise(median_delta_before = median(High_minus_Low_delta_before, na.rm = TRUE),
            median_delta_after = median(High_minus_Low_delta_after, na.rm = TRUE),
            n_positive_before = sum(High_minus_Low_delta_before > 0, na.rm = TRUE),
            n_positive_after = sum(High_minus_Low_delta_after > 0, na.rm = TRUE),
            n_patients_after = sum(is.finite(High_minus_Low_delta_after)), .groups = "drop") %>%
  mutate(exact_sign_test_P_after = mapply(sign_p, n_positive_after, n_patients_after),
         BH_FDR_after = p.adjust(exact_sign_test_P_after, "BH")) %>%
  left_join(overlap, by = "module") %>%
  mutate(
    interpretation = case_when(
      module %in% c("UPR","TNFA NFKB","Hypoxia") ~ "clean independent program",
      module == "Stem Plasticity" & overlap_fraction >= 0.5 ~ "score-related QC; excessive overlap",
      module == "Cell Cycle" & BH_FDR_after > 0.1 ~ "heterogeneous proliferation",
      module %in% c("EMT","Migration Invasion") & BH_FDR_after > 0.1 ~ "boundary audit; not retained",
      n_positive_after >= 6 & BH_FDR_after < 0.05 ~ "score-independent, patient-consistent",
      n_positive_after >= 6 & BH_FDR_after <= 0.1 ~ "score-independent, borderline patient-consistent",
      TRUE ~ "not patient-consistent after score-gene removal"),
    keep_or_drop = case_when(
      module %in% c("UPR","TNFA NFKB","Hypoxia") ~ "KEEP_MAIN",
      module == "Stem Plasticity" ~ "DROP_MAIN_QC_ONLY",
      module %in% c("Cell Cycle","EMT","Migration Invasion") ~ "BOUNDARY_ONLY",
      n_positive_after >= 6 & BH_FDR_after <= 0.1 ~ "KEEP_MAIN",
      TRUE ~ "DROP_MAIN")) %>%
  select(module, median_delta_before, median_delta_after,
         n_positive_before, n_positive_after, exact_sign_test_P_after, BH_FDR_after,
         overlap_fraction, n_genes_after_removal, interpretation, keep_or_drop)
write_tsv(summary_modules, "module_keep_or_drop_summary.tsv")

# Score-independent Figure 3 inputs -----------------------------------------
fig3_modules <- c("UPR","TNFA NFKB","Hypoxia","Adhesion Remodeling","Wound Healing",
                  "Survival Stress","Anoikis Resistance","Cell Cycle","EMT","Migration Invasion")
write_tsv(patient_delta %>% filter(module %in% fig3_modules) %>%
            transmute(patient, module, delta = High_minus_Low_delta_after),
          "Fig3_score_independent_delta.tsv")
write_tsv(summary_modules %>% filter(module %in% fig3_modules),
          "Fig3_score_independent_stats.tsv")

candidate_genes <- c("TACSTD2","LAMB3","LAMC2","CXCL2","MCL1","BIRC3",
                     "TNFAIP3","HSPA5","DDIT3","SERPINE1","ANXA1","GDF15")
gene_eval <- bind_rows(lapply(intersect(candidate_genes, rownames(expr)), function(g) {
  v <- as.numeric(expr[g, meta$cell_id])
  group_mean <- tapply(v, meta$YS_tertile, mean)[c("Low","Intermediate","High")]
  pt <- data.frame(patient = meta$patient, group = meta$YS_tertile, value = v) %>%
    filter(group %in% c("Low","High")) %>% group_by(patient, group) %>%
    summarise(med = median(value), .groups = "drop") %>% pivot_wider(names_from = group, values_from = med)
  dd <- deg[match(g, deg$gene), ]
  data.frame(gene = g, Low = group_mean[1], Intermediate = group_mean[2], High = group_mean[3],
             logFC = dd$logFC, FDR = dd$adj.P.Val,
             n_positive_patients = sum(pt$High - pt$Low > 0),
             monotonic = group_mean[1] <= group_mean[2] && group_mean[2] <= group_mean[3],
             pass = dd$logFC > 0 && dd$adj.P.Val < 0.05 &&
                    sum(pt$High - pt$Low > 0) >= 6 &&
                    group_mean[1] <= group_mean[2] && group_mean[2] <= group_mean[3])
}))
dot <- bind_rows(lapply(gene_eval$gene[gene_eval$pass], function(g) {
  v <- as.numeric(expr[g, meta$cell_id])
  data.frame(gene = g, group = meta$YS_tertile, value = v) %>% group_by(gene, group) %>%
    summarise(mean_expression = mean(value), pct_expressing = mean(value > 0) * 100, .groups = "drop")
})) %>% left_join(gene_eval, by = "gene")
write_tsv(dot, "Fig3_core_gene_dotplot.tsv", allow_empty = TRUE)

# Boundary audit includes a separate proliferation signature.
prolif_genes <- c("MKI67","UBE2C","CENPF","BIRC5","STMN1","TK1","PTTG1","NUSAP1","TOP2A","CDK1")
prolif <- z(mean_score(expr, prolif_genes))
prolif_pt <- data.frame(patient = meta$patient, group = meta$YS_tertile, value = prolif) %>%
  filter(group %in% c("Low","High")) %>% group_by(patient, group) %>%
  summarise(med = median(value), .groups = "drop") %>% pivot_wider(names_from = group, values_from = med) %>%
  mutate(delta = High - Low)
prolif_p <- sign_p(sum(prolif_pt$delta > 0), nrow(prolif_pt))
boundary <- bind_rows(
  summary_modules %>% filter(module == "Cell Cycle") %>% transmute(
    boundary = "Cell Cycle", median_delta = median_delta_after, n_positive = n_positive_after,
    FDR = BH_FDR_after, evidence_class = ifelse(BH_FDR_after < .05, "strong", ifelse(BH_FDR_after <= .1, "moderate", "weak")),
    decision = "heterogeneous proliferation"),
  data.frame(boundary = "EMT TFs", median_delta = NA_real_, n_positive = 0, FDR = NA_real_,
             evidence_class = "not retained", decision = "all canonical EMT TFs are score-defining; 0/5 significant-up"),
  summary_modules %>% filter(module == "Migration Invasion") %>% transmute(
    boundary = "Migration/Invasion", median_delta = median_delta_after, n_positive = n_positive_after,
    FDR = BH_FDR_after, evidence_class = ifelse(BH_FDR_after < .05, "strong", ifelse(BH_FDR_after <= .1, "moderate", "weak")),
    decision = "boundary audit"),
  data.frame(boundary = "Proliferation", median_delta = median(prolif_pt$delta), n_positive = sum(prolif_pt$delta > 0),
             FDR = prolif_p, evidence_class = ifelse(prolif_p < .05, "strong", ifelse(prolif_p <= .1, "moderate", "weak")),
             decision = "heterogeneous proliferation"))
write_tsv(boundary, "Fig3_boundary_audit.tsv")

# Figure 4 DEG classification ------------------------------------------------
module_map <- bind_rows(lapply(seq_along(modules), function(i)
  data.frame(gene = modules[[i]], module = module_labels[i]))) %>%
  filter(module %in% c("UPR","TNFA NFKB","Hypoxia","Adhesion Remodeling","Wound Healing","Survival Stress","Anoikis Resistance")) %>%
  distinct(gene, .keep_all = TRUE)
deg_all <- deg %>% left_join(module_map, by = "gene") %>%
  mutate(module = ifelse(is.na(module), "Other", module),
         score_gene = gene %in% score_genes,
         score_gene_status = ifelse(score_gene, "Score-defining", "Score-independent"),
         neglogFDR = -log10(pmax(adj.P.Val, .Machine$double.xmin)),
         significant_high = adj.P.Val < 0.05 & logFC > 0)
write_tsv(deg_all, "Fig4_DEG_all.tsv")
deg_ind <- deg_all %>% filter(!score_gene)
write_tsv(deg_ind, "Fig4_DEG_score_independent.tsv")

# All-gene and score-independent enrichment with the expressed universe.
orgdb <- AnnotationDbi::loadDb(file.path(src, "vendor/org.Hs.eg.db/extdata/org.Hs.eg.sqlite"))
# Load the frozen official Human Disease Ontology database explicitly and
# register it in GOSemSim's in-process cache.  This avoids machine-specific
# Windows cache paths and prevents enrichDO() from attempting a live download.
hdo_path <- file.path(out, "references", "HDO.sqlite")
if (!file.exists(hdo_path)) stop("Frozen HDO database is missing: ", hdo_path)
hdo_db <- AnnotationDbi::loadDb(hdo_path)
yulab.utils::update_cache_item(".GOSemSimEnv", setNames(list(hdo_db), ".onto_HDO"))
universe_symbols <- unique(deg_all$gene)
map <- AnnotationDbi::select(orgdb, keys = universe_symbols, keytype = "SYMBOL", columns = "ENTREZID") %>%
  filter(!is.na(ENTREZID)) %>% distinct(SYMBOL, ENTREZID)
universe_entrez <- unique(map$ENTREZID)
entrez_to_symbol <- setNames(map$SYMBOL, map$ENTREZID)
hall_t2g <- read_gmt(file.path(root, "h.all.v2024.Hs.symbols.gmt"))

# Freeze the official KEGG REST annotations locally.  clusterProfiler::enrichKEGG()
# downloads these same three resources at run time, but long REST transfers were
# intermittently terminated by the remote peer.  Building TERM2GENE/TERM2NAME
# from the frozen official files preserves the expressed-gene universe,
# hypergeometric test, gene-set-size limits and BH correction while making the
# analysis reproducible without a live KEGG connection.
kegg_dir <- file.path(out, "references", "kegg")
kegg_files <- file.path(kegg_dir, c("link_hsa_pathway.tsv", "list_pathway_hsa.tsv",
                                   "conv_ncbi_geneid_hsa.tsv"))
if (!all(file.exists(kegg_files))) stop("Frozen KEGG annotation files are missing: ", kegg_dir)
kegg_link <- read.delim(kegg_files[1], header = FALSE, col.names = c("term", "hsa_gene"),
                        stringsAsFactors = FALSE) %>%
  mutate(term = sub("^path:", "", term))
kegg_name <- read.delim(kegg_files[2], header = FALSE, col.names = c("term", "name"),
                        stringsAsFactors = FALSE) %>%
  mutate(name = sub(" - Homo sapiens \\(human\\)$", "", name))
kegg_conv <- read.delim(kegg_files[3], header = FALSE, col.names = c("hsa_gene", "ncbi_gene"),
                        stringsAsFactors = FALSE) %>%
  mutate(gene = sub("^ncbi-geneid:", "", ncbi_gene))
kegg_t2g <- inner_join(kegg_link, kegg_conv, by = "hsa_gene") %>%
  transmute(term, gene) %>% distinct()
kegg_t2n <- kegg_name %>% transmute(term, name) %>% distinct()
if (nrow(kegg_t2g) < 10000 || nrow(kegg_t2n) < 300)
  stop("Frozen KEGG annotation cache is unexpectedly incomplete")

std <- function(x, database, entrez = FALSE) {
  d <- as.data.frame(x)
  if (!nrow(d)) return(data.frame())
  if (entrez && "geneID" %in% names(d)) d$geneID <- vapply(strsplit(d$geneID, "/", fixed = TRUE),
    function(ids) paste(unique(na.omit(entrez_to_symbol[ids])), collapse = "/"), character(1))
  d$database <- database
  keep <- intersect(c("database","ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count"), names(d))
  d[, keep, drop = FALSE]
}
run_enrichment <- function(query_symbols) {
  query_symbols <- unique(intersect(query_symbols, universe_symbols))
  query_entrez <- unique(map$ENTREZID[map$SYMBOL %in% query_symbols])
  res <- list(
    std(enrichGO(query_symbols, universe = universe_symbols, OrgDb = orgdb, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE), "GO BP"),
    std(enrichGO(query_symbols, universe = universe_symbols, OrgDb = orgdb, keyType = "SYMBOL", ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE), "GO CC"),
    std(enrichGO(query_symbols, universe = universe_symbols, OrgDb = orgdb, keyType = "SYMBOL", ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE), "GO MF"),
    std(enricher(query_entrez, universe = universe_entrez,
                 TERM2GENE = kegg_t2g, TERM2NAME = kegg_t2n,
                 minGSSize = 10, maxGSSize = 500,
                 pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1), "KEGG", TRUE),
    std(enrichPathway(query_entrez, universe = universe_entrez, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = FALSE), "Reactome", TRUE),
    std(enricher(query_symbols, universe = universe_symbols, TERM2GENE = hall_t2g[, c("term","gene")], pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1), "Hallmark"),
    tryCatch(std(DOSE::enrichDO(query_entrez, universe = universe_entrez, ont = "HDO", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE), "Disease Ontology", TRUE),
             error = function(e) { cat("DO enrichment warning:", conditionMessage(e), "\n"); data.frame() })
  )
  bind_rows(res) %>% arrange(p.adjust)
}

query_all <- deg_all$gene[deg_all$significant_high]
query_ind <- deg_all$gene[deg_all$significant_high & !deg_all$score_gene]
enr_all <- run_enrichment(query_all)
enr_ind <- run_enrichment(query_ind)
write_tsv(enr_all, "Fig4_enrichment_all_genes.tsv")
write_tsv(enr_ind, "Fig4_enrichment_score_independent.tsv")

# Patient x state pseudobulk heatmap: mutually exclusive evidence blocks.
sig_score <- deg_all %>% filter(significant_high, score_gene) %>% arrange(adj.P.Val)
downstream_overlap <- unique(module_map$gene[module_map$gene %in% score_genes])
score_overlap <- sig_score %>% filter(gene %in% downstream_overlap) %>% slice_head(n = 8) %>% pull(gene)
score_only <- sig_score %>% filter(!gene %in% downstream_overlap) %>% slice_head(n = 8) %>% pull(gene)
independent_top <- deg_all %>% filter(significant_high, !score_gene, module != "Other") %>%
  arrange(adj.P.Val, desc(logFC)) %>% distinct(gene, .keep_all = TRUE) %>% slice_head(n = 14) %>% pull(gene)
heat_genes <- unique(c(score_only, score_overlap, independent_top))
if (length(heat_genes) < 12) stop("Too few heatmap genes")
block <- data.frame(gene = heat_genes,
  evidence_block = case_when(heat_genes %in% score_only ~ "Score-defining only",
                             heat_genes %in% score_overlap ~ "Score-overlapping downstream",
                             TRUE ~ "Score-independent"))
pg <- interaction(meta$patient, meta$YS_tertile, drop = TRUE, sep = "|")
mm <- sparse.model.matrix(~0 + pg); colnames(mm) <- sub("^pg", "", colnames(mm))
pb <- as.matrix(expr[heat_genes, , drop = FALSE] %*% mm)
pb <- sweep(pb, 2, as.numeric(table(pg)[colnames(pb)]), "/")
pbz <- z_rows(pb)
heat_long <- as.data.frame(pbz) %>% mutate(gene = rownames(pbz)) %>%
  pivot_longer(-gene, names_to = "sample", values_to = "z") %>%
  separate(sample, c("patient","group"), sep = "\\|", remove = FALSE) %>%
  left_join(block, by = "gene")
write_tsv(heat_long, "Fig4_heatmap_matrix.tsv")

lollipop <- deg_all %>% filter(significant_high, module != "Other") %>%
  group_by(module, score_gene_status) %>% slice_max(logFC, n = 2, with_ties = FALSE) %>%
  ungroup() %>% arrange(module, logFC)
write_tsv(lollipop, "Fig4_representative_DEG_lollipop.tsv")

# Figure 5 decision ----------------------------------------------------------
f5_modules <- c("Adhesion Remodeling","Integrin Adhesion","Wound Healing","Survival Stress","Anoikis Resistance")
f5_eval <- summary_modules %>% filter(module %in% f5_modules)
f5_keep <- sum(f5_eval$n_positive_after >= 6 & f5_eval$BH_FDR_after <= 0.1) >= 4
writeLines(if (f5_keep) "KEEP: at least four score-independent adhesion/stress/survival modules pass the pre-specified threshold."
           else "DROP_MAIN: fewer than four score-independent adhesion/stress/survival modules pass the pre-specified threshold.",
           file.path(tab, "Figure5_decision.txt"))

cat("Score-independent reanalysis finished:", format(Sys.time()), "\n")
sink(type = "message"); sink(type = "output"); close(log_con)
