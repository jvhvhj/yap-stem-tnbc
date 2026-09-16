# Program146 derivation

## Purpose

Program146 is a score-independent transcriptional program associated with the YAP–Stem High state in Wu TNBC malignant epithelial cells. Original-source provenance is resolved. The portable implementation is supplied with static checks; its full biological execution has not been performed.

## Discovery cohort

Wu 2021 / GSE176078: 10,836 malignant epithelial cells from eight TNBC patients: CID3963, CID4465, CID4495, CID44971, CID44991, CID4513, CID4515 and CID4523. The script matches exact cell IDs and preserves patient order encountered in the original metadata; it does not impose a newly sorted order.

## High versus Other definition

High is the existing `YS_tertile == "High"`. Other comprises all non-High cells. No state assignment, signature scoring, tertile calculation or malignant-cell inference is performed.

## Score-gene exclusion

All 38 distinct YAP17 and Stem21 defining genes are removed from the pseudobulk testing matrix **before expression filtering and formal testing**. The exact signature symbols and order are retained. The formal tested universe and Program146 contain none of these score genes.

## Patient-blocked pseudobulk

Integer raw RNA counts are summed by patient × state, producing 16 nonempty samples. Other is the reference level, followed by High. Patient is the biological unit; cells are not independent biological replicates.

## Differential expression

The original filter is `edgeR::filterByExpr(DGEList(pb_decoupled), design=design, min.count=5)`, with `design = ~ patient + state`. DESeq2 uses the same design, `contrast=c("state", "High", "Other")`, `alpha=0.05` and `independentFiltering=FALSE`.

The original `lfcShrink(type="normal")` display branch is retained, including its documented fallback if unavailable. **Primary selection uses unshrunk log2FC**, never the display-shrinkage column. No apeglm substitution is made.

## Program146 selection

All requirements must hold:

- BH FDR < 0.05;
- unshrunk High-versus-Other log2FC ≥ 0.5;
- positive direction in at least six of eight patients;
- membership in the score-excluded, expression-filtered testing universe.

Patient direction is the sign of `log2(DESeq2-normalized pseudobulk count + 1)` in High minus Other. Positive means `>0`, negative means `<0`, and zero is separate. Original `na.rm=TRUE` summaries are preserved; NA does not count as positive and the denominator is not silently reduced from eight.

## Result

Surviving historical outputs contain 17,597 tested genes, 362 BH FDR < 0.05 genes and 146 selected genes. Exact membership matches [Program146.tsv](../../config/signatures/Program146.tsv). These existing outputs were checked during packaging, not regenerated. The reference list is used only for a post-derivation check, never to select or force genes into the result.

## Scripts

| Script | Role |
| --- | --- |
| [July 27 original](original_archive/01_wu_pseudobulk_and_reconciliation_original_2026-07-27.R) | Corrected original derivation, exact-byte provenance copy; do not execute directly |
| [01_derive_program146.R](scripts/01_derive_program146.R) | Portable original calculations; parse/static checks and input-only preflight; full run not performed |
| [August 8 regeneration](original_archive/Figure3_I_effect_consistency_regeneration_2026-08-08.R) | Later regeneration/verification, not original discovery; provenance only |
| [verify_patient_pseudobulk_effects.R](scripts/verify_patient_pseudobulk_effects.R) | Existing downstream verification; historical input-path wiring requires review |
| [patient_gene_effect_consistency.R](scripts/patient_gene_effect_consistency.R) | Existing gene-effect figure source; final display/assembly is separate |
| [program_differential_expression_landscape.R](scripts/program_differential_expression_landscape.R) | Volcano-source family; exact current Fig.4B export remains a content-match candidate |

Current Fig.4A uses this derivation as the scientific basis of its workflow illustration; Fig.4B displays DE results. Current Fig.4E/F are Hallmark profiling/leading-edge architecture, not the historical B/C labels in some source filenames. [Figure-code mapping](../../docs/figure_code_map.tsv) distinguishes scientific source links from unresolved final render/composition provenance.

## Representative genes

ATF3, PPP1R15A, MCL1 and ZFP36L1 are verified members, listed in [Program146_representative_genes.tsv](../../config/signatures/Program146_representative_genes.tsv). CDKN1A, PERP and S100A10 are **not** members. Their use in a separate pathway or DE context must not imply Program146 membership. No new representatives were selected.

## Reproduction

Run from the repository root. Choose a NEW output directory outside the repository whose parent already exists. Existing directories, even empty, are rejected. There is no overwrite or resume option.

```bash
Rscript analysis/04_program146/scripts/01_derive_program146.R \
  --counts /path/to/inputs/Wu2021_TNBC_malignant.rds \
  --cell-metadata /path/to/inputs/Wu2021_continuous_state_cell_metadata.csv \
  --historical-score-independent /path/to/inputs/Fig4_DEG_score_independent.tsv \
  --yap-signature config/signatures/YAP17.tsv \
  --stem-signature config/signatures/Stem21.tsv \
  --output-dir /path/to/reproduction/program146_new_run \
  --dry-run
```

`--counts` is the original **Seurat RDS with the raw RNA/counts layer**, not a normalized matrix or low-dimensional embedding. No alternative loader or normalization is introduced. Metadata is CSV containing `cell_id`, `patient`, `YS_tertile`, `YAP_score`, `Stemness_score` and `YAP_Stem_Score`. IDs and definitions are validated without filtering or reordering.

The explicit additional input is the original High–Low table `0713_score_independent_rebuild/tables/Fig4_DEG_score_independent.tsv`, containing `gene` and `adj.P.Val`. It is retained because the original reconciliation builds its gene union and historical-142 comparison from this table. This historical set does not determine the High–Other criteria. Retaining it avoids rewriting the original reconciliation.

`--dry-run` checks small input tables, fixed definitions, file existence and output safety. It does not load the Seurat object, import analysis packages, create output files, aggregate counts or fit models. `--help` needs no inputs. After environment/input review and authorization, removing `--dry-run` executes the derivation in the isolated directory. This full run has not been validated during packaging.

The scientific runtime requires Seurat, Matrix, edgeR and DESeq2. [Historical environment reconstruction](../../environment/README.md) remains incomplete; no current package installation is presented as an exact historical environment.

On Windows, use a valid native UTF-8 locale for non-ASCII paths. A shell-inherited `LC_ALL=C.UTF-8` is not recognized by this Windows R runtime and can break normalized path lookup. Packaging preflight used the native Windows UTF-8 locale after clearing that unsupported override in the test subprocess only; the script does not change locale itself. Capture the actual locale with `sessionInfo()` when full historical-environment reproduction is later undertaken.

## Outputs

| File | Content |
| --- | --- |
| `pseudobulk_counts.tsv.gz` | Gene identifiers plus 16 raw-count columns |
| `pseudobulk_samples.tsv` | 16 patient/state rows, cell counts and library summaries |
| `pseudobulk_design.tsv` | Original design matrix |
| `score_gene_exclusion.tsv` | One audit row per score gene, 38 rows |
| `score_gene_exclusion_summary.tsv` | Original five-stage exclusion summary |
| `Program146_DE_all.tsv` | Full DE table, including primary unshrunk and original display-shrunken LFC |
| `Program146_DE_FDR05.tsv` | All FDR < 0.05 genes |
| `patient_log2_normalized_difference.tsv.gz` | Original gene × eight-patient differences |
| `patient_direction_summary.tsv` | Positive/negative/zero counts and original descriptive summaries |
| `gene_selection_audit.tsv` | Original reconciliation and selection criteria |
| `threshold_cascade.tsv` | Original count cascade |
| `Program146_genes.txt` | Derived program; not copied from the reference |
| `Program146_genes_historical_alias.txt` | Original second membership write, an explicit identical alias |
| `historical142_comparison_genes.txt` | Historical comparison only |
| `recurrent_downregulated_genes.txt` | Original downregulated counterpart, not Program146 |
| `representative_gene_results.tsv` | Original MCL1/PPP1R15A/ATF3 fields |
| `input_manifest.tsv`, `sessionInfo.txt`, `shrinkage_note.txt` | Input/runtime records |
| `derivation_reference_checks.tsv` | Post-derivation count/membership comparison |

The cascade's original “After expression filter” row separately evaluates the raw gene space for description. Formal testing still excludes the 38 genes before filtering. Its labels must not be interpreted as a different computational order. A discrepancy in derived membership/counts stops the script after recording checks; there is no tuning or corrective reselection.

## Notes

The portable copy retains the primary statistical expressions and replaces machine-specific paths with explicit arguments. Optional module summaries, historical comparison plots, diagnostic/final plots, enrichment and LOPO reporting are omitted because they do not feed Program146 selection. The historical reconciliation input, display-shrinkage branch, factors, seed, missing-value handling and selection expressions are retained. No figures are generated.

The original wrapper's nonzero historical return code and incomplete native intermediate archive are documented in the [archive README](original_archive/README.md). Source provenance is resolved; clean full-run equivalence is not claimed. The repository as a whole remains not ready for release.
