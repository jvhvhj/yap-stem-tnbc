# Post-neoadjuvant residual-TNBC clinical association

This module applies Program146 to public WTS reads from Seo et al. (2025), BioProject [PRJNA1256162](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1256162). The primary subset contains 48 Neo HR-negative/HER2-negative patients, D=34 and ND=14, with RCB II/III. Each patient contributes one fresh-frozen residual-tumor sample collected at curative surgery after NAC and one paired-end WTS run.

The original study selected extreme outcomes: distant metastasis within two years of surgery for D, and no recurrence for ND. The source reports no disease progression or breast-cancer-specific death during seven years of follow-up for ND. This is not an unselected survival cohort or prospective clinical validation.

## Inputs

`metadata/Seo_TNBC48_run_manifest.tsv` supplies public patient/sample IDs and accessions. Obtain paired reads from SRA; raw reads and Salmon output directories are not distributed here. `metadata/Program146_clinical_mapping.tsv` preserves all 146 members, of which 141 have measurable exact GENCODE v36 mappings. Five absent exact symbols are reported separately. No alias rescue or replacement was used.

Processed patient scores, model coefficients and nested-model result tables are not included in this code repository. Patient IDs in the input manifest are the source study's deidentified IDs, not personal identifiers.

## Reproduction sequence

The supplied portable scripts implement the documented workflow. They are not claimed to be byte-identical original execution scripts, and complete numerical rerun equivalence has not been established during release preparation.

1. Validate the 48-run manifest with `scripts/01_prepare_run_manifest.py`. Download paired reads using the listed SRA accessions.
2. Obtain GENCODE v36 transcript FASTA, GTF and the complete matching GRCh38 genome FASTA. `scripts/02_salmon_quantification.sh index` builds the full-genome decoy-aware index with k=31 and `--gencode`; `quantify` uses Salmon 1.10.3 with `-l A --seqBias --gcBias`.
3. `scripts/03_gencode_v36_aggregation.py mapping` creates the exact transcript-to-gene table from the GTF. Its `aggregate` command sums transcript TPM and estimated read counts by exact gene ID. Supply Salmon directories named by public patient ID.
4. `scripts/04_score_program146.py` uses the supplied 141/146 mapping: log2(TPM+1), gene-wise z across 48 patients, mean of 141 gene z scores, then score-wise z across the same 48. Both standard deviations use n−1.
5. `scripts/05_merge_clinical_metadata.py` performs exact-ID merging and defines the 45-patient operative HR-negative/HER2-negative sensitivity subset. It does not re-standardize scores.
6. `scripts/06_firth_models.R` fits unadjusted, RCB-adjusted primary and operative-subtype sensitivity models in R 4.3.3 / logistf 1.26.1, with profile penalized-likelihood intervals.
7. `scripts/07_nested_model_comparison.R` compares the primary model with RCB-only using `anova.logistf(..., method="nested")`.
8. `scripts/08_export_publication_results.py` copies already-produced tables into semantic publication files without analysis.

Use `--help` for Python arguments. For example, from this module:

```bash
python scripts/01_prepare_run_manifest.py --manifest metadata/Seo_TNBC48_run_manifest.tsv --output work/read_manifest.tsv
python scripts/03_gencode_v36_aggregation.py mapping --gtf references/gencode.v36.annotation.gtf.gz --output work/tx2gene.tsv
python scripts/03_gencode_v36_aggregation.py aggregate --manifest metadata/Seo_TNBC48_run_manifest.tsv --mapping work/tx2gene.tsv --quant-dir work/salmon --output-dir work/gene_expression
python scripts/04_score_program146.py --tpm work/gene_expression/gene_TPM.tsv.gz --mapping metadata/Program146_clinical_mapping.tsv --manifest metadata/Seo_TNBC48_run_manifest.tsv --output work/scores.tsv
python scripts/05_merge_clinical_metadata.py --manifest metadata/Seo_TNBC48_run_manifest.tsv --scores work/scores.tsv --output work/clinical.tsv
Rscript scripts/06_firth_models.R work/clinical.tsv work/models
Rscript scripts/07_nested_model_comparison.R work/models/clinical_models.rds work/models/incremental.txt
```

## Interpretation

Existing results show an association with the early distant-recurrence phenotype in an independent post-neoadjuvant residual-TNBC extreme-outcome cohort after accounting for RCB class. The nested comparison is chi-square 17.95559, df=1, P=2.261198e-5. This does not establish a validated clinical biomarker, absolute recurrence-risk prediction, causation, or performance in an unselected clinical population. There was no outcome-guided gene selection, gene reweighting or cutoff optimization.

The project Salmon reprocessing is distinct from the original study's STAR/RSEM quantification. See [clinical environments and references](../../environment/clinical_analysis.md).
