# Independent TNBC validation

## ARTEMIS six-program and Program146 associations

`scripts/artemis_functional_program_associations.R` is the recovered numerical
source for the six-program continuous associations and Program146 transfer.
It preserves the original scientific calculations. Only the input root and
output directory are configurable:

```bash
AHIPPO_ANALYSIS_ROOT=/path/to/analysis_workspace \
AHIPPO_RESULTS_DIR=/path/to/new_results \
Rscript analysis/06_independent_tnbc_validation/scripts/artemis_functional_program_associations.R
```

This script requires the existing workspace input paths declared in its `paths`
vector, including the cancer-cell score cache, h5ad and discovery module/gene
tables. It does not run from the small public repository alone. Do not execute it
against the original frozen output directory.

For six-program associations, each patient's Joint axis and each program score
are separately OLS-residualized against log1p(nCount_RNA), log1p(nFeature_RNA),
S and G2M scores. The association is Spearman correlation of ordinary residuals.
There is no additional z-standardization or reconstruction of Joint from
separately adjusted YAP and Stem scores. This differs from the primary
YAP17–Stem21 model, which uses untransformed nFeature_RNA. S3 lists these models
separately. The broad six-program memberships are not the compact spatial
module memberships.

## Purpose

Author-annotated malignant cells, patient identity and the discovery signatures.

## Related manuscript content

Fig. 5; orthogonal transfer in Fig. 4. Figure-level links are used where the final panel selection remains to be confirmed.

## Dataset

ARTEMIS; GSE180286. Patient identifiers and the original inclusion criteria are preserved. Cells are not substituted for biological replicates.

## Inputs

See [the input manifest](../../config/input_manifest.tsv) and the file reads in each script. Large expression objects are not distributed in this repository. The historical input/output filename contracts have been retained inside an isolated analysis workspace; source folders in this repository are organized by scientific question.

## study-defined definitions

Use the exact [YAP17, Stem21 and Program146 files](../../config/signatures/README.md). No gene substitution, new state threshold or patient exclusion is introduced by this source snapshot. Each module retains its original normalization/standardization scope rather than imposing a new global formula.

## Scripts

The numbering below indexes files, not an automatically executable pipeline. Dependencies and source gaps must be resolved before a full execution sequence can be certified.

| Order | Script | Purpose | Main output |
| --- | --- | --- | --- |
| 1 | [export_cell_cycle_genes.R](scripts/export_cell_cycle_genes.R) | Export cell-cycle reference genes | Source-defined files |
| 2 | [prepare_artemis_and_signatures.py](scripts/prepare_artemis_and_signatures.py) | Prepare ARTEMIS data and study-defined signatures | Source-defined files |
| 3 | [artemis_validation.py](scripts/artemis_validation.py) | ARTEMIS validation analyses | Source-defined files |
| 4 | [gse180286_closure.R](scripts/gse180286_closure.R) | GSE180286 completion analysis | Source-defined files |
| 5 | [prepare_gse180286_transfer_sources.R](scripts/prepare_gse180286_transfer_sources.R) | Prepare GSE180286 transfer source tables | Source-defined files |
| 6 | [orthogonal_program_transfer.R](scripts/orthogonal_program_transfer.R) | Figure 4 orthogonal program transfer | Source-defined files |
| 7 | [figure5_data_closure.R](scripts/figure5_data_closure.R) | Figure 5 final data completion | Source-defined files |
| 8 | [render_figure5.R](scripts/render_figure5.R) | Figure 5 final render | Source-defined files |
| 9 | [artemis_robustness_models.py](scripts/artemis_robustness_models.py) | Patient-level external scoring robustness | Source-defined files |

## Reproduction

The current snapshot supports source inspection and static tests. It does **not** yet support an end-to-end reproduction claim. Do not point these scripts at the original analysis archive. Review [the outstanding source and input requirements](../../docs/reproducibility.md) first.

```bash
python tests/smoke/check_repository.py --rscript Rscript
```

Each recovered script retains its original arguments and scientific parameters. Its argument parser is the execution contract; a new unified biological runner has deliberately not been substituted while source gaps remain.

## Key parameters

Random seeds, eligibility conditions, scoring expressions, model formulas and native missing-value handling remain in the recovered scripts. The script provenance index records portability and non-scientific documentation changes. Dataset-specific technical adjustments must not be replaced by a common newly invented covariate model.

## Outputs

Result tables and rendering outputs are named in the individual scripts. No result has been regenerated during repository preparation. Existing results are used only for source/membership/denominator checks.

## Software

See [environment documentation](../../environment/README.md), recorded package versions and script dependency lists. Historical environments are distinguished from the environment used for static checks.

## Notes

This is a recoverable source snapshot, not a certified full rerun. Missing original sources and unresolved publication selection are listed explicitly. Native diagnostic/status fields remain in code because changing them could alter downstream joins; internal project reports are not included as publication documents.
