# Malignant-cell and CNA sensitivity analysis

## Purpose

The fixed 10,836 malignant epithelial cells and same-patient reference cells; returned RNA-CNA results.

## Related manuscript content

Fig. 2 and Supplementary Fig. S1. Figure-level links are used where the final panel selection remains to be confirmed.

## Dataset

GSE176078. Patient identifiers and the original inclusion criteria are preserved. Cells are not substituted for biological replicates.

## Inputs

See [the input manifest](../../config/input_manifest.tsv) and the file reads in each script. Large expression objects are not distributed in this repository. The historical input/output filename contracts have been retained inside an isolated analysis workspace; source folders in this repository are organized by scientific question.

## Frozen definitions

Use the exact [YAP17, Stem21 and Program146 files](../../config/signatures/README.md). No gene substitution, new state threshold or patient exclusion is introduced by this source snapshot. Each module retains its original normalization/standardization scope rather than imposing a new global formula.

## Scripts

The numbering below indexes files, not an automatically executable pipeline. Dependencies and source gaps must be resolved before a full execution sequence can be certified.

| Order | Script | Purpose | Main output |
| --- | --- | --- | --- |
| 0 | [extract_wu_malignant_population.R](scripts/extract_wu_malignant_population.R) | Derive the fixed malignant analytical population from the source atlas annotation | `cell_ids.tsv`, `per_patient_counts.tsv`, `extraction_summary.tsv` |
| 1 | [run_copykat_patient.R](scripts/copykat/run_copykat_patient.R) | Reference-subsampled CopyKAT inference | Source-defined files |
| 2 | [run_scevan_patient.R](scripts/scevan/run_scevan_patient.R) | Anchored patient-wise SCEVAN inference | Source-defined files |
| 3 | [run_infercnv_patient.R](scripts/infercnv/run_infercnv_patient.R) | Continuous inferCNV burden and native HMM status | Source-defined files |
| 4 | [integrate_caller_results.py](scripts/integrate_caller_results.py) | Exact cell-ID integration and patient-unit sensitivity | Source-defined files |
| 5 | [recover_program_scores.py](scripts/recover_program_scores.py) | Existing Program146 score recovery and sensitivity | Source-defined files |
| 6 | [module_caller_sensitivity.py](scripts/module_caller_sensitivity.py) | Six fixed modules within caller-supported subsets | Source-defined files |
| 7 | [recover_mcl1_expression.R](scripts/recover_mcl1_expression.R) | Read existing normalized expression by exact cell ID | Source-defined files |

### Derivation entry (step 0)

`extract_wu_malignant_population.R` is the public entry to the analytical
population and precedes every caller in this module and every score in
`analysis/03_yap_stem_continuum`. It is **population extraction, not a
classification**: the malignant identity is the source object's own annotation,
and the script applies exactly two filters.

```
full Wu atlas        100,064 cells / 26 patients
  subtype == "TNBC"   42,512 cells / 10 donors
  celltype_major == "Cancer Epithelial"
                      10,836 cells /  8 donors   <- the frozen analytical population
```

Donors whose TNBC subset contains no Cancer Epithelial cells drop out as a
consequence of the annotation, not by a separate exclusion rule. CopyKAT, SCEVAN
and inferCNV results are **not** used to define inclusion; they are supporting
sensitivity evidence computed downstream.

```bash
Rscript analysis/02_malignant_cell_and_cnv/scripts/extract_wu_malignant_population.R \
  --atlas-metadata <Wu2021_full_atlas_umap_metadata.csv> \
  --output-dir <NEW_output_directory>
```

`--dry-run` validates the input table only. Add `--atlas-object <Wu2021_full_atlas.rds>`
to also write the subset Seurat object. The extraction has been verified to
reproduce the frozen 10,836-cell `cell_id` set exactly.

## Reproduction

The current snapshot supports source inspection and static tests. It does **not** yet support an end-to-end reproduction claim. Do not point these scripts at the original analysis archive. Review [the outstanding source and input requirements](../../docs/reproducibility.md) first.

```bash
python tests/smoke/check_repository.py --rscript Rscript
```

Each recovered script retains its original arguments and scientific parameters. Its argument parser is the execution contract; a new unified biological runner has deliberately not been substituted while source gaps remain.

## Key parameters

Random seeds, eligibility conditions, scoring expressions, model formulas and native missing-value handling remain in the recovered scripts. The script provenance index records path relocation only. Dataset-specific technical adjustments must not be replaced by a common newly invented covariate model.

## Outputs

Result tables and rendering outputs are named in the individual scripts. No result has been regenerated during repository preparation. Existing results are used only for source/membership/denominator checks.

## Software

See [environment documentation](../../environment/README.md), recorded package versions and script dependency lists. Historical environments are distinguished from the environment used for static checks.

## Notes

This is a recoverable source snapshot, not a certified full rerun. Missing original sources and unresolved publication selection are listed explicitly. Native diagnostic/status fields remain in code because changing them could alter downstream joins; internal project reports are not included as publication documents.
