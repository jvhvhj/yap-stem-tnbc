# Supplementary robustness analyses

## Purpose

Patient-level effects and source tables retained in the manuscript.

## Related manuscript content

Supplementary Fig. S1–S4. Figure-level links are used where the final panel selection remains to be confirmed.

## Dataset

GSE176078; ARTEMIS; GSE210616; BSW2. Patient identifiers and the original inclusion criteria are preserved. Cells are not substituted for biological replicates.

## Inputs

See [the input manifest](../../config/input_manifest.tsv) and the file reads in each script. Large expression objects are not distributed in this repository. The historical input/output filename contracts have been retained inside an isolated analysis workspace; source folders in this repository are organized by scientific question.

## study-defined definitions

Use the exact [YAP17, Stem21 and Program146 files](../../config/signatures/README.md). No gene substitution, new state threshold or patient exclusion is introduced by this source snapshot. Each module retains its original normalization/standardization scope rather than imposing a new global formula.

## Scripts

The numbering below indexes files, not an automatically executable pipeline. Dependencies and source gaps must be resolved before a full execution sequence can be certified.

| Order | Script | Purpose | Main output |
| --- | --- | --- | --- |
| 1 | [effect_size_and_shrinkage_sensitivity.R](scripts/effect_size_and_shrinkage_sensitivity.R) | Effect-size threshold and LFC-shrinkage sensitivity | Source-defined files |
| 2 | [render_s2a_patient_robustness.R](scripts/render_s2a_patient_robustness.R) | Render S2A patient-resolved robustness | Source-defined files |
| 3 | [render_s2b_matched_null.R](scripts/render_s2b_matched_null.R) | Render S2B matched-null calibration | Source-defined files |
| 4 | [render_s2c_state_definition.R](scripts/render_s2c_state_definition.R) | Render S2C state-definition robustness | Source-defined files |
| 5 | [render_s2d_program_stability.R](scripts/render_s2d_program_stability.R) | Render S2D Program146 stability | Source-defined files |
| 6 | [render_s3a_adjustment_displacement.R](scripts/render_s3a_adjustment_displacement.R) | Render S3A adjustment displacement | Source-defined files |
| 7 | [render_s3b_patient_distribution.R](scripts/render_s3b_patient_distribution.R) | Render S3B patient-level distribution | Source-defined files |
| 8 | [render_s3c_state_sensitivity.R](scripts/render_s3c_state_sensitivity.R) | Render S3C state sensitivity | Source-defined files |
| 9 | [render_s3d_leave_one_state_out.R](scripts/render_s3d_leave_one_state_out.R) | Render S3D leave-one-cancer-state-out robustness | Source-defined files |
| 10 | [state_omission_and_spatial_summaries.py](scripts/state_omission_and_spatial_summaries.py) | Leave-one-cancer-state-out and paired-section summaries | Source-defined files |

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
