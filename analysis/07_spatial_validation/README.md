# Spatial recurrence and replication

## Purpose

Gene-level spatial count matrices, coordinates, section/patient registries and fixed score definitions.

## Related manuscript content

Fig. 6 and Supplementary Fig. S4. Figure-level links are used where the final panel selection remains to be confirmed.

## Dataset

GSE210616; BSW2. Patient identifiers and the original inclusion criteria are preserved. Cells are not substituted for biological replicates.

## Inputs

See [the input manifest](../../config/input_manifest.tsv) and the file reads in each script. Large expression objects are not distributed in this repository. The historical input/output filename contracts have been retained inside an isolated analysis workspace; source folders in this repository are organized by scientific question.

## Frozen definitions

Use the exact [YAP17, Stem21 and Program146 files](../../config/signatures/README.md). No gene substitution, new state threshold or patient exclusion is introduced by this source snapshot. Each module retains its original normalization/standardization scope rather than imposing a new global formula.

## Scripts

The numbering below indexes files, not an automatically executable pipeline. Dependencies and source gaps must be resolved before a full execution sequence can be certified.

| Order | Script | Purpose | Main output |
| --- | --- | --- | --- |
| 1 | [compute_spatial_scores_worker.R](scripts/compute_spatial_scores_worker.R) | Compute section-level YAP-Stem and Program146 spot scores | Source-defined files |
| 2 | [assemble_spatial_coupling_results.R](scripts/assemble_spatial_coupling_results.R) | Assemble patient-recurrent spatial coupling results | Source-defined files |
| 3 | [compute_spatial_enrichment_and_neighborhood_statistics.R](scripts/compute_spatial_enrichment_and_neighborhood_statistics.R) | Compute cohort recurrence and Lee L neighborhood statistics | Source-defined files |
| 4 | [score_spatial_modules.R](scripts/score_spatial_modules.R) | Score predefined spatial functional modules | Source-defined files |
| 5 | [assemble_spatial_module_results.R](scripts/assemble_spatial_module_results.R) | Assemble patient-level spatial module results | Source-defined files |
| 6 | [bsw2_spatial_replication.R](scripts/bsw2_spatial_replication.R) | Independent BSW2 spatial replication | Source-defined files |
| 7 | [render_spatial_robustness_panels.py](scripts/render_spatial_robustness_panels.py) | Figure 6B/C visual render | Source-defined files |
| 8 | [render_spatial_enrichment_panel.R](scripts/render_spatial_enrichment_panel.R) | Figure 6D spatial enrichment render | Source-defined files |
| 9 | [render_neighborhood_coorganization_panel.R](scripts/render_neighborhood_coorganization_panel.R) | Figure 6E neighborhood co-organization render | Source-defined files |
| 10 | [render_patient_recurrent_functional_coupling.R](scripts/render_patient_recurrent_functional_coupling.R) | Figure 6F patient-recurrent functional coupling render | Source-defined files |
| 11 | [render_independent_spatial_replication.R](scripts/render_independent_spatial_replication.R) | Figure 6G independent spatial replication render | Source-defined files |
| 12 | [render_spatial_score_maps.R](scripts/render_spatial_score_maps.R) | Figure 6A left spatial map display reset | Source-defined files |

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
