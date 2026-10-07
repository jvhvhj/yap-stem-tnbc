# YAP–Stem transcriptional continuum

## Purpose

YAP17 and Stem21 gene definitions, RNA/data expression and patient-resolved scores.

## Related manuscript content

Fig. 3 and Supplementary Fig. S2. Figure-level links are used where the final panel selection remains to be confirmed.

## Dataset

GSE176078. Patient identifiers and the original inclusion criteria are preserved. Cells are not substituted for biological replicates.

## Inputs

See [the input manifest](../../config/input_manifest.tsv) and the file reads in each script. Large expression objects are not distributed in this repository. The historical input/output filename contracts have been retained inside an isolated analysis workspace; source folders in this repository are organized by scientific question.

## study-defined definitions

Use the exact [YAP17, Stem21 and Program146 files](../../config/signatures/README.md). No gene substitution, new state threshold or patient exclusion is introduced by this source snapshot. Each module retains its original normalization/standardization scope rather than imposing a new global formula.

## Scripts

The numbering below indexes files, not an automatically executable pipeline. Dependencies and source gaps must be resolved before a full execution sequence can be certified.

| Order | Script | Purpose | Main output |
| --- | --- | --- | --- |
| 0a | [score_yap17_stem21.R](scripts/score_yap17_stem21.R) | YAP17 and Stem21 component scores | `component_scores.tsv`, `signature_coverage.tsv` |
| 0b | [derive_joint_axis_and_tertiles.R](scripts/derive_joint_axis_and_tertiles.R) | Joint axis and Low/Intermediate/High state assignment | `joint_state.tsv`, `state_counts.tsv` |
| 0c | [validate_joint_state_reproduction.R](scripts/validate_joint_state_reproduction.R) | Validate the derivation against the study-defined columns | `state_reproduction_validation.tsv`, `tertile_rule_equivalence.tsv` |
| 1 | [state_definition_concordance_panel.R](scripts/state_definition_concordance_panel.R) | Figure 3E concordance label correction | Source-defined files |
| 2 | [wu_robustness_models.R](scripts/wu_robustness_models.R) | Patient-level scoring and covariate robustness | Source-defined files |

### Derivation entry (steps 0a–0c)

These three scripts are the public chain that turns the study-defined malignant
population into the state variable every downstream analysis reads. Run them from
the repository root, into a new output directory each time.

```bash
Rscript analysis/03_yap_stem_continuum/scripts/score_yap17_stem21.R \
  --counts <Wu_malignant_Seurat_RDS> --output-dir <NEW_dir>

Rscript analysis/03_yap_stem_continuum/scripts/derive_joint_axis_and_tertiles.R \
  --component-scores <NEW_dir>/component_scores.tsv --output-dir <NEW_dir2>

Rscript analysis/03_yap_stem_continuum/scripts/validate_joint_state_reproduction.R \
  --derived <NEW_dir2>/joint_state.tsv \
  --frozen  <frozen_cell_state_metadata.csv> --output-dir <NEW_dir3>
```

Rules fixed before any result was inspected:

- Component score = unweighted arithmetic mean of the **measurable** members of
  the signature in the log-normalised RNA `data` layer. Absent genes stay absent;
  no alias substitution and no gene re-selection.
- `Joint = z(YAP17 mean) + z(Stem21 mean)`, standardised across the study-defined
  10,836-cell population with the sample SD (cohort-level, not per-patient).
- Low / Intermediate / High = equal thirds of the Joint value.
- Validated: 17/17 YAP17 and 21/21 Stem21 genes measurable; High 3,612 /
  Other 7,224; maximum absolute Joint difference from the study-defined column
  `2.0e-14`; **0 state mismatches** across all 10,836 cells.

On the tertile rule: n = 10,836 is exactly divisible by 3 and no tie spans a
boundary, so five standard equal-thirds rules each reproduce the study-defined state
with zero mismatches. That is a property of this population, not evidence about
which code path the original authors ran; `validate_joint_state_reproduction.R`
reports every rule it tested rather than hiding the residual ambiguity.

The study-defined validation outputs are not distributed in this code repository. The Joint axis is
score-defining and must not be used as an outcome of itself.

## Reproduction

The current snapshot supports source inspection and static tests. It does **not** yet support an end-to-end reproduction claim. Do not point these scripts at the original analysis archive. Review [the outstanding source and input requirements](../../docs/reproducibility.md) first.

```bash
python tests/smoke/check_repository.py --rscript Rscript
```

Each recovered script retains its original arguments and scientific parameters. Its argument parser is the execution contract; a new unified biological runner has deliberately not been substituted while source gaps remain.

## Key parameters

The saved Wu forest is an extended-adjustment endpoint: its covariates include Hypoxia and UPR in addition to technical variables. It must not be labelled technical-only or substituted for the separate matched-null technical-adjustment endpoint. The source-object cohort retains CID3963 with its documented receptor-metadata discrepancy; this does not establish independently reconfirmed clinical receptor homogeneity.

Random seeds, eligibility conditions, scoring expressions, model formulas and native missing-value handling remain in the recovered scripts. The script provenance index records path relocation only. Dataset-specific technical adjustments must not be replaced by a common newly invented covariate model.

## Outputs

Result tables and rendering outputs are named in the individual scripts. No result has been regenerated during repository preparation. Existing results are used only for source/membership/denominator checks.

## Software

See [environment documentation](../../environment/README.md), recorded package versions and script dependency lists. Historical environments are distinguished from the environment used for static checks.

## Notes

This is a recoverable source snapshot, not a certified full rerun. Missing original sources and unresolved publication selection are listed explicitly. Native diagnostic/status fields remain in code because changing them could alter downstream joins; internal project reports are not included as publication documents.
