# Data preparation

## Purpose

Raw and processed gene expression, original author cell annotations and exact cell identifiers.

## Related manuscript content

Fig. 1–2. Figure-level links are used where the final panel selection remains to be confirmed.

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
| 1 | [pseudobulk_volcano_panel.R](scripts/pseudobulk_volcano_panel.R) | Figure 2 volcano threshold panel | Source-defined files |
| 2 | [prepare_cna_inputs.R](scripts/prepare_cna_inputs.R) | CNA matrix and identity preparation | Source-defined files |

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
