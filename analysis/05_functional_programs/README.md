# Functional programs

## Purpose

study-defined score-gene-excluded module definitions, existing module scores and external Hallmark resources.

## Related manuscript content

Fig. 4 and Supplementary Fig. S1. Figure-level links are used where the final panel selection remains to be confirmed.

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
| 0 | [build_hallmark_ranking.R](scripts/build_hallmark_ranking.R) | Figure 4E/F GSEA ranking input from the score-excluded Wald table | `Figure4_E_GSEA_v2_ranking.tsv` |
| 1 | [score_independent_reanalysis.R](scripts/score_independent_reanalysis.R) | score-gene-excluded discovery reanalysis | Source-defined files |
| 2 | [finalize_companion_tables.R](scripts/finalize_companion_tables.R) | Finalize score-gene-excluded companion tables | Source-defined files |
| 3 | [build_program_landscape_sources.py](scripts/build_program_landscape_sources.py) | Build program-landscape and lineage source tables | Source-defined files |
| 4 | [lineage_and_external_context.py](scripts/lineage_and_external_context.py) | Lineage and external-cohort context summaries | Source-defined files |
| 5 | [hallmark_gsea_and_leading_edge.R](scripts/hallmark_gsea_and_leading_edge.R) | Hallmark GSEA and leading-edge completion | Source-defined files |
| 6 | [functional_evidence_landscape.R](scripts/functional_evidence_landscape.R) | Figure 4 integrated patient functional evidence | Source-defined files |
| 7 | [leading_edge_architecture.R](scripts/leading_edge_architecture.R) | Figure 4 leading-edge architecture | Source-defined files |
| 8 | [external_hallmark_analysis.R](scripts/external_hallmark_analysis.R) | study-defined Hallmark enrichment and patient-level recurrence | Source-defined files |

### Figure 4E/F chain (previously unresolved)

The upstream step that feeds the Hallmark enrichment was missing. It is now
shipped as step 0 and the chain is closed:

```
score-excluded 17,597-gene DESeq2 Wald table        (Fig. 3 I effect-consistency source;
                                                     design ~ patient + state, High vs Other)
  -> build_hallmark_ranking.R                       deterministic projection; no model refit,
                                                     no gene filtering, no statistic recomputed
  -> Figure4_E_GSEA_v2_ranking.tsv                  17,597 ranked genes
  -> external_hallmark_analysis.R /
     hallmark_gsea_and_leading_edge.R               fgseaMultilevel over the 50 Hallmark sets
  -> Figure4_E_GSEA_v2_results.tsv                  saved enrichment result
  -> Figure 4E / 4F  and  leading_edge_architecture.R
```

All 38 YAP17/Stem21 genes were removed before gene filtering and testing, so the
testable universe is score-gene-excluded. The script asserts the
`score_definition_genes_excluded_before_testing` flag is TRUE for every retained
row rather than trusting it, and refuses to run if the source table has anything
other than 17,597 genes.

Verified against the saved result: 17,597 genes, all flagged score-gene-excluded;
`ranking_statistic` identical to `Wald_stat` (maximum absolute difference 0); 50
Hallmark pathways in the saved enrichment table. The recovered public ranking
reproduces the saved ranking exactly up to the internal order of one tied pair
(`LA16c-395F10.2`, `DEFB124`, identical statistic and FDR), which does not affect
the GSEA ranking.

```bash
Rscript analysis/05_functional_programs/scripts/build_hallmark_ranking.R \
  --source-table <Figure3_I_effect_consistency_source.tsv> --output-dir <NEW_output_directory>
```

`--dry-run` performs validation only. The final render and panel composition of
Figure 4E/4F remain a manual step and are not certified by this chain.

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
