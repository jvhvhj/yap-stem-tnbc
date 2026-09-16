# Analyses and manuscript content

Each entry names the script a reader should open first and the manuscript content
it supports. Scripts live in the module directories; each module README gives the
inputs, the frozen definitions and the current reproduction boundary.

Machine-readable detail — scientific question, dataset, authoritative script and
status per unit of analysis — is in [`docs/figure_code_map.tsv`](../docs/figure_code_map.tsv).

| Analysis | Description | Related figure(s) | Code entry |
| --- | --- | --- | --- |
| Malignant epithelial population | Fixed 10,836 malignant epithelial cells from 8 donors, selected by the source atlas's own TNBC and Cancer Epithelial annotation | Fig. 2 | [`02_malignant_cell_and_cnv/scripts/extract_wu_malignant_population.R`](02_malignant_cell_and_cnv/scripts/extract_wu_malignant_population.R) |
| CopyKAT | Per-patient RNA-CNA classification (sensitivity evidence; does not define inclusion) | Fig. 2; Suppl. Fig. S1 | [`02_malignant_cell_and_cnv/scripts/copykat/run_copykat_patient.R`](02_malignant_cell_and_cnv/scripts/copykat/run_copykat_patient.R) |
| SCEVAN | Per-patient RNA-CNA classification, anchored mode | Fig. 2; Suppl. Fig. S1 | [`02_malignant_cell_and_cnv/scripts/scevan/run_scevan_patient.R`](02_malignant_cell_and_cnv/scripts/scevan/run_scevan_patient.R) |
| inferCNV | Descriptive CNA burden per patient and per cell | Fig. 2; Suppl. Fig. S1 | [`02_malignant_cell_and_cnv/scripts/infercnv/run_infercnv_patient.R`](02_malignant_cell_and_cnv/scripts/infercnv/run_infercnv_patient.R) |
| Multi-caller integration | Caller agreement, patient-unit sensitivity, CNA-burden restriction subsets | Fig. 2; Suppl. Fig. S1 | [`02_malignant_cell_and_cnv/scripts/integrate_caller_results.py`](02_malignant_cell_and_cnv/scripts/integrate_caller_results.py) |
| YAP17 / Stem21 scoring | Unweighted mean of the measurable signature genes in the log-normalised RNA data layer | Fig. 3 | [`03_yap_stem_continuum/scripts/score_yap17_stem21.R`](03_yap_stem_continuum/scripts/score_yap17_stem21.R) |
| Joint axis and state assignment | Joint = sum of the two component z scores; equal thirds define Low / Intermediate / High | Fig. 3 | [`03_yap_stem_continuum/scripts/derive_joint_axis_and_tertiles.R`](03_yap_stem_continuum/scripts/derive_joint_axis_and_tertiles.R) |
| State-assignment validation | Reproduce the frozen Joint column and state labels; report every equal-thirds rule tested | Fig. 3 | [`03_yap_stem_continuum/scripts/validate_joint_state_reproduction.R`](03_yap_stem_continuum/scripts/validate_joint_state_reproduction.R) |
| Patient-resolved association | Within-patient Spearman correlation and bootstrap intervals | Fig. 3; Suppl. Fig. S2 | [`03_yap_stem_continuum/scripts/wu_robustness_models.R`](03_yap_stem_continuum/scripts/wu_robustness_models.R) |
| State-definition concordance | Agreement between pooled-, within-patient- and patient-standardised-tertile definitions | Fig. 3 | [`03_yap_stem_continuum/scripts/state_definition_concordance_panel.R`](03_yap_stem_continuum/scripts/state_definition_concordance_panel.R) |
| Program146 derivation | Patient-blocked High-versus-Other pseudobulk DE; unshrunk log2FC ≥ 0.5 and ≥ 6/8 patient direction over 17,597 score-excluded genes | Fig. 4A/B | [`04_program146/scripts/01_derive_program146.R`](04_program146/scripts/01_derive_program146.R) |
| Pseudobulk verification | Independent patient-level verification of the pseudobulk effects | Fig. 4A/B | [`04_program146/scripts/verify_patient_pseudobulk_effects.R`](04_program146/scripts/verify_patient_pseudobulk_effects.R) |
| Functional-program landscape | Six programs (four Hallmark plus two custom) across within-patient deciles | Fig. 4C/D | [`05_functional_programs/scripts/functional_evidence_landscape.R`](05_functional_programs/scripts/functional_evidence_landscape.R) |
| Hallmark ranking input | Score-excluded 17,597-gene DESeq2 Wald ranking used as the GSEA ranking | Fig. 4E/F | [`05_functional_programs/scripts/build_hallmark_ranking.R`](05_functional_programs/scripts/build_hallmark_ranking.R) |
| Hallmark enrichment | All 50 Hallmark sets over the score-excluded Wald ranking | Fig. 4E/F | [`05_functional_programs/scripts/external_hallmark_analysis.R`](05_functional_programs/scripts/external_hallmark_analysis.R) |
| Leading-edge architecture | Gene–pathway membership for the displayed Hallmark pathways | Fig. 4F | [`05_functional_programs/scripts/leading_edge_architecture.R`](05_functional_programs/scripts/leading_edge_architecture.R) |
| GSE180286 orthogonal transfer | Fixed-signature transfer and ranked-signature enrichment | Fig. 4G/H | [`06_independent_tnbc_validation/scripts/gse180286_closure.R`](06_independent_tnbc_validation/scripts/gse180286_closure.R) |
| ARTEMIS validation | 78-patient association, 77-patient program transfer, state occupancy and within-state coupling | Fig. 5 | [`06_independent_tnbc_validation/scripts/artemis_validation.py`](06_independent_tnbc_validation/scripts/artemis_validation.py) |
| Cancer-state robustness | State occupancy versus within-state coupling; leave-one-state-out influence | Fig. 5; Suppl. Fig. S4 | [`supplementary/scripts/render_s3c_state_sensitivity.R`](supplementary/scripts/render_s3c_state_sensitivity.R) |
| GSE210616 spatial recurrence | 22 patients / 43 sections; spot scores, purity- and depth-adjusted patient recurrence, neighbourhood co-organization | Fig. 6A–F | [`07_spatial_validation/scripts/compute_spatial_scores_worker.R`](07_spatial_validation/scripts/compute_spatial_scores_worker.R) — see [spatial reproduction](../docs/spatial_reproduction.md) |
| BSW2 independent replication | Nine-patient Visium replication under a separate normalisation implementation | Fig. 6F/G | [`07_spatial_validation/scripts/bsw2_spatial_replication.R`](07_spatial_validation/scripts/bsw2_spatial_replication.R) — see [spatial reproduction](../docs/spatial_reproduction.md) |
| Matched-null calibration | Observed versus gene-number/expression/detection-matched random signatures | Suppl. Fig. S2B | [`supplementary/scripts/render_s2b_matched_null.R`](supplementary/scripts/render_s2b_matched_null.R) — renders from a frozen null table; the generator is not included |
| Scoring-method sensitivity | Five predefined scoring approaches across cohorts | Suppl. Fig. S2A | [`supplementary/scripts/render_s2a_patient_robustness.R`](supplementary/scripts/render_s2a_patient_robustness.R) |
| Effect-size and shrinkage sensitivity | Threshold and apeglm-shrinkage sensitivity of Program146 membership | Suppl. Fig. S2D | [`supplementary/scripts/effect_size_and_shrinkage_sensitivity.R`](supplementary/scripts/effect_size_and_shrinkage_sensitivity.R) |
| Hippo/YAP–TNBC schematic | Conceptual overview of Hippo/YAP–TAZ signalling and its TNBC context | Fig. 1 | Not applicable — conceptual figure, no analysis code |

## Derivation chain

The population-to-state chain is public and verified end to end:

```
source atlas annotation
  -> extract_wu_malignant_population.R      10,836 cells / 8 donors
  -> score_yap17_stem21.R                   YAP17 17/17, Stem21 21/21 measurable
  -> derive_joint_axis_and_tertiles.R       Joint = z(YAP) + z(Stem); High 3,612 / Other 7,224
  -> validate_joint_state_reproduction.R    max |delta Joint| 2.0e-14; 0 state mismatches
  -> downstream analyses
```

## Boundaries

- **Execution linkage for the RNA-CNA callers is unresolved.** The CopyKAT, SCEVAN
  and inferCNV scripts are recovered outbound copies; no execution-time copy or
  checksum was returned from the workstation that ran them. They are supporting
  sensitivity evidence and are not part of the population definition.
- **Some rows are content matches, not certified render sources.** Where the
  approved panel version could not be tied to one file, the row status says so.
  A filename or a timestamp is not treated as proof of final placement.
- **Supporting bulk and cross-cancer cohorts are out of code-release scope** —
  see [`data/README.md`](../data/README.md).

See [the reproduction guide](../docs/reproducibility.md) before running anything.

### Datasets

See [`data/README.md`](../data/README.md) and [`config/datasets.tsv`](../config/datasets.tsv).

### Fixed definitions

[`config/signatures/`](../config/signatures/README.md) holds YAP17, Stem21 and
Program146. These are definitions, not claims that every external dataset
measures every member.
