# Clinical Source Data

| File | Records | Content |
|---|---:|---|
| SourceData_Seo_clinical_patient_scores.tsv | 48 | Public patient ID, group, subtypes, RCB class, raw and standardized Program146 scores, operative sensitivity inclusion |
| SourceData_Seo_clinical_model_results.tsv | 7 | All available unadjusted, primary RCB-adjusted and operative-subtype sensitivity coefficients |
| SourceData_Seo_program146_mapping.tsv | 146 | Exact GENCODE v36 mapping, 141 measurable genes and five absent symbols |
| SourceData_Seo_incremental_model.tsv | 1 | Nested penalized likelihood-ratio comparison of RCB-only versus RCB plus Program146 |

Patient/model scores and estimates are copied from existing formal outputs, not recalculated for publication. No final figure-panel mapping is asserted by these semantic filenames. See the [clinical module](../../analysis/08_clinical_residual_tnbc/README.md).
