# Supplementary tables

The S1–S5 publication workbooks are recorded in `table_code_map.tsv`. S1 contains cohort/sample structure, S2 gene memberships and scoring definitions, S3 ARTEMIS patient-level results, S4 spatial results, and S5 the Seo residual-TNBC clinical extension. S2–S5 use the `PUBLICATION_CLEAN.xlsx` suffix.

S2 separates the broad six-program memberships used for discovery/ARTEMIS continuous associations from the compact study-defined modules used for spatial analyses. The five spatial modules exclude YAP17/Stem21 genes and retain Program146 overlap. Compact-module High–Low sensitivities and broad-program continuous correlations are distinct estimands.

S3 documents endpoint-specific adjustment: primary YAP17–Stem21 associations use log1p counts and raw detected-gene counts; Program146 and the six-program associations use log1p for both. All models also include S and G2M scores. For six-program associations, Joint and each program are separately residualized, then correlated with Spearman's rho.

S4 separates patient-level spatial program effects from cohort summaries. Spatial program scores average within-section gene-wise z scores of log1p sctransform-corrected UMI counts. Primary adjusted associations use average-rank residualization for TumorPurity and log1p raw depth, followed by Pearson residual correlation. BSW2 coverage refers to the nine analyzed sections.

S5 distinguishes source-provided D/ND outcome definitions from this study's D=1/ND=0 coding and RCB-II model reference category.

`source_data/SourceData_Seo_*.tsv` contains 48 patient scores, all seven available formal model-coefficient records, 146 mapping records (141 measurable) and one nested-model comparison. No clinical data are added to S3/S4. No supplementary figure is added or renumbered. Figure-panel mapping remains to be supplied when manuscript placement is decided.

`Supplementary_Table_S1_FINAL.xlsx` and `Supplementary_Table_S2_FINAL.xlsx` are compatibility copies of S1 and S2, not separate scientific versions.

Additional cohort/sample inventories, gene-mapping records, robustness/null results and cross-context validation retain their original estimates and analysis scope.

Internal version records and checksums are not included in the publication workbooks.
