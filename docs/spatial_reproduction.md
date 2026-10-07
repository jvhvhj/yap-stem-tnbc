# Spatial reproduction: sources and input preparation

This page defines the expected on-disk layout and the entry scripts for the two
spatial cohorts. It does not restate the Methods. Every path below is the one the
shipped script actually reads; nothing here is inferred from a filename.

The spatial module has two routes:

| Route | Cohort | Role in the study |
| --- | --- | --- |
| A | GSE210616 | Primary spatial recurrence cohort (22 patients / 43 sections) |
| B | BSW2 (Zenodo) | Independent spatial replication cohort (9 patients) |

Both routes consume the same study-defined definitions and share the scoring rule; they
differ in source, normalisation implementation and coordinate universe.

study-defined inputs for both routes:

- `config/signatures/YAP17.tsv` (17 genes), `config/signatures/Stem21.tsv` (21 genes)
- `config/signatures/Program146.tsv` (146 genes)
- the authorised alias map `CYR61→CCN1`, `CTGF→CCN2`, `HNRNPU-AS1→HNRNPU`

Program146 is disjoint from the 38 score-defining genes; the scripts assert this
rather than assuming it.

---

## Route A — GSE210616

### Public source

GEO series [GSE210616](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210616).
Per-GSM processed files are sufficient; raw FASTQ was not deposited by the
original authors for privacy reasons.

| Item | Value |
| --- | --- |
| Patients / sections | 22 / 43 |
| Per-GSM expression | `filtered_feature_bc_matrix.h5` (36,601 features; raw UMI counts) |
| Per-GSM coordinates | `tissue_positions_list.csv.gz` |
| Per-GSM image (context only) | `tissue_hires_image.png.gz` |
| Processing | Space Ranger 1.1.0, GRCh38 |
| Spot denominator | 56,567 in-tissue spots over the 43 sections |

The documented section/patient registries used by the scripts are the documented
project registries, not a re-derivation of the paper's sample table.

### Expected prepared layout

Working directory = an isolated reproduction workspace (not the repository root).

```
.runtime/Figure6_S0Q2_geo/<sample_id>/filtered_feature_bc_matrix.h5
.runtime/Figure6_S1_inputs/section_registry.tsv
.runtime/Figure6_S1_inputs/patient_registry.tsv
.runtime/Figure6_S1_inputs/score_gene_list.tsv
.runtime/Figure6_S1_inputs/146_gene_program.txt
.runtime/Figure6_S1_inputs/Q3_ESTIMATE_spot_scores.tsv
.runtime/Figure6_S1_inputs/source_hashes.tsv
```

- `<sample_id>` is the GSM identifier from `section_registry.tsv`; exactly one
  `filtered_feature_bc_matrix.h5` per sample directory.
- The worker asserts the H5 spot count equals the registry's
  `n_expression_spots`, so a mismatched download fails loudly rather than
  silently subsetting.
- `score_gene_list.tsv` and `146_gene_program.txt` carry the study-defined memberships
  listed above; the worker re-asserts 17 / 21 / 146 and the zero-overlap rule.
- The ESTIMATE-derived spot covariates (`Q3_ESTIMATE_spot_scores.tsv`) are used
  only as a study-defined adjustment variable. Their code origin and permission
  are still under review — see the open items below.

Environment variable `SPATIAL_WORK_ROOT` may be used to point at this workspace;
the manifest row `GSE210616_GEO_FILES` records the same requirement.

### Scoring

`analysis/07_spatial_validation/scripts/compute_spatial_scores_worker.R <worker_id> <n_workers>`

Per section: raw UMI → sctransform `vst` (seeded) → `log1p` corrected counts →
YAP17 and Stem21 means over the measurable members of each signature →
cohort-level z of each mean → `Joint = z(YAP) + z(Stem)`. Program146 is scored as
the mean of per-gene row z scores, with the 142-gene variant computed alongside.

The sctransform implementation used here is source-loaded from
`.runtime/Figure6_S0Q3_sctransform_src` (R sources plus a compiled shared
library); it is not the installed CRAN package. Reproduce that source tree from
its own upstream before running.

### Analysis entry

1. `compute_spatial_scores_worker.R <worker_id> <n_workers>` — writes one chunk
   per section with a per-chunk QC status.
2. `assemble_spatial_coupling_results.R` — joins chunks to the registries and the
   ESTIMATE covariates, and produces the patient- and section-level association
   tables.
3. `compute_spatial_enrichment_and_neighborhood_statistics.R` — cohort recurrence
   and Lee L neighbourhood statistics on the study-defined first-order Visium graph.
4. `score_spatial_modules.R` then `assemble_spatial_module_results.R` — the
   patient-by-functional-module association table.

Rendering scripts (`render_spatial_score_maps.R`, `render_spatial_robustness_panels.py`,
`render_spatial_enrichment_panel.R`, `render_neighborhood_coorganization_panel.R`,
`render_patient_recurrent_functional_coupling.R`) read the assembled tables.

---

## Route B — BSW2 independent replication

### Public source

Zenodo record [10.5281/zenodo.15247064](https://doi.org/10.5281/zenodo.15247064)
holds the BSW2 Visium archive `10x.visium.tar.gz`. The later version of the same
deposit, [10.5281/zenodo.15252874](https://doi.org/10.5281/zenodo.15252874),
carries the identical BSW2 archive (same MD5) and additionally deposits unrelated
data that were not used here.

| Item | Value |
| --- | --- |
| Archive | `10x.visium.tar.gz` |
| MD5 | `dfc29e45270167d126dddb41094aa966` |
| SHA256 | `48a4a92ce5708e7d4b341abc8e02ee7e79ea259b4f6c145674d48290c5bcf19f` |
| Size | 2,295,846,195 bytes |
| Technology | 10x Visium, TNBC FFPE |
| Patients analysed | 9 (sAA1, sAA2, sAA6, sAA8, sEA1, sEA2, sEA5, sEA6, sEA7) |
| Patients in archive but excluded | sAA9 — present in the archive, not in the analysis sample |

Verify the checksum before unpacking. The archive is not redistributed here and
is not a GitHub asset.

### Expected prepared layout

```
Figure6G0_BSW2_archive_audit/Figure6G0_BSW2_section_registry.tsv
Figure6G0_BSW2_archive_audit/audit/extracted_objects/10x.visium/<section_id>/outs/filtered_feature_bc_matrix.h5
Figure6G0_BSW2_archive_audit/audit/extracted_objects/10x.visium/<section_id>/outs/spatial/tissue_positions.csv
Figure6G0_BSW2_archive_audit/audit/extracted_objects/10x.visium/<section_id>/outs/spatial/tissue_hires_image.png
```

`<section_id>` is the archive member name and equals the patient identifier;
this cohort contributes one section per patient. The registry is the documented
section table and is required — the script reads it before touching any H5.
If the extracted tree is placed elsewhere, use `BSW2_WORK_ROOT`.

### Scoring

`analysis/07_spatial_validation/scripts/bsw2_spatial_replication.R`

Normalisation is a deliberate second implementation, not a copy of route A:
`log2(raw UMI / library size × 6000 + 1)`, formula-matched to
`Giotto::normalizeGiotto` (`norm_methods = "standard"`, `library_size_norm = TRUE`,
`scalefactor = 6000`, `log_norm = TRUE`, `log_offset = 1`, `logbase = 2`) with a
recorded deterministic-equivalence check. YAP17/Stem21 means, the Joint axis and
the Program146 score then follow the same rule as route A.

The script asserts `setequal()` between the H5 barcodes and the `in_tissue == 1`
coordinate barcodes, and asserts the spot count against the registry, so a
mis-unpacked or partially extracted section fails rather than degrading quietly.

### Analysis entry

`bsw2_spatial_replication.R`, run from the workspace root. It writes the raw and
depth-/featureness-adjusted patient associations used by the Figure 6G render
(`render_independent_spatial_replication.R`) and the supplementary spatial table.

---

## Open items on both routes

These are recorded, not solved, and they bound what "reproduction" can currently mean:

1. **Prepared-input construction is not fully scripted.** The workspace trees
   above are the contract the workers read; the step that assembles them from a
   fresh download is not yet a shipped, verified script.
2. **ESTIMATE-style covariates.** The code origin and redistribution permission of
   the tumour-composition scoring function used for `Q3_ESTIMATE_spot_scores.tsv`
   are still under review, so that worker is withheld.
3. **Source-loaded sctransform.** Route A loads a compiled implementation from
   `.runtime/`, outside this repository. Obtain it from its own upstream.
4. **Cohort-size wording.** The BSW2 article reports 9 BSW2 spatial patients; the
   Zenodo archive lists 10 samples. The nine analysed here are the ones the
   registry marks as included, with sAA9 explicitly excluded.

No spatial score, patient effect or significance test was recomputed during
repository preparation. The result tables referenced above are existing outputs.

See [the reproduction guide](reproducibility.md) for the repository-wide
boundary and [the input manifest](../config/input_manifest.tsv) for machine-readable
rows.
