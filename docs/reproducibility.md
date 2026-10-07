# Reproduction guide and current limitations

## What this snapshot provides

The repository contains recovered scripts arranged by scientific question, exact study-defined signature files, historical Linux environment records, input descriptions, and static checks. Archived originals were not modified. Source and release checksums are recorded in `script_provenance.tsv`; this is not a substitute for reproducing results.

The population-to-state derivation is public and verified. The chain a reader
should follow is:

```
source atlas annotation
  -> analysis/02_malignant_cell_and_cnv/scripts/extract_wu_malignant_population.R
       100,064 cells / 26 patients -> TNBC 42,512 / 10 -> malignant 10,836 / 8
  -> analysis/03_yap_stem_continuum/scripts/score_yap17_stem21.R
       YAP17 17/17 and Stem21 21/21 measurable in the RNA data layer
  -> analysis/03_yap_stem_continuum/scripts/derive_joint_axis_and_tertiles.R
       Joint = z(YAP) + z(Stem); High 3,612 / Other 7,224
  -> analysis/03_yap_stem_continuum/scripts/validate_joint_state_reproduction.R
       max |delta Joint| 2.0e-14; 0 state mismatches over 10,836 cells
  -> downstream analyses
```

study-defined validation outputs for the last two steps are kept in
`analysis/03_yap_stem_continuum/tables/`. The Figure 4E/F Hallmark chain is also
closed: `analysis/05_functional_programs/scripts/build_hallmark_ranking.R` builds
the ranking input from the saved score-excluded 17,597-gene Wald table, and the
enrichment is produced by `external_hallmark_analysis.R`. See the module READMEs.

## What remains unresolved

1. **Program146 execution equivalence, not source identity.** The corrected July 27 original source is recovered and archived. The August 8 regeneration also matches its recorded checksum and is separately identified. This provenance blocker is resolved. The portable original module has static and input-only checks; the full biological run and historical environment are not verified. Legacy downstream scripts still refer to old temporary paths and require input-contract review.
2. **Successful CNA execution linkage.** CopyKAT, SCEVAN and inferCNV scripts exist locally. SCEVAN/inferCNV outbound script copies match the outbound transfer archive. The returned patient bundles contain results and session records, not executed-script copies or execution-time checksums. It is not yet established whether workstation-only edits occurred.
3. **Incomplete native inferCNV archives.** Two returned compressed patient archives terminate prematurely. Complete small final burden/status tables exist separately for all eight patients. The intact tables do not prove that every native output and log in those archives was recovered.
4. **Spatial input construction.** The score/summary workers require prepared spatial inputs and source-loaded preprocessing components. The sources, expected layout and entry scripts are documented in `spatial_reproduction.md`; the step that assembles a fresh public download into those workspaces is not yet scripted. The ESTIMATE implementation also needs explicit code-origin/permission review.
5. **Manuscript selection.** Author-supplied approved PDFs establish current content/numbering, including current Fig.4A/B derivation/DE and E/F Hallmark/leading-edge content. Some exact render versions and manual compositions remain unresolved. Filenames and dates cannot close these issues. Breast-subtype, cross-cancer and P6 analyses are not automatically included without confirmation of reporting scope.
6. **Release metadata.** The code licence is settled: MIT, Copyright (c) 2026 Yuting Zhang. `CITATION.cff` records the software author. The manuscript's complete author list, article citation and DOI are still pending and are deliberately not invented; correspondence contacts are provided only as authorized by the authors.

These limitations prevent a claim of “one-command reproduction from GEO download to all spatial/manuscript results.” Documentation alone cannot close missing code, input construction, execution linkage or redistribution permissions.

## Directory and path contract

Scientific modules in `analysis/` have descriptive names. To avoid unverified changes to scientific logic, recovered scripts retain their historical relative input and output filenames inside a **separate analysis workspace**. This is a transitional compatibility contract, not the final distribution layout for deposited data.

Absolute private paths in selected release copies have been replaced with workspace-relative paths. Script directories are not intended to be the working directory. Some scripts also accept `AHIPPO_YAP_ROOT` or `AHIPPO_ROOT`; use only an isolated reproduction workspace, never the original archive. Runtime scratch inputs formerly held in temporary directories now use `.runtime/`; R libraries use the configured R library search path and optional `.software/r-library`.

Rendering helpers that call external tools use `AHIPPO_PYTHON` and `AHIPPO_PDFTOPPM`, falling back to executable lookup. An input manifest is provided, but it is not yet an exhaustive, machine-verified dependency graph. Do not start a full run until the unresolved items above have been closed.

## Permitted test now

From the repository root:

```bash
python tests/smoke/check_repository.py --rscript Rscript
```

The test checks exact membership cardinalities and exclusion, source hashes, configuration readability, Python syntax and R syntax. It does not import scientific Python modules, load expression objects, run R scripts, calculate a score, execute a caller, draw a figure or test a biological hypothesis.

## Conditional CNA invocation contract

The recovered caller parsers use these arguments, with a separately prepared workspace as the working directory:

```bash
# Do not execute until the successful-run scripts/parameters are confirmed.
Rscript /path/to/repository/analysis/02_malignant_cell_and_cnv/scripts/scevan/run_scevan_patient.R CID3963 anchored
Rscript /path/to/repository/analysis/02_malignant_cell_and_cnv/scripts/infercnv/run_infercnv_patient.R CID3963 /path/to/hg38_gene_order.tsv
Rscript /path/to/repository/analysis/02_malignant_cell_and_cnv/scripts/copykat/run_copykat_patient.R CID3963 /path/to/workspace 12 1000 123
```

The CopyKAT argument values above describe the recovered parser's parameter slots/defaults, not a certified invocation of every completed patient. Recover each patient's actual parameter record before reproducing that result. No unified caller dispatcher is supplied while this linkage is unresolved.

## How to close the remaining chain

1. Recover the execution-time CopyKAT, SCEVAN and inferCNV scripts, patient invocation records and successful-run session files from the analysis workstation, together with the two incomplete inferCNV archives if intact copies still exist.
2. Review the recovered original Program146 module, explicit input dependencies and static equivalence records. Recover missing input-construction/environment records, then separately authorize isolated full-run comparison; source recovery does not certify execution.
3. Identify the final manuscript panel and supplementary-table versions; do not choose by timestamp alone.
4. Deposit eligible derived inputs or provide verified construction code and exact public download manifests. Check data/resource redistribution terms separately from code licensing.
5. Resolve third-party code attribution for the withheld ESTIMATE-style scoring function.
6. Complete module-specific input validation and numerical comparisons, then run a clean-environment, end-to-end reproduction outside the original working project. Only that successful run can support a stronger reproduction claim.

The Program146 CLI is run from the repository root with a fresh external output directory. Its exact archive copies intentionally preserve historical absolute paths and must not be executed directly; they are exceptions to the portable-script path policy. See the module README for input-only preflight.

No biological analysis has been rerun as part of this preparation. Missing results are not replaced with reconstructed values.
