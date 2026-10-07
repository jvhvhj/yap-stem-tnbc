# Data

This repository distributes **code and study-defined gene definitions, not primary data**.
Raw sequencing files, expression matrices, Seurat/AnnData objects, spatial image
archives and native CNA matrices are not included and must not be uploaded to
GitHub. Large public datasets are not duplicated here — download them from the
sources below.

## Cohorts analysed in this repository

| Dataset | Source | Accession / URL | Role in the study |
| --- | --- | --- | --- |
| Wu breast cancer atlas | GEO | [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) | Discovery cohort: 8 donors, 10,836 malignant epithelial cells |
| ARTEMIS (Yan et al.) | NCBI BioProject | [PRJNA1041570](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1041570) — *"Chemotherapy Response in Triple-negative breast cancer defined by single cell sequencing"* | Independent TNBC validation (78-patient association; 77-patient program transfer) |
| Orthogonal transfer cohort | GEO | [GSE180286](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180286) | Single-patient orthogonal transfer; not multi-patient replication |
| Primary spatial cohort | GEO | [GSE210616](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210616) | Spatial recurrence: 22 patients, 43 sections — see [spatial reproduction](../docs/spatial_reproduction.md) |
| Independent spatial cohort (BSW2) | Zenodo | [10.5281/zenodo.15247064](https://doi.org/10.5281/zenodo.15247064) — archive `10x.visium.tar.gz`, MD5 `dfc29e45270167d126dddb41094aa966` | Nine-patient spatial replication — see [spatial reproduction](../docs/spatial_reproduction.md) |
| BSW2, later deposit version | Zenodo | [10.5281/zenodo.15252874](https://doi.org/10.5281/zenodo.15252874) | Same BSW archive checksum as the record above; additionally deposits data not used here |

## Reported supporting cohorts (no code in this release)

The manuscript reports these as supporting evidence. **No analysis code for them
is included in this repository**, and their prepared objects are not distributed
here. They are listed so the data provenance is not ambiguous.

| Dataset | Source | Accession / URL | Reported role |
| --- | --- | --- | --- |
| Chen breast atlas | CELLxGENE | dataset [`de5416ef-bbfa-41f3-99d4-01a3554173f5`](https://cellxgene.cziscience.com/collections/de5416ef-bbfa-41f3-99d4-01a3554173f5) | Supporting subtype-level association |
| SCAN-B | GEO | [GSE96058](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058) | Supporting population-scale bulk association |
| METABRIC | cBioPortal (EGAS00000000083) | [brca_metabric](https://www.cbioportal.org/study/summary?id=brca_metabric) | Supporting population-scale bulk association |
| LUAD | GEO | [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907) | Cross-cancer specificity / generic-program attenuation |
| CRC | GEO | [GSE132465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465) | Cross-cancer specificity / generic-program attenuation |
| HNSCC | GEO | [GSE103322](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322) | Cross-cancer specificity / generic-program attenuation |

## Reference resources

| Resource | Source | Use |
| --- | --- | --- |
| MSigDB Hallmark gene sets | [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb) — `h.all.v2024.Hs.symbols.gmt` | Functional-program and enrichment definitions |
| HGNC gene nomenclature | [HGNC download archive](https://www.genenames.org/download/archive/) | Approved symbol/alias mapping |
| UCSC hg38 refGene | [UCSC downloads](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/) | Gene order for inferCNV |

## What is included here

Small files needed to inspect the analysis definitions:

- `config/signatures/` — YAP17 (17 genes), Stem21 (21 genes), Program146 (146 genes)
- `config/datasets.tsv`, `config/input_manifest.tsv`
- `analysis/02_malignant_cell_and_cnv/config/patient_inputs.tsv` — anonymised
  patient-count configuration

## Status of each cohort's input chain

For each cohort, an independent reader needs to know whether the prepared
analysis object can be rebuilt from the public download.

| Cohort | Prepared input included | Preparation documented | State |
| --- | --- | --- | --- |
| GSE176078 (Wu) | No | Yes — [population extraction](../analysis/02_malignant_cell_and_cnv/scripts/extract_wu_malignant_population.R), [scoring](../analysis/03_yap_stem_continuum/scripts/score_yap17_stem21.R) and [state assignment](../analysis/03_yap_stem_continuum/scripts/derive_joint_axis_and_tertiles.R) are all scripted and verified | Reproducible from the atlas metadata given the source object |
| ARTEMIS / PRJNA1041570 | No | Partly — the preparation script exists, but the download-to-prepared-object transformation is not fully recorded | Requires a recorded construction chain |
| GSE180286 | No | Partly — a preparation script exists, but the sample-to-patient mapping is unresolved | Requires mapping completion |
| GSE210616 | No | Yes — layout and entry scripts are defined in [spatial reproduction](../docs/spatial_reproduction.md); the download-to-workspace assembly step is not yet scripted | Partially closed |
| BSW2 | No | Yes — source, checksum, expected layout and entry script are defined in [spatial reproduction](../docs/spatial_reproduction.md) | Partially closed |
| Chen / SCAN-B / METABRIC / LUAD / CRC / HNSCC | No | No | Out of code-release scope |

## Notes

- An accession is an access pointer. It is not proof that a prepared analysis
  object is identical to the public file, and it does not by itself grant
  permission to redistribute a derived object.
- Local derived object names used during the analysis must **not** be presented
  as a public deposit.
- Prepared objects need either an appropriate public deposit with access
  conditions or a verified construction script. Neither a placeholder directory
  nor a source-accession link closes that requirement.
- Third-party gene sets, annotation resources and software carry their own
  licences. See [`docs/third_party_code_and_attribution.md`](../docs/third_party_code_and_attribution.md).
