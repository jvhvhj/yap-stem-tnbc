# YAP–Stem transcriptional continuum in triple-negative breast cancer

## About

Analysis code and study-defined gene definitions for a manuscript investigating a YAP–Stem transcriptional continuum in malignant epithelial cells from triple-negative breast cancer (TNBC).

**Code author and maintainer: Yuting Zhang**

This local source snapshot is **not yet a complete, end-to-end reproducible release**. Some original execution records, input-construction steps and final figure compositions remain unresolved. The [reproducibility documentation](docs/reproducibility.md) states these limitations explicitly. No analysis results were recalculated during repository preparation.

## Study overview

- YAP17 and Stem21 scores define the transcriptional axes.
- Patient-level analyses characterize the recurrence and heterogeneity of their association.
- Program146 is a score-gene-excluded, patient-consistent transcriptional program.
- Independent single-cell and spatial datasets provide complementary tests of the study-defined framework.
- RNA-derived CNA callers provide malignant-identity sensitivity evidence without redefining the primary cell population.

The study does not establish a universal pan-cancer mechanism, a clinical biomarker or a drug-response predictor.

## Repository structure

| Directory | Contents |
| --- | --- |
| `analysis/` | Scientific-question modules and recovered analysis/rendering scripts |
| `config/` | Exact signatures and input descriptions |
| `data/` | Data access and redistribution policy |
| `supp/` | Supplementary table descriptions and source mapping |
| `environment/` | Historical software records and dependencies |
| `docs/` | Reproduction boundaries, figure/code mapping and attribution |
| `tests/` | Static checks; not a rerun of biological results |

## Analysis index

See [analyses and related manuscript content](analysis/README.md).

The [Program146 module](analysis/04_program146/README.md) now includes the recovered July 27 original, a separately identified August 8 regeneration source, and a portable derivation with input-only preflight. Program146 source provenance is resolved. Full numerical execution and the complete repository's input and environment dependencies remain unverified.

## Data availability

The clinical extension uses Seo et al. (2025), NCBI SRA BioProject [PRJNA1256162](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1256162): 48 post-NAC residual-TNBC patients (D=34, ND=14), with a 45-patient operative HR-negative/HER2-negative sensitivity subset. This is an independent retrospective extreme-outcome cohort. Existing results support association with the early distant-recurrence phenotype after accounting for RCB class, not absolute recurrence-risk prediction or a validated clinical biomarker.

The [clinical workflow](analysis/08_clinical_residual_tnbc/README.md) is SRA reads → Salmon → exact GENCODE v36 transcript-to-gene aggregation → Program146 scoring → clinical merge → Firth logistic regression. Existing processed clinical tables are supplied without rerunning models. Clinical software records are documented separately from single-cell/CNA environments.

See [dataset access](docs/data_availability.md) and the [input manifest](config/input_manifest.tsv). Raw sequencing files, large expression matrices, Seurat objects, spatial image archives and native CNA matrices are not included. A public accession does not establish permission to redistribute every derived object or third-party resource.

## Reproducibility

Begin with [the reproduction guide](docs/reproducibility.md) and [software environments](environment/README.md). The supported check in this snapshot is:

```bash
python tests/smoke/check_repository.py --rscript Rscript
```

This loads configuration/signature files and parses scripts. It neither executes a scientific analysis nor demonstrates numerical reproduction of the manuscript. No universal run-all command is provided while the missing-source and input-construction issues remain unresolved.

## Citation

See [`CITATION.cff`](CITATION.cff). It gives the repository as a software work and its code author. A machine-readable article citation and DOI will be added once the manuscript citation is finalized; none is asserted here.

## Correspondence

- Tingming Liang: tmliang@njnu.edu.cn
- Li Guo: lguo@njupt.edu.cn

These are authorized correspondence contacts, not a statement of the complete author list or author order.

## License and third-party software

Released under the [MIT License](LICENSE) — Copyright (c) 2026 Yuting Zhang.

The licence covers the author-written analysis and rendering code in this repository only. Third-party software and data retain their own terms. Their authors are not authors of this manuscript merely because their resources are used. See [third-party software and source attribution](docs/third_party_code_and_attribution.md).
