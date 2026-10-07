# Data availability

| Dataset/resource | Public source | Required representation | Analysis module |
| --- | --- | --- | --- |
| Wu breast cancer atlas | [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) | Author cell annotations, gene-level expression, fixed TNBC cell identities | Data preparation, continuum, Program146 and CNA |
| ARTEMIS / Yan | NCBI BioProject [PRJNA1041570](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1041570) — *"Chemotherapy Response in Triple-negative breast cancer defined by single cell sequencing"*; [study code and data-access index](https://github.com/navinlabcode/tnbc-chemo) | Author-annotated single-cell object, patient IDs and cancer-state labels | Independent TNBC validation |
| Orthogonal transfer cohort | [GSE180286](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180286) | Gene expression and exact primary/LN sample-to-patient mapping | Orthogonal program transfer; not multi-patient replication |
| Primary spatial cohort | [GSE210616](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210616) | Count matrices, Visium coordinates, section/patient registry and public H&E context | Spatial validation |
| Independent spatial cohort (BSW2) | Zenodo [10.5281/zenodo.15247064](https://doi.org/10.5281/zenodo.15247064) — archive `10x.visium.tar.gz`, MD5 `dfc29e45270167d126dddb41094aa966`. The later record [10.5281/zenodo.15252874](https://doi.org/10.5281/zenodo.15252874) carries the same BSW archive checksum | Nine included patient samples, 10x expression and coordinates | Independent spatial replication |
| Hallmark gene sets | [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb) | `h.all.v2024.Hs.symbols.gmt`; use the recorded resource version | External pathway profiling |
| Gene nomenclature | [HGNC](https://www.genenames.org/download/archive/) | Recorded HGNC symbol/alias table | Approved unambiguous mappings only |
| Genomic order | [UCSC hg38 refGene](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/) | Recorded gene-order file and generation provenance | inferCNV |
| Seo residual-TNBC clinical cohort | NCBI SRA BioProject [PRJNA1256162](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1256162); [Seo et al. 2025](https://doi.org/10.1016/j.xcrm.2025.102164) | 48 deidentified patients, 48 fresh-frozen post-NAC residual-tumor samples, 48 paired-end WTS runs; Neo HR-negative/HER2-negative, D34 / ND14; operative-subtype sensitivity n45 | `analysis/08_clinical_residual_tnbc/` |

The Seo analysis addresses clinical association in a retrospective extreme-outcome sample, not prospective validation or prediction in an unselected survival cohort. Supplementary workbooks and processed result tables are not distributed in this code repository; raw reads and raw Salmon directories are also excluded. Gene definitions, input metadata and software records required to interpret the code remain included.

The detailed existing filenames are in `config/input_manifest.tsv`. Large data and restricted/third-party resource files are not redistributed here. Accession pages are access pointers, not proof that the prepared analysis objects are identical to public files.

The expected input layout, checksum and analysis entry for both spatial cohorts are in [`spatial_reproduction.md`](spatial_reproduction.md). What remains open is the assembly step that turns a fresh public download into those prepared workspaces; current local derived object names must not be described as a public deposit. No download URL is invented.
