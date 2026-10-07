# Documented clinical analysis environments

Quantification and clinical modelling have separate environment records. No single unified environment is asserted.

| Component | Version/reference | Role |
|---|---|---|
| Salmon | 1.10.3 | Paired-end WTS quantification |
| Genome | GRCh38 | Complete genome decoys |
| Annotation | GENCODE v36 | Transcripts and exact transcript-to-gene mapping |
| R | 4.3.3 | Formal clinical model runtime |
| logistf | 1.26.1 | Firth models and nested comparison |
| Operating system for R models | Ubuntu 24.04.2 LTS | Documented Linux runtime |

The R/logistf versions are recovered from the formal clinical environment record. The Salmon index and quantification options follow the documented project specification: full-genome decoy-aware index, k=31, `--gencode`; quantification `-l A --seqBias --gcBias`. The public wrapper is a portable implementation of that specification, not a recovered original execution log. No package installation was used to infer historical versions.

## References

- Seo ES et al. Cell Reports Medicine (2025). https://doi.org/10.1016/j.xcrm.2025.102164
- Patro R et al. Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods (2017). https://doi.org/10.1038/nmeth.4197
- Frankish A et al. GENCODE reference annotation for the human and mouse genomes. Nucleic Acids Research (2019). https://doi.org/10.1093/nar/gky955
- Heinze G, Schemper M. A solution to the problem of separation in logistic regression. Statistics in Medicine (2002). https://doi.org/10.1002/sim.1047
- Heinze G, Ploner M, Jiricka L, Steiner G (2025). logistf: Firth's Bias-Reduced Logistic Regression. R package 1.26.1. https://georgheinze.r-universe.dev/logistf/citation.html
- Tirosh I et al. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. Science (2016). https://doi.org/10.1126/science.aad0501
