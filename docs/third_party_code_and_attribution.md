# Third-party software and code attribution

Using a package, adapting a function, and redistributing a package's source are different activities. Package citations do not by themselves establish permission to copy source code, and the manuscript's eventual code license cannot replace third-party terms.

## Findings from the available sources

| Component | Evidence | Treatment in this snapshot |
| --- | --- | --- |
| CopyKAT | Analysis scripts call `copykat::copykat`; the project stores package-related results and wrappers | Wrapper code is recovered; exact successful-run package revision remains to be documented. No CopyKAT package source is bundled |
| SCEVAN | Calls `SCEVAN::pipelineCNA`. An upstream source snapshot exists locally at commit `5a49b88ac9445eeffcebb95404e3190992faac04` | No package source is bundled. Local metadata notes an Encoding-field addition to DESCRIPTION; preserve that modification record if redistributing the patched source later |
| yaGST | Local third-party dependency source, version 2017.08.25 | Not bundled. Upstream DESCRIPTION contains both a license placeholder and `GPL (>= 3)`; confirm the applicable terms before redistributing a local copy |
| inferCNV | Analysis wrapper calls package APIs; successful sessions record infercnv 1.18.1 | Wrapper retained. Package source/native matrices not bundled |
| ESTIMATE-style spatial scoring | Local worker defines `estimate_scores_fast` and loads ESTIMATE resources; an explicit attribution trail for that function is not recovered | Worker withheld pending code-origin, algorithm equivalence and permission review. This is not a finding of plagiarism |
| sctransform | Spatial preprocessing source-loads an external source tree and compiled component | Do not present those components as newly authored code; record exact version/source and platform build before packaging |
| Seurat, DESeq2, edgeR, fgsea, plotting and Python libraries | Imported libraries in project wrappers | Record software citations and versions. Do not copy installed package trees into this repository |
| Yan study repository | Consulted for organization of analysis documentation and supplementary material | No scientific source or supplementary file was copied from that reference during this repository preparation |

The available text/source scan is not an exhaustive originality certification. Unattributed copied fragments cannot reliably be ruled out by keyword search. The authors should review the retained scripts and identify any functions adapted from tutorials, repositories or collaborators that are not yet documented.

## What to record for any copied or adapted implementation

- Original repository/file, author and stable commit or release.
- Applicable license and required copyright/NOTICE files.
- Whether the code is used unchanged or modified; describe the changes.
- Which local script uses it and why.
- Whether redistribution is authorized. Where unclear, withhold the source pending clarification rather than deleting its attribution.

These requirements apply separately to code and external data/gene-set resources. Scientific names or public author contacts in upstream attribution must not be mistaken for manuscript authorship.

## Sources inspected

- [SCEVAN source DESCRIPTION at the recorded commit](https://github.com/AntonioDeFalco/SCEVAN/blob/5a49b88ac9445eeffcebb95404e3190992faac04/DESCRIPTION).
- [yaGST DESCRIPTION](https://github.com/miccec/yaGST/blob/master/DESCRIPTION).
- [CopyKAT](https://github.com/navinlabcode/copykat).
- [inferCNV](https://github.com/broadinstitute/infercnv).
- [Reference study's analysis organization](https://github.com/navinlabcode/tnbc-chemo).
- [GitHub guidance on repository licensing](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/licensing-a-repository).

No conclusion about license compatibility is asserted here. The final repository license remains an author/institutional decision after this review.
