# Workstation CNA environment

## Recovered execution records

- Operating system: Ubuntu 24.04.2 LTS, x86_64 Linux.
- R: 4.3.3.
- SCEVAN: 1.0.3.
- infercnv: 1.18.1.
- Matrix: 1.6-5 in these successful-run session records.
- Patient-wise execution; the returned SCEVAN/inferCNV status records report seed 20260728 and 12 cores.

These values come from the returned patient-session files, not the current Windows installation. The native inferCNV archives for two patients are truncated after readable early members. Session text is recoverable, but the archives must not be described as complete native-result deposits.

## Inputs

Eight patients contribute 10,836 fixed malignant epithelial cells and 28,771 same-patient reference cells, totaling 39,607 cells in the CNA input-preparation manifest. CopyKAT's reference-subsampled input is a different caller-specific object: do not claim that its sampled run contains all 39,607 cells.

Patient-level raw counts, exact cell-ID annotations, reference groups and the recorded hg38 gene-order file are required. Raw counts and native large caller objects are not included in GitHub.

## Recovered script parameters

The SCEVAN source uses anchored reference cells, `SUBCLONES=FALSE`, `beta_vega=0.5`, `ClonalCN=TRUE`, `FIXED_NORMAL_CELLS=FALSE`, human genes, 12 cores and seed 20260728. Optional reference subsampling is controlled in the source; its actual successful-run setting must be recovered from execution records.

The inferCNV source uses cutoff 0.1, clustered groups, denoising and 12 cores. Continuous burden is mean absolute deviation from neutral signal 1. The descriptive high-CNA threshold is strictly above the same-patient reference burden's 95th percentile. A separate native i6 HMM run is retained as state-level evidence; it is not converted into an invented malignant-cell binary call.

The recovered CopyKAT script uses human symbol IDs, `ngene.chr=5`, `win.size=25`, `KS.cut=0.1`, Euclidean distance, and `genome="hg20"`. Reference cap, seed and core count are caller arguments. Those parser defaults do not substitute for each successful patient's parameter log.

## Interpretation of returned calls

CopyKAT `not.defined` and SCEVAN `filtered` are not evaluable, not negative. inferCNV-high is descriptive RNA-CNA support, not a validated malignant classification. A reference high-CNA fraction near 5% follows from the threshold construction and is not 95% specificity. CID3963 HMM non-testability under the patient-specific reference limitation is not a negative biological result.

## Missing execution linkage

The exact scripts sent to the workstation are locally recoverable. No returned execution-time source copy or checksum has yet demonstrated whether those scripts were changed on the workstation. Recovering that link, as described in [the reproduction guide](../docs/reproducibility.md), is a precondition for any full-rerun claim. No destructive RAM-disk cleanup wrapper has been made a publication runner.
