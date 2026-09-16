# Software environments

Historical environments and the environment used to inspect code are different records. `package_versions.tsv` lists package versions recovered from successful Linux patient-run session files. `R_sessionInfo.txt` reproduces representative historical session records, with provenance identified by archive/member name.

## Confirmed workstation records

The returned SCEVAN/inferCNV sessions report R 4.3.3 on Ubuntu 24.04.2 LTS. SCEVAN 1.0.3 and infercnv 1.18.1 are recorded in those sessions. The older local Windows records describe earlier attempts and must not replace these successful-run versions.

CopyKAT's exact successful-run package version has not been established from a patient-specific runtime record. A current installed package version or a Windows plotting-session version would not resolve that question.

## Dependency installation

The dependency list in `script_package_dependencies.tsv` is extracted from source imports. It is not an installation lockfile. Missing versions are not filled with currently available releases.

The project was not shown to have been developed under `renv`; no reconstructed lockfile is presented as an original one. Python environment reconstruction remains incomplete. `requirements.txt` documents this rather than pretending that an unpinned install recreates the historical environment.

For the static checks only, use Python 3 and base R. For biological reproduction, recover the missing historical environments and package revisions first. Do not install vendored third-party sources without checking their provenance and terms.

See [workstation details](workstation_cna_environment.md) and [third-party source attribution](../docs/third_party_code_and_attribution.md).
