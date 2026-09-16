#!/usr/bin/env python
# Purpose: Prepare ARTEMIS data and frozen signatures
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 5; orthogonal transfer in Fig. 4.
"""Stages 0-2: Yan h5ad audit, frozen signatures, overlap audit, and subset extraction.

The source h5ad is always opened read-only.  Only the cancer-cell rows and the
predefined genes required by the frozen analysis are materialized.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import logging
import math
import os
import sys
from pathlib import Path

import anndata as ad
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(os.environ.get("AHIPPO_YAP_ROOT", Path.cwd()))
OUT = ROOT / "0716_yan2026_validation"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
REPORTS = OUT / "reports"
INPUTS = OUT / "inputs"
LOGS = OUT / "logs"
for directory in (TABLES, FIGURES, REPORTS, INPUTS, LOGS):
    directory.mkdir(parents=True, exist_ok=True)

LOG_FILE = LOGS / "01_prepare_yan_and_frozen_signatures.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[logging.FileHandler(LOG_FILE, encoding="utf-8"), logging.StreamHandler(sys.stdout)],
)
LOG = logging.getLogger("prepare")

H5AD = Path(os.environ.get("YAN_H5AD", ROOT / "scouter" / "af8c4fce-4c63-4671-b339-91a383cf36f6.h5ad"))
FROZEN_SCORE = ROOT / "0703_rebuild" / "continuous_state_main_figures" / "core_score_gene_lists.csv"
FROZEN_MODULES = ROOT / "0713_score_independent_rebuild" / "tables" / "module_gene_overlap.tsv"
FROZEN_REP = ROOT / "0713_score_independent_rebuild" / "tables" / "Figure4_v3_1_panelD_lollipop_plot_data.tsv"
FROZEN_FORMULA = ROOT / "0710_final_evidence_rebuild" / "tables_final" / "Fig2_score_definition_used.txt"
CELL_CYCLE = INPUTS / "seurat_cell_cycle_genes_updated_2019.tsv"
HGNC = INPUTS / "hgnc_complete_set_2026-07-16.txt"

REQUIRED = [H5AD, FROZEN_SCORE, FROZEN_MODULES, FROZEN_REP, FROZEN_FORMULA, CELL_CYCLE, HGNC]


def sha256(path: Path, chunk: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, na_rep="NA")
    LOG.info("Wrote %s (%s rows)", path, len(df))


def value_summary(series: pd.Series, max_examples: int = 12) -> str:
    vals = series.astype("string").fillna("<NA>").value_counts(dropna=False)
    return " | ".join(f"{idx}:{int(count)}" for idx, count in vals.head(max_examples).items())


def fail_block(reason: str, details: list[str]) -> None:
    report = REPORTS / "BLOCKED_stage0_yan2026.md"
    report.write_text(
        "# ARTEMIS input validation failure\n\n"
        f"**Reason:** {reason}\n\n"
        + "\n".join(f"- {x}" for x in details)
        + "\n\nNo formal score or outcome analysis was run.\n",
        encoding="utf-8",
    )
    LOG.error("Input validation failed: %s", reason)
    raise SystemExit(2)


def contiguous_runs(rows: np.ndarray) -> list[tuple[int, int, int, int]]:
    """Return (output_start, output_end_exclusive, h5_start, h5_end_exclusive)."""
    breaks = np.flatnonzero(np.diff(rows) > 1)
    starts = np.r_[0, breaks + 1]
    ends = np.r_[breaks + 1, len(rows)]
    return [(int(s), int(e), int(rows[s]), int(rows[e - 1]) + 1) for s, e in zip(starts, ends)]


def extract_csr_rows_columns(
    h5_path: Path, rows: np.ndarray, columns: np.ndarray, group: str = "X"
) -> sparse.csr_matrix:
    """Read selected columns for selected CSR rows with bounded memory.

    The public anndata backed sparse indexer scans very slowly for a two-axis
    fancy selection on this 907-million-nonzero object.  This routine reads the
    CSR arrays for 97 contiguous cancer-cell runs, filters columns in memory,
    and never materializes the full matrix.
    """
    rows = np.asarray(rows, dtype=np.int64)
    columns = np.asarray(columns, dtype=np.int64)
    order = np.argsort(columns)
    sorted_cols = columns[order]
    n_out = len(rows)
    out_indptr = np.zeros(n_out + 1, dtype=np.int64)
    out_indices_parts: list[np.ndarray] = []
    out_data_parts: list[np.ndarray] = []
    cursor = 0

    with h5py.File(h5_path, "r") as handle:
        grp = handle[group]
        shape = tuple(int(x) for x in grp.attrs["shape"])
        if shape[0] <= rows.max() or shape[1] <= columns.max():
            raise IndexError("Requested row/column exceeds h5ad matrix shape")
        indptr = grp["indptr"][:]
        data_ds = grp["data"]
        index_ds = grp["indices"]
        for run_no, (out_s, out_e, h5_s, h5_e) in enumerate(contiguous_runs(rows), 1):
            p0, p1 = int(indptr[h5_s]), int(indptr[h5_e])
            run_data = data_ds[p0:p1]
            run_cols = index_ds[p0:p1]
            for local_row, global_row in enumerate(range(h5_s, h5_e)):
                d0 = int(indptr[global_row] - p0)
                d1 = int(indptr[global_row + 1] - p0)
                cols_row = run_cols[d0:d1]
                vals_row = run_data[d0:d1]
                pos = np.searchsorted(sorted_cols, cols_row)
                keep = (pos < len(sorted_cols)) & (sorted_cols[np.minimum(pos, len(sorted_cols) - 1)] == cols_row)
                if np.any(keep):
                    # ``order`` maps sorted-column positions back to the
                    # requested/output column order.
                    mapped = order[pos[keep]].astype(np.int32, copy=False)
                    idx_order = np.argsort(mapped)
                    out_indices_parts.append(mapped[idx_order])
                    out_data_parts.append(vals_row[keep][idx_order].astype(np.float32, copy=False))
                    cursor += int(np.sum(keep))
                out_indptr[out_s + local_row + 1] = cursor
            if run_no % 10 == 0 or run_no == len(contiguous_runs(rows)):
                LOG.info("Extracted cancer-cell CSR run %s/%s", run_no, len(contiguous_runs(rows)))

    out_indices = np.concatenate(out_indices_parts) if out_indices_parts else np.array([], dtype=np.int32)
    out_data = np.concatenate(out_data_parts) if out_data_parts else np.array([], dtype=np.float32)
    return sparse.csr_matrix((out_data, out_indices, out_indptr), shape=(n_out, len(columns)))


def zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    sd = np.nanstd(values, ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(values)
    return (values - np.nanmean(values)) / sd


def mean_columns(matrix: sparse.csr_matrix, positions: list[int]) -> np.ndarray:
    if not positions:
        return np.full(matrix.shape[0], np.nan)
    return np.asarray(matrix[:, positions].mean(axis=1)).ravel()


def main() -> None:
    missing = [str(p) for p in REQUIRED if not p.exists()]
    if missing:
        fail_block("Required input file is missing", missing)

    manifest = []
    for path in REQUIRED:
        manifest.append(
            {
                "input": path.name,
                "path": str(path),
                "bytes": path.stat().st_size,
                "modified": pd.Timestamp(path.stat().st_mtime, unit="s").isoformat(),
                "sha256": sha256(path) if path.stat().st_size < 50_000_000 else "not_computed_large_input",
                "role": "read_only_source",
            }
        )
    write_tsv(pd.DataFrame(manifest), INPUTS / "stage0_2_input_manifest.tsv")

    LOG.info("Opening Yan h5ad in backed mode: %s", H5AD)
    adata = ad.read_h5ad(H5AD, backed="r")
    obs = adata.obs.copy()
    var = adata.var.copy()
    n_cells, n_genes = adata.shape

    obs_summary = pd.DataFrame(
        {
            "column": obs.columns,
            "dtype": [str(obs[c].dtype) for c in obs.columns],
            "n_unique_including_na": [int(obs[c].nunique(dropna=False)) for c in obs.columns],
            "n_missing": [int(obs[c].isna().sum()) for c in obs.columns],
            "value_summary": [value_summary(obs[c]) for c in obs.columns],
        }
    )
    write_tsv(obs_summary, TABLES / "yan2026_obs_column_summary.tsv")
    write_tsv(
        pd.DataFrame(
            {
                "field": var.columns,
                "dtype": [str(var[c].dtype) for c in var.columns],
                "n_unique_including_na": [int(var[c].nunique(dropna=False)) for c in var.columns],
                "n_missing": [int(var[c].isna().sum()) for c in var.columns],
            }
        ),
        TABLES / "yan2026_var_fields.tsv",
    )

    required_obs = ["donor_id", "pCR_status", "author_cell_type", "cell_type", "cell_state", "nCount_RNA", "nFeature_RNA"]
    missing_obs = [c for c in required_obs if c not in obs.columns]
    if missing_obs:
        fail_block("Required patient/response/cancer/technical fields are ambiguous or missing", missing_obs)

    cancer_author = obs["author_cell_type"].astype(str).eq("Tumor").to_numpy()
    cancer_standard = obs["cell_type"].astype(str).eq("abnormal cell").to_numpy()
    label_discordance = int(np.sum(cancer_author != cancer_standard))
    patient = obs["donor_id"].astype(str)
    response = obs["pCR_status"].astype(str)
    patient_n = int(patient.nunique())
    cancer_n = int(cancer_author.sum())
    cancer_patients = int(patient[cancer_author].nunique())

    cancer_counts = (
        pd.DataFrame({"patient": patient, "is_cancer": cancer_author, "response": response})
        .groupby("patient", sort=True)
        .agg(total_cells=("is_cancer", "size"), cancer_cells=("is_cancer", "sum"), pCR_status=("response", "first"))
        .reset_index()
    )
    cancer_counts["has_cancer_cells"] = cancer_counts["cancer_cells"] > 0
    write_tsv(cancer_counts, TABLES / "yan2026_cancer_cell_counts_by_patient_verified.tsv")

    response_counts = (
        cancer_counts.groupby("pCR_status", dropna=False)
        .agg(total_patients=("patient", "size"), patients_with_cancer=("has_cancer_cells", "sum"), cancer_cells=("cancer_cells", "sum"))
        .reset_index()
    )
    write_tsv(response_counts, TABLES / "yan2026_response_counts_verified.tsv")

    with h5py.File(H5AD, "r") as handle:
        layer_keys = list(handle["layers"].keys()) if "layers" in handle else []
        x_encoding = handle["X"].attrs.get("encoding-type", "unknown")
        raw_encoding = handle["raw"]["X"].attrs.get("encoding-type", "unknown")
        h5_keys = list(handle.keys())

    sample_idx = np.flatnonzero(cancer_author)[:5]
    sample_x = adata.X[sample_idx, :]
    sample_raw = adata.raw.X[sample_idx, :]
    if hasattr(sample_x, "to_memory"):
        sample_x = sample_x.to_memory()
    if hasattr(sample_raw, "to_memory"):
        sample_raw = sample_raw.to_memory()
    raw_integer = bool(np.allclose(sample_raw.data, np.round(sample_raw.data)))
    x_noninteger = bool(not np.allclose(sample_x.data, np.round(sample_x.data)))
    x_library = np.array([np.expm1(sample_x.getrow(i).data).sum() for i in range(sample_x.shape[0])])
    x_log1p_10k = bool(np.all((x_library > 9500) & (x_library < 10500)))
    raw_sums = np.asarray(sample_raw.sum(axis=1)).ravel()
    ncount_sample = obs.iloc[sample_idx]["nCount_RNA"].to_numpy(float)
    raw_matches_ncount = bool(np.all(np.abs(raw_sums - ncount_sample) / np.maximum(ncount_sample, 1) < 0.02))

    audit_rows = [
        ("object", "local_cells", n_cells, 427_823, n_cells == 427_823, "local h5ad"),
        ("object", "paper_cells", n_cells, 427_857, n_cells == 427_857, "paper; local difference=-34"),
        ("object", "local_minus_paper_cells", n_cells - 427_857, -34, n_cells - 427_857 == -34, "local vs paper"),
        ("object", "genes", n_genes, 32_354, n_genes == 32_354, "local h5ad"),
        ("patient", "total_patients", patient_n, 101, patient_n == 101, "donor_id"),
        ("cancer", "cancer_cells", cancer_n, 49_275, cancer_n == 49_275, "author_cell_type=Tumor"),
        ("cancer", "patients_with_cancer", cancer_patients, 97, cancer_patients == 97, "donor_id among Tumor"),
        ("cancer", "label_discordance", label_discordance, 0, label_discordance == 0, "Tumor vs abnormal cell"),
        ("expression", "X_encoding", x_encoding, "csr_matrix", str(x_encoding) == "csr_matrix", "h5py attrs"),
        ("expression", "X_dtype", str(adata.X.dtype), "floating log-normalized", x_noninteger, "sampled Tumor rows"),
        ("expression", "X_log1p_library_approximately_10000", x_log1p_10k, True, x_log1p_10k, "sum(expm1(X))"),
        ("expression", "raw_encoding", raw_encoding, "csr_matrix", str(raw_encoding) == "csr_matrix", "h5py attrs"),
        ("expression", "raw_integer_counts", raw_integer, True, raw_integer, "sampled raw.X"),
        ("expression", "raw_matches_nCount_within_2pct", raw_matches_ncount, True, raw_matches_ncount, "sampled raw.X"),
        ("expression", "layers", ",".join(layer_keys) if layer_keys else "none", "none", len(layer_keys) == 0, "h5py layers"),
        ("metadata", "patient_field", "donor_id", "unambiguous", patient_n == 101, "obs"),
        ("metadata", "response_field", "pCR_status", "pCR/RD/Excluded", set(response.unique()) == {"pCR", "RD", "Excluded"}, "obs"),
        ("metadata", "cancer_field", "author_cell_type=Tumor", "cross-checked", label_discordance == 0, "obs"),
        ("metadata", "sample_batch_field", "orig.ident", "7 source batches", obs["orig.ident"].nunique() == 7, "obs"),
    ]
    audit = pd.DataFrame(audit_rows, columns=["section", "item", "observed", "expected", "pass", "source"])
    write_tsv(audit, TABLES / "yan2026_metadata_audit.tsv")

    blocking_items = audit.loc[
        audit["item"].isin(
            [
                "total_patients",
                "cancer_cells",
                "patients_with_cancer",
                "label_discordance",
                "X_log1p_library_approximately_10000",
                "raw_integer_counts",
                "response_field",
            ]
        )
        & ~audit["pass"].astype(bool)
    ]
    if len(blocking_items):
        fail_block("Stage 0 gate failed", blocking_items.to_string(index=False).splitlines())

    structure_report = f"""# Yan 2026 h5ad structure audit

## Gate decision

**PASS.** The expression layer, cancer-cell label, patient identifier and response field are sufficiently defined for formal validation.

## Object structure

- File: `{H5AD}`
- Local object: **{n_cells:,} cells x {n_genes:,} genes**.
- Paper total: **427,857 cells**; the CELLxGENE-curated local object has **34 fewer cells**. The local file does not retain a row-level exclusion ledger, so the exact identities of the 34 cells cannot be reconstructed. This is reported as a curation/version difference, not silently reconciled.
- Top-level h5ad keys: `{', '.join(h5_keys)}`.
- `X`: CSR, `{adata.X.dtype}`; sampled cells satisfy `sum(expm1(X)) ~= 10,000`, establishing log1p library-size-normalized expression.
- `raw.X`: CSR, `{adata.raw.X.dtype}`; sampled values are integer counts and agree with `nCount_RNA` within 2%.
- `layers`: none.
- `obsm`: `{', '.join(adata.obsm.keys())}`; no spatial coordinates are present.
- `uns`: `{', '.join(adata.uns.keys())}`.

## Metadata mapping

- Patient: `donor_id` ({patient_n} patients).
- Source/batch: `orig.ident` ({obs['orig.ident'].nunique()} values); this is not substituted for patient.
- Response: `pCR_status` with pCR, RD and Excluded.
- Cancer cells: `author_cell_type == Tumor`, exactly cross-checked against `cell_type == abnormal cell`.
- Cancer-cell total: {cancer_n:,}; patients with cancer cells: {cancer_patients}.
- Metaprogram/dominant state: `cell_state`; M01-M13 labels are present as categorical assignments, but continuous author module scores are not stored.

## Formal-analysis read policy

The complete 427k x 32k matrix is never materialized.  The validation reads the 49,275 cancer-cell rows and only the frozen signature, representative-gene, and technical-covariate columns from `X`.
"""
    (REPORTS / "yan2026_h5ad_structure.md").write_text(structure_report, encoding="utf-8")
    LOG.info("Input integrity checks completed successfully")

    # ------------------------------------------------------------------
    # Stage 1: frozen definitions and HGNC-aware feature resolution.
    # ------------------------------------------------------------------
    score = pd.read_csv(FROZEN_SCORE)
    score["module"] = score["axis"].replace({"YAP": "YAP clean", "Stemness": "Stemness clean"})
    score["signature_type"] = "score"
    score["source"] = str(FROZEN_SCORE)

    overlap_src = pd.read_csv(FROZEN_MODULES, sep="\t")
    requested_modules = [
        "UPR",
        "TNFA NFKB",
        "Hypoxia",
        "Adhesion Remodeling",
        "Wound Healing",
        "Survival Stress",
        "Anoikis Resistance",
        "Integrin Adhesion",
    ]
    module_sets: dict[str, list[str]] = {}
    module_original: dict[str, list[str]] = {}
    for module in requested_modules:
        row = overlap_src.loc[overlap_src["module"].eq(module)].iloc[0]
        retained = [g for g in str(row["retained_genes"]).split(",") if g and g != "nan"]
        removed = [g for g in str(row["removed_genes"]).split(",") if g and g != "nan"]
        module_sets[module] = retained
        module_original[module] = removed + retained

    rep = pd.read_csv(FROZEN_REP, sep="\t")
    representative_genes = rep["gene"].dropna().astype(str).drop_duplicates().tolist()
    score_genes = set(score["gene"].astype(str))

    shared_gene_modules: dict[str, set[str]] = {}
    for module, genes in module_sets.items():
        shared_gene_modules[module] = {
            gene
            for gene in genes
            if sum(gene in other_genes for other_genes in module_sets.values()) > 1
        }
    module_unique = {module: [g for g in genes if g not in shared_gene_modules[module]] for module, genes in module_sets.items()}

    manifest_rows: list[dict] = []
    for _, row in score.iterrows():
        manifest_rows.append(
            {
                "signature_type": "score",
                "module": row["module"],
                "gene": row["gene"],
                "in_yap_stem_score": True,
                "score_independent_gene": False,
                "pairwise_unique_gene": False,
                "definition": "arithmetic mean of log-normalized expression",
                "source": str(FROZEN_SCORE),
            }
        )
    for module in requested_modules:
        for gene in module_original[module]:
            manifest_rows.append(
                {
                    "signature_type": "module",
                    "module": module,
                    "gene": gene,
                    "in_yap_stem_score": gene in score_genes,
                    "score_independent_gene": gene in module_sets[module],
                    "pairwise_unique_gene": gene in module_unique[module],
                    "definition": "v3.1 frozen module; score-independent version removes YAP/Stem score genes",
                    "source": str(FROZEN_MODULES),
                }
            )
    for gene in representative_genes:
        manifest_rows.append(
            {
                "signature_type": "representative_gene",
                "module": rep.loc[rep["gene"].eq(gene), "module"].iloc[0],
                "gene": gene,
                "in_yap_stem_score": gene in score_genes,
                "score_independent_gene": gene not in score_genes,
                "pairwise_unique_gene": False,
                "definition": "Figure 4 v3.1 panel D frozen plotting gene",
                "source": str(FROZEN_REP),
            }
        )
    manifest_df = pd.DataFrame(manifest_rows).drop_duplicates(["signature_type", "module", "gene"])

    formula_text = FROZEN_FORMULA.read_text(encoding="utf-8", errors="replace").strip().replace("\n", " | ")
    formula_rows = pd.DataFrame(
        [
            {"signature_type": "formula", "module": "Joint YAP-Stem axis", "gene": "NA", "in_yap_stem_score": True,
             "score_independent_gene": False, "pairwise_unique_gene": False,
             "definition": "z(YAP clean score) + z(Stemness clean score); pooled cancer-cell z scores; pooled tertiles",
             "source": f"{FROZEN_FORMULA}; {formula_text}"}
        ]
    )
    manifest_df = pd.concat([manifest_df, formula_rows], ignore_index=True)
    write_tsv(manifest_df, TABLES / "frozen_signature_manifest.tsv")

    hgnc = pd.read_csv(HGNC, sep="\t", dtype=str, low_memory=False).fillna("")
    hgnc_by_symbol = hgnc.set_index("symbol", drop=False)
    yan_symbols = var["gene_symbols"].astype(str).to_numpy()
    symbol_to_col: dict[str, int] = {}
    for idx, symbol in enumerate(yan_symbols):
        symbol_to_col.setdefault(symbol, idx)

    all_frozen = sorted(set(manifest_df.loc[manifest_df["gene"].ne("NA"), "gene"].astype(str)))
    cc = pd.read_csv(CELL_CYCLE, sep="\t")
    technical_genes = cc["gene"].dropna().astype(str).drop_duplicates().tolist()
    requested_genes = all_frozen + [g for g in technical_genes if g not in all_frozen]

    resolution_rows = []
    resolved_symbol: dict[str, str | None] = {}
    for gene in requested_genes:
        approved = gene
        alias_symbols = ""
        previous_symbols = ""
        if gene in hgnc_by_symbol.index:
            rec = hgnc_by_symbol.loc[gene]
            if isinstance(rec, pd.DataFrame):
                rec = rec.iloc[0]
            approved = str(rec.get("symbol", gene))
            alias_symbols = str(rec.get("alias_symbol", ""))
            previous_symbols = str(rec.get("prev_symbol", ""))
        candidates = [gene, approved]
        candidates += [x for x in alias_symbols.split("|") if x]
        candidates += [x for x in previous_symbols.split("|") if x]
        hit = next((x for x in candidates if x in symbol_to_col), None)
        resolved_symbol[gene] = hit
        resolution_rows.append(
            {
                "requested_gene": gene,
                "hgnc_approved_symbol": approved,
                "gene_alias": alias_symbols,
                "previous_symbol": previous_symbols,
                "yan_resolved_symbol": hit if hit else "NA",
                "present": hit is not None,
                "yan_var_index": symbol_to_col.get(hit, pd.NA) if hit else pd.NA,
                "role": "technical_cell_cycle" if gene in technical_genes and gene not in all_frozen else "frozen_analysis",
            }
        )
    resolution = pd.DataFrame(resolution_rows)
    write_tsv(resolution, TABLES / "yan2026_gene_symbol_resolution.tsv")

    frozen_missing = resolution.loc[resolution["role"].eq("frozen_analysis") & ~resolution["present"], "requested_gene"].tolist()
    if frozen_missing:
        LOG.warning("Frozen genes absent after HGNC resolution: %s", ", ".join(frozen_missing))

    extracted_requested = [g for g in requested_genes if resolved_symbol[g] is not None]
    extracted_symbols = [resolved_symbol[g] for g in extracted_requested]
    extracted_cols = np.array([symbol_to_col[s] for s in extracted_symbols], dtype=np.int64)
    cancer_rows = np.flatnonzero(cancer_author)
    LOG.info("Extracting %s genes for %s cancer cells from X", len(extracted_cols), len(cancer_rows))
    expr = extract_csr_rows_columns(H5AD, cancer_rows, extracted_cols, "X")
    sparse.save_npz(INPUTS / "yan_cancer_frozen_gene_expression.npz", expr, compressed=True)
    write_tsv(
        pd.DataFrame(
            {
                "matrix_column": np.arange(len(extracted_requested)),
                "requested_gene": extracted_requested,
                "resolved_symbol": extracted_symbols,
                "h5ad_var_index": extracted_cols,
            }
        ),
        INPUTS / "yan_cancer_frozen_gene_expression_columns.tsv",
    )
    position = {gene: idx for idx, gene in enumerate(extracted_requested)}

    cancer_obs = obs.iloc[cancer_rows].copy()
    cell_meta = pd.DataFrame(
        {
            "cell_id": cancer_obs.index.astype(str),
            "patient": cancer_obs["donor_id"].astype(str).to_numpy(),
            "response": cancer_obs["pCR_status"].astype(str).to_numpy(),
            "batch": cancer_obs["orig.ident"].astype(str).to_numpy(),
            "cell_state": cancer_obs["cell_state"].astype(str).to_numpy(),
            "nCount_RNA": cancer_obs["nCount_RNA"].to_numpy(float),
            "nFeature_RNA": cancer_obs["nFeature_RNA"].to_numpy(float),
        }
    )

    yap_genes = score.loc[score["axis"].eq("YAP"), "gene"].astype(str).tolist()
    stem_genes = score.loc[score["axis"].eq("Stemness"), "gene"].astype(str).tolist()
    cell_meta["YAP_score"] = mean_columns(expr, [position[g] for g in yap_genes if g in position])
    cell_meta["Stemness_score"] = mean_columns(expr, [position[g] for g in stem_genes if g in position])
    cell_meta["YAP_Stem_axis"] = zscore(cell_meta["YAP_score"].to_numpy()) + zscore(cell_meta["Stemness_score"].to_numpy())
    q1, q2 = np.quantile(cell_meta["YAP_Stem_axis"], [1 / 3, 2 / 3])
    cell_meta["YS_state"] = pd.cut(
        cell_meta["YAP_Stem_axis"], bins=[-np.inf, q1, q2, np.inf], labels=["Low", "Intermediate", "High"], include_lowest=True
    ).astype(str)

    for module, genes in module_sets.items():
        cell_meta[f"module__{module}"] = mean_columns(expr, [position[g] for g in genes if g in position])
        unique_genes = module_unique[module]
        cell_meta[f"module_unique__{module}"] = mean_columns(expr, [position[g] for g in unique_genes if g in position])

    for phase in ["S", "G2M"]:
        genes = cc.loc[cc["phase"].eq(phase), "gene"].astype(str).tolist()
        cell_meta[f"{phase}_score"] = mean_columns(expr, [position[g] for g in genes if g in position])
    for gene in representative_genes:
        cell_meta[f"gene__{gene}"] = expr[:, position[gene]].toarray().ravel() if gene in position else np.nan

    cell_meta.to_csv(INPUTS / "yan_cancer_cell_scores.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")
    LOG.info("Wrote cancer-cell score cache with %s rows and %s columns", len(cell_meta), cell_meta.shape[1])

    # Coverage: overall and patient-resolved detection for every frozen gene.
    frozen_gene_order = sorted(set(all_frozen))
    cov_rows = []
    for gene in frozen_gene_order:
        resolved = resolved_symbol[gene]
        approved_row = resolution.loc[resolution["requested_gene"].eq(gene)].iloc[0]
        if resolved is None or gene not in position:
            cov_rows.append(
                {
                    "gene": gene,
                    "hgnc_approved_symbol": approved_row["hgnc_approved_symbol"],
                    "gene_alias": approved_row["gene_alias"],
                    "yan_resolved_symbol": "NA",
                    "present": False,
                    "scope": "overall",
                    "patient": "ALL",
                    "n_cells": len(cell_meta),
                    "n_detected": 0,
                    "detection_rate": 0.0,
                }
            )
            continue
        detected = expr[:, position[gene]].toarray().ravel() > 0
        cov_rows.append(
            {
                "gene": gene,
                "hgnc_approved_symbol": approved_row["hgnc_approved_symbol"],
                "gene_alias": approved_row["gene_alias"],
                "yan_resolved_symbol": resolved,
                "present": True,
                "scope": "overall",
                "patient": "ALL",
                "n_cells": len(cell_meta),
                "n_detected": int(detected.sum()),
                "detection_rate": float(detected.mean()),
            }
        )
        for pt, idx in cell_meta.groupby("patient", sort=True).indices.items():
            di = detected[np.asarray(idx)]
            cov_rows.append(
                {
                    "gene": gene,
                    "hgnc_approved_symbol": approved_row["hgnc_approved_symbol"],
                    "gene_alias": approved_row["gene_alias"],
                    "yan_resolved_symbol": resolved,
                    "present": True,
                    "scope": "patient",
                    "patient": pt,
                    "n_cells": len(di),
                    "n_detected": int(di.sum()),
                    "detection_rate": float(di.mean()),
                }
            )
    coverage = pd.DataFrame(cov_rows)
    write_tsv(coverage, TABLES / "yan2026_signature_gene_coverage.tsv")

    # Stage 2 overlap audit on score-independent module lists.
    overlap_rows = []
    for m1 in requested_modules:
        for m2 in requested_modules:
            s1, s2 = set(module_sets[m1]), set(module_sets[m2])
            shared = sorted(s1 & s2)
            union = s1 | s2
            overlap_rows.append(
                {
                    "module_1": m1,
                    "module_2": m2,
                    "n_genes_1": len(s1),
                    "n_genes_2": len(s2),
                    "shared_gene_count": len(shared),
                    "jaccard_index": len(shared) / len(union) if union else np.nan,
                    "shared_genes": ",".join(shared),
                }
            )
    pairwise = pd.DataFrame(overlap_rows)
    write_tsv(pairwise, TABLES / "module_pairwise_overlap.tsv")

    score_independent_lists = []
    unique_lists = []
    for module in requested_modules:
        for gene in module_sets[module]:
            score_independent_lists.append({"module": module, "gene": gene, "source_version": "YAP_Stem_score_genes_removed"})
        if module_unique[module]:
            for gene in module_unique[module]:
                unique_lists.append({"module": module, "gene": gene, "status": "retained_after_pairwise_shared_gene_removal"})
        else:
            unique_lists.append({"module": module, "gene": "NA", "status": "empty_after_pairwise_shared_gene_removal"})
    write_tsv(pd.DataFrame(score_independent_lists), TABLES / "module_score_independent_gene_lists.tsv")
    write_tsv(pd.DataFrame(unique_lists), TABLES / "module_pairwise_shared_removed_gene_lists.tsv")

    mat = pairwise.pivot(index="module_1", columns="module_2", values="jaccard_index").loc[requested_modules, requested_modules]
    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, ax = plt.subplots(figsize=(7.4, 6.4))
    im = ax.imshow(mat.to_numpy(), vmin=0, vmax=1, cmap="Reds")
    ax.set_xticks(range(len(requested_modules)), requested_modules, rotation=45, ha="right")
    ax.set_yticks(range(len(requested_modules)), requested_modules)
    for i in range(len(requested_modules)):
        for j in range(len(requested_modules)):
            value = mat.iloc[i, j]
            ax.text(j, i, f"{value:.2f}", ha="center", va="center", color="white" if value > 0.55 else "#202428", fontsize=7)
    ax.set_title("Frozen score-independent module gene overlap", weight="bold", pad=10)
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cb.set_label("Jaccard index")
    fig.tight_layout()
    fig.savefig(FIGURES / "module_gene_overlap_heatmap.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / "module_gene_overlap_heatmap.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES / "module_gene_overlap_heatmap.svg", bbox_inches="tight")
    plt.close(fig)

    prep_qc = pd.DataFrame(
        [
            {"stage": 0, "check": "h5ad_backed_read", "status": "PASS", "details": "read_h5ad(backed='r'); full matrix not materialized"},
            {"stage": 0, "check": "expression_layer", "status": "PASS", "details": "X log1p normalized; raw.X integer counts"},
            {"stage": 0, "check": "cancer_label", "status": "PASS", "details": f"Tumor equals abnormal cell; discordance={label_discordance}"},
            {"stage": 0, "check": "patient_response_fields", "status": "PASS", "details": "donor_id and pCR_status verified"},
            {"stage": 1, "check": "frozen_definition_source", "status": "PASS", "details": "read from v3.1 source tables; no Yan-driven selection"},
            {"stage": 1, "check": "frozen_gene_coverage", "status": "PASS" if not frozen_missing else "WARN", "details": f"missing={','.join(frozen_missing) if frozen_missing else 'none'}"},
            {"stage": 2, "check": "pairwise_unique_modules", "status": "WARN" if any(len(v) == 0 for v in module_unique.values()) else "PASS", "details": "; ".join(f"{k}:{len(v)}" for k, v in module_unique.items())},
        ]
    )
    write_tsv(prep_qc, TABLES / "stage0_2_QC.tsv")

    adata.file.close()
    LOG.info("Stages 0-2 completed successfully")


if __name__ == "__main__":
    main()
