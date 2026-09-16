#!/usr/bin/env python
# Purpose: Lineage and external-cohort context summaries
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 4 and Supplementary Fig. S1.
from __future__ import annotations

import hashlib
import math
import os
from pathlib import Path

# Use the operating system environment without a machine-specific WINDIR override.
os.environ.setdefault(
    "MPLCONFIGDIR",
    r".\0723_Core_program_gene_closure_and_multicellular_context_gate\tmp\matplotlib",
)

import anndata as ad
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import statsmodels.api as sm
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy import stats


ROOT = Path(r".")
OUT = ROOT / "0723_Core_program_gene_closure_and_multicellular_context_gate"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
REPORTS = OUT / "reports"
PLOTDATA = OUT / "plotting_data"
LOGS = OUT / "logs"
for d in (TABLES, FIGURES, REPORTS, PLOTDATA, LOGS):
    d.mkdir(parents=True, exist_ok=True)

H5AD = ROOT / "scouter" / "af8c4fce-4c63-4671-b339-91a383cf36f6.h5ad"
YAN_SCORES = ROOT / "0716_yan2026_validation" / "inputs" / "yan_cancer_cell_scores.tsv.gz"
WU_LINEAGE = TABLES / "wu_lineage_cell_scores.tsv.gz"
WU_MODULES = (
    ROOT
    / "0713_score_independent_rebuild"
    / "tables"
    / "module_scores_before_after_removing_score_genes.tsv"
)
LINEAGE_MANIFEST = TABLES / "lineage_signature_manifest.tsv"
LANDSCAPE_SOURCE = (
    ROOT
    / "0720_StageA3_scientific_and_visual_closure_before_v6"
    / "tables"
    / "Figure3_module_joint_decile_plotting_data.tsv"
)
DESEQ_RESULT = TABLES / "locked_DESeq2_High_vs_Low_full.tsv"
FIG5_FROZEN_PLOTDATA = (
    ROOT
    / "0720_module_scale_harmonization_and_affected_panel_rebuild"
    / "tables"
    / "plotting_data"
    / "Figure5_final_v4_1_plotting_data.tsv"
)

MODULES = [
    "UPR",
    "TNFA NFKB",
    "Hypoxia",
    "Survival Stress",
    "Adhesion Remodeling",
    "Wound Healing",
]
MODULE_DISPLAY = {
    "UPR": "UPR",
    "TNFA NFKB": "TNFα–NF-κB",
    "Hypoxia": "Hypoxia",
    "Survival Stress": "Survival stress",
    "Adhesion Remodeling": "Adhesion remodeling",
    "Wound Healing": "Wound healing",
}

LINEAGE_KEYS = {
    "Basal/myoepithelial": "Basal_myoepithelial",
    "Luminal progenitor/secretory luminal": "Luminal_progenitor_secretory_luminal",
    "Mature luminal": "Mature_luminal",
    "Hormone-responsive luminal": "Hormone_responsive_luminal",
}

STATE_MAP = [
    ("Mac-ECM", "mac-ECM", "Mye", "primary_myeloid"),
    ("Mac-angio", "mac-angio", "Mye", "primary_myeloid"),
    ("Mac-CXCL", "mac-M2-CXCL", "Mye", "primary_myeloid"),
    ("Mac-CCL", "mac-M1-CCL", "Mye", "primary_myeloid"),
    ("Mac-IFN", "mac-IFN", "Mye", "primary_myeloid"),
    ("Mac-lip-C1Q", "mac-lip-C1Q", "Mye", "primary_myeloid"),
    ("CAF", "CAFs", "Fibro", "primary_stromal"),
    ("Fibro-matrix", "fibro-matrix", "Fibro", "primary_stromal"),
    ("Fibro-prematrix", "fibro-prematrix", "Fibro", "primary_stromal"),
    ("TEC", "TEC", "Endo", "primary_stromal"),
    ("Endo-prolif", "endo-prolif", "Endo", "primary_stromal"),
    ("Peri-immune", "peri-immune", "Peri", "primary_stromal"),
    ("VSMC-contra", "VSMC-contra", "Peri", "primary_stromal"),
    ("VSMC-synth", "VSMC-synth", "Peri", "primary_stromal"),
    ("CD8-Texh", "CD8-TEXH", "T", "secondary_TNK"),
    ("CD8-TIFN", "CD8-TIFN", "T", "secondary_TNK"),
    ("CD4-Treg", "CD4-TREG", "T", "secondary_TNK"),
    ("CD4-TIFN", "CD4-TIFN", "T", "secondary_TNK"),
    ("NK-CD16high", "NK-CD16high", "T", "secondary_TNK"),
    ("NK-CD16low", "NK-CD16low", "T", "secondary_TNK"),
]

FEATURE_MAP = {
    "Joint axis median": "joint_axis_median",
    "High-state fraction": "high_state_fraction",
    "UPR": "UPR",
    "TNFα–NF-κB": "TNFA NFKB",
    "Hypoxia": "Hypoxia",
    "Survival stress": "Survival Stress",
    "Adhesion remodeling": "Adhesion Remodeling",
    "Wound healing": "Wound Healing",
}

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 7.0,
        "axes.titlesize": 8.0,
        "axes.labelsize": 7.0,
        "xtick.labelsize": 5.8,
        "ytick.labelsize": 5.8,
        "axes.linewidth": 0.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False, na_rep="NA")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def zscore(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    sd = np.nanstd(x, ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(x)
    return (x - np.nanmean(x)) / sd


def bh(p: pd.Series | np.ndarray) -> np.ndarray:
    x = np.asarray(p, dtype=float)
    out = np.full(len(x), np.nan)
    ok = np.isfinite(x)
    if not ok.any():
        return out
    vals = x[ok]
    order = np.argsort(vals)
    ranked = vals[order] * len(vals) / np.arange(1, len(vals) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    temp = np.empty_like(ranked)
    temp[order] = np.minimum(ranked, 1)
    out[ok] = temp
    return out


def safe_spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 5 or np.nanstd(x[ok]) == 0 or np.nanstd(y[ok]) == 0:
        return np.nan, np.nan
    r = stats.spearmanr(x[ok], y[ok])
    return float(r.statistic), float(r.pvalue)


def state_effect(values: np.ndarray, states: np.ndarray) -> float:
    high = values[states == "High"]
    low = values[states == "Low"]
    if len(high) == 0 or len(low) == 0:
        return np.nan
    return float(np.nanmedian(high) - np.nanmedian(low))


def fast_huber_beta(y: np.ndarray, x: np.ndarray, max_iter: int = 20) -> float:
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    ok = np.isfinite(y) & np.all(np.isfinite(x), axis=1)
    y = y[ok]
    x = x[ok]
    if len(y) <= x.shape[1] + 2:
        return np.nan
    try:
        beta = np.linalg.lstsq(x, y, rcond=None)[0]
    except np.linalg.LinAlgError:
        return np.nan
    for _ in range(max_iter):
        resid = y - x @ beta
        scale = np.median(np.abs(resid - np.median(resid))) / 0.6745
        if not np.isfinite(scale) or scale < 1e-8:
            break
        u = resid / (1.345 * scale)
        w = np.ones_like(u)
        mask = np.abs(u) > 1
        w[mask] = 1 / np.abs(u[mask])
        xw = x * np.sqrt(w)[:, None]
        yw = y * np.sqrt(w)
        try:
            new = np.linalg.lstsq(xw, yw, rcond=None)[0]
        except np.linalg.LinAlgError:
            break
        if np.max(np.abs(new - beta)) < 1e-8:
            beta = new
            break
        beta = new
    return float(beta[1])


def fixed_weight_beta(y: np.ndarray, x: np.ndarray, weights: np.ndarray) -> float:
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    weights = np.asarray(weights, dtype=float)
    ok = np.isfinite(y) & np.all(np.isfinite(x), axis=1) & np.isfinite(weights) & (weights > 0)
    if ok.sum() <= x.shape[1] + 1:
        return np.nan
    xx = x[ok]
    yy = y[ok]
    ww = weights[ok]
    a = np.einsum("i,ij,ik->jk", ww, xx, xx)
    b = np.einsum("i,ij,i->j", ww, xx, yy)
    ridge = np.eye(a.shape[0]) * 1e-10
    try:
        return float(np.linalg.solve(a + ridge, b)[1])
    except np.linalg.LinAlgError:
        return np.nan


def fixed_weight_bootstrap_betas(
    y: np.ndarray,
    x: np.ndarray,
    robust_weights: np.ndarray,
    rng: np.random.Generator,
    n_boot: int = 1000,
) -> np.ndarray:
    n = len(y)
    counts = rng.multinomial(n, np.full(n, 1 / n), size=n_boot).astype(float)
    weights = counts * robust_weights[None, :]
    a = np.einsum("bi,ij,ik->bjk", weights, x, x)
    b = np.einsum("bi,ij,i->bj", weights, x, y)
    ridge = np.eye(a.shape[1])[None, :, :] * 1e-10
    out = np.full(n_boot, np.nan)
    try:
        coef = np.linalg.solve(a + ridge, b[:, :, None])[:, :, 0]
        out[:] = coef[:, 1]
    except np.linalg.LinAlgError:
        for i in range(n_boot):
            out[i] = fixed_weight_beta(y, x, weights[i])
    return out


def extract_yan_lineage() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    manifest = pd.read_csv(LINEAGE_MANIFEST, sep="\t")
    yan = pd.read_csv(YAN_SCORES, sep="\t")
    a = ad.read_h5ad(H5AD, backed="r")
    cancer_mask = a.obs["author_cell_type"].astype(str).eq("Tumor").to_numpy()
    idx = np.flatnonzero(cancer_mask)
    cell_ids = a.obs_names[idx].astype(str)
    if len(idx) != 49275:
        raise RuntimeError(f"Expected 49,275 Yan cancer cells, found {len(idx)}")
    if set(cell_ids) != set(yan["cell_id"].astype(str)):
        raise RuntimeError("Yan score cache and author-labeled Tumor cells do not match")

    genes = list(dict.fromkeys(manifest["gene"].astype(str)))
    symbols = a.var["gene_symbols"].astype(str).to_numpy()
    symbol_to_var_index: dict[str, int] = {}
    for i, symbol in enumerate(symbols):
        if symbol not in symbol_to_var_index:
            symbol_to_var_index[symbol] = i
    present = [g for g in genes if g in symbol_to_var_index]
    selected_var_indices = np.array([symbol_to_var_index[g] for g in present], dtype=int)
    gene_index = {g: i for i, g in enumerate(present)}

    # The full X matrix is a 907-million-entry CSR object. Reading a backed
    # fancy row×column slice is slow on Windows, so scan only the 97 contiguous
    # author-labeled Tumor row runs and immediately collapse the fixed gene
    # subset. This never materializes the full matrix.
    runs: list[tuple[int, int, int, int]] = []
    run_start = 0
    for j in range(1, len(idx) + 1):
        if j == len(idx) or idx[j] != idx[j - 1] + 1:
            runs.append((int(idx[run_start]), int(idx[j - 1]), run_start, j))
            run_start = j
    selected_matrix = np.zeros((len(idx), len(present)), dtype=np.float32)
    x_group = a.file["X"]
    n_vars = a.n_vars
    for start, end, out_start, out_end in runs:
        ptr = np.asarray(x_group["indptr"][start : end + 2], dtype=np.int64)
        base = int(ptr[0])
        stop = int(ptr[-1])
        local_ptr = ptr - base
        data = np.asarray(x_group["data"][base:stop])
        indices = np.asarray(x_group["indices"][base:stop], dtype=np.int32)
        block = sp.csr_matrix(
            (data, indices, local_ptr),
            shape=(end - start + 1, n_vars),
        )
        selected_matrix[out_start:out_end] = block[:, selected_var_indices].toarray()
    out = pd.DataFrame({"cell_id": cell_ids})
    coverage = []
    for (signature, version), d in manifest.groupby(["signature", "version"], sort=False):
        defined = list(dict.fromkeys(d["gene"].astype(str)))
        have = [g for g in defined if g in gene_index]
        missing = [g for g in defined if g not in gene_index]
        cols = [gene_index[g] for g in have]
        if len(cols) < 2:
            score = np.full(len(idx), np.nan)
        else:
            score = np.mean(selected_matrix[:, cols], axis=1)
        key = f"lineage__{LINEAGE_KEYS[signature]}__{version}"
        out[key] = zscore(score)
        coverage.append(
            {
                "cohort": "Yan2026",
                "signature": signature,
                "version": version,
                "n_defined": len(defined),
                "n_present": len(have),
                "coverage_fraction": len(have) / len(defined),
                "present_genes": ",".join(have),
                "missing_genes": ",".join(missing),
                "score_algorithm": "arithmetic mean h5ad X log1p-normalized expression followed by cohort z-standardization",
            }
        )
    a.file.close()
    out = yan.merge(out, on="cell_id", how="left", validate="one_to_one")
    out.to_csv(
        TABLES / "yan_lineage_cell_scores.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    return out, pd.DataFrame(coverage), manifest


def prepare_wu() -> pd.DataFrame:
    lineage = pd.read_csv(WU_LINEAGE, sep="\t")
    modules = pd.read_csv(WU_MODULES, sep="\t")
    modules = modules.loc[modules["module"].isin(MODULES)].copy()
    wide = modules.pivot(index="cell_id", columns="module", values="score_after")
    wide.columns = [f"module__{c}" for c in wide.columns]
    out = lineage.merge(wide.reset_index(), on="cell_id", how="left", validate="one_to_one")
    return out


def lineage_analysis_one(cohort: str, d: pd.DataFrame) -> pd.DataFrame:
    if cohort == "Wu2021":
        patient_col, axis_col, state_col = "patient", "YAP_Stem_axis", "YS_state"
    else:
        patient_col, axis_col, state_col = "patient", "YAP_Stem_axis", "YS_state"

    for module in MODULES:
        col = f"module__{module}"
        d[col] = zscore(d[col].to_numpy(float))

    records: list[dict[str, object]] = []
    for version in ("raw_predefined", "core_overlap_removed"):
        lineage_cols = [
            f"lineage__{LINEAGE_KEYS[sig]}__{version}" for sig in LINEAGE_KEYS
        ]
        for patient, g in d.groupby(patient_col, sort=True):
            states = g[state_col].astype(str).to_numpy()
            axis = g[axis_col].to_numpy(float)
            for sig, col in zip(LINEAGE_KEYS, lineage_cols):
                rho, p = safe_spearman(g[col].to_numpy(float), axis)
                records.append(
                    {
                        "record_level": "patient",
                        "analysis": "lineage_joint_correlation",
                        "cohort": cohort,
                        "version": version,
                        "patient": patient,
                        "signature": sig,
                        "module": "",
                        "adjustment": "none",
                        "effect": rho,
                        "P_value": p,
                        "n_cells": len(g),
                    }
                )
                records.append(
                    {
                        "record_level": "patient",
                        "analysis": "lineage_High_minus_Low",
                        "cohort": cohort,
                        "version": version,
                        "patient": patient,
                        "signature": sig,
                        "module": "",
                        "adjustment": "none",
                        "effect": state_effect(g[col].to_numpy(float), states),
                        "P_value": np.nan,
                        "n_cells": len(g),
                    }
                )

            xlin = g[lineage_cols].to_numpy(float)
            xlin = np.column_stack(
                [np.ones(len(g))] + [zscore(xlin[:, j]) for j in range(xlin.shape[1])]
            )
            for module in MODULES:
                y = g[f"module__{module}"].to_numpy(float)
                raw = state_effect(y, states)
                ok = np.isfinite(y) & np.all(np.isfinite(xlin), axis=1)
                if ok.sum() > xlin.shape[1] + 5:
                    beta = np.linalg.lstsq(xlin[ok], y[ok], rcond=None)[0]
                    resid = np.full(len(g), np.nan)
                    resid[ok] = y[ok] - xlin[ok] @ beta
                    adjusted = state_effect(resid, states)
                else:
                    adjusted = np.nan
                records.extend(
                    [
                        {
                            "record_level": "patient",
                            "analysis": "six_core_High_minus_Low",
                            "cohort": cohort,
                            "version": version,
                            "patient": patient,
                            "signature": "four_lineage_signatures",
                            "module": module,
                            "adjustment": "before_lineage_adjustment",
                            "effect": raw,
                            "P_value": np.nan,
                            "n_cells": len(g),
                        },
                        {
                            "record_level": "patient",
                            "analysis": "six_core_High_minus_Low",
                            "cohort": cohort,
                            "version": version,
                            "patient": patient,
                            "signature": "four_lineage_signatures",
                            "module": module,
                            "adjustment": "after_lineage_adjustment",
                            "effect": adjusted,
                            "P_value": np.nan,
                            "n_cells": len(g),
                        },
                    ]
                )

    patient_records = pd.DataFrame(records)
    summaries = []
    group_cols = ["analysis", "cohort", "version", "signature", "module", "adjustment"]
    for keys, g in patient_records.groupby(group_cols, dropna=False, sort=False):
        vals = g["effect"].to_numpy(float)
        vals = vals[np.isfinite(vals)]
        n = len(vals)
        npos = int((vals > 0).sum())
        sign_p = stats.binomtest(npos, n, 0.5).pvalue if n else np.nan
        summaries.append(
            {
                "record_level": "summary",
                **dict(zip(group_cols, keys)),
                "patient": "ALL",
                "effect": float(np.median(vals)) if n else np.nan,
                "P_value": sign_p,
                "n_cells": int(g["n_cells"].sum()),
                "n_patients": n,
                "n_positive": npos,
                "IQR_low": float(np.quantile(vals, 0.25)) if n else np.nan,
                "IQR_high": float(np.quantile(vals, 0.75)) if n else np.nan,
            }
        )
    summary = pd.DataFrame(summaries)
    summary["BH_FDR"] = np.nan
    mask = summary["analysis"].eq("six_core_High_minus_Low")
    for _, idx in summary.loc[mask].groupby(
        ["cohort", "version", "adjustment"], sort=False
    ).groups.items():
        summary.loc[idx, "BH_FDR"] = bh(summary.loc[idx, "P_value"])
    patient_records["n_patients"] = np.nan
    patient_records["n_positive"] = np.nan
    patient_records["IQR_low"] = np.nan
    patient_records["IQR_high"] = np.nan
    patient_records["BH_FDR"] = np.nan
    return pd.concat([patient_records, summary], ignore_index=True, sort=False)


def restrict_and_resummarize_lineage(
    lineage: pd.DataFrame, eligible_yan_patients: set[str]
) -> pd.DataFrame:
    patient_records = lineage.loc[lineage["record_level"].eq("patient")].copy()
    patient_records = patient_records.loc[
        patient_records["cohort"].ne("Yan2026")
        | patient_records["patient"].astype(str).isin(eligible_yan_patients)
    ].copy()
    summaries = []
    group_cols = ["analysis", "cohort", "version", "signature", "module", "adjustment"]
    for keys, g in patient_records.groupby(group_cols, dropna=False, sort=False):
        vals = g["effect"].to_numpy(float)
        vals = vals[np.isfinite(vals)]
        n = len(vals)
        npos = int((vals > 0).sum())
        sign_p = stats.binomtest(npos, n, 0.5).pvalue if n else np.nan
        summaries.append(
            {
                "record_level": "summary",
                **dict(zip(group_cols, keys)),
                "patient": "ALL",
                "effect": float(np.median(vals)) if n else np.nan,
                "P_value": sign_p,
                "n_cells": int(g["n_cells"].sum()),
                "n_patients": n,
                "n_positive": npos,
                "IQR_low": float(np.quantile(vals, 0.25)) if n else np.nan,
                "IQR_high": float(np.quantile(vals, 0.75)) if n else np.nan,
                "BH_FDR": np.nan,
            }
        )
    summary = pd.DataFrame(summaries)
    mask = summary["analysis"].eq("six_core_High_minus_Low")
    for _, idx in summary.loc[mask].groupby(
        ["cohort", "version", "adjustment"], sort=False
    ).groups.items():
        summary.loc[idx, "BH_FDR"] = bh(summary.loc[idx, "P_value"])
    for col in ("n_patients", "n_positive", "IQR_low", "IQR_high", "BH_FDR"):
        patient_records[col] = np.nan
    out = pd.concat([patient_records, summary], ignore_index=True, sort=False)
    out["interpretation_limit"] = (
        "association/confounding audit only; no cell-of-origin or lineage-tracing claim"
    )
    return out


def build_tme_proportions(obs: pd.DataFrame, eligible_patients: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    all_patients = sorted(obs["donor_id"].astype(str).unique())
    total_by_patient = obs.groupby("donor_id", observed=True).size().to_dict()
    parent_counts = obs.groupby(["donor_id", "author_cell_type"], observed=True).size().to_dict()
    state_counts = obs.groupby(["donor_id", "cell_state"], observed=True).size().to_dict()
    cohort_state_counts = obs["cell_state"].astype(str).value_counts().to_dict()
    patient_batch = (
        obs.groupby("donor_id", observed=True)["orig.ident"]
        .agg(lambda x: "|".join(sorted(set(map(str, x)))))
        .to_dict()
    )
    rows = []
    for display, author_state, parent, family in STATE_MAP:
        for patient in all_patients:
            parent_n = int(parent_counts.get((patient, parent), 0))
            state_n = int(state_counts.get((patient, author_state), 0))
            if parent_n < 30:
                status = "NT_parent_lt30"
                prop = np.nan
            elif state_n == 0:
                status = "NT_missing_state_not_zero"
                prop = np.nan
            else:
                status = "testable"
                prop = state_n / parent_n
            rows.append(
                {
                    "patient": patient,
                    "state": display,
                    "author_cell_state": author_state,
                    "parent_compartment": parent,
                    "family": family,
                    "state_cells": state_n,
                    "parent_cells": parent_n,
                    "total_cells": int(total_by_patient.get(patient, 0)),
                    "proportion": prop,
                    "testability_status": status,
                    "batch": patient_batch.get(patient, ""),
                    "cancer_feature_eligible": patient in eligible_patients,
                }
            )
    prop = pd.DataFrame(rows)
    summary = []
    for (state, author_state, parent, family), g in prop.groupby(
        ["state", "author_cell_state", "parent_compartment", "family"], sort=False
    ):
        cohort_total = int(cohort_state_counts.get(author_state, 0))
        n_testable_all = int(g["testability_status"].eq("testable").sum())
        n_testable_eligible = int(
            (g["testability_status"].eq("testable") & g["cancer_feature_eligible"]).sum()
        )
        state_pass = cohort_total >= 100 and n_testable_eligible >= 20
        summary.append(
            {
                "state": state,
                "author_cell_state": author_state,
                "parent_compartment": parent,
                "family": family,
                "cohort_state_cells": cohort_total,
                "n_testable_all_patients": n_testable_all,
                "n_testable_cancer_feature_patients": n_testable_eligible,
                "cohort_total_ge100": cohort_total >= 100,
                "eligible_testable_ge20": n_testable_eligible >= 20,
                "state_testability_pass": state_pass,
            }
        )
    return prop, pd.DataFrame(summary)


def build_cancer_features(yan: pd.DataFrame) -> pd.DataFrame:
    counts = yan["patient"].value_counts()
    eligible = counts[counts >= 50].index
    d = yan.loc[yan["patient"].isin(eligible)].copy()
    rows = []
    for patient, g in d.groupby("patient", sort=True):
        row = {
            "patient": patient,
            "cancer_cells": len(g),
            "batch": "|".join(sorted(set(g["batch"].astype(str)))),
            "joint_axis_median": float(g["YAP_Stem_axis"].median()),
            "high_state_fraction": float(g["YS_state"].astype(str).eq("High").mean()),
        }
        for module in MODULES:
            row[module] = float(g[f"module__{module}"].median())
        rows.append(row)
    out = pd.DataFrame(rows)
    if len(out) != 78:
        raise RuntimeError(f"Expected 78 Yan patients with >=50 cancer cells, found {len(out)}")
    for module in MODULES:
        out[module] = zscore(out[module].to_numpy(float))
    return out


def add_compositional_coordinates(prop: pd.DataFrame) -> pd.DataFrame:
    out = prop.copy()
    out["logit_proportion"] = np.log(
        (out["state_cells"] + 0.5) / (out["parent_cells"] - out["state_cells"] + 0.5)
    )
    out["clr_coordinate"] = np.nan
    map_df = pd.DataFrame(
        STATE_MAP, columns=["state", "author_cell_state", "parent_compartment", "family"]
    )
    for (patient, parent), idx in out.groupby(["patient", "parent_compartment"]).groups.items():
        part = out.loc[idx].copy()
        other = max(
            float(part["parent_cells"].iloc[0] - part["state_cells"].sum()),
            0.0,
        )
        logs = np.log(np.r_[part["state_cells"].to_numpy(float) + 0.5, other + 0.5])
        center = logs.mean()
        out.loc[idx, "clr_coordinate"] = np.log(
            part["state_cells"].to_numpy(float) + 0.5
        ) - center
    out.loc[~out["testability_status"].eq("testable"), ["logit_proportion", "clr_coordinate"]] = np.nan
    return out


def robust_associations(
    prop: pd.DataFrame, feasibility: pd.DataFrame, features: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    p = add_compositional_coordinates(prop)
    p = p.merge(features, on="patient", how="inner", suffixes=("", "_cancer"))
    passing_states = set(
        feasibility.loc[feasibility["state_testability_pass"], "state"].astype(str)
    )
    rows = []
    fit_cache: dict[tuple[str, str], dict[str, object]] = {}
    for state in [x[0] for x in STATE_MAP if x[0] in passing_states]:
        s = p.loc[p["state"].eq(state) & p["testability_status"].eq("testable")].copy()
        family = s["family"].iloc[0]
        for feature_display, feature_col in FEATURE_MAP.items():
            cols = [
                "proportion",
                "logit_proportion",
                "clr_coordinate",
                feature_col,
                "cancer_cells",
                "parent_cells",
                "total_cells",
            ]
            g = s.dropna(subset=cols).copy()
            if len(g) < 20:
                continue
            raw_rho, raw_p = safe_spearman(
                g["proportion"].to_numpy(float), g[feature_col].to_numpy(float)
            )
            logit_rho, logit_p = safe_spearman(
                g["logit_proportion"].to_numpy(float), g[feature_col].to_numpy(float)
            )
            clr_rho, clr_p = safe_spearman(
                g["clr_coordinate"].to_numpy(float), g[feature_col].to_numpy(float)
            )
            x = np.column_stack(
                [
                    zscore(g[feature_col].to_numpy(float)),
                    zscore(np.log10(g["cancer_cells"].to_numpy(float) + 1)),
                    zscore(np.log10(g["parent_cells"].to_numpy(float) + 1)),
                    zscore(np.log10(g["total_cells"].to_numpy(float) + 1)),
                ]
            )
            x = sm.add_constant(x, has_constant="add")
            y = g["logit_proportion"].to_numpy(float)
            try:
                fit = sm.RLM(y, x, M=sm.robust.norms.HuberT()).fit(maxiter=100)
                beta = float(fit.params[1])
                se = float(fit.bse[1])
                pval = float(fit.pvalues[1])
                robust_weights = np.asarray(fit.weights, dtype=float)
            except Exception:
                beta, se, pval = np.nan, np.nan, np.nan
                robust_weights = np.ones(len(y), dtype=float)
            rows.append(
                {
                    "state": state,
                    "family": family,
                    "cancer_feature": feature_display,
                    "feature_column": feature_col,
                    "n_patients": len(g),
                    "raw_spearman_rho": raw_rho,
                    "raw_spearman_P": raw_p,
                    "logit_spearman_rho": logit_rho,
                    "logit_spearman_P": logit_p,
                    "CLR_spearman_rho": clr_rho,
                    "CLR_spearman_P": clr_p,
                    "adjusted_robust_beta": beta,
                    "adjusted_robust_SE": se,
                    "adjusted_robust_P": pval,
                    "adjusted_outcome": "logit(state/parent proportion)",
                    "adjusted_predictor_scale": "per 1 SD cancer feature",
                    "technical_covariates": "log10 cancer cells; log10 parent cells; log10 total cells",
                }
            )
            fit_cache[(state, feature_display)] = {
                "data": g,
                "x": x,
                "y": y,
                "full_beta": beta,
                "robust_weights": robust_weights,
            }
    assoc = pd.DataFrame(rows)
    assoc["raw_BH_FDR_within_family"] = np.nan
    assoc["adjusted_BH_FDR_within_family"] = np.nan
    for _, idx in assoc.groupby("family", sort=False).groups.items():
        assoc.loc[idx, "raw_BH_FDR_within_family"] = bh(
            assoc.loc[idx, "raw_spearman_P"]
        )
        assoc.loc[idx, "adjusted_BH_FDR_within_family"] = bh(
            assoc.loc[idx, "adjusted_robust_P"]
        )
    assoc["raw_adjusted_direction_concordant"] = (
        np.sign(assoc["raw_spearman_rho"]) == np.sign(assoc["adjusted_robust_beta"])
    )

    stability_rows = []
    for r in assoc.itertuples(index=False):
        cache = fit_cache[(r.state, r.cancer_feature)]
        g = cache["data"].reset_index(drop=True)
        x = np.asarray(cache["x"], dtype=float)
        y = np.asarray(cache["y"], dtype=float)
        beta = float(cache["full_beta"])
        robust_weights = np.asarray(cache["robust_weights"], dtype=float)
        seed = int(
            hashlib.sha256(f"{r.state}|{r.cancer_feature}".encode()).hexdigest()[:8], 16
        )
        rng = np.random.default_rng(seed)
        boot = fixed_weight_bootstrap_betas(y, x, robust_weights, rng, n_boot=1000)
        finite = boot[np.isfinite(boot)]
        boot_stability = (
            float(np.mean(np.sign(finite) == np.sign(beta))) if len(finite) else np.nan
        )
        boot_low = float(np.quantile(finite, 0.025)) if len(finite) else np.nan
        boot_high = float(np.quantile(finite, 0.975)) if len(finite) else np.nan

        lopo = []
        lopo_patients = []
        for i, patient in enumerate(g["patient"].astype(str)):
            keep = np.arange(len(g)) != i
            lopo.append(fixed_weight_beta(y[keep], x[keep], robust_weights[keep]))
            lopo_patients.append(patient)
        lopo = np.asarray(lopo, dtype=float)
        valid_lopo = lopo[np.isfinite(lopo)]
        lopo_stable = bool(
            len(valid_lopo)
            and np.all(np.sign(valid_lopo) == np.sign(beta))
        )
        if len(valid_lopo):
            changes = np.abs(valid_lopo - beta)
            max_i = int(np.nanargmax(changes))
            max_change = float(changes[max_i])
            influence_ratio = max_change / max(abs(beta), 1e-8)
            influential_patient = lopo_patients[np.flatnonzero(np.isfinite(lopo))[max_i]]
        else:
            max_change, influence_ratio, influential_patient = np.nan, np.nan, ""
        no_single_dominance = bool(
            np.isfinite(influence_ratio) and influence_ratio <= 1.0 and lopo_stable
        )

        batch_values = g["batch"].astype(str).to_numpy()
        unique_batches = sorted(set(batch_values))
        batch_betas = []
        for batch in unique_batches:
            keep = batch_values != batch
            if keep.sum() >= 20:
                batch_betas.append(
                    fixed_weight_beta(y[keep], x[keep], robust_weights[keep])
                )
        batch_betas = np.asarray(batch_betas, dtype=float)
        valid_batch = batch_betas[np.isfinite(batch_betas)]
        batch_feasible = len(valid_batch) >= 3
        batch_stable = (
            bool(np.all(np.sign(valid_batch) == np.sign(beta)))
            if batch_feasible
            else np.nan
        )

        stability_rows.append(
            {
                "state": r.state,
                "family": r.family,
                "cancer_feature": r.cancer_feature,
                "n_patients": r.n_patients,
                "full_adjusted_beta": beta,
                "bootstrap_iterations": 1000,
                "bootstrap_direction_stability": boot_stability,
                "bootstrap_CI_low": boot_low,
                "bootstrap_CI_high": boot_high,
                "LOPO_direction_stable": lopo_stable,
                "LOPO_beta_min": float(np.nanmin(valid_lopo)) if len(valid_lopo) else np.nan,
                "LOPO_beta_max": float(np.nanmax(valid_lopo)) if len(valid_lopo) else np.nan,
                "max_single_patient_abs_beta_change": max_change,
                "max_single_patient_influence_ratio": influence_ratio,
                "most_influential_patient": influential_patient,
                "no_single_patient_dominance": no_single_dominance,
                "leave_batch_out_feasible": batch_feasible,
                "leave_batch_out_direction_stable": batch_stable,
                "n_leave_batch_out_estimates": len(valid_batch),
            }
        )
    stability = pd.DataFrame(stability_rows)
    gate = assoc.merge(
        stability[
            [
                "state",
                "cancer_feature",
                "bootstrap_direction_stability",
                "LOPO_direction_stable",
                "no_single_patient_dominance",
            ]
        ],
        on=["state", "cancer_feature"],
        how="left",
        validate="one_to_one",
    )
    gate["adequate_coverage"] = gate["n_patients"] >= 20
    gate["primary_state"] = gate["family"].isin(["primary_myeloid", "primary_stromal"])
    gate["primary_edge_gate_pass"] = (
        gate["primary_state"]
        & (gate["adjusted_BH_FDR_within_family"] < 0.05)
        & gate["raw_adjusted_direction_concordant"]
        & (gate["bootstrap_direction_stability"] >= 0.80)
        & gate["LOPO_direction_stable"].fillna(False)
        & gate["adequate_coverage"]
        & gate["no_single_patient_dominance"].fillna(False)
    )
    assoc = assoc.merge(
        gate[
            [
                "state",
                "cancer_feature",
                "primary_edge_gate_pass",
            ]
        ],
        on=["state", "cancer_feature"],
        how="left",
        validate="one_to_one",
    )
    return assoc, stability


def plot_ed6(lineage: pd.DataFrame) -> None:
    landscape = pd.read_csv(LANDSCAPE_SOURCE, sep="\t")
    trends = landscape.loc[landscape["record_type"].eq("across_patient_summary")].copy()
    if len(trends) != 60:
        raise RuntimeError(f"Expected 60 across-patient trend rows, found {len(trends)}")
    de = pd.read_csv(DESEQ_RESULT, sep="\t")
    lin = lineage.loc[
        lineage["record_level"].eq("summary")
        & lineage["analysis"].eq("six_core_High_minus_Low")
        & lineage["version"].eq("core_overlap_removed")
    ].copy()

    fig = plt.figure(figsize=(180 / 25.4, 175 / 25.4))
    outer = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.04, 0.96],
        width_ratios=[1.05, 0.95],
        left=0.07,
        right=0.97,
        bottom=0.075,
        top=0.875,
        wspace=0.22,
        hspace=0.58,
    )
    top = outer[0, :].subgridspec(2, 3, wspace=0.25, hspace=0.38)
    colours = {
        "Stress-adaptive": "#5A7897",
        "Adhesion-remodeling": "#B76850",
    }
    for i, module in enumerate(MODULES):
        ax = fig.add_subplot(top[i // 3, i % 3])
        g = trends.loc[trends["module"].eq(module)].sort_values("joint_decile")
        branch = str(g["branch"].iloc[0])
        ax.fill_between(
            g["joint_decile"],
            g["q1"],
            g["q3"],
            color=colours[branch],
            alpha=0.16,
            lw=0,
        )
        ax.plot(
            g["joint_decile"],
            g["median_score_after"],
            color=colours[branch],
            lw=1.25,
        )
        ax.scatter(
            g["joint_decile"],
            g["median_score_after"],
            color=colours[branch],
            s=7,
            zorder=3,
        )
        ax.axhline(0, color="#C2C6C9", lw=0.4)
        ax.set_title(MODULE_DISPLAY[module], loc="left", fontweight="bold", pad=2)
        ax.set_xticks([1, 5, 10])
        if i // 3 == 1:
            ax.set_xlabel("Joint decile")
        if i % 3 == 0:
            ax.set_ylabel("Cohort-z score")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.spines["left"].set_linewidth(0.45)
        ax.spines["bottom"].set_linewidth(0.45)
    fig.text(0.015, 0.925, "a", fontsize=10, fontweight="bold", va="top")
    fig.text(
        0.07,
        0.925,
        "Original program trends moved from Figure 3",
        fontsize=8.2,
        fontweight="bold",
        va="top",
    )

    ax_ma = fig.add_subplot(outer[1, 0])
    sig = de["BH_FDR"].lt(0.05)
    ax_ma.scatter(
        np.log10(de.loc[~sig, "baseMean"] + 1),
        de.loc[~sig, "log2FC"],
        s=2.0,
        c="#B8BEC3",
        alpha=0.35,
        rasterized=True,
        linewidths=0,
    )
    ax_ma.scatter(
        np.log10(de.loc[sig, "baseMean"] + 1),
        de.loc[sig, "log2FC"],
        s=2.3,
        c="#B53A43",
        alpha=0.55,
        rasterized=True,
        linewidths=0,
    )
    ax_ma.axhline(0, color="#777777", lw=0.5)
    for gene in ("MCL1", "PPP1R15A", "ATF3"):
        r = de.loc[de["gene"].eq(gene)]
        if len(r):
            x = float(np.log10(r["baseMean"].iloc[0] + 1))
            y = float(r["log2FC"].iloc[0])
            ax_ma.scatter([x], [y], s=18, c="#7A1621", edgecolors="white", linewidths=0.4)
            ax_ma.annotate(
                gene,
                (x, y),
                xytext=(3, 3),
                textcoords="offset points",
                fontsize=6.0,
                fontstyle="italic",
            )
    ax_ma.set_title(
        "b  Locked DESeq2 support",
        loc="left",
        fontweight="bold",
        pad=17,
    )
    ax_ma.text(
        0.0,
        1.015,
        "score genes excluded before testing",
        transform=ax_ma.transAxes,
        fontsize=5.6,
        color="#555555",
        va="bottom",
    )
    ax_ma.set_xlabel(r"$\log_{10}(\mathrm{baseMean}+1)$")
    ax_ma.set_ylabel(r"$\log_{2}$ fold change, High vs Low")
    for side in ("top", "right"):
        ax_ma.spines[side].set_visible(False)

    ax_lin = fig.add_subplot(outer[1, 1])
    cols = [
        ("Wu2021", "before_lineage_adjustment", "Wu before"),
        ("Wu2021", "after_lineage_adjustment", "Wu adjusted"),
        ("Yan2026", "before_lineage_adjustment", "Yan before"),
        ("Yan2026", "after_lineage_adjustment", "Yan adjusted"),
    ]
    matrix = np.full((6, 4), np.nan)
    npos = np.full((6, 4), np.nan)
    for i, module in enumerate(MODULES):
        for j, (cohort, adj, _) in enumerate(cols):
            r = lin.loc[
                lin["cohort"].eq(cohort)
                & lin["module"].eq(module)
                & lin["adjustment"].eq(adj)
            ]
            if len(r):
                matrix[i, j] = float(r["effect"].iloc[0])
                npos[i, j] = float(r["n_positive"].iloc[0])
    lim = max(0.8, float(np.nanquantile(np.abs(matrix), 0.95)))
    cmap = LinearSegmentedColormap.from_list(
        "effect", ["#235789", "#E6EDF1", "#F6F4EF", "#F4C3A8", "#B53A43"]
    )
    im = ax_lin.imshow(
        matrix,
        aspect="auto",
        cmap=cmap,
        norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim),
    )
    ax_lin.set_yticks(range(6))
    ax_lin.set_yticklabels([MODULE_DISPLAY[m] for m in MODULES])
    ax_lin.set_xticks(range(4))
    ax_lin.set_xticklabels([x[2] for x in cols], rotation=35, ha="right")
    for i in range(6):
        for j in range(4):
            if np.isfinite(npos[i, j]):
                rr = lin.loc[
                    lin["cohort"].eq(cols[j][0])
                    & lin["module"].eq(MODULES[i])
                    & lin["adjustment"].eq(cols[j][1])
                ]
                denominator = int(rr["n_patients"].iloc[0]) if len(rr) else (8 if j < 2 else 78)
                ax_lin.text(
                    j,
                    i,
                    f"{int(npos[i,j])}/{denominator}",
                    ha="center",
                    va="center",
                    fontsize=5.3,
                    color="black" if abs(matrix[i, j]) < lim * 0.65 else "white",
                )
    ax_lin.set_title(
        "c  Six-core effects before/after lineage adjustment",
        loc="left",
        fontweight="bold",
        pad=17,
    )
    ax_lin.text(
        0,
        1.015,
        "tiles: median patient High–Low effect; labels: positive patients",
        transform=ax_lin.transAxes,
        fontsize=5.4,
        color="#555555",
        va="bottom",
    )
    ax_lin.tick_params(length=0)
    for side in ax_lin.spines.values():
        side.set_visible(False)
    cbar = fig.colorbar(im, ax=ax_lin, fraction=0.045, pad=0.04)
    cbar.set_label("Median effect (cohort SD)", fontsize=5.5)
    cbar.ax.tick_params(labelsize=5.2, length=2)

    fig.suptitle(
        "Extended Data 6 | Program continuity, locked gene support and epithelial-lineage audit",
        x=0.07,
        y=0.982,
        ha="left",
        fontsize=8.8,
        fontweight="bold",
    )
    fig.savefig(FIGURES / "ExtendedData6_storyboard.pdf")
    fig.savefig(FIGURES / "ExtendedData6_storyboard.png", dpi=400)
    plt.close(fig)

    write_tsv(trends, PLOTDATA / "ED6_panelA_original_program_trends.tsv")
    write_tsv(
        lin,
        PLOTDATA / "ED6_panelC_lineage_adjustment_summary.tsv",
    )


def plot_conditional_multicellular(
    assoc: pd.DataFrame,
    stability: pd.DataFrame,
    gate_edges: pd.DataFrame,
    prop: pd.DataFrame,
    features: pd.DataFrame,
) -> None:
    all_states = [
        s
        for s, _, _, family in STATE_MAP
        if family in ("primary_myeloid", "primary_stromal")
        and s in set(assoc["state"])
    ]
    all_features = list(FEATURE_MAP)
    feature_short = [
        "Joint",
        "High fraction",
        "UPR",
        "TNF-NFkB",
        "Hypoxia",
        "Survival",
        "Adhesion",
        "Wound",
    ]
    raw = assoc.pivot(
        index="state", columns="cancer_feature", values="raw_spearman_rho"
    ).reindex(index=all_states, columns=all_features)
    adj = assoc.pivot(
        index="state", columns="cancer_feature", values="adjusted_robust_beta"
    ).reindex(index=all_states, columns=all_features)
    stab = stability.pivot(
        index="state", columns="cancer_feature", values="bootstrap_direction_stability"
    ).reindex(index=all_states, columns=all_features)
    gated = {
        (r.state, r.cancer_feature) for r in gate_edges.itertuples(index=False)
    }
    cmap = LinearSegmentedColormap.from_list(
        "assoc", ["#235789", "#DDE7EC", "#F7F5F0", "#F2BA9A", "#B53A43"]
    )

    pcomp = add_compositional_coordinates(prop).merge(features, on="patient", how="inner")
    selected_states = list(dict.fromkeys(gate_edges["state"].astype(str)))
    order = features.sort_values("joint_axis_median")["patient"].astype(str).tolist()
    landscape_rows = [
        ("Joint axis median", "feature", "joint_axis_median"),
        ("High-state fraction", "feature", "high_state_fraction"),
        ("UPR", "feature", "UPR"),
        ("TNFα–NF-κB", "feature", "TNFA NFKB"),
    ] + [(s, "state", s) for s in selected_states]
    landscape = np.full((len(landscape_rows), len(order)), np.nan)
    landscape_long = []
    feature_index = features.set_index("patient")
    for i, (label, kind, key) in enumerate(landscape_rows):
        if kind == "feature":
            vals = np.array(
                [
                    feature_index.loc[p, key] if p in feature_index.index else np.nan
                    for p in order
                ],
                dtype=float,
            )
        else:
            state_index = (
                pcomp.loc[
                    pcomp["state"].eq(key) & pcomp["testability_status"].eq("testable")
                ]
                .set_index("patient")["logit_proportion"]
            )
            vals = np.array(
                [state_index.get(p, np.nan) for p in order],
                dtype=float,
            )
        landscape[i] = zscore(vals)
        for patient, value, zval in zip(order, vals, landscape[i]):
            landscape_long.append(
                {
                    "patient": patient,
                    "patient_order": order.index(patient) + 1,
                    "row": label,
                    "row_type": kind,
                    "raw_value": value,
                    "row_z": zval,
                }
            )
    write_tsv(
        pd.DataFrame(landscape_long),
        PLOTDATA / "ED_multicellular_patient_landscape.tsv",
    )
    write_tsv(gate_edges, PLOTDATA / "ED_multicellular_network_edges.tsv")

    fig = plt.figure(figsize=(180 / 25.4, 185 / 25.4))
    gs = fig.add_gridspec(
        3,
        4,
        height_ratios=[1.00, 0.58, 0.92],
        left=0.12,
        right=0.975,
        bottom=0.13,
        top=0.91,
        hspace=0.68,
        wspace=0.48,
    )
    ax_raw = fig.add_subplot(gs[0, 0])
    ax_adj = fig.add_subplot(gs[0, 1])
    ax_net = fig.add_subplot(gs[0, 2:4])
    ax_land = fig.add_subplot(gs[1, :])
    forest_gs = gs[2, 0:2].subgridspec(1, 2, wspace=0.13)
    ax_fr = fig.add_subplot(forest_gs[0, 0])
    ax_fa = fig.add_subplot(forest_gs[0, 1], sharey=ax_fr)
    ax_stab = fig.add_subplot(gs[2, 2:4])

    for ax, mat, title, lim, show_y in [
        (ax_raw, raw, "a  Raw patient Spearman", 0.60, True),
        (ax_adj, adj, "Adjusted robust β", 1.00, False),
    ]:
        im = ax.imshow(
            mat.to_numpy(float),
            aspect="auto",
            cmap=cmap,
            norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim),
        )
        ax.set_title(title, loc="left", fontweight="bold", pad=3)
        ax.set_xticks(range(len(all_features)))
        ax.set_xticklabels(feature_short, rotation=45, ha="right")
        ax.set_yticks(range(len(all_states)))
        ax.set_yticklabels(all_states if show_y else [])
        ax.tick_params(length=0)
        for i, state in enumerate(all_states):
            for j, feature in enumerate(all_features):
                if (state, feature) in gated:
                    ax.add_patch(
                        plt.Rectangle(
                            (j - 0.48, i - 0.48),
                            0.96,
                            0.96,
                            fill=False,
                            edgecolor="#2B2B2B",
                            linewidth=0.9,
                        )
                    )
        cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
        cb.ax.tick_params(labelsize=5, length=2)

    ax_net.set_title("b  Gate-passing patient-level network", loc="left", fontweight="bold")
    ax_net.set_xlim(0, 1)
    ax_net.set_ylim(0, 1)
    ax_net.axis("off")
    net_features = list(dict.fromkeys(gate_edges["cancer_feature"].astype(str)))
    net_states = list(dict.fromkeys(gate_edges["state"].astype(str)))
    network_feature_label = {
        "Joint axis median": "Joint axis",
        "High-state fraction": "High fraction",
        "TNFα–NF-κB": "TNF-NFkB",
    }
    fy = np.linspace(0.82, 0.18, len(net_features))
    sy = np.linspace(0.68, 0.32, len(net_states))
    fpos = {f: (0.20, y) for f, y in zip(net_features, fy)}
    spos = {s: (0.88, y) for s, y in zip(net_states, sy)}
    for r in gate_edges.itertuples(index=False):
        x1, y1 = fpos[r.cancer_feature]
        x2, y2 = spos[r.state]
        colour = "#B53A43" if r.adjusted_robust_beta > 0 else "#235789"
        ax_net.plot(
            [x1, x2],
            [y1, y2],
            color=colour,
            lw=0.8 + 2.2 * abs(r.adjusted_robust_beta),
            alpha=0.72,
            zorder=1,
        )
    for label, (x, y) in fpos.items():
        ax_net.scatter([x], [y], s=90, c="#E7DDD3", edgecolors="#6E6257", linewidths=0.6, zorder=2)
        ax_net.text(
            x - 0.035,
            y,
            network_feature_label.get(label, label),
            ha="right",
            va="center",
            fontsize=5.8,
        )
    for label, (x, y) in spos.items():
        ax_net.scatter([x], [y], s=95, c="#D8E5EC", edgecolors="#47677F", linewidths=0.6, zorder=2)
        ax_net.text(x + 0.035, y, label, ha="left", va="center", fontsize=6.0)
    ax_net.text(
        0.50,
        0.04,
        "5 robust edges; red positive, blue negative; width = |adjusted beta|",
        ha="center",
        fontsize=5.0,
        color="#555555",
    )

    land_cmap = LinearSegmentedColormap.from_list(
        "landscape", ["#235789", "#E1E9ED", "#F6F4EF", "#F3BFA2", "#B53A43"]
    )
    im_land = ax_land.imshow(
        landscape,
        aspect="auto",
        cmap=land_cmap,
        norm=TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5),
        interpolation="nearest",
    )
    ax_land.set_title(
        "c  Patient-resolved cancer-program/TME landscape",
        loc="left",
        fontweight="bold",
        pad=3,
    )
    ax_land.set_yticks(range(len(landscape_rows)))
    ax_land.set_yticklabels([x[0] for x in landscape_rows])
    ax_land.set_xticks([0, 19, 39, 59, 77])
    ax_land.set_xticklabels([order[i] for i in [0, 19, 39, 59, 77]])
    ax_land.set_xlabel("78 patients ordered by median Joint axis")
    ax_land.tick_params(length=0)
    for side in ax_land.spines.values():
        side.set_visible(False)
    cb = fig.colorbar(im_land, ax=ax_land, fraction=0.012, pad=0.015)
    cb.set_label("row z-score", fontsize=5.2)
    cb.ax.tick_params(labelsize=5, length=2)

    forest = gate_edges.merge(
        stability,
        on=["state", "family", "cancer_feature", "n_patients"],
        how="left",
        suffixes=("", "_stab"),
    ).copy()
    forest_feature_short = {
        feature: short for feature, short in zip(all_features, feature_short)
    }
    forest_state_short = {
        "Mac-IFN": "IFN",
        "Mac-CCL": "CCL",
    }
    forest["label"] = (
        forest["state"].map(forest_state_short).fillna(forest["state"])
        + " x "
        + forest["cancer_feature"].map(forest_feature_short).fillna(forest["cancer_feature"])
    )
    forest = forest.sort_values(["state", "cancer_feature"])
    ypos = np.arange(len(forest))
    zcrit = stats.norm.ppf(0.975)
    n = forest["n_patients"].to_numpy(float)
    rho = np.clip(forest["raw_spearman_rho"].to_numpy(float), -0.999, 0.999)
    rz = np.arctanh(rho)
    rse = 1 / np.sqrt(np.maximum(n - 3, 1))
    rlo = np.tanh(rz - zcrit * rse)
    rhi = np.tanh(rz + zcrit * rse)
    ax_fr.errorbar(
        rho,
        ypos,
        xerr=np.vstack([rho - rlo, rhi - rho]),
        fmt="o",
        color="#575757",
        ecolor="#A7ACB0",
        markersize=3.2,
        elinewidth=0.8,
        capsize=1.5,
    )
    ax_fr.axvline(0, color="#9A9A9A", lw=0.55)
    ax_fr.set_yticks(ypos)
    ax_fr.set_yticklabels(forest["label"])
    ax_fr.set_xlabel("Raw Spearman rho (95% CI)")
    ax_fr.set_title("d  Raw/adjusted effect forest", loc="left", fontweight="bold", pad=3)
    ax_fa.errorbar(
        forest["adjusted_robust_beta"],
        ypos,
        xerr=np.vstack(
            [
                forest["adjusted_robust_beta"] - forest["bootstrap_CI_low"],
                forest["bootstrap_CI_high"] - forest["adjusted_robust_beta"],
            ]
        ),
        fmt="o",
        color="#7A1621",
        ecolor="#D49A9D",
        markersize=3.2,
        elinewidth=0.8,
        capsize=1.5,
    )
    ax_fa.axvline(0, color="#9A9A9A", lw=0.55)
    ax_fa.set_xlabel("Adjusted robust beta\n(1000-bootstrap 95% CI)")
    ax_fa.tick_params(labelleft=False)
    for ax in (ax_fr, ax_fa):
        ax.invert_yaxis()
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    im_stab = ax_stab.imshow(
        stab.to_numpy(float),
        aspect="auto",
        cmap="Greens",
        vmin=0.5,
        vmax=1.0,
    )
    ax_stab.set_title("e  Network stability matrix", loc="left", fontweight="bold", pad=3)
    ax_stab.set_xticks(range(len(all_features)))
    ax_stab.set_xticklabels(feature_short, rotation=45, ha="right")
    ax_stab.set_yticks(range(len(all_states)))
    ax_stab.set_yticklabels(all_states)
    ax_stab.tick_params(length=0)
    for i, state in enumerate(all_states):
        for j, feature in enumerate(all_features):
            if (state, feature) in gated:
                ax_stab.text(j, i, "●", ha="center", va="center", fontsize=5.5, color="#151515")
    cb = fig.colorbar(im_stab, ax=ax_stab, fraction=0.025, pad=0.02)
    cb.set_label("bootstrap direction stability", fontsize=5.2)
    cb.ax.tick_params(labelsize=5, length=2)

    fig.suptitle(
        "Extended Data | Gated multicellular context",
        x=0.075,
        y=0.985,
        ha="left",
        fontsize=9.0,
        fontweight="bold",
    )
    fig.savefig(FIGURES / "ExtendedData_multicellular_context_storyboard.pdf")
    fig.savefig(FIGURES / "ExtendedData_multicellular_context_storyboard.png", dpi=400)
    plt.close(fig)

    # Non-destructive Figure 5 proposal using frozen Figure 5 evidence plus gated TME panels.
    f5 = pd.read_csv(FIG5_FROZEN_PLOTDATA, sep="\t")
    forest78 = f5.loc[f5["v4_panel"].eq("5B")].sort_values("rho").copy()
    modules = f5.loc[
        f5["v4_panel"].eq("5D") & f5["version"].eq("technical_adjusted")
    ].sort_values("module_order")
    genes = f5.loc[
        f5["v4_panel"].eq("5E") & f5["model"].eq("technical_adjusted_continuous")
    ].copy()
    fig = plt.figure(figsize=(180 / 25.4, 145 / 25.4))
    gs = fig.add_gridspec(
        2,
        4,
        width_ratios=[0.70, 1.15, 1.05, 1.10],
        height_ratios=[1.0, 1.0],
        left=0.105,
        right=0.96,
        bottom=0.09,
        top=0.91,
        wspace=0.44,
        hspace=0.42,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1:3])
    c_gs = gs[1, 0:2].subgridspec(1, 2, width_ratios=[0.70, 0.30], wspace=0.34)
    ax_c = fig.add_subplot(c_gs[0, 0])
    ax_cg = fig.add_subplot(c_gs[0, 1])
    ax_d = fig.add_subplot(gs[0, 3])
    ax_e = fig.add_subplot(gs[1, 2:4])

    ax_a.axis("off")
    ax_a.set_title("a  Yan cohort", loc="left", fontweight="bold")
    cohort_lines = [
        ("All patients", 101),
        ("Cancer cells present", 97),
        (">=50 cancer cells", 78),
        ("Cancer cells", 49275),
    ]
    for i, (label, value) in enumerate(cohort_lines):
        y = 0.80 - i * 0.19
        ax_a.text(0.02, y, f"{value:,}", fontsize=11 if i in (2, 3) else 9, fontweight="bold")
        ax_a.text(0.02, y - 0.08, label, fontsize=5.7, color="#555555")

    yy = np.arange(len(forest78))
    ax_b.hlines(
        yy,
        forest78["bootstrap_ci_low"],
        forest78["bootstrap_ci_high"],
        color="#C5C9CC",
        lw=0.45,
    )
    ax_b.scatter(forest78["rho"], yy, s=5.5, c="#486D8C", edgecolors="none")
    ax_b.axvline(0, color="#8E8E8E", lw=0.55)
    ax_b.set_ylim(-1, len(forest78))
    ax_b.set_yticks([])
    ax_b.set_xlabel("Within-patient YAP–Stem Spearman ρ")
    ax_b.set_title("b  Locked 78-patient forest", loc="left", fontweight="bold")
    ax_b.text(
        0.99,
        0.03,
        "74/78 positive",
        transform=ax_b.transAxes,
        ha="right",
        fontsize=6.2,
        fontweight="bold",
    )
    for side in ("top", "right", "left"):
        ax_b.spines[side].set_visible(False)

    ymod = np.arange(len(modules))
    ax_c.scatter(
        modules["median_delta_z"],
        ymod,
        s=28,
        c=["#5A7897" if b == "Stress-adaptive" else "#B76850" for b in modules["branch"]],
        zorder=3,
    )
    for i, r in enumerate(modules.itertuples(index=False)):
        ax_c.text(
            r.median_delta_z + 0.03,
            i,
            f"{int(r.n_positive_harmonized_z)}/{int(r.n_patients_testable)}",
            va="center",
            fontsize=5.4,
        )
    ax_c.set_yticks(ymod)
    module_short = {
        "UPR": "UPR",
        "TNFA NFKB": "TNF-NFkB",
        "Hypoxia": "Hypoxia",
        "Survival Stress": "Survival",
        "Adhesion Remodeling": "Adhesion",
        "Wound Healing": "Wound",
    }
    ax_c.set_yticklabels([module_short.get(m, m) for m in modules["module"]])
    ax_c.invert_yaxis()
    ax_c.set_xlabel("Technical-adjusted High–Low effect (cohort SD)")
    ax_c.set_title("c  Six-core replication", loc="left", fontweight="bold")
    for side in ("top", "right"):
        ax_c.spines[side].set_visible(False)

    genes = genes.sort_values("median_high_minus_low", ascending=False)
    ygene = np.arange(len(genes))
    ax_cg.scatter(
        genes["median_high_minus_low"],
        ygene,
        s=24,
        c="#7A1621",
        edgecolors="white",
        linewidths=0.4,
        zorder=3,
    )
    ax_cg.axvline(0, color="#9A9A9A", lw=0.5)
    ax_cg.set_yticks(ygene)
    ax_cg.set_yticklabels(genes["gene"])
    for tick in ax_cg.get_yticklabels():
        tick.set_fontstyle("italic")
    ax_cg.invert_yaxis()
    ax_cg.set_xlabel("Adjusted effect")
    ax_cg.set_title("Core genes", loc="left", fontweight="bold")
    for side in ("top", "right"):
        ax_cg.spines[side].set_visible(False)

    ax_d.set_title("d  Selected TME network", loc="left", fontweight="bold")
    ax_d.set_xlim(0, 1)
    ax_d.set_ylim(0, 1)
    ax_d.axis("off")
    for r in gate_edges.itertuples(index=False):
        f_idx = net_features.index(r.cancer_feature)
        s_idx = net_states.index(r.state)
        y1 = 0.83 - f_idx * (0.60 / max(len(net_features) - 1, 1))
        y2 = 0.65 - s_idx * (0.30 / max(len(net_states) - 1, 1))
        colour = "#B53A43" if r.adjusted_robust_beta > 0 else "#235789"
        ax_d.plot([0.24, 0.75], [y1, y2], color=colour, lw=0.8 + 2 * abs(r.adjusted_robust_beta), alpha=0.7)
    for i, f in enumerate(net_features):
        y = 0.83 - i * (0.60 / max(len(net_features) - 1, 1))
        ax_d.scatter([0.24], [y], s=42, c="#E7DDD3", edgecolors="#6E6257", linewidths=0.5)
        ax_d.text(
            0.17,
            y,
            network_feature_label.get(f, f),
            ha="right",
            va="center",
            fontsize=4.9,
        )
    for i, s in enumerate(net_states):
        y = 0.65 - i * (0.30 / max(len(net_states) - 1, 1))
        ax_d.scatter([0.75], [y], s=46, c="#D8E5EC", edgecolors="#47677F", linewidths=0.5)
        ax_d.text(0.82, y, s, ha="left", va="center", fontsize=5.2)

    im = ax_e.imshow(
        landscape,
        aspect="auto",
        cmap=land_cmap,
        norm=TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5),
        interpolation="nearest",
    )
    ax_e.set_title("e  Patient-resolved Joint/TME landscape", loc="left", fontweight="bold")
    ax_e.set_yticks(range(len(landscape_rows)))
    ax_e.set_yticklabels([x[0] for x in landscape_rows])
    ax_e.set_xticks([0, 25, 51, 77])
    ax_e.set_xticklabels([order[i] for i in [0, 25, 51, 77]])
    ax_e.set_xlabel("78 patients ordered by median Joint axis")
    ax_e.tick_params(length=0)
    for side in ax_e.spines.values():
        side.set_visible(False)
    cb = fig.colorbar(im, ax=ax_e, fraction=0.020, pad=0.012)
    cb.ax.tick_params(labelsize=5, length=2)
    cb.set_label("row z", fontsize=5.2)

    fig.suptitle(
        "Figure 5 TME integration option | locked validation plus gated multicellular context",
        x=0.105,
        y=0.982,
        ha="left",
        fontsize=8.7,
        fontweight="bold",
    )
    fig.text(
        0.98,
        0.02,
        "Proposal only • current Figure 5 not overwritten • complete M01–M13 remains in ED8",
        ha="right",
        fontsize=5.3,
        color="#555555",
    )
    fig.savefig(FIGURES / "Figure5_TME_integration_option.pdf")
    fig.savefig(FIGURES / "Figure5_TME_integration_option.png", dpi=400)
    plt.close(fig)


def main() -> None:
    lineage_result_path = TABLES / "lineage_confounding_results.tsv"
    lineage_coverage_path = TABLES / "lineage_gene_coverage_Wu_Yan.tsv"
    yan = pd.read_csv(YAN_SCORES, sep="\t")
    features = build_cancer_features(yan)
    eligible_yan_patients = set(features["patient"].astype(str))
    lineage_cache_valid = False
    if lineage_coverage_path.exists():
        coverage_existing = pd.read_csv(lineage_coverage_path, sep="\t")
        yan_coverage = coverage_existing.loc[coverage_existing["cohort"].eq("Yan2026")]
        lineage_cache_valid = (
            len(yan_coverage) == 8
            and yan_coverage["n_present"].ge(2).all()
            and lineage_result_path.exists()
        )
    if lineage_cache_valid:
        lineage = pd.read_csv(lineage_result_path, sep="\t")
    else:
        yan_with_lineage, coverage_yan, manifest = extract_yan_lineage()
        coverage_wu = pd.read_csv(TABLES / "lineage_gene_coverage_wu.tsv", sep="\t")
        write_tsv(
            pd.concat([coverage_wu, coverage_yan], ignore_index=True),
            lineage_coverage_path,
        )
        wu = prepare_wu()
        lineage = pd.concat(
            [
                lineage_analysis_one("Wu2021", wu),
                lineage_analysis_one("Yan2026", yan_with_lineage),
            ],
            ignore_index=True,
            sort=False,
        )
        lineage["interpretation_limit"] = (
            "association/confounding audit only; no cell-of-origin or lineage-tracing claim"
        )
    lineage = restrict_and_resummarize_lineage(lineage, eligible_yan_patients)
    write_tsv(lineage, lineage_result_path)
    a = ad.read_h5ad(H5AD, backed="r")
    obs = a.obs[
        ["donor_id", "orig.ident", "author_cell_type", "cell_state", "pCR_status"]
    ].copy()
    a.file.close()
    obs["donor_id"] = obs["donor_id"].astype(str)
    obs["orig.ident"] = obs["orig.ident"].astype(str)
    obs["author_cell_type"] = obs["author_cell_type"].astype(str)
    obs["cell_state"] = obs["cell_state"].astype(str)
    prop_path = TABLES / "Yan_TME_state_proportions.tsv"
    feasibility_path = TABLES / "Yan_TME_state_testability.tsv"
    feature_path = TABLES / "Yan_cancer_features_for_TME_gate.tsv"
    if prop_path.exists() and feasibility_path.exists() and feature_path.exists():
        prop = pd.read_csv(prop_path, sep="\t")
        feasibility = pd.read_csv(feasibility_path, sep="\t")
        features = pd.read_csv(feature_path, sep="\t")
    else:
        prop, feasibility = build_tme_proportions(
            obs, set(features["patient"].astype(str))
        )
        write_tsv(prop, prop_path)
        write_tsv(feasibility, feasibility_path)
        write_tsv(features, feature_path)
    assoc, stability = robust_associations(prop, feasibility, features)
    write_tsv(assoc, TABLES / "Yan_TME_raw_adjusted_associations.tsv")
    write_tsv(stability, TABLES / "Yan_TME_edge_stability.tsv")

    gate_edges = assoc.loc[assoc["primary_edge_gate_pass"].fillna(False)].copy()
    n_states_passing = gate_edges["state"].nunique()
    stage6_pass = n_states_passing >= 2
    gate_summary = pd.DataFrame(
        [
            {
                "metric": "primary_edges_passing_all_criteria",
                "value": len(gate_edges),
                "threshold": "descriptive",
                "status": "PASS" if len(gate_edges) else "NONE",
            },
            {
                "metric": "distinct_primary_myeloid_or_stromal_states_passing",
                "value": n_states_passing,
                "threshold": ">=2",
                "status": "PASS" if stage6_pass else "FAIL",
            },
            {
                "metric": "Stage6_conditional_figure_gate",
                "value": "GENERATE" if stage6_pass else "DO_NOT_GENERATE",
                "threshold": ">=2 primary states",
                "status": "PASS" if stage6_pass else "FAIL",
            },
        ]
    )
    write_tsv(gate_summary, TABLES / "Yan_TME_primary_gate_summary.tsv")
    write_tsv(gate_edges, TABLES / "Yan_TME_gate_passing_edges.tsv")

    plot_ed6(lineage)
    if stage6_pass:
        plot_conditional_multicellular(
            assoc, stability, gate_edges, prop, features
        )

    primary_feas = feasibility.loc[
        feasibility["family"].isin(["primary_myeloid", "primary_stromal"])
    ]
    n_primary_testable = int(primary_feas["state_testability_pass"].sum())
    best = assoc.sort_values("adjusted_BH_FDR_within_family").head(10)
    best_lines = "\n".join(
        f"- {r.state} × {r.cancer_feature}: raw ρ={r.raw_spearman_rho:.3f}, "
        f"adjusted β={r.adjusted_robust_beta:.3f}, adjusted FDR={r.adjusted_BH_FDR_within_family:.3g}, "
        f"gate={'PASS' if r.primary_edge_gate_pass else 'FAIL'}"
        for r in best.itertuples(index=False)
    )
    report = f"""# Yan TME multicellular context feasibility and gate

## Frozen inputs and scope

- Author-provided `author_cell_type` and `cell_state` labels were used without reclustering or renaming.
- Cancer features were frozen Yan scores from the 78 patients with at least 50 author-labeled cancer cells.
- State abundance was defined as state cells divided by the same parent-compartment cells per patient.
- A parent compartment required at least 30 cells; a missing state was retained as `NT`, never converted to zero.

## Testability

- {n_primary_testable}/{len(primary_feas)} prespecified primary myeloid/stromal states passed both cohort-total >=100 cells and >=20 testable cancer-feature patients.
- All source labels and the display-name mapping are retained in `Yan_TME_state_testability.tsv`.

## Association and stability gate

Primary edges required all of the following: adjusted within-family FDR <0.05; raw/adjusted direction concordance; 1,000-bootstrap direction stability >=0.80; stable leave-one-patient-out direction; adequate coverage; and no single-patient dominance. Single-patient dominance was prespecified as a maximum leave-one-out beta change greater than the magnitude of the full beta or any direction reversal.

- Passing primary edges: {len(gate_edges)}
- Distinct passing primary states: {n_states_passing}
- Stage 6 conditional figure gate: **{'PASS' if stage6_pass else 'FAIL'}**

## Strongest adjusted associations (reported regardless of direction or significance)

{best_lines}

## Decision

{"At least two primary states passed. The gated Extended Data multicellular context storyboard and a non-destructive Figure 5 option were generated." if stage6_pass else "Fewer than two primary myeloid/stromal states passed the complete stability gate. This is a null/boundary multicellular-context result; no decorative network and no Figure 5 integration option were generated."}

## Limitations

- The analysis is patient-level association, not ligand-receptor causality, spatial co-localization, cell-cell communication or lineage tracing.
- Parent-normalized proportions remain compositional. Logit and CLR coordinates are reported as sensitivities, not as independent biological endpoints.
- Batch leave-out is reported where at least three valid leave-batch estimates remained; it was not substituted for the prespecified primary gate.
- No pCR/RD network mining, NMF, ecotype, spatial niche, classifier, Cox or survival analysis was performed.
"""
    (REPORTS / "Yan_TME_feasibility_report.md").write_text(report, encoding="utf-8")

    provenance_files = [
        H5AD,
        YAN_SCORES,
        WU_LINEAGE,
        WU_MODULES,
        LINEAGE_MANIFEST,
        LANDSCAPE_SOURCE,
        DESEQ_RESULT,
        FIG5_FROZEN_PLOTDATA,
    ]
    provenance = pd.DataFrame(
        [
            {
                "source_file": str(p),
                "sha256": sha256(p),
                "size_bytes": p.stat().st_size,
                "role": {
                    str(H5AD): "Yan author labels and lineage-gene subset",
                    str(YAN_SCORES): "frozen Yan cancer features",
                    str(WU_LINEAGE): "Wu frozen lineage scores",
                    str(WU_MODULES): "Wu six-core module scores",
                    str(LINEAGE_MANIFEST): "predefined lineage signatures",
                    str(LANDSCAPE_SOURCE): "approved 480-row program landscape",
                    str(DESEQ_RESULT): "locked score-gene-excluded DESeq2 support",
                    str(FIG5_FROZEN_PLOTDATA): "frozen Figure 5 retained evidence for non-destructive option",
                }[str(p)],
            }
            for p in provenance_files
        ]
    )
    write_tsv(provenance, TABLES / "provenance_hashes.tsv")

    lineage_summary_n = int(
        (
            lineage["record_level"].eq("summary")
            & lineage["analysis"].eq("six_core_High_minus_Low")
        ).sum()
    )
    qc = pd.DataFrame(
        [
            {
                "check": "Yan_author_cells",
                "observed": len(obs),
                "expected": 427823,
                "status": "PASS" if len(obs) == 427823 else "FAIL",
            },
            {
                "check": "Yan_cancer_cells",
                "observed": len(yan),
                "expected": 49275,
                "status": "PASS" if len(yan) == 49275 else "FAIL",
            },
            {
                "check": "Yan_cancer_feature_patients",
                "observed": len(features),
                "expected": 78,
                "status": "PASS" if len(features) == 78 else "FAIL",
            },
            {
                "check": "TME_prespecified_states",
                "observed": len(feasibility),
                "expected": 20,
                "status": "PASS" if len(feasibility) == 20 else "FAIL",
            },
            {
                "check": "lineage_summary_records",
                "observed": lineage_summary_n,
                "expected": 48,
                "status": "PASS" if lineage_summary_n == 48 else "FAIL",
            },
            {
                "check": "stage6_conditional_output_rule",
                "observed": "generated" if stage6_pass else "not_generated",
                "expected": "generated only when >=2 primary states pass",
                "status": "PASS",
            },
        ]
    )
    write_tsv(qc, TABLES / "numerical_QC.tsv")
    (REPORTS / "limitations.md").write_text(
        "# Limitations\n\n"
        "- Program landscapes are descriptive and introduce no new trend P values.\n"
        "- DESeq2 is a locked supportive analysis; it does not replace or redefine frozen main-figure statistics.\n"
        "- Epithelial lineage scores are association covariates, not cell-of-origin or lineage-tracing evidence.\n"
        "- Yan TME analyses use author labels and patient-level proportions only; no reclustering, CellChat, NMF, ecotype, spatial or response-network analysis was run.\n"
        "- Missing fine states are NT, not zero; this preserves the gate but limits complete compositional coverage.\n",
        encoding="utf-8",
    )
    (LOGS / "03_lineage_and_yan_tme_gate.log").write_text(
        "PASS\n"
        f"Lineage rows: {len(lineage)}\n"
        f"TME association rows: {len(assoc)}\n"
        f"Gate edges: {len(gate_edges)}\n"
        f"Distinct primary states passing: {n_states_passing}\n"
        f"Stage6: {'GENERATED' if stage6_pass else 'NOT_GENERATED'}\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
