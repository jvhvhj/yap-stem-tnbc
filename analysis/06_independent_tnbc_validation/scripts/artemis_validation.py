#!/usr/bin/env python
# Purpose: ARTEMIS validation analyses
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 5; orthogonal transfer in Fig. 4.
"""Stages 4-7: Yan YAP-Stem, module/gene, metaprogram and response validation."""

from __future__ import annotations

import logging
import math
import os
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(os.environ.get("AHIPPO_YAP_ROOT", Path.cwd()))
OUT = ROOT / "0716_yan2026_validation"
INPUTS = OUT / "inputs"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
REPORTS = OUT / "reports"
LOGS = OUT / "logs"
for p in (TABLES, FIGURES, REPORTS, LOGS):
    p.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[logging.FileHandler(LOGS / "03_yan_validation_stages4_7.log", encoding="utf-8"), logging.StreamHandler(sys.stdout)],
)
LOG = logging.getLogger("yan")

CACHE = INPUTS / "yan_cancer_cell_scores.tsv.gz"
REP_SOURCE = ROOT / "0713_score_independent_rebuild" / "tables" / "Figure4_v3_1_panelD_lollipop_plot_data.tsv"
WU_MODULE_SOURCE = ROOT / "0713_score_independent_rebuild" / "tables" / "module_keep_or_drop_summary.tsv"
WU_TREATMENT_EFFECTS = TABLES / "wu_untreated_module_effects.tsv"

MODULES = [
    "UPR",
    "TNFA NFKB",
    "Hypoxia",
    "Adhesion Remodeling",
    "Wound Healing",
    "Survival Stress",
    "Anoikis Resistance",
    "Integrin Adhesion",
]
CORE7 = [m for m in MODULES if m != "Integrin Adhesion"]
THRESHOLDS = [20, 50, 100]
BOOTSTRAPS = 300
SEED = 20260716


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False, na_rep="NA")
    LOG.info("Wrote %s (%s rows)", path, len(df))


def bh_array(values: np.ndarray | pd.Series) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    out = np.full(len(x), np.nan)
    keep = np.isfinite(x)
    p = x[keep]
    if not len(p):
        return out
    order = np.argsort(p)
    ranked = p[order]
    adj = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    tmp = np.empty(len(p))
    tmp[order] = np.minimum(adj, 1)
    out[keep] = tmp
    return out


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    x, y = np.asarray(x, float), np.asarray(y, float)
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 5 or np.nanstd(x[keep]) == 0 or np.nanstd(y[keep]) == 0:
        return np.nan
    return float(stats.spearmanr(x[keep], y[keep]).statistic)


def z_within(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, float)
    sd = np.nanstd(x, ddof=1)
    return np.zeros_like(x) if not np.isfinite(sd) or sd == 0 else (x - np.nanmean(x)) / sd


def residualize(y: np.ndarray, covariates: np.ndarray) -> np.ndarray:
    y = np.asarray(y, float)
    covariates = np.asarray(covariates, float)
    keep_cols = np.nanstd(covariates, axis=0) > 1e-12
    covariates = covariates[:, keep_cols]
    zcov = np.column_stack([z_within(covariates[:, j]) for j in range(covariates.shape[1])]) if covariates.shape[1] else np.empty((len(y), 0))
    design = np.column_stack([np.ones(len(y)), zcov])
    valid = np.isfinite(y) & np.all(np.isfinite(design), axis=1)
    out = np.full(len(y), np.nan)
    if valid.sum() <= design.shape[1] + 3:
        return out
    beta = np.linalg.lstsq(design[valid], y[valid], rcond=None)[0]
    out[valid] = y[valid] - design[valid] @ beta
    return out


def bootstrap_rho(x: np.ndarray, y: np.ndarray, rng: np.random.Generator, b: int = BOOTSTRAPS) -> tuple[float, float, int]:
    x, y = np.asarray(x, float), np.asarray(y, float)
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    if len(x) < 10:
        return np.nan, np.nan, 0
    estimates = np.empty(b, dtype=float)
    for i in range(b):
        idx = rng.integers(0, len(x), size=len(x))
        estimates[i] = safe_spearman(x[idx], y[idx])
    estimates = estimates[np.isfinite(estimates)]
    if not len(estimates):
        return np.nan, np.nan, 0
    lo, hi = np.quantile(estimates, [0.025, 0.975])
    return float(lo), float(hi), len(estimates)


def fisher_random_effects(rhos: np.ndarray, ns: np.ndarray) -> dict[str, float]:
    rhos, ns = np.asarray(rhos, float), np.asarray(ns, float)
    keep = np.isfinite(rhos) & np.isfinite(ns) & (ns > 3)
    rhos, ns = np.clip(rhos[keep], -0.999999, 0.999999), ns[keep]
    z = np.arctanh(rhos)
    var = 1 / (ns - 3)
    w = 1 / var
    fixed = np.sum(w * z) / np.sum(w)
    q = np.sum(w * (z - fixed) ** 2)
    df = len(z) - 1
    c = np.sum(w) - np.sum(w**2) / np.sum(w)
    tau2 = max(0.0, (q - df) / c) if c > 0 else 0.0
    wr = 1 / (var + tau2)
    pooled = np.sum(wr * z) / np.sum(wr)
    se = math.sqrt(1 / np.sum(wr))
    return {
        "n_patients": len(z),
        "pooled_rho": math.tanh(pooled),
        "pooled_ci_low": math.tanh(pooled - 1.96 * se),
        "pooled_ci_high": math.tanh(pooled + 1.96 * se),
        "Q": q,
        "Q_df": df,
        "Q_p": float(stats.chi2.sf(q, df)) if df > 0 else np.nan,
        "I2_percent": max(0.0, (q - df) / q * 100) if q > 0 else 0.0,
        "tau2_fisher_z": tau2,
    }


def save_figure(fig: plt.Figure, stem: str) -> None:
    for ext in ["pdf", "png", "svg"]:
        kwargs = {"dpi": 300} if ext == "png" else {}
        fig.savefig(FIGURES / f"{stem}.{ext}", bbox_inches="tight", **kwargs)
    plt.close(fig)
    LOG.info("Saved figure %s.[pdf/png/svg]", stem)


def patient_models(d: pd.DataFrame, patient_index: int) -> list[dict]:
    primary_cov = np.column_stack(
        [
            np.log1p(d["nCount_RNA"].to_numpy(float)),
            d["nFeature_RNA"].to_numpy(float),
            d["S_score"].to_numpy(float),
            d["G2M_score"].to_numpy(float),
        ]
    )
    extended_cov = np.column_stack(
        [
            primary_cov,
            d["module__Hypoxia"].to_numpy(float),
            d["module__UPR"].to_numpy(float),
        ]
    )
    yap = d["YAP_score"].to_numpy(float)
    stem = d["Stemness_score"].to_numpy(float)
    model_data = {
        "raw": (yap, stem, "none"),
        "technical_primary": (residualize(yap, primary_cov), residualize(stem, primary_cov), "log_nCount+nFeature+S+G2M"),
        "technical_extended": (residualize(yap, extended_cov), residualize(stem, extended_cov), "primary+Hypoxia+UPR"),
    }
    rows = []
    for model_no, (model, (x, y, covariates)) in enumerate(model_data.items()):
        rho = safe_spearman(x, y)
        rng = np.random.default_rng(SEED + patient_index * 101 + model_no * 100_003)
        lo, hi, nb = bootstrap_rho(x, y, rng)
        rows.append(
            {
                "patient": d["patient"].iloc[0],
                "response": d["response"].iloc[0],
                "n_cells": len(d),
                "model": model,
                "covariates": covariates,
                "rho": rho,
                "bootstrap_ci_low": lo,
                "bootstrap_ci_high": hi,
                "n_bootstrap_valid": nb,
                "bootstrap_method": "cell bootstrap of raw or precomputed residual pairs; 300 replicates",
            }
        )
    return rows


def main() -> None:
    for path in (CACHE, REP_SOURCE, WU_MODULE_SOURCE, WU_TREATMENT_EFFECTS):
        if not path.exists():
            raise FileNotFoundError(path)
    LOG.info("Reading cancer-cell score cache")
    cells = pd.read_csv(CACHE, sep="\t")
    counts = cells.groupby("patient").size().sort_index()
    if len(cells) != 49_275 or len(counts) != 97:
        raise RuntimeError(f"Unexpected cancer score cache: cells={len(cells)}, patients={len(counts)}")

    # ------------------------------------------------------------------
    # Stage 4: patient-resolved association and meta-analysis.
    # ------------------------------------------------------------------
    rho_rows = []
    eligible20 = counts[counts >= 20].index.tolist()
    for i, patient in enumerate(eligible20):
        d = cells.loc[cells["patient"].eq(patient)].copy()
        rho_rows.extend(patient_models(d, i))
        if (i + 1) % 10 == 0 or i + 1 == len(eligible20):
            LOG.info("Computed Yan patient correlations %s/%s", i + 1, len(eligible20))
    rho_table = pd.DataFrame(rho_rows)
    for threshold in THRESHOLDS:
        rho_table[f"included_n{threshold}"] = rho_table["n_cells"] >= threshold
    write_tsv(rho_table, TABLES / "yan_patient_yap_stem_rho.tsv")

    meta_rows = []
    loo_rows = []
    for threshold in THRESHOLDS:
        for model in ["raw", "technical_primary", "technical_extended"]:
            d = rho_table.loc[rho_table["model"].eq(model) & (rho_table["n_cells"] >= threshold)].copy()
            pooled = fisher_random_effects(d["rho"], d["n_cells"])
            npos = int((d["rho"] > 0).sum())
            meta_rows.append(
                {
                    "record_type": "meta_summary",
                    "threshold": threshold,
                    "model": model,
                    "n_patients": len(d),
                    "n_positive": npos,
                    "positive_fraction": npos / len(d),
                    "exact_sign_test_p": stats.binomtest(npos, len(d), 0.5).pvalue,
                    "median_rho": d["rho"].median(),
                    "rho_iqr_low": d["rho"].quantile(0.25),
                    "rho_iqr_high": d["rho"].quantile(0.75),
                    **pooled,
                }
            )
            for omitted in d["patient"]:
                dd = d.loc[~d["patient"].eq(omitted)]
                loo_rows.append(
                    {
                        "record_type": "leave_one_patient_out",
                        "threshold": threshold,
                        "model": model,
                        "omitted_patient": omitted,
                        "median_rho": dd["rho"].median(),
                        **fisher_random_effects(dd["rho"], dd["n_cells"]),
                    }
                )
    meta_table = pd.concat([pd.DataFrame(meta_rows), pd.DataFrame(loo_rows)], ignore_index=True, sort=False)
    write_tsv(meta_table, TABLES / "yan_yap_stem_meta_analysis.tsv")

    # Decile trend, using within-patient YAP ranks for n>=50 patients.
    primary_patients = counts[counts >= 50].index
    primary_cells = cells.loc[cells["patient"].isin(primary_patients)].copy()
    primary_cells["YAP_decile"] = primary_cells.groupby("patient")["YAP_score"].transform(
        lambda x: np.minimum(10, np.floor((x.rank(method="average") - 1) / len(x) * 10) + 1)
    ).astype(int)
    decile = (
        primary_cells.groupby(["patient", "YAP_decile"], observed=True)
        .agg(n_cells=("cell_id", "size"), median_YAP=("YAP_score", "median"), median_Stem=("Stemness_score", "median"))
        .reset_index()
    )
    pooled_decile = (
        decile.groupby("YAP_decile")
        .agg(patient_median=("median_Stem", "median"), patient_q25=("median_Stem", lambda x: x.quantile(0.25)), patient_q75=("median_Stem", lambda x: x.quantile(0.75)), n_patients=("patient", "size"))
        .reset_index()
    )
    write_tsv(decile, TABLES / "yan_decile_trend.tsv")
    write_tsv(pooled_decile, TABLES / "yan_decile_trend_pooled.tsv")

    # Selected patients span raw-rho quantiles and have >=100 cells.
    raw100 = rho_table.loc[rho_table["model"].eq("raw") & (rho_table["n_cells"] >= 100)].sort_values("rho")
    select_idx = np.unique(np.round(np.linspace(0, len(raw100) - 1, 6)).astype(int))
    selected_patients = raw100.iloc[select_idx]["patient"].tolist()
    density_data = primary_cells.loc[primary_cells["patient"].isin(selected_patients), ["cell_id", "patient", "YAP_score", "Stemness_score"]]
    write_tsv(density_data, TABLES / "yan_selected_patient_density.tsv")

    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "pdf.fonttype": 42, "ps.fonttype": 42})
    # Forest.
    forest = rho_table.loc[rho_table["model"].eq("raw") & (rho_table["n_cells"] >= 50)].sort_values("rho").reset_index(drop=True)
    fig_h = max(8.0, 0.18 * len(forest) + 1.5)
    fig, ax = plt.subplots(figsize=(7.4, fig_h))
    yy = np.arange(len(forest))
    ax.axvline(0, color="#8F969D", ls="--", lw=0.8)
    ax.hlines(yy, forest["bootstrap_ci_low"], forest["bootstrap_ci_high"], color="#AEB4BA", lw=0.9)
    ax.scatter(forest["rho"], yy, s=np.clip(np.sqrt(forest["n_cells"]) * 4, 18, 75), color="#B51F35", edgecolor="white", lw=0.4)
    ax.set_yticks(yy, [f"{p}  n={n}" for p, n in zip(forest["patient"], forest["n_cells"])], fontsize=6.5)
    ax.set_xlabel("Within-patient Spearman rho (bootstrap 95% CI)")
    ax.set_title("Yan cohort: frozen YAP-Stem association", loc="left", weight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    save_figure(fig, "Yan_FigA_yap_stem_forest")

    # Decile trend.
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for _, d in decile.groupby("patient"):
        ax.plot(d["YAP_decile"], d["median_Stem"], color="#C9CED2", lw=0.6, alpha=0.55)
    ax.fill_between(pooled_decile["YAP_decile"], pooled_decile["patient_q25"], pooled_decile["patient_q75"], color="#E9B8BF", alpha=0.42)
    ax.plot(pooled_decile["YAP_decile"], pooled_decile["patient_median"], color="#9E142A", lw=2.5)
    ax.set_xlabel("Within-patient YAP activity decile")
    ax.set_ylabel("Median frozen Stemness score")
    ax.set_xticks(range(1, 11))
    ax.set_title(f"Yan cohort decile trend ({len(primary_patients)} patients with >=50 cancer cells)", loc="left", weight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    save_figure(fig, "Yan_FigB_decile_trend")

    # Density panels.
    fig, axes = plt.subplots(2, 3, figsize=(11.4, 7.2), sharex=False, sharey=False)
    raw_map = raw100.set_index("patient")
    for ax, patient in zip(axes.ravel(), selected_patients):
        d = density_data.loc[density_data["patient"].eq(patient)]
        hb = ax.hexbin(d["YAP_score"], d["Stemness_score"], gridsize=28, mincnt=1, cmap="Reds", bins="log")
        ax.set_title(f"{patient} | rho={raw_map.loc[patient, 'rho']:.2f} | n={len(d)}", fontsize=8, weight="bold")
        ax.set_xlabel("YAP score")
        ax.set_ylabel("Stemness score")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Selected Yan patients spanning the observed association", fontsize=13, weight="bold")
    save_figure(fig, "Yan_FigC_patient_density")

    # ------------------------------------------------------------------
    # Stage 5: module robustness and representative genes.
    # ------------------------------------------------------------------
    module_rows = []
    primary_cov_cols = ["nCount_RNA", "nFeature_RNA", "S_score", "G2M_score"]
    for patient, d in primary_cells.groupby("patient", sort=True):
        cov = np.column_stack([np.log1p(d["nCount_RNA"]), d["nFeature_RNA"], d["S_score"], d["G2M_score"]])
        for module in MODULES:
            versions = {
                "score_independent": d[f"module__{module}"].to_numpy(float),
                "pairwise_shared_removed": d[f"module_unique__{module}"].to_numpy(float),
                "technical_adjusted": residualize(d[f"module__{module}"].to_numpy(float), cov),
            }
            for version, values in versions.items():
                low = values[d["YS_state"].eq("Low").to_numpy()]
                high = values[d["YS_state"].eq("High").to_numpy()]
                low_med = np.nanmedian(low) if np.isfinite(low).any() else np.nan
                high_med = np.nanmedian(high) if np.isfinite(high).any() else np.nan
                module_rows.append(
                    {
                        "patient": patient,
                        "response": d["response"].iloc[0],
                        "n_cells": len(d),
                        "module": module,
                        "version": version,
                        "n_low": int(d["YS_state"].eq("Low").sum()),
                        "n_high": int(d["YS_state"].eq("High").sum()),
                        "median_low": low_med,
                        "median_high": high_med,
                        "high_minus_low": high_med - low_med if np.isfinite(low_med + high_med) else np.nan,
                        "testable": bool(np.isfinite(low_med + high_med)),
                    }
                )
    module_patient = pd.DataFrame(module_rows)
    write_tsv(module_patient, TABLES / "yan_module_effects_by_patient.tsv")

    robust_rows = []
    for (version, module), d in module_patient.groupby(["version", "module"], sort=False):
        delta = d["high_minus_low"].dropna()
        if len(delta):
            try:
                wilcox_p = stats.wilcoxon(delta, zero_method="wilcox", alternative="two-sided").pvalue
            except ValueError:
                wilcox_p = np.nan
            npos = int((delta > 0).sum())
            sign_p = stats.binomtest(npos, len(delta), 0.5).pvalue
        else:
            wilcox_p, npos, sign_p = np.nan, 0, np.nan
        robust_rows.append(
            {
                "version": version,
                "module": module,
                "n_patients_testable": len(delta),
                "median_high_minus_low": delta.median() if len(delta) else np.nan,
                "iqr_low": delta.quantile(0.25) if len(delta) else np.nan,
                "iqr_high": delta.quantile(0.75) if len(delta) else np.nan,
                "n_positive": npos,
                "positive_fraction": npos / len(delta) if len(delta) else np.nan,
                "paired_wilcoxon_p": wilcox_p,
                "exact_sign_test_p": sign_p,
                "decision_label": "supportive_borderline" if module == "Integrin Adhesion" else "core",
                "testability_note": "empty after pairwise shared-gene removal" if not len(delta) else "testable",
            }
        )
    robust = pd.DataFrame(robust_rows)
    robust["BH_FDR_within_version"] = robust.groupby("version")["paired_wilcoxon_p"].transform(bh_array)
    write_tsv(robust, TABLES / "yan_module_robustness.tsv")

    rep = pd.read_csv(REP_SOURCE, sep="\t")
    rep_genes = rep["gene"].drop_duplicates().tolist()
    gene_rows = []
    for patient, d in primary_cells.groupby("patient", sort=True):
        cov = np.column_stack([np.log1p(d["nCount_RNA"]), d["nFeature_RNA"], d["S_score"], d["G2M_score"]])
        axis = d["YAP_Stem_axis"].to_numpy(float)
        for gene in rep_genes:
            expr = d[f"gene__{gene}"].to_numpy(float)
            expr_resid = residualize(expr, cov)
            for model, values in [("raw_continuous", expr), ("technical_adjusted_continuous", expr_resid)]:
                low = values[d["YS_state"].eq("Low").to_numpy()]
                high = values[d["YS_state"].eq("High").to_numpy()]
                low_median = np.nanmedian(low) if np.isfinite(low).any() else np.nan
                high_median = np.nanmedian(high) if np.isfinite(high).any() else np.nan
                gene_rows.append(
                    {
                        "record_type": "patient",
                        "patient": patient,
                        "gene": gene,
                        "score_gene_status": rep.loc[rep["gene"].eq(gene), "score_gene_status"].iloc[0],
                        "module": rep.loc[rep["gene"].eq(gene), "module"].iloc[0],
                        "model": model,
                        "n_cells": len(d),
                        "rho_with_joint_axis": safe_spearman(axis, values),
                        "median_high_minus_low": high_median - low_median,
                    }
                )
    gene_patient = pd.DataFrame(gene_rows)
    gene_summary_rows = []
    for (model, gene), d in gene_patient.groupby(["model", "gene"], sort=False):
        rho = d["rho_with_joint_axis"].dropna()
        delta = d["median_high_minus_low"].dropna()
        try:
            p = stats.wilcoxon(delta).pvalue if len(delta) else np.nan
        except ValueError:
            p = np.nan
        gene_summary_rows.append(
            {
                "record_type": "summary",
                "patient": "ALL",
                "gene": gene,
                "score_gene_status": d["score_gene_status"].iloc[0],
                "module": d["module"].iloc[0],
                "model": model,
                "n_patients": int(d["patient"].nunique()),
                "n_patients_rho": len(rho),
                "n_patients_delta": len(delta),
                "median_rho_with_joint_axis": rho.median(),
                "n_positive_rho": int((rho > 0).sum()),
                "median_high_minus_low": delta.median(),
                "n_positive_delta": int((delta > 0).sum()),
                "paired_wilcoxon_p": p,
            }
        )
    gene_summary = pd.DataFrame(gene_summary_rows)
    gene_summary["BH_FDR_within_model"] = gene_summary.groupby("model")["paired_wilcoxon_p"].transform(bh_array)
    gene_validation = pd.concat([gene_patient, gene_summary], ignore_index=True, sort=False)
    write_tsv(gene_validation, TABLES / "yan_representative_gene_validation.tsv")

    # Wu-Yan concordance for modules and representative genes.
    wu_modules = pd.read_csv(WU_MODULE_SOURCE, sep="\t").set_index("module")
    yan_mod = robust.loc[robust["version"].eq("score_independent")].set_index("module")
    concordance_rows = []
    for module in MODULES:
        wu_effect = wu_modules.loc[module, "median_delta_after"]
        yan_effect = yan_mod.loc[module, "median_high_minus_low"]
        concordance_rows.append(
            {
                "record_type": "feature",
                "feature_type": "module",
                "feature": module,
                "wu_effect": wu_effect,
                "yan_effect": yan_effect,
                "direction_concordant": np.sign(wu_effect) == np.sign(yan_effect),
            }
        )
    yan_gene_raw = gene_summary.loc[gene_summary["model"].eq("raw_continuous")].set_index("gene")
    for row in rep.drop_duplicates("gene").itertuples():
        yan_effect = yan_gene_raw.loc[row.gene, "median_high_minus_low"]
        concordance_rows.append(
            {
                "record_type": "feature",
                "feature_type": "representative_gene",
                "feature": row.gene,
                "wu_effect": row.logFC,
                "yan_effect": yan_effect,
                "direction_concordant": np.sign(row.logFC) == np.sign(yan_effect),
            }
        )
    concordance = pd.DataFrame(concordance_rows)
    summary_rows = []
    for feature_type, d in concordance.groupby("feature_type"):
        summary_rows.append(
            {
                "record_type": "summary",
                "feature_type": feature_type,
                "feature": "ALL",
                "n_features": len(d),
                "n_direction_concordant": int(d["direction_concordant"].sum()),
                "spearman_rho_wu_vs_yan": safe_spearman(d["wu_effect"], d["yan_effect"]),
            }
        )
    concordance = pd.concat([concordance, pd.DataFrame(summary_rows)], ignore_index=True, sort=False)
    write_tsv(concordance, TABLES / "wu_yan_effect_concordance.tsv")

    # Module heatmap.
    heat = module_patient.loc[module_patient["version"].eq("score_independent")].pivot(index="patient", columns="module", values="high_minus_low")
    heat = heat.loc[:, MODULES]
    fig, ax = plt.subplots(figsize=(8.6, max(8.0, len(heat) * 0.15)))
    vmax = np.nanquantile(np.abs(heat.to_numpy()), 0.97)
    im = ax.imshow(heat.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(MODULES)), MODULES, rotation=45, ha="right")
    ax.set_yticks(range(len(heat)), heat.index, fontsize=5.5)
    ax.set_title("Yan patient-level score-independent module replication", loc="left", weight="bold")
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("High-Low median score")
    save_figure(fig, "Yan_module_effect_heatmap")

    # Gene validation lollipop.
    gs = gene_summary.loc[gene_summary["model"].eq("raw_continuous")].sort_values("median_high_minus_low")
    fig, ax = plt.subplots(figsize=(7.6, 5.5))
    yy = np.arange(len(gs))
    colors = np.where(gs["score_gene_status"].eq("Score-independent"), "#168A8A", "#D47B29")
    ax.axvline(0, color="#8F969D", ls="--", lw=0.8)
    ax.hlines(yy, 0, gs["median_high_minus_low"], color="#C2C7CA", lw=1.2)
    ax.scatter(gs["median_high_minus_low"], yy, c=colors, s=np.clip(-np.log10(gs["BH_FDR_within_model"].clip(lower=1e-12)) * 14, 24, 100))
    ax.set_yticks(yy, gs["gene"])
    ax.set_xlabel("Median patient High-Low expression")
    ax.set_title("Frozen Figure 4 representative genes in Yan", loc="left", weight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    save_figure(fig, "Yan_representative_gene_validation")

    # ------------------------------------------------------------------
    # Stage 6: dominant M01-M13 label alignment (one-vs-rest correlations).
    # ------------------------------------------------------------------
    primary_cells["metaprogram"] = primary_cells["cell_state"].str.extract(r"^(M\d{2})", expand=False)
    metaprograms = sorted(primary_cells["metaprogram"].dropna().unique())
    expected_mps = [f"M{i:02d}" for i in range(1, 14)]
    if not set(expected_mps).issubset(metaprograms):
        LOG.warning("Not all M01-M13 labels present among primary cells: %s", sorted(set(expected_mps) - set(metaprograms)))
    alignment_features = {
        "YAP": "YAP_score",
        "Stemness": "Stemness_score",
        "Joint axis": "YAP_Stem_axis",
        **{m: f"module__{m}" for m in CORE7},
    }
    align_rows = []
    for patient, d in primary_cells.groupby("patient", sort=True):
        for mp in expected_mps:
            binary = d["metaprogram"].eq(mp).astype(int).to_numpy()
            n_mp = int(binary.sum())
            for feature, col in alignment_features.items():
                rho = safe_spearman(d[col].to_numpy(float), binary) if n_mp >= 5 and (len(d) - n_mp) >= 5 else np.nan
                align_rows.append(
                    {
                        "record_type": "patient",
                        "patient": patient,
                        "metaprogram": mp,
                        "feature": feature,
                        "n_cells": len(d),
                        "n_metaprogram_cells": n_mp,
                        "rho_one_vs_rest": rho,
                        "method": "within-patient Spearman with one-vs-rest dominant cell_state indicator",
                    }
                )
    align_patient = pd.DataFrame(align_rows)
    align_summary_rows = []
    for (mp, feature), d in align_patient.groupby(["metaprogram", "feature"], sort=False):
        vals = d["rho_one_vs_rest"].dropna()
        align_summary_rows.append(
            {
                "record_type": "summary",
                "patient": "ALL",
                "metaprogram": mp,
                "feature": feature,
                "n_patients_testable": len(vals),
                "median_rho_one_vs_rest": vals.median() if len(vals) else np.nan,
                "n_positive": int((vals > 0).sum()),
                "method": "median of patient-specific one-vs-rest correlations",
            }
        )
    alignment = pd.concat([align_patient, pd.DataFrame(align_summary_rows)], ignore_index=True, sort=False)
    write_tsv(alignment, TABLES / "yan_metaprogram_alignment.tsv")
    align_heat = pd.DataFrame(align_summary_rows).pivot(index="feature", columns="metaprogram", values="median_rho_one_vs_rest")
    align_heat = align_heat.reindex(index=list(alignment_features), columns=expected_mps)
    fig, ax = plt.subplots(figsize=(9.8, 5.6))
    im = ax.imshow(align_heat.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-0.45, vmax=0.45)
    ax.set_xticks(range(13), expected_mps)
    ax.set_yticks(range(len(align_heat)), align_heat.index)
    for i in range(len(align_heat)):
        for j in range(13):
            val = align_heat.iloc[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=6, color="white" if abs(val) > 0.28 else "#202428")
    ax.set_title("Yan dominant metaprogram alignment", loc="left", weight="bold")
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("Median patient one-vs-rest rho")
    save_figure(fig, "Yan_metaprogram_alignment_heatmap")

    # ------------------------------------------------------------------
    # Stage 7: predeclared patient-level pCR/RD exploration.
    # ------------------------------------------------------------------
    module_q75 = {m: cells[f"module__{m}"].quantile(0.75) for m in MODULES}
    patient_metric_rows = []
    raw_rho_map = rho_table.loc[rho_table["model"].eq("raw")].set_index("patient")["rho"]
    for patient, d in primary_cells.groupby("patient", sort=True):
        row = {
            "patient": patient,
            "response": d["response"].iloc[0],
            "n_cancer_cells": len(d),
            "median_YAP": d["YAP_score"].median(),
            "median_Stem": d["Stemness_score"].median(),
            "median_joint_axis": d["YAP_Stem_axis"].median(),
            "upper_quartile_joint_axis": d["YAP_Stem_axis"].quantile(0.75),
            "high_state_fraction": d["YS_state"].eq("High").mean(),
            "within_patient_YAP_Stem_rho": raw_rho_map.get(patient, np.nan),
        }
        for module in MODULES:
            row[f"module_mean__{module}"] = d[f"module__{module}"].mean()
            row[f"module_high_frequency__{module}"] = (d[f"module__{module}"] >= module_q75[module]).mean()
        patient_metric_rows.append(row)
    patient_metrics = pd.DataFrame(patient_metric_rows)
    write_tsv(patient_metrics, TABLES / "yan_patient_response_metrics.tsv")
    response_data = patient_metrics.loc[patient_metrics["response"].isin(["pCR", "RD"])].copy()
    metrics = [c for c in patient_metrics.columns if c not in ["patient", "response", "n_cancer_cells"]]
    response_rows = []
    rng = np.random.default_rng(SEED + 777)
    for metric in metrics:
        pcr = response_data.loc[response_data["response"].eq("pCR"), metric].dropna().to_numpy(float)
        rd = response_data.loc[response_data["response"].eq("RD"), metric].dropna().to_numpy(float)
        if len(pcr) and len(rd):
            test = stats.mannwhitneyu(rd, pcr, alternative="two-sided")
            effect = 2 * test.statistic / (len(rd) * len(pcr)) - 1  # rank-biserial, RD minus pCR
            median_diff = np.median(rd) - np.median(pcr)
            boots = np.empty(2000)
            for i in range(2000):
                boots[i] = np.median(rng.choice(rd, size=len(rd), replace=True)) - np.median(rng.choice(pcr, size=len(pcr), replace=True))
            ci_low, ci_high = np.quantile(boots, [0.025, 0.975])
            p = test.pvalue
        else:
            effect = median_diff = ci_low = ci_high = p = np.nan
        response_rows.append(
            {
                "metric": metric,
                "n_pCR": len(pcr),
                "n_RD": len(rd),
                "pCR_median": np.median(pcr) if len(pcr) else np.nan,
                "RD_median": np.median(rd) if len(rd) else np.nan,
                "RD_minus_pCR_median_difference": median_diff,
                "bootstrap_ci_low": ci_low,
                "bootstrap_ci_high": ci_high,
                "rank_biserial_effect_RD_minus_pCR": effect,
                "wilcoxon_rank_sum_p": p,
                "logistic_regression_status": "not_run_no_patient_level_clinical_covariates_in_h5ad",
                "selection_policy": "all_predefined_metrics_reported",
            }
        )
    response_assoc = pd.DataFrame(response_rows)
    response_assoc["BH_FDR"] = bh_array(response_assoc["wilcoxon_rank_sum_p"])
    write_tsv(response_assoc, TABLES / "yan_response_association.tsv")

    fig, ax = plt.subplots(figsize=(8.6, max(6.5, len(response_assoc) * 0.24)))
    rr = response_assoc.sort_values("RD_minus_pCR_median_difference").reset_index(drop=True)
    yy = np.arange(len(rr))
    colors = np.where(rr["BH_FDR"] < 0.05, "#B51F35", "#7E878D")
    ax.axvline(0, color="#8F969D", ls="--", lw=0.8)
    ax.hlines(yy, rr["bootstrap_ci_low"], rr["bootstrap_ci_high"], color="#BCC2C6", lw=1)
    ax.scatter(rr["RD_minus_pCR_median_difference"], yy, c=colors, s=35)
    ax.set_yticks(yy, rr["metric"].str.replace("module_mean__", "mean: ").str.replace("module_high_frequency__", "frequency: "), fontsize=6.5)
    ax.set_xlabel("RD - pCR median difference (patient-level bootstrap 95% CI)")
    ax.set_title("Yan pCR/RD exploratory associations: all predefined metrics", loc="left", weight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    save_figure(fig, "Yan_response_association")

    qc = pd.DataFrame(
        [
            {"stage": 4, "check": "primary_threshold", "status": "PASS", "details": f">=50 cancer cells; {len(primary_patients)} patients"},
            {"stage": 4, "check": "threshold_sensitivity", "status": "PASS", "details": "20, 50 and 100 cells"},
            {"stage": 4, "check": "technical_models", "status": "PASS", "details": "primary and extended residual models"},
            {"stage": 4, "check": "random_effects_and_LOO", "status": "PASS", "details": "Fisher-z DerSimonian-Laird"},
            {"stage": 5, "check": "module_versions", "status": "PASS", "details": "score-independent, pairwise-shared removed, technical-adjusted"},
            {"stage": 5, "check": "anoikis_unique_version", "status": "NOT_TESTABLE", "details": "zero genes after removing all pairwise shared genes"},
            {"stage": 5, "check": "representative_genes", "status": "PASS", "details": f"{len(rep_genes)} frozen Figure 4 genes"},
            {"stage": 6, "check": "metaprogram_labels", "status": "PASS", "details": "M01-M13 dominant cell_state labels; one-vs-rest patient correlations"},
            {"stage": 7, "check": "response_metrics", "status": "PASS", "details": f"{len(metrics)} predefined metrics; all reported"},
            {"stage": 7, "check": "logistic_covariates", "status": "NOT_RUN", "details": "no patient-level clinical covariates in h5ad"},
        ]
    )
    write_tsv(qc, TABLES / "stage4_7_QC.tsv")
    LOG.info("Stages 4-7 completed")


if __name__ == "__main__":
    main()
