#!/usr/bin/env python
# Purpose: Patient-level external scoring robustness
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 5; orthogonal transfer in Fig. 4.
"""Yan conventional/patient-aware feature effects using the frozen 182-gene cache."""

from __future__ import annotations

import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats


ROOT = Path(os.environ.get("AHIPPO_YAP_ROOT", Path.cwd()))
OUT = ROOT / "0717_methodological_framework_benchmark"
TABLES = OUT / "tables"
LOGS = OUT / "logs"
TABLES.mkdir(parents=True, exist_ok=True)
LOGS.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[logging.FileHandler(LOGS / "03_yan_framework_effects.log", encoding="utf-8"), logging.StreamHandler(sys.stdout)],
)
LOG = logging.getLogger("yan_effects")

CACHE = ROOT / "0716_yan2026_validation" / "inputs" / "yan_cancer_cell_scores.tsv.gz"
EXPR = ROOT / "0716_yan2026_validation" / "inputs" / "yan_cancer_frozen_gene_expression.npz"
COLS = ROOT / "0716_yan2026_validation" / "inputs" / "yan_cancer_frozen_gene_expression_columns.tsv"
MANIFEST = ROOT / "0716_yan2026_validation" / "tables" / "frozen_signature_manifest.tsv"
MODULES = ["UPR", "TNFA NFKB", "Hypoxia", "Adhesion Remodeling", "Wound Healing", "Survival Stress", "Anoikis Resistance", "Integrin Adhesion"]


def write_tsv(frame: pd.DataFrame, name: str) -> None:
    path = TABLES / name
    frame.to_csv(path, sep="\t", index=False, na_rep="NA")
    LOG.info("Wrote %s (%s rows)", path, len(frame))


def bh(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    out = np.full(len(values), np.nan)
    keep = np.isfinite(values)
    p = values[keep]
    if len(p):
        order = np.argsort(p)
        ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        temp = np.empty(len(p))
        temp[order] = np.minimum(ranked, 1)
        out[keep] = temp
    return out


def zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    sd = np.nanstd(values, ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(values)
    return (values - np.nanmean(values)) / sd


def residualize(values: np.ndarray, covariates: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    covariates = np.asarray(covariates, dtype=float)
    keep = np.nanstd(covariates, axis=0, ddof=1) > 1e-12
    covariates = covariates[:, keep]
    if covariates.shape[1]:
        covariates = np.column_stack([zscore(covariates[:, j]) for j in range(covariates.shape[1])])
    design = np.column_stack([np.ones(len(values)), covariates])
    q, _ = np.linalg.qr(design, mode="reduced")
    return values - q @ (q.T @ values)


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 5 or np.std(x[keep]) == 0 or np.std(y[keep]) == 0:
        return np.nan
    return float(stats.spearmanr(x[keep], y[keep]).statistic)


def safe_wilcoxon(x: np.ndarray, y: np.ndarray | None = None) -> float:
    try:
        if y is None:
            return float(stats.wilcoxon(x).pvalue)
        return float(stats.mannwhitneyu(x, y, alternative="two-sided").pvalue)
    except ValueError:
        return np.nan


def main() -> None:
    for path in (CACHE, EXPR, COLS, MANIFEST):
        if not path.exists():
            raise FileNotFoundError(path)
    cache = pd.read_csv(CACHE, sep="\t")
    expr = sparse.load_npz(EXPR).tocsr()
    columns = pd.read_csv(COLS, sep="\t")
    manifest = pd.read_csv(MANIFEST, sep="\t")
    if expr.shape[0] != len(cache):
        raise RuntimeError("Frozen Yan expression cache is not row-aligned")
    gene_to_col = dict(zip(columns["requested_gene"].astype(str), columns["matrix_column"].astype(int)))
    eligible_patients = cache.groupby("patient", sort=True).size()
    eligible_patients = eligible_patients[eligible_patients >= 50].index.astype(str)
    eligible_mask = cache["patient"].astype(str).isin(eligible_patients).to_numpy()

    module_manifest = manifest.loc[manifest["signature_type"].eq("module") & manifest["module"].isin(MODULES)].copy()
    module_sets: dict[tuple[str, str], list[str]] = {}
    for module in MODULES:
        part = module_manifest.loc[module_manifest["module"].eq(module)]
        module_sets[(module, "original")] = part["gene"].astype(str).drop_duplicates().tolist()
        score_independent = part["score_independent_gene"].astype(str).str.lower().eq("true")
        pairwise_unique = part["pairwise_unique_gene"].astype(str).str.lower().eq("true")
        module_sets[(module, "score_independent")] = part.loc[score_independent, "gene"].astype(str).drop_duplicates().tolist()
        module_sets[(module, "pairwise_unique")] = part.loc[pairwise_unique, "gene"].astype(str).drop_duplicates().tolist()

    rep_manifest = manifest.loc[manifest["signature_type"].eq("representative_gene")].drop_duplicates("gene")
    rep_genes = [g for g in rep_manifest["gene"].astype(str) if g in gene_to_col]

    def feature_score(genes: list[str]) -> np.ndarray:
        positions = [gene_to_col[g] for g in genes if g in gene_to_col]
        if not positions:
            return np.full(len(cache), np.nan)
        return np.asarray(expr[:, positions].mean(axis=1)).ravel()

    vectors: dict[tuple[str, str, str], np.ndarray] = {}
    count_rows = []
    for (module, version), genes in module_sets.items():
        vectors[("module", module, version)] = feature_score(genes)
        count_rows.append({"feature_type": "module", "feature": module, "version": version, "n_genes": len(genes), "n_genes_present": sum(g in gene_to_col for g in genes)})
    for gene in rep_genes:
        vectors[("representative_gene", gene, "gene")] = np.asarray(expr[:, gene_to_col[gene]].toarray()).ravel()
    write_tsv(pd.DataFrame(count_rows), "yan_framework_feature_gene_counts.tsv")

    joint = cache["YAP_Stem_axis"].to_numpy(float)
    global_group = np.where(joint >= np.nanmedian(joint), "High", "Low")
    yap = cache["YAP_score"].to_numpy(float)
    stem = cache["Stemness_score"].to_numpy(float)
    hh_group = np.where((yap >= np.nanmedian(yap)) & (stem >= np.nanmedian(stem)), "HH", "Other")
    state = cache["YS_state"].astype(str).to_numpy()

    rows: list[dict] = []
    for (feature_type, feature, version), values in vectors.items():
        if not np.isfinite(values).any():
            continue
        for workflow, group, high_label, low_label in (
            ("conventional_global_median", global_group, "High", "Low"),
            ("conventional_HH_other", hh_group, "HH", "Other"),
        ):
            high = values[eligible_mask & (group == high_label) & np.isfinite(values)]
            low = values[eligible_mask & (group == low_label) & np.isfinite(values)]
            rows.append(
                {
                    "record_type": "pooled", "cohort": "Yan2026", "feature_type": feature_type,
                    "feature": feature, "workflow": workflow, "version": version, "patient": "POOLED",
                    "n_cells": len(high) + len(low), "n_high": len(high), "n_low": len(low),
                    "effect": float(np.median(high) - np.median(low)), "p_value": safe_wilcoxon(high, low),
                    "n_patients": np.nan, "n_positive": np.nan, "BH_FDR": np.nan,
                }
            )

        for workflow in ("patient_raw", "technical_adjusted"):
            patient_effects = []
            for patient in eligible_patients:
                idx = np.flatnonzero(cache["patient"].astype(str).to_numpy() == patient)
                val = values[idx]
                if workflow == "technical_adjusted":
                    covariates = np.column_stack(
                        [
                            np.log1p(cache.loc[idx, "nCount_RNA"].to_numpy(float)),
                            cache.loc[idx, "nFeature_RNA"].to_numpy(float),
                            cache.loc[idx, "S_score"].to_numpy(float),
                            cache.loc[idx, "G2M_score"].to_numpy(float),
                        ]
                    )
                    val = residualize(val, covariates)
                high = val[(state[idx] == "High") & np.isfinite(val)]
                low = val[(state[idx] == "Low") & np.isfinite(val)]
                effect = float(np.median(high) - np.median(low)) if len(high) and len(low) else np.nan
                patient_effects.append(effect)
                rows.append(
                    {
                        "record_type": "patient", "cohort": "Yan2026", "feature_type": feature_type,
                        "feature": feature, "workflow": workflow, "version": version, "patient": patient,
                        "n_cells": len(idx), "n_high": len(high), "n_low": len(low), "effect": effect,
                        "p_value": np.nan, "n_patients": np.nan, "n_positive": np.nan, "BH_FDR": np.nan,
                    }
                )
            effects = np.asarray(patient_effects, dtype=float)
            effects = effects[np.isfinite(effects)]
            rows.append(
                {
                    "record_type": "summary", "cohort": "Yan2026", "feature_type": feature_type,
                    "feature": feature, "workflow": workflow, "version": version, "patient": "ALL",
                    "n_cells": int(eligible_mask.sum()), "n_high": int(np.sum(eligible_mask & (state == "High"))),
                    "n_low": int(np.sum(eligible_mask & (state == "Low"))),
                    "effect": float(np.median(effects)) if len(effects) else np.nan,
                    "p_value": safe_wilcoxon(effects), "n_patients": len(effects),
                    "n_positive": int((effects > 0).sum()), "BH_FDR": np.nan,
                }
            )

    effect_df = pd.DataFrame(rows)
    summary_mask = effect_df["record_type"].eq("summary")
    for _, idx in effect_df.loc[summary_mask].groupby(["workflow", "feature_type", "version"]).groups.items():
        effect_df.loc[idx, "BH_FDR"] = bh(effect_df.loc[idx, "p_value"].to_numpy(float))
    write_tsv(effect_df, "framework_feature_effects_yan.tsv")

    # Continuous patient-wise module correlations for raw and adjusted axes.
    correlation_rows = []
    for version in ("score_independent", "pairwise_unique"):
        for module in MODULES:
            values = vectors[("module", module, version)]
            if not np.isfinite(values).any():
                continue
            for patient in eligible_patients:
                idx = np.flatnonzero(cache["patient"].astype(str).to_numpy() == patient)
                covariates = np.column_stack(
                    [
                        np.log1p(cache.loc[idx, "nCount_RNA"].to_numpy(float)),
                        cache.loc[idx, "nFeature_RNA"].to_numpy(float),
                        cache.loc[idx, "S_score"].to_numpy(float),
                        cache.loc[idx, "G2M_score"].to_numpy(float),
                    ]
                )
                for model in ("raw", "technical_adjusted"):
                    y = yap[idx].copy(); s = stem[idx].copy(); val = values[idx].copy()
                    if model == "technical_adjusted":
                        y = residualize(y, covariates); s = residualize(s, covariates); val = residualize(val, covariates)
                    correlation_rows.append(
                        {
                            "cohort": "Yan2026", "patient": patient, "module": module,
                            "version": version, "model": model, "n_cells": len(idx),
                            "rho_with_joint": spearman(zscore(y) + zscore(s), val),
                        }
                    )
    write_tsv(pd.DataFrame(correlation_rows), "framework_module_correlations_yan.tsv")

    technical_rows = []
    score_vars = {"YAP": yap, "Stemness": stem, "Joint": joint}
    tech_vars = {
        "log_nCount": np.log1p(cache["nCount_RNA"].to_numpy(float)),
        "nFeature": cache["nFeature_RNA"].to_numpy(float),
        "S_score": cache["S_score"].to_numpy(float),
        "G2M_score": cache["G2M_score"].to_numpy(float),
    }
    for scope in ["POOLED", *eligible_patients]:
        idx = np.flatnonzero(eligible_mask) if scope == "POOLED" else np.flatnonzero(cache["patient"].astype(str).to_numpy() == scope)
        for score_name, score_values in score_vars.items():
            for covariate, cov_values in tech_vars.items():
                technical_rows.append(
                    {
                        "cohort": "Yan2026", "scope": "pooled" if scope == "POOLED" else "patient",
                        "patient": scope, "score": score_name, "covariate": covariate,
                        "n_cells": len(idx), "spearman_rho": spearman(score_values[idx], cov_values[idx]),
                    }
                )
    write_tsv(pd.DataFrame(technical_rows), "score_technical_association_yan.tsv")
    LOG.info("Yan framework effects completed")


if __name__ == "__main__":
    main()
