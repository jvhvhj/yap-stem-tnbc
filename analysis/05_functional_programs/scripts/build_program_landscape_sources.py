#!/usr/bin/env python
# Purpose: Build program-landscape and lineage source tables
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 4 and Supplementary Fig. S1.
from __future__ import annotations

import hashlib
import os
import re
from pathlib import Path

# Use the operating system environment without a machine-specific WINDIR override.
os.environ.setdefault(
    "MPLCONFIGDIR",
    r".\0723_Core_program_gene_closure_and_multicellular_context_gate\tmp\matplotlib",
)

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle


ROOT = Path(r".")
OUT = ROOT / "0723_Core_program_gene_closure_and_multicellular_context_gate"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
PLOTDATA = OUT / "plotting_data"
LOGS = OUT / "logs"

LANDSCAPE_SOURCE = (
    ROOT
    / "0720_StageA3_scientific_and_visual_closure_before_v6"
    / "tables"
    / "Figure3_module_joint_decile_plotting_data.tsv"
)
WU_SUPP = ROOT / "0716_yan2026_validation" / "inputs" / "Wu2021_Supplementary_Tables.xlsx"
HPA = ROOT / "0706_defensive_extension" / "references" / "proteinatlas" / "proteinatlas.tsv"
GENE_LIST_DIR = ROOT / "0719_frozen_score_definition_audit" / "signature_gene_lists"

for d in (TABLES, FIGURES, PLOTDATA, LOGS):
    d.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 7.0,
        "axes.titlesize": 8.0,
        "axes.labelsize": 7.0,
        "xtick.labelsize": 5.8,
        "ytick.labelsize": 6.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.5,
    }
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False, na_rep="NA")


def extract_target_ncpm(value: object, target: str) -> float:
    if not isinstance(value, str):
        return np.nan
    for item in value.split(";"):
        if item.startswith(target + ":"):
            try:
                return float(item.split(":", 1)[1].strip())
            except ValueError:
                return np.nan
    return np.nan


def read_frozen_union() -> tuple[set[str], set[str]]:
    score = set()
    for name in ("YAP_clean__clean_score.tsv", "Stemness_clean__clean_score.tsv"):
        d = pd.read_csv(GENE_LIST_DIR / name, sep="\t")
        score.update(d["gene"].astype(str))
    core = set()
    for name in (
        "UPR__score_independent.tsv",
        "TNFA_NFKB__score_independent.tsv",
        "Hypoxia__score_independent.tsv",
        "Adhesion_Remodeling__score_independent.tsv",
        "Wound_Healing__score_independent.tsv",
        "Survival_Stress__score_independent.tsv",
    ):
        d = pd.read_csv(GENE_LIST_DIR / name, sep="\t")
        core.update(d["gene"].astype(str))
    return score, core


def build_lineage_manifest() -> pd.DataFrame:
    score_genes, core_genes = read_frozen_union()
    rows: list[dict[str, object]] = []

    wu = pd.read_excel(WU_SUPP, sheet_name="Supplementary Table 10", header=3)
    wu = wu.loc[
        wu["major lineage"].eq("Normal Epithelial")
        & (pd.to_numeric(wu["p_val_adj"], errors="coerce") < 0.05)
        & (pd.to_numeric(wu["pct.1"], errors="coerce") >= 0.25)
        & (pd.to_numeric(wu["avg_logFC"], errors="coerce") > 0)
    ].copy()
    wu["avg_logFC"] = pd.to_numeric(wu["avg_logFC"], errors="coerce")
    wu["pct.1"] = pd.to_numeric(wu["pct.1"], errors="coerce")

    wu_map = {
        "Basal/myoepithelial": "Myoepithelial",
        "Luminal progenitor/secretory luminal": "Luminal Progenitors",
        "Mature luminal": "Mature Luminal",
    }
    for signature, cluster in wu_map.items():
        selected = (
            wu.loc[wu["cluster"].eq(cluster)]
            .sort_values(["avg_logFC", "pct.1", "gene"], ascending=[False, False, True])
            .drop_duplicates("gene")
            .head(20)
        )
        for rank, r in enumerate(selected.itertuples(index=False), start=1):
            rows.append(
                {
                    "signature": signature,
                    "gene": str(r.gene),
                    "rank": rank,
                    "source_dataset": "Wu et al. 2021 Supplementary Table 10",
                    "source_file": str(WU_SUPP),
                    "source_field": f"Normal Epithelial / {cluster}",
                    "selection_rule": (
                        "published within-lineage marker; adjusted P<0.05; pct.1>=0.25; "
                        "avg_logFC>0; top 20 by avg_logFC then pct.1"
                    ),
                    "source_effect": float(r.avg_logFC),
                    "source_detection": float(getattr(r, "_2", np.nan))
                    if False
                    else float(r._asdict().get("pct.1", np.nan)),
                }
            )

    hpa_cols = [
        "Gene",
        "RNA single cell type specificity",
        "RNA single cell type specific nCPM",
    ]
    hpa = pd.read_csv(HPA, sep="\t", usecols=hpa_cols, low_memory=False)
    target = "Breast hormone-responsive cells"
    hpa["target_nCPM"] = hpa["RNA single cell type specific nCPM"].map(
        lambda x: extract_target_ncpm(x, target)
    )
    selected = (
        hpa.loc[
            hpa["target_nCPM"].notna()
            & hpa["RNA single cell type specificity"].isin(
                ["Cell type enriched", "Group enriched", "Cell type enhanced"]
            )
        ]
        .sort_values(["target_nCPM", "Gene"], ascending=[False, True])
        .drop_duplicates("Gene")
        .head(20)
    )
    for rank, r in enumerate(selected.itertuples(index=False), start=1):
        rows.append(
            {
                "signature": "Hormone-responsive luminal",
                "gene": str(r.Gene),
                "rank": rank,
                "source_dataset": "Human Protein Atlas single-cell type specificity",
                "source_file": str(HPA),
                "source_field": target,
                "selection_rule": (
                    "published cell-type enriched/group-enriched/enhanced entry; "
                    "top 20 by target single-cell nCPM"
                ),
                "source_effect": float(r.target_nCPM),
                "source_detection": np.nan,
            }
        )

    raw = pd.DataFrame(rows)
    raw["version"] = "raw_predefined"
    raw["excluded_overlap"] = False
    raw["overlap_class"] = np.select(
        [
            raw["gene"].isin(score_genes) & raw["gene"].isin(core_genes),
            raw["gene"].isin(score_genes),
            raw["gene"].isin(core_genes),
        ],
        ["score_and_core", "score_gene", "six_core_program"],
        default="none",
    )
    clean = raw.loc[raw["overlap_class"].eq("none")].copy()
    clean["version"] = "core_overlap_removed"
    manifest = pd.concat([raw, clean], ignore_index=True)
    manifest["source_file_sha256"] = manifest["source_file"].map(
        {str(WU_SUPP): sha256(WU_SUPP), str(HPA): sha256(HPA)}
    )
    manifest["definition_frozen_before_outcome_analysis"] = True
    manifest["notes"] = np.where(
        manifest["version"].eq("raw_predefined"),
        "Primary lineage signature for confounding audit.",
        "Sensitivity version excludes every YAP/Stem score gene and every gene in the six frozen core programs.",
    )
    return manifest


def plot_landscape() -> pd.DataFrame:
    d = pd.read_csv(LANDSCAPE_SOURCE, sep="\t")
    d = d.loc[d["record_type"].eq("patient_decile")].copy()
    if len(d) != 480:
        raise RuntimeError(f"Expected 480 patient-decile rows, found {len(d)}")
    patients = list(dict.fromkeys(d["patient"].astype(str)))
    if len(patients) != 8:
        raise RuntimeError(f"Expected 8 patients, found {len(patients)}")
    modules = [
        "UPR",
        "TNFA NFKB",
        "Hypoxia",
        "Survival Stress",
        "Adhesion Remodeling",
        "Wound Healing",
    ]
    display = {
        "UPR": "UPR",
        "TNFA NFKB": "TNFα–NF-κB",
        "Hypoxia": "Hypoxia",
        "Survival Stress": "Survival stress",
        "Adhesion Remodeling": "Adhesion remodeling",
        "Wound Healing": "Wound healing",
    }
    branch = {
        "UPR": "Stress-adaptive",
        "TNFA NFKB": "Stress-adaptive",
        "Hypoxia": "Stress-adaptive",
        "Survival Stress": "Stress-adaptive",
        "Adhesion Remodeling": "Adhesion-remodeling",
        "Wound Healing": "Adhesion-remodeling",
    }
    expected = pd.MultiIndex.from_product(
        [modules, patients, range(1, 11)], names=["module", "patient", "joint_decile"]
    )
    x = d.set_index(["module", "patient", "joint_decile"]).reindex(expected)
    if x["median_score_after"].isna().any():
        raise RuntimeError("Landscape contains missing module×patient×decile cells")
    matrix = x["median_score_after"].to_numpy().reshape(len(modules), len(patients) * 10)

    fig = plt.figure(figsize=(180 / 25.4, 105 / 25.4), constrained_layout=False)
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[0.028, 0.972],
        height_ratios=[0.17, 0.83],
        left=0.19,
        right=0.97,
        bottom=0.20,
        top=0.83,
        wspace=0.03,
        hspace=0.05,
    )
    ax_branch = fig.add_subplot(gs[1, 0])
    ax_top = fig.add_subplot(gs[0, 1])
    ax = fig.add_subplot(gs[1, 1])

    cmap = LinearSegmentedColormap.from_list(
        "cohort_z", ["#235789", "#DCE6ED", "#F6F4EF", "#F5C4A8", "#B53A43"]
    )
    norm = TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5)
    im = ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)
    ax.set_yticks(np.arange(len(modules)))
    ax.set_yticklabels([display[m] for m in modules])
    ax.tick_params(axis="y", length=0, pad=5)
    ax.set_xticks([])
    for i in range(1, 8):
        ax.axvline(i * 10 - 0.5, color="white", lw=1.5)
    for y in [3.5]:
        ax.axhline(y, color="white", lw=1.6)
    for side in ax.spines.values():
        side.set_visible(False)

    ax_branch.set_xlim(0, 1)
    ax_branch.set_ylim(5.5, -0.5)
    ax_branch.add_patch(Rectangle((0.05, -0.48), 0.68, 3.96, color="#778FA7", lw=0))
    ax_branch.add_patch(Rectangle((0.05, 3.52), 0.68, 1.96, color="#C17B63", lw=0))
    ax_branch.axis("off")

    ax_top.set_xlim(-0.5, 79.5)
    ax_top.set_ylim(0, 1)
    for i, patient in enumerate(patients):
        ax_top.text(i * 10 + 4.5, 0.72, patient, ha="center", va="center", fontsize=6.3)
        for decile in (1, 5, 10):
            ax_top.text(
                i * 10 + decile - 1,
                0.18,
                str(decile),
                ha="center",
                va="center",
                fontsize=5.0,
                color="#555555",
            )
    for i in range(1, 8):
        ax_top.axvline(i * 10 - 0.5, color="#BFC4C8", lw=0.45)
    ax_top.axis("off")

    fig.text(
        0.19,
        0.93,
        "Score-independent programs form a patient-resolved continuous landscape",
        ha="left",
        va="top",
        fontsize=9.2,
        fontweight="bold",
    )
    fig.text(
        0.19,
        0.885,
        "6 programs × 8 patients × 10 Joint-axis deciles; descriptive medians only",
        ha="left",
        va="top",
        fontsize=6.3,
        color="#555555",
    )
    fig.add_artist(
        Rectangle(
            (0.735, 0.875),
            0.010,
            0.018,
            transform=fig.transFigure,
            facecolor="#778FA7",
            edgecolor="none",
        )
    )
    fig.text(0.749, 0.884, "Stress-adaptive", fontsize=5.8, va="center", color="#344B61")
    fig.add_artist(
        Rectangle(
            (0.835, 0.875),
            0.010,
            0.018,
            transform=fig.transFigure,
            facecolor="#C17B63",
            edgecolor="none",
        )
    )
    fig.text(
        0.849,
        0.884,
        "Adhesion-remodeling",
        fontsize=5.8,
        va="center",
        color="#7D4939",
    )
    arrow_ax = fig.add_axes([0.19, 0.065, 0.36, 0.065])
    arrow_ax.annotate(
        "",
        xy=(0.98, 0.58),
        xytext=(0.02, 0.58),
        arrowprops=dict(arrowstyle="-|>", lw=0.85, color="#444444"),
    )
    arrow_ax.text(
        0.50,
        0.02,
        "Joint YAP–Stem axis: low → high within each patient",
        ha="center",
        va="bottom",
        fontsize=5.8,
        color="#444444",
    )
    arrow_ax.axis("off")

    cax = fig.add_axes([0.70, 0.075, 0.22, 0.025])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", ticks=[-2, 0, 2])
    cb.set_label("Median module score (shared cohort z-scale)", fontsize=5.8, labelpad=1)
    cb.ax.tick_params(labelsize=5.3, length=2, pad=1)
    cb.outline.set_linewidth(0.4)

    fig.text(0.018, 0.955, "a", fontsize=10.5, fontweight="bold", va="top")
    pdf = FIGURES / "Figure3_program_landscape_storyboard.pdf"
    png = FIGURES / "Figure3_program_landscape_storyboard.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=400)
    plt.close(fig)

    plot = d.copy()
    plot["module_display"] = plot["module"].map(display)
    plot["branch_frozen"] = plot["module"].map(branch)
    plot["patient_order"] = plot["patient"].map({p: i + 1 for i, p in enumerate(patients)})
    plot["module_order"] = plot["module"].map({m: i + 1 for i, m in enumerate(modules)})
    plot["colour_scale_vmin"] = -2.5
    plot["colour_scale_vmax"] = 2.5
    plot["statistical_role"] = "descriptive_only_no_new_trend_P"
    return plot


def main() -> None:
    lineage = build_lineage_manifest()
    write_tsv(lineage, TABLES / "lineage_signature_manifest.tsv")
    landscape = plot_landscape()
    write_tsv(landscape, PLOTDATA / "Figure3_program_landscape_plotting_data.tsv")

    provenance = pd.DataFrame(
        [
            {
                "artifact": "Figure3_program_landscape_storyboard",
                "source_file": str(LANDSCAPE_SOURCE),
                "source_sha256": sha256(LANDSCAPE_SOURCE),
                "rows_used": len(landscape),
                "definition": "approved patient×decile×six-program dataset; no recomputation",
            },
            {
                "artifact": "lineage_signature_manifest_Wu",
                "source_file": str(WU_SUPP),
                "source_sha256": sha256(WU_SUPP),
                "rows_used": int(lineage["source_file"].eq(str(WU_SUPP)).sum()),
                "definition": "deterministic published-marker selection before outcome analysis",
            },
            {
                "artifact": "lineage_signature_manifest_HPA",
                "source_file": str(HPA),
                "source_sha256": sha256(HPA),
                "rows_used": int(lineage["source_file"].eq(str(HPA)).sum()),
                "definition": "deterministic published single-cell specificity selection before outcome analysis",
            },
        ]
    )
    write_tsv(provenance, TABLES / "stage1_lineage_provenance.tsv")
    (LOGS / "01_program_landscape_and_lineage_manifest.log").write_text(
        "PASS\n"
        f"Landscape rows: {len(landscape)}\n"
        f"Lineage manifest rows: {len(lineage)}\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
