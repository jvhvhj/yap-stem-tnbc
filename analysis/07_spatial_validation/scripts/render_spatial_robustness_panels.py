# Purpose: Figure 6B/C visual render
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 6 and Supplementary Fig. S4.
from __future__ import annotations

import hashlib
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import pandas as pd


ROOT = Path(r".")
OUT = ROOT / "Figure6BC_visual_reset"
SOURCE = OUT / "Figure6BC_visual_reset_source.tsv"
B_PDF = OUT / "Figure6B_adjustment_landscape.pdf"
B_PNG = OUT / "Figure6B_adjustment_landscape.png"
C_PDF = OUT / "Figure6C_section_agreement_landscape.pdf"
C_PNG = OUT / "Figure6C_section_agreement_landscape.png"
BC_PDF = OUT / "Figure6BC_combined_preview.pdf"
BC_PNG = OUT / "Figure6BC_combined_preview.png"
QC = OUT / "Figure6BC_visual_reset_QC.md"
SEMANTIC = OUT / "Figure6BC_semantic_readout.md"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


OUT.mkdir(parents=True, exist_ok=True)
src = pd.read_csv(SOURCE, sep="\t", dtype=str, keep_default_na=False)
b_data = src.loc[src["record_type"].eq("PATIENT_ADJUSTMENT_LANDSCAPE")].copy()
c_data = src.loc[src["record_type"].eq("PAIRED_SECTION_AGREEMENT")].copy()
b_summary = src.loc[(src["panel"].eq("B")) & (src["record_type"].eq("PANEL_SUMMARY"))].iloc[0]
c_summary = src.loc[(src["panel"].eq("C")) & (src["record_type"].eq("PANEL_SUMMARY"))].iloc[0]

for column in ["display_order", "rho_raw", "rho_adjusted", "delta_rho"]:
    b_data[column] = pd.to_numeric(b_data[column], errors="raise")
for column in ["display_order", "section_1_rho", "section_2_rho", "pair_mean_rho", "pair_difference_rho"]:
    c_data[column] = pd.to_numeric(c_data[column], errors="raise")

b_data = b_data.sort_values("display_order", kind="stable")
c_data = c_data.sort_values("display_order", kind="stable")

assert b_data.shape[0] == 22
assert c_data.shape[0] == 21
assert b_data["rho_adjusted"].gt(0).all()
assert c_data["section_1_rho"].gt(0).all() and c_data["section_2_rho"].gt(0).all()
assert np.allclose(b_data["delta_rho"], b_data["rho_adjusted"] - b_data["rho_raw"], rtol=0, atol=1e-14)
assert np.allclose(c_data["pair_mean_rho"], (c_data["section_1_rho"] + c_data["section_2_rho"]) / 2, rtol=0, atol=1e-14)
assert np.allclose(c_data["pair_difference_rho"], c_data["section_2_rho"] - c_data["section_1_rho"], rtol=0, atol=1e-14)

cohort_n = int(float(b_summary["cohort_n"]))
positive_n = int(float(b_summary["positive_n"]))
cohort_median = float(b_summary["cohort_median_adjusted_rho"])
cohort_p = float(b_summary["cohort_p"])
paired_n = int(float(c_summary["paired_n"]))
paired_positive_n = int(float(c_summary["paired_both_positive_n"]))
median_pair_difference = float(c_summary["median_pair_difference_rho"])

BLUE = "#2F6DB3"
BLUE_LIGHT = "#AFC8E4"
CORAL = "#D7605C"
CORAL_DARK = "#8A3936"
ZERO = "#F7F7F7"
DARK = "#26323A"
MUTED = "#69777F"
GRID = "#D1D8DC"
SEQUENTIAL = LinearSegmentedColormap.from_list("positive_rho", [ZERO, "#F2C4C0", CORAL], N=256)

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 6.5,
    "axes.labelsize": 6.5,
    "xtick.labelsize": 5.9,
    "ytick.labelsize": 5.9,
    "axes.linewidth": 0.55,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.facecolor": "white",
    "savefig.transparent": False,
    "figure.facecolor": "white",
})


def panel_header(ax, letter, title, statistics):
    ax.axis("off")
    ax.text(0.0, 0.96, letter, ha="left", va="top", fontsize=10.5, fontweight="bold", color="#111111")
    ax.text(0.075, 0.94, title, ha="left", va="top", fontsize=7.8, fontweight="bold", color=DARK, linespacing=1.05)
    ax.text(0.075, 0.04, statistics, ha="left", va="bottom", fontsize=6.0, color=MUTED, linespacing=1.12)


def clean_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(DARK)
    ax.spines["bottom"].set_color(DARK)
    ax.tick_params(length=2.0, width=0.5, colors=DARK, pad=1.5)


def draw_b(fig, slot):
    spec = slot.subgridspec(2, 1, height_ratios=[0.22, 0.78], hspace=0.015)
    header = fig.add_subplot(spec[0, 0])
    panel_header(
        header,
        "B",
        "Patient-level robustness of spatial\nYAP-Stem programme coupling",
        f"n = {cohort_n} · {positive_n}/{cohort_n} adjusted positive\n"
        f"median adjusted ρ = {cohort_median:.3f} · P = 2.38 × 10$^{{-7}}$",
    )

    ax = fig.add_subplot(spec[1, 0])
    x = b_data["rho_adjusted"].to_numpy(float)
    y = b_data["delta_rho"].to_numpy(float)
    strengthened = y > 0
    colors = np.where(strengthened, CORAL, BLUE)
    ax.scatter(x, y, s=24, c=colors, edgecolors="white", linewidths=0.55, zorder=3)
    ax.axhline(0, color="#7E8A91", lw=0.75, ls=(0, (2.4, 2.4)), zorder=0)
    ax.axvline(cohort_median, color=CORAL, lw=0.75, ls=(0, (3.2, 2.2)), alpha=0.85, zorder=0)
    ax.text(
        cohort_median + 0.004, -0.229, f"median {cohort_median:.3f}",
        ha="left", va="bottom", rotation=90, fontsize=5.3, color=CORAL_DARK,
    )
    ax.text(0.252, 0.028, "strengthened", ha="left", va="center", fontsize=5.5, color=CORAL_DARK)
    ax.text(0.252, -0.227, "attenuated", ha="left", va="bottom", fontsize=5.5, color=BLUE)

    x_min, x_max = 0.245, 0.575
    y_min, y_max = -0.235, 0.045
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    for value in x:
        ax.plot([value, value], [y_min, y_min + 0.006], color=MUTED, lw=0.45, alpha=0.55, clip_on=True)
    for value in y:
        ax.plot([x_max - 0.006, x_max], [value, value], color=MUTED, lw=0.45, alpha=0.55, clip_on=True)

    offsets = {
        "maximum_adjusted_rho": (-7, 7, "right", "bottom"),
        "minimum_adjusted_rho": (7, -5, "left", "top"),
        "largest_positive_delta": (7, 6, "left", "bottom"),
        "largest_negative_delta": (7, 6, "left", "bottom"),
    }
    for _, row in b_data.loc[b_data["display_label"].ne("")].iterrows():
        first_rule = row["deterministic_label_reason"].split(";")[0]
        dx, dy, ha, va = offsets[first_rule]
        ax.annotate(
            row["patient_id"], (row["rho_adjusted"], row["delta_rho"]),
            xytext=(dx, dy), textcoords="offset points", ha=ha, va=va,
            fontsize=5.8, color=DARK,
            arrowprops=dict(arrowstyle="-", color=GRID, lw=0.45, shrinkA=1.5, shrinkB=2.0),
        )

    ax.set_xlabel("Adjusted spatial Spearman ρ", labelpad=2.0, color=DARK)
    ax.set_ylabel("Δρ (adjusted - raw)", labelpad=2.0, color=DARK)
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    clean_axes(ax)
    return ax


def draw_c(fig, slot):
    spec = slot.subgridspec(2, 1, height_ratios=[0.20, 0.80], hspace=0.015)
    header = fig.add_subplot(spec[0, 0])
    panel_header(
        header,
        "C",
        "Paired-section recurrence of spatial\nYAP-Stem programme coupling",
        f"{paired_positive_n}/{paired_n} paired patients · both sections positive",
    )

    body = spec[1, 0].subgridspec(1, 2, width_ratios=[0.79, 0.21], wspace=0.08)
    ax = fig.add_subplot(body[0, 0])
    x = c_data["pair_mean_rho"].to_numpy(float)
    y = c_data["pair_difference_rho"].to_numpy(float)
    section2_stronger = y > 0
    colors = np.where(section2_stronger, CORAL, BLUE)
    ax.scatter(x, y, s=25, marker="D", c=colors, edgecolors="white", linewidths=0.55, zorder=3)
    ax.axhline(0, color="#7E8A91", lw=0.75, ls=(0, (2.4, 2.4)), zorder=0)
    ax.axhline(median_pair_difference, color=MUTED, lw=0.65, ls=(0, (3.0, 2.4)), alpha=0.75, zorder=0)
    ax.text(0.565, median_pair_difference - 0.005, f"median difference {median_pair_difference:.3f}",
            ha="right", va="top", fontsize=5.3, color=MUTED)
    ax.text(0.252, 0.215, "Section 2 stronger", ha="left", va="top", fontsize=5.4, color=CORAL_DARK)
    ax.text(0.252, -0.150, "Section 1 stronger", ha="left", va="bottom", fontsize=5.4, color=BLUE)

    x_min, x_max = 0.245, 0.575
    y_min, y_max = -0.160, 0.225
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    for value in x:
        ax.plot([value, value], [y_min, y_min + 0.008], color=MUTED, lw=0.45, alpha=0.50, clip_on=True)
    for value in y:
        ax.plot([x_max - 0.007, x_max], [value, value], color=MUTED, lw=0.45, alpha=0.50, clip_on=True)

    ax.set_xlabel("Mean spatial Spearman ρ", labelpad=2.0, color=DARK)
    ax.set_ylabel("Section 2 - Section 1 ρ", labelpad=2.0, color=DARK)
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    clean_axes(ax)

    strip = fig.add_subplot(body[0, 1])
    matrix = c_data[["section_1_rho", "section_2_rho"]].to_numpy(float)
    strip.imshow(matrix, cmap=SEQUENTIAL, norm=Normalize(0, 0.67), aspect="auto", interpolation="nearest")
    strip.set_xticks([0, 1])
    strip.set_xticklabels(["S1", "S2"], fontsize=5.4)
    strip.xaxis.tick_top()
    strip.set_yticks(np.arange(c_data.shape[0]))
    strip.set_yticklabels(c_data["patient_id"].tolist(), fontsize=5.0)
    strip.tick_params(axis="x", length=0, pad=2.0, colors=DARK)
    strip.tick_params(axis="y", length=0, pad=1.8, colors=DARK)
    strip.set_title("Original ρ", fontsize=5.7, fontweight="bold", color=DARK, pad=11)
    for spine in strip.spines.values():
        spine.set_visible(False)
    for boundary in np.arange(-0.5, c_data.shape[0] + 0.5, 1):
        strip.axhline(boundary, color="white", lw=0.30)
    strip.axvline(0.5, color="white", lw=0.35)
    return ax, strip


def save_individual(panel, pdf_path, png_path):
    fig = plt.figure(figsize=(90 / 25.4, 78 / 25.4), facecolor="white")
    grid = fig.add_gridspec(1, 1, left=0.12, right=0.97, top=0.97, bottom=0.12)
    if panel == "B":
        draw_b(fig, grid[0, 0])
    else:
        draw_c(fig, grid[0, 0])
    fig.savefig(pdf_path, format="pdf", dpi=600, metadata={"Title": f"Figure 6{panel} visual reset"})
    fig.savefig(png_path, format="png", dpi=600)
    plt.close(fig)


save_individual("B", B_PDF, B_PNG)
save_individual("C", C_PDF, C_PNG)

combined = plt.figure(figsize=(180 / 25.4, 80 / 25.4), facecolor="white")
combined_grid = combined.add_gridspec(1, 2, left=0.045, right=0.985, top=0.97, bottom=0.12, wspace=0.18)
draw_b(combined, combined_grid[0, 0])
draw_c(combined, combined_grid[0, 1])
combined.savefig(BC_PDF, format="pdf", dpi=600, metadata={"Title": "Figure 6B-C visual reset preview"})
combined.savefig(BC_PNG, format="png", dpi=600)
plt.close(combined)

semantic_text = f"""# Figure 6B/C visual reset semantic readout

## Design references

- Reused the audited common-quantitative-axis, patient-as-replication-unit, nested-repeat and unequal-information-area principles from the existing Figure 6 design-reference registry.
- No new literature search was required because both approved geometries have defensible precedents: an effect-retention/shift landscape for B and a descriptive Bland-Altman-style agreement landscape for C.
- Neither panel uses forest, lollipop or dumbbell grammar.

## Figure 6B

- x position is the existing TumorPurity+nCount-adjusted spatial Spearman rho.
- y position is the deterministic display transformation delta rho = adjusted rho - raw rho.
- Blue indicates attenuation (delta rho < 0); coral indicates strengthening (delta rho > 0).
- The horizontal zero line means no adjustment-induced change. The vertical coral line is the frozen cohort median adjusted rho ({cohort_median:.3f}).
- Marginal rugs show the 22 observed x and y values without adding density estimation.
- Labels are reproducible and restricted to four predeclared rules: largest positive delta, largest negative delta, maximum adjusted rho and minimum adjusted rho. The labelled patients are {', '.join(b_data.loc[b_data['display_label'].ne(''), 'patient_id'])}.

## Figure 6C

- x position is pair_mean_rho = (Section 1 rho + Section 2 rho) / 2.
- y position is pair_difference_rho = Section 2 rho - Section 1 rho.
- Coral diamonds indicate Section 2 stronger; blue diamonds indicate Section 1 stronger.
- The horizontal zero line is exact section agreement. The secondary dashed line is the descriptive median pair difference ({median_pair_difference:.3f}); no limits of agreement or new inferential test were added.
- The narrow two-column strip displays the two original positive rho values on one shared 0-to-0.67 scale. It is subordinate and contains no printed cell values.
- Sections remain nested within {paired_n} patients and are not treated as independent biological samples.

## Claim boundary

The reset visualizes adjustment retention and repeat-section agreement using deterministic transformations of frozen values. It does not recompute scores, correlations or P values, and it does not add new biological inference.

Standalone spatial robustness panels; scientific interpretation follows the source results.
"""
SEMANTIC.write_text(semantic_text, encoding="utf-8")

qc_text = f"""# Figure 6B/C visual reset QC

## Numerical non-modification

- Source patient rows: {b_data.shape[0]}/22
- Source paired-patient rows: {c_data.shape[0]}/21
- Adjusted rho positive: {int(b_data['rho_adjusted'].gt(0).sum())}/22
- Both sections positive: {int((c_data['section_1_rho'].gt(0) & c_data['section_2_rho'].gt(0)).sum())}/21
- Frozen median adjusted rho: {cohort_median:.15f}
- Frozen cohort P: {cohort_p:.15g}
- B adjusted rho range: {b_data['rho_adjusted'].min():.15f} to {b_data['rho_adjusted'].max():.15f}
- B delta rho range: {b_data['delta_rho'].min():.15f} to {b_data['delta_rho'].max():.15f}
- C pair mean range: {c_data['pair_mean_rho'].min():.15f} to {c_data['pair_mean_rho'].max():.15f}
- C pair difference range: {c_data['pair_difference_rho'].min():.15f} to {c_data['pair_difference_rho'].max():.15f}
- C median pair difference: {median_pair_difference:.15f}
- New correlation, P value, regression or limits of agreement: NO
- Figure 6A modified: NO
- Figure 6D started: NO
- Source SHA-256: {sha256(SOURCE)}
- Figure 6B PDF SHA-256: {sha256(B_PDF)}
- Figure 6C PDF SHA-256: {sha256(C_PDF)}
- Combined PDF SHA-256: {sha256(BC_PDF)}

## Information-density success gate

- B: every adjusted patient positive visible within seconds: PENDING_VISUAL_REVIEW
- B: adjusted range, attenuation/strengthening and cohort median visible: PENDING_VISUAL_REVIEW
- C: both-section positivity visible: PENDING_VISUAL_REVIEW
- C: average coupling, section difference and discrepancy frequency visible: PENDING_VISUAL_REVIEW
- B and C use distinct visual grammars: PENDING_VISUAL_REVIEW
- Text clipping: PENDING_VISUAL_REVIEW
- Text overlap: PENDING_VISUAL_REVIEW
- True-size readability: PENDING_VISUAL_REVIEW

Standalone spatial robustness panels; scientific interpretation follows the source results.
"""
QC.write_text(qc_text, encoding="utf-8")

print(f"WROTE {B_PDF}")
print(f"WROTE {C_PDF}")
print(f"WROTE {BC_PDF}")
print(f"PATIENTS={b_data.shape[0]} PAIRED={c_data.shape[0]}")
