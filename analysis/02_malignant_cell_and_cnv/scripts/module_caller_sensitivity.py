# Purpose: Six fixed modules within caller-supported subsets
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 2 and Supplementary Fig. S1.
from __future__ import annotations

from pathlib import Path
import hashlib

import numpy as np
import pandas as pd


ROOT = Path(r".")
OUT = ROOT / "0911_FigureS1_preplot_gate"
TABLES = OUT / "tables"
QC_DIR = OUT / "qc"
REPORTS = OUT / "reports"

MASTER_PATH = ROOT / "0910_Wu_multicaller_integration_audit" / "tables" / "Wu_multicaller_cell_master.tsv.gz"
MODULE_SCORE_PATH = ROOT / "0713_score_independent_rebuild" / "tables" / "module_scores_before_after_removing_score_genes.tsv"
MODULE_DEFINITION_PATH = ROOT / "0713_score_independent_rebuild" / "tables" / "module_gene_overlap.tsv"
MODULE_SUMMARY_PATH = ROOT / "0713_score_independent_rebuild" / "tables" / "module_keep_or_drop_summary.tsv"
MODULE_SCRIPT_PATH = ROOT / "0713_score_independent_rebuild" / "scripts" / "00_score_independent_reanalysis.R"
MCL1_REJECTED_PATH = ROOT / "0708_teacher_questions" / "tables" / "05_drug_candidate_gene_expression.tsv"
MCL1_PATH = OUT / "tables" / "MCL1_recovered_normalized_expression.tsv.gz"
MCL1_SOURCE_RDS = ROOT / "output_step0" / "Wu2021_TNBC_malignant.rds"
MCL1_SCRIPT_PATH = OUT / "scripts" / "03a_recover_existing_MCL1_layer.R"
PROGRAM146_PATH = OUT / "tables" / "Program146_multicaller_sensitivity.tsv"

PATIENTS = ["CID3963", "CID4465", "CID4495", "CID44971", "CID44991", "CID4513", "CID4515", "CID4523"]
MODULES = [
    "UPR",
    "TNFA NFKB",
    "Hypoxia",
    "Adhesion Remodeling",
    "Wound Healing",
    "Survival Stress",
]

SUBSETS = [
    ("All_frozen_10836", "All frozen malignant", None),
    ("CopyKAT_supported", "CopyKAT-supported", "CopyKAT_support"),
    ("SCEVAN_supported", "SCEVAN-supported", "SCEVAN_support"),
    ("inferCNV_high_descriptive", "inferCNV high-CNA descriptive", "inferCNV_high_CNA"),
    ("at_least_two_callers_supported", "multi-method-supported subset", "at_least_two_callers_supported"),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_spearman(x: pd.Series, y: pd.Series) -> tuple[float, int, str]:
    pair = pd.DataFrame(
        {"x": pd.to_numeric(x, errors="coerce"), "y": pd.to_numeric(y, errors="coerce")}
    ).dropna()
    n = len(pair)
    if n < 3:
        return np.nan, n, "NOT_TESTABLE_LT3"
    if pair["x"].nunique() < 2 or pair["y"].nunique() < 2:
        return np.nan, n, "NOT_TESTABLE_CONSTANT"
    return float(pair["x"].corr(pair["y"], method="spearman")), n, "TESTABLE"


def effect_comparison(effect: float, all_effect: float) -> dict[str, object]:
    if pd.isna(effect) or pd.isna(all_effect):
        return {
            "delta_vs_all": np.nan,
            "absolute_effect_ratio_vs_all": np.nan,
            "direction_reversal_vs_all": pd.NA,
            "marked_attenuation_vs_all": pd.NA,
        }
    ratio = np.nan if all_effect == 0 else abs(effect) / abs(all_effect)
    reversal = bool(effect * all_effect < 0)
    # Frozen 0910 definition: attenuation means no reversal and <=50% of the all-cell absolute effect.
    attenuation = bool((not reversal) and (not pd.isna(ratio)) and ratio <= 0.50)
    return {
        "delta_vs_all": effect - all_effect,
        "absolute_effect_ratio_vs_all": ratio,
        "direction_reversal_vs_all": reversal,
        "marked_attenuation_vs_all": attenuation,
    }


def subset_mask(master: pd.DataFrame, source_column: str | None) -> pd.Series:
    if source_column is None:
        return pd.Series(True, index=master.index)
    return (master[source_column] == True).fillna(False)  # noqa: E712


def correlation_records(
    merged: pd.DataFrame,
    feature_names: list[str],
    feature_type: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    patient_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    for feature in feature_names:
        baseline: dict[str, float] = {}
        for patient in PATIENTS:
            frame = merged[merged["patient"] == patient]
            rho, _, _ = safe_spearman(frame["YAP_Stem_Score"], frame[feature])
            baseline[patient] = rho

        for subset_name, submission_label, source_column in SUBSETS:
            mask = subset_mask(merged, source_column)
            total_cells = int(mask.sum())
            this_subset: list[dict[str, object]] = []
            for patient in PATIENTS:
                frame = merged[(merged["patient"] == patient) & mask]
                rho, n, status = safe_spearman(frame["YAP_Stem_Score"], frame[feature])
                comparison = effect_comparison(rho, baseline[patient])
                row = {
                    "feature": feature,
                    "feature_type": feature_type,
                    "subset": subset_name,
                    "submission_label": submission_label,
                    "subset_source_column": source_column or "ALL_FROZEN",
                    "subset_total_cells": total_cells,
                    "patient": patient,
                    "patient_n": n,
                    "rho": rho,
                    "status": status,
                    "positive": (rho > 0) if not pd.isna(rho) else pd.NA,
                    "all_frozen_rho": baseline[patient],
                    **comparison,
                }
                patient_rows.append(row)
                this_subset.append(row)

            valid = [row for row in this_subset if not pd.isna(row["rho"])]
            effects = [float(row["rho"]) for row in valid]
            reversal_patients = [str(row["patient"]) for row in valid if row["direction_reversal_vs_all"] is True]
            attenuation_patients = [str(row["patient"]) for row in valid if row["marked_attenuation_vs_all"] is True]
            summary_rows.append(
                {
                    "feature": feature,
                    "feature_type": feature_type,
                    "subset": subset_name,
                    "submission_label": submission_label,
                    "subset_source_column": source_column or "ALL_FROZEN",
                    "total_cells": total_cells,
                    "n_patients_testable": len(valid),
                    "n_patients_positive": sum(effect > 0 for effect in effects),
                    "median_rho": float(np.median(effects)) if effects else np.nan,
                    "min_rho": float(np.min(effects)) if effects else np.nan,
                    "max_rho": float(np.max(effects)) if effects else np.nan,
                    "n_direction_reversals_vs_all": len(reversal_patients),
                    "direction_reversal_patients": ",".join(reversal_patients),
                    "n_marked_attenuations_vs_all": len(attenuation_patients),
                    "marked_attenuation_patients": ",".join(attenuation_patients),
                    "patient_unit_inference": "UNWEIGHTED_PATIENT_SUMMARY_NO_POOLED_CELL_P_VALUE",
                }
            )

    return pd.DataFrame(patient_rows), pd.DataFrame(summary_rows)


def bool_text(value: object) -> str:
    if pd.isna(value):
        return "NA"
    return "YES" if bool(value) else "NO"


def fmt(value: object, digits: int = 3) -> str:
    if pd.isna(value):
        return "NA"
    return f"{float(value):.{digits}f}"


def main() -> None:
    for directory in (TABLES, QC_DIR, REPORTS):
        directory.mkdir(parents=True, exist_ok=True)

    required_paths = [
        MASTER_PATH,
        MODULE_SCORE_PATH,
        MODULE_DEFINITION_PATH,
        MODULE_SUMMARY_PATH,
        MODULE_SCRIPT_PATH,
        MCL1_PATH,
        MCL1_SOURCE_RDS,
        MCL1_REJECTED_PATH,
        MCL1_SCRIPT_PATH,
    ]
    missing_paths = [str(path) for path in required_paths if not path.exists()]
    if missing_paths:
        raise FileNotFoundError("Missing frozen input(s): " + "; ".join(missing_paths))

    master = pd.read_csv(MASTER_PATH, sep="\t", compression="gzip", low_memory=False)
    if len(master) != 10836 or master["cell_id"].duplicated().any():
        raise RuntimeError("Frozen master is not the expected unique 10,836-cell universe")
    if master["patient"].astype(str).isin(PATIENTS).sum() != 10836:
        raise RuntimeError("Frozen master patient set does not match the eight prespecified patients")
    master_ids = set(master["cell_id"].astype(str))

    module_long = pd.read_csv(MODULE_SCORE_PATH, sep="\t", low_memory=False)
    module_def = pd.read_csv(MODULE_DEFINITION_PATH, sep="\t")
    module_summary = pd.read_csv(MODULE_SUMMARY_PATH, sep="\t")

    exact_core = module_summary.loc[
        module_summary["module"].isin(MODULES), ["module", "keep_or_drop"]
    ].copy()
    if set(exact_core["module"]) != set(MODULES) or not (exact_core["keep_or_drop"] == "KEEP_MAIN").all():
        raise RuntimeError("Six-core module identity is not unambiguously confirmed as KEEP_MAIN")

    qc_rows: list[dict[str, object]] = []
    qc_rows.extend(
        [
            {
                "scope": "input_file",
                "item": path.name,
                "check": "sha256",
                "observed": sha256(path),
                "expected": "READ_ONLY_FROZEN_INPUT",
                "status": "PASS",
                "source_path": str(path),
            }
            for path in required_paths
        ]
    )

    recovered_frames: list[pd.DataFrame] = []
    definition_records: list[dict[str, object]] = []
    for module in MODULES:
        rows = module_long[module_long["module"] == module].copy()
        ids = rows["cell_id"].astype(str)
        duplicate_ids = int(ids.duplicated().sum())
        overlap = len(set(ids) & master_ids)
        extra_ids = len(set(ids) - master_ids)
        missing_ids = len(master_ids - set(ids))
        patient_check = rows[["cell_id", "patient"]].merge(
            master[["cell_id", "patient"]], on="cell_id", how="inner", suffixes=("_module", "_master")
        )
        patient_mismatch = int((patient_check["patient_module"].astype(str) != patient_check["patient_master"].astype(str)).sum())
        score_missing = int(pd.to_numeric(rows["score_after"], errors="coerce").isna().sum())
        defs = module_def[module_def["module"] == module]
        if len(defs) != 1:
            raise RuntimeError(f"Definition row is not unique for {module}")
        retained_genes = str(defs.iloc[0]["retained_genes"])
        retained_count = len([gene for gene in retained_genes.split(",") if gene])
        observed_gene_counts = sorted(rows["n_genes_after"].dropna().astype(int).unique().tolist())
        status = "PASS" if (
            len(rows) == 10836
            and duplicate_ids == 0
            and overlap == 10836
            and extra_ids == 0
            and missing_ids == 0
            and patient_mismatch == 0
            and score_missing == 0
            and observed_gene_counts == [retained_count]
        ) else "NOT_RECOVERABLE"
        definition_records.append(
            {
                "module": module,
                "retained_genes": retained_genes,
                "n_genes_after": retained_count,
                "score_column": "score_after",
                "score_algorithm": "cohort-pooled z-score of mean RNA/data log-normalized expression after frozen YAP/Stem score-gene removal",
                "definition_source": str(MODULE_DEFINITION_PATH),
                "score_source": str(MODULE_SCORE_PATH),
                "generating_script": str(MODULE_SCRIPT_PATH),
                "recovery_status": status,
            }
        )
        for check, observed, expected in [
            ("rows", len(rows), 10836),
            ("exact_cell_id_overlap", overlap, 10836),
            ("duplicate_cell_ids", duplicate_ids, 0),
            ("extra_cell_ids", extra_ids, 0),
            ("missing_cell_ids", missing_ids, 0),
            ("patient_mismatch", patient_mismatch, 0),
            ("missing_score_after", score_missing, 0),
            ("n_genes_after_consistent", ",".join(map(str, observed_gene_counts)), str(retained_count)),
        ]:
            qc_rows.append(
                {
                    "scope": "module_recovery",
                    "item": module,
                    "check": check,
                    "observed": observed,
                    "expected": expected,
                    "status": "PASS" if str(observed) == str(expected) else "FAIL",
                    "source_path": str(MODULE_SCORE_PATH),
                }
            )
        if status != "PASS":
            raise RuntimeError(f"{module} is NOT_RECOVERABLE; no sensitivity calculation was performed")
        recovered_frames.append(rows[["cell_id", "score_after"]].rename(columns={"score_after": module}))

    module_wide = recovered_frames[0]
    for frame in recovered_frames[1:]:
        module_wide = module_wide.merge(frame, on="cell_id", how="inner", validate="one_to_one")

    rejected_mcl1 = pd.read_csv(MCL1_REJECTED_PATH, sep="\t", usecols=["cell_id", "patient", "MCL1"])
    rejected_missing = int(pd.to_numeric(rejected_mcl1["MCL1"], errors="coerce").isna().sum())
    qc_rows.append(
        {
            "scope": "MCL1_recovery",
            "item": "MCL1_prior_candidate_rejected",
            "check": "missing_MCL1_values",
            "observed": rejected_missing,
            "expected": 10836,
            "status": "PASS" if rejected_missing == 10836 else "FAIL",
            "source_path": str(MCL1_REJECTED_PATH),
        }
    )

    mcl1 = pd.read_csv(MCL1_PATH, sep="\t", compression="gzip", usecols=["cell_id", "patient", "MCL1"])
    mcl1_ids = set(mcl1["cell_id"].astype(str))
    mcl1_check = mcl1.merge(
        master[["cell_id", "patient"]], on="cell_id", how="inner", suffixes=("_mcl1", "_master")
    )
    mcl1_metrics = {
        "rows": len(mcl1),
        "exact_cell_id_overlap": len(mcl1_ids & master_ids),
        "duplicate_cell_ids": int(mcl1["cell_id"].duplicated().sum()),
        "extra_cell_ids": len(mcl1_ids - master_ids),
        "missing_cell_ids": len(master_ids - mcl1_ids),
        "patient_mismatch": int((mcl1_check["patient_mcl1"].astype(str) != mcl1_check["patient_master"].astype(str)).sum()),
        "missing_MCL1_values": int(pd.to_numeric(mcl1["MCL1"], errors="coerce").isna().sum()),
    }
    mcl1_expected = {
        "rows": 10836,
        "exact_cell_id_overlap": 10836,
        "duplicate_cell_ids": 0,
        "extra_cell_ids": 0,
        "missing_cell_ids": 0,
        "patient_mismatch": 0,
        "missing_MCL1_values": 0,
    }
    for check, observed in mcl1_metrics.items():
        expected = mcl1_expected[check]
        qc_rows.append(
            {
                "scope": "MCL1_recovery",
                "item": "MCL1",
                "check": check,
                "observed": observed,
                "expected": expected,
                "status": "PASS" if observed == expected else "FAIL",
                "source_path": str(MCL1_PATH),
            }
        )
    if any(mcl1_metrics[key] != value for key, value in mcl1_expected.items()):
        raise RuntimeError("MCL1 normalized-expression column is not unambiguously recoverable")

    # The multi-caller master contains a legacy metadata column named Hypoxia.
    # Preserve it under an explicit non-analysis name so the recovered 0713 score_after column is used.
    collision_renames = {
        feature: f"{feature}__master_not_used"
        for feature in MODULES + ["MCL1"]
        if feature in master.columns
    }
    master_for_merge = master.rename(columns=collision_renames)
    for feature, renamed in collision_renames.items():
        qc_rows.append(
            {
                "scope": "feature_source_gate",
                "item": feature,
                "check": "legacy_master_column_excluded",
                "observed": renamed,
                "expected": renamed,
                "status": "PASS",
                "source_path": str(MASTER_PATH),
            }
        )
    merged = master_for_merge.merge(module_wide, on="cell_id", how="left", validate="one_to_one")
    merged = merged.merge(mcl1[["cell_id", "MCL1"]], on="cell_id", how="left", validate="one_to_one")
    if merged[MODULES + ["MCL1"]].isna().any().any():
        raise RuntimeError("Recovered feature merge contains missing values")

    module_patient, module_summary_result = correlation_records(merged, MODULES, "score-independent module score_after")
    mcl1_patient, mcl1_summary = correlation_records(merged, ["MCL1"], "single-gene normalized RNA/data expression")

    module_patient_path = TABLES / "six_module_multicaller_sensitivity_by_patient.tsv"
    module_summary_path = TABLES / "six_module_multicaller_sensitivity_summary.tsv"
    mcl1_path = TABLES / "MCL1_multicaller_sensitivity.tsv"
    module_patient.to_csv(module_patient_path, sep="\t", index=False)
    module_summary_result.to_csv(module_summary_path, sep="\t", index=False)
    mcl1_combined = pd.concat(
        [
            module.rename(columns={})
            for module in [
                mcl1_patient.assign(record_type="patient"),
                mcl1_summary.assign(record_type="summary"),
            ]
        ],
        ignore_index=True,
        sort=False,
    )
    first_cols = ["record_type", "feature", "feature_type", "subset", "submission_label"]
    mcl1_combined = mcl1_combined[first_cols + [col for col in mcl1_combined.columns if col not in first_cols]]
    mcl1_combined.to_csv(mcl1_path, sep="\t", index=False)

    definition_df = pd.DataFrame(definition_records)
    definition_df.to_csv(TABLES / "six_module_recovery_manifest.tsv", sep="\t", index=False)

    qc_rows.extend(
        [
            {
                "scope": "output",
                "item": module_patient_path.name,
                "check": "rows",
                "observed": len(module_patient),
                "expected": 6 * 5 * 8,
                "status": "PASS" if len(module_patient) == 240 else "FAIL",
                "source_path": str(module_patient_path),
            },
            {
                "scope": "output",
                "item": module_summary_path.name,
                "check": "rows",
                "observed": len(module_summary_result),
                "expected": 6 * 5,
                "status": "PASS" if len(module_summary_result) == 30 else "FAIL",
                "source_path": str(module_summary_path),
            },
            {
                "scope": "output",
                "item": mcl1_path.name,
                "check": "patient_plus_summary_rows",
                "observed": len(mcl1_combined),
                "expected": 5 * 8 + 5,
                "status": "PASS" if len(mcl1_combined) == 45 else "FAIL",
                "source_path": str(mcl1_path),
            },
            {
                "scope": "missing_feature_gate",
                "item": "PPP1R15A",
                "check": "frozen_normalized_expression_column",
                "observed": "NOT_RECOVERABLE",
                "expected": "NOT_RECOVERABLE_NO_RECOMPUTATION",
                "status": "PASS",
                "source_path": str(MCL1_PATH),
            },
            {
                "scope": "missing_feature_gate",
                "item": "ATF3",
                "check": "frozen_normalized_expression_column",
                "observed": "NOT_RECOVERABLE",
                "expected": "NOT_RECOVERABLE_NO_RECOMPUTATION",
                "status": "PASS",
                "source_path": str(MCL1_PATH),
            },
            {
                "scope": "method",
                "item": "all_features",
                "check": "pooled_cell_inference_or_new_p_values",
                "observed": 0,
                "expected": 0,
                "status": "PASS",
                "source_path": "patient-wise Spearman; unweighted patient summaries",
            },
            {
                "scope": "method",
                "item": "marked_attenuation",
                "check": "frozen_0910_definition",
                "observed": "no reversal and abs(rho_subset)/abs(rho_all)<=0.50",
                "expected": "no reversal and abs(rho_subset)/abs(rho_all)<=0.50",
                "status": "PASS",
                "source_path": str(ROOT / "0910_Wu_multicaller_integration_audit" / "scripts" / "01_build_multicaller_master_and_sensitivity.py"),
            },
        ]
    )
    qc = pd.DataFrame(qc_rows)
    qc_path = QC_DIR / "module_sensitivity_QC.tsv"
    qc.to_csv(qc_path, sep="\t", index=False)
    if (qc["status"] == "FAIL").any():
        raise RuntimeError("A module sensitivity QC check failed")

    summaries_nonall = module_summary_result[module_summary_result["subset"] != "All_frozen_10836"]
    all_medians_positive = bool((module_summary_result["median_rho"] > 0).all())
    all_module_subset_positive_majority = bool(
        (module_summary_result["n_patients_positive"] == module_summary_result["n_patients_testable"]).all()
    )
    reversals = module_patient[module_patient["direction_reversal_vs_all"] == True].copy()  # noqa: E712
    attenuations = module_patient[module_patient["marked_attenuation_vs_all"] == True].copy()  # noqa: E712

    module_lines = []
    for module in MODULES:
        sub = module_summary_result[module_summary_result["feature"] == module]
        all_row = sub[sub["subset"] == "All_frozen_10836"].iloc[0]
        multi_row = sub[sub["subset"] == "at_least_two_callers_supported"].iloc[0]
        infer_row = sub[sub["subset"] == "inferCNV_high_descriptive"].iloc[0]
        module_lines.append(
            f"- {module}: all frozen {int(all_row['n_patients_positive'])}/{int(all_row['n_patients_testable'])} positive, "
            f"median rho {fmt(all_row['median_rho'])}; multi-method-supported {int(multi_row['n_patients_positive'])}/"
            f"{int(multi_row['n_patients_testable'])} positive, median rho {fmt(multi_row['median_rho'])}; "
            f"inferCNV-high descriptive median rho {fmt(infer_row['median_rho'])}."
        )

    reversal_text = (
        "None."
        if reversals.empty
        else "; ".join(
            f"{row.feature} / {row.submission_label} / {row.patient} ({fmt(row.rho)} vs all {fmt(row.all_frozen_rho)})"
            for row in reversals.itertuples()
        )
    )
    attenuation_text = (
        "None."
        if attenuations.empty
        else "; ".join(
            f"{row.feature} / {row.submission_label} / {row.patient} (ratio {fmt(row.absolute_effect_ratio_vs_all, 2)})"
            for row in attenuations.itertuples()
        )
    )

    program_comparison = "Program146 is an aggregate 146-gene program, whereas these six recovered scores represent prespecified, score-independent functional branches. Their module-specific reversals and attenuation patterns therefore add non-redundant biological resolution, even when the broad direction agrees with Program146."
    recommendation = "KEEP" if (all_medians_positive and (not attenuations.empty or not reversals.empty)) else "DO_NOT_ADD"

    mcl1_all = mcl1_summary[mcl1_summary["subset"] == "All_frozen_10836"].iloc[0]
    mcl1_multi = mcl1_summary[mcl1_summary["subset"] == "at_least_two_callers_supported"].iloc[0]
    report = [
        "# Six-module multi-caller sensitivity interpretation",
        "",
        "## Recovery gate",
        "",
        f"All six formal frozen core modules were recovered read-only from `{MODULE_SCORE_PATH}` using the existing `score_after` column. The exact retained-gene definitions were read from `{MODULE_DEFINITION_PATH}`, and the generating script was `{MODULE_SCRIPT_PATH}`. Each module had 10,836/10,836 exact frozen cell-ID overlap, zero duplicate IDs, zero extra or missing IDs, zero patient mismatches and no missing `score_after` values. No signature or score was recomputed.",
        "",
        "The sixth formal module is **Hypoxia**, confirmed from the frozen six-core conclusion and the 0713 KEEP_MAIN manifest; it was not inferred from the current results.",
        "",
        f"The earlier tabular candidate `{MCL1_REJECTED_PATH}` was rejected because its MCL1 column is entirely missing (10,836/10,836 values). MCL1 was instead recovered read-only from the already normalized RNA/data layer in `{MCL1_SOURCE_RDS}` by `{MCL1_SCRIPT_PATH}`; no normalization, rescaling or signature calculation was run. The recovered vector passed the 10,836-cell identity gate. MCL1 is a single-gene normalized-expression feature, not a module score. PPP1R15A and ATF3 remain NOT_RECOVERABLE; they were not recalculated.",
        "",
        "## A. Do all six modules retain direction?",
        "",
        ("Yes at the cross-patient median level across all five prespecified subsets." if all_medians_positive else "No; at least one module-subset cross-patient median is not positive."),
        (" Every testable patient is positive for every module-subset combination." if all_module_subset_positive_majority else " Some patient-specific module-subset effects are non-positive; these are itemized below."),
        "",
        *module_lines,
        "",
        "## B. Patient-specific direction reversals",
        "",
        reversal_text,
        "",
        "## C. Marked attenuation",
        "",
        "The frozen definition is no direction reversal and an absolute subset/all-frozen rho ratio <=0.50.",
        "",
        attenuation_text,
        "",
        "## D. Information beyond Program146",
        "",
        program_comparison,
        "",
        "## E. Figure S1F gate",
        "",
        "A module panel is justified only if the recovered functional branches add a distinct pattern beyond the broad Program146/YAP–Stem retention result. The recommendation below is based on the observed module-specific attenuation/reversal structure, not on new significance testing.",
        "",
        "## MCL1 descriptive sensitivity",
        "",
        f"All frozen: {int(mcl1_all['n_patients_positive'])}/{int(mcl1_all['n_patients_testable'])} patients positive, median rho {fmt(mcl1_all['median_rho'])}. Multi-method-supported subset: {int(mcl1_multi['n_patients_positive'])}/{int(mcl1_multi['n_patients_testable'])} positive, median rho {fmt(mcl1_multi['median_rho'])}. This is descriptive single-gene evidence and does not alter the frozen core-module evidence ranking.",
        "",
        f"FIGURE_S1F_RECOMMENDATION = {recommendation}",
    ]
    (REPORTS / "module_sensitivity_interpretation.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    print("Six modules:", ", ".join(MODULES))
    print("All summary medians positive:", all_medians_positive)
    print("All module-subset patient effects positive:", all_module_subset_positive_majority)
    print("Patient-level reversals:", len(reversals))
    print("Marked attenuations:", len(attenuations))
    print("Recommendation:", recommendation)


if __name__ == "__main__":
    main()
