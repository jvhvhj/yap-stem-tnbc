# Purpose: Existing Program146 score recovery and sensitivity
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 2 and Supplementary Fig. S1.
from __future__ import annotations

import gzip
import hashlib
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(r".")
OUT = ROOT / "0911_FigureS1_preplot_gate"
TABLES = OUT / "tables"
REPORTS = OUT / "reports"

PATIENTS = [
    "CID3963",
    "CID4465",
    "CID4495",
    "CID44971",
    "CID44991",
    "CID4513",
    "CID4515",
    "CID4523",
]

MASTER_PATH = ROOT / "0910_Wu_multicaller_integration_audit" / "tables" / "Wu_multicaller_cell_master.tsv.gz"
PATIENT_SUMMARY_PATH = ROOT / "0910_Wu_multicaller_integration_audit" / "tables" / "Wu_multicaller_patient_summary.tsv"
PROGRAM_PATH = ROOT / "P6_0_pan_cancer_program_specificity_audit" / "03_adjusted_coupling" / "Wu_cell_scores_for_P6C.tsv.gz"
PROGRAM_GENERATOR = ROOT / "P6_0_pan_cancer_program_specificity_audit" / "scripts" / "02_prepare_wu_cell_scores.R"
PROGRAM_GENE_LIST = ROOT / "0728_sensitivity_and_cnv_closure" / "02_direction_gate_sensitivity" / "gate_6of8_primary_genes.txt"
PROGRAM_REPRO_QC = ROOT / "P6_0_pan_cancer_program_specificity_audit" / "03_adjusted_coupling" / "Wu_score_reproduction_QC.tsv"
MODULE_PATH = ROOT / "0713_score_independent_rebuild" / "tables" / "module_scores_before_after_removing_score_genes.tsv"
MODULE_GENERATOR = ROOT / "0713_score_independent_rebuild" / "scripts" / "00_score_independent_reanalysis.R"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def timestamp(path: Path, attr: str) -> str:
    value = getattr(path.stat(), attr)
    return datetime.fromtimestamp(value).astimezone().isoformat(timespec="seconds")


def patient_from_cell_id(series: pd.Series) -> pd.Series:
    return series.astype("string").str.extract(r"^(CID\d+)_", expand=False)


def id_audit(
    frame: pd.DataFrame,
    *,
    id_col: str,
    patient_col: str | None,
    master: pd.DataFrame,
) -> dict:
    ids = frame[id_col].astype("string")
    unique_ids = set(ids.dropna().astype(str))
    master_ids = set(master["cell_id"].astype(str))
    n_patient_mismatch = 0
    if patient_col and patient_col in frame.columns:
        probe = frame[[id_col, patient_col]].drop_duplicates(id_col).copy()
        expected = master[["cell_id", "patient"]].rename(columns={"cell_id": id_col, "patient": "expected_patient"})
        probe = probe.merge(expected, on=id_col, how="inner")
        n_patient_mismatch = int((probe[patient_col].astype(str) != probe["expected_patient"].astype(str)).sum())
    else:
        derived = patient_from_cell_id(ids)
        expected_map = master.set_index("cell_id")["patient"]
        expected = ids.map(expected_map)
        n_patient_mismatch = int(((derived != expected) & derived.notna() & expected.notna()).sum())
    return {
        "n_rows": int(len(frame)),
        "n_unique_cell_ids": int(len(unique_ids)),
        "exact_id_overlap_with_frozen_10836": int(len(unique_ids & master_ids)),
        "frozen_ids_missing_from_candidate": int(len(master_ids - unique_ids)),
        "extra_ids_not_in_frozen_10836": int(len(unique_ids - master_ids)),
        "duplicate_id_rows": int(ids.dropna().duplicated().sum()),
        "patient_mismatch_rows": n_patient_mismatch,
    }


def candidate_row(
    path: Path,
    frame: pd.DataFrame,
    *,
    artifact_role: str,
    id_col: str,
    patient_col: str | None,
    requested_columns: list[str],
    definition_status: str,
    provenance: str,
    master: pd.DataFrame,
) -> dict:
    audit = id_audit(frame, id_col=id_col, patient_col=patient_col, master=master)
    return {
        "artifact_role": artifact_role,
        "file_path": str(path),
        "created_time": timestamp(path, "st_ctime"),
        "modified_time": timestamp(path, "st_mtime"),
        "file_size_bytes": path.stat().st_size,
        "columns": "|".join(map(str, frame.columns)),
        "cell_id_field": id_col,
        "patient_field": patient_col or "",
        "requested_score_or_result_columns": "|".join(requested_columns),
        **audit,
        "definition_confirmation": definition_status,
        "provenance_evidence": provenance,
        "file_sha256": sha256(path),
    }


def safe_spearman(x: pd.Series, y: pd.Series) -> tuple[float, int, str]:
    pair = pd.DataFrame({"x": pd.to_numeric(x, errors="coerce"), "y": pd.to_numeric(y, errors="coerce")}).dropna()
    n = len(pair)
    if n < 3:
        return np.nan, n, "NOT_TESTABLE_LT3"
    if pair["x"].nunique() < 2 or pair["y"].nunique() < 2:
        return np.nan, n, "NOT_TESTABLE_CONSTANT"
    return float(spearmanr(pair["x"], pair["y"]).statistic), n, "TESTABLE"


def comparison(effect: float, all_effect: float) -> dict:
    if pd.isna(effect) or pd.isna(all_effect):
        return {
            "delta_vs_all": np.nan,
            "absolute_effect_ratio_vs_all": np.nan,
            "direction_reversal_vs_all": np.nan,
            "marked_attenuation_vs_all": np.nan,
        }
    ratio = np.nan if all_effect == 0 else abs(effect) / abs(all_effect)
    reversal = bool(effect * all_effect < 0)
    attenuation = bool((not reversal) and (not pd.isna(ratio)) and ratio <= 0.50)
    return {
        "delta_vs_all": effect - all_effect,
        "absolute_effect_ratio_vs_all": ratio,
        "direction_reversal_vs_all": reversal,
        "marked_attenuation_vs_all": attenuation,
    }


def parse_counts(text: str) -> dict[str, int]:
    out: dict[str, int] = {}
    for item in str(text).split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        out[key.strip()] = int(float(value.strip()))
    return out


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)

    master = pd.read_csv(MASTER_PATH, sep="\t", low_memory=False)
    if len(master) != 10836 or master["cell_id"].duplicated().any():
        raise RuntimeError("Frozen multi-caller master is not the expected unique 10,836-cell table")
    if list(dict.fromkeys(master["patient"].astype(str))) != [
        "CID4465", "CID4495", "CID44971", "CID44991", "CID4513", "CID4515", "CID4523", "CID3963"
    ]:
        # Row order is not used for plotting; patient identities still have to be exact.
        if set(master["patient"].astype(str)) != set(PATIENTS):
            raise RuntimeError("Frozen master patient set differs from the fixed eight-patient set")

    program = pd.read_csv(PROGRAM_PATH, sep="\t", low_memory=False)
    module = pd.read_csv(MODULE_PATH, sep="\t", low_memory=False)
    f2 = pd.read_csv(
        ROOT / "0727_Figure2_complete_patient_resolved_revision_and_CNV_extension_plan" / "plotting_data" / "Figure2_panelAB_cell_scores.tsv",
        sep="\t",
        low_memory=False,
    )
    adhesion = pd.read_csv(
        ROOT / "adhesion_stemness_link" / "tables" / "adhesion_stemness_module_scores_cell_metadata.csv",
        low_memory=False,
    )
    lineage = pd.read_csv(
        ROOT / "0723_Core_program_gene_closure_and_multicellular_context_gate" / "tables" / "wu_lineage_cell_scores.tsv.gz",
        sep="\t",
        low_memory=False,
    )
    gene_expression = pd.read_csv(
        ROOT / "0708_teacher_questions" / "tables" / "05_drug_candidate_gene_expression.tsv",
        sep="\t",
        low_memory=False,
    )

    program_cols = [c for c in program.columns if re.search(r"Program146|program_score|UPR|Hypoxia", c, re.I)]
    module_cols = sorted(module["module"].dropna().astype(str).unique().tolist())
    candidates = [
        candidate_row(
            PROGRAM_PATH,
            program,
            artifact_role="Wu per-cell frozen axis plus Program146 artifact",
            id_col="cell_id",
            patient_col="patient_id",
            requested_columns=program_cols,
            definition_status="CONFIRMED_FOR_PROGRAM146_ONLY",
            provenance=(
                f"Generated by {PROGRAM_GENERATOR}; reads {PROGRAM_GENE_LIST}; enforces 146/146 feature coverage; "
                f"Program146_score is arithmetic mean of frozen log-normalized RNA/data expression. "
                "generic_UPR and generic_hypoxia are adjustment covariates and are not the frozen score-independent modules."
            ),
            master=master,
        ),
        candidate_row(
            MODULE_PATH,
            module,
            artifact_role="Wu per-cell long-format before/after score-gene-removal module scores",
            id_col="cell_id",
            patient_col="patient",
            requested_columns=module_cols,
            definition_status="CONFIRMED_FOR_0713_FROZEN_SCORE_INDEPENDENT_MODULES",
            provenance=f"Generated by {MODULE_GENERATOR}; score_after is the frozen score-gene-removed value.",
            master=master,
        ),
        candidate_row(
            ROOT / "0727_Figure2_complete_patient_resolved_revision_and_CNV_extension_plan" / "plotting_data" / "Figure2_panelAB_cell_scores.tsv",
            f2,
            artifact_role="Wu per-cell Figure 2 frozen YAP/Stem/Joint-axis plotting data",
            id_col="cell_id",
            patient_col="patient",
            requested_columns=[c for c in f2.columns if "score" in c.lower() or "axis" in c.lower()],
            definition_status="CONFIRMED_FOR_YAP17_STEM21_JOINT_ONLY",
            provenance="Figure 2 frozen plotting data; contains no Program146 or requested core-module/gene column.",
            master=master,
        ),
        candidate_row(
            ROOT / "adhesion_stemness_link" / "tables" / "adhesion_stemness_module_scores_cell_metadata.csv",
            adhesion,
            artifact_role="Pre-0713 Wu per-cell adhesion/stemness exploratory module table",
            id_col="cell_id",
            patient_col="patient",
            requested_columns=[c for c in adhesion.columns if "score" in c.lower()],
            definition_status="NOT_CONFIRMED_AS_FROZEN_SCORE_INDEPENDENT_DEFINITION",
            provenance="Modified before the 0713 score-independent freeze; column names and definitions do not establish equivalence.",
            master=master,
        ),
        candidate_row(
            ROOT / "0723_Core_program_gene_closure_and_multicellular_context_gate" / "tables" / "wu_lineage_cell_scores.tsv.gz",
            lineage,
            artifact_role="Wu per-cell lineage program audit",
            id_col="cell_id",
            patient_col="patient",
            requested_columns=[c for c in lineage.columns if "lineage__" in c],
            definition_status="OUT_OF_SCOPE_LINEAGE_PROGRAMS",
            provenance="Lineage-score artifact; contains no Program146 or requested frozen six-core module columns.",
            master=master,
        ),
        candidate_row(
            ROOT / "0708_teacher_questions" / "tables" / "05_drug_candidate_gene_expression.tsv",
            gene_expression,
            artifact_role="Wu per-cell normalized expression for predefined candidate genes",
            id_col="cell_id",
            patient_col="patient",
            requested_columns=[c for c in ["MCL1", "PPP1R15A", "ATF3"] if c in gene_expression.columns],
            definition_status="CONFIRMED_FOR_MCL1_NORMALIZED_RNA_DATA_EXPRESSION_ONLY",
            provenance=(
                "Generated by 0708_teacher_questions/05_drug_exploration.R from RNA/data log-normalized expression; "
                "the script extracts the single MCL1 row without aggregation or rescaling. PPP1R15A and ATF3 are absent."
            ),
            master=master,
        ),
    ]

    # The RDS audit was read-only and performed with base R because the object is a data.frame.
    rds = ROOT / "0630" / "meta_from_robustness.rds"
    candidates.append(
        {
            "artifact_role": "Pre-0713 Wu per-cell robustness metadata RDS",
            "file_path": str(rds),
            "created_time": timestamp(rds, "st_ctime"),
            "modified_time": timestamp(rds, "st_mtime"),
            "file_size_bytes": rds.stat().st_size,
            "columns": "orig.ident|nCount_RNA|nFeature_RNA|X|percent.mito|subtype|celltype_subset|celltype_minor|celltype_major|HH_program_score|YAP_score|Stem_score|P53_PATHWAY_score|UPR_score|TNFA_NFKB_score|G2M_CHECKPOINT_score|HH_group|YAP_z|Stem_z|HH_add_score|HH_product_score|HH_cutoff_0.5|HH_cutoff_0.67|HH_cutoff_0.75|patient_id",
            "cell_id_field": "rownames",
            "patient_field": "patient_id",
            "requested_score_or_result_columns": "HH_program_score|UPR_score|TNFA_NFKB_score",
            "n_rows": 10836,
            "n_unique_cell_ids": 10836,
            "exact_id_overlap_with_frozen_10836": 10836,
            "frozen_ids_missing_from_candidate": 0,
            "extra_ids_not_in_frozen_10836": 0,
            "duplicate_id_rows": 0,
            "patient_mismatch_rows": 0,
            "definition_confirmation": "NOT_CONFIRMED_AS_PROGRAM146_OR_0713_SCORE_INDEPENDENT_MODULES",
            "provenance_evidence": "Modified 0630, before the 0713 module freeze and 0728 Program146 definition; names alone are insufficient.",
            "file_sha256": sha256(rds),
        }
    )

    candidate_df = pd.DataFrame(candidates)
    candidate_df.to_csv(TABLES / "frozen_score_recovery_candidate_audit.tsv", sep="\t", index=False)

    # Whole-project filename/path inventory. This is intentionally read-only and excludes installed package trees
    # from candidate adjudication while retaining their counts in the search summary.
    tokens = [
        "program146", "strict_v2", "program_score", "mcl1", "ppp1r15a", "atf3", "upr", "tnfa", "nfkb",
        "adhesion", "wound", "survival", "hypoxia", "module_scores", "cell_scores",
    ]
    permitted_suffixes = {".csv", ".tsv", ".txt", ".rds", ".rdata", ".rda", ".parquet", ".feather", ".gz"}
    inventory_rows = []
    # `candidate_paths.txt` is the result of a read-only `rg --files <project> | rg -i <keywords>`
    # whole-project search. Reusing it avoids a second slow traversal through bundled package trees.
    candidate_path_file = OUT / "tmp" / "candidate_paths.txt"
    candidate_paths = [Path(line.strip()) for line in candidate_path_file.read_text(encoding="utf-8-sig").splitlines() if line.strip()]
    examined = len(candidate_paths)
    for path in candidate_paths:
        if not path.is_file() or path.suffix.lower() not in permitted_suffixes:
            continue
        low = str(path).lower()
        hits = [token for token in tokens if token in low]
        if not hits:
            continue
        inventory_rows.append(
            {
                "file_path": str(path),
                "extension": path.suffix.lower(),
                "file_size_bytes": path.stat().st_size,
                "created_time": timestamp(path, "st_ctime"),
                "modified_time": timestamp(path, "st_mtime"),
                "path_keyword_hits": "|".join(hits),
                "adjudication_scope": (
                    "INSTALLED_PACKAGE_OR_EXTERNAL_CACHE_NOT_PROJECT_RESULT"
                    if any(part.lower() in {"python_env", "r_library", "00_provenance", "source_archive"} for part in path.parts)
                    else "PROJECT_ARTIFACT_NAME_MATCH"
                ),
            }
        )
    inventory = pd.DataFrame(inventory_rows).sort_values(["adjudication_scope", "file_path"])
    inventory.to_csv(TABLES / "project_score_artifact_search_inventory.tsv", sep="\t", index=False)

    # Program146 must pass source, gene-list, cell-ID, and patient-ID gates.
    expected_program_hash = "371b2d0b6c4c99aebe452746ae0e165a4425b0c32ea640a58f5e0230719584b2"
    program_gene_hash = sha256(PROGRAM_GENE_LIST)
    p_audit = id_audit(program, id_col="cell_id", patient_col="patient_id", master=master)
    reproduction_qc = pd.read_csv(PROGRAM_REPRO_QC, sep="\t")
    program_coverage_ok = bool(
        ((reproduction_qc["item"] == "Program146_coverage") & (pd.to_numeric(reproduction_qc["observed"], errors="coerce") == 146) & (reproduction_qc["status"] == "VERIFIED")).any()
    )
    program_recoverable = bool(
        program_gene_hash == expected_program_hash
        and program_coverage_ok
        and p_audit["exact_id_overlap_with_frozen_10836"] == 10836
        and p_audit["duplicate_id_rows"] == 0
        and p_audit["patient_mismatch_rows"] == 0
        and p_audit["frozen_ids_missing_from_candidate"] == 0
        and p_audit["extra_ids_not_in_frozen_10836"] == 0
        and "Program146_score" in program.columns
    )

    gate = pd.DataFrame(
        [
            {"check": "Program146 gene-list SHA256", "observed": program_gene_hash, "expected": expected_program_hash, "status": "PASS" if program_gene_hash == expected_program_hash else "FAIL"},
            {"check": "Program146 feature coverage", "observed": 146 if program_coverage_ok else "not verified", "expected": 146, "status": "PASS" if program_coverage_ok else "FAIL"},
            {"check": "Frozen exact cell-ID overlap", "observed": p_audit["exact_id_overlap_with_frozen_10836"], "expected": 10836, "status": "PASS" if p_audit["exact_id_overlap_with_frozen_10836"] == 10836 else "FAIL"},
            {"check": "Duplicate cell IDs", "observed": p_audit["duplicate_id_rows"], "expected": 0, "status": "PASS" if p_audit["duplicate_id_rows"] == 0 else "FAIL"},
            {"check": "Patient mismatches", "observed": p_audit["patient_mismatch_rows"], "expected": 0, "status": "PASS" if p_audit["patient_mismatch_rows"] == 0 else "FAIL"},
            {"check": "PROGRAM146_EXISTING_SCORE_RECOVERABLE", "observed": "YES" if program_recoverable else "NO", "expected": "YES only if all gates pass", "status": "PASS" if program_recoverable else "HOLD"},
        ]
    )
    gate.to_csv(TABLES / "Program146_recovery_gate.tsv", sep="\t", index=False)

    if not program_recoverable:
        raise RuntimeError("PROGRAM146_EXISTING_SCORE_RECOVERABLE = NO; no sensitivity table was generated")

    merged = master.merge(
        program[["cell_id", "patient_id", "Joint_axis", "Program146_score"]],
        on="cell_id",
        how="left",
        validate="one_to_one",
    )
    if merged["Program146_score"].isna().any():
        raise RuntimeError("Recovered Program146 has missing values after exact-ID merge")
    if (merged["patient"].astype(str) != merged["patient_id"].astype(str)).any():
        raise RuntimeError("Patient mismatch after Program146 merge")
    max_joint_diff = float(np.max(np.abs(merged["YAP_Stem_Score"] - merged["Joint_axis"])))

    subset_specs = [
        ("All_frozen_10836", "All frozen malignant epithelial cells", pd.Series(True, index=merged.index)),
        ("CopyKAT_supported", "CopyKAT-supported subset", merged["CopyKAT_support"] == True),
        ("SCEVAN_supported", "SCEVAN-supported subset", merged["SCEVAN_support"] == True),
        ("inferCNV_high_descriptive", "inferCNV high-CNA descriptive subset", merged["inferCNV_high_CNA"] == True),
        ("at_least_two_callers_supported", "multi-method-supported subset", merged["at_least_two_callers_supported"] == True),
    ]
    subset_specs = [(sid, label, mask.fillna(False)) for sid, label, mask in subset_specs]

    all_effect = {}
    for patient in PATIENTS:
        frame = merged[merged["patient"] == patient]
        all_effect[patient] = safe_spearman(frame["YAP_Stem_Score"], frame["Program146_score"])[0]

    rows = []
    for subset_id, display_label, mask in subset_specs:
        subset_rows = []
        for patient in PATIENTS:
            patient_all = merged[merged["patient"] == patient]
            frame = merged[(merged["patient"] == patient) & mask]
            rho, n_pair, status = safe_spearman(frame["YAP_Stem_Score"], frame["Program146_score"])
            row = {
                "record_type": "patient",
                "feature": "Program146",
                "frozen_column": "Program146_score",
                "source_subset_field": subset_id,
                "submission_facing_subset_label": display_label,
                "patient": patient,
                "n_frozen_patient_cells": len(patient_all),
                "n_subset_cells": len(frame),
                "n_complete_pairs": n_pair,
                "effect_type": "within_patient_Spearman_frozen_Joint_axis_vs_recovered_Program146",
                "effect": rho,
                "status": status,
                "positive": bool(rho > 0) if not pd.isna(rho) else np.nan,
                **comparison(rho, all_effect[patient]),
                "n_patients_testable": np.nan,
                "n_patients_positive": np.nan,
                "median_effect": np.nan,
                "min_effect": np.nan,
                "max_effect": np.nan,
                "n_direction_reversals": np.nan,
                "n_marked_attenuations": np.nan,
            }
            rows.append(row)
            subset_rows.append(row)
        valid = [row["effect"] for row in subset_rows if not pd.isna(row["effect"])]
        rows.append(
            {
                "record_type": "cross_patient_summary",
                "feature": "Program146",
                "frozen_column": "Program146_score",
                "source_subset_field": subset_id,
                "submission_facing_subset_label": display_label,
                "patient": "ALL_8_PATIENTS",
                "n_frozen_patient_cells": len(merged),
                "n_subset_cells": int(mask.sum()),
                "n_complete_pairs": int(sum(row["n_complete_pairs"] for row in subset_rows)),
                "effect_type": "patient_is_statistical_unit; unweighted_summary_of_patient_effects",
                "effect": np.nan,
                "status": "SUMMARY",
                "positive": np.nan,
                "delta_vs_all": np.nan,
                "absolute_effect_ratio_vs_all": np.nan,
                "direction_reversal_vs_all": np.nan,
                "marked_attenuation_vs_all": np.nan,
                "n_patients_testable": len(valid),
                "n_patients_positive": int(sum(x > 0 for x in valid)),
                "median_effect": float(np.median(valid)),
                "min_effect": float(np.min(valid)),
                "max_effect": float(np.max(valid)),
                "n_direction_reversals": int(sum(row["direction_reversal_vs_all"] is True for row in subset_rows)),
                "n_marked_attenuations": int(sum(row["marked_attenuation_vs_all"] is True for row in subset_rows)),
            }
        )
    sensitivity = pd.DataFrame(rows)
    sensitivity.to_csv(TABLES / "Program146_multicaller_sensitivity.tsv", sep="\t", index=False, na_rep="NA")

    # Prepare S1A source data strictly from the frozen patient summary.
    summary = pd.read_csv(PATIENT_SUMMARY_PATH, sep="\t")
    caller_map = {
        "CopyKAT": ("CopyKAT", [("aneuploid", "support", "support"), ("diploid", "non-support", "non_support"), ("not.defined", "not evaluable", "not_evaluable")]),
        "SCEVAN": ("SCEVAN", [("tumor", "support", "support"), ("normal", "non-support", "non_support"), ("filtered", "not evaluable", "not_evaluable")]),
        "inferCNV-high-descriptive": ("inferCNV", [("True", "high CNA burden (descriptive)", "infer_high"), ("False", "below threshold", "infer_below")]),
    }
    hmm = summary[summary["caller"] == "inferCNV-HMM-native-status"].set_index("patient")
    s1_rows = []
    for patient_index, patient in enumerate(PATIENTS, start=1):
        for method_index, (source_caller, (method_label, categories)) in enumerate(caller_map.items(), start=1):
            hit = summary[(summary["patient"] == patient) & (summary["caller"] == source_caller)]
            if len(hit) != 1:
                raise RuntimeError(f"Expected one {source_caller} summary row for {patient}")
            source = hit.iloc[0]
            counts = parse_counts(source["raw_class_counts"])
            method_total = 0
            for category_index, (source_category, display_category, semantic_class) in enumerate(categories, start=1):
                count = int(counts.get(source_category, 0))
                method_total += count
                s1_rows.append(
                    {
                        "patient_order": patient_index,
                        "patient": patient,
                        "method_order": method_index,
                        "method": method_label,
                        "source_caller_field": source_caller,
                        "category_order": category_index,
                        "source_category": source_category,
                        "display_category": display_category,
                        "semantic_class": semantic_class,
                        "count": count,
                        "fraction_of_frozen_patient_cells": count / int(source["n_frozen_cells"]),
                        "n_frozen_cells": int(source["n_frozen_cells"]),
                        "n_evaluable_source": int(source["n_evaluable"]),
                        "n_supported_source": source["n_supported"],
                        "n_not_evaluable_source": int(source["n_NA"]),
                        "source_native_status": source["native_status"],
                        "source_interpretation": source["interpretation"],
                        "inferCNV_HMM_patient_status": hmm.loc[patient, "native_status"],
                        "inferCNV_HMM_patient_reason": hmm.loc[patient, "interpretation"],
                        "source_file": str(PATIENT_SUMMARY_PATH),
                    }
                )
            if method_total != int(source["n_frozen_cells"]):
                raise RuntimeError(f"Composition does not sum to patient total for {patient}/{method_label}")
    s1a = pd.DataFrame(s1_rows)
    s1a.to_csv(TABLES / "FigureS1A_source_data.tsv", sep="\t", index=False)

    feature_status = pd.DataFrame(
        [
            ["Program146", "YES", "Program146_score", str(PROGRAM_PATH), "Exact frozen 146-gene definition and 10,836/10,836 cell-ID map verified."],
            ["UPR", "YES", "score_after", str(MODULE_PATH), "Recoverable from long-format module=UPR; exact 0713 score-independent definition."],
            ["TNFA-NFKB", "YES", "score_after", str(MODULE_PATH), "Recoverable from long-format module=TNFA-NFKB; exact 0713 score-independent definition."],
            ["Hypoxia", "YES", "score_after", str(MODULE_PATH), "Recoverable from long-format module=Hypoxia; excluded from Figure S1 by instruction."],
            ["Adhesion Remodeling", "YES", "score_after", str(MODULE_PATH), "Recoverable from exact 0713 long-format module artifact."],
            ["Wound Healing", "YES", "score_after", str(MODULE_PATH), "Recoverable from exact 0713 long-format module artifact."],
            ["Survival Stress", "YES", "score_after", str(MODULE_PATH), "Recoverable from exact 0713 long-format module artifact."],
            ["MCL1", "YES_EXPRESSION_ONLY", "MCL1", str(ROOT / "0708_teacher_questions" / "tables" / "05_drug_candidate_gene_expression.tsv"), "Exact-ID-mappable RNA/data log-normalized single-gene expression artifact; not a module score and not used for Program146 sensitivity."],
            ["PPP1R15A", "NO", "", "", "No unambiguous precomputed frozen per-cell result column located; raw expression is not substituted."],
            ["ATF3", "NO", "", "", "No unambiguous precomputed frozen per-cell result column located; raw expression is not substituted."],
        ],
        columns=["feature", "existing_frozen_per_cell_result_recoverable", "value_field", "source_path", "adjudication"],
    )
    feature_status.to_csv(TABLES / "frozen_score_feature_recovery_status.tsv", sep="\t", index=False)

    summary_rows = sensitivity[sensitivity["record_type"] == "cross_patient_summary"].copy()
    lines = [
        "# Frozen-score recovery audit",
        "",
        "## Gate decision",
        "",
        f"`PROGRAM146_EXISTING_SCORE_RECOVERABLE = {'YES' if program_recoverable else 'NO'}`",
        "",
        "This is a recovery and identity-mapping audit only. No biological score, gene expression value, caller output, master-table field, threshold, or signature was recomputed or overwritten.",
        "",
        "## Why Program146 passes",
        "",
        f"- Candidate: `{PROGRAM_PATH}`.",
        f"- Generator: `{PROGRAM_GENERATOR}`.",
        f"- Frozen gene list: `{PROGRAM_GENE_LIST}`; SHA256 `{program_gene_hash}` (expected frozen hash matched).",
        "- The generator explicitly enforces 146/146 gene coverage and computes the arithmetic mean from the frozen RNA/data layer; the saved reproduction QC reports Program146 coverage as VERIFIED.",
        f"- Exact cell-ID overlap: {p_audit['exact_id_overlap_with_frozen_10836']}/10,836; duplicate IDs: {p_audit['duplicate_id_rows']}; patient mismatches: {p_audit['patient_mismatch_rows']}; extra IDs: {p_audit['extra_ids_not_in_frozen_10836']}.",
        f"- Frozen Joint-axis reproduction check after ID merge: maximum absolute difference {max_joint_diff:.3g}.",
        "",
        "## Important non-equivalences",
        "",
        "- `HH_program_score` in `0630/meta_from_robustness.rds` predates the 0728 strict_v2 Program146 definition and is not treated as Program146.",
        "- `generic_UPR` and `generic_hypoxia` in the P6 artifact are adjustment covariates, not the frozen 0713 score-independent UPR/Hypoxia modules.",
        "- The exact score-independent modules are recoverable from `0713_score_independent_rebuild/tables/module_scores_before_after_removing_score_genes.tsv` using `score_after` and the frozen module label. They were not used to add panels in this gate.",
        "- MCL1 is recoverable only as an existing RNA/data log-normalized single-gene expression column in `0708_teacher_questions/tables/05_drug_candidate_gene_expression.tsv`; its generator and exact 10,836-cell mapping are auditable. It is not a module score and was not used here.",
        "- PPP1R15A and ATF3 have no unambiguous precomputed Wu per-cell result column. Existing expression matrices were not read to create replacement values.",
        "",
        "## Program146 multi-caller sensitivity (numeric only)",
        "",
        "The recovered score was joined in memory by exact cell ID to the unchanged multi-caller master. Effects are within-patient Spearman correlations between the frozen Joint axis and recovered Program146 score. Cross-patient summaries are unweighted; patients, not cells, are the summary unit.",
        "",
    ]
    for _, row in summary_rows.iterrows():
        lines.append(
            f"- {row['submission_facing_subset_label']}: {int(row['n_subset_cells'])} cells; "
            f"{int(row['n_patients_positive'])}/{int(row['n_patients_testable'])} positive; "
            f"median rho={row['median_effect']:.3f} (range {row['min_effect']:.3f} to {row['max_effect']:.3f}); "
            f"direction reversals={int(row['n_direction_reversals'])}."
        )
    lines += [
        "",
        "## Search scope",
        "",
        f"- Candidate paths returned by the whole-project read-only keyword search: {examined}.",
        f"- Path/name keyword matches retained in the inventory: {len(inventory)}.",
        "- Installed package trees and external source caches were inventory-labelled and excluded from biological-result adjudication.",
        "- Detailed columns, dates, ID overlap, duplicates, patient mismatches, hashes, and provenance decisions are in `frozen_score_recovery_candidate_audit.tsv`.",
    ]
    (REPORTS / "frozen_score_recovery_audit.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
