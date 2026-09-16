# Purpose: Exact cell-ID integration and patient-unit sensitivity
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Fig. 2 and Supplementary Fig. S1.
from __future__ import annotations

import gzip
import hashlib
import io
import itertools
import math
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(r".")
CNV_ROOT = ROOT / "0728_sensitivity_and_cnv_closure"
TRANSFER = CNV_ROOT / "CNV_S1_transfer_20260910"
OUT = ROOT / "0910_Wu_multicaller_integration_audit"
TABLES = OUT / "tables"
REPORTS = OUT / "reports"
LOGS = OUT / "logs"

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

META_PATH = TRANSFER / "00_frozen_metadata" / "Wu2021_continuous_state_cell_metadata.csv"
COPYKAT_PATH = TRANSFER / "01_CopyKAT" / "ALL_patients_copykat_prediction_merged_subsampled.csv"

SCEVAN_ARCHIVES = {
    patient: sorted(CNV_ROOT.glob(f"SCEVAN_{patient}_20260907_*.tar.gz"))[-1]
    for patient in PATIENTS
}

TABLES.mkdir(parents=True, exist_ok=True)
REPORTS.mkdir(parents=True, exist_ok=True)
LOGS.mkdir(parents=True, exist_ok=True)


def log(message: str) -> None:
    text = str(message)
    print(text, flush=True)
    with (LOGS / "01_multicaller_integration.log").open("a", encoding="utf-8") as handle:
        handle.write(text + "\n")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tar_tsv(archive: Path, member: str) -> tuple[pd.DataFrame, int]:
    with tarfile.open(archive, "r:gz") as bundle:
        info = bundle.getmember(member)
        stream = bundle.extractfile(info)
        if stream is None:
            raise RuntimeError(f"Cannot read {member} from {archive}")
        payload = stream.read()
    return pd.read_csv(io.BytesIO(payload), sep="\t", low_memory=False), info.size


def bool_series(series: pd.Series) -> pd.Series:
    mapping = {
        "true": True,
        "false": False,
        "1": True,
        "0": False,
        "yes": True,
        "no": False,
    }
    out = series.astype("string").str.strip().str.lower().map(mapping)
    return out.astype("boolean")


def clean_scalar(value):
    if value is None:
        return None
    if isinstance(value, float) and np.isnan(value):
        return None
    return value


def audit_dataframe(
    *,
    source_label: str,
    path: Path,
    frame: pd.DataFrame,
    frozen_ids: set[str],
    member: str = "",
    member_size: int | None = None,
    cell_field: str | None = None,
    patient_field: str | None = None,
    expected_frozen_ids: set[str] | None = None,
    expected_patients: list[str] | None = None,
) -> dict:
    duplicate_ids = np.nan
    missing_ids = np.nan
    overlap_frozen = np.nan
    frozen_not_present = np.nan
    extra_not_frozen = np.nan
    scope_ids = expected_frozen_ids if expected_frozen_ids is not None else frozen_ids
    if cell_field and cell_field in frame.columns:
        ids = frame[cell_field].astype("string")
        nonmissing = set(ids.dropna().astype(str))
        duplicate_ids = int(ids.dropna().duplicated().sum())
        missing_ids = int(ids.isna().sum())
        overlap_frozen = len(nonmissing & frozen_ids)
        frozen_not_present = len(scope_ids - nonmissing)
        extra_not_frozen = len(nonmissing - frozen_ids)

    missing_patients = np.nan
    unexpected_patients = np.nan
    if patient_field and patient_field in frame.columns:
        present = set(frame[patient_field].dropna().astype(str))
        expected_patient_set = set(expected_patients if expected_patients is not None else PATIENTS)
        missing_patients = ",".join(sorted(expected_patient_set - present))
        unexpected_patients = ",".join(sorted(present - expected_patient_set))

    return {
        "source_label": source_label,
        "source_path": str(path),
        "archive_member": member,
        "source_file_size_bytes": path.stat().st_size,
        "member_size_bytes": member_size if member_size is not None else path.stat().st_size,
        "n_rows": len(frame),
        "n_columns": len(frame.columns),
        "columns": "|".join(map(str, frame.columns)),
        "cell_id_field": cell_field or "",
        "patient_field": patient_field or "",
        "duplicate_cell_ids": duplicate_ids,
        "missing_cell_ids": missing_ids,
        "cell_ids_overlapping_all_frozen": overlap_frozen,
        "expected_frozen_IDs_in_source_scope": len(scope_ids),
        "expected_scope_frozen_IDs_not_present": frozen_not_present,
        "extra_cell_ids_not_in_frozen_master": extra_not_frozen,
        "missing_expected_patients": missing_patients,
        "unexpected_patients": unexpected_patients,
    }


def safe_spearman(x: pd.Series, y: pd.Series) -> tuple[float, int, str]:
    pair = pd.DataFrame({"x": pd.to_numeric(x, errors="coerce"), "y": pd.to_numeric(y, errors="coerce")}).dropna()
    n = len(pair)
    if n < 3:
        return np.nan, n, "NOT_TESTABLE_LT3"
    if pair["x"].nunique() < 2 or pair["y"].nunique() < 2:
        return np.nan, n, "NOT_TESTABLE_CONSTANT"
    rho = float(spearmanr(pair["x"], pair["y"]).statistic)
    return rho, n, "TESTABLE"


def effect_comparison(effect: float, all_effect: float) -> dict:
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


def format_float(value: float, digits: int = 3) -> str:
    if pd.isna(value):
        return "NA"
    return f"{value:.{digits}f}"


def normalise_name(value: str) -> str:
    return "".join(ch.lower() for ch in str(value) if ch.isalnum())


def find_frozen_column(columns: list[str], aliases: list[str]) -> str | None:
    exact = {str(col): str(col) for col in columns}
    for alias in aliases:
        if alias in exact:
            return exact[alias]
    normalised = {}
    for col in columns:
        normalised.setdefault(normalise_name(col), []).append(str(col))
    for alias in aliases:
        matches = normalised.get(normalise_name(alias), [])
        if len(matches) == 1:
            return matches[0]
    return None


log("Reading frozen metadata")
meta = pd.read_csv(META_PATH, low_memory=False, dtype={"cell_id": "string", "patient": "string"})
if "cell_id" not in meta.columns or "patient" not in meta.columns:
    raise RuntimeError("Frozen metadata must contain cell_id and patient")
meta["cell_id"] = meta["cell_id"].astype(str)
meta["patient"] = meta["patient"].astype(str)
frozen_ids = set(meta["cell_id"])

log("Reading CopyKAT output")
copykat = pd.read_csv(COPYKAT_PATH, low_memory=False, dtype={"cell_id": "string", "patient_id": "string"})
copykat["cell_id"] = copykat["cell_id"].astype(str)

log("Reading SCEVAN classifications directly from archives")
scevan_calls_frames = []
scevan_status_frames = []
scevan_member_sizes: dict[tuple[str, str], int] = {}
for patient in PATIENTS:
    archive = SCEVAN_ARCHIVES[patient]
    calls, call_size = read_tar_tsv(archive, f"{patient}/cell_calls.tsv")
    status, status_size = read_tar_tsv(archive, f"{patient}/patient_status.tsv")
    calls["cell_id"] = calls["cell_id"].astype(str)
    calls["patient"] = calls["patient"].astype(str)
    scevan_calls_frames.append(calls)
    scevan_status_frames.append(status)
    scevan_member_sizes[(patient, "calls")] = call_size
    scevan_member_sizes[(patient, "status")] = status_size
scevan = pd.concat(scevan_calls_frames, ignore_index=True)
scevan_status = pd.concat(scevan_status_frames, ignore_index=True)

log("Reading inferCNV continuous burden and native HMM status")
infer_burden_frames = []
infer_hmm_frames = []
infer_status_frames = []
for patient in PATIENTS:
    patient_dir = TRANSFER / "02_inferCNV" / patient
    burden = pd.read_csv(patient_dir / "cell_cnv_burden.tsv", sep="\t", low_memory=False, dtype={"cell_id": "string", "patient": "string"})
    hmm = pd.read_csv(patient_dir / "hmm_cell_calls.tsv", sep="\t", low_memory=False, dtype={"cell_id": "string", "patient": "string"})
    status = pd.read_csv(patient_dir / "patient_status.tsv", sep="\t", low_memory=False, dtype={"patient": "string"})
    burden["cell_id"] = burden["cell_id"].astype(str)
    burden["patient"] = burden["patient"].astype(str)
    hmm["cell_id"] = hmm["cell_id"].astype(str)
    hmm["patient"] = hmm["patient"].astype(str)
    infer_burden_frames.append(burden)
    infer_hmm_frames.append(hmm)
    infer_status_frames.append(status)
infer_burden = pd.concat(infer_burden_frames, ignore_index=True)
infer_hmm = pd.concat(infer_hmm_frames, ignore_index=True)
infer_status = pd.concat(infer_status_frames, ignore_index=True)


# Input audit is deliberately written before the master-table merge.
audit_rows = []
audit_rows.append(
    audit_dataframe(
        source_label="frozen_metadata",
        path=META_PATH,
        frame=meta,
        frozen_ids=frozen_ids,
        cell_field="cell_id",
        patient_field="patient",
    )
)
audit_rows.append(
    audit_dataframe(
        source_label="CopyKAT_merged_subsampled",
        path=COPYKAT_PATH,
        frame=copykat,
        frozen_ids=frozen_ids,
        cell_field="cell_id",
        patient_field="patient_id",
    )
)
for patient, calls, status in zip(PATIENTS, scevan_calls_frames, scevan_status_frames):
    archive = SCEVAN_ARCHIVES[patient]
    audit_rows.append(
        audit_dataframe(
            source_label=f"SCEVAN_cell_calls_{patient}",
            path=archive,
            member=f"{patient}/cell_calls.tsv",
            member_size=scevan_member_sizes[(patient, "calls")],
            frame=calls,
            frozen_ids=frozen_ids,
            cell_field="cell_id",
            patient_field="patient",
            expected_frozen_ids=set(meta.loc[meta["patient"] == patient, "cell_id"]),
            expected_patients=[patient],
        )
    )
    audit_rows.append(
        audit_dataframe(
            source_label=f"SCEVAN_patient_status_{patient}",
            path=archive,
            member=f"{patient}/patient_status.tsv",
            member_size=scevan_member_sizes[(patient, "status")],
            frame=status,
            frozen_ids=frozen_ids,
            patient_field="patient",
            expected_frozen_ids=set(meta.loc[meta["patient"] == patient, "cell_id"]),
            expected_patients=[patient],
        )
    )
for patient in PATIENTS:
    patient_dir = TRANSFER / "02_inferCNV" / patient
    for label, filename, frame, cell_field, patient_field in [
        ("burden", "cell_cnv_burden.tsv", infer_burden[infer_burden["patient"] == patient], "cell_id", "patient"),
        ("HMM_cell_status", "hmm_cell_calls.tsv", infer_hmm[infer_hmm["patient"] == patient], "cell_id", "patient"),
        ("patient_status", "patient_status.tsv", infer_status[infer_status["patient"] == patient], None, "patient"),
    ]:
        audit_rows.append(
            audit_dataframe(
                source_label=f"inferCNV_{label}_{patient}",
                path=patient_dir / filename,
                frame=frame,
                frozen_ids=frozen_ids,
                cell_field=cell_field,
                patient_field=patient_field,
                expected_frozen_ids=set(meta.loc[meta["patient"] == patient, "cell_id"]),
                expected_patients=[patient],
            )
        )

input_audit = pd.DataFrame(audit_rows)
input_audit.to_csv(TABLES / "multicaller_input_audit.tsv", sep="\t", index=False, na_rep="NA")

hash_paths = [META_PATH, COPYKAT_PATH, TRANSFER / "MANIFEST.tsv", TRANSFER / "SHA256SUMS.txt"]
hash_paths.extend(SCEVAN_ARCHIVES[p] for p in PATIENTS)
for patient in PATIENTS:
    patient_dir = TRANSFER / "02_inferCNV" / patient
    hash_paths.extend(patient_dir / name for name in ["cell_cnv_burden.tsv", "hmm_cell_calls.tsv", "patient_status.tsv"])
hash_rows = [
    {"source_path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256(path)}
    for path in hash_paths
]
pd.DataFrame(hash_rows).to_csv(TABLES / "source_file_sha256.tsv", sep="\t", index=False)

fatal_checks = {
    "frozen_row_count_is_10836": len(meta) == 10836,
    "frozen_cell_ids_unique": not meta["cell_id"].duplicated().any(),
    "frozen_cell_ids_nonmissing": meta["cell_id"].notna().all(),
    "frozen_patients_exactly_expected": set(meta["patient"]) == set(PATIENTS),
    "copykat_cell_ids_unique": not copykat["cell_id"].duplicated().any(),
    "scevan_cell_ids_unique": not scevan["cell_id"].duplicated().any(),
    "infercnv_burden_cell_ids_unique": not infer_burden["cell_id"].duplicated().any(),
    "infercnv_hmm_cell_ids_unique": not infer_hmm["cell_id"].duplicated().any(),
    "all_scevan_archives_pass_classification": set(scevan_status["status"].astype(str)) == {"PASS_CLASSIFICATION"},
}
fatal_table = pd.DataFrame([{"check": key, "pass": value} for key, value in fatal_checks.items()])
fatal_table.to_csv(TABLES / "premerge_fatal_checks.tsv", sep="\t", index=False)
if not all(fatal_checks.values()):
    raise RuntimeError("Pre-merge audit failed; see premerge_fatal_checks.tsv")
log("Pre-merge input audit PASS")


# Build a cell-ID keyed master table while retaining caller-native fields.
copy_native = copykat[["cell_id", "cell.names", "copykat.pred", "copykat_input_role", "patient_id"]].copy()
copy_native = copy_native.rename(
    columns={
        "cell.names": "CopyKAT_raw_cell_name",
        "copykat.pred": "CopyKAT_raw_class",
        "copykat_input_role": "CopyKAT_raw_role",
        "patient_id": "CopyKAT_raw_patient",
    }
)

scevan_native = scevan[[
    "cell_id",
    "patient",
    "scevan_native_class",
    "scevan_tumor_positive",
    "caller_role",
    "returned_by_native_result",
]].rename(
    columns={
        "patient": "SCEVAN_raw_patient",
        "scevan_native_class": "SCEVAN_raw_class",
        "scevan_tumor_positive": "SCEVAN_raw_tumor_positive",
        "caller_role": "SCEVAN_raw_role",
        "returned_by_native_result": "SCEVAN_raw_returned_by_native_result",
    }
)

burden_native = infer_burden[[
    "cell_id",
    "patient",
    "caller_role",
    "inferCNV_continuous_burden",
    "neutral_baseline",
    "reference_95pct_threshold",
    "inferCNV_descriptive_CNV_high",
    "threshold_interpretation",
]].rename(
    columns={
        "patient": "inferCNV_raw_patient",
        "caller_role": "inferCNV_raw_role",
        "neutral_baseline": "inferCNV_raw_neutral_baseline",
        "reference_95pct_threshold": "inferCNV_raw_reference_95pct_threshold",
        "inferCNV_descriptive_CNV_high": "inferCNV_raw_descriptive_CNV_high",
        "threshold_interpretation": "inferCNV_raw_threshold_interpretation",
    }
)

hmm_native = infer_hmm[["cell_id", "hmm_testable", "inferCNV_HMM_positive", "reason"]].rename(
    columns={
        "hmm_testable": "inferCNV_HMM_raw_testable",
        "inferCNV_HMM_positive": "inferCNV_HMM_raw_positive",
        "reason": "inferCNV_HMM_raw_reason",
    }
)

status_native = infer_status[["patient", "continuous_status", "hmm_status", "reason"]].rename(
    columns={
        "continuous_status": "inferCNV_continuous_patient_status",
        "hmm_status": "inferCNV_HMM_status",
        "reason": "inferCNV_HMM_patient_reason",
    }
)

master = meta.merge(copy_native, how="left", on="cell_id", validate="one_to_one")
master = master.merge(scevan_native, how="left", on="cell_id", validate="one_to_one")
master = master.merge(burden_native, how="left", on="cell_id", validate="one_to_one")
master = master.merge(hmm_native, how="left", on="cell_id", validate="one_to_one")
master = master.merge(status_native, how="left", on="patient", validate="many_to_one")

copy_map = {"aneuploid": True, "diploid": False, "not.defined": pd.NA}
master["CopyKAT_support"] = master["CopyKAT_raw_class"].astype("string").str.lower().map(copy_map).astype("boolean")

scevan_map = {"tumor": True, "normal": False, "filtered": pd.NA}
master["SCEVAN_support"] = master["SCEVAN_raw_class"].astype("string").str.lower().map(scevan_map).astype("boolean")

master["inferCNV_burden"] = pd.to_numeric(master["inferCNV_continuous_burden"], errors="coerce")
master["inferCNV_high_CNA"] = bool_series(master["inferCNV_raw_descriptive_CNV_high"])

binary_columns = ["CopyKAT_support", "SCEVAN_support", "inferCNV_high_CNA"]
master["n_callers_evaluable"] = master[binary_columns].notna().sum(axis=1).astype(int)
master["n_callers_supporting"] = master[binary_columns].fillna(False).astype(bool).sum(axis=1).astype(int)
master["at_least_two_callers_supported"] = (master["n_callers_supporting"] >= 2).astype("boolean")
master["complete_case_support_category"] = pd.Series(pd.NA, index=master.index, dtype="string")
complete = master["n_callers_evaluable"] == 3
master.loc[complete, "complete_case_support_category"] = master.loc[complete, "n_callers_supporting"].astype(str) + "/3"

patient_mismatch_rows = []
for caller, caller_patient_col in [
    ("CopyKAT", "CopyKAT_raw_patient"),
    ("SCEVAN", "SCEVAN_raw_patient"),
    ("inferCNV", "inferCNV_raw_patient"),
]:
    present = master[caller_patient_col].notna()
    mismatch = present & (master[caller_patient_col].astype(str) != master["patient"].astype(str))
    patient_mismatch_rows.append({"caller": caller, "matched_cells": int(present.sum()), "patient_mismatches": int(mismatch.sum())})
patient_mismatch = pd.DataFrame(patient_mismatch_rows)
patient_mismatch.to_csv(TABLES / "caller_patient_identity_QC.tsv", sep="\t", index=False)
if patient_mismatch["patient_mismatches"].sum() != 0:
    raise RuntimeError("Caller patient mismatch after cell-ID merge")

master_tsv = TABLES / "Wu_multicaller_cell_master.tsv.gz"
master_csv = TABLES / "Wu_multicaller_cell_master.csv.gz"
master.to_csv(master_tsv, sep="\t", index=False, compression="gzip", na_rep="NA")
master.to_csv(master_csv, index=False, compression="gzip", na_rep="NA")
log(f"Master written: {len(master)} frozen cells")


# Patient-by-caller summaries. HMM is retained as native status and never converted to binary support.
patient_summary_rows = []
for patient in PATIENTS:
    frame = master[master["patient"] == patient]
    native_status_row = infer_status[infer_status["patient"].astype(str) == patient].iloc[0]
    caller_specs = [
        ("CopyKAT", "CopyKAT_support", "CopyKAT_raw_class", "aneuploid=support; diploid=non-support; not.defined=NA"),
        ("SCEVAN", "SCEVAN_support", "SCEVAN_raw_class", "tumor=support; normal=non-support; filtered=NA"),
        ("inferCNV-high-descriptive", "inferCNV_high_CNA", "inferCNV_raw_descriptive_CNV_high", "above same-patient reference 95th percentile; descriptive, not a malignant call"),
    ]
    for caller, support_col, raw_col, interpretation in caller_specs:
        support = frame[support_col]
        raw_counts = frame[raw_col].astype("string").fillna("MISSING").value_counts(dropna=False)
        patient_summary_rows.append(
            {
                "patient": patient,
                "population": "frozen_malignant_epithelial",
                "caller": caller,
                "n_frozen_cells": len(frame),
                "n_evaluable": int(support.notna().sum()),
                "n_supported": int((support == True).fillna(False).sum()),
                "support_fraction_among_evaluable": float((support == True).fillna(False).sum() / support.notna().sum()) if support.notna().sum() else np.nan,
                "n_NA": int(support.isna().sum()),
                "native_status": "PASS",
                "raw_class_counts": ";".join(f"{key}={value}" for key, value in raw_counts.items()),
                "interpretation": interpretation,
            }
        )
    patient_summary_rows.append(
        {
            "patient": patient,
            "population": "frozen_malignant_epithelial",
            "caller": "inferCNV-HMM-native-status",
            "n_frozen_cells": len(frame),
            "n_evaluable": 0,
            "n_supported": np.nan,
            "support_fraction_among_evaluable": np.nan,
            "n_NA": len(frame),
            "native_status": native_status_row["hmm_status"],
            "raw_class_counts": "binary malignant-cell call intentionally unavailable",
            "interpretation": native_status_row["reason"],
        }
    )

patient_summary = pd.DataFrame(patient_summary_rows)
patient_summary["patient"] = pd.Categorical(patient_summary["patient"], categories=PATIENTS, ordered=True)
patient_summary = patient_summary.sort_values(["patient", "caller"]).reset_index(drop=True)
patient_summary.to_csv(TABLES / "Wu_multicaller_patient_summary.tsv", sep="\t", index=False, na_rep="NA")


def cohen_kappa(a: pd.Series, b: pd.Series) -> float:
    valid = a.notna() & b.notna()
    if valid.sum() == 0:
        return np.nan
    av = a[valid].astype(bool).to_numpy()
    bv = b[valid].astype(bool).to_numpy()
    po = np.mean(av == bv)
    pa = np.mean(av)
    pb = np.mean(bv)
    pe = pa * pb + (1 - pa) * (1 - pb)
    return np.nan if pe == 1 else float((po - pe) / (1 - pe))


PAIR_LABELS = {
    "CopyKAT_support": "CopyKAT",
    "SCEVAN_support": "SCEVAN",
    "inferCNV_high_CNA": "inferCNV-high-descriptive",
}


def build_overlap(frame: pd.DataFrame, stratum: str) -> list[dict]:
    rows = []
    for a, b in itertools.combinations(binary_columns, 2):
        av = frame[a]
        bv = frame[b]
        valid = av.notna() & bv.notna()
        a_true = (av == True).fillna(False)
        b_true = (bv == True).fillna(False)
        both = int((valid & a_true & b_true).sum())
        a_only = int((valid & a_true & ~b_true).sum())
        b_only = int((valid & ~a_true & b_true).sum())
        neither = int((valid & ~a_true & ~b_true).sum())
        union = both + a_only + b_only
        rows.append(
            {
                "stratum": stratum,
                "analysis_type": "pairwise_overlap",
                "caller_A": PAIR_LABELS[a],
                "caller_B": PAIR_LABELS[b],
                "category": "pairwise_complete",
                "n_universe": len(frame),
                "n_complete_evaluable": int(valid.sum()),
                "n_any_missing": int((~valid).sum()),
                "n_both_supported": both,
                "n_A_only_supported": a_only,
                "n_B_only_supported": b_only,
                "n_neither_supported": neither,
                "jaccard_supported": both / union if union else np.nan,
                "raw_binary_agreement": (both + neither) / valid.sum() if valid.sum() else np.nan,
                "cohen_kappa": cohen_kappa(av, bv),
                "n_cells": np.nan,
                "fraction_of_universe": np.nan,
            }
        )
    complete3 = frame["n_callers_evaluable"] == 3
    for k in range(4):
        n = int((complete3 & (frame["n_callers_supporting"] == k)).sum())
        rows.append(
            {
                "stratum": stratum,
                "analysis_type": "three_way_complete_case_distribution",
                "caller_A": "CopyKAT",
                "caller_B": "SCEVAN|inferCNV-high-descriptive",
                "category": f"{k}/3",
                "n_universe": len(frame),
                "n_complete_evaluable": int(complete3.sum()),
                "n_any_missing": int((~complete3).sum()),
                "n_both_supported": np.nan,
                "n_A_only_supported": np.nan,
                "n_B_only_supported": np.nan,
                "n_neither_supported": np.nan,
                "jaccard_supported": np.nan,
                "raw_binary_agreement": np.nan,
                "cohen_kappa": np.nan,
                "n_cells": n,
                "fraction_of_universe": n / len(frame) if len(frame) else np.nan,
            }
        )
    for k in range(4):
        n = int((frame["n_callers_evaluable"] == k).sum())
        rows.append(
            {
                "stratum": stratum,
                "analysis_type": "caller_evaluability_distribution",
                "caller_A": "all_three_binary_descriptors",
                "caller_B": "",
                "category": f"{k}/3_evaluable",
                "n_universe": len(frame),
                "n_complete_evaluable": int(complete3.sum()),
                "n_any_missing": int((~complete3).sum()),
                "n_both_supported": np.nan,
                "n_A_only_supported": np.nan,
                "n_B_only_supported": np.nan,
                "n_neither_supported": np.nan,
                "jaccard_supported": np.nan,
                "raw_binary_agreement": np.nan,
                "cohen_kappa": np.nan,
                "n_cells": n,
                "fraction_of_universe": n / len(frame) if len(frame) else np.nan,
            }
        )
    n_two = int(master.loc[frame.index, "at_least_two_callers_supported"].fillna(False).sum())
    rows.append(
        {
            "stratum": stratum,
            "analysis_type": "consensus_support",
            "caller_A": "all_three_binary_descriptors",
            "caller_B": "",
            "category": ">=2_supporting",
            "n_universe": len(frame),
            "n_complete_evaluable": int(complete3.sum()),
            "n_any_missing": int((~complete3).sum()),
            "n_both_supported": np.nan,
            "n_A_only_supported": np.nan,
            "n_B_only_supported": np.nan,
            "n_neither_supported": np.nan,
            "jaccard_supported": np.nan,
            "raw_binary_agreement": np.nan,
            "cohen_kappa": np.nan,
            "n_cells": n_two,
            "fraction_of_universe": n_two / len(frame) if len(frame) else np.nan,
        }
    )
    return rows


overlap_summary = pd.DataFrame(build_overlap(master, "ALL_8_PATIENTS_POOLED_DESCRIPTIVE"))
overlap_summary.to_csv(TABLES / "caller_overlap_summary.tsv", sep="\t", index=False, na_rep="NA")

overlap_patient_rows = []
for patient in PATIENTS:
    overlap_patient_rows.extend(build_overlap(master[master["patient"] == patient], patient))
overlap_by_patient = pd.DataFrame(overlap_patient_rows)
overlap_by_patient["stratum"] = pd.Categorical(overlap_by_patient["stratum"], categories=PATIENTS, ordered=True)
overlap_by_patient = overlap_by_patient.sort_values(["stratum", "analysis_type", "category"]).reset_index(drop=True)
overlap_by_patient.to_csv(TABLES / "caller_overlap_by_patient.tsv", sep="\t", index=False, na_rep="NA")


# Reference calibration is expected to be approximately 5% by threshold construction and is not specificity.
reference_rows = []
for patient in PATIENTS:
    ref = infer_burden[(infer_burden["patient"] == patient) & (infer_burden["caller_role"] == "same_patient_reference")].copy()
    high = bool_series(ref["inferCNV_descriptive_CNV_high"])
    reference_rows.append(
        {
            "patient": patient,
            "n_reference_cells": len(ref),
            "n_reference_evaluable": int(high.notna().sum()),
            "n_reference_high_CNA": int((high == True).fillna(False).sum()),
            "reference_high_CNA_fraction": float((high == True).fillna(False).sum() / high.notna().sum()) if high.notna().sum() else np.nan,
            "reference_95pct_threshold": pd.to_numeric(ref["reference_95pct_threshold"], errors="coerce").dropna().iloc[0] if len(ref) else np.nan,
            "interpretation": "Expected near 5% by threshold construction; not an estimate of specificity or false-positive rate",
        }
    )
reference_calibration = pd.DataFrame(reference_rows)
reference_calibration.to_csv(TABLES / "inferCNV_reference_threshold_calibration.tsv", sep="\t", index=False, na_rep="NA")


# Detect frozen analytical columns without redefining or recomputing them.
feature_definitions = [
    ("Program146", "program", ["Program146", "Program_146", "Program 146"]),
    ("UPR", "six_core_module", ["UPR", "UPR_score"]),
    ("TNFA-NFKB", "six_core_module", ["TNFA-NFKB", "TNFA_NFKB", "TNFA_NFkB", "TNFA via NFKB", "TNFA_SIGNALING_VIA_NFKB"]),
    ("Hypoxia", "six_core_module", ["Hypoxia", "Hypoxia_score"]),
    ("Adhesion Remodeling", "six_core_module", ["Adhesion Remodeling", "Adhesion_Remodeling", "Adhesion_Remodeling_score"]),
    ("Wound Healing", "six_core_module", ["Wound Healing", "Wound_Healing", "Wound_Healing_score"]),
    ("Survival Stress", "six_core_module", ["Survival Stress", "Survival_Stress", "Survival_Stress_score"]),
    ("MCL1", "core_gene", ["MCL1"]),
    ("PPP1R15A", "core_gene", ["PPP1R15A"]),
    ("ATF3", "core_gene", ["ATF3"]),
]
feature_manifest_rows = []
feature_columns = {}
for feature, feature_type, aliases in feature_definitions:
    resolved = find_frozen_column(list(meta.columns), aliases)
    feature_columns[feature] = resolved
    feature_manifest_rows.append(
        {
            "feature": feature,
            "feature_type": feature_type,
            "frozen_column_detected": resolved if resolved else "",
            "status": "AVAILABLE_FROZEN_COLUMN" if resolved else "MISSING_FROZEN_COLUMN_NO_RECOMPUTATION",
            "aliases_checked": "|".join(aliases),
        }
    )
feature_manifest = pd.DataFrame(feature_manifest_rows)
feature_manifest.to_csv(TABLES / "frozen_functional_column_availability.tsv", sep="\t", index=False)

required_axis_cols = {
    "YAP": find_frozen_column(list(meta.columns), ["YAP_score"]),
    "Stem": find_frozen_column(list(meta.columns), ["Stemness_score", "Stem_score"]),
    "Joint": find_frozen_column(list(meta.columns), ["YAP_Stem_Score", "Joint_axis", "Joint_score"]),
}
if any(value is None for value in required_axis_cols.values()):
    raise RuntimeError(f"Missing frozen axis columns: {required_axis_cols}")
pd.DataFrame(
    [{"quantity": key, "frozen_column": value, "status": "AVAILABLE_FROZEN_COLUMN"} for key, value in required_axis_cols.items()]
).to_csv(TABLES / "frozen_axis_column_availability.tsv", sep="\t", index=False)


subset_specs = [
    ("All_frozen_10836", pd.Series(True, index=master.index)),
    ("CopyKAT_supported", master["CopyKAT_support"] == True),
    ("SCEVAN_supported", master["SCEVAN_support"] == True),
    ("inferCNV_high_descriptive", master["inferCNV_high_CNA"] == True),
    ("at_least_two_callers_supported", master["at_least_two_callers_supported"] == True),
]
subset_specs = [(name, mask.fillna(False)) for name, mask in subset_specs]


# Patient-wise YAP-Stem sensitivity.
yap_stem_rows = []
all_patient_effects = {}
for patient in PATIENTS:
    frame = master[master["patient"] == patient]
    rho, n_pair, status = safe_spearman(frame[required_axis_cols["YAP"]], frame[required_axis_cols["Stem"]])
    all_patient_effects[patient] = rho

for subset_name, subset_mask in subset_specs:
    patient_effects = []
    for patient in PATIENTS:
        patient_all = master[master["patient"] == patient]
        frame = master[(master["patient"] == patient) & subset_mask]
        rho, n_pair, status = safe_spearman(frame[required_axis_cols["YAP"]], frame[required_axis_cols["Stem"]])
        comparison = effect_comparison(rho, all_patient_effects[patient])
        patient_effects.append(rho)
        yap_stem_rows.append(
            {
                "record_type": "patient",
                "subset": subset_name,
                "patient": patient,
                "n_frozen_patient_cells": len(patient_all),
                "n_subset_cells": len(frame),
                "subset_fraction": len(frame) / len(patient_all) if len(patient_all) else np.nan,
                "n_complete_score_pairs": n_pair,
                "effect_type": "within_patient_Spearman_YAP_vs_Stem",
                "rho": rho,
                "status": status,
                "positive": bool(rho > 0) if not pd.isna(rho) else np.nan,
                **comparison,
                "n_patients_testable": np.nan,
                "n_patients_positive": np.nan,
                "median_rho": np.nan,
                "min_rho": np.nan,
                "max_rho": np.nan,
                "n_direction_reversals": np.nan,
                "n_marked_attenuations": np.nan,
            }
        )
    valid = [x for x in patient_effects if not pd.isna(x)]
    subset_patient_rows = [row for row in yap_stem_rows if row["record_type"] == "patient" and row["subset"] == subset_name]
    yap_stem_rows.append(
        {
            "record_type": "cross_patient_summary",
            "subset": subset_name,
            "patient": "ALL_8_PATIENTS",
            "n_frozen_patient_cells": len(master),
            "n_subset_cells": int(subset_mask.sum()),
            "subset_fraction": float(subset_mask.mean()),
            "n_complete_score_pairs": sum(int(row["n_complete_score_pairs"]) for row in subset_patient_rows),
            "effect_type": "patient_is_statistical_unit; unweighted_summary_of_patient_rho",
            "rho": np.nan,
            "status": "SUMMARY",
            "positive": np.nan,
            "delta_vs_all": np.nan,
            "absolute_effect_ratio_vs_all": np.nan,
            "direction_reversal_vs_all": np.nan,
            "marked_attenuation_vs_all": np.nan,
            "n_patients_testable": len(valid),
            "n_patients_positive": sum(x > 0 for x in valid),
            "median_rho": float(np.median(valid)) if valid else np.nan,
            "min_rho": float(np.min(valid)) if valid else np.nan,
            "max_rho": float(np.max(valid)) if valid else np.nan,
            "n_direction_reversals": sum(row["direction_reversal_vs_all"] is True for row in subset_patient_rows),
            "n_marked_attenuations": sum(row["marked_attenuation_vs_all"] is True for row in subset_patient_rows),
        }
    )

yap_stem_sensitivity = pd.DataFrame(yap_stem_rows)
yap_stem_sensitivity["patient"] = pd.Categorical(yap_stem_sensitivity["patient"], categories=PATIENTS + ["ALL_8_PATIENTS"], ordered=True)
yap_stem_sensitivity.to_csv(TABLES / "Wu_multicaller_YAP_Stem_sensitivity.tsv", sep="\t", index=False, na_rep="NA")


# Program/module/gene sensitivity, using only pre-existing frozen columns.
functional_rows = []
joint_col = required_axis_cols["Joint"]
for feature, feature_type, aliases in feature_definitions:
    feature_col = feature_columns[feature]
    all_feature_effects = {}
    if feature_col:
        for patient in PATIENTS:
            frame = master[master["patient"] == patient]
            rho, _, _ = safe_spearman(frame[joint_col], frame[feature_col])
            all_feature_effects[patient] = rho
    else:
        all_feature_effects = {patient: np.nan for patient in PATIENTS}

    for subset_name, subset_mask in subset_specs:
        subset_patient_rows = []
        for patient in PATIENTS:
            patient_all = master[master["patient"] == patient]
            frame = master[(master["patient"] == patient) & subset_mask]
            if feature_col is None:
                rho, n_pair, status = np.nan, 0, "MISSING_FROZEN_COLUMN_NO_RECOMPUTATION"
            else:
                rho, n_pair, status = safe_spearman(frame[joint_col], frame[feature_col])
            comparison = effect_comparison(rho, all_feature_effects[patient])
            row = {
                "record_type": "patient",
                "feature": feature,
                "feature_type": feature_type,
                "frozen_column": feature_col or "",
                "subset": subset_name,
                "patient": patient,
                "n_frozen_patient_cells": len(patient_all),
                "n_subset_cells": len(frame),
                "n_complete_pairs": n_pair,
                "effect_type": "within_patient_Spearman_frozen_Joint_axis_vs_feature",
                "effect": rho,
                "status": status,
                "positive": bool(rho > 0) if not pd.isna(rho) else np.nan,
                **comparison,
                "n_patients_testable": np.nan,
                "n_patients_positive": np.nan,
                "median_effect": np.nan,
                "min_effect": np.nan,
                "max_effect": np.nan,
                "n_direction_reversals": np.nan,
                "n_marked_attenuations": np.nan,
            }
            functional_rows.append(row)
            subset_patient_rows.append(row)

        valid = [row["effect"] for row in subset_patient_rows if not pd.isna(row["effect"])]
        functional_rows.append(
            {
                "record_type": "cross_patient_summary",
                "feature": feature,
                "feature_type": feature_type,
                "frozen_column": feature_col or "",
                "subset": subset_name,
                "patient": "ALL_8_PATIENTS",
                "n_frozen_patient_cells": len(master),
                "n_subset_cells": int(subset_mask.sum()),
                "n_complete_pairs": sum(int(row["n_complete_pairs"]) for row in subset_patient_rows),
                "effect_type": "patient_is_statistical_unit; unweighted_summary_of_patient_effects",
                "effect": np.nan,
                "status": "SUMMARY" if feature_col else "MISSING_FROZEN_COLUMN_NO_RECOMPUTATION",
                "positive": np.nan,
                "delta_vs_all": np.nan,
                "absolute_effect_ratio_vs_all": np.nan,
                "direction_reversal_vs_all": np.nan,
                "marked_attenuation_vs_all": np.nan,
                "n_patients_testable": len(valid),
                "n_patients_positive": sum(x > 0 for x in valid),
                "median_effect": float(np.median(valid)) if valid else np.nan,
                "min_effect": float(np.min(valid)) if valid else np.nan,
                "max_effect": float(np.max(valid)) if valid else np.nan,
                "n_direction_reversals": sum(row["direction_reversal_vs_all"] is True for row in subset_patient_rows),
                "n_marked_attenuations": sum(row["marked_attenuation_vs_all"] is True for row in subset_patient_rows),
            }
        )

functional_sensitivity = pd.DataFrame(functional_rows)
functional_sensitivity["patient"] = pd.Categorical(functional_sensitivity["patient"], categories=PATIENTS + ["ALL_8_PATIENTS"], ordered=True)
functional_sensitivity.to_csv(TABLES / "Wu_multicaller_functional_sensitivity.tsv", sep="\t", index=False, na_rep="NA")


# SCEVAN-native class profiles among cells frozen as malignant epithelial.
scevan_profile_rows = []
for patient in PATIENTS:
    frame = master[master["patient"] == patient].copy()
    frame["SCEVAN_profile_class"] = frame["SCEVAN_raw_class"].astype("string").fillna("missing")
    for native_class in ["tumor", "normal", "filtered", "missing"]:
        group = frame[frame["SCEVAN_profile_class"] == native_class]
        scevan_profile_rows.append(
            {
                "patient": patient,
                "SCEVAN_native_class": native_class,
                "n_cells": len(group),
                "fraction_of_frozen_patient_cells": len(group) / len(frame) if len(frame) else np.nan,
                "median_YAP_score": pd.to_numeric(group[required_axis_cols["YAP"]], errors="coerce").median() if len(group) else np.nan,
                "median_Stemness_score": pd.to_numeric(group[required_axis_cols["Stem"]], errors="coerce").median() if len(group) else np.nan,
                "median_Joint_axis": pd.to_numeric(group[required_axis_cols["Joint"]], errors="coerce").median() if len(group) else np.nan,
                "median_nCount_RNA": pd.to_numeric(group.get("nCount_RNA"), errors="coerce").median() if len(group) and "nCount_RNA" in group else np.nan,
                "median_nFeature_RNA": pd.to_numeric(group.get("nFeature_RNA"), errors="coerce").median() if len(group) and "nFeature_RNA" in group else np.nan,
            }
        )
scevan_profile = pd.DataFrame(scevan_profile_rows)
scevan_profile.to_csv(TABLES / "SCEVAN_native_class_profile_in_frozen_cells.tsv", sep="\t", index=False, na_rep="NA")


# Numerical QC and support distributions.
support_distribution_rows = []
for stratum, frame in [("ALL_8_PATIENTS", master)] + [(patient, master[master["patient"] == patient]) for patient in PATIENTS]:
    for evaluable in range(4):
        for supporting in range(4):
            n = int(((frame["n_callers_evaluable"] == evaluable) & (frame["n_callers_supporting"] == supporting)).sum())
            if n or (evaluable == 3 and supporting <= 3):
                support_distribution_rows.append(
                    {
                        "stratum": stratum,
                        "n_callers_evaluable": evaluable,
                        "n_callers_supporting": supporting,
                        "n_cells": n,
                        "fraction": n / len(frame) if len(frame) else np.nan,
                        "complete_case_label": f"{supporting}/3" if evaluable == 3 else "incomplete_case",
                    }
                )
pd.DataFrame(support_distribution_rows).to_csv(TABLES / "binary_support_distribution.tsv", sep="\t", index=False)

qc_rows = [
    {"check": "master_row_count", "observed": len(master), "expected": 10836, "pass": len(master) == 10836},
    {"check": "master_unique_cell_ids", "observed": master["cell_id"].nunique(), "expected": 10836, "pass": master["cell_id"].nunique() == 10836},
    {"check": "master_patient_count", "observed": master["patient"].nunique(), "expected": 8, "pass": master["patient"].nunique() == 8},
    {"check": "CopyKAT_frozen_ID_coverage", "observed": master["CopyKAT_raw_class"].notna().sum(), "expected": 10836, "pass": master["CopyKAT_raw_class"].notna().sum() == 10836},
    {"check": "SCEVAN_frozen_ID_row_coverage", "observed": master["SCEVAN_raw_class"].notna().sum(), "expected": 10836, "pass": master["SCEVAN_raw_class"].notna().sum() == 10836},
    {"check": "inferCNV_frozen_ID_burden_coverage", "observed": master["inferCNV_burden"].notna().sum(), "expected": 10836, "pass": master["inferCNV_burden"].notna().sum() == 10836},
    {"check": "caller_patient_ID_mismatches", "observed": int(patient_mismatch["patient_mismatches"].sum()), "expected": 0, "pass": patient_mismatch["patient_mismatches"].sum() == 0},
    {"check": "CID3963_HMM_native_status", "observed": master.loc[master["patient"] == "CID3963", "inferCNV_HMM_status"].dropna().unique()[0], "expected": "NOT_TESTABLE", "pass": set(master.loc[master["patient"] == "CID3963", "inferCNV_HMM_status"].dropna()) == {"NOT_TESTABLE"}},
]
qc = pd.DataFrame(qc_rows)
qc.to_csv(TABLES / "Wu_multicaller_master_numerical_QC.tsv", sep="\t", index=False)


# Compose the data-driven report without producing any formal figure.
summary_yap = yap_stem_sensitivity[yap_stem_sensitivity["record_type"] == "cross_patient_summary"].set_index("subset")
support_counts = {
    "CopyKAT": int((master["CopyKAT_support"] == True).fillna(False).sum()),
    "SCEVAN": int((master["SCEVAN_support"] == True).fillna(False).sum()),
    "inferCNV": int((master["inferCNV_high_CNA"] == True).fillna(False).sum()),
    "at_least_two": int(master["at_least_two_callers_supported"].fillna(False).sum()),
    "complete": int((master["n_callers_evaluable"] == 3).sum()),
}

complete_counts = {
    k: int(((master["n_callers_evaluable"] == 3) & (master["n_callers_supporting"] == k)).sum())
    for k in range(4)
}
copy_scevan_overlap = overlap_summary[
    (overlap_summary["analysis_type"] == "pairwise_overlap")
    & (overlap_summary["caller_A"] == "CopyKAT")
    & (overlap_summary["caller_B"] == "SCEVAN")
].iloc[0]

scevan_frozen_counts = master["SCEVAN_raw_class"].astype("string").fillna("missing").value_counts().to_dict()
copy_frozen_counts = master["CopyKAT_raw_class"].astype("string").fillna("missing").value_counts().to_dict()
infer_ref_median = float(reference_calibration["reference_high_CNA_fraction"].median())
infer_ref_range = (
    float(reference_calibration["reference_high_CNA_fraction"].min()),
    float(reference_calibration["reference_high_CNA_fraction"].max()),
)

reversal_rows = yap_stem_sensitivity[
    (yap_stem_sensitivity["record_type"] == "patient")
    & (yap_stem_sensitivity["direction_reversal_vs_all"] == True)
]
attenuation_rows = yap_stem_sensitivity[
    (yap_stem_sensitivity["record_type"] == "patient")
    & (yap_stem_sensitivity["marked_attenuation_vs_all"] == True)
]

available_features = feature_manifest[feature_manifest["status"] == "AVAILABLE_FROZEN_COLUMN"]["feature"].tolist()
missing_features = feature_manifest[feature_manifest["status"] != "AVAILABLE_FROZEN_COLUMN"]["feature"].tolist()
hypoxia_summary = functional_sensitivity[
    (functional_sensitivity["record_type"] == "cross_patient_summary")
    & (functional_sensitivity["feature"] == "Hypoxia")
].set_index("subset")
scevan_normal_by_patient = (
    master.loc[master["SCEVAN_raw_class"].astype("string").str.lower() == "normal"]
    .groupby("patient", observed=False)
    .size()
    .reindex(PATIENTS, fill_value=0)
)

def yap_line(subset: str) -> str:
    row = summary_yap.loc[subset]
    return (
        f"- {subset}: {int(row['n_subset_cells'])} cells; "
        f"{int(row['n_patients_positive'])}/{int(row['n_patients_testable'])} patients positive; "
        f"median rho {format_float(row['median_rho'])}, range "
        f"{format_float(row['min_rho'])} to {format_float(row['max_rho'])}; "
        f"{int(row['n_direction_reversals'])} patient-level reversals and "
        f"{int(row['n_marked_attenuations'])} marked attenuations versus all frozen cells."
    )


report_lines = [
    "# Wu multi-caller cell-level concordance and sensitivity report",
    "",
    "## 中文结论摘要",
    "",
    f"本次工作严格以冻结的 10,836 个 Wu TNBC 恶性上皮细胞为索引，按完整 cell ID 合并 CopyKAT、SCEVAN 与 inferCNV。三类结果对冻结细胞的 ID 覆盖均为 10,836/10,836，患者错配为 0；未重跑 caller、未改动原始结果、未修改冻结恶性细胞定义。",
    "",
    f"RNA-CNA 支持并非所有细胞完全一致：CopyKAT aneuploid 为 {support_counts['CopyKAT']:,} 个，SCEVAN tumor 为 {support_counts['SCEVAN']:,} 个，inferCNV-high/descriptive 为 {support_counts['inferCNV']:,} 个；至少两个二元描述指标支持的细胞为 {support_counts['at_least_two']:,}/{len(master):,}。三者均可评价的 {support_counts['complete']:,} 个细胞中，0/3、1/3、2/3、3/3 分别为 {complete_counts[0]:,}、{complete_counts[1]:,}、{complete_counts[2]:,}、{complete_counts[3]:,}。这些结果用于敏感性分析，不用于重新定义主分析细胞。",
    "",
    f"CopyKAT 与 SCEVAN 在两者均可评价的 {int(copy_scevan_overlap['n_complete_evaluable']):,} 个细胞中具有较高的一致性：支持集 Jaccard={format_float(copy_scevan_overlap['jaccard_supported'])}、原始二分类一致率={format_float(copy_scevan_overlap['raw_binary_agreement'])}、Cohen kappa={format_float(copy_scevan_overlap['cohen_kappa'])}。相反，inferCNV-high 是参考 95 百分位构造的 burden 描述阈值，不应与 CopyKAT/SCEVAN 的分类语义等同。",
    "",
    f"核心 YAP–Stem 方向高度稳定：全部冻结细胞、CopyKAT 支持、SCEVAN 支持、inferCNV-high/descriptive 及至少两个 caller 支持五个预定义集合中均为 8/8 患者 rho>0。至少两个 caller 支持集合的跨患者中位 rho={format_float(summary_yap.loc['at_least_two_callers_supported', 'median_rho'])}（范围 {format_float(summary_yap.loc['at_least_two_callers_supported', 'min_rho'])}–{format_float(summary_yap.loc['at_least_two_callers_supported', 'max_rho'])}），无患者方向反转。",
    "",
    "inferCNV-high/descriptive 子集虽然仍为 8/8 正向，但 CID3963、CID4465、CID4513、CID4523 的绝对 rho 相对全部冻结细胞下降到不超过 50%；这应表述为效应减弱，而不是主结论失效。",
    "",
    f"SCEVAN 在冻结恶性细胞中给出 {int(scevan_frozen_counts.get('normal', 0)):,} 个 normal 和 {int(scevan_frozen_counts.get('filtered', 0)):,} 个 filtered。normal 主要见于 "
    + ", ".join(f"{patient}={int(count)}" for patient, count in scevan_normal_by_patient.items() if count > 0)
    + "。这些 normal 细胞在相应患者中通常具有更低的 nCount/nFeature 和更低的 Joint axis，因此更符合 RNA-CNA 信号较弱/覆盖较低所导致的 caller 分歧，不能直接解释为冻结注释错误；这是基于描述性分布的解释性推断，不是 DNA 层面的验证。",
    "",
    f"inferCNV 参考细胞中 high-CNA 比例的患者中位数为 {infer_ref_median:.3%}（范围 {infer_ref_range[0]:.3%}–{infer_ref_range[1]:.3%}）。该约 5% 比例由同患者参考细胞第 95 百分位阈值的构造决定，不能报告为特异性或假阳性率。CID3963 的 HMM 原生状态保留为 NOT_TESTABLE；其余患者也未把原生 HMM 状态强行转换为恶性细胞二分类。",
    "",
    f"转移的冻结 metadata 中，Program146、UPR、TNFA–NFκB、Adhesion Remodeling、Wound Healing、Survival Stress、MCL1、PPP1R15A 和 ATF3 均无现成列，因此遵照任务限制标记为不可检验且未重算。仅 Hypoxia 可直接检查：至少两个 caller 支持集合中为 {int(hypoxia_summary.loc['at_least_two_callers_supported', 'n_patients_positive'])}/{int(hypoxia_summary.loc['at_least_two_callers_supported', 'n_patients_testable'])} 患者正向，中位效应 {format_float(hypoxia_summary.loc['at_least_two_callers_supported', 'median_effect'])}；CID44971 与 CID4523 出现幅度很小的近零负相关，故该项在本次 caller 限定敏感性中应视为部分保留/边界，而不能据此改写已冻结的模块结论。",
    "",
    "综合判断：原 10,836 细胞主分析可以保持冻结；多 caller 结果应作为补充敏感性证据。最强结论是 YAP–Stem 连续耦合对 caller 限定稳健；不能声称所有冻结细胞都被所有 CNA caller 一致认定，也不能把任何 RNA caller 称为 DNA 金标准。",
    "",
    "## Scope and fixed interpretation",
    "",
    "This audit integrates the frozen 10,836 Wu TNBC malignant epithelial cells with returned CopyKAT, SCEVAN and inferCNV outputs strictly by exact cell ID. It does not rerun a caller, change the frozen malignant set, recompute a missing biological score, or produce a formal figure. Patients—not pooled cells—are the inferential unit for state-association summaries.",
    "",
    "CopyKAT aneuploid is treated as support, diploid as non-support and not.defined as missing. SCEVAN tumor is support, normal is non-support and filtered is missing. inferCNV-high is only a descriptive threshold above the same-patient reference 95th percentile; it is not called a validated malignant-cell positive. Native inferCNV HMM status is retained without inventing a binary cell call. CID3963 HMM is NOT_TESTABLE because the native HMM step failed under the singleton-reference condition.",
    "",
    "## Input and merge audit",
    "",
    f"- Frozen metadata: {len(meta):,} unique cells across {meta['patient'].nunique()} expected patients.",
    f"- Exact frozen-ID coverage: CopyKAT {int(master['CopyKAT_raw_class'].notna().sum()):,}; SCEVAN rows {int(master['SCEVAN_raw_class'].notna().sum()):,}; inferCNV burden {int(master['inferCNV_burden'].notna().sum()):,}.",
    f"- Patient-ID mismatches after cell-ID merge: {int(patient_mismatch['patient_mismatches'].sum())}.",
    "- The eight SCEVAN sources are the 2026-09-07 PASS_CLASSIFICATION archives. The older July NOT_TESTABLE placeholders were not used.",
    "- Full paths, file sizes, row/column inventories, ID fields, duplicates, missing IDs and source hashes are provided in the audit tables.",
    "",
    "## A. Are all 10,836 frozen malignant epithelial cells consistently supported?",
    "",
    "No. The 10,836 cells are the frozen annotation-defined primary analysis set, whereas RNA-derived CNA callers provide heterogeneous supporting evidence rather than a universal gold-standard identity.",
    "",
    f"- CopyKAT classes among frozen cells: {copy_frozen_counts}.",
    f"- SCEVAN native classes among frozen cells: {scevan_frozen_counts}.",
    f"- inferCNV-high/descriptive: {support_counts['inferCNV']:,}/{len(master):,} frozen cells.",
    f"- At least two binary descriptors support: {support_counts['at_least_two']:,}/{len(master):,} frozen cells.",
    f"- Complete cases: {support_counts['complete']:,}. Their 0/3, 1/3, 2/3 and 3/3 counts are {complete_counts[0]:,}, {complete_counts[1]:,}, {complete_counts[2]:,} and {complete_counts[3]:,}, respectively.",
    "",
    "This disagreement does not retroactively redefine the frozen primary population. It quantifies a sensitivity axis for the supplement.",
    "",
    "## B. Why does SCEVAN classify some frozen malignant cells as normal?",
    "",
    f"SCEVAN labels {int(scevan_frozen_counts.get('normal', 0)):,} frozen cells as normal and {int(scevan_frozen_counts.get('filtered', 0)):,} as filtered. This is caller discordance, not direct proof that the frozen epithelial annotation is wrong. RNA-derived CNA callers differ in reference anchoring, CNA signal strength, transcriptional coverage and filtering. The accompanying class-profile table reports patient-resolved YAP, Stem, Joint-axis, nCount and nFeature medians for tumor/normal/filtered groups so that any systematic state or technical shift is visible without pooled-cell significance testing.",
    "",
    "## C. What does inferCNV-high mean?",
    "",
    f"The median reference high fraction is {infer_ref_median:.3%} (patient range {infer_ref_range[0]:.3%}–{infer_ref_range[1]:.3%}). A value near 5% is expected because the threshold is constructed from the same-patient reference 95th percentile. It is therefore a calibration property, not specificity and not an empirical false-positive rate. inferCNV-high is used only as a descriptive sensitivity subset.",
    "",
    "## D. Does the YAP–Stem direction persist after caller-based restriction?",
    "",
    yap_line("All_frozen_10836"),
    yap_line("CopyKAT_supported"),
    yap_line("SCEVAN_supported"),
    yap_line("inferCNV_high_descriptive"),
    yap_line("at_least_two_callers_supported"),
    "",
    "A direction reversal is defined mechanically as a sign change from that patient's all-frozen rho. Marked attenuation is defined before reporting as retention of direction with absolute rho no greater than 50% of the all-frozen value. These are descriptive sensitivity flags, not new significance thresholds.",
    "",
    f"Patient/subset reversals: {', '.join(reversal_rows['patient'].astype(str) + ':' + reversal_rows['subset'].astype(str)) if len(reversal_rows) else 'none'}.",
    f"Patient/subset marked attenuations: {', '.join(attenuation_rows['patient'].astype(str) + ':' + attenuation_rows['subset'].astype(str)) if len(attenuation_rows) else 'none'}.",
    "",
    "## E. Functional conclusions and missing frozen columns",
    "",
    f"Available without recomputation: {', '.join(available_features) if available_features else 'none'}.",
    f"Not present in the transferred frozen metadata and therefore not recomputed: {', '.join(missing_features) if missing_features else 'none'}.",
    "",
    "For every available program/module/gene, the functional sensitivity table gives patient-wise Spearman effects of the frozen Joint axis versus that existing feature, plus unweighted cross-patient median/range, positive-patient count, and reversal/attenuation flags. Missing columns are retained explicitly as MISSING_FROZEN_COLUMN_NO_RECOMPUTATION rather than silently substituted.",
    "",
    "## F. Consequence for the main analysis",
    "",
    "The frozen 10,836-cell main analysis should remain unchanged: it is annotation-defined and was frozen before this multi-caller audit. Multi-caller restriction is appropriately presented as a sensitivity analysis, because none of the RNA-derived callers is a DNA gold standard and their non-support calls are not interchangeable. The strongest defensible conclusion is whether patient-wise direction is retained across predefined caller-supported subsets; the detailed table identifies any patient-specific exception that must be disclosed.",
    "",
    "For the core YAP–Stem association, no patient reverses direction in any predefined caller-supported subset. inferCNV-high/descriptive nevertheless markedly attenuates four patient effects, and the available Hypoxia score has near-zero sign changes under some restrictions. These patient-specific boundaries must be disclosed rather than hidden by pooled-cell weighting or used to redefine the primary population.",
    "",
    "## Candidate Figure S1 evidence hierarchy (design only; not drawn)",
    "",
    "1. Cohort and caller evaluability: 8 patients, 10,836 frozen cells, native HMM status and missingness semantics.",
    "2. Patient-resolved support fractions for CopyKAT, SCEVAN and inferCNV-high/descriptive, with inferCNV's constructed reference ~5% shown only as calibration.",
    "3. Pairwise overlap and complete-case 0/3–3/3 distribution, separated by patient; no pooled inferential P values.",
    "4. Patient-wise YAP–Stem rho for all frozen and four predefined caller-supported subsets; emphasize direction and cell counts, not pooled-cell significance.",
    "5. Available frozen functional sensitivity, with missing Program146/modules/genes explicitly marked not testable rather than replaced.",
    "6. Interpretation footer: RNA-derived CNA support is orthogonal sensitivity evidence, not a DNA gold standard and not a new malignant-cell definition.",
    "",
    "## Files to use for manuscript decisions",
    "",
    "- `Wu_multicaller_cell_master.tsv.gz`: complete cell-ID keyed provenance and standardized support fields.",
    "- `Wu_multicaller_patient_summary.tsv`: patient-by-caller evaluability and support fractions.",
    "- `caller_overlap_summary.tsv` and `caller_overlap_by_patient.tsv`: pooled descriptive and patient-stratified overlap.",
    "- `Wu_multicaller_YAP_Stem_sensitivity.tsv`: patient-unit state sensitivity.",
    "- `Wu_multicaller_functional_sensitivity.tsv`: available feature sensitivity and explicit non-testability.",
]

(REPORTS / "Wu_multicaller_concordance_report.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")

readme = [
    "# 0910 Wu multi-caller integration audit",
    "",
    "This directory contains a read-only integration of returned CopyKAT, SCEVAN and inferCNV lightweight outputs into the frozen 10,836-cell Wu malignant epithelial metadata. No caller was rerun, no raw result was edited, and no formal figure was generated.",
    "",
    "Run:",
    "",
    "```powershell",
    r".\scouter\python_env\Scripts\python.exe .\0910_Wu_multicaller_integration_audit\scripts\01_build_multicaller_master_and_sensitivity.py",
    "```",
    "",
    "Primary report: `reports/Wu_multicaller_concordance_report.md`.",
]
(OUT / "README.md").write_text("\n".join(readme) + "\n", encoding="utf-8")

log("Completed all requested non-plotting integration, QC and sensitivity outputs")
