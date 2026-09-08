#!/usr/bin/env python3
"""Recode primary and sensitivity cohorts for causal-forest analyses.

Workflow position:
    primary_causalML_after_eda.csv   -> primary_single_variable.csv
    sensitive_causalML_after_eda.csv -> sensitive_single_variable.csv

The same recoding function is applied independently to both cohorts. Both the
single-variable and combined-variable causal-forest scripts consume these
outputs; the combined R workflow builds paired treatments internally.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd


PRIMARY_INPUT_PATH = Path(os.environ.get(
    "PRIMARY_INPUT_PATH",
    "/home/rmhiund/causal_analysis/EDA/primary_causalML_after_eda.csv",
))
SENSITIVITY_INPUT_PATH = Path(os.environ.get(
    "SENSITIVITY_INPUT_PATH",
    "/home/rmhiund/causal_analysis/EDA/sensitive_causalML_after_eda.csv",
))
OUTPUT_DIR = Path(os.environ.get(
    "OUTPUT_DIR",
    "/home/rmhiund/causal_analysis/Cohort",
))

RECODE_MAPS = {
    "mental": {0.0: 0.0, 1.0: 1.0, 2.0: 2.0, 3.0: 2.0},
    "alcohol": {0.0: 0.0, 1.0: 0.0, 2.0: 1.0, 3.0: 2.0, 4.0: 2.0},
    "diet": {0.0: 0.0, 1.0: 0.0, 2.0: 1.0, 3.0: 2.0},
}
DOMAIN_COLUMNS = {
    "sleep": ("sleep_category_0.0", "sleep_category_2.0", "sleep_t1_mod"),
    "smoking": ("smoking_category_0.0", "smoking_category_2.0", "smoking_t1_mod"),
    "alcohol": ("alcohol_score_0.0", "alcohol_score_2.0", "alcohol_t1_mod"),
    "diet": ("diet_score_0.0", "diet_score_2.0", "diet_t1_mod"),
    "physical": ("physical_category_0.0", "physical_category_2.0", "pa_t1_mod"),
    "mental": ("mental_score_0.0", "mental_score_2.0", "mental_t1_mod"),
}


def modify_t1(t0: pd.Series, t1: pd.Series) -> pd.Series:
    """Collapse extreme transitions while retaining exact moderate moves."""
    result = t1.astype(float).copy()
    result[(t0 == 0) & t1.isin([1, 2])] = 1.0
    result[(t0 == 2) & t1.isin([0, 1])] = 1.0
    result[t0.isna() | t1.isna()] = np.nan
    return result


def require_columns(frame: pd.DataFrame, path: Path) -> None:
    required = {"eid", "age", "sex", "bmi", "CVD_outcome", "cvd_timing_group"}
    for t0, t2, _ in DOMAIN_COLUMNS.values():
        required.update([t0, t2])
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {missing}")


def prepare_cohort(path: Path, expected_groups: set[str]) -> pd.DataFrame:
    frame = pd.read_csv(path, low_memory=False)
    require_columns(frame, path)
    if frame["eid"].duplicated().any():
        raise ValueError(f"{path} contains duplicated eid values")
    unexpected = set(frame["cvd_timing_group"].dropna().unique()) - expected_groups
    if unexpected:
        raise ValueError(f"{path} contains unexpected timing groups: {sorted(unexpected)}")

    frame[["mental_score_0.0", "mental_score_2.0"]] = frame[
        ["mental_score_0.0", "mental_score_2.0"]
    ].replace(RECODE_MAPS["mental"])
    frame[["alcohol_score_0.0", "alcohol_score_2.0"]] = frame[
        ["alcohol_score_0.0", "alcohol_score_2.0"]
    ].replace(RECODE_MAPS["alcohol"])
    frame[["diet_score_0.0", "diet_score_2.0"]] = frame[
        ["diet_score_0.0", "diet_score_2.0"]
    ].replace(RECODE_MAPS["diet"])

    for _, (t0, t2, modified) in DOMAIN_COLUMNS.items():
        frame[modified] = modify_t1(frame[t0], frame[t2])

    modelling_columns = ["age", "sex", "bmi", "CVD_outcome"]
    for t0, _, modified in DOMAIN_COLUMNS.values():
        modelling_columns.extend([t0, modified])
    frame[modelling_columns] = frame[modelling_columns].apply(
        pd.to_numeric, errors="raise"
    )
    if frame[modelling_columns].isna().any().any():
        missing = frame[modelling_columns].isna().sum()
        raise ValueError(
            "causal-forest variables still contain missing values:\n"
            + missing[missing > 0].to_string()
        )
    if not set(frame["CVD_outcome"].unique()).issubset({0, 1, 0.0, 1.0}):
        raise ValueError("CVD_outcome must be coded 0/1")
    return frame


def transition_table(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    healthier_when_higher = {"sleep", "diet", "physical"}
    for domain, (t0, _, modified) in DOMAIN_COLUMNS.items():
        for (baseline, destination), subset in frame.groupby([t0, modified]):
            if baseline == destination:
                transition = "maintained"
            else:
                healthier = (
                    destination > baseline
                    if domain in healthier_when_higher
                    else destination < baseline
                )
                transition = "improved" if healthier else "deteriorated"
            rows.append({
                "domain": domain,
                "baseline": baseline,
                "destination": destination,
                "transition": transition,
                "n": len(subset),
                "events": int(subset["CVD_outcome"].sum()),
                "event_rate_percent": subset["CVD_outcome"].mean() * 100,
            })
    return pd.DataFrame(rows).sort_values(["domain", "baseline", "destination"])


def main() -> None:
    primary_groups = {
        "2_CVD_between_baseline_and_imaging",
        "3_CVD_after_imaging",
        "4_no_CVD",
    }
    sensitivity_groups = {"3_CVD_after_imaging", "4_no_CVD"}
    primary = prepare_cohort(PRIMARY_INPUT_PATH, primary_groups)
    sensitivity = prepare_cohort(SENSITIVITY_INPUT_PATH, sensitivity_groups)

    if not set(sensitivity["eid"]).issubset(set(primary["eid"])):
        raise ValueError("sensitivity cohort is not a subset of the primary cohort")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    primary_path = OUTPUT_DIR / "primary_single_variable.csv"
    sensitivity_path = OUTPUT_DIR / "sensitive_single_variable.csv"
    transition_path = OUTPUT_DIR / "primary_transition_summary.csv"

    primary.to_csv(primary_path, index=False)
    sensitivity.to_csv(sensitivity_path, index=False)
    transition_table(primary).to_csv(transition_path, index=False)

    print(f"primary: {len(primary):,} -> {primary_path}")
    print(f"sensitivity: {len(sensitivity):,} -> {sensitivity_path}")
    print(f"transition summary -> {transition_path}")


if __name__ == "__main__":
    main()
