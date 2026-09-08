#!/usr/bin/env python3
"""Assign CVD timing relative to baseline and the imaging visit.

Workflow position:
    longitudinal_after_eda.csv
        + imaging_visit_date.tsv
        + bhf_all_individuals_plus_cvd_events.csv
        -> primary_causalML_after_eda.csv
        -> sensitive_causalML_after_eda.csv

Primary cohort:
  participants without pre-baseline CVD whose post-baseline status is either
  CVD between baseline and imaging, CVD after imaging, or no CVD.

Sensitivity cohort:
  participants without pre-baseline CVD whose status is CVD after imaging or
  no CVD. This removes outcomes occurring before the exposure period ended.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd


LONGITUDINAL_EDA_PATH = Path(os.environ.get(
    "LONGITUDINAL_EDA_PATH",
    "/home/rmhiund/causal_analysis/EDA/longitudinal_after_eda.csv",
))
IMAGING_DATE_PATH = Path(os.environ.get(
    "IMAGING_DATE_PATH",
    "/home/rmhiund/causal_analysis/Dataset/imaging_visit_date.tsv",
))
EVENT_SOURCE_PATH = Path(os.environ.get(
    "EVENT_SOURCE_PATH",
    "/home/rmhiund/causal_analysis/Dataset/bhf_all_individuals_plus_cvd_events.csv",
))
OUTPUT_DIR = Path(os.environ.get(
    "OUTPUT_DIR",
    "/home/rmhiund/causal_analysis/EDA",
))


def read_table(path: Path) -> pd.DataFrame:
    separator = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=separator, low_memory=False)


def require_columns(frame: pd.DataFrame, columns: list[str], name: str) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f"{name} is missing required columns: {missing}")


def coerce_binary(series: pd.Series, name: str) -> pd.Series:
    """Convert common logical encodings to nullable 0/1 integers."""
    mapping = {
        True: 1, False: 0, 1: 1, 0: 0, 1.0: 1, 0.0: 0,
        "true": 1, "false": 0, "1": 1, "0": 0,
        "yes": 1, "no": 0,
    }
    normalised = series.map(
        lambda value: mapping.get(value.strip().lower(), np.nan)
        if isinstance(value, str) else mapping.get(value, np.nan)
    )
    invalid = series.notna() & normalised.isna()
    if invalid.any():
        examples = series.loc[invalid].astype(str).unique()[:10]
        raise ValueError(f"{name} contains unsupported values: {examples.tolist()}")
    return normalised.astype("Int64")


def main() -> None:
    longitudinal = read_table(LONGITUDINAL_EDA_PATH)
    imaging = read_table(IMAGING_DATE_PATH)
    events = read_table(EVENT_SOURCE_PATH)

    require_columns(longitudinal, ["eid"], "longitudinal EDA dataset")
    require_columns(imaging, ["eid", "53-2.0"], "imaging-date dataset")
    require_columns(
        events,
        [
            "eid", "defined_baseline_date", "def_CVD_BEFORE",
            "def_CVD_AFTER", "def_CVD_AFTER_days_from_baseline",
        ],
        "CVD event dataset",
    )

    for name, frame in {
        "longitudinal": longitudinal,
        "imaging": imaging,
        "events": events,
    }.items():
        if frame["eid"].duplicated().any():
            raise ValueError(f"{name} dataset contains duplicated eid values")

    imaging = imaging[["eid", "53-2.0"]].copy()
    imaging["imaging_date"] = pd.to_datetime(imaging.pop("53-2.0"), errors="coerce")

    event_columns = [
        "eid", "defined_baseline_date", "def_CVD_BEFORE", "def_CVD_AFTER",
        "def_CVD_AFTER_days_from_baseline",
    ]
    if "def_CVD_AF_HF_AFTER" in events.columns:
        event_columns.append("def_CVD_AF_HF_AFTER")
    events = events[event_columns].copy()
    events["defined_baseline_date"] = pd.to_datetime(
        events["defined_baseline_date"], errors="coerce"
    )
    events["def_CVD_BEFORE"] = coerce_binary(
        events["def_CVD_BEFORE"], "def_CVD_BEFORE"
    )
    events["def_CVD_AFTER"] = coerce_binary(
        events["def_CVD_AFTER"], "def_CVD_AFTER"
    )
    events["def_CVD_AFTER_days_from_baseline"] = pd.to_numeric(
        events["def_CVD_AFTER_days_from_baseline"], errors="coerce"
    )
    events["CVD_event_date_after_baseline"] = (
        events["defined_baseline_date"]
        + pd.to_timedelta(events["def_CVD_AFTER_days_from_baseline"], unit="D")
    )

    timing = events.merge(imaging, on="eid", how="left", validate="one_to_one")
    combined = longitudinal.merge(
        timing, on="eid", how="left", validate="one_to_one", suffixes=("", "_event")
    )

    # Prefer the authoritative event table for the primary outcome, while
    # retaining a mismatch audit if the EDA input already carried the outcome.
    if "CVD_outcome" in combined.columns:
        existing = pd.to_numeric(combined["CVD_outcome"], errors="coerce")
        mismatch = existing.notna() & combined["def_CVD_AFTER"].notna() & (
            existing != combined["def_CVD_AFTER"].astype(float)
        )
        if mismatch.any():
            raise ValueError(
                f"CVD outcome disagrees with the event table for {mismatch.sum()} participants"
            )
    combined["CVD_outcome"] = combined["def_CVD_AFTER"].astype("Int64")

    combined["cvd_timing_group"] = pd.NA
    before = combined["def_CVD_BEFORE"] == 1
    after = combined["def_CVD_AFTER"] == 1
    no_after = combined["def_CVD_AFTER"] == 0
    dated_event = combined["CVD_event_date_after_baseline"].notna()
    dated_imaging = combined["imaging_date"].notna()

    combined.loc[before & after, "cvd_timing_group"] = "0_CVD_always"
    combined.loc[before & no_after, "cvd_timing_group"] = "1_CVD_before_baseline"
    combined.loc[
        ~before & after & dated_event & dated_imaging
        & (combined["CVD_event_date_after_baseline"] < combined["imaging_date"]),
        "cvd_timing_group",
    ] = "2_CVD_between_baseline_and_imaging"
    combined.loc[
        ~before & after & dated_event & dated_imaging
        & (combined["CVD_event_date_after_baseline"] >= combined["imaging_date"]),
        "cvd_timing_group",
    ] = "3_CVD_after_imaging"
    combined.loc[
        ~before & no_after & dated_imaging,
        "cvd_timing_group",
    ] = "4_no_CVD"

    primary_groups = {
        "2_CVD_between_baseline_and_imaging",
        "3_CVD_after_imaging",
        "4_no_CVD",
    }
    sensitivity_groups = {"3_CVD_after_imaging", "4_no_CVD"}
    primary = combined[combined["cvd_timing_group"].isin(primary_groups)].copy()
    sensitivity = combined[
        combined["cvd_timing_group"].isin(sensitivity_groups)
    ].copy()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    primary_path = OUTPUT_DIR / "primary_causalML_after_eda.csv"
    sensitivity_path = OUTPUT_DIR / "sensitive_causalML_after_eda.csv"
    audit_path = OUTPUT_DIR / "cvd_timing_audit.csv"

    primary.to_csv(primary_path, index=False)
    sensitivity.to_csv(sensitivity_path, index=False)
    (
        combined["cvd_timing_group"]
        .value_counts(dropna=False)
        .rename_axis("cvd_timing_group")
        .reset_index(name="n")
        .to_csv(audit_path, index=False)
    )

    print("CVD timing groups:")
    print(combined["cvd_timing_group"].value_counts(dropna=False).to_string())
    print(f"primary cohort: {len(primary):,} -> {primary_path}")
    print(f"sensitivity cohort: {len(sensitivity):,} -> {sensitivity_path}")
    print(f"timing audit -> {audit_path}")


if __name__ == "__main__":
    main()
