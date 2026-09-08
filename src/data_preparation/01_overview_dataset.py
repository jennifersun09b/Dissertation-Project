#!/usr/bin/env python3
"""Inspect the source datasets used by the dissertation workflow.

This file is deliberately descriptive only: it does not merge cohorts, assign
CVD timing groups, filter participants, or write analysis datasets. The imaging
visit and CVD timing logic lives in ``04_imaging_visit_separation.py``.

Environment variables may override every input path:
  BASELINE_SOURCE_PATH, OUTCOME_SOURCE_PATH, LIFESTYLE_SOURCE_PATH,
  EVENT_SOURCE_PATH, IMAGING_DATE_PATH
"""

import os
from pathlib import Path

import pandas as pd


BASELINE_SOURCE_PATH = Path(os.environ.get(
    "BASELINE_SOURCE_PATH",
    "/home/rmhiund/causal_analysis/Dataset/group_1_clean.tsv",
))
OUTCOME_SOURCE_PATH = Path(os.environ.get(
    "OUTCOME_SOURCE_PATH",
    "/home/rmhiund/causal_analysis/Dataset/group_1_outcomes_df_without_qc_df.tsv",
))
LIFESTYLE_SOURCE_PATH = Path(os.environ.get(
    "LIFESTYLE_SOURCE_PATH",
    "/myriadfs/projects/ICS_UKB/projects/kasia_project/causal_analysis/ukb_tabular_data_causal_analysis.tsv",
))
EVENT_SOURCE_PATH = Path(os.environ.get(
    "EVENT_SOURCE_PATH",
    "/home/rmhiund/causal_analysis/Dataset/bhf_all_individuals_plus_cvd_events.csv",
))
IMAGING_DATE_PATH = Path(os.environ.get(
    "IMAGING_DATE_PATH",
    "/home/rmhiund/causal_analysis/Dataset/imaging_visit_date.tsv",
))


def load_table(path: Path) -> pd.DataFrame:
    """Load a CSV/TSV using the delimiter implied by its suffix."""
    separator = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=separator, low_memory=False)


def describe_table(name: str, frame: pd.DataFrame) -> None:
    """Print a compact, reproducible dataset overview."""
    print(f"\n{'=' * 80}\n{name}\n{'=' * 80}")
    print(f"shape: {frame.shape[0]:,} rows x {frame.shape[1]:,} columns")
    print("first columns:", list(frame.columns[:30]))
    print("\ndtypes:")
    print(frame.dtypes.value_counts().to_string())
    missing = (frame.isna().mean() * 100).sort_values(ascending=False)
    print("\ncolumns with the most missing data (%):")
    print(missing.head(20).round(2).to_string())
    if "eid" in frame.columns:
        print(f"\neid unique: {frame['eid'].nunique(dropna=True):,}")
        print(f"eid duplicated rows: {frame['eid'].duplicated().sum():,}")


def main() -> None:
    sources = {
        "Baseline population and covariates": BASELINE_SOURCE_PATH,
        "Outcome table": OUTCOME_SOURCE_PATH,
        "Raw repeated lifestyle assessments": LIFESTYLE_SOURCE_PATH,
        "CVD event timing table": EVENT_SOURCE_PATH,
        "Imaging visit dates": IMAGING_DATE_PATH,
    }

    for name, path in sources.items():
        print(f"\nLoading {name}: {path}")
        frame = load_table(path)
        describe_table(name, frame)

        if name == "Baseline population and covariates":
            # Notes retained from the original notebook:
            # - eid is the participant identifier.
            # - genetic_sex: 0 = female, 1 = male.
            # - merged_bmi is the merged baseline BMI measure.
            # - CRP, creatinine, HDL, total cholesterol, HbA1c, SBP, DBP and
            #   LDL are candidate baseline risk-prediction covariates.
            if "sepsis" in frame.columns:
                print("\nsepsis value counts:")
                print(frame["sepsis"].value_counts(dropna=False).to_string())

        if name == "Outcome table":
            # CVD_AF_HF denotes the composite cardiovascular disease, atrial
            # fibrillation, or heart-failure endpoint. Keep the primary CVD
            # endpoint and this composite endpoint as separate outcomes.
            for outcome in ("def_CVD_AFTER", "def_CVD_AF_HF_AFTER"):
                if outcome in frame.columns:
                    print(f"\n{outcome} value counts:")
                    print(frame[outcome].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
