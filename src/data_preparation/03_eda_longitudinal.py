#!/usr/bin/env python3
"""Clean the longitudinal dataset without changing lifestyle exposure classes.

Workflow position:
    longitudinal.csv -> longitudinal_after_eda.csv

The six lifestyle variables at baseline (T0) and imaging (T2) define the causal
exposures. They must remain exact observed categories; imputing them with MICE
can create fractional transition values. Therefore this script:
  - requires all T0/T2 lifestyle categories to be observed;
  - imputes only age, BMI, and sex when needed;
  - retains ``eid`` for the later imaging/CVD timing merge;
  - does not perform a train/test split, because this is cohort preparation,
    not predictive-model evaluation.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer


INPUT_PATH = Path(os.environ.get(
    "INPUT_PATH",
    "/home/rmhiund/causal_analysis/data_merge/longitudinal.csv",
))
OUTPUT_PATH = Path(os.environ.get(
    "OUTPUT_PATH",
    "/home/rmhiund/causal_analysis/EDA/longitudinal_after_eda.csv",
))

EXPOSURE_COLUMNS = [
    "sleep_category_0.0", "sleep_category_2.0",
    "smoking_category_0.0", "smoking_category_2.0",
    "alcohol_score_0.0", "alcohol_score_2.0",
    "diet_score_0.0", "diet_score_2.0",
    "physical_category_0.0", "physical_category_2.0",
    "mental_score_0.0", "mental_score_2.0",
]
IDENTIFIER_AND_COVARIATES = [
    "eid", "genetic_sex", "age_defined_baseline", "merged_bmi",
]


def require_columns(frame: pd.DataFrame, columns: list[str]) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f"longitudinal dataset is missing required columns: {missing}")


def validate_categories(frame: pd.DataFrame) -> None:
    """Ensure exposure values are integral observed categories."""
    for column in EXPOSURE_COLUMNS:
        numeric = pd.to_numeric(frame[column], errors="coerce")
        invalid_numeric = frame[column].notna() & numeric.isna()
        if invalid_numeric.any():
            raise ValueError(f"{column} contains non-numeric category values")
        non_integral = numeric.notna() & ~np.isclose(numeric, np.round(numeric))
        if non_integral.any():
            examples = numeric.loc[non_integral].unique()[:10]
            raise ValueError(
                f"{column} contains fractional categories: {examples.tolist()}"
            )
        frame[column] = numeric


def main() -> None:
    data = pd.read_csv(INPUT_PATH, low_memory=False)
    require_columns(data, IDENTIFIER_AND_COVARIATES + EXPOSURE_COLUMNS)

    if data["eid"].duplicated().any():
        raise ValueError("longitudinal dataset contains duplicated eid values")

    data = data[IDENTIFIER_AND_COVARIATES + EXPOSURE_COLUMNS].copy()
    validate_categories(data)

    missing_exposure = data[EXPOSURE_COLUMNS].isna().any(axis=1)
    if missing_exposure.any():
        raise ValueError(
            f"{missing_exposure.sum()} rows have a missing T0/T2 lifestyle category. "
            "The longitudinal merge should require complete observed exposures."
        )

    data = data.rename(columns={
        "genetic_sex": "sex",
        "age_defined_baseline": "age",
        "merged_bmi": "bmi",
    })

    data[["age", "bmi"]] = SimpleImputer(strategy="median").fit_transform(
        data[["age", "bmi"]]
    )
    data[["sex"]] = SimpleImputer(strategy="most_frequent").fit_transform(
        data[["sex"]]
    )

    for column in ["age", "bmi", "sex"] + EXPOSURE_COLUMNS:
        data[column] = pd.to_numeric(data[column], errors="raise")

    if data[["age", "bmi", "sex"] + EXPOSURE_COLUMNS].isna().any().any():
        raise ValueError("EDA output still contains missing modelling variables")

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(OUTPUT_PATH, index=False)

    print(f"input: {INPUT_PATH}")
    print(f"rows: {len(data):,}; columns: {data.shape[1]}")
    print("missing values after cleaning:")
    print(data.isna().sum().to_string())
    print(f"wrote: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
