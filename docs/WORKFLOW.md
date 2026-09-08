# Workflow runbook

Execution order, exact commands, and environment variables for reproducing the
dissertation analyses. The overview and diagrams live in the top-level
[README](../README.md); this file is the operational reference.

All paths below use the HPC layout under `/home/rmhiund/causal_analysis/`.
Override any path with the environment variables listed at the end.

## 1. Dataset overview

```bash
python src/data_preparation/01_overview_dataset.py
```

Descriptive only. Inspects the raw baseline, outcome, repeated-lifestyle,
CVD-event, and imaging-date tables and prints shapes, dtypes, missingness, and
key value counts. It does not merge, filter, or write an analysis dataset.

## 2. Data merge

```bash
python src/data_preparation/02_data_merge.py
```

Computes the six lifestyle scores (sleep, smoking, alcohol, diet, physical
activity, mental health) at baseline (instance 0.0, T0) and the imaging visit
(instance 2.0, T2), merges them with baseline covariates and both outcomes,
and writes:

- `data_merge/baseline.csv`: every merged participant, with `def_CVD_AFTER`
  (primary endpoint) and `def_CVD_AF_HF_AFTER` (composite of CVD, atrial
  fibrillation, heart failure). Input to the predictive model.
- `data_merge/longitudinal.csv`: participants with complete observed lifestyle
  exposures at both T0 and T2. Input to the causal branch.

## 3. Longitudinal EDA

```bash
python src/data_preparation/03_eda_longitudinal.py
```

Writes `EDA/longitudinal_after_eda.csv`. Lifestyle exposure categories are not
imputed because fractional or synthetic categories would invalidate the
transition definitions. Only age, BMI, and sex are imputed when required.
`eid` is retained for the next merge.

## 4. Imaging-visit and CVD-timing separation

```bash
python src/data_preparation/04_imaging_visit_separation.py
```

Merges `longitudinal_after_eda.csv` with the imaging-date and CVD-event files,
assigns every participant a CVD timing group relative to T0 and T2, and
writes:

- `EDA/primary_causalML_after_eda.csv`: CVD between baseline and imaging, CVD
  after imaging, and no CVD. Pre-baseline CVD is excluded.
- `EDA/sensitive_causalML_after_eda.csv`: CVD after imaging and no CVD only.
  Cleaner temporality, because every outcome follows the T2 exposure.
- `EDA/cvd_timing_audit.csv`: counts for every timing group, including
  excluded and unclassified participants.

## 5. Cohort recoding

```bash
python src/data_preparation/05_cohort_recoding.py
```

Applies the same lifestyle recoding to both cohorts (collapsing sparse
categories and building the modified T2 transition variables `*_t1_mod`) and
writes:

- `Cohort/primary_single_variable.csv`
- `Cohort/sensitive_single_variable.csv`
- `Cohort/primary_transition_summary.csv`

Both R workflows read these files. The combined-variable workflow builds its
paired treatments internally, so no separate combined CSV is needed.

## 6A. Baseline predictive model

`src/predictive_model/predictive_model.py` contains preprocessing and
modelling. Default input is `data_merge/baseline.csv`; the outcome is chosen
with `OUTCOME`.

```bash
# Primary endpoint
OUTCOME=def_CVD_AFTER python src/predictive_model/predictive_model.py

# Sensitivity composite endpoint
OUTCOME=def_CVD_AF_HF_AFTER \
RESULT_DIR=/home/rmhiund/causal_analysis/predictive_model/results_composite \
python src/predictive_model/predictive_model.py
```

Set `RUN_HEAVY=1` for the full bootstrap and tuning budget. The logistic
pipeline is exported to `predictive_model.joblib` for the web app.

## 6B. Causal forests

Run `single_variable.R` and `combined_variable.R` once per cohort by setting
`DATA_PATH`, `ANALYSIS_LABEL`, and `OUTPUT_DIR`.

```bash
# Primary cohort
DATA_PATH=/home/rmhiund/causal_analysis/Cohort/primary_single_variable.csv \
ANALYSIS_LABEL=primary \
OUTPUT_DIR=/home/rmhiund/causal_analysis/results/single_primary \
Rscript src/causal_forest/single_variable.R

DATA_PATH=/home/rmhiund/causal_analysis/Cohort/primary_single_variable.csv \
ANALYSIS_LABEL=primary \
OUTPUT_DIR=/home/rmhiund/causal_analysis/results/combined_primary \
Rscript src/causal_forest/combined_variable.R

# Sensitivity cohort
DATA_PATH=/home/rmhiund/causal_analysis/Cohort/sensitive_single_variable.csv \
ANALYSIS_LABEL=sensitivity \
OUTPUT_DIR=/home/rmhiund/causal_analysis/results/single_sensitivity \
Rscript src/causal_forest/single_variable.R

DATA_PATH=/home/rmhiund/causal_analysis/Cohort/sensitive_single_variable.csv \
ANALYSIS_LABEL=sensitivity \
OUTPUT_DIR=/home/rmhiund/causal_analysis/results/combined_sensitivity \
Rscript src/causal_forest/combined_variable.R
```

Both scripts write the dissertation CSV/RDS outputs and, via `reticulate`,
the joblib tables used by the web app (`single_variable_ate.joblib`,
`combined_variable_ate.joblib`).

## 6C. Additional causal analyses

- `src/causal_forest/unmeasured_confounding.R`: unmeasured-confounding
  sensitivity analysis (E-values and quantitative bias analysis). Refits the
  same forests as `single_variable.R` with the same seed and settings. Run
  downstream of a recoded cohort with `DATA_PATH` and `OUTPUT_DIR`; restrict
  with `CONTRAST_DOMAINS`, and speed up with `N_BOOT` and `NUM_TREES`.
- `src/causal_forest/baseline_reference.R`: the baseline causal reference
  model, one causal forest per lifestyle domain at baseline. It is not part
  of the predictive branch in 6A.

## 7. Web app

```bash
cd cvd_webapp
pip install -r requirements.txt
streamlit run app.py
```

Copy the three joblib artefacts from 6A and 6B into `cvd_webapp/model/`.

## Environment variables

| Variable | Meaning |
| --- | --- |
| `BASELINE_SOURCE_PATH` | raw baseline covariate table |
| `OUTCOME_SOURCE_PATH` | raw outcome table |
| `LIFESTYLE_SOURCE_PATH` | raw repeated lifestyle assessments |
| `EVENT_SOURCE_PATH` | CVD event-timing file |
| `IMAGING_DATE_PATH` | imaging visit dates |
| `INPUT_PATH`, `OUTPUT_PATH`, `OUTPUT_DIR` | per-step input file and output location |
| `DATA_PATH` | analysis dataset for the modelling scripts |
| `RESULT_DIR` | output folder for `predictive_model.py` |
| `OUTCOME` | `def_CVD_AFTER` (default) or `def_CVD_AF_HF_AFTER` |
| `ANALYSIS_LABEL` | `primary` or `sensitivity`, used in R output names |
| `RUN_HEAVY`, `N_JOBS` | compute budget for the predictive model |
