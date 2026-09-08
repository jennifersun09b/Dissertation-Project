# Cardiovascular Risk Prediction and Longitudinal Lifestyle Change in UK Biobank

**Predictive and Causal Machine-Learning Models** — MSc Health Data Science dissertation, Institute of Health Informatics, UCL (Sep 2026).

The full dissertation is in [`docs/dissertation.pdf`](docs/dissertation.pdf). The deployed results web application lives in [`cvd_webapp/`](cvd_webapp/) (also maintained at [jennifersun09b/cvd-webapp](https://github.com/jennifersun09b/cvd-webapp)).

## Code overview


This repository contains the full analysis pipeline for the dissertation. The
workflow moves from raw UK Biobank source tables to two modelling branches, a
causal-forest branch built on the longitudinal lifestyle-change cohort and a
predictive-modelling branch built on the baseline cohort.

## Pipeline at a glance

```
overview_dataset.py            (descriptive only, no outputs)
        |
data_merge_new.py
        |-------------------------------------------.
        v                                            v
data_merge/longitudinal.csv                 data_merge/baseline.csv
        |                                            |
EDA_longitudinal.py                          predictive_model.py
        |                                    (preprocessing + modelling,
EDA/longitudinal_after_eda.csv                run once per outcome)
        |
imaging_visit_separation.py
  + imaging_visit_date.tsv
  + bhf_all_individuals_plus_cvd_events.csv
        |
        |-- EDA/primary_causalML_after_eda.csv
        |-- EDA/sensitive_causalML_after_eda.csv
        |-- EDA/cvd_timing_audit.csv
        |
cohort.py   (lifestyle recoding, applied to both cohorts)
        |
        |-- Cohort/primary_single_variable.csv
        |-- Cohort/sensitive_single_variable.csv
        |-- Cohort/primary_transition_summary.csv
        |
single_variable3.R  and  combined_variable6.R
  (each run once per cohort: primary and sensitivity)
```

## Step 1 — Dataset overview

`overview_dataset.py` inspects the five raw source tables (baseline
covariates, outcomes, repeated lifestyle assessments, CVD event timing, and
imaging visit dates). It is deliberately descriptive only. It prints shapes,
dtypes, missingness, and key value counts, and it does not merge, filter, or
write any analysis dataset. The imaging-visit and CVD-timing separation logic
lives in `imaging_visit_separation.py` (step 4).

```bash
python overview_dataset.py
```

## Step 2 — Data merge

`data_merge_new.py` computes the six lifestyle scores (sleep, smoking,
alcohol, diet, physical activity, mental health) at baseline (instance 0.0,
T0) and at the imaging visit (instance 2.0, T2), merges them with the
baseline covariates and both outcomes, and writes two datasets.

- `data_merge/baseline.csv` holds every merged participant with baseline
  covariates, the T0/T2 lifestyle variables, and both outcome columns,
  `def_CVD_AFTER` (primary endpoint) and `def_CVD_AF_HF_AFTER` (sensitivity
  composite of CVD, atrial fibrillation, and heart failure). This file feeds
  the predictive-modelling branch.
- `data_merge/longitudinal.csv` keeps only participants with complete
  observed lifestyle exposures at both T0 and T2. This file feeds the
  causal-forest branch.

```bash
python data_merge_new.py
```

## Step 3 — Longitudinal EDA

`EDA_longitudinal.py` cleans `longitudinal.csv` and writes
`EDA/longitudinal_after_eda.csv`. Lifestyle exposure categories are never
imputed, because fractional or synthetic categories would invalidate the
transition definitions used later. Only age, BMI, and sex are imputed when
required, and `eid` is retained for the imaging/CVD-timing merge in step 4.

```bash
python EDA_longitudinal.py
```

## Step 4 — Imaging-visit and CVD-timing separation

`imaging_visit_separation.py` merges `longitudinal_after_eda.csv` with the
imaging-visit-date file and the CVD event-timing file, assigns every
participant a CVD timing group relative to baseline (T0) and the imaging
visit (T2), and splits the cohort.

- `EDA/primary_causalML_after_eda.csv` — the primary cohort: no pre-baseline
  CVD, keeping participants with CVD between baseline and imaging, CVD after
  imaging, or no CVD.
- `EDA/sensitive_causalML_after_eda.csv` — the sensitivity cohort: CVD after
  imaging or no CVD only. This definition has cleaner temporality because
  every outcome occurs after the T2 exposure measurement.
- `EDA/cvd_timing_audit.csv` — counts for every timing group, including
  excluded and unclassified participants.

```bash
python imaging_visit_separation.py
```

## Step 5 — Cohort recoding

`cohort.py` applies the same lifestyle recoding independently to the primary
and sensitivity cohorts (collapsing sparse categories and building the
modified T2 transition variables, `*_t1_mod`) and writes:

- `Cohort/primary_single_variable.csv`
- `Cohort/sensitive_single_variable.csv`
- `Cohort/primary_transition_summary.csv`

Both recoded files serve the single-variable and the combined-variable
causal-forest scripts. The combined-variable workflow constructs its paired
treatments internally, so no separate combined CSV is needed.

```bash
python cohort.py
```

## Step 6A — Causal-forest modelling (longitudinal branch)

`single_variable3.R` estimates per-domain causal forests and
`combined_variable6.R` estimates paired-treatment causal forests. Each script
runs once per cohort, controlled by `DATA_PATH`, `ANALYSIS_LABEL`, and
`OUTPUT_DIR`.

```bash
# Primary cohort
DATA_PATH=.../Cohort/primary_single_variable.csv ANALYSIS_LABEL=primary \
  OUTPUT_DIR=.../results/single_primary Rscript single_variable3.R
DATA_PATH=.../Cohort/primary_single_variable.csv ANALYSIS_LABEL=primary \
  OUTPUT_DIR=.../results/combined_primary Rscript combined_variable6.R

# Sensitivity cohort
DATA_PATH=.../Cohort/sensitive_single_variable.csv ANALYSIS_LABEL=sensitivity \
  OUTPUT_DIR=.../results/single_sensitivity Rscript single_variable3.R
DATA_PATH=.../Cohort/sensitive_single_variable.csv ANALYSIS_LABEL=sensitivity \
  OUTPUT_DIR=.../results/combined_sensitivity Rscript combined_variable6.R
```

## Step 6B — Baseline predictive modelling

`predictive_model.py` contains both the preprocessing and the modelling for
the baseline branch. Its default input is `data_merge/baseline.csv`, and the
outcome is selected with the `OUTCOME` environment variable.

```bash
# Primary endpoint
OUTCOME=def_CVD_AFTER python predictive_model.py

# Sensitivity composite endpoint
OUTCOME=def_CVD_AF_HF_AFTER RESULT_DIR=.../results_composite python predictive_model.py
```

## Additional analyses

- `unmeasured_confounder_refined.R` — unmeasured-confounding sensitivity
  analysis, run downstream of a recoded causal cohort.
- `baseline.R` — a legacy baseline causal-forest analysis. It is not part of
  the baseline predictive-model branch in step 6B.
- `cvd_webapp/` — the results web application, which consumes the model
  artifacts exported by the modelling scripts.

## Path configuration

Every Python script reads its input and output locations from environment
variables (`BASELINE_SOURCE_PATH`, `OUTCOME_SOURCE_PATH`,
`LIFESTYLE_SOURCE_PATH`, `EVENT_SOURCE_PATH`, `IMAGING_DATE_PATH`,
`INPUT_PATH`, `OUTPUT_PATH`, `OUTPUT_DIR`, `DATA_PATH`, `RESULT_DIR`,
`OUTCOME`), falling back to the HPC paths under
`/home/rmhiund/causal_analysis/`. Override them locally as needed, for
example:

```bash
OUTPUT_DIR=./data_merge python data_merge_new.py
```
