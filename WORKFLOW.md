# Dissertation code workflow

The code is organised into the following execution order.

## 1. Dataset overview

`overview_dataset.py` is descriptive only. It inspects the raw baseline,
outcome, repeated-lifestyle, CVD-event, and imaging-date datasets. It does not
create an analysis cohort.

## 2. Data merge

`data_merge_new.py` calculates the lifestyle scores and writes:

- `data_merge/baseline.csv`: all baseline records, including both
  `def_CVD_AFTER` (primary endpoint) and `def_CVD_AF_HF_AFTER` (sensitivity
  composite endpoint). This is the input to `predictive_model.py`.
- `data_merge/longitudinal.csv`: participants with complete observed T0 and T2
  lifestyle exposures. This is the input to `EDA_longitudinal.py`.

## 3. Longitudinal EDA

`EDA_longitudinal.py` writes `EDA/longitudinal_after_eda.csv`.

Lifestyle exposure categories are not imputed because fractional or synthetic
categories would invalidate transition definitions. Only age, BMI, and sex are
imputed when required. `eid` is retained for the next merge.

## 4. Imaging-visit and CVD timing separation

`imaging_visit_separation.py` combines `longitudinal_after_eda.csv` with the
imaging-date and CVD-event files. It writes:

- `EDA/primary_causalML_after_eda.csv`: CVD between baseline and imaging, CVD
  after imaging, and no CVD; pre-baseline CVD is excluded.
- `EDA/sensitive_causalML_after_eda.csv`: CVD after imaging and no CVD only.
- `EDA/cvd_timing_audit.csv`: counts for every timing group, including excluded
  and unclassified participants.

The sensitivity definition has cleaner temporality because the CVD outcome
occurs after the T2 exposure measurement.

## 5. Cohort recoding

`cohort.py` independently applies the same lifestyle recoding to both files and
writes:

- `Cohort/primary_single_variable.csv`
- `Cohort/sensitive_single_variable.csv`
- `Cohort/primary_transition_summary.csv`

Both R causal-forest workflows use these recoded files. The combined-variable
workflow constructs paired treatments internally, so a separate combined CSV
is unnecessary.

## 6A. Causal-forest modelling

Run `single_variable3.R` and `combined_variable6.R` once for each cohort by
setting `DATA_PATH`, `ANALYSIS_LABEL`, and `OUTPUT_DIR`.

Primary example:

```bash
DATA_PATH=/home/rmhiund/causal_analysis/Cohort/primary_single_variable.csv \
ANALYSIS_LABEL=primary \
OUTPUT_DIR=/home/rmhiund/causal_analysis/results/single_primary \
Rscript single_variable3.R

DATA_PATH=/home/rmhiund/causal_analysis/Cohort/primary_single_variable.csv \
ANALYSIS_LABEL=primary \
OUTPUT_DIR=/home/rmhiund/causal_analysis/results/combined_primary \
Rscript combined_variable6.R
```

Sensitivity example:

```bash
DATA_PATH=/home/rmhiund/causal_analysis/Cohort/sensitive_single_variable.csv \
ANALYSIS_LABEL=sensitivity \
OUTPUT_DIR=/home/rmhiund/causal_analysis/results/single_sensitivity \
Rscript single_variable3.R

DATA_PATH=/home/rmhiund/causal_analysis/Cohort/sensitive_single_variable.csv \
ANALYSIS_LABEL=sensitivity \
OUTPUT_DIR=/home/rmhiund/causal_analysis/results/combined_sensitivity \
Rscript combined_variable6.R
```

## 6B. Baseline predictive modelling

`predictive_model.py` contains both preprocessing and modelling. Its default
input is `data_merge/baseline.csv`.

Primary endpoint:

```bash
OUTCOME=def_CVD_AFTER python predictive_model.py
```

Sensitivity composite endpoint:

```bash
OUTCOME=def_CVD_AF_HF_AFTER \
RESULT_DIR=/home/rmhiund/causal_analysis/predictive_model/results_composite \
python predictive_model.py
```

## Additional analysis files

- `unmeasured_confounder_refined.R` is a sensitivity analysis downstream of a
  recoded causal cohort.
- `baseline.R` is a separate legacy baseline causal-forest analysis. It is not
  part of the baseline predictive-model branch described above.

