# Causal and Predictive Machine Learning for Cardiovascular Risk in UK Biobank

Code for an MSc Health Data Science dissertation (UCL, 2026) that combines two
questions in one cohort: how well can lifestyle and demographic information
predict incident cardiovascular disease (CVD), and does changing a lifestyle
factor causally alter that risk. The pipeline runs from raw UK Biobank tables
to a baseline prediction model, causal-forest estimates of lifestyle change,
external validation, and a Streamlit app that presents both.

## Key results

- **Prediction.** Logistic regression, random forest, XGBoost and LightGBM were
  compared on 458,840 participants using age, sex, BMI and six lifestyle
  scores. All reached AUC 0.714 to 0.718 with good calibration; logistic
  regression was kept for interpretability (AUC 0.718, 95% CI 0.712 to 0.724).
- **Causal effects.** In 71,428 participants with lifestyle measured at two
  visits, doubly robust causal forests estimated 24 single-domain lifestyle
  transitions and multi-arm forests estimated 108 joint contrasts. Most showed
  no effect, with confidence intervals tight enough to exclude absolute risk
  differences above roughly one percentage point.
- **Reverse causation identified by design.** The one transition surviving
  multiple-testing correction (reducing alcohol, +0.98 pp) disappeared when
  only CVD events after the second lifestyle measurement were counted and
  recurred in a baseline-level reference model, consistent with illness
  prompting the change rather than the change affecting risk.
- **External validation.** Re-estimating the transitions in NHEFS (US cohort)
  gave directional agreement for 8 of 11 testable contrasts, with overlapping
  intervals throughout.

## How it works

![Project structure](docs/figures/project_structure.svg)

A shared preparation stage builds the analysis cohorts, then two independent
branches produce model artefacts that the web app loads.

**Data preparation** (`src/data_preparation/`, run in order). Scores six
lifestyle domains (sleep, smoking, alcohol, diet, physical activity, mental
health) at baseline and imaging visits from raw touchscreen fields, merges
covariates and outcomes, classifies each participant's CVD timing relative to
the two visits, and recodes exposures to 0/1/2 levels with explicit
transition variables. Lifestyle categories are never imputed.

**Predictive model** (`src/predictive_model/`). scikit-learn pipelines with
spline terms for age and BMI, one-hot categories and standardised scores;
randomised hyperparameter search; evaluation by AUC with DeLong tests,
calibration slope and intercept, scaled Brier score, decision curves and
bootstrap optimism correction. The final pipeline is exported with joblib.

**Causal forests** (`src/causal_forest/`, R with `grf`). Each transition is a
binary treatment (moved vs stayed) under target-trial rules, with honest
forests, doubly robust average treatment effects, overlap weighting and
trimming, heterogeneity tests (BLP, RATE), subgroup effects and seed-stability
checks. `combined_variable.R` uses multi-arm forests for paired changes.
`unmeasured_confounding.R` adds E-values and quantitative bias analysis.
`baseline_reference.R` estimates baseline-level effects as a temporally
unambiguous reference. Results are exported to joblib via reticulate.

**Web app** (`cvd_webapp/`). A Streamlit questionnaire reproduces the UK
Biobank items, computes absolute risk from the predictive pipeline, looks up
the causal effect of each available lifestyle improvement, and shows every
number with its confidence interval. Predictive and causal results are kept
separate. See [how the app computes risk](docs/figures/risk_calculation.svg)
and [what is new in the design](docs/figures/innovation.svg).

## Repository layout

```
src/data_preparation/    01_overview → 02_data_merge → 03_eda → 04_cvd_timing → 05_cohort_recoding
src/predictive_model/    predictive_model.py
src/causal_forest/       single_variable.R · combined_variable.R · unmeasured_confounding.R · baseline_reference.R
cvd_webapp/              Streamlit app, model loading, questionnaire scoring, recommendations, Dockerfile
docs/                    PROJECT_DETAILS.md (full methods) · WORKFLOW.md (commands) · figures/
```

## Getting started

```bash
pip install -r requirements.txt
Rscript src/causal_forest/install_packages.R

# run the app with the committed model artefacts
cd cvd_webapp && pip install -r requirements.txt && streamlit run app.py
```

Pipeline scripts read paths from environment variables (defaults point to the
HPC layout). Step-by-step commands are in [`docs/WORKFLOW.md`](docs/WORKFLOW.md);
design notes, the pipeline flowchart and the file-rename history are in
[`docs/PROJECT_DETAILS.md`](docs/PROJECT_DETAILS.md).

UK Biobank data cannot be redistributed. No participant-level data, derived
cohorts or result tables are stored here; only the three small model artefacts
the app needs are committed.

**Stack:** Python (pandas, scikit-learn, XGBoost, LightGBM, joblib, Plotly, Streamlit), R (grf, tidyverse, reticulate), Docker.
