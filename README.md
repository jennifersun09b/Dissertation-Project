# Causal and Predictive Machine Learning for Cardiovascular Risk in UK Biobank

A predictive model of incident cardiovascular disease (CVD) and a causal
analysis of longitudinal lifestyle change, estimated in the same UK Biobank
cohort and delivered through an interactive risk tool. MSc Health Data Science
dissertation, University College London, 2026.

Established risk calculators identify who is at risk but say little about
whether changing a lifestyle factor would lower that risk. This project
addresses both questions with complementary methods: a calibrated prediction
model on 458,840 participants, and doubly robust causal forests on 71,428
participants with lifestyle measured at two visits. A timing-aware design
fits three temporal models of the same cohort so that reverse causation is
detected by construction. Estimates were externally validated in an independent
US cohort, and both result streams feed a Streamlit application.

UK Biobank data cannot be redistributed; the repository contains code, documentation and the three model artefacts the app requires.

![Project structure](docs/figures/project_structure.svg)

## Results

- **Prediction.** Logistic regression, random forest, XGBoost and LightGBM
  reached AUC 0.714 to 0.718 with good calibration from age, sex, BMI and six
  lifestyle scores. Logistic regression was retained for interpretability
  (AUC 0.718, 95% CI 0.712 to 0.724).
- **Causal effects.** 24 single-domain lifestyle transitions and 108 joint
  contrasts showed little evidence of effect. Confidence intervals excluded
  absolute risk differences above roughly one percentage point.
- **Reverse causation.** The only transition surviving multiple-testing
  correction (reducing alcohol, +0.98 pp) disappeared when only events after
  the second measurement were counted and recurred in a baseline-level
  reference model, consistent with illness prompting the change.
- **External validation.** Directional agreement for 8 of 11 testable
  contrasts in NHEFS, with overlapping intervals throughout.

## Methods

| Component | Approach |
| --- | --- |
| `src/data_preparation/` | Five ordered scripts: score six lifestyle domains at two visits from raw touchscreen fields, merge covariates and outcomes, classify CVD timing relative to both visits, recode exposures to 0/1/2 levels with explicit transition variables. Lifestyle categories are never imputed. |
| `src/predictive_model/` | scikit-learn pipelines with spline terms, randomised hyperparameter search; AUC with DeLong tests, calibration slope and intercept, scaled Brier score, decision curves, bootstrap optimism correction. |
| `src/causal_forest/` | R `grf`. Each transition is a binary treatment under target-trial rules: honest forests, doubly robust ATEs, overlap weighting and trimming, heterogeneity tests (BLP, RATE), subgroup and seed-stability checks; multi-arm forests for paired changes; E-values and quantitative bias analysis; a baseline-level reference model. |
| `cvd_webapp/` | Streamlit questionnaire reproducing the UK Biobank items; absolute risk from the predictive pipeline, causal effect of each available improvement from the forest tables, every estimate shown with its interval. Predictive and causal results are kept separate. |

Further diagrams: [how the app computes risk](docs/figures/risk_calculation.svg),
[what is new in the design](docs/figures/innovation.svg).

## Getting started

```bash
pip install -r requirements.txt
Rscript src/causal_forest/install_packages.R
cd cvd_webapp && pip install -r requirements.txt && streamlit run app.py
```

Pipeline scripts read paths from environment variables. Commands for every
step are in [`docs/WORKFLOW.md`](docs/WORKFLOW.md); design notes, the pipeline
flowchart and the file history are in [`docs/PROJECT_DETAILS.md`](docs/PROJECT_DETAILS.md).

**Stack:** Python (pandas, scikit-learn, XGBoost, LightGBM, joblib, Plotly, Streamlit), R (grf, tidyverse, reticulate), Docker.
