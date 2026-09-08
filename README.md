# Causal and Predictive ML for Cardiovascular Risk in UK Biobank

![Python](https://img.shields.io/badge/Python-3.10+-3776AB?logo=python&logoColor=white)
![R](https://img.shields.io/badge/R-grf%20%7C%20tidyverse-276DC3?logo=r&logoColor=white)
![scikit-learn](https://img.shields.io/badge/scikit--learn-pipelines-F7931E?logo=scikit-learn&logoColor=white)
![Streamlit](https://img.shields.io/badge/Streamlit-app-FF4B4B?logo=streamlit&logoColor=white)
![Docker](https://img.shields.io/badge/Docker-Render%20blueprint-2496ED?logo=docker&logoColor=white)

**Can machine learning say both *who* is at risk of heart disease and *which lifestyle change* would actually lower it?**
MSc Health Data Science dissertation (UCL, 2026) answering both questions in one 500,000-participant cohort,
end to end: data engineering on raw UK Biobank tables, predictive and causal modelling, external validation,
and a deployed decision-support app.

## Results at a glance

| | |
| ---: | --- |
| **458,840** | participants in the baseline CVD prediction model, four model families compared |
| **AUC 0.718** (95% CI 0.712 to 0.724) | from self-reported lifestyle and demographics only, no blood tests or clinical history |
| **71,428** | participants with lifestyle measured at two visits, used for causal analysis |
| **24 + 108** | single and joint lifestyle-change effects estimated with doubly robust causal forests |
| **8 of 11** | causal contrasts replicated in direction in an independent US cohort (NHEFS) |
| **1 app** | Streamlit tool that ships both models to end users with full uncertainty |

**Headline finding.** Most short-horizon lifestyle changes had no detectable effect, and the confidence
intervals were tight enough to rule out changes above about one percentage point of absolute risk. The single
result that survived multiple-testing correction, higher risk after cutting alcohol, was traced to reverse
causation by a timing-aware design that fitted three temporal models of the same cohort. Not finding an
effect, and proving why an apparent one was spurious, is the substance of the work.

## What this project demonstrates

| Capability | Evidence in this repository |
| --- | --- |
| **Causal inference at scale** | Target-trial emulation, honest doubly robust causal forests, multi-arm forests for joint treatments, overlap weighting and trimming, E-values and quantitative bias analysis. [`src/causal_forest/`](src/causal_forest/) |
| **Predictive modelling done properly** | Logistic regression vs random forest, XGBoost and LightGBM with tuned pipelines, calibration slope and intercept, scaled Brier, DeLong tests, decision curves, bootstrap optimism correction, subgroup performance. [`predictive_model.py`](src/predictive_model/predictive_model.py) |
| **Data engineering on biobank data** | Reproducible five-step pipeline turning raw touchscreen fields into analysis cohorts, configured by environment variables and run on an HPC cluster. [`src/data_preparation/`](src/data_preparation/) |
| **Epidemiological study design** | Primary, strict post-exposure and baseline-reference temporal models; two outcome definitions; sensitivity analyses for timing, overlap and unmeasured confounding; external validation. [`docs/figures/innovation.svg`](docs/figures/innovation.svg) |
| **Shipping models to users** | scikit-learn pipeline and R results exported as versioned joblib artefacts (R to Python via reticulate), loaded by a Streamlit app with Docker and Render deployment. [`cvd_webapp/`](cvd_webapp/) |
| **Communicating uncertainty** | Every estimate shown with its interval, forest plots for effects, predictive and causal numbers kept separate and explained in plain language. [`docs/figures/risk_calculation.svg`](docs/figures/risk_calculation.svg) |

## How it fits together

![Project structure](docs/figures/project_structure.svg)

One preparation stage feeds two independent modelling branches. Branch A fits the baseline risk model on
every participant. Branch B fits causal forests on the longitudinal cohort. Each branch exports an artefact,
and the web app consumes both. Further diagrams: [how the app computes risk](docs/figures/risk_calculation.svg)
and [what is new in the design](docs/figures/innovation.svg).

## Tech stack

**Python:** pandas, NumPy, scikit-learn, XGBoost, LightGBM, joblib, matplotlib, Plotly, Streamlit.
**R:** grf (generalized random forests), tidyverse, reticulate.
**Infrastructure:** UCL HPC cluster, Docker, Render, Git.

## Repository map

```
src/data_preparation/    01–05: overview → merge & score → EDA → CVD timing → cohort recoding
src/predictive_model/    baseline CVD risk model, evaluation, artefact export
src/causal_forest/       single-variable, combined-variable, confounding sensitivity, baseline reference
cvd_webapp/              Streamlit app, model loading, questionnaire scoring, recommendations
docs/                    PROJECT_DETAILS.md (full methods), WORKFLOW.md (runbook), figures/
```

## Run the app locally

```bash
cd cvd_webapp
pip install -r requirements.txt
streamlit run app.py
```

Full pipeline commands, environment variables and design notes are in
[`docs/PROJECT_DETAILS.md`](docs/PROJECT_DETAILS.md) and [`docs/WORKFLOW.md`](docs/WORKFLOW.md).
UK Biobank data cannot be redistributed, so no participant data or result tables are stored here.

---

Jennifer Sun · MSc Health Data Science, University College London · [github.com/jennifersun09b](https://github.com/jennifersun09b)
