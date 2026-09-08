# CVD Lifestyle Risk Explorer

Streamlit front end for the dissertation's predictive (Stream 1) and causal
(Stream 2) CVD analyses. Three parts: an overview, a questionnaire covering
the nine modelling variables, and a results page with predicted risk, causal
effect estimates, recommendations, and references.

## Run

The predictive model was serialized with scikit-learn 1.0.2, so the webpage
must run with Python 3.9 and the package versions in `requirements.txt`. The
existing `.venv` uses Python 3.14 and must not be used for this model.

```bash
cd cvd_webapp
conda create -n cvd-webapp python=3.9.18 -y
conda activate cvd-webapp
pip install -r requirements.txt
streamlit run app.py
```

The app starts without model files and shows a "model not loaded" notice in
each results panel until they are supplied.

## Public deployment

The included `Dockerfile` keeps Python 3.9 and scikit-learn 1.0.2 aligned with
the predictive model. To publish with Render:

1. Put this directory in a private GitHub repository. Keep all three files in
   `model/` tracked; they are required at runtime.
2. In Render, select **New > Blueprint** and connect the repository.
3. Render reads `render.yaml`, builds the container, and gives the service a
   public `onrender.com` URL.

The free Render plan may sleep after inactivity, so the first request after a
quiet period can take longer. Use a paid always-on instance for a dissertation
demonstration where immediate startup matters.

Do not commit `.venv`; `.dockerignore` excludes it from the image. The app does
not require API keys or other deployment secrets.

To test the deployment image locally when Docker is installed:

```bash
docker build -t cvd-webapp .
docker run --rm -p 8501:8501 cvd-webapp
```

Then open `http://localhost:8501`.

## Model files you supply (place in `cvd_webapp/model/`)

### 1. `predictive_model.joblib` or `logistic.joblib`

A fitted sklearn `Pipeline` (preprocessing included) whose `predict_proba`
accepts a DataFrame with exactly these 9 raw columns:

```
age, bmi, alcohol_score_0.0, diet_score_0.0, mental_score_0.0,
sleep_category_0.0, smoking_category_0.0, physical_category_0.0, sex
```

Raw codings: sleep/smoking/physical 0-2, alcohol 0-4, diet 0-3, mental 0-3,
sex 0 = female / 1 = male.

The loader first looks for `predictive_model.joblib`, then falls back to the
existing `model/logistic.joblib`. Both are loaded with the matching
scikit-learn 1.0.2 runtime pinned in `requirements.txt`.

### 2. `single_variable_ate.joblib`

A pandas DataFrame (or a dict containing one) with one row per lifestyle
transition. Layout of `single_variable/primary_ate_table1.csv` works as-is:

- `domain` — one of `sleep, diet, pa, mental, alcohol, smoking`
- `transition` — string like `0->1` (or separate `baseline`/`destination` int columns)
- `ate`/`estimate`, `ate_low`/`low`, `ate_high`/`high` — risk difference as a
  fraction (the app multiplies by 100 to show percentage points)
- optional `sig` — True/False; if absent, significance is derived from the CI

Transitions use the collapsed causal 0/1/2 coding.

### 3. `combined_variable_ate.joblib`

A pandas DataFrame in the layout of
`causalml/two_variable_results_primary_1/multiarm_ate.csv`:

- `pair` (e.g. `sleep_diet`), `model_type` (`improvement`/`deterioration`/`mixed`),
  `arm`, `arm_label` (e.g. `sleep improves only`), `estimate`, `low`, `high`,
  optional `p`, `p_fdr`

The app currently displays the `improvement` scenario for pairs where the
user has room to improve both factors.

Creating the two table joblibs from the existing CSVs is one line each, e.g.:

```python
import joblib, pandas as pd
joblib.dump(pd.read_csv("single_variable/primary_ate_table1.csv"), "cvd_webapp/model/single_variable_ate.joblib")
joblib.dump(pd.read_csv("causalml/two_variable_results_primary_1/multiarm_ate.csv"), "cvd_webapp/model/combined_variable_ate.joblib")
```

Alternatively the analyses can write the joblib directly, skipping the CSV step:
`causalml/single_variable_results_primary_1/single_variable3_joblib.R` and
`causalml/two_variable_results_primary_1/combined_variable6_joblib.R` are copies
of the CSV-writing scripts whose only output is one `.joblib` each. Set
`JOBLIB_PATH` to the file in `model/` and they land in place. Each payload is a
dict of every analysis table, with the app-facing one first — which is the entry
`model_io` picks up.

After adding or replacing files, use the sidebar's **Reload model files**
button (or restart the app).

## Where things live

- `app.py` — page flow and charts
- `theme.py` — design tokens, stylesheet, and HTML components (hero, stat
  cards, advice cards, chart legend). Edit visual styling here.
- `.streamlit/config.toml` — Streamlit theme (colours, fonts, radii) and
  `toolbarMode = "minimal"` to hide the dev toolbar. Tokens are duplicated in
  `theme.py`; change both together.
- `questions.py` — questionnaire wording and scoring (raw scores + causal
  0/1/2 recodes, reproducing `data_merge_new.ipynb` and
  `cohort/build_table1_causal_cohorts.py`). Edit wording here.
- `model_io.py` — joblib loading + schema validation
- `recommendations.py` — advice text, transition logic, references, disclaimer
- `COHORT_EVENT_RATE` in `app.py` (5.98%) is the Stream 1 primary cohort event
  rate used to contextualize the user's predicted risk; change it if the
  supplied predictive model was trained on a different cohort/outcome.

## Front-end notes

The page is patient-facing, so each section leads with one short line and the
methodological detail (cohort, C-statistic, pp definition, null-effect caveats,
references, disclaimer) sits in `st.expander` drawers rather than on the
surface. The three stages — overview, questionnaire, results — render one at a
time instead of stacking.

The questionnaire is a bordered `st.container`, **not** an `st.form`: widgets
inside a form do not rerun the script when they change, so a follow-up question
could not appear in response to the answer above it (the alcohol frequency
question is only shown to current drinkers). The card styling therefore hangs
off a `:has()` rule in `theme.py` anchored to the `.qsec.first` marker that
`theme.form_section(first=True)` emits. If you add another bordered container to
this page, check that rule still selects only the questionnaire.

Forest-plot colour encodes **polarity**, not category: blue below zero (risk
down), red above (risk up), gray for intervals crossing zero. Marker shape is a
second, non-colour channel for the same distinction (filled = statistically
clear, hollow = no clear effect), so the chart stays readable for colour-vision
deficiency and in print. Blue/red is used rather than green/red both because
green/red fails CVD separation and because green is the UI's chrome colour.
