#!/usr/bin/env python3
"""
Predictive modelling (improved version) - HPC batch script.

Non-interactive translation of model_building_improved.ipynb, meant to be run as
a submitted job (e.g. SGE qsub) rather than in a notebook:
  * uses the matplotlib "Agg" backend (no display needed)
  * saves every figure to disk instead of calling plt.show()
  * writes all tables/plots/models into RESULT_DIR
  * prints progress to stdout (captured in the job's .o log)

Environment variables (all optional):
  RUN_HEAVY = 1   -> bigger tuning search + more bootstrap repeats (final run)
  N_JOBS          -> parallel cores; defaults to $NSLOTS (SGE) or -1 (all cores)
  DATA_PATH       -> override the input CSV
  RESULT_DIR      -> override the output folder
  OUTCOME         -> outcome column; defaults to def_CVD_AFTER
  PREDICTIVE_JOBLIB_PATH -> exact app artifact destination

Primary analyses: CVD outcome
Sensitive analyses: CVD/HF/AF outcome
"""

import os
import tempfile
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")            # no interactive display on the compute nodes
import matplotlib.pyplot as plt

import joblib

from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer, SimpleImputer
from sklearn.preprocessing import StandardScaler, OneHotEncoder, SplineTransformer
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.base import clone
from sklearn.model_selection import (train_test_split, GridSearchCV,
                                     RandomizedSearchCV, StratifiedKFold)
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (roc_auc_score, roc_curve, precision_recall_curve,
                             average_precision_score, brier_score_loss)
from sklearn.calibration import calibration_curve
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------
RANDOM_STATE = 42

RUN_HEAVY = os.environ.get("RUN_HEAVY", "0") == "1"
N_JOBS    = int(os.environ.get("N_JOBS", os.environ.get("NSLOTS", "-1")))

DATA_PATH  = os.environ.get(
    "DATA_PATH", "/home/rmhiund/causal_analysis/data_merge/baseline.csv")
RESULT_DIR = os.environ.get(
    "RESULT_DIR", "/home/rmhiund/causal_analysis/predictive_model/results_improved_primary")
MODEL_DIR  = os.path.join(RESULT_DIR, "models")

OUTCOME = os.environ.get("OUTCOME", "def_CVD_AFTER")
PREDICTIVE_JOBLIB_PATH = os.environ.get(
    "PREDICTIVE_JOBLIB_PATH",
    os.path.join(MODEL_DIR, "predictive_model.joblib"),
)

ESTABLISHED_NUM = ["age", "bmi"]
SCORE_COLS      = ["alcohol_score_0.0", "diet_score_0.0", "mental_score_0.0"]
CATEGORY_COLS   = ["sleep_category_0.0", "smoking_category_0.0", "physical_category_0.0"]
BINARY_COLS     = ["sex"]
PREDICTORS      = ESTABLISHED_NUM + SCORE_COLS + CATEGORY_COLS + BINARY_COLS

N_BOOT = 1000 if RUN_HEAVY else 300
N_ITER = 25 if RUN_HEAVY else 8


def log(*args):
    """Print immediately so progress shows up in the job log while running."""
    print(*args, flush=True)


def export_app_model(model):
    """Atomically write and smoke-test the pipeline used by cvd_webapp."""
    destination = os.path.abspath(PREDICTIVE_JOBLIB_PATH)
    os.makedirs(os.path.dirname(destination), exist_ok=True)
    fd, temporary = tempfile.mkstemp(
        prefix=os.path.basename(destination) + ".tmp-",
        dir=os.path.dirname(destination),
    )
    os.close(fd)
    try:
        joblib.dump(model, temporary, compress=3)
        restored = joblib.load(temporary)
        probe = pd.DataFrame([{column: 1.0 for column in PREDICTORS}])
        probe.loc[0, ["age", "bmi"]] = [55.0, 26.0]
        probabilities = np.asarray(restored.predict_proba(probe))
        if probabilities.shape != (1, 2) or not np.isfinite(probabilities).all():
            raise RuntimeError(
                f"reloaded predictive model returned invalid shape {probabilities.shape}"
            )
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    log("wrote and verified app model:", destination)


# ----------------------------------------------------------------------------
# Preprocessing
# ----------------------------------------------------------------------------
def build_preprocessors():
    # Logistic model: splines for age/bmi, scaling for scores, one-hot for
    # categories, most-frequent fill for sex. Everything is fitted inside the
    # pipeline (train data only) so nothing leaks from the test set.
    preprocessor = ColumnTransformer(transformers=[
        ("age_bmi", Pipeline([
            ("impute", IterativeImputer(random_state=RANDOM_STATE, max_iter=10)),
            ("spline", SplineTransformer(degree=3, n_knots=4, include_bias=False))]),
         ESTABLISHED_NUM),
        ("scores", Pipeline([
            ("impute", IterativeImputer(random_state=RANDOM_STATE, max_iter=10)),
            ("scale", StandardScaler())]), SCORE_COLS),
        ("cats", Pipeline([
            ("impute", SimpleImputer(strategy="most_frequent")),
            ("ohe", OneHotEncoder(drop="first"))]), CATEGORY_COLS),
        ("sex", SimpleImputer(strategy="most_frequent"), BINARY_COLS),
    ])

    # Tree models: impute only, no scaling.
    preprocessor_tree = ColumnTransformer(transformers=[
        ("num", IterativeImputer(random_state=RANDOM_STATE, max_iter=10),
         ESTABLISHED_NUM + SCORE_COLS),
        ("cats", SimpleImputer(strategy="most_frequent"), CATEGORY_COLS),
        ("sex", SimpleImputer(strategy="most_frequent"), BINARY_COLS),
    ])
    return preprocessor, preprocessor_tree


def build_models(preprocessor, preprocessor_tree):
    logistic = Pipeline([
        ("preprocessor", preprocessor),
        ("classifier", LogisticRegression(max_iter=4000, random_state=RANDOM_STATE))])

    random_forest = Pipeline([
        ("preprocessor", preprocessor_tree),
        ("classifier", RandomForestClassifier(n_jobs=N_JOBS, random_state=RANDOM_STATE))])

    xgboost = Pipeline([
        ("preprocessor", preprocessor_tree),
        ("classifier", XGBClassifier(eval_metric="logloss", n_jobs=N_JOBS,
                                     random_state=RANDOM_STATE))])

    lightgbm = Pipeline([
        ("preprocessor", preprocessor_tree),
        ("classifier", LGBMClassifier(n_jobs=N_JOBS, random_state=RANDOM_STATE,
                                      verbose=-1))])

    return {"logistic": logistic, "random_forest": random_forest,
            "xgboost": xgboost, "lightgbm": lightgbm}


PARAM_GRIDS = {
    "logistic": {
        "classifier__penalty": ["l2"],
        "classifier__C": [0.01, 0.1, 1, 10, 100],
        "classifier__solver": ["lbfgs"]},
    "random_forest": {
        "classifier__n_estimators": [300, 500, 800],
        "classifier__min_samples_leaf": [20, 50, 100],
        "classifier__max_features": ["sqrt", 0.5]},
    "xgboost": {
        "classifier__n_estimators": [400, 600, 800],
        "classifier__learning_rate": [0.03, 0.05, 0.1],
        "classifier__max_depth": [3, 4, 6],
        "classifier__subsample": [0.8, 1.0],
        "classifier__colsample_bytree": [0.8, 1.0]},
    "lightgbm": {
        "classifier__n_estimators": [400, 600, 800],
        "classifier__learning_rate": [0.03, 0.05, 0.1],
        "classifier__num_leaves": [31, 63, 127],
        "classifier__subsample": [0.8, 1.0],
        "classifier__colsample_bytree": [0.8, 1.0]},
}


def tune_models(models, X_train, y_train):
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
    best_models = {}
    for name, model in models.items():
        log("tuning:", name)
        if name == "logistic":
            search = GridSearchCV(model, PARAM_GRIDS[name], cv=cv,
                                  scoring="roc_auc", n_jobs=N_JOBS, refit=True)
        else:
            search = RandomizedSearchCV(model, PARAM_GRIDS[name], n_iter=N_ITER, cv=3,
                                        scoring="roc_auc", n_jobs=N_JOBS,
                                        random_state=RANDOM_STATE, refit=True)
        search.fit(X_train, y_train)
        best_models[name] = search.best_estimator_
        log("  best cv roc_auc:", round(search.best_score_, 4),
            " params:", search.best_params_)
    return best_models


# ----------------------------------------------------------------------------
# Metric helpers
# ----------------------------------------------------------------------------
def logit(p, eps=1e-6):
    p = np.clip(p, eps, 1 - eps)
    return np.log(p / (1 - p))


def scaled_brier(y, p):
    brier = brier_score_loss(y, p)
    prev = np.mean(y)
    return brier, 1 - brier / (prev * (1 - prev))


def calibration_slope_intercept(y, p):
    lp = logit(p).reshape(-1, 1)
    try:
        m = LogisticRegression(penalty=None, solver="lbfgs", max_iter=1000).fit(lp, y)
    except (ValueError, TypeError):
        m = LogisticRegression(penalty="none", solver="lbfgs", max_iter=1000).fit(lp, y)
    slope = float(m.coef_[0][0])
    try:
        import statsmodels.api as sm
        glm = sm.GLM(np.asarray(y), np.ones((len(y), 1)),
                     family=sm.families.Binomial(), offset=logit(p)).fit()
        intercept = float(glm.params[0])
    except Exception:
        intercept = float(logit(np.mean(y)) - logit(np.mean(p)))
    return slope, intercept


def bootstrap_ci(y, p, metric_fn, n_boot=1000, seed=RANDOM_STATE):
    y = np.asarray(y); p = np.asarray(p)
    rng = np.random.default_rng(seed)
    stats = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) == 2:
            stats.append(metric_fn(y[idx], p[idx]))
    lo, hi = np.percentile(stats, [2.5, 97.5])
    return metric_fn(y, p), lo, hi


# ----------------------------------------------------------------------------
# DeLong test
# ----------------------------------------------------------------------------
def _midrank(x):
    order = np.argsort(x)
    z = x[order]
    n = len(x)
    t = np.zeros(n)
    i = 0
    while i < n:
        j = i
        while j < n and z[j] == z[i]:
            j += 1
        t[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    out = np.empty(n)
    out[order] = t
    return out


def delong_test(y_true, prob_a, prob_b):
    from scipy import stats
    y_true = np.asarray(y_true)
    order = (-y_true).argsort(kind="mergesort")   # positives first
    m = int(y_true.sum())
    preds = np.vstack((prob_a, prob_b))[:, order]
    n = preds.shape[1] - m
    pos, neg = preds[:, :m], preds[:, m:]
    tx = np.vstack([_midrank(pos[r]) for r in range(2)])
    ty = np.vstack([_midrank(neg[r]) for r in range(2)])
    tz = np.vstack([_midrank(preds[r]) for r in range(2)])
    aucs = tz[:, :m].sum(axis=1) / m / n - (m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx) / n
    v10 = 1.0 - (tz[:, m:] - ty) / m
    cov = np.cov(v01) / m + np.cov(v10) / n
    l = np.array([[1.0, -1.0]])
    var = float(l @ cov @ l.T)
    z = (aucs[0] - aucs[1]) / np.sqrt(var + 1e-12)
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    return float(aucs[0]), float(aucs[1]), float(p)


# ----------------------------------------------------------------------------
# Decision curve analysis
# ----------------------------------------------------------------------------
def net_benefit(y_true, y_prob, thresholds):
    y_true = np.asarray(y_true)
    n = len(y_true)
    prev = y_true.mean()
    nb_model, nb_all = [], []
    for t in thresholds:
        pred_pos = y_prob >= t
        tp = np.sum(pred_pos & (y_true == 1))
        fp = np.sum(pred_pos & (y_true == 0))
        w = t / (1 - t)
        nb_model.append(tp / n - (fp / n) * w)
        nb_all.append(prev - (1 - prev) * w)
    return np.array(nb_model), np.array(nb_all)


# ----------------------------------------------------------------------------
# Internal validation (optional, heavy)
# ----------------------------------------------------------------------------
def bootstrap_optimism(estimator, X, y, n_boot=200, seed=RANDOM_STATE):
    rng = np.random.default_rng(seed)
    Xv = X.reset_index(drop=True); yv = y.reset_index(drop=True)
    n = len(Xv)
    base = clone(estimator).fit(Xv, yv)
    apparent = roc_auc_score(yv, base.predict_proba(Xv)[:, 1])
    opt = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        Xb, yb = Xv.iloc[idx], yv.iloc[idx]
        if yb.nunique() < 2:
            continue
        m = clone(estimator).fit(Xb, yb)
        auc_boot = roc_auc_score(yb, m.predict_proba(Xb)[:, 1])
        auc_orig = roc_auc_score(yv, m.predict_proba(Xv)[:, 1])
        opt.append(auc_boot - auc_orig)
    optimism = float(np.mean(opt))
    return apparent, optimism, apparent - optimism


# ----------------------------------------------------------------------------
# Main pipeline
# ----------------------------------------------------------------------------
def main():
    os.makedirs(RESULT_DIR, exist_ok=True)
    os.makedirs(MODEL_DIR, exist_ok=True)
    log(f"RUN_HEAVY={RUN_HEAVY}  N_JOBS={N_JOBS}  N_BOOT={N_BOOT}  N_ITER={N_ITER}")
    log("DATA_PATH :", DATA_PATH)
    log("RESULT_DIR:", RESULT_DIR)
    log("OUTCOME   :", OUTCOME)
    log("APP MODEL :", PREDICTIVE_JOBLIB_PATH)

    # --- load ---
    df = pd.read_csv(DATA_PATH, sep=",")
    log("loaded:", df.shape)

    rename_map = {
        "genetic_sex": "sex",
        "age_defined_baseline": "age",
        "merged_bmi": "bmi",
    }
    df = df.rename(columns={
        old: new for old, new in rename_map.items()
        if old in df.columns and new not in df.columns
    })
    required_columns = set(PREDICTORS + [OUTCOME])
    missing_columns = sorted(required_columns - set(df.columns))
    if missing_columns:
        raise ValueError(f"input dataset is missing required columns: {missing_columns}")

    y_all = pd.to_numeric(df[OUTCOME], errors="raise").astype(int)
    if not set(y_all.unique()).issubset({0, 1}):
        raise ValueError(f"{OUTCOME} must be coded 0/1")
    log("N=%d  events=%d  prevalence=%.2f%%"
        % (len(df), int(y_all.sum()), y_all.mean() * 100))
    log("missing %:\n", ((df.isnull().sum() / len(df)) * 100).round(3).to_string())

    # --- split ---
    X = df[PREDICTORS].copy()
    Y = df[OUTCOME].astype(int)
    X_train, X_test, y_train, y_test = train_test_split(
        X, Y, test_size=0.30, stratify=Y, random_state=RANDOM_STATE)
    log("train:", X_train.shape, " test:", X_test.shape)
    log("train prevalence: %.2f%%   test prevalence: %.2f%%"
        % (y_train.mean() * 100, y_test.mean() * 100))

    # --- build + tune ---
    preprocessor, preprocessor_tree = build_preprocessors()
    models = build_models(preprocessor, preprocessor_tree)
    best_models = tune_models(models, X_train, y_train)

    # --- evaluate ---
    metrics, probs = [], {}
    for name, model in best_models.items():
        y_prob = model.predict_proba(X_test)[:, 1]
        probs[name] = y_prob
        auc, auc_lo, auc_hi = bootstrap_ci(y_test, y_prob, roc_auc_score, n_boot=N_BOOT)
        brier, brier_lo, brier_hi = bootstrap_ci(y_test, y_prob, brier_score_loss, n_boot=N_BOOT)
        _, sbrier = scaled_brier(y_test, y_prob)
        slope, intercept = calibration_slope_intercept(y_test, y_prob)
        metrics.append({
            "model": name,
            "c_statistic_auc": auc, "auc_ci_low": auc_lo, "auc_ci_high": auc_hi,
            "average_precision": average_precision_score(y_test, y_prob),
            "brier": brier, "brier_ci_low": brier_lo, "brier_ci_high": brier_hi,
            "scaled_brier": sbrier,
            "calibration_slope": slope, "calibration_intercept": intercept})

    metrics_df = pd.DataFrame(metrics).set_index("model")
    metrics_df.to_csv(os.path.join(RESULT_DIR, "model_metrics.csv"))
    log("metrics:\n", metrics_df.round(4).to_string())

    # --- ROC ---
    plt.figure(figsize=(6, 6))
    for name, y_prob in probs.items():
        fpr, tpr, _ = roc_curve(y_test, y_prob)
        plt.plot(fpr, tpr, label=f"{name} (AUC={roc_auc_score(y_test, y_prob):.3f})")
    plt.plot([0, 1], [0, 1], "k--", label="No skill")
    plt.xlabel("False Positive Rate"); plt.ylabel("True Positive Rate")
    plt.title("ROC Curve"); plt.legend(fontsize=8)
    plt.savefig(os.path.join(RESULT_DIR, "roc_curve.png"), dpi=150, bbox_inches="tight")
    plt.close()

    # --- PR ---
    plt.figure(figsize=(6, 6))
    for name, y_prob in probs.items():
        precision, recall, _ = precision_recall_curve(y_test, y_prob)
        plt.plot(recall, precision,
                 label=f"{name} (AP={average_precision_score(y_test, y_prob):.3f})")
    plt.axhline(y_test.mean(), ls="--", c="k", label=f"Prevalence={y_test.mean():.3f}")
    plt.xlabel("Recall"); plt.ylabel("Precision")
    plt.title("Precision-Recall Curve"); plt.legend(fontsize=8)
    plt.savefig(os.path.join(RESULT_DIR, "pr_curve.png"), dpi=150, bbox_inches="tight")
    plt.close()

    # --- calibration ---
    plt.figure(figsize=(6, 6))
    calib_rows = []                       # collect the binned points for each model
    for name, y_prob in probs.items():
        frac_pos, mean_pred = calibration_curve(y_test, y_prob, n_bins=10, strategy="quantile")
        plt.plot(mean_pred, frac_pos, marker="o", label=name)
        for b, (mp, fp) in enumerate(zip(mean_pred, frac_pos), start=1):
            calib_rows.append({"model": name, "bin": b,
                               "mean_predicted_risk": mp,
                               "observed_event_frequency": fp})
    plt.plot([0, 1], [0, 1], "k--", label="Perfect calibration")
    plt.xlabel("Mean predicted risk"); plt.ylabel("Observed event frequency")
    plt.title("Calibration Plot"); plt.legend(fontsize=8)
    plt.savefig(os.path.join(RESULT_DIR, "calibration_plot.png"), dpi=150, bbox_inches="tight")
    plt.close()

    # binned calibration points (reliability curve) for every model -> CSV
    calib_df = pd.DataFrame(calib_rows)
    calib_df.to_csv(os.path.join(RESULT_DIR, "calibration_curve.csv"), index=False)
    log("saved calibration_curve.csv (%d rows)" % len(calib_df))

    # raw per-patient test-set predictions: lets you rebuild ANY calibration
    # curve later (binned / flexible / primary-vs-sensitive) without the joblib
    pred_df = pd.DataFrame({"y_true": y_test.to_numpy()})
    for name, y_prob in probs.items():
        pred_df[name + "_prob"] = y_prob
    pred_df.to_csv(os.path.join(RESULT_DIR, "test_predictions.csv"), index=False)
    log("saved test_predictions.csv (%d rows)" % len(pred_df))

    # --- DeLong ---
    ref = "logistic"
    rows = []
    for name in best_models:
        if name == ref:
            continue
        a, b, p = delong_test(y_test, probs[ref], probs[name])
        rows.append({"comparison": f"{ref} vs {name}", "auc_logistic": a,
                     "auc_other": b, "difference": a - b, "p_value": p})
    delong_df = pd.DataFrame(rows)
    delong_df.to_csv(os.path.join(RESULT_DIR, "delong_tests.csv"), index=False)
    log("delong:\n", delong_df.round(4).to_string())

    # --- decision curve ---
    thresholds = np.linspace(0.01, 0.5, 60)
    plt.figure(figsize=(7, 6))
    nb_all = None
    for name, y_prob in probs.items():
        nb_model, nb_all = net_benefit(y_test, y_prob, thresholds)
        plt.plot(thresholds, nb_model, label=name)
    plt.plot(thresholds, nb_all, "k--", label="Treat all")
    plt.axhline(0, color="gray", lw=1, label="Treat none")
    plt.ylim(-0.02, None)
    plt.xlabel("Threshold probability"); plt.ylabel("Net benefit")
    plt.title("Decision Curve Analysis"); plt.legend(fontsize=8)
    plt.savefig(os.path.join(RESULT_DIR, "decision_curve.png"), dpi=150, bbox_inches="tight")
    plt.close()

    # --- interpretation: logistic odds ratios ---
    log_model = best_models["logistic"]
    # sklearn 1.0.2 imputers lack get_feature_names_out, but they don't change
    # columns. Build names in ColumnTransformer order: splines, scores, one-hot, sex.
    pre = log_model.named_steps["preprocessor"]
    spline = pre.named_transformers_["age_bmi"].named_steps["spline"]
    ohe    = pre.named_transformers_["cats"].named_steps["ohe"]
    feat_names = (list(spline.get_feature_names_out(ESTABLISHED_NUM))
                  + SCORE_COLS
                  + list(ohe.get_feature_names_out(CATEGORY_COLS))
                  + BINARY_COLS)
    coefs = log_model.named_steps["classifier"].coef_[0]
    or_df = pd.DataFrame({"feature": feat_names, "coef": coefs, "odds_ratio": np.exp(coefs)})
    or_df = or_df.reindex(np.argsort(-np.abs(coefs))).reset_index(drop=True)
    or_df.to_csv(os.path.join(RESULT_DIR, "logistic_odds_ratios.csv"), index=False)

    # --- interpretation: tree feature importances ---
    importances = {}
    for name in ["random_forest", "xgboost", "lightgbm"]:
        mdl = best_models[name]
        # tree preprocessor only imputes (no one-hot), so names are the columns in order
        fn = ESTABLISHED_NUM + SCORE_COLS + CATEGORY_COLS + BINARY_COLS
        importances[name] = pd.Series(mdl.named_steps["classifier"].feature_importances_, index=fn)
    imp_df = pd.DataFrame(importances).fillna(0.0)
    imp_df.to_csv(os.path.join(RESULT_DIR, "tree_feature_importance.csv"))

    # --- subgroup performance ---
    y_arr = y_test.reset_index(drop=True)
    sex_arr = X_test["sex"].reset_index(drop=True)
    age_arr = X_test["age"].reset_index(drop=True)

    def subgroup_row(mask, label):
        mask = np.asarray(mask)
        if len(np.unique(y_arr[mask])) < 2:
            return None
        row = {"subgroup": label, "n": int(mask.sum()), "events": int(y_arr[mask].sum())}
        for name in probs:
            row[name + "_auc"] = roc_auc_score(y_arr[mask], probs[name][mask])
        return row

    sg_rows = [subgroup_row(sex_arr == 0, "sex=0"),
               subgroup_row(sex_arr == 1, "sex=1"),
               subgroup_row(age_arr < 60, "age<60"),
               subgroup_row(age_arr >= 60, "age>=60")]
    subgroup_df = pd.DataFrame([r for r in sg_rows if r is not None])
    subgroup_df.to_csv(os.path.join(RESULT_DIR, "subgroup_auc.csv"), index=False)
    log("subgroups:\n", subgroup_df.round(4).to_string())

    # --- save models ---
    for name, model in best_models.items():
        joblib.dump(model, os.path.join(MODEL_DIR, f"{name}.joblib"))
    log("saved models to:", MODEL_DIR)
    export_app_model(best_models["logistic"])

    # --- optional internal validation ---
    if RUN_HEAVY:
        app, opt, corrected = bootstrap_optimism(best_models["logistic"],
                                                 X_train, y_train, n_boot=200)
        log("logistic apparent AUC=%.4f  optimism=%.4f  corrected AUC=%.4f"
            % (app, opt, corrected))

    log("DONE.")


if __name__ == "__main__":
    main()
