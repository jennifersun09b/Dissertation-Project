# Loading layer for the three user-supplied joblib files.
#
# Expected files (drop them into cvd_webapp/model/):
#   predictive_model.joblib      fitted sklearn Pipeline with predict_proba,
#                                taking a DataFrame of the 9 raw columns:
#                                age, bmi, alcohol_score_0.0, diet_score_0.0,
#                                mental_score_0.0, sleep_category_0.0,
#                                smoking_category_0.0, physical_category_0.0, sex
#   single_variable_ate.joblib   precomputed single-variable transition table
#                                (DataFrame, or dict containing one) with per-
#                                transition ATE + 95% CI, as exported from the
#                                grf analysis (primary_ate_table1.csv layout)
#   combined_variable_ate.joblib precomputed two-variable multi-arm table
#                                (multiarm_ate.csv layout: pair, model_type,
#                                arm_label, estimate, low, high, p, p_fdr)
#
# Every loader returns (payload, error_message); exactly one is None. The app
# renders a "model not loaded" notice from the message rather than crashing.

from pathlib import Path

import joblib
import pandas as pd

MODEL_DIR = Path(__file__).resolve().parent / "model"

PREDICTIVE_FILE = MODEL_DIR / "predictive_model.joblib"
SINGLE_ATE_FILE = MODEL_DIR / "single_variable_ate.joblib"
COMBINED_ATE_FILE = MODEL_DIR / "combined_variable_ate.joblib"

PREDICTOR_COLUMNS = [
    "age", "bmi", "alcohol_score_0.0", "diet_score_0.0", "mental_score_0.0",
    "sleep_category_0.0", "smoking_category_0.0", "physical_category_0.0", "sex",
]

# Column aliases accepted for the single-variable table
_SINGLE_ALIASES = {
    "estimate": ["estimate", "ate"],
    "low": ["low", "ate_low", "ci_low", "conf_low"],
    "high": ["high", "ate_high", "ci_high", "conf_high"],
    "domain": ["domain"],
    "transition": ["transition", "contrast", "contrast_id"],
}


def _missing(path):
    return (
        f"`{path.name}` not found in `{MODEL_DIR}`. "
        "Add the file and reload the page."
    )


def load_predictive():
    """Load the predictive pipeline and smoke-test predict_proba on one row."""
    if not PREDICTIVE_FILE.exists():
        return None, _missing(PREDICTIVE_FILE)
    try:
        model = joblib.load(PREDICTIVE_FILE)
    except Exception as exc:  # old-sklearn pickles raise many exception types
        return None, (
            f"Could not deserialize `{PREDICTIVE_FILE.name}`: {exc}. "
            "Joblib files pickled under an older scikit-learn may need to be "
            "re-exported with the version pinned in requirements.txt."
        )
    if not hasattr(model, "predict_proba"):
        return None, f"`{PREDICTIVE_FILE.name}` has no predict_proba method."
    try:
        probe = pd.DataFrame([{c: 1.0 for c in PREDICTOR_COLUMNS}])
        probe["age"], probe["bmi"] = 55.0, 26.0
        model.predict_proba(probe)
    except Exception as exc:
        return None, (
            f"`{PREDICTIVE_FILE.name}` loaded but failed a test prediction on "
            f"the expected 9 columns ({', '.join(PREDICTOR_COLUMNS)}): {exc}"
        )
    return model, None


def _as_dataframe(obj, filename):
    """Normalize supported causal-table payloads to a DataFrame."""
    if isinstance(obj, pd.DataFrame):
        return obj
    if isinstance(obj, dict):
        if obj.get("__format__") == "cvd_table_v1":
            records = obj.get("records")
            columns = obj.get("columns")
            if not isinstance(records, list) or not isinstance(columns, list):
                raise ValueError(
                    f"`{filename}` has an invalid cvd_table_v1 payload"
                )
            return pd.DataFrame.from_records(records, columns=columns)
        for value in obj.values():
            if isinstance(value, pd.DataFrame):
                return value
        try:
            return pd.DataFrame(obj)
        except Exception:
            pass
    raise ValueError(f"`{filename}` is a {type(obj).__name__}, expected a DataFrame")


def _resolve(df, role):
    for name in _SINGLE_ALIASES[role]:
        if name in df.columns:
            return name
    return None


def _parse_transition(df):
    """Add integer baseline/destination columns parsed from the table."""
    if {"baseline", "destination"}.issubset(df.columns):
        return df
    col = _resolve(df, "transition")
    if col is None:
        raise ValueError("no transition/baseline/destination columns found")
    parts = (
        df[col].astype(str)
        .str.extract(r"(\d)\s*(?:->|→|_)\s*(\d)")
        .astype(float)
    )
    df["baseline"], df["destination"] = parts[0], parts[1]
    if df["baseline"].isna().any():
        raise ValueError(f"could not parse transitions from column `{col}`")
    return df


def load_single_ate():
    """Load and normalize the single-variable transition-effect table."""
    if not SINGLE_ATE_FILE.exists():
        return None, _missing(SINGLE_ATE_FILE)
    try:
        df = _as_dataframe(joblib.load(SINGLE_ATE_FILE), SINGLE_ATE_FILE.name)
        df = df.copy()
        for role in ("estimate", "low", "high"):
            col = _resolve(df, role)
            if col is None:
                raise ValueError(f"no `{role}` column (accepted: {_SINGLE_ALIASES[role]})")
            df[role] = pd.to_numeric(df[col], errors="coerce")
        if "domain" not in df.columns:
            raise ValueError("no `domain` column")
        df = _parse_transition(df)
        df["baseline"] = df["baseline"].astype(int)
        df["destination"] = df["destination"].astype(int)
        if "sig" in df.columns:
            df["sig"] = df["sig"].astype(str).str.lower().isin(["true", "1", "yes"])
        else:
            df["sig"] = (df["low"] > 0) | (df["high"] < 0)
        return df, None
    except Exception as exc:
        return None, f"`{SINGLE_ATE_FILE.name}` could not be read: {exc}"


def load_combined_ate():
    """Load the two-variable multi-arm table (multiarm_ate.csv layout)."""
    if not COMBINED_ATE_FILE.exists():
        return None, _missing(COMBINED_ATE_FILE)
    try:
        df = _as_dataframe(joblib.load(COMBINED_ATE_FILE), COMBINED_ATE_FILE.name)
        df = df.copy()
        required = {"pair", "model_type", "arm_label", "estimate", "low", "high"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"missing columns: {sorted(missing)}")
        for col in ("estimate", "low", "high", "p", "p_fdr"):
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")
        return df, None
    except Exception as exc:
        return None, f"`{COMBINED_ATE_FILE.name}` could not be read: {exc}"
