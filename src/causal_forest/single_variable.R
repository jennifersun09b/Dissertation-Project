# Run with: Rscript single_variable.R

# =============================================================================
# Refined single-variable causal-forest analysis
#
# This integrated version preserves every dissertation CSV/RDS output and also
# writes a validated Python joblib artifact for downstream application use.
#
# Main design
#   - Exact baseline-specific transitions: stayed at baseline versus moved to a
#     prespecified destination value.
#   - Binary causal forest for each transition.
#   - Full eligible-sample ATE is the main estimand for the selected cohort.
#
# Sensitivity estimands
#   - Overlap-weighted ATE.
#   - ATE after propensity trimming at 0.01-0.99, 0.05-0.95 and 0.10-0.90.
#
# Heterogeneity evidence
#   - test_calibration() differential forest prediction.
#   - best_linear_projection() for prespecified sex, age and BMI groups.
#   - RATE/AUTOC using benefit/harm-oriented CATE priorities.
#   - Subgroup ATEs with arm, event and propensity support.
#   - Participant bootstrap of direct subgroup differences.
#   - Seed-stability evaluation.
#
# CATE distribution outputs
#   - Per-participant out-of-bag CATEs with standard errors.
#   - Exact histogram counts and a fine quantile grid, so the distribution can
#     be plotted without reloading the per-participant file.
#   - Adaptive CATE quartiles, carrying explicit flags for estimates that exceed
#     or fall below what the stratum's event rate can arithmetically support,
#     plus a propensity-trimmed variant of the same quartiles.
#
# Every CSV carries at least one row. safe_write() substitutes an explicit
# status/note row when a table is empty, so a missing result is always visible as
# a statement rather than as a file containing only column names.
#
# Notes on two things that previously produced empty or unusable output
#   - tidy_inference() now strips S3 class attributes before coercing. grf and
#     lmtest return objects classed as "coeftest" and
#     "rank_average_treatment_effect"; as.data.frame() dispatches on class, so
#     those attributes stopped it reaching as.data.frame.matrix()/.list() and it
#     fell through to as.data.frame.default(), which errors. The surrounding
#     tryCatch swallowed the error, so the calibration, best-linear-projection
#     and RATE outputs were silently written as header-only files.
#   - The subgroup bootstrap now resamples the forest's doubly robust scores by
#     default instead of refitting a forest per replicate. The score bootstrap
#     is conditional on the fitted nuisance functions and is thousands of times
#     cheaper, so it can be left on. Set BOOT_REFIT=1 for the full refit.
#
# No arbitrary sample-size or event-count exclusions are applied. Every
# observed transition is attempted. Sparse, zero-event and failed models are
# retained in status/support outputs and must be interpreted accordingly.
#
# The default input is the primary cohort. To run the sensitivity cohort, use
# the same script with a different DATA_PATH, ANALYSIS_LABEL and OUTPUT_DIR.
# This keeps all modelling decisions identical.
#
# Example sensitive run
#   DATA_PATH=/path/sensitive_single_variable.csv \
#   OUTCOME_NAME=CVD_outcome ANALYSIS_LABEL=sensitive \
#   OUTPUT_DIR=/path/results_sensitive Rscript single_variable_sensitive_grf_refined.R
#
# Example focused bootstrap/stability run
#   CONTRAST_DOMAINS=diet,smoking N_BOOT=200 \
#   STABILITY_SEEDS=42,101,202 Rscript single_variable_sensitive_grf_refined.R
# =============================================================================

suppressPackageStartupMessages({
  library(grf)
  library(tidyverse)
})

options(
  width = 10000,
  tibble.width = Inf,
  tibble.print_max = Inf,
  tibble.print_min = Inf
)

# =============================================================================
# 1. Configuration
# =============================================================================

env_chr <- function(name, default = "") {
  value <- Sys.getenv(name, unset = default)
  if (!nzchar(value)) default else value
}

env_int <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (is.na(value)) default else value
}

env_bool <- function(name, default = FALSE) {
  value <- tolower(Sys.getenv(name, unset = ifelse(default, "1", "0")))
  value %in% c("1", "true", "yes", "y")
}

split_csv <- function(value) {
  value <- trimws(value)
  if (!nzchar(value)) return(character())
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

DATA_PATH <- env_chr(
  "DATA_PATH",
  "/home/rmhiund/causal_analysis/Cohort/primary_single_variable.csv"
)
OUTCOME_NAME <- env_chr("OUTCOME_NAME", "CVD_outcome")
ANALYSIS_LABEL <- env_chr("ANALYSIS_LABEL", "primary")
TEMPORAL_DESIGN <- env_chr(
  "TEMPORAL_DESIGN",
  "broad T0-to-T2 transition analysis; temporal interpretation is exploratory"
)
OUTPUT_DIR <- env_chr(
  "OUTPUT_DIR",
  file.path(getwd(), paste0("single_variable_results_", ANALYSIS_LABEL))
)

# Additional app-facing output; original dissertation outputs remain unchanged.
JOBLIB_PATH <- env_chr("JOBLIB_PATH", file.path(OUTPUT_DIR, "single_variable_ate.joblib"))
JOBLIB_BRIDGE <- tolower(env_chr("JOBLIB_BRIDGE", "auto"))
JOBLIB_COMPRESS <- env_int("JOBLIB_COMPRESS", 3)
PYTHON_BIN <- env_chr("PYTHON_BIN", "python3")
if (!JOBLIB_BRIDGE %in% c("auto", "reticulate", "python")) {
  stop("JOBLIB_BRIDGE must be auto, reticulate or python; got: ", JOBLIB_BRIDGE)
}

RANDOM_SEED <- env_int("RANDOM_SEED", 42)
NUM_TREES <- env_int("NUM_TREES", 5000)
BOOT_TREES <- env_int("BOOT_TREES", 1500)
NUM_THREADS <- env_int("NUM_THREADS", 1)
N_BOOT <- env_int("N_BOOT", 1000)
BOOT_REFIT <- env_bool("BOOT_REFIT", FALSE)
TUNE_PARAMETERS <- env_chr("TUNE_PARAMETERS", "all")
SAVE_FORESTS <- env_bool("SAVE_FORESTS", FALSE)
EXPORT_INDIVIDUAL_CATE <- env_bool("EXPORT_INDIVIDUAL_CATE", TRUE)
CATE_HIST_BINS <- env_int("CATE_HIST_BINS", 60)
STABILITY_SEEDS <- suppressWarnings(as.integer(split_csv(
  env_chr("STABILITY_SEEDS", as.character(RANDOM_SEED))
)))
STABILITY_SEEDS <- unique(STABILITY_SEEDS[!is.na(STABILITY_SEEDS)])
if (length(STABILITY_SEEDS) == 0) STABILITY_SEEDS <- RANDOM_SEED

CONTRAST_DOMAINS_RAW <- env_chr("CONTRAST_DOMAINS", "all")
EXTRA_COVARIATES <- split_csv(env_chr("EXTRA_COVARIATES", ""))

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(JOBLIB_PATH), recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
  flush.console()
}

log_msg("Analysis: ", ANALYSIS_LABEL)
log_msg("Data: ", DATA_PATH)
log_msg("Outcome: ", OUTCOME_NAME)
log_msg("CSV output directory: ", OUTPUT_DIR)
log_msg("Joblib output: ", JOBLIB_PATH, " (bridge: ", JOBLIB_BRIDGE, ")")
log_msg("Trees: ", NUM_TREES, "; bootstrap repetitions: ", N_BOOT)
log_msg(
  "Bootstrap mode: ",
  if (BOOT_REFIT) "full forest refit per replicate (slow)" else "doubly robust score resampling"
)
log_msg("Stability seeds: ", paste(STABILITY_SEEDS, collapse = ","))
if (length(STABILITY_SEEDS) < 2) {
  log_msg(
    "NOTE: only one stability seed requested, so seed_stability.csv will record ",
    "a single fit and cannot demonstrate stability. Pass ",
    "STABILITY_SEEDS=42,43,44 to test it."
  )
}
log_msg(
  "CATE distribution export: histogram (", CATE_HIST_BINS,
  " bins) and quantile grid always written; per-participant file ",
  if (EXPORT_INDIVIDUAL_CATE) {
    "ON (one row per participant per contrast, expect roughly 50-100 MB)"
  } else {
    "OFF (set EXPORT_INDIVIDUAL_CATE=1 to write it)"
  }
)

# =============================================================================
# 2. Transition registry and covariates
# =============================================================================

# health_direction describes the clinical direction of each transition. It is
# NOT the direction of the raw score: for mental health and alcohol a higher
# score is worse, so an increase is a deterioration. The interpretation_flag
# column below records the score convention for each domain.
contrasts <- tribble(
  ~domain, ~t0_col, ~t1_col, ~baseline, ~destination, ~raw_direction, ~health_direction, ~interpretation_flag,
  "sleep", "sleep_category_0.0", "sleep_t1_mod", 0, 1, "increase", "improvement", "standard transition",
  "sleep", "sleep_category_0.0", "sleep_t1_mod", 1, 2, "increase", "improvement", "standard transition",
  "sleep", "sleep_category_0.0", "sleep_t1_mod", 1, 0, "decrease", "deterioration", "standard transition",
  "sleep", "sleep_category_0.0", "sleep_t1_mod", 2, 1, "decrease", "deterioration", "standard transition",

  "diet", "diet_score_0.0", "diet_t1_mod", 0, 1, "increase", "improvement", "standard transition",
  "diet", "diet_score_0.0", "diet_t1_mod", 1, 2, "increase", "improvement", "standard transition",
  "diet", "diet_score_0.0", "diet_t1_mod", 1, 0, "decrease", "deterioration", "standard transition",
  "diet", "diet_score_0.0", "diet_t1_mod", 2, 1, "decrease", "deterioration", "standard transition",

  "pa", "physical_category_0.0", "pa_t1_mod", 0, 1, "increase", "improvement", "standard transition",
  "pa", "physical_category_0.0", "pa_t1_mod", 1, 2, "increase", "improvement", "standard transition",
  "pa", "physical_category_0.0", "pa_t1_mod", 1, 0, "decrease", "deterioration", "standard transition",
  "pa", "physical_category_0.0", "pa_t1_mod", 2, 1, "decrease", "deterioration", "standard transition",

  "mental", "mental_score_0.0", "mental_t1_mod", 0, 1, "increase", "deterioration", "higher score = worse mental health",
  "mental", "mental_score_0.0", "mental_t1_mod", 1, 2, "increase", "deterioration", "higher score = worse mental health",
  "mental", "mental_score_0.0", "mental_t1_mod", 1, 0, "decrease", "improvement", "higher score = worse mental health",
  "mental", "mental_score_0.0", "mental_t1_mod", 2, 1, "decrease", "improvement", "higher score = worse mental health",

  "alcohol", "alcohol_score_0.0", "alcohol_t1_mod", 0, 1, "increase", "deterioration", "higher score = greater alcohol burden",
  "alcohol", "alcohol_score_0.0", "alcohol_t1_mod", 1, 2, "increase", "deterioration", "higher score = greater alcohol burden",
  "alcohol", "alcohol_score_0.0", "alcohol_t1_mod", 1, 0, "decrease", "improvement", "lower score = lower alcohol burden",
  "alcohol", "alcohol_score_0.0", "alcohol_t1_mod", 2, 1, "decrease", "improvement", "lower score = lower alcohol burden",

  "smoking", "smoking_category_0.0", "smoking_t1_mod", 0, 1, "increase", "deterioration", "never to previous: possible initiation/recoding",
  "smoking", "smoking_category_0.0", "smoking_t1_mod", 1, 2, "increase", "deterioration", "previous to current: relapse",
  "smoking", "smoking_category_0.0", "smoking_t1_mod", 1, 0, "decrease", "not_causally_interpretable", "previous to never: likely reporting inconsistency",
  "smoking", "smoking_category_0.0", "smoking_t1_mod", 2, 1, "decrease", "improvement", "current to previous: cessation"
) |>
  mutate(
    contrast_id = paste(domain, baseline, destination, sep = "_"),
    transition = paste0(baseline, " -> ", destination),
    reference = paste0(baseline, " -> ", baseline),
    active = paste0(baseline, " -> ", destination),
    .before = 1
  )

valid_domains <- unique(contrasts$domain)
selected_domains <- if (tolower(CONTRAST_DOMAINS_RAW) == "all") {
  valid_domains
} else {
  requested <- split_csv(CONTRAST_DOMAINS_RAW)
  unknown <- setdiff(requested, valid_domains)
  if (length(unknown) > 0) stop("Unknown CONTRAST_DOMAINS: ", paste(unknown, collapse = ", "))
  requested
}
contrasts <- filter(contrasts, domain %in% selected_domains)

base_covariates <- c(
  "age", "sex", "bmi",
  "sleep_category_0.0", "smoking_category_0.0", "alcohol_score_0.0",
  "diet_score_0.0", "physical_category_0.0", "mental_score_0.0"
)
covariate_names <- unique(c(base_covariates, EXTRA_COVARIATES))

# =============================================================================
# 3. Load and validate data
# =============================================================================

df <- readr::read_csv(DATA_PATH, show_col_types = FALSE, progress = FALSE)

required_cols <- unique(c(
  covariate_names,
  OUTCOME_NAME,
  contrasts$t0_col,
  contrasts$t1_col
))
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

non_numeric <- covariate_names[
  !vapply(df[covariate_names], is.numeric, logical(1))
]
if (length(non_numeric) > 0) {
  stop("Encode these covariates numerically before GRF: ", paste(non_numeric, collapse = ", "))
}

X_all <- df |>
  select(all_of(covariate_names)) |>
  as.matrix()
Y_all <- as.numeric(df[[OUTCOME_NAME]])

if (!all(na.omit(unique(Y_all)) %in% c(0, 1))) {
  stop(OUTCOME_NAME, " must be coded 0/1.")
}

cohort_summary <- tibble(
  analysis_label = ANALYSIS_LABEL,
  temporal_design = TEMPORAL_DESIGN,
  n_loaded = nrow(df),
  events = sum(Y_all == 1, na.rm = TRUE),
  non_events = sum(Y_all == 0, na.rm = TRUE),
  missing_outcome = sum(is.na(Y_all)),
  complete_covariates = sum(complete.cases(X_all)),
  event_rate = mean(Y_all, na.rm = TRUE)
)

# =============================================================================
# 4. General helpers
# =============================================================================

normalise_name <- function(x) {
  x <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  gsub("^_+|_+$", "", x)
}

find_column <- function(names_vector, candidates) {
  for (candidate in candidates) {
    hit <- which(names_vector == candidate)
    if (length(hit) > 0) return(hit[1])
  }
  NA_integer_
}

# grf and lmtest return objects that carry explicit S3 class attributes:
# test_calibration() and best_linear_projection() return "coeftest" matrices and
# rank_average_treatment_effect() returns a classed list. as.data.frame()
# dispatches on class, so those attributes prevent it reaching
# as.data.frame.matrix() or as.data.frame.list(); it falls through to
# as.data.frame.default(), which errors on them. unclass() first is what makes
# these objects tidyable, and skipping it is why the calibration, BLP and RATE
# CSVs were previously written with a header and no rows.
to_inference_frame <- function(object) {
  if (is.null(object)) return(NULL)
  bare <- unclass(object)

  if (is.matrix(bare)) {
    row_names <- rownames(bare)
    dat <- as.data.frame(bare, check.names = FALSE, stringsAsFactors = FALSE)
    if (!is.null(row_names)) rownames(dat) <- row_names
    return(dat)
  }

  if (is.list(bare)) {
    # Keep scalar fields only. This drops the TOC data frame returned by
    # rank_average_treatment_effect() and any per-observation vectors, which
    # would otherwise leave the fields with unequal lengths.
    keep <- vapply(bare, function(x) is.atomic(x) && length(x) == 1, logical(1))
    if (!any(keep)) return(NULL)
    return(as.data.frame(bare[keep], check.names = FALSE, stringsAsFactors = FALSE))
  }

  if (is.atomic(bare) && !is.null(names(bare))) {
    return(as.data.frame(as.list(bare), check.names = FALSE, stringsAsFactors = FALSE))
  }

  NULL
}

tidy_inference <- function(object, source = NA_character_) {
  if (is.null(object)) {
    return(tibble(source = source, error_message = "estimator returned NULL"))
  }

  if (is.atomic(object) && !is.matrix(object) && !is.null(names(object))) {
    clean <- normalise_name(names(object))
    estimate_idx <- find_column(clean, c("estimate", "coefficient", "coef"))
    se_idx <- find_column(clean, c("std_error", "std_err", "se"))
    p_idx <- find_column(clean, c("p_value", "p", "pr_t", "pr_z"))
    if (!is.na(estimate_idx)) {
      estimate <- as.numeric(object[[estimate_idx]])
      std_error <- if (!is.na(se_idx)) as.numeric(object[[se_idx]]) else NA_real_
      statistic <- if (is.finite(std_error) && std_error > 0) estimate / std_error else NA_real_
      p_value <- if (!is.na(p_idx)) {
        as.numeric(object[[p_idx]])
      } else if (is.finite(statistic)) {
        2 * pnorm(abs(statistic), lower.tail = FALSE)
      } else {
        NA_real_
      }
      return(tibble(
        source = source,
        term = source,
        estimate = estimate,
        std_error = std_error,
        statistic = statistic,
        p_value = p_value
      ))
    }
  }

  dat <- tryCatch(to_inference_frame(object), error = function(e) NULL)
  if (is.null(dat) || nrow(dat) == 0) {
    # Return a visible diagnostic rather than a zero-row tibble. A zero-row
    # return is what safe_bind() silently drops, and it is how the calibration,
    # BLP and RATE files came to be written with a header and nothing else.
    return(tibble(
      source = source,
      error_message = paste0(
        "could not coerce an object of class ",
        paste(class(object), collapse = "/"),
        " into a tidy inference table"
      )
    ))
  }
  clean_names <- normalise_name(names(dat))
  names(dat) <- clean_names

  estimate_col <- find_column(clean_names, c("estimate", "coefficient", "coef"))
  se_col <- find_column(clean_names, c("std_error", "std_err", "se"))
  stat_col <- find_column(clean_names, c("t_value", "z_value", "statistic", "z", "t"))
  p_col <- find_column(clean_names, c("pr_t", "pr_z", "p_value", "p"))
  target_col <- find_column(clean_names, "target")

  terms <- rownames(dat)
  if (is.null(terms) || all(terms == as.character(seq_len(nrow(dat))))) {
    terms <- if (!is.na(target_col)) {
      as.character(dat[[target_col]])
    } else {
      paste0("row_", seq_len(nrow(dat)))
    }
  }

  n_rows <- nrow(dat)
  estimate <- if (!is.na(estimate_col)) as.numeric(dat[[estimate_col]]) else rep(NA_real_, n_rows)
  std_error <- if (!is.na(se_col)) as.numeric(dat[[se_col]]) else rep(NA_real_, n_rows)
  statistic <- if (!is.na(stat_col)) as.numeric(dat[[stat_col]]) else rep(NA_real_, n_rows)
  p_value <- if (!is.na(p_col)) as.numeric(dat[[p_col]]) else rep(NA_real_, n_rows)

  # RATE/AUTOC reports only an estimate and a standard error, so derive the test
  # statistic and two-sided p-value when the object does not carry them. Without
  # this the RATE family has no p_value column and its FDR step is skipped.
  needs_stat <- !is.finite(statistic) & is.finite(estimate) &
    is.finite(std_error) & std_error > 0
  statistic[needs_stat] <- estimate[needs_stat] / std_error[needs_stat]

  needs_p <- !is.finite(p_value) & is.finite(statistic)
  p_value[needs_p] <- 2 * pnorm(abs(statistic[needs_p]), lower.tail = FALSE)

  tibble(
    source = source,
    term = terms,
    estimate = estimate,
    std_error = std_error,
    statistic = statistic,
    p_value = p_value
  )
}

extract_ate <- function(object) {
  if (is.atomic(object) && !is.null(names(object))) {
    clean <- normalise_name(names(object))
    estimate_idx <- find_column(clean, c("estimate", "coefficient", "coef"))
    se_idx <- find_column(clean, c("std_error", "std_err", "se"))
    if (!is.na(estimate_idx) && !is.na(se_idx)) {
      estimate <- as.numeric(object[[estimate_idx]])
      std_error <- as.numeric(object[[se_idx]])
      z <- estimate / std_error
      return(tibble(
        estimate = estimate,
        std_error = std_error,
        low = estimate - 1.96 * std_error,
        high = estimate + 1.96 * std_error,
        z = z,
        p_value = 2 * pnorm(abs(z), lower.tail = FALSE)
      ))
    }
  }
  dat <- as.data.frame(object)
  if (!all(c("estimate", "std.err") %in% names(dat))) {
    stop("ATE output did not contain estimate and std.err.")
  }
  estimate <- as.numeric(dat$estimate[1])
  std_error <- as.numeric(dat$std.err[1])
  z <- estimate / std_error
  tibble(
    estimate = estimate,
    std_error = std_error,
    low = estimate - 1.96 * std_error,
    high = estimate + 1.96 * std_error,
    z = z,
    p_value = 2 * pnorm(abs(z), lower.tail = FALSE)
  )
}

prefix_result <- function(dat, spec) {
  if (nrow(dat) == 0) return(dat)
  dat |>
    mutate(
      analysis_label = ANALYSIS_LABEL,
      contrast_id = spec$contrast_id,
      domain = spec$domain,
      transition = spec$transition,
      baseline = spec$baseline,
      destination = spec$destination,
      health_direction = spec$health_direction,
      interpretation_flag = spec$interpretation_flag,
      .before = 1
    )
}

safe_write <- function(dat, filename, note = "no rows produced by this run") {
  # safe_bind() returns a one-column, zero-row tibble when every contrast
  # contributed nothing, which write_csv happily writes as a header-only file.
  # Checking nrow as well as ncol is what guarantees every CSV carries at least
  # one interpretable row rather than a bare column name.
  if (is.null(dat) || ncol(dat) == 0 || nrow(dat) == 0) {
    dat <- tibble(
      analysis_label = ANALYSIS_LABEL,
      status = "no_rows",
      note = note
    )
  }
  readr::write_csv(dat, file.path(OUTPUT_DIR, filename))
}

# --- joblib bridge -----------------------------------------------------------
#
# joblib is a Python, pickle-based format, so there is no native R writer. Two
# bridges are attempted, in this order:
#
#   1. reticulate, when the package is installed and its Python has joblib and
#      pandas. The table is converted in memory; no intermediate files exist.
#   2. an external interpreter (PYTHON_BIN, default python3). The table is
#      staged temporarily as CSV, converted to built-in Python records, and
#      written with joblib. The staging directory is then deleted.
#
# Both paths reload and validate the artifact before atomically replacing the
# destination. Either way the only persistent file is the .joblib.

joblib_dump_reticulate <- function(table, path, compress) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return("the reticulate package is not installed")
  }
  tryCatch({
    joblib <- reticulate::import("joblib", convert = FALSE)
    json <- reticulate::import("json", convert = FALSE)
    payload <- reticulate::r_to_py(as.data.frame(table))
    wire <- reticulate::import_builtins(convert = FALSE)$dict()
    reticulate::py_set_item(wire, "__format__", "cvd_table_v1")
    reticulate::py_set_item(
      wire, "columns", reticulate::r_to_py(as.list(names(table)))
    )
    reticulate::py_set_item(
      wire, "records", json$loads(payload$to_json(orient = "records"))
    )
    joblib$dump(wire, path, compress = as.integer(compress))

    restored <- joblib$load(path)
    restored_format <- reticulate::py_to_r(restored$get("__format__"))
    restored_records <- reticulate::py_to_r(restored$get("records"))
    if (!identical(restored_format, "cvd_table_v1")) {
      stop("reloaded joblib has the wrong wire-format marker")
    }
    if (length(restored_records) != nrow(table)) {
      stop("reloaded joblib row count does not match the exported table")
    }
    NULL
  }, error = function(e) conditionMessage(e))
}

joblib_dump_python <- function(table, path, compress, python_bin) {
  if (!nzchar(python_bin)) python_bin <- "python3"
  if (!nzchar(Sys.which(python_bin)[[1]]) && !file.exists(python_bin)) {
    return(paste0("no executable Python found at '", python_bin, "'"))
  }

  staging <- file.path(tempdir(), paste0("joblib_staging_", Sys.getpid()))
  unlink(staging, recursive = TRUE, force = TRUE)
  dir.create(staging, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(staging, recursive = TRUE, force = TRUE), add = TRUE)

  table_path <- file.path(staging, "single_variable_ate.csv")
  readr::write_csv(as.data.frame(table), table_path, na = "")

  # true_values/false_values restore the logical columns that write_csv emitted
  # as TRUE/FALSE; pandas would otherwise read them back as strings.
  script_path <- file.path(staging, "dump.py")
  writeLines(c(
    "import json, os, sys",
    "import pandas as pd",
    "import joblib",
    "source, out, compress = sys.argv[1], sys.argv[2], int(sys.argv[3])",
    "required = {'domain', 'baseline', 'destination', 'ate', 'ate_low', 'ate_high'}",
    "payload = pd.read_csv(",
    "    source, true_values=['TRUE'], false_values=['FALSE'], low_memory=False",
    ")",
    "missing = required.difference(payload.columns)",
    "if missing:",
    "    raise ValueError(f'missing app columns before dump: {sorted(missing)}')",
    "if payload.empty:",
    "    raise ValueError('single-variable payload has no estimable rows')",
    "wire = {",
    "    '__format__': 'cvd_table_v1',",
    "    'columns': payload.columns.tolist(),",
    "    'records': json.loads(payload.to_json(orient='records')),",
    "}",
    "joblib.dump(wire, out, compress=compress)",
    "restored = joblib.load(out)",
    "if restored.get('__format__') != 'cvd_table_v1':",
    "    raise ValueError('reloaded payload has the wrong wire-format marker')",
    "restored_df = pd.DataFrame.from_records(",
    "    restored['records'], columns=restored['columns']",
    ")",
    "if restored_df.shape != payload.shape:",
    "    raise ValueError(f'reloaded shape {restored_df.shape} != source shape {payload.shape}')",
    "missing = required.difference(restored_df.columns)",
    "if missing:",
    "    raise ValueError(f'missing app columns after reload: {sorted(missing)}')"
  ), script_path)

  output <- suppressWarnings(system2(
    python_bin,
    shQuote(c(script_path, table_path, path, as.character(compress))),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status")
  if (!is.null(status) && status != 0) {
    return(paste(c(paste0(python_bin, " exited with status ", status), output),
                 collapse = "\n"))
  }
  NULL
}

write_joblib <- function(table, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  problems <- character()
  temporary_path <- paste0(path, ".tmp-", Sys.getpid())
  unlink(temporary_path, force = TRUE)
  on.exit(unlink(temporary_path, force = TRUE), add = TRUE)

  if (JOBLIB_BRIDGE %in% c("auto", "reticulate")) {
    failure <- joblib_dump_reticulate(table, temporary_path, JOBLIB_COMPRESS)
    if (is.null(failure)) {
      if (!file.rename(temporary_path, path)) {
        stop("validated joblib could not replace destination: ", path)
      }
      log_msg("Wrote and verified ", path, " via reticulate.")
      return(invisible(path))
    }
    problems <- c(problems, paste0("reticulate bridge: ", failure))
    unlink(temporary_path, force = TRUE)
  }

  if (JOBLIB_BRIDGE %in% c("auto", "python")) {
    failure <- joblib_dump_python(
      table, temporary_path, JOBLIB_COMPRESS, PYTHON_BIN
    )
    if (is.null(failure)) {
      if (!file.rename(temporary_path, path)) {
        stop("validated joblib could not replace destination: ", path)
      }
      log_msg("Wrote and verified ", path, " via ", PYTHON_BIN, ".")
      return(invisible(path))
    }
    problems <- c(problems, paste0("external Python bridge: ", failure))
  }

  stop(
    "Could not write ", path, ". joblib needs Python, so install either\n",
    "  - the R package reticulate, with joblib and pandas in its Python, or\n",
    "  - a python3 on PATH (or set PYTHON_BIN) with joblib and pandas.\n",
    "Attempts:\n", paste(problems, collapse = "\n")
  )
}

safe_bind <- function(results, field) {
  if (length(results) == 0) return(tibble())
  bind_rows(lapply(results, function(x) x[[field]]), .id = "result_id")
}

add_fdr_subset <- function(dat, subset_idx, p_col = "p_value", output_col = "p_fdr") {
  if (nrow(dat) == 0 || !p_col %in% names(dat)) return(dat)
  dat[[output_col]] <- NA_real_
  valid <- subset_idx & is.finite(dat[[p_col]])
  dat[[output_col]][valid] <- p.adjust(dat[[p_col]][valid], method = "BH")
  dat
}

weighted_mean_safe <- function(x, w) {
  valid <- is.finite(x) & is.finite(w) & w >= 0
  if (!any(valid) || sum(w[valid]) <= 0) return(NA_real_)
  sum(w[valid] * x[valid]) / sum(w[valid])
}

weighted_var_safe <- function(x, w) {
  valid <- is.finite(x) & is.finite(w) & w >= 0
  if (sum(valid) < 2 || sum(w[valid]) <= 0) return(NA_real_)
  mu <- weighted_mean_safe(x[valid], w[valid])
  sum(w[valid] * (x[valid] - mu)^2) / sum(w[valid])
}

effective_sample_size <- function(w) {
  w <- w[is.finite(w) & w >= 0]
  if (length(w) == 0 || sum(w^2) == 0) return(NA_real_)
  sum(w)^2 / sum(w^2)
}

calculate_balance <- function(X, W, weights, weighting, spec) {
  bind_rows(lapply(colnames(X), function(variable) {
    x <- X[, variable]
    mean_1 <- weighted_mean_safe(x[W == 1], weights[W == 1])
    mean_0 <- weighted_mean_safe(x[W == 0], weights[W == 0])
    var_1 <- weighted_var_safe(x[W == 1], weights[W == 1])
    var_0 <- weighted_var_safe(x[W == 0], weights[W == 0])
    pooled_sd <- sqrt((var_1 + var_0) / 2)
    tibble(
      weighting = weighting,
      variable = variable,
      mean_moved = mean_1,
      mean_stayed = mean_0,
      smd = if (is.finite(pooled_sd) && pooled_sd > 0) (mean_1 - mean_0) / pooled_sd else NA_real_
    )
  })) |>
    prefix_result(spec)
}

make_subgroups <- function(X) {
  list(
    female_sex_0 = which(X[, "sex"] == 0),
    male_sex_1 = which(X[, "sex"] == 1),
    age_lt_60 = which(X[, "age"] < 60),
    age_ge_60 = which(X[, "age"] >= 60),
    bmi_lt_25 = which(X[, "bmi"] < 25),
    bmi_25_to_lt_30 = which(X[, "bmi"] >= 25 & X[, "bmi"] < 30),
    bmi_ge_30 = which(X[, "bmi"] >= 30)
  )
}

# =============================================================================
# 5. Contrast-level estimators
# =============================================================================

estimate_ate_row <- function(forest, subset, estimand, X, Y, W, spec, target = "all") {
  idx <- which(subset)
  if (length(idx) == 0) {
    return(prefix_result(tibble(
      estimand = estimand,
      n = 0,
      n_stayed = 0,
      n_moved = 0,
      events = 0,
      error_message = "empty estimand subset"
    ), spec))
  }

  result <- tryCatch({
    raw <- if (target == "overlap") {
      average_treatment_effect(forest, target.sample = "overlap")
    } else {
      average_treatment_effect(forest, subset = subset, target.sample = "all")
    }
    extract_ate(raw)
  }, error = function(e) tibble(error_message = e$message))

  result |>
    mutate(
      estimand = estimand,
      n = length(idx),
      n_stayed = sum(W[idx] == 0),
      n_moved = sum(W[idx] == 1),
      events = sum(Y[idx] == 1),
      events_stayed = sum(Y[idx][W[idx] == 0] == 1),
      events_moved = sum(Y[idx][W[idx] == 1] == 1),
      .before = 1
    ) |>
    prefix_result(spec)
}

propensity_summary <- function(propensity, W, spec) {
  bind_rows(lapply(c(0, 1), function(arm) {
    values <- propensity[W == arm]
    tibble(
      arm = arm,
      arm_label = ifelse(arm == 1, "moved", "stayed"),
      n = length(values),
      minimum = min(values, na.rm = TRUE),
      q01 = quantile(values, 0.01, na.rm = TRUE),
      q05 = quantile(values, 0.05, na.rm = TRUE),
      median = median(values, na.rm = TRUE),
      q95 = quantile(values, 0.95, na.rm = TRUE),
      q99 = quantile(values, 0.99, na.rm = TRUE),
      maximum = max(values, na.rm = TRUE)
    )
  })) |>
    prefix_result(spec)
}

subgroup_support <- function(subgroup_data, Y, W, propensity, spec) {
  groups <- make_subgroups(subgroup_data)
  bind_rows(lapply(names(groups), function(group_name) {
    group_idx <- groups[[group_name]]
    bind_rows(lapply(c(0, 1), function(arm) {
      idx <- group_idx[W[group_idx] == arm]
      p <- propensity[idx]
      tibble(
        subgroup = group_name,
        arm = arm,
        arm_label = ifelse(arm == 1, "moved", "stayed"),
        n = length(idx),
        events = sum(Y[idx] == 1),
        event_rate = if (length(idx) > 0) mean(Y[idx]) else NA_real_,
        propensity_q01 = if (length(p) > 0) quantile(p, 0.01, na.rm = TRUE) else NA_real_,
        propensity_q05 = if (length(p) > 0) quantile(p, 0.05, na.rm = TRUE) else NA_real_,
        propensity_median = if (length(p) > 0) median(p, na.rm = TRUE) else NA_real_,
        propensity_q95 = if (length(p) > 0) quantile(p, 0.95, na.rm = TRUE) else NA_real_,
        support_flag = case_when(
          length(idx) == 0 ~ "empty subgroup-arm cell",
          sum(Y[idx] == 1) == 0 ~ "zero events; estimate may be unstable",
          TRUE ~ "events observed"
        )
      )
    }))
  })) |>
    prefix_result(spec)
}

subgroup_ates <- function(forest, subgroup_data, propensity, spec) {
  groups <- make_subgroups(subgroup_data)
  estimands <- list(
    full_sample = rep(TRUE, nrow(subgroup_data)),
    trim_005_095 = propensity > 0.05 & propensity < 0.95
  )

  bind_rows(lapply(names(estimands), function(estimand_name) {
    estimand_mask <- estimands[[estimand_name]]
    bind_rows(lapply(names(groups), function(group_name) {
      idx <- intersect(groups[[group_name]], which(estimand_mask))
      if (length(idx) == 0) {
        return(tibble(
          estimand = estimand_name,
          subgroup = group_name,
          n = 0,
          error_message = "empty subgroup"
        ))
      }
      tryCatch({
        estimate <- extract_ate(average_treatment_effect(forest, subset = idx))
        estimate |>
          mutate(
            estimand = estimand_name,
            subgroup = group_name,
            n = length(idx),
            inference = "subgroup ATE; differences tested with BLP or bootstrap",
            .before = 1
          )
      }, error = function(e) tibble(
        estimand = estimand_name,
        subgroup = group_name,
        n = length(idx),
        error_message = e$message
      ))
    }))
  })) |>
    prefix_result(spec)
}

SUBGROUP_CONTRASTS <- list(
  list(high = "male_sex_1", low = "female_sex_0", label = "male - female"),
  list(high = "age_ge_60", low = "age_lt_60", label = "age >=60 - age <60"),
  list(high = "bmi_ge_30", low = "bmi_lt_25", label = "BMI >=30 - BMI <25")
)

summarise_bootstrap_draws <- function(values, modifier, observed, method) {
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    return(tibble(
      modifier = modifier,
      method = method,
      n_success = 0L,
      estimate = observed,
      bootstrap_mean = NA_real_,
      bootstrap_se = NA_real_,
      low = NA_real_,
      high = NA_real_,
      p_value = NA_real_
    ))
  }
  p_value <- 2 * min(mean(values <= 0), mean(values >= 0))
  tibble(
    modifier = modifier,
    method = method,
    n_success = length(values),
    estimate = observed,
    bootstrap_mean = mean(values),
    bootstrap_se = sd(values),
    low = as.numeric(quantile(values, 0.025, names = FALSE)),
    high = as.numeric(quantile(values, 0.975, names = FALSE)),
    p_value = min(1, p_value)
  )
}

# Default path: resample participants and recompute each subgroup difference from
# the forest's doubly robust scores. average_treatment_effect(target.sample =
# "all") is the mean of those scores over the requested subset, so a subgroup
# difference is a difference of score means and the bootstrap is exact for that
# quantity, conditional on the fitted nuisance functions. No forest is refitted,
# which is why this can be left on by default.
bootstrap_scores <- function(forest, subgroup_data, spec) {
  scores <- tryCatch(as.numeric(get_scores(forest)), error = function(e) NULL)
  if (is.null(scores) || length(scores) != nrow(subgroup_data)) {
    return(prefix_result(tibble(
      modifier = NA_character_,
      method = "doubly_robust_score_bootstrap",
      n_success = 0L,
      error_message = "doubly robust scores unavailable"
    ), spec))
  }

  groups <- make_subgroups(subgroup_data)
  n <- length(scores)
  membership <- matrix(
    FALSE, nrow = n, ncol = length(groups),
    dimnames = list(NULL, names(groups))
  )
  for (group_name in names(groups)) {
    membership[groups[[group_name]], group_name] <- TRUE
  }

  difference_vector <- function(idx) {
    s <- scores[idx]
    member <- membership[idx, , drop = FALSE]
    group_mean <- function(group_name) {
      values <- s[member[, group_name]]
      if (length(values) == 0) NA_real_ else mean(values)
    }
    vapply(SUBGROUP_CONTRASTS, function(cmp) {
      group_mean(cmp$high) - group_mean(cmp$low)
    }, numeric(1))
  }

  observed <- difference_vector(seq_len(n))

  set.seed(RANDOM_SEED + 8000)
  draws <- matrix(NA_real_, nrow = N_BOOT, ncol = length(SUBGROUP_CONTRASTS))
  for (b in seq_len(N_BOOT)) {
    draws[b, ] <- difference_vector(sample.int(n, n, replace = TRUE))
  }

  bind_rows(lapply(seq_along(SUBGROUP_CONTRASTS), function(k) {
    summarise_bootstrap_draws(
      draws[, k], SUBGROUP_CONTRASTS[[k]]$label, observed[k],
      "doubly_robust_score_bootstrap"
    )
  })) |>
    mutate(
      n_bootstrap = N_BOOT,
      inference = paste(
        "participant bootstrap of doubly robust scores;",
        "conditional on the fitted nuisance functions"
      )
    ) |>
    prefix_result(spec)
}

# Opt-in path (BOOT_REFIT=1): refit a complete forest per replicate. This also
# propagates uncertainty in the nuisance functions, at roughly N_BOOT times the
# cost of the main fit.
bootstrap_refit <- function(X, subgroup_data, Y, W, spec) {
  n <- nrow(X)
  set.seed(RANDOM_SEED + 8000)

  draws <- bind_rows(lapply(seq_len(N_BOOT), function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    X_b <- X[idx, , drop = FALSE]
    subgroup_b <- subgroup_data[idx, , drop = FALSE]
    Y_b <- Y[idx]
    W_b <- W[idx]
    if (length(unique(W_b)) < 2 || length(unique(Y_b)) < 2) return(tibble())

    forest_b <- tryCatch(
      causal_forest(
        X = X_b,
        Y = Y_b,
        W = W_b,
        num.trees = BOOT_TREES,
        honesty = TRUE,
        ci.group.size = 2,
        num.threads = NUM_THREADS,
        seed = RANDOM_SEED + b
      ),
      error = function(e) NULL
    )
    if (is.null(forest_b)) return(tibble())

    groups <- make_subgroups(subgroup_b)
    get_estimate <- function(group_name) {
      group_idx <- groups[[group_name]]
      if (length(group_idx) == 0) return(NA_real_)
      tryCatch(
        extract_ate(average_treatment_effect(forest_b, subset = group_idx))$estimate[1],
        error = function(e) NA_real_
      )
    }

    values <- vapply(SUBGROUP_CONTRASTS, function(cmp) {
      get_estimate(cmp$high) - get_estimate(cmp$low)
    }, numeric(1))
    tibble(bootstrap = b, comparison = seq_along(values), difference = values)
  }))

  if (nrow(draws) == 0) {
    return(prefix_result(tibble(
      modifier = NA_character_,
      method = "forest_refit_bootstrap",
      n_success = 0L,
      error_message = "no bootstrap replicate produced an estimable forest"
    ), spec))
  }

  bind_rows(lapply(seq_along(SUBGROUP_CONTRASTS), function(k) {
    values <- draws$difference[draws$comparison == k]
    summarise_bootstrap_draws(
      values, SUBGROUP_CONTRASTS[[k]]$label, mean(values[is.finite(values)]),
      "forest_refit_bootstrap"
    )
  })) |>
    mutate(
      n_bootstrap = N_BOOT,
      inference = "participant bootstrap; complete forest refitted per replicate"
    ) |>
    prefix_result(spec)
}

bootstrap_subgroup_differences <- function(forest, X, subgroup_data, Y, W, spec) {
  if (N_BOOT <= 0) {
    return(prefix_result(tibble(
      modifier = NA_character_,
      method = if (BOOT_REFIT) "forest_refit_bootstrap" else "doubly_robust_score_bootstrap",
      n_success = 0L,
      n_bootstrap = 0L,
      inference = "not run: N_BOOT = 0"
    ), spec))
  }
  if (BOOT_REFIT) {
    bootstrap_refit(X, subgroup_data, Y, W, spec)
  } else {
    bootstrap_scores(forest, subgroup_data, spec)
  }
}

empty_result <- function(spec, status, reason, cohort_flow, arm_support) {
  list(
    status = prefix_result(tibble(status = status, reason = reason), spec),
    cohort_flow = cohort_flow,
    arm_support = arm_support,
    propensity = tibble(),
    weight_summary = tibble(),
    balance = tibble(),
    ate = tibble(),
    calibration = tibble(),
    blp = tibble(),
    rate = tibble(),
    cate_summary = tibble(),
    cate_quartiles = tibble(),
    # These three must be listed here as well as in run_transition(). safe_bind()
    # pulls one named field out of every result, so a field that is absent for
    # non-estimable contrasts is silently skipped rather than filled with NA.
    cate_individual = tibble(),
    cate_histogram = tibble(),
    cate_quantiles = tibble(),
    variable_importance = tibble(),
    subgroup_support = tibble(),
    subgroup_ate = tibble(),
    bootstrap_difference = tibble(),
    stability = tibble(),
    forest = NULL
  )
}

# =============================================================================
# 6. Run one exact transition
# =============================================================================

run_transition <- function(spec) {
  t0 <- df[[spec$t0_col]]
  t1 <- df[[spec$t1_col]]

  stayed <- !is.na(t0) & !is.na(t1) & t0 == spec$baseline & t1 == spec$baseline
  moved <- !is.na(t0) & !is.na(t1) & t0 == spec$baseline & t1 == spec$destination
  transition_eligible <- stayed | moved
  complete <- transition_eligible & complete.cases(X_all) & !is.na(Y_all)

  X_complete <- X_all[complete, , drop = FALSE]
  Y <- Y_all[complete]
  W <- as.integer(moved[complete])

  varies <- apply(X_complete, 2, function(x) length(unique(x[is.finite(x)])) > 1)
  dropped_covariates <- colnames(X_complete)[!varies]
  X <- X_complete[, varies, drop = FALSE]
  subgroup_data <- X_complete[, c("age", "sex", "bmi"), drop = FALSE]

  cohort_flow <- prefix_result(tibble(
    n_loaded = nrow(df),
    n_with_transition_values = sum(!is.na(t0) & !is.na(t1)),
    n_transition_eligible = sum(transition_eligible),
    n_complete_analysis = sum(complete),
    n_excluded_missing_covariate_or_outcome = sum(transition_eligible) - sum(complete),
    dropped_constant_covariates = paste(dropped_covariates, collapse = ",")
  ), spec)

  arm_support <- prefix_result(tibble(
    arm = c(0, 1),
    arm_label = c("stayed", "moved"),
    n = c(sum(W == 0), sum(W == 1)),
    events = c(sum(Y[W == 0] == 1), sum(Y[W == 1] == 1)),
    non_events = c(sum(Y[W == 0] == 0), sum(Y[W == 1] == 0)),
    event_rate = c(mean(Y[W == 0]), mean(Y[W == 1])),
    support_flag = c(
      ifelse(sum(Y[W == 0] == 1) == 0, "zero events; model may fail", "events observed"),
      ifelse(sum(Y[W == 1] == 1) == 0, "zero events; model may fail", "events observed")
    )
  ), spec)

  if (nrow(X) == 0 || length(unique(W)) < 2) {
    return(empty_result(spec, "not_estimable", "reference or moved arm absent", cohort_flow, arm_support))
  }
  if (length(unique(Y)) < 2) {
    return(empty_result(spec, "not_estimable", "outcome has no variation", cohort_flow, arm_support))
  }

  set.seed(RANDOM_SEED)
  forest <- tryCatch(
    causal_forest(
      X = X,
      Y = Y,
      W = W,
      num.trees = NUM_TREES,
      honesty = TRUE,
      ci.group.size = 2,
      tune.parameters = TUNE_PARAMETERS,
      num.threads = NUM_THREADS,
      seed = RANDOM_SEED
    ),
    error = function(e) e
  )

  if (inherits(forest, "error")) {
    return(empty_result(spec, "fit_failed", forest$message, cohort_flow, arm_support))
  }

  propensity <- as.numeric(forest$W.hat)
  trim_001 <- propensity > 0.01 & propensity < 0.99
  trim_005 <- propensity > 0.05 & propensity < 0.95
  trim_010 <- propensity > 0.10 & propensity < 0.90

  propensity_table <- propensity_summary(propensity, W, spec)
  overlap_weights <- ifelse(W == 1, 1 - propensity, propensity)
  weight_summary <- prefix_result(bind_rows(lapply(c(0, 1), function(arm) {
    tibble(
      weighting = "overlap",
      arm = arm,
      n = sum(W == arm),
      weight_sum = sum(overlap_weights[W == arm]),
      effective_sample_size = effective_sample_size(overlap_weights[W == arm])
    )
  })), spec)

  balance <- bind_rows(
    calculate_balance(X, W, rep(1, length(W)), "unweighted", spec),
    calculate_balance(X, W, overlap_weights, "overlap_weighted", spec)
  )

  ate <- bind_rows(
    estimate_ate_row(forest, rep(TRUE, length(W)), "full_eligible_sample", X, Y, W, spec),
    estimate_ate_row(forest, rep(TRUE, length(W)), "overlap_weighted", X, Y, W, spec, target = "overlap"),
    estimate_ate_row(forest, trim_001, "trim_001_099", X, Y, W, spec),
    estimate_ate_row(forest, trim_005, "trim_005_095", X, Y, W, spec),
    estimate_ate_row(forest, trim_010, "trim_010_090", X, Y, W, spec)
  ) |>
    mutate(
      retained_fraction = n / length(W),
      inference = case_when(
        estimand == "full_eligible_sample" ~ "full eligible-sample ATE",
        TRUE ~ "sensitivity estimand"
      )
    )

  calibration <- tryCatch(
    prefix_result(tidy_inference(test_calibration(forest), "test_calibration"), spec),
    error = function(e) prefix_result(tibble(error_message = e$message), spec)
  )

  modifier_matrix <- cbind(
    male_sex_1 = as.numeric(subgroup_data[, "sex"] == 1),
    age_ge_60 = as.numeric(subgroup_data[, "age"] >= 60),
    bmi_25_to_lt_30 = as.numeric(
      subgroup_data[, "bmi"] >= 25 & subgroup_data[, "bmi"] < 30
    ),
    bmi_ge_30 = as.numeric(subgroup_data[, "bmi"] >= 30)
  )
  modifier_varies <- apply(
    modifier_matrix, 2,
    function(x) length(unique(x[is.finite(x)])) > 1
  )
  modifier_matrix <- modifier_matrix[, modifier_varies, drop = FALSE]

  blp <- if (ncol(modifier_matrix) == 0) {
    prefix_result(tibble(error_message = "prespecified modifiers had no variation"), spec)
  } else {
    tryCatch(
      prefix_result(tidy_inference(
        best_linear_projection(forest, A = modifier_matrix),
        "best_linear_projection"
      ), spec),
      error = function(e) prefix_result(tibble(error_message = e$message), spec)
    )
  }

  prediction <- tryCatch(
    predict(forest, estimate.variance = TRUE),
    error = function(e) NULL
  )
  cate <- if (!is.null(prediction)) as.numeric(prediction$predictions) else rep(NA_real_, nrow(X))
  cate_se <- if (!is.null(prediction)) sqrt(as.numeric(prediction$variance.estimates)) else rep(NA_real_, nrow(X))

  priorities <- if (spec$health_direction == "improvement") {
    -cate
  } else if (spec$health_direction == "deterioration") {
    cate
  } else {
    abs(cate)
  }

  rate <- tryCatch(
    prefix_result(tidy_inference(
      rank_average_treatment_effect(
        forest,
        priorities = priorities,
        target = "AUTOC"
      ),
      "rank_average_treatment_effect"
    ), spec),
    error = function(e) prefix_result(tibble(error_message = e$message), spec)
  )

  cate_summary <- prefix_result(tibble(
    n = sum(is.finite(cate)),
    mean_cate = mean(cate, na.rm = TRUE),
    sd_cate = sd(cate, na.rm = TRUE),
    min_cate = min(cate, na.rm = TRUE),
    max_cate = max(cate, na.rm = TRUE),
    mean_individual_se = mean(cate_se, na.rm = TRUE),
    inference = "descriptive out-of-bag CATE distribution; not individual causal proof"
  ), spec)

  # ---------------------------------------------------------------------------
  # CATE distribution exports.
  #
  # Previously only the six summary statistics above were written, which is not
  # enough to draw the distribution. The histogram and quantile grid below make
  # the shape reproducible from a small file; the per-participant file is written
  # for anything that needs the raw draws (it is large, so it is switchable).
  # ---------------------------------------------------------------------------
  finite_cate <- cate[is.finite(cate)]

  cate_histogram <- tibble()
  cate_quantiles <- tibble()
  if (length(finite_cate) >= 2 && diff(range(finite_cate)) > 0) {
    breaks <- seq(min(finite_cate), max(finite_cate), length.out = CATE_HIST_BINS + 1)
    bin_index <- cut(finite_cate, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    counts <- tabulate(bin_index, nbins = CATE_HIST_BINS)
    bin_width <- diff(breaks)
    cate_histogram <- prefix_result(tibble(
      bin = seq_len(CATE_HIST_BINS),
      bin_low = breaks[-length(breaks)],
      bin_high = breaks[-1],
      bin_mid = (breaks[-length(breaks)] + breaks[-1]) / 2,
      count = counts,
      proportion = counts / length(finite_cate),
      density = counts / (length(finite_cate) * bin_width),
      n_total = length(finite_cate)
    ), spec)

    probs <- c(0.001, 0.005, seq(0.01, 0.99, by = 0.01), 0.995, 0.999)
    cate_quantiles <- prefix_result(tibble(
      probability = probs,
      cate = as.numeric(quantile(finite_cate, probs, names = FALSE)),
      n_total = length(finite_cate)
    ), spec)
  }

  cate_individual <- tibble()
  if (EXPORT_INDIVIDUAL_CATE) {
    # Deliberately NOT prefix_result(): this table has one row per participant
    # per contrast (of the order of 700k rows in total), and the eight repeated
    # label columns would roughly double the file for no information. contrast_id
    # joins it back to every other output.
    cate_individual <- tibble(
      contrast_id = spec$contrast_id,
      row_id = seq_along(cate),
      moved = as.integer(W),
      cvd_event = as.integer(Y),
      propensity = as.numeric(propensity),
      cate = cate,
      cate_se = cate_se
    )
  }

  # ---------------------------------------------------------------------------
  # Adaptive CATE quartiles.
  #
  # These are exploratory and can be arithmetically impossible: quartiles are
  # formed by ranking on out-of-bag predictions and then averaging doubly robust
  # scores inside those same data-defined bins, so extreme inverse-propensity
  # weights pile up in the tail quartiles and can push a risk difference past the
  # stratum's entire event rate. The flags below make that visible in the file
  # instead of leaving it to be discovered downstream, and the trimmed variant
  # reuses the fitted forest (no refit) to show how much of it is propensity tail.
  # ---------------------------------------------------------------------------
  stratum_event_rate <- mean(Y, na.rm = TRUE)

  quartile_rows <- function(mask, estimand_name) {
    # !is.na(mask) matters: the trimming mask is built from propensity, and an NA
    # there would otherwise reach a subscripted assignment, which errors in R.
    eligible <- is.finite(cate) & !is.na(mask) & mask
    if (sum(eligible) < 4) return(tibble())
    quartile <- rep(NA_integer_, length(cate))
    quartile[eligible] <- ntile(cate[eligible], 4)
    bind_rows(lapply(1:4, function(q) {
      idx <- which(quartile == q)
      tryCatch({
        extract_ate(average_treatment_effect(forest, subset = idx)) |>
          mutate(
            estimand = estimand_name,
            quartile = q,
            n = length(idx),
            inference = "exploratory adaptive CATE quartile",
            .before = 1
          )
      }, error = function(e) tibble(
        estimand = estimand_name,
        quartile = q,
        n = length(idx),
        error_message = e$message
      ))
    }))
  }

  cate_quartiles <- bind_rows(
    quartile_rows(rep(TRUE, length(cate)), "full_sample"),
    quartile_rows(propensity >= 0.05 & propensity <= 0.95, "trim_005_095")
  )
  if (nrow(cate_quartiles) > 0) {
    # Guard: if every quartile errored there is no estimate column to test.
    if (!"estimate" %in% names(cate_quartiles)) {
      cate_quartiles$estimate <- NA_real_
    }
    cate_quartiles <- cate_quartiles |>
      mutate(
        stratum_event_rate = stratum_event_rate,
        exceeds_stratum_risk = is.finite(estimate) &
          abs(estimate) > stratum_event_rate,
        arithmetically_impossible = is.finite(estimate) &
          estimate < -stratum_event_rate
      ) |>
      prefix_result(spec)
  }

  variable_importance_table <- tryCatch(
    prefix_result(tibble(
      variable = colnames(X),
      importance = as.numeric(variable_importance(forest)),
      inference = "splitting importance; not a causal effect or formal modifier test"
    ) |>
      arrange(desc(importance)), spec),
    error = function(e) prefix_result(tibble(error_message = e$message), spec)
  )

  subgroup_support_table <- subgroup_support(
    subgroup_data, Y, W, propensity, spec
  )
  subgroup_ate_table <- subgroup_ates(
    forest, subgroup_data, propensity, spec
  )
  bootstrap_difference <- bootstrap_subgroup_differences(
    forest, X, subgroup_data, Y, W, spec
  )

  stability <- bind_rows(lapply(STABILITY_SEEDS, function(stability_seed) {
    stability_forest <- if (stability_seed == RANDOM_SEED) {
      forest
    } else {
      tryCatch(
        causal_forest(
          X = X,
          Y = Y,
          W = W,
          num.trees = NUM_TREES,
          honesty = TRUE,
          ci.group.size = 2,
          tune.parameters = TUNE_PARAMETERS,
          num.threads = NUM_THREADS,
          seed = stability_seed
        ),
        error = function(e) NULL
      )
    }
    if (is.null(stability_forest)) {
      return(prefix_result(tibble(seed = stability_seed, status = "fit_failed"), spec))
    }

    stability_cate <- as.numeric(predict(stability_forest)$predictions)
    stability_ate <- tryCatch(
      extract_ate(average_treatment_effect(stability_forest))$estimate[1],
      error = function(e) NA_real_
    )
    stability_calibration <- tryCatch(
      tidy_inference(test_calibration(stability_forest), "test_calibration"),
      error = function(e) tibble()
    )
    differential <- if (
      nrow(stability_calibration) > 0 &&
      "term" %in% names(stability_calibration)
    ) {
      stability_calibration |>
        filter(grepl("differential", term, ignore.case = TRUE))
    } else {
      tibble()
    }

    prefix_result(tibble(
      seed = stability_seed,
      n_seeds_requested = length(STABILITY_SEEDS),
      is_reference_seed = stability_seed == RANDOM_SEED,
      status = "fitted",
      ate = stability_ate,
      cate_sd = sd(stability_cate, na.rm = TRUE),
      cate_correlation_with_reference_seed = cor(cate, stability_cate, use = "complete.obs"),
      differential_calibration_estimate = if (nrow(differential) > 0) differential$estimate[1] else NA_real_,
      differential_calibration_p = if (nrow(differential) > 0) differential$p_value[1] else NA_real_,
      inference = if (length(STABILITY_SEEDS) < 2) {
        "single seed: records the reference fit only, does not demonstrate stability"
      } else {
        "refit under an independent seed; compare ate and cate_sd across seeds"
      }
    ), spec)
  }))

  list(
    status = prefix_result(tibble(status = "fitted", reason = NA_character_), spec),
    cohort_flow = cohort_flow,
    arm_support = arm_support,
    propensity = propensity_table,
    weight_summary = weight_summary,
    balance = balance,
    ate = ate,
    calibration = calibration,
    blp = blp,
    rate = rate,
    cate_summary = cate_summary,
    cate_quartiles = cate_quartiles,
    cate_individual = cate_individual,
    cate_histogram = cate_histogram,
    cate_quantiles = cate_quantiles,
    variable_importance = variable_importance_table,
    subgroup_support = subgroup_support_table,
    subgroup_ate = subgroup_ate_table,
    bootstrap_difference = bootstrap_difference,
    stability = stability,
    forest = forest
  )
}

# =============================================================================
# 7. Run all transitions without arbitrary support thresholds
# =============================================================================

all_results <- list()
for (i in seq_len(nrow(contrasts))) {
  spec <- contrasts[i, ]
  log_msg(
    "Contrast ", spec$contrast_id, ": ", spec$transition,
    " (", spec$health_direction, ")"
  )
  all_results[[spec$contrast_id]] <- run_transition(spec)
}

# =============================================================================
# 8. Consolidate outputs and control multiplicity
# =============================================================================

model_status <- safe_bind(all_results, "status")
contrast_flow <- safe_bind(all_results, "cohort_flow")
arm_support_table <- safe_bind(all_results, "arm_support")
propensity_table <- safe_bind(all_results, "propensity")
weight_table <- safe_bind(all_results, "weight_summary")
balance_table <- safe_bind(all_results, "balance")
ate_table <- safe_bind(all_results, "ate")
calibration_table <- safe_bind(all_results, "calibration")
blp_table <- safe_bind(all_results, "blp")
rate_table <- safe_bind(all_results, "rate")
cate_summary_table <- safe_bind(all_results, "cate_summary")
cate_quartile_table <- safe_bind(all_results, "cate_quartiles")
cate_hist_table <- safe_bind(all_results, "cate_histogram")
cate_quantile_table <- safe_bind(all_results, "cate_quantiles")
cate_individual_table <- safe_bind(all_results, "cate_individual")
importance_table <- safe_bind(all_results, "variable_importance")
subgroup_support_table <- safe_bind(all_results, "subgroup_support")
subgroup_ate_table <- safe_bind(all_results, "subgroup_ate")
bootstrap_table <- safe_bind(all_results, "bootstrap_difference")
stability_table <- safe_bind(all_results, "stability")

if (nrow(ate_table) > 0 && all(c("estimand", "p_value") %in% names(ate_table))) {
  ate_table <- add_fdr_subset(
    ate_table,
    ate_table$estimand == "full_eligible_sample",
    "p_value", "p_fdr_full_sample_family"
  )
}

if (nrow(calibration_table) > 0 && all(c("term", "p_value") %in% names(calibration_table))) {
  calibration_table <- add_fdr_subset(
    calibration_table,
    grepl("differential", calibration_table$term, ignore.case = TRUE),
    "p_value", "p_fdr_differential_family"
  )
}

if (nrow(blp_table) > 0 && all(c("term", "p_value") %in% names(blp_table))) {
  blp_table <- add_fdr_subset(
    blp_table,
    !grepl("intercept", blp_table$term, ignore.case = TRUE),
    "p_value", "p_fdr_modifier_family"
  )
}

if (nrow(rate_table) > 0 && "p_value" %in% names(rate_table)) {
  rate_table <- add_fdr_subset(
    rate_table,
    rep(TRUE, nrow(rate_table)),
    "p_value", "p_fdr_rate_family"
  )
}

if (nrow(bootstrap_table) > 0 && "p_value" %in% names(bootstrap_table)) {
  bootstrap_table <- add_fdr_subset(
    bootstrap_table,
    rep(TRUE, nrow(bootstrap_table)),
    "p_value", "p_fdr_subgroup_family"
  )
}

# =============================================================================
# 9. Save dissertation-ready outputs
# =============================================================================

safe_write(cohort_summary, "cohort_summary.csv")
safe_write(contrasts, "transition_registry.csv")
safe_write(model_status, "model_status.csv")
safe_write(contrast_flow, "contrast_participant_flow.csv")
safe_write(arm_support_table, "arm_counts_events.csv")
safe_write(propensity_table, "propensity_summary.csv")
safe_write(weight_table, "overlap_weight_effective_sample_size.csv")
safe_write(balance_table, "covariate_balance.csv")
safe_write(ate_table, "ate_full_overlap_and_trimmed_estimands.csv")
safe_write(calibration_table, "heterogeneity_test_calibration.csv",
  "test_calibration() produced no rows for any contrast")
safe_write(blp_table, "heterogeneity_best_linear_projection.csv",
  "best_linear_projection() produced no rows for any contrast")
safe_write(rate_table, "heterogeneity_rate_autoc.csv",
  "rank_average_treatment_effect() produced no rows for any contrast")
safe_write(cate_summary_table, "cate_summary_descriptive.csv")
safe_write(cate_hist_table, "cate_distribution_histogram.csv",
  "no contrast produced a CATE distribution with more than one distinct value")
safe_write(cate_quantile_table, "cate_distribution_quantiles.csv",
  "no contrast produced a CATE distribution with more than one distinct value")
safe_write(cate_quartile_table, "cate_quartiles_exploratory.csv",
  "no contrast had at least four participants with finite CATEs")
if (EXPORT_INDIVIDUAL_CATE) {
  safe_write(cate_individual_table, "cate_individual_predictions.csv",
    "no contrast produced per-participant CATEs")
}
safe_write(importance_table, "variable_importance_descriptive.csv")
safe_write(subgroup_support_table, "subgroup_arm_event_overlap_support.csv")
safe_write(subgroup_ate_table, "subgroup_ate.csv")
safe_write(bootstrap_table, "bootstrap_subgroup_differences.csv",
  "subgroup bootstrap did not run; check N_BOOT in run_configuration.csv")
safe_write(stability_table, "seed_stability.csv",
  "no forest was refitted; check STABILITY_SEEDS in run_configuration.csv")

if (SAVE_FORESTS) {
  saveRDS(all_results, file.path(OUTPUT_DIR, "all_single_variable_forests.rds"))
}

run_configuration <- tibble(
  setting = c(
    "DATA_PATH", "OUTCOME_NAME", "ANALYSIS_LABEL", "TEMPORAL_DESIGN",
    "OUTPUT_DIR", "RANDOM_SEED", "NUM_TREES", "BOOT_TREES",
    "NUM_THREADS", "N_BOOT", "BOOT_REFIT", "TUNE_PARAMETERS", "SAVE_FORESTS",
    "EXPORT_INDIVIDUAL_CATE", "CATE_HIST_BINS",
    "STABILITY_SEEDS", "CONTRAST_DOMAINS", "EXTRA_COVARIATES",
    "JOBLIB_PATH", "JOBLIB_BRIDGE", "JOBLIB_COMPRESS", "PYTHON_BIN"
  ),
  value = c(
    DATA_PATH, OUTCOME_NAME, ANALYSIS_LABEL, TEMPORAL_DESIGN,
    OUTPUT_DIR, RANDOM_SEED, NUM_TREES, BOOT_TREES,
    NUM_THREADS, N_BOOT, BOOT_REFIT, TUNE_PARAMETERS, SAVE_FORESTS,
    EXPORT_INDIVIDUAL_CATE, CATE_HIST_BINS,
    paste(STABILITY_SEEDS, collapse = ","), CONTRAST_DOMAINS_RAW,
    paste(EXTRA_COVARIATES, collapse = ","),
    JOBLIB_PATH, JOBLIB_BRIDGE, JOBLIB_COMPRESS, PYTHON_BIN
  )
)
safe_write(run_configuration, "run_configuration.csv")
capture.output(sessionInfo(), file = file.path(OUTPUT_DIR, "session_info.txt"))


# App-facing joblib payload
column_or_na <- function(dat, name, fill = NA_real_) {
  if (name %in% names(dat)) dat[[name]] else rep(fill, nrow(dat))
}

# The table the web app reads: one row per transition, using the full eligible
# sample estimand for the selected cohort. ate/ate_low/ate_high are
# risk differences as fractions; the app multiplies by 100 for percentage points.
build_webapp_ate <- function(dat) {
  if (is.null(dat) || nrow(dat) == 0 ||
      !all(c("estimand", "domain") %in% names(dat))) {
    return(tibble())
  }
  keep <- !is.na(dat$estimand) & dat$estimand == "full_eligible_sample"
  primary <- dat[keep, , drop = FALSE]
  if (nrow(primary) == 0) return(tibble())

  estimate <- as.numeric(column_or_na(primary, "estimate"))
  low <- as.numeric(column_or_na(primary, "low"))
  high <- as.numeric(column_or_na(primary, "high"))

  tibble(
    analysis_label = ANALYSIS_LABEL,
    domain = as.character(column_or_na(primary, "domain", NA_character_)),
    transition = as.character(column_or_na(primary, "transition", NA_character_)),
    baseline = as.integer(column_or_na(primary, "baseline")),
    destination = as.integer(column_or_na(primary, "destination")),
    health = as.character(column_or_na(primary, "health_direction", NA_character_)),
    contrast_id = as.character(column_or_na(primary, "contrast_id", NA_character_)),
    n = column_or_na(primary, "n"),
    n_stayed = column_or_na(primary, "n_stayed"),
    n_moved = column_or_na(primary, "n_moved"),
    events = column_or_na(primary, "events"),
    ate = estimate,
    ate_low = low,
    ate_high = high,
    # Same rule the app falls back to when sig is absent: the 95% interval does
    # not cross zero.
    sig = is.finite(low) & is.finite(high) & (low > 0 | high < 0),
    p_value = column_or_na(primary, "p_value"),
    p_fdr = column_or_na(primary, "p_fdr_full_sample_family"),
    interpretation_flag = as.character(
      column_or_na(primary, "interpretation_flag", NA_character_)
    ),
    estimand = "full_eligible_sample"
  )
}

validate_webapp_ate <- function(dat) {
  required <- c(
    "domain", "transition", "baseline", "destination",
    "ate", "ate_low", "ate_high"
  )
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0) {
    stop("App-facing ATE table is missing columns: ", paste(missing, collapse = ", "))
  }
  if (nrow(dat) == 0) {
    stop("No transition effects were available for the app joblib.")
  }
  if (any(is.na(dat$domain) | !nzchar(trimws(dat$domain)))) {
    stop("App-facing ATE table contains a missing or blank domain.")
  }
  numeric_columns <- c("ate", "ate_low", "ate_high")
  if (any(!vapply(dat[numeric_columns], is.numeric, logical(1)))) {
    stop("ate, ate_low, and ate_high must be numeric.")
  }
  if (any(!is.finite(as.matrix(dat[numeric_columns])))) {
    stop("App-facing ATE table contains non-finite estimates or intervals.")
  }
  if (any(dat$ate_low > dat$ate | dat$ate > dat$ate_high)) {
    stop("At least one estimate lies outside its confidence interval.")
  }
  if (any(duplicated(dat[c("domain", "baseline", "destination")]))) {
    stop("Duplicate domain/baseline/destination rows found in app-facing ATE table.")
  }
  dat |>
    arrange(domain, baseline, destination)
}

app_ate <- validate_webapp_ate(build_webapp_ate(ate_table))
log_msg(
  "App payload: ", nrow(app_ate), " rows; ", ncol(app_ate),
  " columns; domains=", paste(sort(unique(app_ate$domain)), collapse = ",")
)
write_joblib(app_ate, JOBLIB_PATH)
log_msg("Completed ", length(all_results), " transition analyses.")
log_msg("CSV results written to: ", OUTPUT_DIR)
log_msg("Joblib written to: ", JOBLIB_PATH)
