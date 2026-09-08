#!/usr/bin/env Rscript

# =============================================================================
# Refined two-variable causal-forest workflow
#
# This integrated version preserves every dissertation CSV/RDS output and also
# writes a validated Python joblib artifact for downstream application use.
#
# Purpose
#   1. Estimate overall paired lifestyle effects with multi-arm causal forests.
#   2. Evaluate treatment-effect heterogeneity for each active arm versus the
#      stable reference with a separate binary causal forest.
#   3. Keep every observed arm. Sparse or unsupported analyses are attempted
#      and reported as failures rather than removed by arbitrary thresholds.
#   4. Support sensitive and strict-outcome runs through environment variables.
#
# Main heterogeneity evidence
#   - test_calibration() differential forest prediction
#   - best_linear_projection() for prespecified effect modifiers
#   - RATE/AUTOC for binary contrasts
#   - subgroup ATEs with arm/event/propensity support summaries
#   - optional full-model bootstrap of subgroup differences
#   - optional seed-stability analysis
#
# Exploratory only
#   - CATE ranges and quartiles
#   - multi-arm variable importance
#   - additive interaction point estimates without bootstrap uncertainty
#
# Example sensitive run
#   DATA_PATH=/path/sensitive.csv OUTCOME_NAME=CVD_outcome \
#   ANALYSIS_LABEL=sensitive Rscript combined_variable4_heterogeneity_refined.R
#
# Example strict post-T2 run
#   DATA_PATH=/path/strict.csv OUTCOME_NAME=CVD_post_T2 \
#   ANALYSIS_LABEL=strict Rscript combined_variable4_heterogeneity_refined.R
#
# Optional focused heterogeneity run
#   HET_PAIRS=diet_smoking,sleep_pa HET_MODEL_TYPES=improvement \
#   N_BOOT=200 STABILITY_SEEDS=42,101,202 Rscript ...
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

env_num <- function(name, default) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
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
OUTPUT_DIR <- env_chr(
  "OUTPUT_DIR",
  file.path(getwd(), paste0("two_variable_results_", ANALYSIS_LABEL))
)

# Additional app-facing output; original dissertation outputs remain unchanged.
JOBLIB_PATH <- env_chr("JOBLIB_PATH", file.path(OUTPUT_DIR, "combined_variable_ate.joblib"))
JOBLIB_BRIDGE <- tolower(env_chr("JOBLIB_BRIDGE", "auto"))
JOBLIB_COMPRESS <- env_int("JOBLIB_COMPRESS", 3)
PYTHON_BIN <- env_chr("PYTHON_BIN", "python3")
if (!JOBLIB_BRIDGE %in% c("auto", "reticulate", "python")) {
  stop("JOBLIB_BRIDGE must be auto, reticulate or python; got: ", JOBLIB_BRIDGE)
}

RANDOM_SEED <- env_int("RANDOM_SEED", 42)
NUM_THREADS <- env_int("NUM_THREADS", 1)
MULTIARM_TREES <- env_int("MULTIARM_TREES", 3000)
HET_TREES <- env_int("HET_TREES", 4000)
BOOT_TREES <- env_int("BOOT_TREES", 1500)
N_BOOT <- env_int("N_BOOT", 0)
RUN_MULTIARM <- env_bool("RUN_MULTIARM", TRUE)
RUN_HETEROGENEITY <- env_bool("RUN_HETEROGENEITY", TRUE)
SAVE_FORESTS <- env_bool("SAVE_FORESTS", FALSE)
ENFORCE_STRUCTURAL_FEASIBILITY <- env_bool(
  "ENFORCE_STRUCTURAL_FEASIBILITY", TRUE
)
HET_TUNE_PARAMETERS <- env_chr("HET_TUNE_PARAMETERS", "all")

HET_PAIRS_RAW <- env_chr("HET_PAIRS", "all")
HET_MODEL_TYPES <- split_csv(env_chr(
  "HET_MODEL_TYPES", "improvement,deterioration,mixed"
))
HET_ACTIVE_ARMS_RAW <- env_chr("HET_ACTIVE_ARMS", "all")
STABILITY_SEEDS <- suppressWarnings(as.integer(split_csv(
  env_chr("STABILITY_SEEDS", as.character(RANDOM_SEED))
)))
STABILITY_SEEDS <- unique(STABILITY_SEEDS[!is.na(STABILITY_SEEDS)])
if (length(STABILITY_SEEDS) == 0) STABILITY_SEEDS <- RANDOM_SEED

EXTRA_COVARIATES <- split_csv(env_chr("EXTRA_COVARIATES", ""))

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(JOBLIB_PATH), recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
  flush.console()
}

log_msg("Analysis label: ", ANALYSIS_LABEL)
log_msg("Data: ", DATA_PATH)
log_msg("Outcome: ", OUTCOME_NAME)
log_msg("CSV output directory: ", OUTPUT_DIR)
log_msg("Joblib output: ", JOBLIB_PATH, " (bridge: ", JOBLIB_BRIDGE, ")")
log_msg("Multi-arm trees: ", MULTIARM_TREES, "; binary trees: ", HET_TREES)
log_msg("Bootstrap repetitions: ", N_BOOT)

# A numerical increase is healthier for mental health, as confirmed for this
# project's score. Alcohol and smoking use the opposite direction.
lifestyle_vars <- list(
  sleep = list(
    t0 = "sleep_category_0.0", t1 = "sleep_t1_mod", health_sign = 1
  ),
  diet = list(
    t0 = "diet_score_0.0", t1 = "diet_t1_mod", health_sign = 1
  ),
  pa = list(
    t0 = "physical_category_0.0", t1 = "pa_t1_mod", health_sign = 1
  ),
  mental = list(
    t0 = "mental_score_0.0", t1 = "mental_t1_mod", health_sign = -1
  ),
  alcohol = list(
    t0 = "alcohol_score_0.0", t1 = "alcohol_t1_mod", health_sign = -1
  ),
  smoking = list(
    t0 = "smoking_category_0.0", t1 = "smoking_t1_mod", health_sign = -1
  )
)

base_covariates <- c(
  "age", "sex", "bmi",
  "sleep_category_0.0", "smoking_category_0.0", "alcohol_score_0.0",
  "diet_score_0.0", "physical_category_0.0", "mental_score_0.0"
)
covariate_names <- unique(c(base_covariates, EXTRA_COVARIATES))

valid_model_types <- c("improvement", "deterioration", "mixed")
bad_model_types <- setdiff(HET_MODEL_TYPES, valid_model_types)
if (length(bad_model_types) > 0) {
  stop("Unknown HET_MODEL_TYPES: ", paste(bad_model_types, collapse = ", "))
}

# =============================================================================
# 2. Load and validate data
# =============================================================================

df <- readr::read_csv(DATA_PATH, show_col_types = FALSE, progress = FALSE)

required_cols <- unique(c(
  covariate_names,
  OUTCOME_NAME,
  unlist(lapply(lifestyle_vars, function(x) c(x$t0, x$t1)))
))
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

non_numeric_covariates <- covariate_names[
  !vapply(df[covariate_names], is.numeric, logical(1))
]
if (length(non_numeric_covariates) > 0) {
  stop(
    "GRF covariates must be numeric. Encode these variables first: ",
    paste(non_numeric_covariates, collapse = ", ")
  )
}

X <- df |>
  select(all_of(covariate_names)) |>
  as.matrix()
Y <- as.numeric(df[[OUTCOME_NAME]])

if (!all(na.omit(unique(Y)) %in% c(0, 1))) {
  stop(OUTCOME_NAME, " must be coded 0/1.")
}

cohort_summary <- tibble(
  analysis_label = ANALYSIS_LABEL,
  n_loaded = nrow(df),
  events = sum(Y == 1, na.rm = TRUE),
  non_events = sum(Y == 0, na.rm = TRUE),
  missing_outcome = sum(is.na(Y)),
  complete_covariates = sum(complete.cases(X)),
  event_rate = mean(Y, na.rm = TRUE)
)
readr::write_csv(cohort_summary, file.path(OUTPUT_DIR, "cohort_summary.csv"))
log_msg("Loaded N=", nrow(df), "; events=", sum(Y == 1, na.rm = TRUE))

# =============================================================================
# 3. General helpers
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

tidy_inference <- function(object, source = NA_character_) {
  if (is.null(object)) return(tibble())

  if (is.atomic(object) && !is.null(names(object))) {
    clean_vector_names <- normalise_name(names(object))
    estimate_idx <- find_column(clean_vector_names, c("estimate", "coefficient", "coef"))
    se_idx <- find_column(clean_vector_names, c("std_error", "std_err", "se"))
    p_idx <- find_column(clean_vector_names, c("p_value", "p", "pr_t", "pr_z"))
    if (!is.na(estimate_idx)) {
      estimate <- as.numeric(object[[estimate_idx]])
      std_error <- if (!is.na(se_idx)) as.numeric(object[[se_idx]]) else NA_real_
      return(tibble(
        source = source,
        term = source,
        estimate = estimate,
        std_error = std_error,
        statistic = if (is.finite(std_error) && std_error > 0) estimate / std_error else NA_real_,
        p_value = if (!is.na(p_idx)) as.numeric(object[[p_idx]]) else if (
          is.finite(std_error) && std_error > 0
        ) 2 * pnorm(abs(estimate / std_error), lower.tail = FALSE) else NA_real_
      ))
    }
  }

  dat <- tryCatch(as.data.frame(object, check.names = FALSE), error = function(e) NULL)
  if (is.null(dat) || nrow(dat) == 0) return(tibble())

  original_names <- names(dat)
  clean_names <- normalise_name(original_names)
  names(dat) <- clean_names

  estimate_col <- find_column(clean_names, c("estimate", "coefficient", "coef"))
  se_col <- find_column(clean_names, c("std_error", "std_err", "se"))
  stat_col <- find_column(clean_names, c("t_value", "z_value", "statistic", "z", "t"))
  p_col <- find_column(clean_names, c(
    "pr_t", "pr_z", "p_value", "p", "pr_t_value", "pr_z_value"
  ))

  terms <- rownames(dat)
  if (is.null(terms) || all(terms == as.character(seq_len(nrow(dat))))) {
    terms <- paste0("row_", seq_len(nrow(dat)))
  }

  tibble(
    source = source,
    term = terms,
    estimate = if (!is.na(estimate_col)) as.numeric(dat[[estimate_col]]) else NA_real_,
    std_error = if (!is.na(se_col)) as.numeric(dat[[se_col]]) else NA_real_,
    statistic = if (!is.na(stat_col)) as.numeric(dat[[stat_col]]) else NA_real_,
    p_value = if (!is.na(p_col)) as.numeric(dat[[p_col]]) else NA_real_
  )
}

extract_ate <- function(object) {
  if (is.atomic(object) && !is.null(names(object))) {
    clean_vector_names <- normalise_name(names(object))
    estimate_idx <- find_column(clean_vector_names, c("estimate", "coefficient", "coef"))
    se_idx <- find_column(clean_vector_names, c("std_error", "std_err", "se"))
    if (!is.na(estimate_idx) && !is.na(se_idx)) {
      estimate <- as.numeric(object[[estimate_idx]])
      std_error <- as.numeric(object[[se_idx]])
      return(tibble(
        contrast = "1 - 0",
        estimate = estimate,
        std.err = std_error,
        low = estimate - 1.96 * std_error,
        high = estimate + 1.96 * std_error,
        z = estimate / std_error,
        p = 2 * pnorm(abs(z), lower.tail = FALSE)
      ))
    }
  }

  dat <- as.data.frame(object)
  if (nrow(dat) == 0) return(tibble())
  if (!"contrast" %in% names(dat)) {
    dat <- tibble::rownames_to_column(dat, "contrast")
  }
  if (!"estimate" %in% names(dat) || !"std.err" %in% names(dat)) {
    stop("ATE output did not contain estimate and std.err columns.")
  }
  dat |>
    mutate(
      estimate = as.numeric(estimate),
      std.err = as.numeric(std.err),
      low = estimate - 1.96 * std.err,
      high = estimate + 1.96 * std.err,
      z = estimate / std.err,
      p = 2 * pnorm(abs(z), lower.tail = FALSE)
    )
}

extract_single_ate <- function(object) {
  dat <- extract_ate(object)
  if (nrow(dat) == 0) return(c(estimate = NA_real_, std_error = NA_real_))
  c(estimate = dat$estimate[1], std_error = dat$std.err[1])
}

add_fdr <- function(dat, p_col = "p_value", output_col = "p_fdr") {
  if (nrow(dat) == 0 || !p_col %in% names(dat)) return(dat)
  p <- dat[[p_col]]
  adjusted <- rep(NA_real_, length(p))
  valid <- is.finite(p)
  adjusted[valid] <- p.adjust(p[valid], method = "BH")
  dat[[output_col]] <- adjusted
  dat
}

safe_bind <- function(results, field) {
  if (length(results) == 0) return(tibble())
  bind_rows(lapply(results, function(x) x[[field]]), .id = "result_id")
}

safe_write <- function(dat, filename) {
  if (is.null(dat)) dat <- tibble()
  if (ncol(dat) == 0) dat <- tibble(note = "no rows")
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
#      written as a CSV under tempdir(), a short generated script reads it and
#      dumps the DataFrame, and the staging directory is deleted.
#
# Both bridges reload and validate the artifact before it atomically replaces
# the destination. Either way the only file left behind is the .joblib.

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

  table_path <- file.path(staging, "combined_variable_ate.csv")
  readr::write_csv(as.data.frame(table), table_path, na = "")

  # true_values/false_values restore the logical columns that write_csv emitted
  # as TRUE/FALSE; pandas would otherwise read them back as strings.
  script_path <- file.path(staging, "dump.py")
  writeLines(c(
    "import json, os, sys",
    "import pandas as pd",
    "import joblib",
    "source, out, compress = sys.argv[1], sys.argv[2], int(sys.argv[3])",
    "required = {'pair', 'model_type', 'arm_label', 'estimate', 'low', 'high'}",
    "payload = pd.read_csv(",
    "    source, true_values=['TRUE'], false_values=['FALSE'], low_memory=False",
    ")",
    "missing = required.difference(payload.columns)",
    "if missing:",
    "    raise ValueError(f'missing app columns before dump: {sorted(missing)}')",
    "if payload.empty:",
    "    raise ValueError('combined-variable payload has no estimable rows')",
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

# =============================================================================
# 4. Treatment construction
# =============================================================================

get_health_change <- function(t0, t1, health_sign) {
  health_delta <- (t1 - t0) * health_sign
  case_when(
    is.na(t0) | is.na(t1) ~ NA_character_,
    health_delta > 0 ~ "improve",
    health_delta < 0 ~ "deteriorate",
    TRUE ~ "stable"
  )
}

get_room_to_move <- function(t0, t1, health_sign) {
  values <- c(t0, t1)
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    return(list(
      can_improve = rep(FALSE, length(t0)),
      can_deteriorate = rep(FALSE, length(t0))
    ))
  }

  observed_range <- range(values)
  can_increase <- !is.na(t0) & t0 < observed_range[2]
  can_decrease <- !is.na(t0) & t0 > observed_range[1]

  if (health_sign == 1) {
    list(can_improve = can_increase, can_deteriorate = can_decrease)
  } else {
    list(can_improve = can_decrease, can_deteriorate = can_increase)
  }
}

make_pair_treatment <- function(df, var1, var2, model_type) {
  v1 <- lifestyle_vars[[var1]]
  v2 <- lifestyle_vars[[var2]]

  a0 <- df[[v1$t0]]
  a1 <- df[[v1$t1]]
  b0 <- df[[v2$t0]]
  b1 <- df[[v2$t1]]

  chg1 <- get_health_change(a0, a1, v1$health_sign)
  chg2 <- get_health_change(b0, b1, v2$health_sign)
  room1 <- get_room_to_move(a0, a1, v1$health_sign)
  room2 <- get_room_to_move(b0, b1, v2$health_sign)

  W <- rep(NA_integer_, nrow(df))

  if (model_type == "improvement") {
    W[chg1 == "stable" & chg2 == "stable"] <- 0
    W[chg1 == "improve" & chg2 == "stable"] <- 1
    W[chg1 == "stable" & chg2 == "improve"] <- 2
    W[chg1 == "improve" & chg2 == "improve"] <- 3
    if (ENFORCE_STRUCTURAL_FEASIBILITY) {
      eligible <- room1$can_improve & room2$can_improve
      W[!eligible] <- NA_integer_
    }
  } else if (model_type == "deterioration") {
    W[chg1 == "stable" & chg2 == "stable"] <- 0
    W[chg1 == "deteriorate" & chg2 == "stable"] <- 1
    W[chg1 == "stable" & chg2 == "deteriorate"] <- 2
    W[chg1 == "deteriorate" & chg2 == "deteriorate"] <- 3
    if (ENFORCE_STRUCTURAL_FEASIBILITY) {
      eligible <- room1$can_deteriorate & room2$can_deteriorate
      W[!eligible] <- NA_integer_
    }
  } else if (model_type == "mixed") {
    W[chg1 == "stable" & chg2 == "stable"] <- 0
    W[chg1 == "improve" & chg2 == "deteriorate"] <- 1
    W[chg1 == "deteriorate" & chg2 == "improve"] <- 2
    if (ENFORCE_STRUCTURAL_FEASIBILITY) {
      eligible <- room1$can_improve & room1$can_deteriorate &
        room2$can_improve & room2$can_deteriorate
      W[!eligible] <- NA_integer_
    }
  } else {
    stop("Unknown model type: ", model_type)
  }

  W
}

make_arm_labels <- function(var1, var2, model_type) {
  reference <- paste(var1, "stable +", var2, "stable")
  switch(
    model_type,
    improvement = c(
      "0" = reference,
      "1" = paste(var1, "improves only"),
      "2" = paste(var2, "improves only"),
      "3" = paste(var1, "+", var2, "both improve")
    ),
    deterioration = c(
      "0" = reference,
      "1" = paste(var1, "deteriorates only"),
      "2" = paste(var2, "deteriorates only"),
      "3" = paste(var1, "+", var2, "both deteriorate")
    ),
    mixed = c(
      "0" = reference,
      "1" = paste(var1, "improves +", var2, "deteriorates"),
      "2" = paste(var1, "deteriorates +", var2, "improves")
    )
  )
}

arm_summary <- function(W, Y, labels, pair_name, model_type) {
  levels_present <- sort(unique(W[!is.na(W)]))
  if (length(levels_present) == 0) return(tibble())

  bind_rows(lapply(levels_present, function(arm_value) {
    idx <- which(W == arm_value)
    tibble(
      analysis_label = ANALYSIS_LABEL,
      pair = pair_name,
      model_type = model_type,
      arm = as.character(arm_value),
      arm_label = unname(labels[as.character(arm_value)]),
      n = length(idx),
      events = sum(Y[idx] == 1, na.rm = TRUE),
      non_events = sum(Y[idx] == 0, na.rm = TRUE),
      event_rate = mean(Y[idx], na.rm = TRUE),
      support_flag = case_when(
        sum(Y[idx] == 1, na.rm = TRUE) == 0 ~
          "zero events: attempted; inference may fail",
        sum(Y[idx] == 0, na.rm = TRUE) == 0 ~
          "zero non-events: attempted; inference may fail",
        TRUE ~ "events and non-events observed"
      )
    )
  }))
}

overlap_summary <- function(propensity, pair_name, model_type, contrast = NA_character_) {
  if (is.null(propensity)) return(tibble())
  if (is.null(dim(propensity))) propensity <- matrix(propensity, ncol = 1)
  probability_names <- colnames(propensity)
  if (is.null(probability_names)) {
    probability_names <- paste0("probability_", seq_len(ncol(propensity)))
  }

  as_tibble(propensity, .name_repair = "minimal") |>
    setNames(probability_names) |>
    pivot_longer(everything(), names_to = "propensity_component", values_to = "propensity") |>
    group_by(propensity_component) |>
    summarise(
      analysis_label = ANALYSIS_LABEL,
      pair = pair_name,
      model_type = model_type,
      contrast = contrast,
      n = sum(is.finite(propensity)),
      minimum = min(propensity, na.rm = TRUE),
      q01 = quantile(propensity, 0.01, na.rm = TRUE),
      q05 = quantile(propensity, 0.05, na.rm = TRUE),
      median = median(propensity, na.rm = TRUE),
      q95 = quantile(propensity, 0.95, na.rm = TRUE),
      q99 = quantile(propensity, 0.99, na.rm = TRUE),
      maximum = max(propensity, na.rm = TRUE),
      .groups = "drop"
    )
}

# =============================================================================
# 5. Multi-arm models for overall paired effects
# =============================================================================

empty_multi_result <- function(pair_name, model_type, status, reason, arm_table = tibble()) {
  list(
    status = tibble(
      analysis_label = ANALYSIS_LABEL,
      pair = pair_name,
      model_type = model_type,
      status = status,
      reason = reason
    ),
    arm_table = arm_table,
    overlap = tibble(),
    ate = tibble(),
    cate_summary = tibble(),
    variable_importance = tibble(),
    interaction = tibble(),
    forest = NULL
  )
}

fit_multiarm_model <- function(var1, var2, model_type) {
  pair_name <- paste(var1, var2, sep = "_")
  labels <- make_arm_labels(var1, var2, model_type)
  W_raw <- make_pair_treatment(df, var1, var2, model_type)

  complete <- !is.na(W_raw) & !is.na(Y) & complete.cases(X)
  W_e <- W_raw[complete]
  X_e <- X[complete, , drop = FALSE]
  Y_e <- Y[complete]

  arm_tab <- arm_summary(W_e, Y_e, labels, pair_name, model_type)
  observed <- sort(unique(W_e))

  if (!(0 %in% observed)) {
    return(empty_multi_result(pair_name, model_type, "not_estimable", "reference arm absent", arm_tab))
  }
  if (length(observed) < 2) {
    return(empty_multi_result(pair_name, model_type, "not_estimable", "no comparison arm", arm_tab))
  }
  if (length(unique(Y_e)) < 2) {
    return(empty_multi_result(pair_name, model_type, "not_estimable", "outcome has no variation", arm_tab))
  }

  levels_ordered <- c(0, setdiff(observed, 0))
  W_fac <- factor(as.character(W_e), levels = as.character(levels_ordered))

  set.seed(RANDOM_SEED)
  forest <- tryCatch(
    multi_arm_causal_forest(
      X = X_e,
      Y = Y_e,
      W = W_fac,
      num.trees = MULTIARM_TREES,
      honesty = TRUE,
      ci.group.size = 2,
      num.threads = NUM_THREADS,
      seed = RANDOM_SEED
    ),
    error = function(e) e
  )

  if (inherits(forest, "error")) {
    return(empty_multi_result(pair_name, model_type, "fit_failed", forest$message, arm_tab))
  }

  ate <- tryCatch(
    extract_ate(average_treatment_effect(forest)) |>
      mutate(
        analysis_label = ANALYSIS_LABEL,
        pair = pair_name,
        model_type = model_type,
        arm = sub(" -.*", "", contrast),
        arm_label = unname(labels[arm]),
        .before = 1
      ),
    error = function(e) tibble(
      analysis_label = ANALYSIS_LABEL,
      pair = pair_name,
      model_type = model_type,
      error_message = e$message
    )
  )

  pred <- tryCatch(predict(forest, estimate.variance = TRUE), error = function(e) NULL)
  cate_summary <- tibble()
  if (!is.null(pred)) {
    contrast_names <- paste(setdiff(levels(W_fac), "0"), "- 0")
    tau <- drop(pred$predictions)
    expected_values <- nrow(X_e) * length(contrast_names)

    if (length(tau) != expected_values) {
      cate_summary <- tibble(
        analysis_label = ANALYSIS_LABEL,
        pair = pair_name,
        model_type = model_type,
        error_message = paste0(
          "Unexpected prediction dimensions: ",
          paste(dim(pred$predictions), collapse = " x "),
          "; expected ", nrow(X_e), " x ", length(contrast_names)
        )
      )
    } else {
      # GRF may return n x contrasts x 1; drop the singleton outcome
      # dimension and restore a stable two-dimensional contrast matrix.
      tau <- matrix(
        tau,
        nrow = nrow(X_e),
        ncol = length(contrast_names),
        dimnames = list(NULL, contrast_names)
      )
      cate_summary <- bind_rows(lapply(seq_along(contrast_names), function(j) {
        arm <- sub(" -.*", "", contrast_names[j])
        values <- tau[, j]
        tibble(
          analysis_label = ANALYSIS_LABEL,
          pair = pair_name,
          model_type = model_type,
          contrast = contrast_names[j],
          arm = arm,
          arm_label = unname(labels[arm]),
          n = sum(is.finite(values)),
          mean_cate = mean(values, na.rm = TRUE),
          sd_cate = sd(values, na.rm = TRUE),
          min_cate = min(values, na.rm = TRUE),
          max_cate = max(values, na.rm = TRUE),
          inference = "descriptive multi-arm CATE distribution"
        )
      }))
    }
  }

  vi <- tryCatch({
    raw <- variable_importance(forest)
    values <- if (is.matrix(raw)) rowMeans(raw) else as.numeric(raw)
    tibble(
      analysis_label = ANALYSIS_LABEL,
      pair = pair_name,
      model_type = model_type,
      variable = colnames(X_e),
      importance = values,
      inference = "splitting importance; not a causal effect or formal modifier test"
    ) |>
      arrange(desc(importance))
  }, error = function(e) tibble(
    analysis_label = ANALYSIS_LABEL,
    pair = pair_name,
    model_type = model_type,
    error_message = e$message
  ))

  interaction <- tibble()
  if (model_type %in% c("improvement", "deterioration") &&
      all(c("arm", "estimate") %in% names(ate))) {
    estimates <- split(ate$estimate, ate$arm)
    if (all(c("1", "2", "3") %in% names(estimates))) {
      interaction <- tibble(
        analysis_label = ANALYSIS_LABEL,
        pair = pair_name,
        model_type = model_type,
        interaction = "arm3 - arm1 - arm2",
        estimate = estimates[["3"]][1] - estimates[["1"]][1] - estimates[["2"]][1],
        inference = "point estimate only; bootstrap the full multi-arm model for inference"
      )
    }
  }

  list(
    status = tibble(
      analysis_label = ANALYSIS_LABEL,
      pair = pair_name,
      model_type = model_type,
      status = "fitted",
      reason = NA_character_
    ),
    arm_table = arm_tab,
    overlap = overlap_summary(forest$W.hat, pair_name, model_type),
    ate = ate,
    cate_summary = cate_summary,
    variable_importance = vi,
    interaction = interaction,
    forest = forest
  )
}

# =============================================================================
# 6. Binary contrasts for formal heterogeneity evaluation
# =============================================================================

make_subgroups <- function(X_binary) {
  list(
    female_sex_0 = which(X_binary[, "sex"] == 0),
    male_sex_1 = which(X_binary[, "sex"] == 1),
    age_lt_60 = which(X_binary[, "age"] < 60),
    age_ge_60 = which(X_binary[, "age"] >= 60)
  )
}

binary_subgroup_support <- function(
    X_b, Y_b, W_b, propensity, pair_name, model_type, contrast, arm_label
) {
  groups <- make_subgroups(X_b)
  bind_rows(lapply(names(groups), function(group_name) {
    group_idx <- groups[[group_name]]
    bind_rows(lapply(c(0, 1), function(treatment_value) {
      idx <- group_idx[W_b[group_idx] == treatment_value]
      p <- propensity[idx]
      tibble(
        analysis_label = ANALYSIS_LABEL,
        pair = pair_name,
        model_type = model_type,
        contrast = contrast,
        arm_label = arm_label,
        subgroup = group_name,
        treatment = treatment_value,
        n = length(idx),
        events = sum(Y_b[idx] == 1, na.rm = TRUE),
        event_rate = if (length(idx) > 0) mean(Y_b[idx], na.rm = TRUE) else NA_real_,
        propensity_q01 = if (length(p) > 0) quantile(p, 0.01, na.rm = TRUE) else NA_real_,
        propensity_q05 = if (length(p) > 0) quantile(p, 0.05, na.rm = TRUE) else NA_real_,
        propensity_median = if (length(p) > 0) median(p, na.rm = TRUE) else NA_real_,
        propensity_q95 = if (length(p) > 0) quantile(p, 0.95, na.rm = TRUE) else NA_real_,
        support_flag = case_when(
          length(idx) == 0 ~ "empty subgroup-treatment cell",
          sum(Y_b[idx] == 1, na.rm = TRUE) == 0 ~ "zero events; inference may be unstable",
          TRUE ~ "events observed"
        )
      )
    }))
  }))
}

estimate_subgroup_ates <- function(forest, X_b, pair_name, model_type, contrast, arm_label) {
  groups <- make_subgroups(X_b)
  bind_rows(lapply(names(groups), function(group_name) {
    idx <- groups[[group_name]]
    if (length(idx) == 0) return(tibble())
    tryCatch({
      est <- extract_single_ate(average_treatment_effect(forest, subset = idx))
      tibble(
        analysis_label = ANALYSIS_LABEL,
        pair = pair_name,
        model_type = model_type,
        contrast = contrast,
        arm_label = arm_label,
        subgroup = group_name,
        n = length(idx),
        estimate = est[["estimate"]],
        std_error = est[["std_error"]],
        low = estimate - 1.96 * std_error,
        high = estimate + 1.96 * std_error,
        inference = "subgroup ATE; compare subgroups using BLP or full bootstrap"
      )
    }, error = function(e) tibble(
      analysis_label = ANALYSIS_LABEL,
      pair = pair_name,
      model_type = model_type,
      contrast = contrast,
      arm_label = arm_label,
      subgroup = group_name,
      n = length(idx),
      error_message = e$message
    ))
  }))
}

bootstrap_subgroup_differences <- function(
    X_b, Y_b, W_b, pair_name, model_type, contrast, arm_label
) {
  if (N_BOOT <= 0) return(tibble())
  n <- nrow(X_b)
  set.seed(RANDOM_SEED + 9000)
  bootstrap_indices <- replicate(N_BOOT, sample.int(n, n, replace = TRUE), simplify = FALSE)

  draws <- bind_rows(lapply(seq_along(bootstrap_indices), function(b) {
    idx <- bootstrap_indices[[b]]
    X_boot <- X_b[idx, , drop = FALSE]
    Y_boot <- Y_b[idx]
    W_boot <- W_b[idx]

    if (length(unique(W_boot)) < 2 || length(unique(Y_boot)) < 2) return(tibble())

    forest_boot <- tryCatch(
      causal_forest(
        X = X_boot,
        Y = Y_boot,
        W = W_boot,
        num.trees = BOOT_TREES,
        honesty = TRUE,
        ci.group.size = 2,
        num.threads = NUM_THREADS,
        seed = RANDOM_SEED + b
      ),
      error = function(e) NULL
    )
    if (is.null(forest_boot)) return(tibble())

    groups <- make_subgroups(X_boot)
    get_est <- function(group_name) {
      group_idx <- groups[[group_name]]
      if (length(group_idx) == 0) return(NA_real_)
      tryCatch(
        extract_single_ate(
          average_treatment_effect(forest_boot, subset = group_idx)
        )[["estimate"]],
        error = function(e) NA_real_
      )
    }

    female <- get_est("female_sex_0")
    male <- get_est("male_sex_1")
    younger <- get_est("age_lt_60")
    older <- get_est("age_ge_60")

    tibble(
      bootstrap = b,
      sex_difference = male - female,
      age_difference = older - younger
    )
  }))

  summarise_draw <- function(values, modifier) {
    values <- values[is.finite(values)]
    if (length(values) == 0) {
      return(tibble(
        modifier = modifier,
        n_success = 0,
        estimate = NA_real_,
        low = NA_real_,
        high = NA_real_,
        p_value = NA_real_
      ))
    }
    p_boot <- 2 * min(mean(values <= 0), mean(values >= 0))
    tibble(
      modifier = modifier,
      n_success = length(values),
      estimate = mean(values),
      low = quantile(values, 0.025),
      high = quantile(values, 0.975),
      p_value = min(1, p_boot)
    )
  }

  bind_rows(
    summarise_draw(draws$sex_difference, "male - female"),
    summarise_draw(draws$age_difference, "age >=60 - age <60")
  ) |>
    mutate(
      analysis_label = ANALYSIS_LABEL,
      pair = pair_name,
      model_type = model_type,
      contrast = contrast,
      arm_label = arm_label,
      inference = "participant bootstrap; complete binary forest refitted",
      .before = 1
    )
}

empty_binary_result <- function(pair_name, model_type, contrast, arm_label, status, reason, support) {
  list(
    status = tibble(
      analysis_label = ANALYSIS_LABEL,
      pair = pair_name,
      model_type = model_type,
      contrast = contrast,
      arm_label = arm_label,
      status = status,
      reason = reason
    ),
    support = support,
    overall_ate = tibble(),
    calibration = tibble(),
    blp = tibble(),
    rate = tibble(),
    subgroup_ate = tibble(),
    subgroup_support = tibble(),
    cate_summary = tibble(),
    cate_quartiles = tibble(),
    variable_importance = tibble(),
    bootstrap_difference = tibble(),
    stability = tibble(),
    forest = NULL
  )
}

fit_binary_contrast <- function(var1, var2, model_type, active_arm) {
  pair_name <- paste(var1, var2, sep = "_")
  contrast <- paste(active_arm, "- 0")
  labels <- make_arm_labels(var1, var2, model_type)
  arm_label <- unname(labels[as.character(active_arm)])
  W_raw <- make_pair_treatment(df, var1, var2, model_type)

  complete <- !is.na(W_raw) & !is.na(Y) & complete.cases(X)
  keep <- complete & W_raw %in% c(0, active_arm)
  X_b <- X[keep, , drop = FALSE]
  Y_b <- Y[keep]
  W_b <- as.integer(W_raw[keep] == active_arm)

  support <- arm_summary(
    W_b, Y_b,
    c("0" = labels[["0"]], "1" = arm_label),
    pair_name, model_type
  ) |>
    mutate(
      original_arm = if_else(arm == "1", as.character(active_arm), "0"),
      contrast = contrast
    )

  observed_original <- unique(W_raw[keep])
  if (!(0 %in% observed_original) || !(active_arm %in% observed_original)) {
    return(empty_binary_result(
      pair_name, model_type, contrast, arm_label,
      "not_estimable", "reference or active arm absent", support
    ))
  }
  if (length(unique(Y_b)) < 2) {
    return(empty_binary_result(
      pair_name, model_type, contrast, arm_label,
      "not_estimable", "outcome has no variation", support
    ))
  }

  set.seed(RANDOM_SEED)
  forest <- tryCatch(
    causal_forest(
      X = X_b,
      Y = Y_b,
      W = W_b,
      num.trees = HET_TREES,
      honesty = TRUE,
      ci.group.size = 2,
      tune.parameters = HET_TUNE_PARAMETERS,
      num.threads = NUM_THREADS,
      seed = RANDOM_SEED
    ),
    error = function(e) e
  )

  if (inherits(forest, "error")) {
    return(empty_binary_result(
      pair_name, model_type, contrast, arm_label,
      "fit_failed", forest$message, support
    ))
  }

  prefix <- function(dat) {
    if (nrow(dat) == 0) return(dat)
    dat |>
      mutate(
        analysis_label = ANALYSIS_LABEL,
        pair = pair_name,
        model_type = model_type,
        contrast = contrast,
        arm_label = arm_label,
        .before = 1
      )
  }

  overall_ate <- tryCatch(
    prefix(extract_ate(average_treatment_effect(forest))),
    error = function(e) prefix(tibble(error_message = e$message))
  )

  calibration <- tryCatch(
    prefix(tidy_inference(test_calibration(forest), "test_calibration")),
    error = function(e) prefix(tibble(error_message = e$message))
  )

  modifier_matrix <- cbind(
    male_sex_1 = as.numeric(X_b[, "sex"] == 1),
    age_ge_60 = as.numeric(X_b[, "age"] >= 60)
  )
  variable_columns <- apply(modifier_matrix, 2, function(x) length(unique(x[is.finite(x)])) > 1)
  modifier_matrix <- modifier_matrix[, variable_columns, drop = FALSE]
  blp <- if (ncol(modifier_matrix) == 0) {
    prefix(tibble(error_message = "effect modifiers had no variation"))
  } else {
    tryCatch(
      prefix(tidy_inference(
        best_linear_projection(forest, A = modifier_matrix),
        "best_linear_projection"
      )),
      error = function(e) prefix(tibble(error_message = e$message))
    )
  }

  tau_oob <- tryCatch(drop(predict(forest)$predictions), error = function(e) rep(NA_real_, nrow(X_b)))
  priorities <- switch(
    model_type,
    improvement = -tau_oob,
    deterioration = tau_oob,
    mixed = abs(tau_oob)
  )
  rate <- tryCatch(
    prefix(tidy_inference(
      rank_average_treatment_effect(
        forest,
        priorities = priorities,
        target = "AUTOC"
      ),
      "rank_average_treatment_effect"
    )),
    error = function(e) prefix(tibble(error_message = e$message))
  )

  subgroup_ate <- estimate_subgroup_ates(
    forest, X_b, pair_name, model_type, contrast, arm_label
  )
  subgroup_support <- binary_subgroup_support(
    X_b, Y_b, W_b, forest$W.hat,
    pair_name, model_type, contrast, arm_label
  )

  cate_summary <- prefix(tibble(
    n = sum(is.finite(tau_oob)),
    mean_cate = mean(tau_oob, na.rm = TRUE),
    sd_cate = sd(tau_oob, na.rm = TRUE),
    min_cate = min(tau_oob, na.rm = TRUE),
    max_cate = max(tau_oob, na.rm = TRUE),
    inference = "descriptive out-of-bag CATE distribution"
  ))

  cate_quartiles <- tibble()
  if (sum(is.finite(tau_oob)) >= 4) {
    quartile <- rep(NA_integer_, length(tau_oob))
    quartile[is.finite(tau_oob)] <- ntile(tau_oob[is.finite(tau_oob)], 4)
    cate_quartiles <- bind_rows(lapply(1:4, function(q) {
      idx <- which(quartile == q)
      tryCatch({
        est <- extract_single_ate(average_treatment_effect(forest, subset = idx))
        prefix(tibble(
          quartile = q,
          n = length(idx),
          estimate = est[["estimate"]],
          std_error = est[["std_error"]],
          low = estimate - 1.96 * std_error,
          high = estimate + 1.96 * std_error,
          inference = "exploratory; quartiles selected using estimated OOB CATE"
        ))
      }, error = function(e) prefix(tibble(
        quartile = q,
        n = length(idx),
        error_message = e$message
      )))
    }))
  }

  vi <- tryCatch(
    prefix(tibble(
      variable = colnames(X_b),
      importance = as.numeric(variable_importance(forest)),
      inference = "splitting importance; not a formal effect-modifier test"
    ) |>
      arrange(desc(importance))),
    error = function(e) prefix(tibble(error_message = e$message))
  )

  bootstrap_difference <- bootstrap_subgroup_differences(
    X_b, Y_b, W_b, pair_name, model_type, contrast, arm_label
  )

  stability <- bind_rows(lapply(STABILITY_SEEDS, function(stability_seed) {
    stability_forest <- if (stability_seed == RANDOM_SEED) {
      forest
    } else {
      tryCatch(
        causal_forest(
          X = X_b,
          Y = Y_b,
          W = W_b,
          num.trees = HET_TREES,
          honesty = TRUE,
          ci.group.size = 2,
          num.threads = NUM_THREADS,
          seed = stability_seed
        ),
        error = function(e) NULL
      )
    }
    if (is.null(stability_forest)) {
      return(prefix(tibble(seed = stability_seed, status = "fit_failed")))
    }
    stability_tau <- drop(predict(stability_forest)$predictions)
    stability_ate <- tryCatch(
      extract_single_ate(average_treatment_effect(stability_forest))[["estimate"]],
      error = function(e) NA_real_
    )
    stability_calibration <- tryCatch(
      tidy_inference(test_calibration(stability_forest), "test_calibration"),
      error = function(e) tibble()
    )
    differential <- if (
      nrow(stability_calibration) > 0 &&
      all(c("term", "estimate", "p_value") %in% names(stability_calibration))
    ) {
      stability_calibration |>
        filter(grepl("differential", term, ignore.case = TRUE))
    } else {
      tibble()
    }
    prefix(tibble(
      seed = stability_seed,
      status = "fitted",
      ate = stability_ate,
      cate_sd = sd(stability_tau, na.rm = TRUE),
      differential_calibration_estimate = if (nrow(differential) > 0) differential$estimate[1] else NA_real_,
      differential_calibration_p = if (nrow(differential) > 0) differential$p_value[1] else NA_real_
    ))
  }))

  list(
    status = prefix(tibble(status = "fitted", reason = NA_character_)),
    support = support,
    overall_ate = overall_ate,
    calibration = calibration,
    blp = blp,
    rate = rate,
    subgroup_ate = subgroup_ate,
    subgroup_support = subgroup_support,
    cate_summary = cate_summary,
    cate_quartiles = cate_quartiles,
    variable_importance = vi,
    bootstrap_difference = bootstrap_difference,
    stability = stability,
    forest = forest
  )
}

# =============================================================================
# 7. Run models
# =============================================================================

pair_list <- combn(names(lifestyle_vars), 2, simplify = FALSE)
pair_names <- vapply(pair_list, function(x) paste(x, collapse = "_"), character(1))

selected_heterogeneity_pairs <- if (tolower(HET_PAIRS_RAW) == "all") {
  pair_names
} else {
  requested <- split_csv(HET_PAIRS_RAW)
  unknown <- setdiff(requested, pair_names)
  if (length(unknown) > 0) {
    stop("Unknown HET_PAIRS: ", paste(unknown, collapse = ", "))
  }
  requested
}

active_arm_filter <- if (tolower(HET_ACTIVE_ARMS_RAW) == "all") {
  NULL
} else {
  suppressWarnings(as.integer(split_csv(HET_ACTIVE_ARMS_RAW)))
}

multiarm_results <- list()
if (RUN_MULTIARM) {
  for (pair in pair_list) {
    for (model_type in valid_model_types) {
      result_id <- paste(pair[1], pair[2], model_type, sep = "_")
      log_msg("Multi-arm: ", result_id)
      multiarm_results[[result_id]] <- fit_multiarm_model(
        pair[1], pair[2], model_type
      )
    }
  }
}

heterogeneity_results <- list()
if (RUN_HETEROGENEITY) {
  for (pair in pair_list) {
    pair_name <- paste(pair, collapse = "_")
    if (!pair_name %in% selected_heterogeneity_pairs) next

    for (model_type in HET_MODEL_TYPES) {
      W_raw <- make_pair_treatment(df, pair[1], pair[2], model_type)
      observed_active <- setdiff(sort(unique(W_raw[!is.na(W_raw)])), 0)
      if (!is.null(active_arm_filter)) {
        observed_active <- intersect(observed_active, active_arm_filter)
      }

      for (active_arm in observed_active) {
        result_id <- paste(pair_name, model_type, paste0("arm", active_arm), sep = "_")
        log_msg("Binary heterogeneity: ", result_id)
        heterogeneity_results[[result_id]] <- fit_binary_contrast(
          pair[1], pair[2], model_type, active_arm
        )
      }
    }
  }
}

# =============================================================================
# 8. Consolidate, correct multiplicity and save
# =============================================================================

multi_status <- safe_bind(multiarm_results, "status")
multi_arms <- safe_bind(multiarm_results, "arm_table")
multi_overlap <- safe_bind(multiarm_results, "overlap")
multi_ate <- safe_bind(multiarm_results, "ate")
multi_cate <- safe_bind(multiarm_results, "cate_summary")
multi_vi <- safe_bind(multiarm_results, "variable_importance")
multi_interaction <- safe_bind(multiarm_results, "interaction")
if ("p" %in% names(multi_ate)) multi_ate <- add_fdr(multi_ate, "p", "p_fdr")

het_status <- safe_bind(heterogeneity_results, "status")
het_support <- safe_bind(heterogeneity_results, "support")
het_ate <- safe_bind(heterogeneity_results, "overall_ate")
het_calibration <- safe_bind(heterogeneity_results, "calibration")
het_blp <- safe_bind(heterogeneity_results, "blp")
het_rate <- safe_bind(heterogeneity_results, "rate")
het_subgroup_ate <- safe_bind(heterogeneity_results, "subgroup_ate")
het_subgroup_support <- safe_bind(heterogeneity_results, "subgroup_support")
het_cate <- safe_bind(heterogeneity_results, "cate_summary")
het_quartiles <- safe_bind(heterogeneity_results, "cate_quartiles")
het_vi <- safe_bind(heterogeneity_results, "variable_importance")
het_bootstrap <- safe_bind(heterogeneity_results, "bootstrap_difference")
het_stability <- safe_bind(heterogeneity_results, "stability")

if ("p" %in% names(het_ate)) het_ate <- add_fdr(het_ate, "p", "p_fdr")

if (
  nrow(het_calibration) > 0 &&
  all(c("term", "p_value") %in% names(het_calibration))
) {
  differential_idx <- grepl("differential", het_calibration$term, ignore.case = TRUE)
  het_calibration$p_fdr <- NA_real_
  het_calibration$p_fdr[differential_idx] <- p.adjust(
    het_calibration$p_value[differential_idx], method = "BH"
  )
}

if (
  nrow(het_blp) > 0 &&
  all(c("term", "p_value") %in% names(het_blp))
) {
  modifier_idx <- !grepl("intercept", het_blp$term, ignore.case = TRUE)
  het_blp$p_fdr <- NA_real_
  het_blp$p_fdr[modifier_idx] <- p.adjust(
    het_blp$p_value[modifier_idx], method = "BH"
  )
}

if (nrow(het_rate) > 0 && "p_value" %in% names(het_rate)) {
  het_rate <- add_fdr(het_rate, "p_value", "p_fdr")
}
if (nrow(het_bootstrap) > 0 && "p_value" %in% names(het_bootstrap)) {
  het_bootstrap <- add_fdr(het_bootstrap, "p_value", "p_fdr")
}

safe_write(multi_status, "multiarm_model_status.csv")
safe_write(multi_arms, "multiarm_arm_support.csv")
safe_write(multi_overlap, "multiarm_overlap.csv")
safe_write(multi_ate, "multiarm_ate.csv")
safe_write(multi_cate, "multiarm_cate_summary_descriptive.csv")
safe_write(multi_vi, "multiarm_variable_importance_descriptive.csv")
safe_write(multi_interaction, "multiarm_interaction_point_estimates.csv")

safe_write(het_status, "heterogeneity_model_status.csv")
safe_write(het_support, "heterogeneity_contrast_support.csv")
safe_write(het_ate, "heterogeneity_binary_ate.csv")
safe_write(het_calibration, "heterogeneity_test_calibration.csv")
safe_write(het_blp, "heterogeneity_best_linear_projection.csv")
safe_write(het_rate, "heterogeneity_rate_autoc.csv")
safe_write(het_subgroup_ate, "heterogeneity_subgroup_ate.csv")
safe_write(het_subgroup_support, "heterogeneity_subgroup_support.csv")
safe_write(het_cate, "heterogeneity_cate_summary_descriptive.csv")
safe_write(het_quartiles, "heterogeneity_cate_quartiles_exploratory.csv")
safe_write(het_vi, "heterogeneity_variable_importance_descriptive.csv")
safe_write(het_bootstrap, "heterogeneity_bootstrap_subgroup_differences.csv")
safe_write(het_stability, "heterogeneity_seed_stability.csv")

if (SAVE_FORESTS) {
  saveRDS(
    list(multiarm = multiarm_results, heterogeneity = heterogeneity_results),
    file.path(OUTPUT_DIR, "all_fitted_forests.rds")
  )
}

capture.output(sessionInfo(), file = file.path(OUTPUT_DIR, "session_info.txt"))

run_config <- tibble(
  setting = c(
    "DATA_PATH", "OUTCOME_NAME", "ANALYSIS_LABEL", "RANDOM_SEED",
    "NUM_THREADS", "MULTIARM_TREES", "HET_TREES", "BOOT_TREES",
    "N_BOOT", "RUN_MULTIARM", "RUN_HETEROGENEITY",
    "ENFORCE_STRUCTURAL_FEASIBILITY", "HET_PAIRS", "HET_MODEL_TYPES",
    "HET_ACTIVE_ARMS", "STABILITY_SEEDS", "HET_TUNE_PARAMETERS",
    "JOBLIB_PATH", "JOBLIB_BRIDGE", "JOBLIB_COMPRESS", "PYTHON_BIN"
  ),
  value = c(
    DATA_PATH, OUTCOME_NAME, ANALYSIS_LABEL, RANDOM_SEED,
    NUM_THREADS, MULTIARM_TREES, HET_TREES, BOOT_TREES,
    N_BOOT, RUN_MULTIARM, RUN_HETEROGENEITY,
    ENFORCE_STRUCTURAL_FEASIBILITY, HET_PAIRS_RAW,
    paste(HET_MODEL_TYPES, collapse = ","), HET_ACTIVE_ARMS_RAW,
    paste(STABILITY_SEEDS, collapse = ","), HET_TUNE_PARAMETERS,
    JOBLIB_PATH, JOBLIB_BRIDGE, JOBLIB_COMPRESS, PYTHON_BIN
  )
)
safe_write(run_config, "run_configuration.csv")


# App-facing joblib payload
# The table the web app reads: the multi-arm ATE table exactly as
# multiarm_ate.csv had it. estimate/low/high are risk differences as fractions;
# the app multiplies by 100 for percentage points.
build_webapp_ate <- function(dat) {
  required <- c("pair", "model_type", "arm_label", "estimate", "low", "high")
  if (is.null(dat) || nrow(dat) == 0 || !all(required %in% names(dat))) {
    return(tibble())
  }
  # Failed analyses carry missing estimates. Keep only complete effect rows and
  # export the stable schema consumed by cvd_webapp/model_io.py.
  dat |>
    filter(is.finite(estimate), is.finite(low), is.finite(high)) |>
    select(any_of(c(
      "pair", "model_type", "arm", "arm_label",
      "estimate", "low", "high", "p", "p_fdr"
    ))) |>
    arrange(pair, factor(model_type, levels = valid_model_types), arm)
}

validate_webapp_ate <- function(dat) {
  required <- c("pair", "model_type", "arm_label", "estimate", "low", "high")
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0) {
    stop("App-facing ATE table is missing columns: ", paste(missing, collapse = ", "))
  }
  if (nrow(dat) == 0) {
    stop("No estimable multi-arm effects were available for the app joblib.")
  }
  if (any(is.na(dat$pair) | !nzchar(trimws(dat$pair)))) {
    stop("App-facing ATE table contains a missing or blank pair.")
  }
  if (any(is.na(dat$arm_label) | !nzchar(trimws(dat$arm_label)))) {
    stop("App-facing ATE table contains a missing or blank arm_label.")
  }
  unknown_types <- setdiff(unique(dat$model_type), valid_model_types)
  if (length(unknown_types) > 0) {
    stop("App-facing ATE table contains unknown model_type values: ",
         paste(unknown_types, collapse = ", "))
  }
  numeric_columns <- c("estimate", "low", "high")
  if (any(!vapply(dat[numeric_columns], is.numeric, logical(1)))) {
    stop("estimate, low, and high must be numeric.")
  }
  if (any(!is.finite(as.matrix(dat[numeric_columns])))) {
    stop("App-facing ATE table contains non-finite estimates or intervals.")
  }
  if (any(dat$low > dat$estimate | dat$estimate > dat$high)) {
    stop("At least one estimate lies outside its confidence interval.")
  }
  duplicate_key <- duplicated(dat[c("pair", "model_type", "arm_label")])
  if (any(duplicate_key)) {
    stop("Duplicate pair/model_type/arm_label rows found in app-facing ATE table.")
  }
  if (!"improvement" %in% dat$model_type) {
    stop("No improvement scenarios were exported; the current app would be empty.")
  }
  dat
}

app_ate <- validate_webapp_ate(build_webapp_ate(multi_ate))
log_msg(
  "App payload: ", nrow(app_ate), " rows; ", ncol(app_ate),
  " columns; ", n_distinct(app_ate$pair), " lifestyle pairs; scenarios=",
  paste(sort(unique(app_ate$model_type)), collapse = ",")
)
write_joblib(app_ate, JOBLIB_PATH)
log_msg("Completed.")
log_msg("Multi-arm models attempted: ", length(multiarm_results))
log_msg("Binary heterogeneity contrasts attempted: ", length(heterogeneity_results))
log_msg("CSV results written to: ", OUTPUT_DIR)
log_msg("Joblib written to: ", JOBLIB_PATH)
