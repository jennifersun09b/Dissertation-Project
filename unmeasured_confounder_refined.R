# Run with: Rscript unmeasured_confounder_refined.R

# =============================================================================
# Unmeasured-confounding sensitivity analysis for the single-variable
# causal-forest models
# =============================================================================
# Plan, in plain words:
#   1. Load the data (the SAME cohort file the primary model used).
#   2. Refit the SAME binary causal forest, with the SAME seed, tree count and
#      tuning settings as single_variable3.R, so the sensitivity analysis is
#      attached to the numbers that appear in the primary results table.
#   3. Put the effect on a risk scale, using doubly robust arm-specific risks
#      instead of a clipped average of the control-arm prediction.
#   4. Run three unmeasured-confounding analyses on that risk scale:
#        a) E-values (VanderWeele & Ding), for the point estimate and for the
#           confidence limit closest to the null.
#        b) Quantitative bias analysis (QBA) for a generic binary confounder U,
#           over the prespecified scenarios AND over a factorial grid.
#        c) Tipping-point analysis: the exact strength of U that would move the
#           estimate to zero, solved in closed form.
#   5. Write every result to CSV, with a status row when something cannot be
#      computed, so a missing number is always visible as a statement.
#
# What changed relative to unmeasured_confounder_running.R, and why
# -----------------------------------------------------------------------------
#   1. DATA_PATH, OUTCOME_NAME, ANALYSIS_LABEL and OUTPUT_DIR are now read from
#      the environment, exactly as in single_variable3.R. The old script
#      hardcoded the sensitive-cohort path, so the primary cohort was never
#      given an E-value. The default here is the PRIMARY cohort; the sensitive
#      run is one environment variable away (example below).
#
#   2. The contrast registry now carries health_direction and
#      interpretation_flag, matching single_variable3.R. The old `direction`
#      column labelled mental-health and alcohol increases as "increase", which
#      is a deterioration, so the sensitivity table could not be joined to the
#      primary table without re-reading the labels by hand.
#
#   3. The forest is now fitted with honesty = TRUE and ci.group.size = 2, the
#      same as the primary script. The old script left these at their defaults.
#
#   4. E-values are computed for ALL FIVE estimands (full sample, overlap
#      weighted, and the three propensity trims), not only for the 0.05-0.95
#      trimmed subset. The primary estimand in single_variable3.R is the FULL
#      eligible sample, so the old script produced E-values for an estimand that
#      is not the headline one.
#
#   5. The reference risk is now the doubly robust AIPW mean of Y(0), not the
#      mean of a clipped mu0 prediction. Because the AIPW arm risks satisfy
#      risk_moved - risk_stayed = the grf ATE exactly, the risk scale and the
#      risk-difference scale can no longer disagree. The identity is checked and
#      the discrepancy is written out as ate_consistency_diff.
#
#   6. The risk-ratio confidence interval is now a real interval. The old file
#      mapped the ATE limits while holding the reference risk fixed and said so
#      in a note. Here the log-RR standard error comes from the influence
#      function of both arm risks, and a participant bootstrap of the same
#      scores is reported alongside it. Nine of the twenty-four transitions
#      previously returned no E-value at all because a mapped limit fell outside
#      (0, 1); that failure mode is gone, and a crude (unadjusted) E-value is
#      also reported as a floor so no row is ever silently blank.
#
#   7. QBA now reports uncertainty. Every prespecified scenario is recomputed on
#      each bootstrap replicate, so the adjusted risk difference has a 95%
#      interval and it is possible to say whether the adjusted estimate still
#      excludes zero. The old output was a point estimate only.
#
#   8. QBA is extended with a factorial grid and with a closed-form tipping
#      point: given a prevalence pattern for U, the exact RR_UY that drives the
#      adjusted risk difference to zero. Seven prespecified scenarios cannot
#      show where the boundary is; the tipping point states it directly.
#
#   9. Benjamini-Hochberg FDR is applied across the 24-transition family, so the
#      E-value table can be read next to the primary q-values. An E-value only
#      carries weight for a transition that survived FDR, and the table now says
#      which those are.
#
#  10. Heterogeneity output (variable importance, calibration, best linear
#      projection, RATE, CATE quartiles, subgroup ATEs) was computed and printed
#      but never saved. It is now written to CSV, using the same unclass() fix
#      as single_variable3.R, without which the calibration, BLP and RATE
#      objects coerce to nothing.
#
#  11. Subgroup E-values are added, for the same seven subgroups the primary
#      script reports, so a subgroup finding can be judged on the same scale as
#      the overall finding.
#
#  12. safe_write(), run_configuration.csv and session_info.txt are added, so
#      every file carries at least one interpretable row and the run is
#      reproducible from its own output.
#
# Example primary run (the default)
#   DATA_PATH=/home/rmhiund/causal_analysis/Cohort/primary_single_variable.csv \
#   OUTCOME_NAME=CVD_outcome ANALYSIS_LABEL=primary \
#   OUTPUT_DIR=/path/unmeasured_confounder_primary Rscript unmeasured_confounder_refined.R
#
# Example sensitive run
#   DATA_PATH=/home/rmhiund/causal_analysis/Cohort/sensitive_single_variable.csv \
#   ANALYSIS_LABEL=sensitive \
#   OUTPUT_DIR=/path/unmeasured_confounder_sensitive Rscript unmeasured_confounder_refined.R
#
# Example fast check (no bootstrap, one domain)
#   CONTRAST_DOMAINS=mental N_BOOT=0 NUM_TREES=1000 Rscript unmeasured_confounder_refined.R
#
# Each comparison is binary:
#     stayed (W = 0) = people who kept their baseline value
#     moved  (W = 1) = people who moved to a given destination value
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
TEMPORAL_DESIGN <- env_chr(
  "TEMPORAL_DESIGN",
  "broad T0-to-T2 transition analysis; temporal interpretation is exploratory"
)
OUTPUT_DIR <- env_chr(
  "OUTPUT_DIR",
  file.path(getwd(), paste0("unmeasured_confounder_results_", ANALYSIS_LABEL))
)

# These four must match single_variable3.R for the sensitivity analysis to refer
# to the same fitted model as the primary result.
RANDOM_SEED <- env_int("RANDOM_SEED", 42)
NUM_TREES <- env_int("NUM_TREES", 5000)
NUM_THREADS <- env_int("NUM_THREADS", 1)
TUNE_PARAMETERS <- env_chr("TUNE_PARAMETERS", "all")

N_BOOT <- env_int("N_BOOT", 1000)
# The AIPW arm-risk scores divide by the estimated propensity. On the untrimmed
# full sample a propensity very close to 0 or 1 can produce an enormous score.
# The floor bounds that division only; it does not change which people are in
# the estimand, and it is inert on any trimmed subset. grf itself does not floor,
# so whenever the floor actually binds the arm risks and the grf ATE separate by
# construction: propensity_summary.csv reports how many people were affected and
# ate_estimands_with_risk_scale.csv reports the resulting discrepancy. The
# default is deliberately small so that the identity holds on almost every row;
# the trimmed estimands are the ones to quote when it does not.
PROPENSITY_FLOOR <- env_num("PROPENSITY_FLOOR", 0.001)
PRIMARY_ESTIMAND <- env_chr("PRIMARY_ESTIMAND", paste0("full_sample_", ANALYSIS_LABEL))
RUN_HETEROGENEITY <- env_bool("RUN_HETEROGENEITY", TRUE)
RUN_SUBGROUP_EVALUES <- env_bool("RUN_SUBGROUP_EVALUES", TRUE)
SAVE_FORESTS <- env_bool("SAVE_FORESTS", FALSE)
AGE_CUT <- env_num("AGE_CUT", 60)

CONTRAST_DOMAINS_RAW <- env_chr("CONTRAST_DOMAINS", "all")
EXTRA_COVARIATES <- split_csv(env_chr("EXTRA_COVARIATES", ""))

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
  flush.console()
}

log_msg("Analysis: ", ANALYSIS_LABEL)
log_msg("Data: ", DATA_PATH)
log_msg("Outcome: ", OUTCOME_NAME)
log_msg("Output: ", OUTPUT_DIR)
log_msg("Trees: ", NUM_TREES, "; tuning: ", TUNE_PARAMETERS, "; seed: ", RANDOM_SEED)
log_msg("Primary estimand for the E-value headline: ", PRIMARY_ESTIMAND)
log_msg(
  "Bootstrap replicates: ", N_BOOT,
  if (N_BOOT <= 0) " (percentile intervals and QBA intervals disabled)" else
    " (participant resampling of the doubly robust arm-risk scores)"
)
log_msg("Propensity floor used inside the AIPW arm-risk scores: ", PROPENSITY_FLOOR)

# =============================================================================
# 2. Transition registry and covariates
# =============================================================================

# health_direction describes the CLINICAL direction of each transition. It is
# not the direction of the raw score: for mental health and alcohol a higher
# score is worse, so an increase is a deterioration. interpretation_flag records
# the score convention. These two columns are copied from single_variable3.R so
# that this file and the primary results table use identical labels.
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
  covariate_names, OUTCOME_NAME, contrasts$t0_col, contrasts$t1_col
))
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

non_numeric <- covariate_names[!vapply(df[covariate_names], is.numeric, logical(1))]
if (length(non_numeric) > 0) {
  stop("Encode these covariates numerically before GRF: ", paste(non_numeric, collapse = ", "))
}

X_all <- df |> select(all_of(covariate_names)) |> as.matrix()
Y_all <- as.numeric(df[[OUTCOME_NAME]])

if (!all(na.omit(unique(Y_all)) %in% c(0, 1))) {
  stop(OUTCOME_NAME, " must be coded 0/1.")
}

cohort_summary <- tibble(
  analysis_label = ANALYSIS_LABEL,
  temporal_design = TEMPORAL_DESIGN,
  data_path = DATA_PATH,
  n_loaded = nrow(df),
  events = sum(Y_all == 1, na.rm = TRUE),
  non_events = sum(Y_all == 0, na.rm = TRUE),
  missing_outcome = sum(is.na(Y_all)),
  complete_covariates = sum(complete.cases(X_all)),
  event_rate = mean(Y_all, na.rm = TRUE)
)

log_msg(
  "Cohort loaded: n = ", cohort_summary$n_loaded,
  "; events = ", cohort_summary$events,
  "; event rate = ", round(cohort_summary$event_rate, 5)
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

# grf and lmtest return objects carrying explicit S3 class attributes:
# test_calibration() and best_linear_projection() return "coeftest" matrices and
# rank_average_treatment_effect() returns a classed list. as.data.frame()
# dispatches on class, so those attributes stop it reaching as.data.frame.matrix()
# or as.data.frame.list(); it falls through to as.data.frame.default(), which
# errors. unclass() first is what makes these objects tidyable. Without it these
# tables coerce to nothing and the CSVs are written with a header and no rows.
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

  dat <- tryCatch(to_inference_frame(object), error = function(e) NULL)
  if (is.null(dat) || nrow(dat) == 0) {
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
    terms <- if (!is.na(target_col)) as.character(dat[[target_col]]) else paste0("row_", seq_len(nrow(dat)))
  }

  n_rows <- nrow(dat)
  estimate <- if (!is.na(estimate_col)) as.numeric(dat[[estimate_col]]) else rep(NA_real_, n_rows)
  std_error <- if (!is.na(se_col)) as.numeric(dat[[se_col]]) else rep(NA_real_, n_rows)
  statistic <- if (!is.na(stat_col)) as.numeric(dat[[stat_col]]) else rep(NA_real_, n_rows)
  p_value <- if (!is.na(p_col)) as.numeric(dat[[p_col]]) else rep(NA_real_, n_rows)

  # RATE/AUTOC reports only an estimate and a standard error, so derive the test
  # statistic and two-sided p-value when the object does not carry them.
  needs_stat <- !is.finite(statistic) & is.finite(estimate) & is.finite(std_error) & std_error > 0
  statistic[needs_stat] <- estimate[needs_stat] / std_error[needs_stat]
  needs_p <- !is.finite(p_value) & is.finite(statistic)
  p_value[needs_p] <- 2 * pnorm(abs(statistic[needs_p]), lower.tail = FALSE)

  tibble(
    source = source, term = terms, estimate = estimate,
    std_error = std_error, statistic = statistic, p_value = p_value
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
        estimate = estimate, std_error = std_error,
        low = estimate - 1.96 * std_error, high = estimate + 1.96 * std_error,
        z = z, p_value = 2 * pnorm(abs(z), lower.tail = FALSE)
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
    estimate = estimate, std_error = std_error,
    low = estimate - 1.96 * std_error, high = estimate + 1.96 * std_error,
    z = z, p_value = 2 * pnorm(abs(z), lower.tail = FALSE)
  )
}

prefix_result <- function(dat, spec) {
  if (is.null(dat) || nrow(dat) == 0) return(tibble())
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

# Every CSV carries at least one row. When a table is empty, an explicit
# status/note row is substituted, so a missing result is visible as a statement
# rather than as a file containing only column names.
safe_write <- function(dat, filename, note = "no rows produced by this run") {
  if (is.null(dat) || ncol(dat) == 0 || nrow(dat) == 0) {
    dat <- tibble(analysis_label = ANALYSIS_LABEL, status = "no_rows", note = note)
  }
  readr::write_csv(dat, file.path(OUTPUT_DIR, filename))
}

safe_bind <- function(results, field) {
  if (length(results) == 0) return(tibble())
  bind_rows(lapply(results, function(x) x[[field]]))
}

add_fdr_subset <- function(dat, subset_idx, p_col = "p_value", output_col = "p_fdr") {
  if (nrow(dat) == 0 || !p_col %in% names(dat)) return(dat)
  dat[[output_col]] <- NA_real_
  valid <- subset_idx & is.finite(dat[[p_col]])
  dat[[output_col]][valid] <- p.adjust(dat[[p_col]][valid], method = "BH")
  dat
}

make_subgroups <- function(X) {
  list(
    female_sex_0 = which(X[, "sex"] == 0),
    male_sex_1 = which(X[, "sex"] == 1),
    age_lt_60 = which(X[, "age"] < AGE_CUT),
    age_ge_60 = which(X[, "age"] >= AGE_CUT),
    bmi_lt_25 = which(X[, "bmi"] < 25),
    bmi_25_to_lt_30 = which(X[, "bmi"] >= 25 & X[, "bmi"] < 30),
    bmi_ge_30 = which(X[, "bmi"] >= 30)
  )
}

# =============================================================================
# 5. E-value helpers
# =============================================================================
# These helpers do not assume the unknown confounder is any particular variable.
# U can represent one omitted factor or the combined effect of several omitted
# factors after adjustment for the measured covariates.

# Standard E-value for a risk ratio. Protective RRs are inverted so the reported
# E-value is always >= 1 and has the usual interpretation.
evalue_from_rr <- function(rr) {
  rr <- as.numeric(rr)
  if (length(rr) != 1 || !is.finite(rr) || rr <= 0) return(NA_real_)
  rr_strength <- if (rr < 1) 1 / rr else rr
  rr_strength + sqrt(rr_strength * (rr_strength - 1))
}

# E-value for the confidence limit closest to the null value RR = 1.
# If the interval already includes 1, its E-value is exactly 1.
evalue_from_ci <- function(rr, low, high) {
  values <- as.numeric(c(rr, low, high))
  if (any(!is.finite(values)) || any(values <= 0) || low > high) return(NA_real_)
  if (low <= 1 && high >= 1) return(1)
  closest_limit <- if (high < 1) 1 / high else low
  evalue_from_rr(closest_limit)
}

evalue_ci_note <- function(rr, low, high) {
  if (any(!is.finite(c(rr, low, high)))) return("CI E-value not calculated")
  if (low <= 1 && high >= 1) {
    "RR interval includes 1, so the CI E-value is 1 by definition"
  } else {
    "CI E-value taken at the confidence limit closest to RR = 1"
  }
}

# =============================================================================
# 6. Doubly robust arm risks
# =============================================================================
# The forest's nuisance estimates satisfy E[Y|X] = mu0(X) + e(X) * tau(X), so
#   mu0(X) = Y.hat - e(X) * tau(X)      mu1(X) = Y.hat + (1 - e(X)) * tau(X).
# The AIPW score for each potential outcome is the plug-in prediction plus an
# inverse-probability-weighted residual correction:
#   score1_i = mu1_i + (W_i / e_i) * (Y_i - mu1_i)
#   score0_i = mu0_i + ((1 - W_i) / (1 - e_i)) * (Y_i - mu0_i)
# Their difference is algebraically identical to the doubly robust score grf
# uses inside average_treatment_effect(), so
#   mean(score1) - mean(score0) = the grf ATE
# holds exactly on any subset and under any weighting. That identity is what
# makes the risk-ratio scale and the risk-difference scale mutually consistent,
# and it is checked and reported as ate_consistency_diff.
#
# The previous approach averaged a clipped mu0 prediction and then added the
# ATE. Clipping breaks the identity, the reference risk is not doubly robust,
# and a mapped confidence limit can leave (0, 1), which is why nine of the
# twenty-four transitions previously produced no E-value.
arm_risk_scores <- function(forest, Y, W) {
  propensity <- as.numeric(forest$W.hat)
  y_hat <- as.numeric(forest$Y.hat)
  tau <- as.numeric(predict(forest)$predictions)

  e_bounded <- pmin(pmax(propensity, PROPENSITY_FLOOR), 1 - PROPENSITY_FLOOR)
  mu1 <- y_hat + (1 - propensity) * tau
  mu0 <- y_hat - propensity * tau

  list(
    propensity = propensity,
    tau = tau,
    mu0 = mu0,
    mu1 = mu1,
    score0 = mu0 + ((1 - W) / (1 - e_bounded)) * (Y - mu0),
    score1 = mu1 + (W / e_bounded) * (Y - mu1),
    n_floored = sum(propensity < PROPENSITY_FLOOR | propensity > 1 - PROPENSITY_FLOOR)
  )
}

weighted_mean_w <- function(x, w) sum(w * x) / sum(w)

# Sandwich standard error of a weighted mean, treating the weights as fixed.
weighted_mean_se <- function(x, w, mu) sqrt(sum(w^2 * (x - mu)^2)) / sum(w)

# Point estimates and delta-method inference for both arm risks, the risk
# difference and the log risk ratio, from the AIPW scores.
#
# The log-RR standard error comes from the influence function
#   psi_i = (score1_i - r1) / r1 - (score0_i - r0) / r0,
# which propagates the uncertainty in BOTH arms and their covariance. The
# previous file held the reference risk fixed, which understates the interval.
risk_scale_from_scores <- function(scores, idx, weights = NULL) {
  s0 <- scores$score0[idx]
  s1 <- scores$score1[idx]
  w <- if (is.null(weights)) rep(1, length(idx)) else weights[idx]

  ok <- is.finite(s0) & is.finite(s1) & is.finite(w) & w > 0
  if (sum(ok) < 2) {
    return(tibble(
      n_used = sum(ok),
      risk_stayed = NA_real_, risk_stayed_se = NA_real_,
      risk_moved = NA_real_, risk_moved_se = NA_real_,
      rd_score = NA_real_, rd_score_se = NA_real_,
      rr = NA_real_, log_rr_se = NA_real_,
      rr_low = NA_real_, rr_high = NA_real_,
      risk_scale_status = "not calculated: fewer than two usable scores"
    ))
  }

  s0 <- s0[ok]; s1 <- s1[ok]; w <- w[ok]
  r0 <- weighted_mean_w(s0, w)
  r1 <- weighted_mean_w(s1, w)
  d <- s1 - s0
  rd <- weighted_mean_w(d, w)

  r0_se <- weighted_mean_se(s0, w, r0)
  r1_se <- weighted_mean_se(s1, w, r1)
  rd_se <- weighted_mean_se(d, w, rd)

  if (!is.finite(r0) || r0 <= 0 || !is.finite(r1) || r1 <= 0) {
    return(tibble(
      n_used = length(w),
      risk_stayed = r0, risk_stayed_se = r0_se,
      risk_moved = r1, risk_moved_se = r1_se,
      rd_score = rd, rd_score_se = rd_se,
      rr = NA_real_, log_rr_se = NA_real_,
      rr_low = NA_real_, rr_high = NA_real_,
      risk_scale_status = paste(
        "RR not calculated: a doubly robust arm risk was not strictly positive",
        "(this happens when a sparse arm has almost no events)"
      )
    ))
  }

  psi <- (s1 - r1) / r1 - (s0 - r0) / r0
  log_rr <- log(r1) - log(r0)
  log_rr_se <- weighted_mean_se(psi, w, 0)

  tibble(
    n_used = length(w),
    risk_stayed = r0, risk_stayed_se = r0_se,
    risk_moved = r1, risk_moved_se = r1_se,
    rd_score = rd, rd_score_se = rd_se,
    rr = exp(log_rr), log_rr_se = log_rr_se,
    rr_low = exp(log_rr - 1.96 * log_rr_se),
    rr_high = exp(log_rr + 1.96 * log_rr_se),
    risk_scale_status = "estimated from doubly robust arm-specific AIPW scores"
  )
}

# Participant bootstrap of the same scores. No forest is refitted: the AIPW
# arm risks are means of fixed per-person scores, so resampling people is exact
# for those quantities conditional on the fitted nuisance functions. This is the
# same argument single_variable3.R uses for its subgroup bootstrap.
bootstrap_arm_risks <- function(scores, idx, weights = NULL, seed_offset = 0) {
  if (N_BOOT <= 0) return(NULL)
  s0 <- scores$score0[idx]
  s1 <- scores$score1[idx]
  w <- if (is.null(weights)) rep(1, length(idx)) else weights[idx]
  ok <- is.finite(s0) & is.finite(s1) & is.finite(w) & w > 0
  if (sum(ok) < 10) return(NULL)
  s0 <- s0[ok]; s1 <- s1[ok]; w <- w[ok]
  n <- length(s0)

  set.seed(RANDOM_SEED + 9000 + seed_offset)
  draws <- matrix(NA_real_, nrow = N_BOOT, ncol = 2, dimnames = list(NULL, c("r0", "r1")))
  for (b in seq_len(N_BOOT)) {
    pick <- sample.int(n, n, replace = TRUE)
    wb <- w[pick]
    sw <- sum(wb)
    draws[b, 1] <- sum(wb * s0[pick]) / sw
    draws[b, 2] <- sum(wb * s1[pick]) / sw
  }
  draws
}

summarise_bootstrap_risks <- function(draws) {
  if (is.null(draws)) {
    return(tibble(
      n_bootstrap = 0L, n_bootstrap_usable = 0L,
      rd_boot_low = NA_real_, rd_boot_high = NA_real_,
      rr_boot_low = NA_real_, rr_boot_high = NA_real_,
      bootstrap_note = if (N_BOOT <= 0) "not run: N_BOOT = 0" else "not run: too few usable scores"
    ))
  }
  rd <- draws[, "r1"] - draws[, "r0"]
  usable_rr <- is.finite(draws[, "r0"]) & draws[, "r0"] > 0 &
    is.finite(draws[, "r1"]) & draws[, "r1"] > 0
  rr <- rep(NA_real_, nrow(draws))
  rr[usable_rr] <- draws[usable_rr, "r1"] / draws[usable_rr, "r0"]

  tibble(
    n_bootstrap = nrow(draws),
    n_bootstrap_usable = sum(usable_rr),
    rd_boot_low = as.numeric(quantile(rd, 0.025, na.rm = TRUE, names = FALSE)),
    rd_boot_high = as.numeric(quantile(rd, 0.975, na.rm = TRUE, names = FALSE)),
    rr_boot_low = if (any(usable_rr)) as.numeric(quantile(rr, 0.025, na.rm = TRUE, names = FALSE)) else NA_real_,
    rr_boot_high = if (any(usable_rr)) as.numeric(quantile(rr, 0.975, na.rm = TRUE, names = FALSE)) else NA_real_,
    bootstrap_note = paste(
      "percentile interval from participant resampling of the doubly robust",
      "arm-risk scores; conditional on the fitted nuisance functions"
    )
  )
}

# Crude, unadjusted arm risks. These are always finite whenever an arm contains
# at least one person, so they give an E-value floor for transitions where the
# doubly robust RR cannot be formed. They are NOT confounder-adjusted and must
# be reported as a descriptive comparison, not as a causal estimate.
crude_risk_scale <- function(Y, W, idx) {
  y <- Y[idx]; w <- W[idx]
  n0 <- sum(w == 0); n1 <- sum(w == 1)
  if (n0 == 0 || n1 == 0) {
    return(tibble(
      crude_risk_stayed = NA_real_, crude_risk_moved = NA_real_,
      crude_rr = NA_real_, crude_evalue = NA_real_,
      crude_note = "one arm was empty in this estimand"
    ))
  }
  r0 <- mean(y[w == 0]); r1 <- mean(y[w == 1])
  rr <- if (r0 > 0) r1 / r0 else NA_real_
  tibble(
    crude_risk_stayed = r0,
    crude_risk_moved = r1,
    crude_rr = rr,
    crude_evalue = evalue_from_rr(rr),
    crude_note = if (is.finite(rr)) {
      "unadjusted arm risks; reported as an E-value floor, not a causal estimate"
    } else {
      "unadjusted RR undefined because the reference arm had no events"
    }
  )
}

# =============================================================================
# 7. Quantitative bias analysis for a generic binary confounder U
# =============================================================================
# Assumptions throughout this section: U is binary; its risk ratio with the
# outcome, RR_UY, is the same in both arms; there is no U-by-treatment
# interaction; and no other bias mechanism is represented. p_u_stayed and
# p_u_moved are residual prevalences AFTER adjustment for the measured
# covariates X, which is the quantity that matters here.
#
# Both directions are included because we do not know whether U is more common
# among participants who stayed or among participants who moved.
qba_scenarios <- tribble(
  ~scenario,                 ~p_u_stayed, ~p_u_moved, ~rr_u_y,
  "no prevalence imbalance",         0.10,        0.10,     2.00,
  "weak: U more common moved",       0.10,        0.15,     1.50,
  "moderate: U more common moved",   0.10,        0.25,     2.00,
  "strong: U more common moved",     0.10,        0.35,     3.00,
  "weak: U more common stayed",      0.15,        0.10,     1.50,
  "moderate: U more common stayed",  0.25,        0.10,     2.00,
  "strong: U more common stayed",    0.35,        0.10,     3.00
)

# Remove the confounding produced by an assumed prevalence imbalance: strip U
# out of each observed arm risk, then standardise both arms to one common U
# prevalence. Vectorised over r0/r1 so the same call serves the point estimate
# and every bootstrap replicate.
qba_adjust <- function(r0, r1, p_u_stayed, p_u_moved, rr_u_y, p_u_standard = NULL) {
  if (is.null(p_u_standard)) p_u_standard <- (p_u_stayed + p_u_moved) / 2
  f0 <- 1 + p_u_stayed * (rr_u_y - 1)
  f1 <- 1 + p_u_moved * (rr_u_y - 1)
  fs <- 1 + p_u_standard * (rr_u_y - 1)
  r0_adj <- r0 / f0 * fs
  r1_adj <- r1 / f1 * fs
  list(
    p_u_standard = p_u_standard,
    risk_stayed_adjusted = r0_adj,
    risk_moved_adjusted = r1_adj,
    rd_adjusted = r1_adj - r0_adj,
    rr_adjusted = r1_adj / r0_adj
  )
}

# Exact tipping point. Standardising both arms multiplies them by the same
# factor, so the adjusted risk difference is zero when
#   r1 / (1 + p1 (R - 1)) = r0 / (1 + p0 (R - 1)),
# which rearranges to
#   R = 1 + (r0 - r1) / (r1 * p0 - r0 * p1).
# No search or simulation is needed. A value above 1 describes a
# risk-increasing U; a value in (0, 1) describes a protective U.
qba_tipping_rr <- function(r0, r1, p_u_stayed, p_u_moved) {
  denominator <- r1 * p_u_stayed - r0 * p_u_moved
  if (!is.finite(denominator) || abs(denominator) < 1e-12) return(NA_real_)
  1 + (r0 - r1) / denominator
}

run_qba_scenarios <- function(r0, r1, boot_draws, spec, scenarios = qba_scenarios) {
  if (!is.finite(r0) || !is.finite(r1)) {
    return(prefix_result(
      scenarios |>
        mutate(
          p_u_standard = (p_u_stayed + p_u_moved) / 2,
          risk_stayed_observed = r0, risk_moved_observed = r1,
          risk_stayed_adjusted = NA_real_, risk_moved_adjusted = NA_real_,
          rd_observed = NA_real_, rd_adjusted = NA_real_,
          rd_adjusted_low = NA_real_, rd_adjusted_high = NA_real_,
          rr_observed = NA_real_, rr_adjusted = NA_real_,
          result = "not calculated: the doubly robust arm risks were not finite",
          inference = "no bootstrap interval because the point estimate is missing"
        ),
      spec
    ))
  }

  rd_observed <- r1 - r0
  rr_observed <- if (r0 > 0) r1 / r0 else NA_real_

  rows <- lapply(seq_len(nrow(scenarios)), function(k) {
    s <- scenarios[k, ]
    point <- qba_adjust(r0, r1, s$p_u_stayed, s$p_u_moved, s$rr_u_y)

    if (is.null(boot_draws)) {
      low <- NA_real_; high <- NA_real_; p_boot <- NA_real_
      inference <- if (N_BOOT <= 0) "point estimate only; N_BOOT = 0" else
        "point estimate only; bootstrap unavailable for this contrast"
    } else {
      boot <- qba_adjust(
        boot_draws[, "r0"], boot_draws[, "r1"],
        s$p_u_stayed, s$p_u_moved, s$rr_u_y
      )
      values <- boot$rd_adjusted
      values <- values[is.finite(values)]
      if (length(values) < 10) {
        low <- NA_real_; high <- NA_real_; p_boot <- NA_real_
        inference <- "point estimate only; too few usable bootstrap replicates"
      } else {
        low <- as.numeric(quantile(values, 0.025, names = FALSE))
        high <- as.numeric(quantile(values, 0.975, names = FALSE))
        p_boot <- min(1, 2 * min(mean(values <= 0), mean(values >= 0)))
        inference <- "percentile interval from the arm-risk score bootstrap"
      }
    }

    tibble(
      scenario = s$scenario,
      p_u_stayed = s$p_u_stayed,
      p_u_moved = s$p_u_moved,
      rr_u_y = s$rr_u_y,
      p_u_standard = point$p_u_standard,
      risk_stayed_observed = r0,
      risk_moved_observed = r1,
      risk_stayed_adjusted = point$risk_stayed_adjusted,
      risk_moved_adjusted = point$risk_moved_adjusted,
      rd_observed = rd_observed,
      rd_adjusted = point$rd_adjusted,
      rd_adjusted_low = low,
      rd_adjusted_high = high,
      rd_adjusted_p_bootstrap = p_boot,
      rr_observed = rr_observed,
      rr_adjusted = point$rr_adjusted,
      adjusted_risks_in_range = all(
        is.finite(c(point$risk_stayed_adjusted, point$risk_moved_adjusted)) &
          c(point$risk_stayed_adjusted, point$risk_moved_adjusted) > 0 &
          c(point$risk_stayed_adjusted, point$risk_moved_adjusted) < 1
      ),
      result = case_when(
        !is.finite(point$rd_adjusted) ~ "not calculated",
        sign(point$rd_adjusted) != sign(rd_observed) ~ "direction reversed",
        abs(point$rd_adjusted) < abs(rd_observed) ~ "moved toward null",
        abs(point$rd_adjusted) > abs(rd_observed) ~ "moved away from null",
        TRUE ~ "unchanged"
      ),
      interval_result = case_when(
        !is.finite(low) | !is.finite(high) ~ "no interval available",
        low <= 0 & high >= 0 ~ "adjusted interval includes zero",
        TRUE ~ "adjusted interval excludes zero"
      ),
      inference = inference
    )
  })

  prefix_result(bind_rows(rows), spec)
}

# Factorial grid. Seven prespecified scenarios show seven points; the grid shows
# the surface those points sit on, which is what a robustness claim needs.
qba_grid_specification <- expand_grid(
  p_u_stayed = c(0.05, 0.10, 0.20, 0.30, 0.50),
  prevalence_difference = c(0.00, 0.02, 0.05, 0.10, 0.15, 0.20),
  rr_u_y = c(1.25, 1.50, 2.00, 3.00, 5.00),
  imbalance_direction = c("U more common moved", "U more common stayed")
) |>
  mutate(
    p_u_moved = if_else(
      imbalance_direction == "U more common moved",
      p_u_stayed + prevalence_difference,
      p_u_stayed - prevalence_difference
    )
  ) |>
  filter(p_u_moved > 0, p_u_moved < 1) |>
  # A zero prevalence difference is the same cell in both directions.
  filter(!(prevalence_difference == 0 & imbalance_direction == "U more common stayed")) |>
  distinct(p_u_stayed, p_u_moved, rr_u_y, .keep_all = TRUE)

run_qba_grid <- function(r0, r1, spec, grid = qba_grid_specification) {
  if (!is.finite(r0) || !is.finite(r1)) {
    return(prefix_result(
      grid |>
        mutate(
          rd_observed = NA_real_, rd_adjusted = NA_real_,
          rr_adjusted = NA_real_, attenuation_ratio = NA_real_,
          result = "not calculated: the doubly robust arm risks were not finite"
        ),
      spec
    ))
  }
  rd_observed <- r1 - r0
  adjusted <- qba_adjust(r0, r1, grid$p_u_stayed, grid$p_u_moved, grid$rr_u_y)

  prefix_result(
    grid |>
      mutate(
        p_u_standard = adjusted$p_u_standard,
        risk_stayed_observed = r0,
        risk_moved_observed = r1,
        rd_observed = rd_observed,
        rd_adjusted = adjusted$rd_adjusted,
        rr_adjusted = adjusted$rr_adjusted,
        attenuation_ratio = adjusted$rd_adjusted / rd_observed,
        result = case_when(
          !is.finite(adjusted$rd_adjusted) ~ "not calculated",
          sign(adjusted$rd_adjusted) != sign(rd_observed) ~ "direction reversed",
          abs(adjusted$rd_adjusted) < abs(rd_observed) ~ "moved toward null",
          abs(adjusted$rd_adjusted) > abs(rd_observed) ~ "moved away from null",
          TRUE ~ "unchanged"
        )
      ),
    spec
  )
}

tipping_point_specification <- expand_grid(
  p_u_stayed = c(0.05, 0.10, 0.20, 0.30, 0.50),
  prevalence_difference = c(0.02, 0.05, 0.10, 0.15, 0.20),
  imbalance_direction = c("U more common moved", "U more common stayed")
) |>
  mutate(
    p_u_moved = if_else(
      imbalance_direction == "U more common moved",
      p_u_stayed + prevalence_difference,
      p_u_stayed - prevalence_difference
    )
  ) |>
  filter(p_u_moved > 0, p_u_moved < 1)

run_tipping_points <- function(r0, r1, spec, grid = tipping_point_specification) {
  if (!is.finite(r0) || !is.finite(r1)) {
    return(prefix_result(
      grid |>
        mutate(
          rr_u_y_tipping = NA_real_,
          tipping_interpretation = "not calculated: the doubly robust arm risks were not finite"
        ),
      spec
    ))
  }

  tipping <- vapply(
    seq_len(nrow(grid)),
    function(k) qba_tipping_rr(r0, r1, grid$p_u_stayed[k], grid$p_u_moved[k]),
    numeric(1)
  )

  prefix_result(
    grid |>
      mutate(
        risk_stayed_observed = r0,
        risk_moved_observed = r1,
        rd_observed = r1 - r0,
        rr_u_y_tipping = tipping,
        tipping_interpretation = case_when(
          !is.finite(tipping) ~ "no finite solution: this prevalence pattern cannot null the estimate",
          tipping > 1 ~ "a risk-increasing U of at least this strength would null the estimate",
          tipping > 0 & tipping < 1 ~ "only a protective U of this strength would null the estimate",
          TRUE ~ "solution is not an admissible risk ratio"
        )
      ),
    spec
  )
}

# =============================================================================
# 8. Contrast-level estimation
# =============================================================================

ESTIMAND_DEFINITIONS <- function() {
  c(
    paste0("full_sample_", ANALYSIS_LABEL),
    "overlap_weighted",
    "trim_001_099",
    "trim_005_095",
    "trim_010_090"
  )
}

empty_result <- function(spec, status, reason, cohort_flow, arm_support) {
  list(
    status = prefix_result(tibble(status = status, reason = reason), spec),
    cohort_flow = cohort_flow,
    arm_support = arm_support,
    propensity = tibble(),
    ate = tibble(),
    evalue = tibble(),
    evalue_subgroup = tibble(),
    qba = tibble(),
    qba_grid = tibble(),
    tipping = tibble(),
    calibration = tibble(),
    blp = tibble(),
    rate = tibble(),
    variable_importance = tibble(),
    cate_quartiles = tibble(),
    subgroup_ate = tibble(),
    subgroup_difference = tibble(),
    forest = NULL
  )
}

propensity_summary <- function(propensity, W, spec, n_floored) {
  bind_rows(lapply(c(0, 1), function(arm) {
    values <- propensity[W == arm]
    tibble(
      arm = arm,
      arm_label = ifelse(arm == 1, "moved", "stayed"),
      n = length(values),
      minimum = min(values, na.rm = TRUE),
      q01 = quantile(values, 0.01, na.rm = TRUE, names = FALSE),
      q05 = quantile(values, 0.05, na.rm = TRUE, names = FALSE),
      median = median(values, na.rm = TRUE),
      q95 = quantile(values, 0.95, na.rm = TRUE, names = FALSE),
      q99 = quantile(values, 0.99, na.rm = TRUE, names = FALSE),
      maximum = max(values, na.rm = TRUE),
      n_outside_propensity_floor = n_floored,
      propensity_floor = PROPENSITY_FLOOR
    )
  })) |>
    prefix_result(spec)
}

run_transition <- function(spec, index) {
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
      ifelse(sum(Y[W == 0] == 1) == 0, "zero events; the risk ratio cannot be formed", "events observed"),
      ifelse(sum(Y[W == 1] == 1) == 0, "zero events; the risk ratio cannot be formed", "events observed")
    )
  ), spec)

  log_msg(
    "---- ", spec$domain, " | ", spec$health_direction, " | ", spec$transition,
    " | stayed ", sum(W == 0), " moved ", sum(W == 1), " events ", sum(Y == 1)
  )

  if (nrow(X) == 0 || length(unique(W)) < 2) {
    return(empty_result(spec, "not_estimable", "reference or moved arm absent", cohort_flow, arm_support))
  }
  if (length(unique(Y)) < 2) {
    return(empty_result(spec, "not_estimable", "outcome has no variation", cohort_flow, arm_support))
  }

  set.seed(RANDOM_SEED)
  forest <- tryCatch(
    causal_forest(
      X = X, Y = Y, W = W,
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

  scores <- arm_risk_scores(forest, Y, W)
  propensity <- scores$propensity
  propensity_table <- propensity_summary(propensity, W, spec, scores$n_floored)

  overlap_weights <- propensity * (1 - propensity)
  estimand_names <- ESTIMAND_DEFINITIONS()
  estimand_subsets <- list(
    rep(TRUE, length(W)),
    rep(TRUE, length(W)),
    propensity > 0.01 & propensity < 0.99,
    propensity > 0.05 & propensity < 0.95,
    propensity > 0.10 & propensity < 0.90
  )
  estimand_weights <- list(NULL, overlap_weights, NULL, NULL, NULL)
  estimand_targets <- c("all", "overlap", "all", "all", "all")

  ate_rows <- list()
  evalue_rows <- list()
  primary_r0 <- NA_real_
  primary_r1 <- NA_real_
  primary_boot <- NULL

  for (k in seq_along(estimand_names)) {
    estimand <- estimand_names[k]
    subset_flag <- estimand_subsets[[k]]
    idx <- which(subset_flag)

    if (length(idx) == 0) {
      ate_rows[[estimand]] <- prefix_result(tibble(
        estimand = estimand, n = 0L, n_stayed = 0L, n_moved = 0L, events = 0L,
        error_message = "empty estimand subset"
      ), spec)
      next
    }

    grf_ate <- tryCatch({
      raw <- if (estimand_targets[k] == "overlap") {
        average_treatment_effect(forest, target.sample = "overlap")
      } else {
        average_treatment_effect(forest, subset = subset_flag, target.sample = "all")
      }
      extract_ate(raw)
    }, error = function(e) tibble(error_message = e$message))

    risk_scale <- risk_scale_from_scores(scores, idx, estimand_weights[[k]])
    boot_draws <- bootstrap_arm_risks(
      scores, idx, estimand_weights[[k]],
      seed_offset = index * 100 + k
    )
    boot_summary <- summarise_bootstrap_risks(boot_draws)
    crude <- crude_risk_scale(Y, W, idx)

    ate_rows[[estimand]] <- grf_ate |>
      mutate(
        estimand = estimand,
        n = length(idx),
        n_stayed = sum(W[idx] == 0),
        n_moved = sum(W[idx] == 1),
        events = sum(Y[idx] == 1),
        events_stayed = sum(Y[idx][W[idx] == 0] == 1),
        events_moved = sum(Y[idx][W[idx] == 1] == 1),
        retained_fraction = length(idx) / length(W),
        .before = 1
      ) |>
      bind_cols(
        risk_scale |> select(rd_score, rd_score_se),
        boot_summary |> select(n_bootstrap, rd_boot_low, rd_boot_high)
      ) |>
      mutate(
        # The AIPW arm risks are algebraically the same doubly robust scores grf
        # averages, so this difference should be numerically zero. A non-zero
        # value means the two paths have diverged and the risk-scale results for
        # that row must not be trusted.
        ate_consistency_diff = if ("estimate" %in% names(grf_ate)) {
          rd_score - grf_ate$estimate[1]
        } else {
          NA_real_
        },
        inference = if (estimand == PRIMARY_ESTIMAND) {
          "primary estimand; Wald interval, with a score bootstrap interval alongside"
        } else {
          "sensitivity estimand"
        }
      ) |>
      prefix_result(spec)

    evalue_rows[[estimand]] <- risk_scale |>
      mutate(
        estimand = estimand,
        n = length(idx),
        n_stayed = sum(W[idx] == 0),
        n_moved = sum(W[idx] == 1),
        events_stayed = sum(Y[idx][W[idx] == 0] == 1),
        events_moved = sum(Y[idx][W[idx] == 1] == 1),
        .before = 1
      ) |>
      bind_cols(boot_summary, crude) |>
      mutate(
        evalue_point = evalue_from_rr(rr),
        evalue_ci = evalue_from_ci(rr, rr_low, rr_high),
        evalue_ci_bootstrap = evalue_from_ci(rr, rr_boot_low, rr_boot_high),
        evalue_ci_basis = evalue_ci_note(rr, rr_low, rr_high),
        evalue_reported = dplyr::coalesce(evalue_point, crude_evalue),
        evalue_reported_source = if_else(
          is.finite(evalue_point), "doubly robust RR", "crude RR fallback"
        ),
        interpretation = paste(
          "the minimum risk-ratio strength an unmeasured confounder would need",
          "with BOTH the transition and CVD, above and beyond the measured",
          "covariates, to move this estimate to the null"
        )
      ) |>
      prefix_result(spec)

    if (estimand == PRIMARY_ESTIMAND) {
      primary_r0 <- risk_scale$risk_stayed[1]
      primary_r1 <- risk_scale$risk_moved[1]
      primary_boot <- boot_draws
    }
  }

  ate <- bind_rows(ate_rows)
  evalue <- bind_rows(evalue_rows)

  # The QBA is anchored to the primary estimand, so that the bias analysis and
  # the headline effect describe the same target population.
  qba <- run_qba_scenarios(primary_r0, primary_r1, primary_boot, spec)
  qba_grid <- run_qba_grid(primary_r0, primary_r1, spec)
  tipping <- run_tipping_points(primary_r0, primary_r1, spec)

  # ---------------------------------------------------------------------------
  # Subgroup E-values, on the same seven subgroups the primary script reports.
  # ---------------------------------------------------------------------------
  evalue_subgroup <- tibble()
  if (RUN_SUBGROUP_EVALUES) {
    groups <- make_subgroups(subgroup_data)
    evalue_subgroup <- bind_rows(lapply(names(groups), function(group_name) {
      idx <- groups[[group_name]]
      if (length(idx) < 2 || length(unique(W[idx])) < 2) {
        return(tibble(
          subgroup = group_name, n = length(idx),
          n_stayed = sum(W[idx] == 0), n_moved = sum(W[idx] == 1),
          events = sum(Y[idx] == 1),
          rr = NA_real_, rr_low = NA_real_, rr_high = NA_real_,
          evalue_point = NA_real_, evalue_ci = NA_real_,
          risk_scale_status = "subgroup contained fewer than two people or only one arm"
        ))
      }
      rs <- risk_scale_from_scores(scores, idx)
      rs |>
        mutate(
          subgroup = group_name,
          n = length(idx),
          n_stayed = sum(W[idx] == 0),
          n_moved = sum(W[idx] == 1),
          events = sum(Y[idx] == 1),
          .before = 1
        ) |>
        mutate(
          evalue_point = evalue_from_rr(rr),
          evalue_ci = evalue_from_ci(rr, rr_low, rr_high),
          evalue_ci_basis = evalue_ci_note(rr, rr_low, rr_high)
        )
    })) |>
      prefix_result(spec)
  }

  # ---------------------------------------------------------------------------
  # Heterogeneity output. These estimators were already being run in the previous
  # version of this file, but only printed to the log, so nothing survived the
  # run. They are written to CSV here.
  # ---------------------------------------------------------------------------
  calibration <- tibble(); blp <- tibble(); rate <- tibble()
  # Named with a _tbl suffix so it does not shadow grf's variable_importance().
  variable_importance_tbl <- tibble(); cate_quartiles <- tibble()
  subgroup_ate <- tibble(); subgroup_difference <- tibble()

  if (RUN_HETEROGENEITY) {
    calibration <- tryCatch(
      prefix_result(tidy_inference(test_calibration(forest), "test_calibration"), spec),
      error = function(e) prefix_result(tibble(source = "test_calibration", error_message = e$message), spec)
    )

    modifier_matrix <- cbind(
      male_sex_1 = as.numeric(subgroup_data[, "sex"] == 1),
      age_ge_60 = as.numeric(subgroup_data[, "age"] >= AGE_CUT),
      bmi_25_to_lt_30 = as.numeric(subgroup_data[, "bmi"] >= 25 & subgroup_data[, "bmi"] < 30),
      bmi_ge_30 = as.numeric(subgroup_data[, "bmi"] >= 30)
    )
    modifier_varies <- apply(modifier_matrix, 2, function(x) length(unique(x[is.finite(x)])) > 1)
    modifier_matrix <- modifier_matrix[, modifier_varies, drop = FALSE]

    blp <- if (ncol(modifier_matrix) == 0) {
      prefix_result(tibble(source = "best_linear_projection",
                           error_message = "prespecified modifiers had no variation"), spec)
    } else {
      tryCatch(
        prefix_result(tidy_inference(
          best_linear_projection(forest, A = modifier_matrix), "best_linear_projection"
        ), spec),
        error = function(e) prefix_result(tibble(source = "best_linear_projection",
                                                 error_message = e$message), spec)
      )
    }

    prediction <- tryCatch(predict(forest, estimate.variance = TRUE), error = function(e) NULL)
    cate <- if (!is.null(prediction)) as.numeric(prediction$predictions) else rep(NA_real_, nrow(X))

    # Priorities follow the clinical direction, matching single_variable3.R:
    # for an improvement the benefit is a NEGATIVE risk difference.
    priorities <- if (spec$health_direction == "improvement") {
      -cate
    } else if (spec$health_direction == "deterioration") {
      cate
    } else {
      abs(cate)
    }

    rate <- tryCatch(
      prefix_result(tidy_inference(
        rank_average_treatment_effect(forest, priorities = priorities, target = "AUTOC"),
        "rank_average_treatment_effect"
      ), spec),
      error = function(e) prefix_result(tibble(source = "rank_average_treatment_effect",
                                               error_message = e$message), spec)
    )

    variable_importance_tbl <- tryCatch(
      prefix_result(tibble(
        variable = colnames(X),
        importance = as.numeric(variable_importance(forest))
      ) |> arrange(desc(importance)), spec),
      error = function(e) prefix_result(tibble(error_message = e$message), spec)
    )

    trim_005 <- propensity > 0.05 & propensity < 0.95
    quartile <- dplyr::ntile(cate, 4)
    cate_quartiles <- prefix_result(bind_rows(lapply(1:4, function(q) {
      idx_flag <- trim_005 & quartile == q & !is.na(quartile)
      if (sum(idx_flag) < 50) {
        return(tibble(
          quartile = q, n = sum(idx_flag),
          estimate = NA_real_, std_error = NA_real_, low = NA_real_, high = NA_real_,
          note = "fewer than 50 people in this quartile after 0.05-0.95 trimming"
        ))
      }
      est <- tryCatch(extract_ate(average_treatment_effect(forest, subset = idx_flag)),
                      error = function(e) tibble(estimate = NA_real_, std_error = NA_real_,
                                                 low = NA_real_, high = NA_real_))
      est |>
        mutate(quartile = q, n = sum(idx_flag), note = "exploratory; quartiles are data-defined", .before = 1)
    })), spec)

    groups <- make_subgroups(subgroup_data)
    subgroup_ate <- prefix_result(bind_rows(lapply(names(groups), function(group_name) {
      idx <- groups[[group_name]]
      if (length(idx) < 50 || length(unique(W[idx])) < 2) {
        return(tibble(
          subgroup = group_name, n = length(idx),
          n_stayed = sum(W[idx] == 0), n_moved = sum(W[idx] == 1),
          events = sum(Y[idx] == 1),
          estimate = NA_real_, std_error = NA_real_, low = NA_real_, high = NA_real_,
          error_message = "subgroup too small or single-armed"
        ))
      }
      est <- tryCatch(extract_ate(average_treatment_effect(forest, subset = idx)),
                      error = function(e) tibble(error_message = e$message))
      est |>
        mutate(
          subgroup = group_name, n = length(idx),
          n_stayed = sum(W[idx] == 0), n_moved = sum(W[idx] == 1),
          events = sum(Y[idx] == 1), .before = 1
        )
    })), spec)

    subgroup_pairs <- list(
      list(label = "sex: male - female", high = "male_sex_1", low = "female_sex_0"),
      list(label = paste0("age: >=", AGE_CUT, " - <", AGE_CUT), high = "age_ge_60", low = "age_lt_60"),
      list(label = "bmi: >=30 - <25", high = "bmi_ge_30", low = "bmi_lt_25")
    )
    subgroup_difference <- prefix_result(bind_rows(lapply(subgroup_pairs, function(cmp) {
      row_high <- subgroup_ate[subgroup_ate$subgroup == cmp$high, ]
      row_low <- subgroup_ate[subgroup_ate$subgroup == cmp$low, ]
      if (nrow(row_high) == 0 || nrow(row_low) == 0 ||
          !"estimate" %in% names(subgroup_ate) ||
          !is.finite(row_high$estimate[1]) || !is.finite(row_low$estimate[1])) {
        return(tibble(comparison = cmp$label, difference = NA_real_, std_error = NA_real_,
                      z = NA_real_, p_value = NA_real_,
                      inference = "not estimable in one or both subgroups"))
      }
      d <- row_high$estimate[1] - row_low$estimate[1]
      se <- sqrt(row_high$std_error[1]^2 + row_low$std_error[1]^2)
      z <- d / se
      tibble(
        comparison = cmp$label, difference = d, std_error = se, z = z,
        p_value = 2 * pnorm(abs(z), lower.tail = FALSE),
        inference = paste(
          "independent-subgroup normal approximation; treats the two subgroup",
          "estimates as independent, which is conservative only if they are"
        )
      )
    })), spec)
  }

  list(
    status = prefix_result(tibble(status = "fitted", reason = NA_character_), spec),
    cohort_flow = cohort_flow,
    arm_support = arm_support,
    propensity = propensity_table,
    ate = ate,
    evalue = evalue,
    evalue_subgroup = evalue_subgroup,
    qba = qba,
    qba_grid = qba_grid,
    tipping = tipping,
    calibration = calibration,
    blp = blp,
    rate = rate,
    variable_importance = variable_importance_tbl,
    cate_quartiles = cate_quartiles,
    subgroup_ate = subgroup_ate,
    subgroup_difference = subgroup_difference,
    forest = if (SAVE_FORESTS) forest else NULL
  )
}

# =============================================================================
# 9. Run every transition
# =============================================================================

all_results <- list()
for (i in seq_len(nrow(contrasts))) {
  spec <- contrasts[i, ]
  all_results[[spec$contrast_id]] <- tryCatch(
    run_transition(spec, i),
    error = function(e) {
      log_msg("FAILED: ", spec$contrast_id, " - ", e$message)
      empty_result(spec, "run_failed", e$message, tibble(), tibble())
    }
  )
}

# =============================================================================
# 10. Collect, apply FDR, and write
# =============================================================================

status_tbl <- safe_bind(all_results, "status")
cohort_flow_tbl <- safe_bind(all_results, "cohort_flow")
arm_support_tbl <- safe_bind(all_results, "arm_support")
propensity_tbl <- safe_bind(all_results, "propensity")
ate_tbl <- safe_bind(all_results, "ate")
evalue_tbl <- safe_bind(all_results, "evalue")
evalue_subgroup_tbl <- safe_bind(all_results, "evalue_subgroup")
qba_tbl <- safe_bind(all_results, "qba")
qba_grid_tbl <- safe_bind(all_results, "qba_grid")
tipping_tbl <- safe_bind(all_results, "tipping")
calibration_tbl <- safe_bind(all_results, "calibration")
blp_tbl <- safe_bind(all_results, "blp")
rate_tbl <- safe_bind(all_results, "rate")
variable_importance_tbl <- safe_bind(all_results, "variable_importance")
cate_quartile_tbl <- safe_bind(all_results, "cate_quartiles")
subgroup_ate_tbl <- safe_bind(all_results, "subgroup_ate")
subgroup_difference_tbl <- safe_bind(all_results, "subgroup_difference")

# Benjamini-Hochberg within the primary-estimand family only. Applying it across
# all five estimands would count the same transition five times. An E-value is
# only worth quoting for a transition whose effect survived this correction, so
# the q-value is carried into the E-value table.
if (nrow(ate_tbl) > 0 && "estimand" %in% names(ate_tbl)) {
  ate_tbl <- add_fdr_subset(
    ate_tbl,
    ate_tbl$estimand == PRIMARY_ESTIMAND,
    output_col = "p_fdr_primary_family"
  )
}

if (nrow(evalue_tbl) > 0 && nrow(ate_tbl) > 0 &&
    all(c("contrast_id", "estimand") %in% names(ate_tbl))) {
  ate_join <- ate_tbl |>
    select(any_of(c("contrast_id", "estimand", "estimate", "std_error",
                    "low", "high", "p_value", "p_fdr_primary_family",
                    "ate_consistency_diff"))) |>
    rename_with(~ paste0("ate_", .x), any_of(c("estimate", "std_error", "low", "high", "p_value")))

  evalue_tbl <- evalue_tbl |>
    left_join(ate_join, by = c("contrast_id", "estimand"))

  # add_fdr_subset() returns the table untouched when no p_value column survived,
  # which happens if every ATE call errored. Creating the column explicitly keeps
  # the mutate below valid instead of failing on a missing name.
  if (!"p_fdr_primary_family" %in% names(evalue_tbl)) {
    evalue_tbl$p_fdr_primary_family <- NA_real_
  }

  evalue_tbl <- evalue_tbl |>
    mutate(
      fdr_supported = case_when(
        !is.finite(p_fdr_primary_family) ~ NA,
        p_fdr_primary_family < 0.05 ~ TRUE,
        TRUE ~ FALSE
      ),
      evalue_relevance = case_when(
        is.na(fdr_supported) ~ "not applicable: this row is a sensitivity estimand",
        fdr_supported ~ "transition survived FDR; the E-value is the quantity to report",
        TRUE ~ paste(
          "transition did not survive FDR, so there is no positive finding for a",
          "confounder to explain away; the E-value is reported for completeness only"
        )
      )
    )
}

safe_write(cohort_summary, "cohort_summary.csv")
safe_write(contrasts, "transition_registry.csv")
safe_write(status_tbl, "contrast_status.csv", "no contrast produced a status row")
safe_write(cohort_flow_tbl, "contrast_participant_flow.csv")
safe_write(arm_support_tbl, "arm_counts_events.csv")
safe_write(propensity_tbl, "propensity_summary.csv")
safe_write(ate_tbl, "ate_estimands_with_risk_scale.csv")
safe_write(evalue_tbl, "evalue_summary.csv", "no transition produced a risk-scale result")
safe_write(evalue_subgroup_tbl, "evalue_subgroups.csv",
           "subgroup E-values were switched off or no subgroup was estimable")
safe_write(qba_tbl, "qba_prespecified_scenarios.csv")
safe_write(qba_grid_tbl, "qba_factorial_grid.csv")
safe_write(tipping_tbl, "qba_tipping_points.csv")
safe_write(calibration_tbl, "heterogeneity_test_calibration.csv",
           "heterogeneity output was switched off (RUN_HETEROGENEITY=0)")
safe_write(blp_tbl, "heterogeneity_best_linear_projection.csv",
           "heterogeneity output was switched off (RUN_HETEROGENEITY=0)")
safe_write(rate_tbl, "heterogeneity_rate_autoc.csv",
           "heterogeneity output was switched off (RUN_HETEROGENEITY=0)")
safe_write(variable_importance_tbl, "variable_importance_descriptive.csv",
           "heterogeneity output was switched off (RUN_HETEROGENEITY=0)")
safe_write(cate_quartile_tbl, "cate_quartiles_exploratory.csv",
           "heterogeneity output was switched off (RUN_HETEROGENEITY=0)")
safe_write(subgroup_ate_tbl, "subgroup_ate.csv",
           "heterogeneity output was switched off (RUN_HETEROGENEITY=0)")
safe_write(subgroup_difference_tbl, "subgroup_difference_tests.csv",
           "heterogeneity output was switched off (RUN_HETEROGENEITY=0)")

# ---------------------------------------------------------------------------
# Headline table: one row per transition, on the primary estimand.
# ---------------------------------------------------------------------------
headline <- tibble()
if (nrow(evalue_tbl) > 0 && "estimand" %in% names(evalue_tbl)) {
  headline <- evalue_tbl |>
    filter(estimand == PRIMARY_ESTIMAND) |>
    select(any_of(c(
      "analysis_label", "contrast_id", "domain", "transition", "health_direction",
      "interpretation_flag", "n", "n_stayed", "n_moved", "events_stayed", "events_moved",
      "ate_estimate", "ate_low", "ate_high", "ate_p_value", "p_fdr_primary_family",
      "fdr_supported", "risk_stayed", "risk_moved", "rr", "rr_low", "rr_high",
      "rr_boot_low", "rr_boot_high", "evalue_point", "evalue_ci",
      "evalue_ci_bootstrap", "evalue_reported", "evalue_reported_source",
      "crude_rr", "crude_evalue", "risk_scale_status", "evalue_relevance"
    )))
}
safe_write(headline, "table_evalue_headline.csv",
           "no transition produced a primary-estimand risk-scale result")

# ---------------------------------------------------------------------------
# QBA reading aid: for each transition, the weakest prespecified scenario that
# reverses the direction, and the smallest tipping-point RR_UY across the grid.
# ---------------------------------------------------------------------------
qba_summary <- tibble()
if (nrow(qba_tbl) > 0 && "result" %in% names(qba_tbl)) {
  reversal <- qba_tbl |>
    filter(result == "direction reversed") |>
    group_by(contrast_id) |>
    slice_min(rr_u_y, n = 1, with_ties = FALSE) |>
    ungroup() |>
    select(contrast_id, weakest_reversing_scenario = scenario,
           weakest_reversing_rr_u_y = rr_u_y)

  crossing_null <- qba_tbl |>
    group_by(contrast_id) |>
    summarise(
      n_scenarios = n(),
      n_direction_reversed = sum(result == "direction reversed", na.rm = TRUE),
      n_toward_null = sum(result == "moved toward null", na.rm = TRUE),
      n_interval_includes_zero = sum(interval_result == "adjusted interval includes zero", na.rm = TRUE),
      .groups = "drop"
    )

  tipping_summary <- if (nrow(tipping_tbl) > 0 && "rr_u_y_tipping" %in% names(tipping_tbl)) {
    tipping_tbl |>
      filter(is.finite(rr_u_y_tipping), rr_u_y_tipping > 1) |>
      group_by(contrast_id) |>
      summarise(
        smallest_tipping_rr_u_y = min(rr_u_y_tipping),
        smallest_tipping_p_u_stayed = p_u_stayed[which.min(rr_u_y_tipping)],
        smallest_tipping_p_u_moved = p_u_moved[which.min(rr_u_y_tipping)],
        .groups = "drop"
      )
  } else {
    tibble(contrast_id = character())
  }

  qba_summary <- qba_tbl |>
    distinct(analysis_label, contrast_id, domain, transition, health_direction) |>
    left_join(crossing_null, by = "contrast_id") |>
    left_join(reversal, by = "contrast_id") |>
    left_join(tipping_summary, by = "contrast_id") |>
    mutate(
      reading = case_when(
        is.na(n_scenarios) ~ "no scenario was calculated",
        n_direction_reversed > 0 ~ "at least one prespecified scenario reverses the direction",
        TRUE ~ "no prespecified scenario reverses the direction"
      )
    )
}
safe_write(qba_summary, "qba_summary_by_transition.csv",
           "no transition produced a QBA scenario")

# ---------------------------------------------------------------------------
# Audit trail
# ---------------------------------------------------------------------------
run_configuration <- tibble(
  setting = c(
    "DATA_PATH", "OUTCOME_NAME", "ANALYSIS_LABEL", "TEMPORAL_DESIGN", "OUTPUT_DIR",
    "RANDOM_SEED", "NUM_TREES", "NUM_THREADS", "TUNE_PARAMETERS",
    "N_BOOT", "PROPENSITY_FLOOR", "PRIMARY_ESTIMAND",
    "RUN_HETEROGENEITY", "RUN_SUBGROUP_EVALUES", "SAVE_FORESTS", "AGE_CUT",
    "CONTRAST_DOMAINS", "EXTRA_COVARIATES", "COVARIATES", "N_CONTRASTS",
    "GRF_VERSION", "R_VERSION"
  ),
  value = c(
    DATA_PATH, OUTCOME_NAME, ANALYSIS_LABEL, TEMPORAL_DESIGN, OUTPUT_DIR,
    as.character(RANDOM_SEED), as.character(NUM_TREES), as.character(NUM_THREADS),
    TUNE_PARAMETERS, as.character(N_BOOT), as.character(PROPENSITY_FLOOR),
    PRIMARY_ESTIMAND, as.character(RUN_HETEROGENEITY), as.character(RUN_SUBGROUP_EVALUES),
    as.character(SAVE_FORESTS), as.character(AGE_CUT),
    CONTRAST_DOMAINS_RAW, paste(EXTRA_COVARIATES, collapse = ","),
    paste(covariate_names, collapse = ","), as.character(nrow(contrasts)),
    as.character(utils::packageVersion("grf")), R.version.string
  )
)
safe_write(run_configuration, "run_configuration.csv")

output_manifest <- tribble(
  ~file, ~contents,
  "cohort_summary.csv", "cohort size, event count and event rate for this run",
  "transition_registry.csv", "the 24 prespecified transitions with their clinical direction labels",
  "contrast_status.csv", "one row per transition: fitted, not_estimable, fit_failed or run_failed",
  "contrast_participant_flow.csv", "eligibility and exclusion counts per transition",
  "arm_counts_events.csv", "arm sizes, event counts and event rates per transition",
  "propensity_summary.csv", "propensity distribution per arm, and how many people hit the propensity floor",
  "ate_estimands_with_risk_scale.csv", "ATE for all five estimands, with the score-based risk difference and the consistency check",
  "evalue_summary.csv", "risk-scale results and E-values for all five estimands",
  "table_evalue_headline.csv", "one row per transition on the primary estimand; the table to put in the dissertation",
  "evalue_subgroups.csv", "E-values within the seven prespecified subgroups",
  "qba_prespecified_scenarios.csv", "the seven prespecified confounder scenarios with bootstrap intervals",
  "qba_factorial_grid.csv", "prevalence x strength grid showing the whole bias surface",
  "qba_tipping_points.csv", "closed-form RR_UY that would drive the estimate to zero, per prevalence pattern",
  "qba_summary_by_transition.csv", "reading aid: weakest reversing scenario and smallest tipping point per transition",
  "heterogeneity_test_calibration.csv", "differential forest prediction test",
  "heterogeneity_best_linear_projection.csv", "BLP on the prespecified sex, age and BMI modifiers",
  "heterogeneity_rate_autoc.csv", "RATE/AUTOC with direction-aware priorities",
  "variable_importance_descriptive.csv", "split-frequency variable importance",
  "cate_quartiles_exploratory.csv", "CATE quartile ATEs after 0.05-0.95 trimming",
  "subgroup_ate.csv", "subgroup ATEs with arm and event support",
  "subgroup_difference_tests.csv", "subgroup difference tests, normal approximation",
  "run_configuration.csv", "every setting this run used",
  "session_info.txt", "R and package versions"
)
safe_write(output_manifest, "output_manifest.csv")

writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "session_info.txt"))

if (SAVE_FORESTS) {
  saveRDS(
    list(
      all_results = all_results,
      qba_scenarios = qba_scenarios,
      qba_grid_specification = qba_grid_specification,
      tipping_point_specification = tipping_point_specification,
      covariate_names = covariate_names,
      run_configuration = run_configuration
    ),
    file.path(OUTPUT_DIR, paste0("unmeasured_confounder_", ANALYSIS_LABEL, ".rds"))
  )
}

# =============================================================================
# 11. Console summary
# =============================================================================

cat("\n================ CONTRAST STATUS ================\n")
print(status_tbl, n = Inf)

if (nrow(headline) > 0) {
  cat("\n================ E-VALUE HEADLINE (primary estimand) ================\n")
  print(headline, n = Inf)

  n_with_evalue <- sum(is.finite(headline$evalue_point))
  cat("\nTransitions with a doubly robust E-value: ", n_with_evalue, " of ", nrow(headline), "\n", sep = "")
  if ("fdr_supported" %in% names(headline)) {
    cat("Transitions surviving FDR at q < 0.05: ",
        sum(headline$fdr_supported %in% TRUE), "\n", sep = "")
  }
}

if (nrow(qba_summary) > 0) {
  cat("\n================ QBA SUMMARY BY TRANSITION ================\n")
  print(qba_summary, n = Inf)
}

log_msg("Done. All output written to ", OUTPUT_DIR)
