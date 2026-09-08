# =============================================================================
# Baseline causal model (Model 2) — one causal forest per lifestyle variable
# =============================================================================
# Plan, in plain words:
#   1. Load the WHOLE preprocessed dataset from the Python notebook.
#   2. Pick the outcome (CVD) and, for each lifestyle domain, define:
#         W (treatment) = 1 if the person is at a "healthy" baseline level, else 0
#         X (confounders) = age, sex, bmi + the OTHER lifestyle domains
#   3. Write ONE function that runs a single causal forest + reads its effect.
#   4. Loop over the six lifestyle domains and run it on each.
#   5. Collect the average effects into one table and save.
#
# This is the temporally-clean BASELINE reference model (project doc §6.2, Model 2):
# baseline lifestyle level is measured before every CVD event, so it has no
# reverse-causation / timing problem. Use it to judge whether the longitudinal
# transition results (Model 1 / Model 3) point in a believable direction.
# =============================================================================

library(grf)
library(tidyverse)

# ----------------------------------------------------------------------------
# 1. Load the whole preprocessed dataset (output of baseline_causal_preprocessing.ipynb)
# ----------------------------------------------------------------------------
df <- read_csv("/home/rmhiund/causal_analysis/grf_model/baseline/baseline_preprocessed.csv")

cat("Total observations:", nrow(df), "\n")
cat("CVD event rate:", round(mean(df$def_CVD_AFTER, na.rm = TRUE), 3), "\n\n")

# ----------------------------------------------------------------------------
# 2. Outcome, the lifestyle domains, and how "healthy" is coded for each
# ----------------------------------------------------------------------------
outcome        <- "def_CVD_AFTER"
base_covars    <- c("age", "sex", "bmi")               # always confounders
lifestyle      <- c("sleep_category_0.0", "smoking_category_0.0", "alcohol_score_0.0",
                    "diet_score_0.0", "physical_category_0.0", "mental_score_0.0")

# which baseline VALUES count as the "healthy" group (W = 1) for each domain.
# (matches the coding used in your single-variable scripts)
healthy_values <- list(
  sleep_category_0.0    = c(2),      # 2 = good sleep            (higher = healthier)
  diet_score_0.0        = c(2, 3),   # higher = better diet      (higher = healthier)
  physical_category_0.0 = c(2),      # 2 = active                (higher = healthier)
  smoking_category_0.0  = c(0),      # 0 = never smoker          (higher = worse)
  alcohol_score_0.0     = c(0, 1),   # lower = less alcohol      (higher = worse)
  mental_score_0.0      = c(0)       # 0 = low mental burden     (higher = worse)
)

Y <- df[[outcome]]

# Helper 1 — ATE inside one subgroup, returns NA if the subgroup is too small.
subgroup_ate <- function(forest, in_group, W_sub, label,
                         min_n = 50, min_each = 20) {
  n1 <- sum(in_group & W_sub == 1)
  n0 <- sum(in_group & W_sub == 0)
  if (sum(in_group) < min_n || n1 < min_each || n0 < min_each) {
    return(tibble(subgroup = label, n = sum(in_group),
                  ate = NA_real_, se = NA_real_, low = NA_real_, high = NA_real_))
  }
  est <- tryCatch(average_treatment_effect(forest, subset = in_group),
                  error = function(e) c(estimate = NA, std.err = NA))
  tibble(subgroup = label, n = sum(in_group),
         ate = unname(est["estimate"]), se = unname(est["std.err"]),
         low = unname(est["estimate"]) - 1.96 * unname(est["std.err"]),
         high = unname(est["estimate"]) + 1.96 * unname(est["std.err"]))
}
# Helper 2 — test whether two (disjoint) subgroups differ.
diff_test <- function(g1, g0, name) {
  d  <- g1$ate - g0$ate
  se <- sqrt(g1$se^2 + g0$se^2)
  tibble(contrast = name, diff = d, se = se, z = d / se, p = 2 * pnorm(-abs(d / se)))
}

# ----------------------------------------------------------------------------
# 3. One function = one causal forest for a given lifestyle domain
# ----------------------------------------------------------------------------
run_forest <- function(domain) {

  # Step 1: treatment W = 1 if baseline level is "healthy" for this domain
  W_all <- as.integer(df[[domain]] %in% healthy_values[[domain]])

  # Step 2: confounders X = age, sex, bmi + every OTHER lifestyle domain
  x_cols <- c(base_covars, setdiff(lifestyle, domain))
  X_all  <- as.matrix(df[, x_cols])

  # Step 3: keep complete rows only
  rows  <- complete.cases(X_all) & !is.na(Y) & !is.na(W_all)
  X_sub <- X_all[rows, , drop = FALSE]
  Y_sub <- Y[rows]
  W_sub <- W_all[rows]

  cat("\n----", domain, "(healthy vs less-healthy baseline) ----\n")
  cat("healthy (W=1):", sum(W_sub == 1), " less-healthy (W=0):", sum(W_sub == 0),
      " events:", sum(Y_sub), "\n")

  n_healthy   <- sum(W_sub == 1); n_less <- sum(W_sub == 0)
  ev_healthy  <- sum(Y_sub[W_sub == 1]); ev_less <- sum(Y_sub[W_sub == 0])
  rate_healthy<- mean(Y_sub[W_sub == 1]); rate_less <- mean(Y_sub[W_sub == 0])

  # Step 4: fit the causal forest
  forest <- causal_forest(X_sub, Y_sub, W_sub,
                          num.trees = 5000, seed = 42,
                          tune.parameters = "all")

  # Step 5: overlap — drop people with propensity too close to 0 or 1
  prop <- forest$W.hat                       # P(healthy | covariates) = propensity
  keep <- prop > 0.05 & prop < 0.95
  cat("propensity range:", round(min(prop), 3), "-", round(max(prop), 3),
      " | kept for overlap:", sum(keep), "(", round(mean(keep) * 100, 1), "%)\n")

  # Step 6: the AVERAGE treatment effect (healthy vs less-healthy), absolute-risk scale
  ate  <- average_treatment_effect(forest, subset = keep)
  low  <- ate["estimate"] - 1.96 * ate["std.err"]
  high <- ate["estimate"] + 1.96 * ate["std.err"]
  cat("ATE:", round(ate["estimate"], 4),
      " 95% CI [", round(low, 4), ",", round(high, 4), "]\n")

  # Step 6b: is the forest well-calibrated / is there real heterogeneity?
  cat("\nCalibration test (omnibus):\n")
  print(tryCatch(test_calibration(forest), error = function(e) "n/a"))

  # ==========================================================================
  # Step 7: CATE — FOUR ways to study how the effect VARIES between people
  # ==========================================================================
  tau     <- predict(forest, estimate.variance = TRUE)
  cate    <- as.numeric(tau$predictions)
  cate_se <- as.numeric(sqrt(tau$variance.estimates))

  # Way 1 — Variable importance
  importance <- tibble(variable   = colnames(X_sub),
                       importance = as.numeric(variable_importance(forest))) |>
    arrange(desc(importance))
  # Way 2 — Best Linear Projection
  blp  <- best_linear_projection(forest, X_sub)
  # Way 3 — RATE (AUTOC)
  rate <- rank_average_treatment_effect(forest, cate, subset = keep)
  # Way 4 — CATE by quartile
  q <- ntile(cate, 4)
  cate_quartiles <- map_dfr(1:4, function(k) {
    idx <- keep & q == k
    if (sum(idx) < 50)
      return(tibble(quartile = k, n = sum(idx),
                    ate = NA_real_, se = NA_real_, low = NA_real_, high = NA_real_))
    est <- average_treatment_effect(forest, subset = idx)
    tibble(quartile = k, n = sum(idx),
           ate = est["estimate"], se = est["std.err"],
           low = est["estimate"] - 1.96 * est["std.err"],
           high = est["estimate"] + 1.96 * est["std.err"])
  })

  cat("\nWay 1 - Variable importance:\n"); print(importance)
  cat("\nWay 2 - BLP:\n");                 print(blp)
  cat("\nWay 3 - RATE:\n");                print(rate)
  cat("\nWay 4 - CATE quartiles:\n");      print(cate_quartiles)

  # Step 7b: subgroup ATEs by sex and age (reuse the SAME forest, no refit)
  age_sub <- X_sub[, "age"]; sex_sub <- X_sub[, "sex"]; age_cut <- 60
  sg_female <- subgroup_ate(forest, keep & sex_sub == 0,       W_sub, "female (sex=0)")
  sg_male   <- subgroup_ate(forest, keep & sex_sub == 1,       W_sub, "male (sex=1)")
  sg_young  <- subgroup_ate(forest, keep & age_sub <  age_cut, W_sub, sprintf("age < %g",  age_cut))
  sg_old    <- subgroup_ate(forest, keep & age_sub >= age_cut, W_sub, sprintf("age >= %g", age_cut))
  subgroups <- bind_rows(sg_female, sg_male, sg_young, sg_old)
  het_test  <- bind_rows(
    diff_test(sg_male, sg_female, "sex: male - female"),
    diff_test(sg_old,  sg_young,  sprintf("age: >=%g - <%g", age_cut, age_cut)))
  cat("\nSubgroup ATEs (sex, age):\n"); print(subgroups)
  cat("\nSubgroup difference tests:\n"); print(het_test)

  # Step 8: hand everything back
  list(
    summary = tibble(
      domain            = domain,
      n_total           = length(Y_sub),
      n_healthy         = n_healthy,
      n_less_healthy    = n_less,
      events_total      = sum(Y_sub),
      events_healthy    = ev_healthy,
      events_less       = ev_less,
      rate_healthy      = rate_healthy,
      rate_less_healthy = rate_less,
      n_overlap         = sum(keep),
      pct_overlap       = mean(keep),
      ate               = ate["estimate"],
      se                = ate["std.err"],
      ate_low           = low,
      ate_high          = high),
    forest = forest, importance = importance, blp = blp, rate = rate,
    cate_quartiles = cate_quartiles, subgroups = subgroups, het_test = het_test,
    cate = cate, cate_se = cate_se)
}

# ----------------------------------------------------------------------------
# 4. Run one forest per lifestyle domain
# ----------------------------------------------------------------------------
all_results <- list()
for (dom in lifestyle) {
  all_results[[dom]] <- run_forest(dom)
}

# ----------------------------------------------------------------------------
# 5. Collect the average effects into one table
# ----------------------------------------------------------------------------
summary_tbl <- map_dfr(all_results, ~ .x$summary)
cat("\n================ BASELINE ATE SUMMARY (Model 2) ================\n")
print(summary_tbl, n = Inf)

subgroup_tbl <- map_dfr(names(all_results),
                        ~ all_results[[.x]]$subgroups |> mutate(domain = .x))
het_tbl      <- map_dfr(names(all_results),
                        ~ all_results[[.x]]$het_test  |> mutate(domain = .x))
# --- Variable importance (Way 1): one row per covariate per domain ---
importance_tbl <- map_dfr(names(all_results),
  ~ all_results[[.x]]$importance |> mutate(domain = .x))

# --- Best Linear Projection (Way 2): tidy the coefficient table per domain ---
blp_tbl <- map_dfr(names(all_results), function(dom) {
  b <- all_results[[dom]]$blp
  tibble(
    domain    = dom,
    term      = rownames(b),
    estimate  = b[, "Estimate"],
    std_error = b[, "Std. Error"],
    t_value   = b[, "t value"],
    p_value   = b[, "Pr(>|t|)"]
  )
})

# --- RATE / AUTOC (Way 3): one row per domain ---
rate_tbl <- map_dfr(names(all_results), function(dom) {
  r <- all_results[[dom]]$rate
  tibble(
    domain   = dom,
    target   = r$target,
    estimate = r$estimate,
    std_err  = r$std.err
  )
})

# --- CATE by quartile (Way 4): four rows per domain ---
cate_quartile_tbl <- map_dfr(names(all_results),
  ~ all_results[[.x]]$cate_quartiles |> mutate(domain = .x))

# --- Per-person CATE (every individual's estimated effect + SE) ---
cate_person_tbl <- map_dfr(names(all_results), function(dom) {
  res <- all_results[[dom]]
  tibble(domain = dom, cate = res$cate, cate_se = res$cate_se)
})

cat("\n============ VARIABLE IMPORTANCE (all domains) ============\n"); print(importance_tbl, n = Inf)
cat("\n============ BLP (all domains) ============\n");                 print(blp_tbl, n = Inf)
cat("\n============ RATE / AUTOC (all domains) ============\n");         print(rate_tbl, n = Inf)
cat("\n============ CATE QUARTILES (all domains) ============\n");       print(cate_quartile_tbl, n = Inf)
cat("\n============ SUBGROUP ATE (age/sex) ============\n"); print(subgroup_tbl, n = Inf)
cat("\n============ SUBGROUP DIFFERENCE TESTS ============\n"); print(het_tbl, n = Inf)

# ----------------------------------------------------------------------------
# 6. Save everything
# ----------------------------------------------------------------------------
outdir <- "results_baseline_causal"
dir.create(outdir, showWarnings = FALSE)

write_csv(summary_tbl,        file.path(outdir, "ate_summary.csv"))         # ATE per domain
write_csv(subgroup_tbl,       file.path(outdir, "subgroup_ate.csv"))        # ATE by sex/age
write_csv(het_tbl,            file.path(outdir, "subgroup_difference_tests.csv"))
write_csv(importance_tbl,     file.path(outdir, "variable_importance.csv")) # Way 1
write_csv(blp_tbl,            file.path(outdir, "blp.csv"))                  # Way 2
write_csv(rate_tbl,           file.path(outdir, "rate_autoc.csv"))          # Way 3
write_csv(cate_quartile_tbl,  file.path(outdir, "cate_quartiles.csv"))      # Way 4
write_csv(cate_person_tbl,    file.path(outdir, "cate_per_person.csv"))     # per-person CATE

# keep the full objects too (forests etc.) for later re-analysis
saveRDS(list(all_results = all_results, summary_tbl = summary_tbl,
             base_covars = base_covars, lifestyle = lifestyle,
             healthy_values = healthy_values),
        file.path(outdir, "baseline_causal_forests.rds"))

cat("\nDone. All metric tables saved to ./", outdir, "/ :\n", sep = "")
for (f in list.files(outdir)) cat("  -", f, "\n")
