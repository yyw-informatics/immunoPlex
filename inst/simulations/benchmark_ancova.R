#!/usr/bin/env Rscript
# benchmark_ancova.R — Benchmark 3: ANCOVA Validation (v3)
#
# ADEMP Structure:
#   Aims:    Does ancova_one()/ancova_fit() correctly recover group effects in
#            the presence of genuine confounders, nonlinearity, outliers, and
#            heteroscedasticity? Does panel-level FDR control hold after BH
#            correction across a 10-analyte panel?
#   DGP:     simulate_ancova_data() with factorial grid over sample size
#            (n/grp 20, 50, 100), effect size (0, 0.3, 0.5, 1.0),
#            confounding type, outlier rate, heteroscedasticity, and outlier
#            direction (symmetric vs. asymmetric upward). Panel sub-grid adds
#            n_analytes=10 with 3 signal analytes.
#   Estimand: Group effect on log-concentration (adjusted for covariates).
#             Partial omega-squared for effect size recovery.
#             Panel FDR, FWER, and sensitivity (multi-analyte).
#   Methods: (1) ancova_one + outlier removal (3x IQR), (2) ancova_one raw,
#            (3) unadjusted t-test/ANOVA, (4) standard parametric ANCOVA,
#            (5) ancova_one + tight fence (1.5x IQR), (6) robust ANCOVA (rlm).
#   Performance: bias, coverage, power, type I error, omega-sq recovery,
#                panel FDR, panel sensitivity, per-analyte effect recovery.
#
# v3 changes (backbone harmonization + multi-analyte panel):
#   - Aligned n_per_group to backbone: c(20, 50, 100) [was c(20, 40, 80)]
#   - Aligned effect_size to backbone: c(0, 0.3, 0.5, 1.0) [was c(0, 0.3, 0.6, 1.0)]
#   - Added 18-cell panel sub-grid: n_analytes=10, n_signal=3
#   - Panel metrics: FDR, FWER, sensitivity via BH correction
#   - Three new panel figures: FDR, sensitivity, effect recovery
#   - Uses "ancova_v3" cache prefix (independent of v2 cache)
#
# Usage:
#   Rscript benchmark_ancova.R [--n_reps N] [--n_cores N] [--cache_dir DIR]
#
# Results are cached as .rds files; re-running tops up to the target replication
# count without re-running existing replications.

# ---- Setup ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# Load immunoPlex: prefer devtools::load_all() for development, fall back to library()
if (requireNamespace("devtools", quietly = TRUE) &&
    file.exists(file.path(getwd(), "DESCRIPTION"))) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(immunoPlex)
}

# Source benchmark helpers (relative to this script's location)
script_dir <- if (exists("script_dir")) script_dir else {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("--file=", "", file_arg)))
  } else {
    "inst/simulations"
  }
}
source(file.path(script_dir, "helpers_benchmark.R"))


# ---- Configuration ----------------------------------------------------------

cli_args <- commandArgs(trailingOnly = TRUE)
parse_cli_arg <- function(flag, default) {
  idx <- which(cli_args == flag)
  if (length(idx) > 0 && idx < length(cli_args)) return(cli_args[idx + 1])
  default
}

N_REPS     <- as.integer(parse_cli_arg("--n_reps",    "500"))
N_CORES    <- as.integer(parse_cli_arg("--n_cores",   "1"))
CACHE_DIR  <- parse_cli_arg("--cache_dir", file.path(script_dir, "cache"))
EXTRA_N    <- parse_cli_arg("--extra_n", NULL)
PHASE_B    <- parse_cli_arg("--phase_b", NULL)
BASE_SEED  <- 20240301L
BENCHMARK  <- "ancova_v3"
VERBOSE    <- TRUE

# Extra-n supplement mode
EXTRA_SIZES  <- if (!is.null(EXTRA_N)) as.integer(strsplit(EXTRA_N, ",")[[1]]) else NULL
EXTRA_OFFSET <- 10000L

# Phase B offset for sub-grid cell IDs
PHASE_B_OFFSET <- 20000L

if (!dir.exists(CACHE_DIR)) dir.create(CACHE_DIR, recursive = TRUE)


# ---- DGP Grid ---------------------------------------------------------------

# Main factorial grid:
#   n_per_group (20, 50, 100) x effect_size (0, 0.3, 0.5, 1.0)
#   x confounder_type (none, balanced, confounding)
#   x outlier_rate (0, 0.05, 0.10) x heteroscedastic (FALSE, TRUE)
# = 3 x 4 x 3 x 3 x 2 = 216 cells

dgp_grid <- make_dgp_grid(
  n_per_group     = if (!is.null(EXTRA_SIZES)) EXTRA_SIZES else c(20L, 50L, 100L),
  effect_size     = c(0, 0.3, 0.5, 1.0),
  confounder_type = c("none", "balanced", "confounding"),
  outlier_rate    = c(0, 0.05, 0.10),
  heteroscedastic = c(FALSE, TRUE)
)
dgp_grid$nonlinear_cov <- FALSE
dgp_grid$outlier_direction <- "symmetric"

# Additional sub-analysis: nonlinear covariate effect (quadratic age)
# Tests robustness of rank-based ANCOVA to model misspecification
nonlinear_grid <- make_dgp_grid(
  n_per_group     = if (!is.null(EXTRA_SIZES)) EXTRA_SIZES else c(50L, 100L),
  effect_size     = c(0, 0.5, 1.0),
  confounder_type = "confounding",
  outlier_rate    = 0,
  heteroscedastic = FALSE
)
nonlinear_grid$nonlinear_cov <- TRUE
nonlinear_grid$outlier_direction <- "symmetric"
nonlinear_grid$cell_id <- nonlinear_grid$cell_id + nrow(dgp_grid)

# Additional sub-analysis: asymmetric (upward-only) outlier injection
# Tests whether outlier removal and robust methods handle directional
# contamination better than symmetric injection (which does not bias
# the group effect estimate).
asymmetric_grid <- make_dgp_grid(
  n_per_group     = if (!is.null(EXTRA_SIZES)) EXTRA_SIZES else c(50L, 100L),
  effect_size     = c(0, 0.3, 1.0),
  confounder_type = "confounding",
  outlier_rate    = c(0.05, 0.10),
  heteroscedastic = FALSE
)
asymmetric_grid$nonlinear_cov <- FALSE
asymmetric_grid$outlier_direction <- "upward"
asymmetric_grid$cell_id <- asymmetric_grid$cell_id + nrow(dgp_grid) + nrow(nonlinear_grid)

dgp_grid <- bind_rows(dgp_grid, nonlinear_grid, asymmetric_grid)
if (!is.null(EXTRA_SIZES)) dgp_grid$cell_id <- dgp_grid$cell_id + EXTRA_OFFSET

if (VERBOSE) {
  message("=== Benchmark 3: ANCOVA Validation (v3) ===")
  message("DGP grid: ", nrow(dgp_grid), " cells x ", N_REPS, " reps x 6 methods = ",
          nrow(dgp_grid) * N_REPS * 6L, " total fits")
  message("  Main grid: 216 cells [n/grp={20,50,100} x eff={0,0.3,0.5,1.0} x ...]")
  message("  Nonlinear sub-grid: ", nrow(nonlinear_grid), " cells")
  message("  Asymmetric outlier sub-grid: ", nrow(asymmetric_grid), " cells")
  message("Cache dir: ", CACHE_DIR)
}


# ---- Method Fitting Functions -----------------------------------------------

#' Fit ancova_one and extract log-scale results
#'
#' Uses the rank-based p-value for inference and the log-validation model
#' for bias/coverage assessment on the original scale.
#'
#' @param d           data.frame of simulated data.
#' @param extra_covs  Character vector of additional covariates (beyond "age"), or NULL.
#' @param outlier_removal Logical.
#' @param method_name Character label.
#' @return One-row data.frame with standard benchmark columns.
fit_ancova_method <- function(d, extra_covs, outlier_removal, method_name) {
  fail_row <- data.frame(
    method = method_name, estimate = NA_real_, se = NA_real_,
    ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
    converged = FALSE, omega_sq = NA_real_, n_removed = 0L,
    stringsAsFactors = FALSE
  )

  tryCatch({
    result <- ancova_one(
      data = d, outcome = "value", group = "group", covariate = "age",
      covariates = extra_covs,
      outlier_removal = outlier_removal,
      log_validate = TRUE
    )

    if (is.null(result$model)) return(fail_row)

    # Log-scale estimate from validation model
    est <- if (!is.null(result$log_effect)) result$log_effect else NA_real_

    log_ci <- c(NA_real_, NA_real_)
    log_se <- NA_real_

    if (!is.null(result$log_model)) {
      group_term <- paste0("group", result$group_levels[2])
      tryCatch({
        ci <- confint(result$log_model, level = 0.95)
        if (group_term %in% rownames(ci)) {
          log_ci <- as.numeric(ci[group_term, ])
        }
        cc <- summary(result$log_model)$coefficients
        if (group_term %in% rownames(cc)) {
          log_se <- as.numeric(cc[group_term, "Std. Error"])
        }
      }, error = function(e) NULL)
    }

    data.frame(
      method    = method_name,
      estimate  = est,
      se        = log_se,
      ci_lo     = log_ci[1],
      ci_hi     = log_ci[2],
      p_value   = result$group_p,   # rank-based p-value
      converged = TRUE,
      omega_sq  = result$omega_sq_partial,
      n_removed = as.integer(result$n_removed),
      stringsAsFactors = FALSE
    )
  }, error = function(e) fail_row)
}


#' Fit ancova_one with a tight IQR fence (1.5x multiplier)
#'
#' Identical to fit_ancova_method with outlier_removal = TRUE, except
#' the IQR multiplier is set to 1.5 (vs. the default 3).
#'
#' @param d           data.frame of simulated data.
#' @param extra_covs  Character vector of additional covariates (beyond "age"), or NULL.
#' @param method_name Character label.
#' @return One-row data.frame with standard benchmark columns.
fit_ancova_tight_fence <- function(d, extra_covs, method_name) {
  fail_row <- data.frame(
    method = method_name, estimate = NA_real_, se = NA_real_,
    ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
    converged = FALSE, omega_sq = NA_real_, n_removed = 0L,
    stringsAsFactors = FALSE
  )

  tryCatch({
    result <- ancova_one(
      data = d, outcome = "value", group = "group", covariate = "age",
      covariates = extra_covs,
      outlier_removal = TRUE,
      iqr_multiplier = 1.5,
      log_validate = TRUE
    )

    if (is.null(result$model)) return(fail_row)

    est <- if (!is.null(result$log_effect)) result$log_effect else NA_real_

    log_ci <- c(NA_real_, NA_real_)
    log_se <- NA_real_

    if (!is.null(result$log_model)) {
      group_term <- paste0("group", result$group_levels[2])
      tryCatch({
        ci <- confint(result$log_model, level = 0.95)
        if (group_term %in% rownames(ci)) {
          log_ci <- as.numeric(ci[group_term, ])
        }
        cc <- summary(result$log_model)$coefficients
        if (group_term %in% rownames(cc)) {
          log_se <- as.numeric(cc[group_term, "Std. Error"])
        }
      }, error = function(e) NULL)
    }

    data.frame(
      method    = method_name,
      estimate  = est,
      se        = log_se,
      ci_lo     = log_ci[1],
      ci_hi     = log_ci[2],
      p_value   = result$group_p,
      converged = TRUE,
      omega_sq  = result$omega_sq_partial,
      n_removed = as.integer(result$n_removed),
      stringsAsFactors = FALSE
    )
  }, error = function(e) fail_row)
}


#' Fit standard lm (unadjusted ANOVA or parametric ANCOVA)
#'
#' @param d          data.frame of simulated data.
#' @param covariates Character vector of covariate names to include, or NULL.
#' @param method_name Character label.
#' @return One-row data.frame with standard benchmark columns.
fit_lm_method <- function(d, covariates, method_name) {
  fail_row <- data.frame(
    method = method_name, estimate = NA_real_, se = NA_real_,
    ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
    converged = FALSE, omega_sq = NA_real_, n_removed = 0L,
    stringsAsFactors = FALSE
  )

  tryCatch({
    # Build formula
    rhs <- if (length(covariates) > 0) {
      paste(c("group", covariates), collapse = " + ")
    } else {
      "group"
    }
    fml <- as.formula(paste("value ~", rhs))

    fit <- lm(fml, data = d)
    cc <- summary(fit)$coefficients

    # Find group coefficient
    group_levels <- levels(d$group)
    group_term <- paste0("group", group_levels[2])

    if (!(group_term %in% rownames(cc))) return(fail_row)

    est  <- as.numeric(cc[group_term, "Estimate"])
    se   <- as.numeric(cc[group_term, "Std. Error"])
    pval <- as.numeric(cc[group_term, "Pr(>|t|)"])

    ci <- confint(fit, level = 0.95)
    ci_vals <- as.numeric(ci[group_term, ])

    # Partial omega-squared via effectsize
    # Note: for one-way designs, effectsize returns Omega2 instead of Omega2_partial
    omega_sq <- tryCatch({
      es <- effectsize::omega_squared(fit, partial = TRUE)
      grp_row <- which(es$Parameter == "group")
      omega_col <- if ("Omega2_partial" %in% names(es)) "Omega2_partial" else "Omega2"
      if (length(grp_row) > 0) es[[omega_col]][grp_row] else NA_real_
    }, error = function(e) NA_real_)

    data.frame(
      method    = method_name,
      estimate  = est,
      se        = se,
      ci_lo     = ci_vals[1],
      ci_hi     = ci_vals[2],
      p_value   = pval,
      converged = TRUE,
      omega_sq  = omega_sq,
      n_removed = 0L,
      stringsAsFactors = FALSE
    )
  }, error = function(e) fail_row)
}


#' Fit robust ANCOVA via MASS::rlm (M-estimator)
#'
#' Uses Huber's M-estimator for robustness to outliers while still adjusting
#' for covariates.  Falls back to OLS if MASS is not available.
#'
#' @param d          data.frame of simulated data.
#' @param covariates Character vector of covariate names to include, or NULL.
#' @param method_name Character label.
#' @return One-row data.frame with standard benchmark columns.
fit_rlm_method <- function(d, covariates, method_name) {
  fail_row <- data.frame(
    method = method_name, estimate = NA_real_, se = NA_real_,
    ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
    converged = FALSE, omega_sq = NA_real_, n_removed = 0L,
    stringsAsFactors = FALSE
  )

  if (!requireNamespace("MASS", quietly = TRUE)) return(fail_row)

  tryCatch({
    # Build formula
    rhs <- if (length(covariates) > 0) {
      paste(c("group", covariates), collapse = " + ")
    } else {
      "group"
    }
    fml <- as.formula(paste("value ~", rhs))

    fit <- MASS::rlm(fml, data = d, method = "M")
    cc <- summary(fit)$coefficients

    # Find group coefficient
    group_levels <- levels(d$group)
    group_term <- paste0("group", group_levels[2])

    if (!(group_term %in% rownames(cc))) return(fail_row)

    est  <- as.numeric(cc[group_term, "Value"])
    se   <- as.numeric(cc[group_term, "Std. Error"])

    # rlm does not provide p-values or df.residual; compute manually
    df_resid <- nrow(d) - length(coef(fit))
    t_stat   <- est / se
    pval     <- 2 * pt(abs(t_stat), df = df_resid, lower.tail = FALSE)

    # Wald CI
    ci_lo <- est - qt(0.975, df = df_resid) * se
    ci_hi <- est + qt(0.975, df = df_resid) * se

    # Partial omega-squared from an equivalent OLS fit (for comparability)
    omega_sq <- tryCatch({
      ols_fit <- lm(fml, data = d)
      es <- effectsize::omega_squared(ols_fit, partial = TRUE)
      grp_row <- which(es$Parameter == "group")
      omega_col <- if ("Omega2_partial" %in% names(es)) "Omega2_partial" else "Omega2"
      if (length(grp_row) > 0) es[[omega_col]][grp_row] else NA_real_
    }, error = function(e) NA_real_)

    data.frame(
      method    = method_name,
      estimate  = est,
      se        = se,
      ci_lo     = ci_lo,
      ci_hi     = ci_hi,
      p_value   = pval,
      converged = TRUE,
      omega_sq  = omega_sq,
      n_removed = 0L,
      stringsAsFactors = FALSE
    )
  }, error = function(e) fail_row)
}


# ---- Asymmetric Outlier Injection -------------------------------------------

#' Inject outliers in a single direction (upward or downward)
#'
#' Replicates the symmetric injection logic of simulate_ancova_data() but
#' forces all shifts to the same sign.
#'
#' @param d              data.frame with a "value" column.
#' @param outlier_rate   Fraction of observations to shift.
#' @param direction      "upward" (+3 SD) or "downward" (-3 SD).
#' @param seed           RNG seed for reproducibility.
#' @return Modified data.frame with outliers injected.
inject_asymmetric_outliers <- function(d, outlier_rate, direction, seed) {
  if (outlier_rate <= 0) return(d)
  n_obs <- nrow(d)
  n_outlier <- round(n_obs * outlier_rate)
  if (n_outlier == 0) return(d)

  set.seed(seed + 2000L)
  outlier_idx <- sort(sample.int(n_obs, n_outlier))
  shift_sign <- if (direction == "upward") 1 else -1
  d$value[outlier_idx] <- d$value[outlier_idx] + shift_sign * 3 * 1.0
  d
}


# ---- Single-Replication Function --------------------------------------------

#' Run one replication for a single DGP cell
#'
#' Generates data via simulate_ancova_data(), applies all 6 methods, and
#' returns a data.frame of results. For single-analyte mode (n_analytes=1),
#' returns 6 rows (one per method). For multi-analyte panel mode, returns
#' 6 x n_analytes rows with additional analyte, is_signal columns.
#'
#' @param seed            RNG seed.
#' @param n_per_group     Sample size per group.
#' @param effect_size     True group effect (log-scale); 0 = null.
#' @param confounder_type "none", "balanced", or "confounding".
#' @param outlier_rate    Fraction of outliers injected (0, 0.05, 0.10).
#' @param heteroscedastic Logical. Heteroscedastic residuals?
#' @param nonlinear_cov   Logical. Quadratic covariate effect?
#' @param outlier_direction Character. "symmetric" (default, +/- 3 SD),
#'   "upward" (+3 SD only), or "downward" (-3 SD only).
#' @param n_analytes      Integer. Number of analytes (1 = single-analyte,
#'   >1 = multi-analyte panel mode). Default 1L.
#' @param n_signal        Integer. Number of signal analytes when n_analytes > 1.
#'   Default 3L.
#' @param covariate_group_cor Numeric or NULL. Explicit covariate-group
#'   correlation passed to simulate_ancova_data(). When NULL (default), the
#'   correlation is derived from confounder_type (0.4 for "confounding", 0
#'   otherwise).
#' @return data.frame with columns: method, estimate, se, ci_lo, ci_hi,
#'   p_value, converged, omega_sq, n_removed, true_effect. Panel mode adds:
#'   analyte, is_signal.
run_one_rep <- function(seed, n_per_group, effect_size, confounder_type,
                        outlier_rate, heteroscedastic, nonlinear_cov,
                        outlier_direction = "symmetric",
                        n_analytes = 1L, n_signal = 3L,
                        covariate_group_cor = NULL) {

  # Map confounder type to DGP correlation (override with explicit value)
  cov_group_cor <- if (!is.null(covariate_group_cor)) {
    covariate_group_cor
  } else if (confounder_type == "confounding") {
    0.4
  } else {
    0
  }

  # For asymmetric outliers, generate clean data then inject post-hoc
  sim_outlier_rate <- if (outlier_direction == "symmetric") outlier_rate else 0

  # --- Multi-analyte panel mode ---
  if (n_analytes > 1L) {
    sim <- simulate_ancova_data(
      n_subjects          = n_per_group * 2L,
      n_analytes          = n_analytes,
      signal_analytes     = if (effect_size != 0) n_signal else 0L,
      group_effects       = effect_size,
      covariate_group_cor = cov_group_cor,
      heteroscedastic     = heteroscedastic,
      outlier_rate        = sim_outlier_rate,
      nonlinear_cov       = nonlinear_cov,
      seed                = seed
    )

    d_full <- as.data.frame(sim$data)
    dgp <- sim$meta$dgp_params
    analyte_names <- dgp$analyte_names
    true_effects  <- dgp$group_effects   # per-analyte vector
    signal_mask   <- dgp$signal_analytes  # per-analyte logical

    extra_covs <- if (confounder_type != "none") "bmi" else NULL
    adj_covs   <- if (confounder_type != "none") c("age", "bmi") else "age"

    all_results <- vector("list", length(analyte_names))

    for (k in seq_along(analyte_names)) {
      a <- analyte_names[k]
      d <- d_full[d_full$cytokine == a, , drop = FALSE]

      if (!is.factor(d$group)) d$group <- factor(d$group)
      if (levels(d$group)[1] != "control") {
        d$group <- relevel(d$group, ref = "control")
      }

      m1 <- fit_ancova_method(d, extra_covs, outlier_removal = TRUE,
                              method_name = "ANCOVA + outlier")
      m2 <- fit_ancova_method(d, extra_covs, outlier_removal = FALSE,
                              method_name = "ANCOVA raw")
      m3 <- fit_lm_method(d, covariates = NULL, method_name = "t-test/ANOVA")
      m4 <- fit_lm_method(d, covariates = adj_covs,
                           method_name = "Parametric ANCOVA")
      m5 <- fit_ancova_tight_fence(d, extra_covs,
                                    method_name = "ANCOVA + tight fence")
      m6 <- fit_rlm_method(d, covariates = adj_covs,
                            method_name = "Robust ANCOVA")

      results_k <- rbind(m1, m2, m3, m4, m5, m6)
      results_k$true_effect <- true_effects[k]
      results_k$analyte     <- a
      results_k$is_signal   <- signal_mask[k]
      all_results[[k]] <- results_k
    }

    return(do.call(rbind, all_results))
  }

  # --- Single-analyte mode (original) ---

  # Simulate data
  sim <- simulate_ancova_data(
    n_subjects          = n_per_group * 2L,
    n_analytes          = 1L,
    signal_analytes     = if (effect_size != 0) 1L else 0L,
    group_effects       = effect_size,
    covariate_group_cor = cov_group_cor,
    heteroscedastic     = heteroscedastic,
    outlier_rate        = sim_outlier_rate,
    nonlinear_cov       = nonlinear_cov,
    seed                = seed
  )

  d <- as.data.frame(sim$data)

  # Inject asymmetric outliers post-hoc if requested
  if (outlier_direction != "symmetric" && outlier_rate > 0) {
    d <- inject_asymmetric_outliers(d, outlier_rate, outlier_direction, seed)
  }

  # Ensure group is factor with control as reference
  if (!is.factor(d$group)) d$group <- factor(d$group)
  if (levels(d$group)[1] != "control") {
    d$group <- relevel(d$group, ref = "control")
  }

  true_effect <- effect_size

  # Covariates for adjusted methods depend on confounder_type:
  #   "none":        ANCOVA uses age only; parametric uses age only
  #   "balanced":    ANCOVA uses age + bmi; parametric uses age + bmi
  #   "confounding": ANCOVA uses age + bmi; parametric uses age + bmi
  extra_covs <- if (confounder_type != "none") "bmi" else NULL
  adj_covs   <- if (confounder_type != "none") c("age", "bmi") else "age"

  # --- Method 1: ANCOVA + outlier removal (3x IQR fence) ---
  m1 <- fit_ancova_method(d, extra_covs, outlier_removal = TRUE,
                          method_name = "ANCOVA + outlier")

  # --- Method 2: ANCOVA raw (no outlier removal) ---
  m2 <- fit_ancova_method(d, extra_covs, outlier_removal = FALSE,
                          method_name = "ANCOVA raw")

  # --- Method 3: Unadjusted t-test/ANOVA (ignoring covariates) ---
  m3 <- fit_lm_method(d, covariates = NULL, method_name = "t-test/ANOVA")

  # --- Method 4: Standard parametric ANCOVA ---
  m4 <- fit_lm_method(d, covariates = adj_covs,
                       method_name = "Parametric ANCOVA")

  # --- Method 5: ANCOVA + tight fence (1.5x IQR) ---
  m5 <- fit_ancova_tight_fence(d, extra_covs, method_name = "ANCOVA + tight fence")

  # --- Method 6: Robust ANCOVA (MASS::rlm) ---
  m6 <- fit_rlm_method(d, covariates = adj_covs,
                        method_name = "Robust ANCOVA")

  # Combine
  results <- rbind(m1, m2, m3, m4, m5, m6)
  results$true_effect <- true_effect
  results
}


# ---- Run Benchmark ----------------------------------------------------------

run_benchmark_ancova <- function(dgp_grid, n_reps, base_seed,
                                 cache_dir, n_cores, verbose) {

  all_summaries <- list()

  for (i in seq_len(nrow(dgp_grid))) {
    row <- dgp_grid[i, ]
    cell_label <- paste0(
      "cell_", row$cell_id,
      " [n/grp=", row$n_per_group,
      "/eff=", row$effect_size,
      "/conf=", row$confounder_type,
      "/out=", row$outlier_rate,
      "/het=", row$heteroscedastic,
      if (row$nonlinear_cov) "/nonlin" else "",
      if (row$outlier_direction != "symmetric") paste0("/", row$outlier_direction) else "",
      "]"
    )

    # Build the per-replication function (closure over DGP parameters)
    rep_fn <- function(seed) {
      run_one_rep(
        seed              = seed,
        n_per_group       = row$n_per_group,
        effect_size       = row$effect_size,
        confounder_type   = row$confounder_type,
        outlier_rate      = row$outlier_rate,
        heteroscedastic   = row$heteroscedastic,
        nonlinear_cov     = row$nonlinear_cov,
        outlier_direction = row$outlier_direction
      )
    }

    cache_file <- file.path(cache_dir,
                            paste0(BENCHMARK, "_cell", row$cell_id, ".rds"))

    raw_results <- run_benchmark_cell(
      rep_fn     = rep_fn,
      n_reps     = n_reps,
      base_seed  = base_seed + (i - 1) * 10000L,
      cache_file = cache_file,
      n_cores    = n_cores,
      cell_label = cell_label,
      verbose    = verbose
    )

    # Aggregate per method within this cell
    methods <- unique(raw_results$method[!is.na(raw_results$method)])
    for (m in methods) {
      m_reps <- raw_results[!is.na(raw_results$method) & raw_results$method == m, ]
      agg <- aggregate_cell(m_reps, truth = row$effect_size)

      # Omega-squared summary
      omega_vals <- m_reps$omega_sq[is.finite(m_reps$omega_sq)]
      agg$mean_omega_sq <- if (length(omega_vals) > 0) mean(omega_vals) else NA_real_
      agg$sd_omega_sq   <- if (length(omega_vals) > 1) sd(omega_vals) else NA_real_

      # Outlier removal summary
      if ("n_removed" %in% names(m_reps)) {
        removed_vals <- m_reps$n_removed[is.finite(m_reps$n_removed)]
        agg$mean_n_removed <- if (length(removed_vals) > 0) mean(removed_vals) else NA_real_
      } else {
        agg$mean_n_removed <- NA_real_
      }

      # Attach DGP labels
      agg$method          <- m
      agg$n_per_group     <- row$n_per_group
      agg$effect_size     <- row$effect_size
      agg$confounder_type <- row$confounder_type
      agg$outlier_rate    <- row$outlier_rate
      agg$heteroscedastic   <- row$heteroscedastic
      agg$nonlinear_cov     <- row$nonlinear_cov
      agg$outlier_direction <- row$outlier_direction
      agg$cell_id           <- row$cell_id

      all_summaries[[length(all_summaries) + 1]] <- agg
    }
  }

  bind_rows(all_summaries)
}


# ---- Phase B6: Covariate-Group Correlation Sub-Grid -------------------------

if (!is.null(PHASE_B)) {

  # Phase B6 sub-grid:
  #   n_per_group = {50, 100} x effect_size = {0, 0.5, 1.0}
  #   x covariate_group_cor = {0.3, 0.6} x outlier_rate = {0, 0.05}
  #   x heteroscedastic = {FALSE, TRUE}
  # = 2 x 3 x 2 x 2 x 2 = 48 cells
  # Fixed: confounder_type = "confounding", nonlinear_cov = FALSE,
  #        outlier_direction = "symmetric"

  phaseB_base <- make_dgp_grid(
    n_per_group     = c(50L, 100L),
    effect_size     = c(0, 0.5, 1.0),
    confounder_type = "confounding",
    outlier_rate    = c(0, 0.05),
    heteroscedastic = c(FALSE, TRUE)
  )
  phaseB_base$nonlinear_cov     <- FALSE
  phaseB_base$outlier_direction <- "symmetric"

  # Expand over covariate_group_cor
  phaseB_grid <- bind_rows(
    lapply(c(0.3, 0.6), function(cgc) {
      g <- phaseB_base
      g$covariate_group_cor <- cgc
      g
    })
  )
  # Assign cell IDs with Phase B offset
  phaseB_grid$cell_id <- seq_len(nrow(phaseB_grid)) + PHASE_B_OFFSET

  if (VERBOSE) {
    message("\n=== Phase B6: Covariate-Group Correlation ===")
    message("DGP grid: ", nrow(phaseB_grid), " cells x ", N_REPS,
            " reps x 6 methods = ",
            nrow(phaseB_grid) * N_REPS * 6L, " total fits")
    message("  covariate_group_cor = {0.3, 0.6}")
    message("Cache dir: ", CACHE_DIR)
  }

  # Run Phase B6 benchmark using same loop structure as main grid
  phaseB_summaries <- list()

  for (i in seq_len(nrow(phaseB_grid))) {
    row <- phaseB_grid[i, ]
    cell_label <- paste0(
      "phaseB_cell_", row$cell_id,
      " [n/grp=", row$n_per_group,
      "/eff=", row$effect_size,
      "/cov_grp_cor=", row$covariate_group_cor,
      "/out=", row$outlier_rate,
      "/het=", row$heteroscedastic, "]"
    )

    rep_fn <- function(seed) {
      run_one_rep(
        seed                = seed,
        n_per_group         = row$n_per_group,
        effect_size         = row$effect_size,
        confounder_type     = row$confounder_type,
        outlier_rate        = row$outlier_rate,
        heteroscedastic     = row$heteroscedastic,
        nonlinear_cov       = row$nonlinear_cov,
        outlier_direction   = row$outlier_direction,
        covariate_group_cor = row$covariate_group_cor
      )
    }

    cache_file <- file.path(CACHE_DIR,
                            paste0(BENCHMARK, "_phaseB_cell", row$cell_id, ".rds"))

    raw_results <- run_benchmark_cell(
      rep_fn     = rep_fn,
      n_reps     = N_REPS,
      base_seed  = BASE_SEED + 5000000L + (i - 1) * 10000L,
      cache_file = cache_file,
      n_cores    = N_CORES,
      cell_label = cell_label,
      verbose    = VERBOSE
    )

    # Aggregate per method within this cell
    methods <- unique(raw_results$method[!is.na(raw_results$method)])
    for (m in methods) {
      m_reps <- raw_results[!is.na(raw_results$method) & raw_results$method == m, ]
      agg <- aggregate_cell(m_reps, truth = row$effect_size)

      # Omega-squared summary
      omega_vals <- m_reps$omega_sq[is.finite(m_reps$omega_sq)]
      agg$mean_omega_sq <- if (length(omega_vals) > 0) mean(omega_vals) else NA_real_
      agg$sd_omega_sq   <- if (length(omega_vals) > 1) sd(omega_vals) else NA_real_

      # Outlier removal summary
      if ("n_removed" %in% names(m_reps)) {
        removed_vals <- m_reps$n_removed[is.finite(m_reps$n_removed)]
        agg$mean_n_removed <- if (length(removed_vals) > 0) mean(removed_vals) else NA_real_
      } else {
        agg$mean_n_removed <- NA_real_
      }

      # Attach DGP labels
      agg$method              <- m
      agg$n_per_group         <- row$n_per_group
      agg$effect_size         <- row$effect_size
      agg$confounder_type     <- row$confounder_type
      agg$outlier_rate        <- row$outlier_rate
      agg$heteroscedastic     <- row$heteroscedastic
      agg$nonlinear_cov       <- row$nonlinear_cov
      agg$outlier_direction   <- row$outlier_direction
      agg$covariate_group_cor <- row$covariate_group_cor
      agg$cell_id             <- row$cell_id

      phaseB_summaries[[length(phaseB_summaries) + 1]] <- agg
    }
  }

  phaseB_summary_df <- bind_rows(phaseB_summaries)

  # Acceptance annotations
  phaseB_null <- phaseB_summary_df %>%
    filter(effect_size == 0) %>%
    annotate_acceptance(truth = 0, type1_bounds = c(0.025, 0.075))

  phaseB_alt_list <- lapply(
    unique(phaseB_summary_df$effect_size[phaseB_summary_df$effect_size != 0]),
    function(es) {
      phaseB_summary_df %>%
        filter(effect_size == es) %>%
        annotate_acceptance(truth = es,
                            bias_ratio_threshold = 0.15,
                            coverage_bounds = c(0.93, 0.97))
    }
  )
  phaseB_alt <- bind_rows(phaseB_alt_list)

  phaseB_annotated <- bind_rows(phaseB_null, phaseB_alt)

  # Save Phase B summary (merge with main summary if it exists)
  summary_file <- file.path(CACHE_DIR, paste0(BENCHMARK, "_summary.rds"))
  if (file.exists(summary_file)) {
    existing <- readRDS(summary_file)
    # Ensure covariate_group_cor column exists in both
    if (!"covariate_group_cor" %in% names(existing)) {
      existing$covariate_group_cor <- NA_real_
    }
    # Remove any previous Phase B rows and append fresh
    existing <- existing[is.na(existing$covariate_group_cor) |
                         !existing$cell_id %in% phaseB_annotated$cell_id, ]
    combined <- bind_rows(existing, phaseB_annotated)
    saveRDS(combined, summary_file)
    if (VERBOSE) message("Phase B summary merged into -> ", summary_file)
  } else {
    saveRDS(phaseB_annotated, summary_file)
    if (VERBOSE) message("Phase B summary saved -> ", summary_file)
  }

  if (VERBOSE) {
    message("\n=== Phase B6: Summary ===")
    phaseB_fmt <- format_summary_table(phaseB_annotated)
    key_b6 <- phaseB_fmt %>%
      filter(outlier_rate == 0, !heteroscedastic) %>%
      select(n_per_group, effect_size, covariate_group_cor, method,
             n_reps, bias_fmt, coverage_fmt, rejection_fmt,
             starts_with("pass_"))
    print(as.data.frame(key_b6), right = FALSE)
    message("\n=== Phase B6 complete ===")
  }

} else {

# ---- Execute ----------------------------------------------------------------

if (VERBOSE) message("\nStarting benchmark execution...")

summary_df <- run_benchmark_ancova(
  dgp_grid  = dgp_grid,
  n_reps    = N_REPS,
  base_seed = BASE_SEED,
  cache_dir = CACHE_DIR,
  n_cores   = N_CORES,
  verbose   = VERBOSE
)

# Acceptance annotations — separate null from alternative scenarios
summary_null <- summary_df %>%
  filter(effect_size == 0) %>%
  annotate_acceptance(truth = 0, type1_bounds = c(0.025, 0.075))

summary_alt <- summary_df %>%
  filter(effect_size != 0) %>%
  annotate_acceptance(truth = 1,   # non-zero: skips type1 check
                      bias_ratio_threshold = 0.15,
                      coverage_bounds = c(0.93, 0.97))

summary_annotated <- bind_rows(summary_null, summary_alt)

# Save final summary (merge with existing when running extra sizes)
summary_file <- file.path(CACHE_DIR, paste0(BENCHMARK, "_summary.rds"))
if (!is.null(EXTRA_SIZES)) {
  merge_with_existing_summary(summary_annotated, summary_file, verbose = VERBOSE)
} else {
  save_main_summary_preserving_supplements(summary_annotated, summary_file,
                                           verbose = VERBOSE)
}


# ---- Omega-Squared Recovery Analysis ----------------------------------------

# Compute correlation between estimated omega-sq and true effect size
# across main grid cells only (excluding null, nonlinear, and asymmetric)
omega_recovery <- summary_annotated %>%
  filter(effect_size > 0, !nonlinear_cov, outlier_direction == "symmetric") %>%
  group_by(method) %>%
  summarise(
    omega_effect_cor = cor(mean_omega_sq, effect_size,
                           use = "pairwise.complete.obs",
                           method = "spearman"),
    pass_omega = omega_effect_cor > 0.85,
    .groups = "drop"
  )

if (VERBOSE) {
  message("\n=== Omega-Squared Recovery (Spearman cor with effect size) ===")
  print(as.data.frame(omega_recovery), right = FALSE)
}

# Save omega recovery
saveRDS(omega_recovery, file.path(CACHE_DIR, paste0(BENCHMARK, "_omega_recovery.rds")))


# ---- Panel Sub-Grid (Multi-Analyte) -----------------------------------------

# 18 cells: n_per_group={20,50,100} x effect_size={0,0.5,1.0}
#           x confounder_type={none,confounding}
# Fixed: outlier_rate=0, heteroscedastic=FALSE (panel simplicity)
N_ANALYTES_PANEL <- 10L
N_SIGNAL_PANEL   <- 3L

panel_grid <- make_dgp_grid(
  n_per_group     = if (!is.null(EXTRA_SIZES)) EXTRA_SIZES else c(20L, 50L, 100L),
  effect_size     = c(0, 0.5, 1.0),
  confounder_type = c("none", "confounding")
)
if (!is.null(EXTRA_SIZES)) panel_grid$cell_id <- panel_grid$cell_id + EXTRA_OFFSET

if (VERBOSE) {
  message("\n=== Panel Extension: ", nrow(panel_grid), " cells x ", N_REPS,
          " reps x 6 methods x ", N_ANALYTES_PANEL, " analytes ===")
}

panel_raw_list <- list()

for (i in seq_len(nrow(panel_grid))) {
  row <- panel_grid[i, ]
  cell_label <- paste0(
    "panel_cell_", row$cell_id,
    " [n/grp=", row$n_per_group,
    "/eff=", row$effect_size,
    "/conf=", row$confounder_type, "]"
  )

  rep_fn <- function(seed) {
    run_one_rep(
      seed              = seed,
      n_per_group       = row$n_per_group,
      effect_size       = row$effect_size,
      confounder_type   = row$confounder_type,
      outlier_rate      = 0,
      heteroscedastic   = FALSE,
      nonlinear_cov     = FALSE,
      outlier_direction = "symmetric",
      n_analytes        = N_ANALYTES_PANEL,
      n_signal          = N_SIGNAL_PANEL
    )
  }

  cache_file <- file.path(CACHE_DIR,
                          paste0(BENCHMARK, "_panel_cell", row$cell_id, ".rds"))

  raw_results <- run_benchmark_cell(
    rep_fn     = rep_fn,
    n_reps     = N_REPS,
    base_seed  = BASE_SEED + 3000000L + (i - 1) * 10000L,
    cache_file = cache_file,
    n_cores    = N_CORES,
    cell_label = cell_label,
    verbose    = VERBOSE
  )

  panel_raw_list[[i]] <- list(raw = raw_results, row = row)
}

# Aggregate panel-level metrics per cell
panel_summaries <- list()

for (i in seq_along(panel_raw_list)) {
  item <- panel_raw_list[[i]]
  row  <- item$row

  panel_agg <- aggregate_panel_cell(item$raw, alpha = 0.05, fdr_method = "BH")

  # Attach DGP labels
  panel_agg$n_per_group     <- row$n_per_group
  panel_agg$effect_size     <- row$effect_size
  panel_agg$confounder_type <- row$confounder_type
  panel_agg$cell_id         <- row$cell_id

  panel_summaries[[i]] <- panel_agg
}

panel_summary_df <- bind_rows(panel_summaries)

# Save panel summary (merge with existing when running extra sizes)
panel_summary_file <- file.path(CACHE_DIR, paste0(BENCHMARK, "_panel_summary.rds"))
if (!is.null(EXTRA_SIZES)) {
  merge_with_existing_summary(panel_summary_df, panel_summary_file, verbose = VERBOSE)
} else {
  saveRDS(panel_summary_df, panel_summary_file)
}
if (VERBOSE) {
  message("\nPanel summary saved -> ", panel_summary_file)
  if ("method" %in% names(panel_summary_df) && nrow(panel_summary_df) > 0) {
    message("\n=== Panel Summary (selected) ===")
    print(as.data.frame(
      panel_summary_df %>%
        select(n_per_group, effect_size, confounder_type, method,
               n_reps, panel_fdr, panel_sensitivity, signal_mean_bias)
    ), right = FALSE)
  } else {
    message("WARNING: Panel aggregation returned no results")
  }
}

# Save per-analyte effect recovery data for standalone plotting
panel_effect_recovery <- list()
for (i in seq_along(panel_raw_list)) {
  item <- panel_raw_list[[i]]
  row  <- item$row
  raw  <- item$raw
  if (row$effect_size == 0) next
  sig_rows <- raw[!raw$error & raw$is_signal & raw$converged, , drop = FALSE]
  if (nrow(sig_rows) == 0) next
  sig_summary <- sig_rows %>%
    group_by(method, analyte) %>%
    summarise(
      mean_estimate  = mean(estimate, na.rm = TRUE),
      sd_estimate    = sd(estimate, na.rm = TRUE),
      true_effect    = first(true_effect),
      n_reps         = n(),
      .groups = "drop"
    ) %>%
    mutate(
      bias = mean_estimate - true_effect,
      n_per_group     = row$n_per_group,
      effect_size     = row$effect_size,
      confounder_type = row$confounder_type
    )
  panel_effect_recovery[[length(panel_effect_recovery) + 1]] <- sig_summary
}
if (length(panel_effect_recovery) > 0) {
  effect_recovery_df <- bind_rows(panel_effect_recovery)
  effect_recovery_file <- file.path(CACHE_DIR,
                                    paste0(BENCHMARK, "_panel_effect_recovery.rds"))
  if (!is.null(EXTRA_SIZES)) {
    merge_with_existing_summary(effect_recovery_df, effect_recovery_file, verbose = VERBOSE)
  } else {
    saveRDS(effect_recovery_df, effect_recovery_file)
  }
  if (VERBOSE) message("Panel effect recovery saved -> ", effect_recovery_file)
}



# Figures are generated by the standalone plot_B3.R script.
# Run: Rscript inst/simulations/plot_B3.R [--cache_dir PATH]


# ---- Formatted Summary Table ------------------------------------------------

summary_formatted <- format_summary_table(summary_annotated)

if (VERBOSE) {
  message("\n=== Benchmark 3: Summary (selected cells) ===\n")

  # Print key scenarios: confounding, no outliers, homoscedastic, symmetric
  key_summary <- summary_formatted %>%
    filter(outlier_rate == 0, !heteroscedastic, !nonlinear_cov,
           outlier_direction == "symmetric") %>%
    select(n_per_group, effect_size, confounder_type, method,
           n_reps, bias_fmt, coverage_fmt, rejection_fmt,
           starts_with("pass_"))

  print(as.data.frame(key_summary), right = FALSE)
}

if (VERBOSE) {
  message("\n=== Benchmark 3 (v3) complete ===")
  if (exists("fig_dir")) message("Figures -> ", fig_dir)
  message("Single-analyte summary -> ", summary_file)
  if (exists("panel_summary_file")) message("Panel summary -> ", panel_summary_file)
}

}  # end if (!is.null(PHASE_B)) / else
