#!/usr/bin/env Rscript
# benchmark_censoring.R — Benchmark 1: Censoring-Aware vs. Naive Methods
#
# ADEMP Structure:
#   Aims:    Quantify bias, coverage, and power differences between
#            censoring-aware and naive modeling across censoring rates and
#            sample sizes.
#   DGP:     simulate_immunoassay() with clean assay process (no plate/lot/
#            replicate effects). Single analyte or 10-analyte panel,
#            cross-sectional.
#   Estimand: Marginal group effect beta_group on log-concentration.
#   Methods: (1) Oracle, (2) Tobit, (3) AFT, (4) Gamma+LOD/2,
#            (5) Gaussian+LOD/2, (6) Naive lm()+LOD/2.
#   Performance: bias, coverage, RMSE, power, type I error, convergence,
#                computation time.
#
# Usage:
#   Rscript benchmark_censoring.R [--n_reps N] [--n_cores N] [--cache_dir DIR]
#   Rscript benchmark_censoring.R --phase_b 1 [--n_reps N] [--n_cores N]
#
# Results are cached as .rds files; re-running tops up to the target
# replication count without re-running existing replications.
#
# Design notes (v3, backbone harmonization + panel extension):
#   - v3 grid: n_subjects={20,40,100,200}, effect_size={0,0.5,1.0} (60 cells).
#     (v2: n_subjects={20,50,100,200}, effect_size={0,0.2,0.5,1.0}, 80 cells.)
#   - Multi-analyte panel extension: n_analytes=10, n_signal=3 with per-analyte
#     fitting and BH FDR correction. Separate panel sub-grid at backbone
#     n={40,100,200} x effect={0,0.5,1.0} x censoring={0,0.30,0.70} (27 cells).
#   - Paired seed design (from v2): seeds are keyed by (n_subjects, effect_size)
#     only, NOT by censoring_target. This ensures identical latent data across
#     censoring levels, isolating the censoring effect and making Oracle
#     results constant within each design group.
#   - Cache versioning via CACHE_VERSION to invalidate stale caches after
#     seed scheme changes.
#   - Sub-analyses: auto-family selection, noise sensitivity (residual_sd=1.5).
#   - Phase B2: informative censoring sub-grid (--phase_b). Treatment group
#     gets lower LOD (lod_ratio < 1), creating differential censoring.
#     54-cell grid: censoring={0.10,0.30,0.50} x n={40,100,200} x
#     effect={0,0.5,1.0} x lod_ratio={0.5,0.75}. Cell IDs offset by 20000.
#   - Figures saved in both PNG and PDF.
#   - Post-hoc Oracle sanity diagnostics verify paired design correctness.
#
# Original notes:
#   - Signal sparsity factor (all vs. 5/20 differ) from the plan is dropped
#     because fit_one() fits each analyte independently, so multi-analyte
#     sparsity does not affect single-analyte estimation.
#   - AFT uses survreg(dist="lognormal") on raw-scale data via fit_one(),
#     correctly modeling log(Y_raw) ~ Normal(Xb, sigma^2).
#   - Gamma+LOD/2 uses raw-scale data with Gamma(link="log") and manual
#     LOD/2 substitution. The log-link coefficient estimates the same
#     estimand (log-scale group effect) as the other methods.
#   - Gaussian+LOD/2 is identical to Naive LOD/2 in bias/RMSE; retained to
#     demonstrate the z-interval vs. t-interval coverage difference.

# ---- Setup ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# Load immunoPlex
if (requireNamespace("devtools", quietly = TRUE) &&
    file.exists(file.path(getwd(), "DESCRIPTION"))) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(immunoPlex)
}

# Source benchmark helpers
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

N_REPS     <- as.integer(parse_cli_arg("--n_reps",    "200"))
N_CORES    <- as.integer(parse_cli_arg("--n_cores",   "1"))
CACHE_DIR  <- parse_cli_arg("--cache_dir", file.path(script_dir, "cache"))
EXTRA_N    <- parse_cli_arg("--extra_n", NULL)
PHASE_B    <- parse_cli_arg("--phase_b", NULL)
PHASE_C    <- parse_cli_arg("--phase_c", NULL)
BASE_SEED  <- 20241001L
BENCHMARK  <- "censoring"
VERBOSE    <- TRUE
# Bump CACHE_VERSION when the DGP or seed scheme changes to invalidate stale
# cache files. v2: paired seed design. v3: backbone harmonization + panel.
CACHE_VERSION <- 3L

# Extra-n supplement mode
EXTRA_SIZES  <- if (!is.null(EXTRA_N)) as.integer(strsplit(EXTRA_N, ",")[[1]]) else NULL
EXTRA_OFFSET <- 10000L

# Phase B2: informative censoring sub-grid offset
PHASE_B_OFFSET <- 20000L
# Phase C: noise sub-analysis fresh-namespace offset
PHASE_C_OFFSET <- 30000L

if (!dir.exists(CACHE_DIR)) dir.create(CACHE_DIR, recursive = TRUE)


# ---- DGP Grid ---------------------------------------------------------------

dgp_grid <- make_dgp_grid(
  censoring_target = c(0, 0.10, 0.30, 0.50, 0.70),
  n_subjects       = if (!is.null(EXTRA_SIZES)) EXTRA_SIZES else c(20L, 40L, 100L, 200L),
  effect_size      = c(0, 0.5, 1.0)
)
if (!is.null(EXTRA_SIZES)) dgp_grid$cell_id <- dgp_grid$cell_id + EXTRA_OFFSET

# Paired seed design: cells sharing (n_subjects, effect_size) get the same
# seed range so that only the LOD differs — the latent biology is identical.
# This eliminates confounding between RNG draws and censoring level, ensuring
# the Oracle produces identical results across censoring rates (as it should).
dgp_grid <- dgp_grid %>%
  group_by(n_subjects, effect_size) %>%
  mutate(design_id = cur_group_id()) %>%
  ungroup()

if (VERBOSE) {
  message("=== Benchmark 1: Censoring-Aware vs. Naive Methods ===")
  message("DGP grid: ", nrow(dgp_grid), " cells x ", N_REPS, " reps = ",
          nrow(dgp_grid) * N_REPS, " total runs")
  message("Design groups (shared seeds): ",
          max(dgp_grid$design_id), " unique (n, effect) combos")
  message("Cache dir: ", CACHE_DIR)
}


# ---- Method: Extract coefficient from survreg (Tobit/AFT) -------------------

.extract_survreg <- function(fit_obj, method_name) {
  fail_row <- data.frame(
    method = method_name, estimate = NA_real_, se = NA_real_,
    ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
    converged = FALSE, time_s = NA_real_, stringsAsFactors = FALSE
  )
  if (!fit_obj$converged || is.null(fit_obj$model)) return(fail_row)

  coef_tab <- tryCatch(summary(fit_obj$model)$table, error = function(e) NULL)
  if (is.null(coef_tab)) return(fail_row)

  group_rows <- grep("^group", rownames(coef_tab))
  if (length(group_rows) == 0) return(fail_row)

  est  <- coef_tab[group_rows[1], "Value"]
  se   <- coef_tab[group_rows[1], "Std. Error"]
  zval <- est / se
  pval <- 2 * stats::pnorm(-abs(zval))

  data.frame(
    method = method_name, estimate = est, se = se,
    ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se,
    p_value = pval, converged = TRUE, time_s = NA_real_,
    stringsAsFactors = FALSE
  )
}


# ---- Method: Extract coefficient from glmmTMB (Gamma/Gaussian) --------------

.extract_glmmtmb <- function(fit_obj, method_name) {
  fail_row <- data.frame(
    method = method_name, estimate = NA_real_, se = NA_real_,
    ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
    converged = FALSE, time_s = NA_real_, stringsAsFactors = FALSE
  )
  if (!fit_obj$converged || is.null(fit_obj$model)) return(fail_row)

  coef_tab <- tryCatch(
    summary(fit_obj$model)$coefficients$cond,
    error = function(e) NULL
  )
  if (is.null(coef_tab)) return(fail_row)

  group_rows <- grep("^group", rownames(coef_tab))
  if (length(group_rows) == 0) return(fail_row)

  est  <- coef_tab[group_rows[1], "Estimate"]
  se   <- coef_tab[group_rows[1], "Std. Error"]
  pval <- coef_tab[group_rows[1], "Pr(>|z|)"]

  data.frame(
    method = method_name, estimate = est, se = se,
    ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se,
    p_value = pval, converged = TRUE, time_s = NA_real_,
    stringsAsFactors = FALSE
  )
}


# ---- Method: Extract coefficient from lm ------------------------------------

.extract_lm <- function(fit, method_name) {
  fail_row <- data.frame(
    method = method_name, estimate = NA_real_, se = NA_real_,
    ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
    converged = FALSE, time_s = NA_real_, stringsAsFactors = FALSE
  )

  coef_tab <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  if (is.null(coef_tab)) return(fail_row)

  group_rows <- grep("^group", rownames(coef_tab))
  if (length(group_rows) == 0) return(fail_row)

  est  <- coef_tab[group_rows[1], "Estimate"]
  se   <- coef_tab[group_rows[1], "Std. Error"]
  pval <- coef_tab[group_rows[1], "Pr(>|t|)"]
  ci   <- tryCatch(
    stats::confint(fit)[grep("^group", rownames(coef_tab))[1], ],
    error = function(e) c(est - 1.96 * se, est + 1.96 * se)
  )

  data.frame(
    method = method_name, estimate = est, se = se,
    ci_lo = ci[1], ci_hi = ci[2],
    p_value = pval, converged = TRUE, time_s = NA_real_,
    stringsAsFactors = FALSE
  )
}


# ---- Single-Replication Function --------------------------------------------

#' Run one replication for a single DGP cell
#'
#' Generates data, applies all 6 methods, returns one row per method.
#' When n_analytes > 1, loops over analytes and returns one row per
#' method x analyte (e.g. 60 rows for 6 methods x 10 analytes).
#'
#' @param seed          RNG seed.
#' @param n_subjects    Sample size.
#' @param effect_size   True group effect (log-scale); 0 = null.
#' @param censoring_target Target censoring fraction (0–0.70).
#' @param residual_sd   Residual SD on log scale (default 1.0).
#' @param n_analytes    Number of analytes (1 = single-analyte, >1 = panel).
#' @param n_signal      Number of signal analytes when n_analytes > 1.
#' @return data.frame. Single-analyte: columns method, estimate, se, ci_lo,
#'   ci_hi, p_value, converged, time_s, true_effect, realized_censoring.
#'   Multi-analyte: additionally analyte, is_signal, with true_effect
#'   set per analyte (effect_size for signal, 0 for null).
run_one_rep <- function(seed, n_subjects, effect_size, censoring_target,
                        residual_sd = 1.0, n_analytes = 1L, n_signal = 3L) {

  # --- Multi-analyte panel path (n_analytes > 1) ---
  if (n_analytes > 1L) {
    sim <- simulate_immunoassay(
      n_subjects      = n_subjects,
      n_timepoints    = 1L,
      n_analytes      = n_analytes,
      design          = "cross_sectional",
      group_levels    = c("control", "treatment"),
      group_effects   = effect_size,
      signal_analytes = if (effect_size != 0) n_signal else 0L,
      effect_direction = "up",
      re_intercept_sd = 0,
      residual_sd     = residual_sd,
      # For zero-censoring panel cells, use tiny lod_quantile to avoid
      # needing per-analyte lod_values (analyte names not known a priori).
      lod_quantile    = if (censoring_target > 0) censoring_target else 0.001,
      seed            = seed
    )

    d_all  <- as.data.frame(sim$data[!is.na(sim$data$value), ])
    tr_all <- as.data.frame(sim$truth[!is.na(sim$data$value), ])

    analyte_names_vec <- unique(d_all$cytokine)
    sig_lookup <- tapply(tr_all$signal_analyte, d_all$cytokine,
                          function(x) as.logical(x[1]))

    all_results <- list()

    for (aname in analyte_names_vec) {
      aidx <- d_all$cytokine == aname
      d  <- d_all[aidx, , drop = FALSE]
      tr <- tr_all[aidx, , drop = FALSE]
      is_sig   <- sig_lookup[[aname]]
      true_eff <- if (is_sig) effect_size else 0
      rcens    <- mean(d$cens_lod, na.rm = TRUE)

      res <- list()

      # Oracle
      oracle_fit <- tryCatch(
        stats::lm(value ~ group,
                  data = data.frame(value = tr$true_latent, group = d$group)),
        error = function(e) NULL)
      res[[1]] <- if (!is.null(oracle_fit)) {
        .extract_lm(oracle_fit, "Oracle")
      } else {
        data.frame(method = "Oracle", estimate = NA_real_, se = NA_real_,
                   ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
                   converged = FALSE, time_s = NA_real_,
                   stringsAsFactors = FALSE)
      }

      # Tobit
      tobit_fit <- tryCatch(
        fit_one(d, family = "tobit", fixed = "group", random = NULL),
        error = function(e) list(converged = FALSE, model = NULL))
      res[[2]] <- .extract_survreg(tobit_fit, "Tobit")

      # AFT
      aft_fit <- tryCatch(
        fit_one(d, family = "aft", fixed = "group", random = NULL),
        error = function(e) list(converged = FALSE, model = NULL))
      res[[3]] <- .extract_survreg(aft_fit, "AFT")

      # Gamma + LOD/2
      d_g <- d
      d_g$value <- d$value_raw
      blod_g <- d_g$cens_lod & !is.na(d_g$cens_lod)
      if (any(blod_g)) d_g$value[blod_g] <- d_g$lod[blod_g] / 2
      d_g$cens_lod  <- FALSE
      d_g$cens_ulod <- FALSE
      gamma_fit <- tryCatch(
        fit_one(d_g, family = "gamma", fixed = "group", random = NULL,
                impute_lod = FALSE),
        error = function(e) list(converged = FALSE, model = NULL))
      res[[4]] <- .extract_glmmtmb(gamma_fit, "Gamma+LOD/2")

      # Gaussian + LOD/2
      d_gs <- d
      blod <- d_gs$cens_lod & !is.na(d_gs$cens_lod)
      if (any(blod)) d_gs$value[blod] <- log(d_gs$lod[blod] / 2)
      d_gs$cens_lod  <- FALSE
      d_gs$cens_ulod <- FALSE
      gauss_fit <- tryCatch(
        fit_one(d_gs, family = "gaussian", fixed = "group", random = NULL),
        error = function(e) list(converged = FALSE, model = NULL))
      res[[5]] <- .extract_glmmtmb(gauss_fit, "Gaussian+LOD/2")

      # Naive lm() + LOD/2
      d_n <- d
      blod_n <- d_n$cens_lod & !is.na(d_n$cens_lod)
      if (any(blod_n)) d_n$value[blod_n] <- log(d_n$lod[blod_n] / 2)
      naive_fit <- tryCatch(
        stats::lm(value ~ group, data = d_n),
        error = function(e) NULL)
      res[[6]] <- if (!is.null(naive_fit)) {
        .extract_lm(naive_fit, "Naive LOD/2")
      } else {
        data.frame(method = "Naive LOD/2", estimate = NA_real_,
                   se = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                   p_value = NA_real_, converged = FALSE, time_s = NA_real_,
                   stringsAsFactors = FALSE)
      }

      a_res <- do.call(rbind, res)
      a_res$analyte            <- aname
      a_res$is_signal          <- is_sig
      a_res$true_effect        <- true_eff
      a_res$realized_censoring <- rcens
      a_res$converged[is.na(a_res$converged)] <- FALSE
      all_results[[length(all_results) + 1]] <- a_res
    }

    return(do.call(rbind, all_results))
  }

  # --- Simulate data (clean assay: no plate/lot/replicate effects) ---
  sim <- simulate_immunoassay(
    n_subjects      = n_subjects,
    n_timepoints    = 1L,
    n_analytes      = 1L,
    design          = "cross_sectional",
    group_levels    = c("control", "treatment"),
    group_effects   = effect_size,
    signal_analytes = if (effect_size != 0) 1L else 0L,
    effect_direction = "up",
    re_intercept_sd = 0,
    residual_sd     = residual_sd,
    # Censoring control
    lod_quantile    = if (censoring_target > 0) censoring_target else NULL,
    lod_values      = if (censoring_target == 0) c("IL-1b" = 1e-10) else NULL,
    seed            = seed
  )

  d <- as.data.frame(sim$data[!is.na(sim$data$value), ])
  tr <- as.data.frame(sim$truth[!is.na(sim$data$value), ])
  realized_cens <- mean(d$cens_lod, na.rm = TRUE)

  results_list <- list()

  # --- Method 1: Oracle (analyze latent uncensored truth) ---
  t0 <- proc.time()["elapsed"]
  d_oracle <- data.frame(
    value = tr$true_latent,
    group = d$group
  )
  oracle_fit <- tryCatch(
    stats::lm(value ~ group, data = d_oracle),
    error = function(e) NULL
  )
  time_oracle <- proc.time()["elapsed"] - t0
  if (!is.null(oracle_fit)) {
    r <- .extract_lm(oracle_fit, "Oracle")
    r$time_s <- time_oracle
    results_list[[1]] <- r
  } else {
    results_list[[1]] <- data.frame(
      method = "Oracle", estimate = NA_real_, se = NA_real_,
      ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
      converged = FALSE, time_s = time_oracle, stringsAsFactors = FALSE
    )
  }

  # --- Method 2: Tobit (censoring-aware, via fit_one) ---
  t0 <- proc.time()["elapsed"]
  tobit_fit <- tryCatch(
    fit_one(d, family = "tobit", fixed = "group", random = NULL),
    error = function(e) list(converged = FALSE, model = NULL)
  )
  time_tobit <- proc.time()["elapsed"] - t0
  r <- .extract_survreg(tobit_fit, "Tobit")
  r$time_s <- time_tobit
  results_list[[2]] <- r

  # --- Method 3: AFT (censoring-aware, lognormal via fit_one) ---
  t0 <- proc.time()["elapsed"]
  aft_fit <- tryCatch(
    fit_one(d, family = "aft", fixed = "group", random = NULL),
    error = function(e) list(converged = FALSE, model = NULL)
  )
  time_aft <- proc.time()["elapsed"] - t0
  r <- .extract_survreg(aft_fit, "AFT")
  r$time_s <- time_aft
  results_list[[3]] <- r

  # --- Method 4: Gamma + LOD/2 (raw-scale, Gamma GLMM, LOD/2 substitution) ---
  t0 <- proc.time()["elapsed"]
  d_gamma <- d
  d_gamma$value <- d$value_raw  # raw-scale response for Gamma
  below_lod_g <- d_gamma$cens_lod & !is.na(d_gamma$cens_lod)
  if (any(below_lod_g)) {
    d_gamma$value[below_lod_g] <- d_gamma$lod[below_lod_g] / 2
  }
  d_gamma$cens_lod  <- FALSE  # treated as observed after substitution
  d_gamma$cens_ulod <- FALSE
  gamma_fit <- tryCatch(
    fit_one(d_gamma, family = "gamma", fixed = "group", random = NULL,
            impute_lod = FALSE),
    error = function(e) list(converged = FALSE, model = NULL)
  )
  time_gamma <- proc.time()["elapsed"] - t0
  r <- .extract_glmmtmb(gamma_fit, "Gamma+LOD/2")
  r$time_s <- time_gamma
  results_list[[4]] <- r

  # --- Method 5: Gaussian + LOD/2 (log-scale, LOD/2 substitution) ---
  # NOTE: Gaussian+LOD/2 produces identical bias and RMSE to Naive LOD/2 (§5.1,
  # §6a). The only difference is CI construction: glmmTMB uses Wald z-intervals
  # while lm() uses t-intervals, causing minor coverage differences at small n.
  # Retained to demonstrate this inferential distinction.
  t0 <- proc.time()["elapsed"]
  d_gauss <- d
  below_lod <- d_gauss$cens_lod & !is.na(d_gauss$cens_lod)
  if (any(below_lod)) {
    d_gauss$value[below_lod] <- log(d_gauss$lod[below_lod] / 2)
  }
  d_gauss$cens_lod  <- FALSE
  d_gauss$cens_ulod <- FALSE
  gauss_fit <- tryCatch(
    fit_one(d_gauss, family = "gaussian", fixed = "group", random = NULL),
    error = function(e) list(converged = FALSE, model = NULL)
  )
  time_gauss <- proc.time()["elapsed"] - t0
  r <- .extract_glmmtmb(gauss_fit, "Gaussian+LOD/2")
  r$time_s <- time_gauss
  results_list[[5]] <- r

  # --- Method 6: Naive lm() + LOD/2 ---
  t0 <- proc.time()["elapsed"]
  d_naive <- d
  below_lod_n <- d_naive$cens_lod & !is.na(d_naive$cens_lod)
  if (any(below_lod_n)) {
    d_naive$value[below_lod_n] <- log(d_naive$lod[below_lod_n] / 2)
  }
  naive_fit <- tryCatch(
    stats::lm(value ~ group, data = d_naive),
    error = function(e) NULL
  )
  time_naive <- proc.time()["elapsed"] - t0
  if (!is.null(naive_fit)) {
    r <- .extract_lm(naive_fit, "Naive LOD/2")
    r$time_s <- time_naive
    results_list[[6]] <- r
  } else {
    results_list[[6]] <- data.frame(
      method = "Naive LOD/2", estimate = NA_real_, se = NA_real_,
      ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
      converged = FALSE, time_s = time_naive, stringsAsFactors = FALSE
    )
  }

  # --- Combine ---
  results <- do.call(rbind, results_list)
  results$true_effect <- effect_size
  results$realized_censoring <- realized_cens
  results
}


# ---- Single-Replication Function: Informative Censoring ---------------------

#' Run one replication with informative (differential) censoring
#'
#' Simulates data normally, then post-hoc lowers the LOD for treatment subjects
#' by a factor of \code{lod_ratio}, creating differential censoring between
#' groups. This mimics real-world batch effects where treatment group samples
#' are assayed with different LODs.
#'
#' @param seed            RNG seed.
#' @param n_subjects      Sample size.
#' @param effect_size     True group effect (log-scale); 0 = null.
#' @param censoring_target Target censoring fraction (0-0.50).
#' @param lod_ratio       Ratio of treatment LOD to control LOD (< 1 means
#'   treatment has a lower LOD, hence less censoring).
#' @param residual_sd     Residual SD on log scale (default 1.0).
#' @return data.frame with columns: method, estimate, se, ci_lo, ci_hi,
#'   p_value, converged, time_s, true_effect, realized_censoring,
#'   realized_censoring_control, realized_censoring_treatment.
run_one_rep_informative <- function(seed, n_subjects, effect_size,
                                    censoring_target, lod_ratio,
                                    residual_sd = 1.0) {

  # --- Simulate data with uniform LOD (set by censoring_target) ---
  sim <- simulate_immunoassay(
    n_subjects      = n_subjects,
    n_timepoints    = 1L,
    n_analytes      = 1L,
    design          = "cross_sectional",
    group_levels    = c("control", "treatment"),
    group_effects   = effect_size,
    signal_analytes = if (effect_size != 0) 1L else 0L,
    effect_direction = "up",
    re_intercept_sd = 0,
    residual_sd     = residual_sd,
    lod_quantile    = if (censoring_target > 0) censoring_target else NULL,
    lod_values      = if (censoring_target == 0) c("IL-1b" = 1e-10) else NULL,
    seed            = seed
  )

  d  <- as.data.frame(sim$data[!is.na(sim$data$value), ])
  tr <- as.data.frame(sim$truth[!is.na(sim$data$value), ])

  # --- Post-hoc: lower LOD for treatment group ---
  is_treatment <- d$group == "treatment"
  lod_control  <- d$lod[!is_treatment][1]  # uniform LOD from simulation
  lod_treat    <- lod_control * lod_ratio

  # Update LOD for treatment subjects
  d$lod[is_treatment] <- lod_treat

  # Re-apply censoring for treatment subjects based on their new (lower) LOD.
  # Values below the new treatment LOD are censored; values between
  # lod_treatment and lod_control that were previously censored are now observed.
  # Work on raw scale since lod is on raw scale.
  raw_values_treat <- tr$true_raw[is_treatment]
  new_cens_treat   <- raw_values_treat < lod_treat

  # Update censoring flag and observed values for treatment subjects

  d$cens_lod[is_treatment] <- new_cens_treat
  # For newly uncensored treatment values, restore the true observed value
  d$value_raw[is_treatment] <- ifelse(
    new_cens_treat,
    lod_treat,         # censored at the new treatment LOD
    raw_values_treat   # uncensored: use true raw value
  )
  d$value[is_treatment] <- log(d$value_raw[is_treatment])

  realized_cens         <- mean(d$cens_lod, na.rm = TRUE)
  realized_cens_control <- mean(d$cens_lod[!is_treatment], na.rm = TRUE)
  realized_cens_treat   <- mean(d$cens_lod[is_treatment], na.rm = TRUE)

  results_list <- list()

  # --- Method 1: Oracle (analyze latent uncensored truth) ---
  t0 <- proc.time()["elapsed"]
  d_oracle <- data.frame(
    value = tr$true_latent,
    group = d$group
  )
  oracle_fit <- tryCatch(
    stats::lm(value ~ group, data = d_oracle),
    error = function(e) NULL
  )
  time_oracle <- proc.time()["elapsed"] - t0
  if (!is.null(oracle_fit)) {
    r <- .extract_lm(oracle_fit, "Oracle")
    r$time_s <- time_oracle
    results_list[[1]] <- r
  } else {
    results_list[[1]] <- data.frame(
      method = "Oracle", estimate = NA_real_, se = NA_real_,
      ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
      converged = FALSE, time_s = time_oracle, stringsAsFactors = FALSE
    )
  }

  # --- Method 2: Tobit (censoring-aware, via fit_one) ---
  t0 <- proc.time()["elapsed"]
  tobit_fit <- tryCatch(
    fit_one(d, family = "tobit", fixed = "group", random = NULL),
    error = function(e) list(converged = FALSE, model = NULL)
  )
  time_tobit <- proc.time()["elapsed"] - t0
  r <- .extract_survreg(tobit_fit, "Tobit")
  r$time_s <- time_tobit
  results_list[[2]] <- r

  # --- Method 3: AFT (censoring-aware, lognormal via fit_one) ---
  t0 <- proc.time()["elapsed"]
  aft_fit <- tryCatch(
    fit_one(d, family = "aft", fixed = "group", random = NULL),
    error = function(e) list(converged = FALSE, model = NULL)
  )
  time_aft <- proc.time()["elapsed"] - t0
  r <- .extract_survreg(aft_fit, "AFT")
  r$time_s <- time_aft
  results_list[[3]] <- r

  # --- Method 4: Gamma + LOD/2 (raw-scale, Gamma GLMM, LOD/2 substitution) ---
  t0 <- proc.time()["elapsed"]
  d_gamma <- d
  d_gamma$value <- d$value_raw  # raw-scale response for Gamma
  below_lod_g <- d_gamma$cens_lod & !is.na(d_gamma$cens_lod)
  if (any(below_lod_g)) {
    d_gamma$value[below_lod_g] <- d_gamma$lod[below_lod_g] / 2
  }
  d_gamma$cens_lod  <- FALSE
  d_gamma$cens_ulod <- FALSE
  gamma_fit <- tryCatch(
    fit_one(d_gamma, family = "gamma", fixed = "group", random = NULL,
            impute_lod = FALSE),
    error = function(e) list(converged = FALSE, model = NULL)
  )
  time_gamma <- proc.time()["elapsed"] - t0
  r <- .extract_glmmtmb(gamma_fit, "Gamma+LOD/2")
  r$time_s <- time_gamma
  results_list[[4]] <- r

  # --- Method 5: Gaussian + LOD/2 (log-scale, LOD/2 substitution) ---
  t0 <- proc.time()["elapsed"]
  d_gauss <- d
  below_lod <- d_gauss$cens_lod & !is.na(d_gauss$cens_lod)
  if (any(below_lod)) {
    d_gauss$value[below_lod] <- log(d_gauss$lod[below_lod] / 2)
  }
  d_gauss$cens_lod  <- FALSE
  d_gauss$cens_ulod <- FALSE
  gauss_fit <- tryCatch(
    fit_one(d_gauss, family = "gaussian", fixed = "group", random = NULL),
    error = function(e) list(converged = FALSE, model = NULL)
  )
  time_gauss <- proc.time()["elapsed"] - t0
  r <- .extract_glmmtmb(gauss_fit, "Gaussian+LOD/2")
  r$time_s <- time_gauss
  results_list[[5]] <- r

  # --- Method 6: Naive lm() + LOD/2 ---
  t0 <- proc.time()["elapsed"]
  d_naive <- d
  below_lod_n <- d_naive$cens_lod & !is.na(d_naive$cens_lod)
  if (any(below_lod_n)) {
    d_naive$value[below_lod_n] <- log(d_naive$lod[below_lod_n] / 2)
  }
  naive_fit <- tryCatch(
    stats::lm(value ~ group, data = d_naive),
    error = function(e) NULL
  )
  time_naive <- proc.time()["elapsed"] - t0
  if (!is.null(naive_fit)) {
    r <- .extract_lm(naive_fit, "Naive LOD/2")
    r$time_s <- time_naive
    results_list[[6]] <- r
  } else {
    results_list[[6]] <- data.frame(
      method = "Naive LOD/2", estimate = NA_real_, se = NA_real_,
      ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
      converged = FALSE, time_s = time_naive, stringsAsFactors = FALSE
    )
  }

  # --- Combine ---
  results <- do.call(rbind, results_list)
  results$true_effect               <- effect_size
  results$realized_censoring        <- realized_cens
  results$realized_censoring_control   <- realized_cens_control
  results$realized_censoring_treatment <- realized_cens_treat
  results
}


# ---- Run Benchmark ----------------------------------------------------------

run_benchmark_censoring <- function(dgp_grid, n_reps, base_seed,
                                    cache_dir, n_cores, verbose) {
  all_summaries <- list()

  for (i in seq_len(nrow(dgp_grid))) {
    row <- dgp_grid[i, ]
    cell_label <- sprintf("cell_%d [cens=%.0f%%/n=%d/eff=%.1f]",
                          row$cell_id, row$censoring_target * 100,
                          row$n_subjects, row$effect_size)

    rep_fn <- function(seed) {
      run_one_rep(
        seed             = seed,
        n_subjects       = row$n_subjects,
        effect_size      = row$effect_size,
        censoring_target = row$censoring_target
      )
    }

    cache_file <- file.path(cache_dir,
                            paste0(BENCHMARK, "_v", CACHE_VERSION,
                                   "_cell", row$cell_id, ".rds"))

    # Paired design: seed depends on design_id (n_subjects x effect_size),
    # NOT on censoring_target. Same latent data across censoring levels.
    raw_results <- run_benchmark_cell(
      rep_fn     = rep_fn,
      n_reps     = n_reps,
      base_seed  = base_seed + (row$design_id - 1) * 10000L,
      cache_file = cache_file,
      n_cores    = n_cores,
      cell_label = cell_label,
      verbose    = verbose
    )

    # Aggregate per method
    methods <- unique(raw_results$method[!is.na(raw_results$method)])
    for (m in methods) {
      m_reps <- raw_results[!is.na(raw_results$method) &
                            raw_results$method == m, ]
      agg <- aggregate_cell(m_reps, truth = row$effect_size)

      # Add timing summary
      agg$mean_time_s <- mean(m_reps$time_s, na.rm = TRUE)

      # Cell identifiers
      agg$method            <- m
      agg$censoring_target  <- row$censoring_target
      agg$n_subjects        <- row$n_subjects
      agg$effect_size       <- row$effect_size
      agg$cell_id           <- row$cell_id
      all_summaries[[length(all_summaries) + 1]] <- agg
    }
  }

  bind_rows(all_summaries)
}


# ---- Phase B2: Informative Censoring Sub-Grid --------------------------------
#
# When --phase_b is passed, run ONLY this sub-grid and exit.
# Investigates differential censoring where treatment group has a lower LOD
# (e.g., due to batch effects), creating informative censoring that may bias
# group comparisons. Key question: does Tobit handle informative censoring
# better than LOD/2 substitution?

if (!is.null(PHASE_C) && !identical(PHASE_C, "noise")) {
  stop("Unknown PHASE_C mode: ", PHASE_C, " (expected 'noise')")
}

if (!is.null(PHASE_C) && identical(PHASE_C, "noise")) {

  # ---- Phase C: Noise Sub-Analysis (residual_sd = 1.5) ----------------------
  # Standalone re-run of the noise sub-analysis with cell_ids offset by
  # PHASE_C_OFFSET so caches are in a fresh namespace and do not collide
  # with the main-grid v3 noise caches.

  if (VERBOSE) message("\n=== Phase C: Noise Sub-Analysis (residual_sd = 1.5) ===")

  noise_grid <- make_dgp_grid(
    censoring_target = c(0, 0.30, 0.70),
    n_subjects       = c(50L, 200L),
    effect_size      = c(0, 0.5)
  )
  noise_grid <- noise_grid %>%
    group_by(n_subjects, effect_size) %>%
    mutate(design_id = cur_group_id()) %>%
    ungroup()
  noise_grid$cell_id <- noise_grid$cell_id + PHASE_C_OFFSET

  if (VERBOSE) {
    message("Noise grid: ", nrow(noise_grid), " cells x ", N_REPS, " reps")
  }

  noise_summaries_c <- list()
  for (i in seq_len(nrow(noise_grid))) {
    row <- noise_grid[i, ]
    cell_label <- sprintf("phaseC_noise_%d [cens=%.0f%%/n=%d/eff=%.1f/sd=1.5]",
                          row$cell_id, row$censoring_target * 100,
                          row$n_subjects, row$effect_size)

    rep_fn <- local({
      ns <- row$n_subjects; es <- row$effect_size; ct <- row$censoring_target
      function(seed) {
        run_one_rep(seed = seed, n_subjects = ns, effect_size = es,
                    censoring_target = ct, residual_sd = 1.5)
      }
    })

    cache_file <- file.path(CACHE_DIR,
                            paste0(BENCHMARK, "_v", CACHE_VERSION,
                                   "_phaseC_noise_cell", row$cell_id, ".rds"))

    raw_results <- run_benchmark_cell(
      rep_fn     = rep_fn,
      n_reps     = N_REPS,
      base_seed  = BASE_SEED + 800000L + (row$design_id - 1) * 10000L,
      cache_file = cache_file,
      n_cores    = N_CORES,
      cell_label = cell_label,
      verbose    = VERBOSE
    )

    methods <- unique(raw_results$method[!is.na(raw_results$method)])
    for (m in methods) {
      m_reps <- raw_results[!is.na(raw_results$method) &
                            raw_results$method == m, ]
      agg <- aggregate_cell(m_reps, truth = row$effect_size)
      agg$mean_time_s       <- mean(m_reps$time_s, na.rm = TRUE)
      agg$method            <- m
      agg$censoring_target  <- row$censoring_target
      agg$n_subjects        <- row$n_subjects
      agg$effect_size       <- row$effect_size
      agg$residual_sd       <- 1.5
      agg$cell_id           <- row$cell_id
      noise_summaries_c[[length(noise_summaries_c) + 1]] <- agg
    }
  }

  noise_summary_c <- bind_rows(noise_summaries_c)
  noise_file_c <- file.path(CACHE_DIR,
                            paste0(BENCHMARK, "_phaseC_noise_summary.rds"))
  saveRDS(noise_summary_c, noise_file_c)
  if (VERBOSE) message("Phase C noise summary saved -> ", noise_file_c)

  # NOTE: Phase C cell_ids are offset by PHASE_C_OFFSET (30000) and therefore
  # do NOT overlap with any prior v3 noise cell_ids. Merging into the shared
  # noise_summary file would *append* rather than *replace*, which was
  # misleading for downstream plot scripts that expect a single canonical
  # source. The standalone file above (censoring_phaseC_noise_summary.rds) is
  # the authoritative Phase C noise summary; plot scripts should read it
  # directly.

  if (VERBOSE) message("Phase C noise mode: exiting.")
  quit(save = "no", status = 0)
}

if (!is.null(PHASE_B)) {

  phaseB_grid <- make_dgp_grid(
    censoring_target = c(0.10, 0.30, 0.50),
    n_subjects       = c(40L, 100L, 200L),
    effect_size      = c(0, 0.5, 1.0),
    lod_ratio        = c(0.5, 0.75)
  )
  phaseB_grid$cell_id <- phaseB_grid$cell_id + PHASE_B_OFFSET

  # Paired seed design: cells sharing (n_subjects, effect_size) get the same
  # seed range so that only LOD/censoring parameters differ.
  phaseB_grid <- phaseB_grid %>%
    group_by(n_subjects, effect_size) %>%
    mutate(design_id = cur_group_id()) %>%
    ungroup()

  if (VERBOSE) {
    message("=== Phase B2: Informative Censoring Sub-Grid ===")
    message("DGP grid: ", nrow(phaseB_grid), " cells x ", N_REPS, " reps = ",
            nrow(phaseB_grid) * N_REPS, " total runs")
    message("Design groups (shared seeds): ",
            max(phaseB_grid$design_id), " unique (n, effect) combos")
    message("Cache dir: ", CACHE_DIR)
  }

  phaseB_summaries <- list()

  for (i in seq_len(nrow(phaseB_grid))) {
    row <- phaseB_grid[i, ]
    cell_label <- sprintf(
      "phaseB_%d [cens=%.0f%%/n=%d/eff=%.1f/lod_ratio=%.2f]",
      row$cell_id, row$censoring_target * 100,
      row$n_subjects, row$effect_size, row$lod_ratio
    )

    # Use local() for safe closure binding in loop
    rep_fn <- local({
      ns <- row$n_subjects
      es <- row$effect_size
      ct <- row$censoring_target
      lr <- row$lod_ratio
      function(seed) {
        run_one_rep_informative(
          seed             = seed,
          n_subjects       = ns,
          effect_size      = es,
          censoring_target = ct,
          lod_ratio        = lr
        )
      }
    })

    cache_file <- file.path(CACHE_DIR,
                            paste0(BENCHMARK, "_v", CACHE_VERSION,
                                   "_phaseB_cell", row$cell_id, ".rds"))

    # Paired design: seed depends on design_id (n_subjects x effect_size),
    # NOT on censoring_target or lod_ratio. Offset by 600000 to avoid
    # collision with main grid (0-based), panel (500000), noise (800000),
    # and auto (900000) sub-grids.
    raw_results <- run_benchmark_cell(
      rep_fn     = rep_fn,
      n_reps     = N_REPS,
      base_seed  = BASE_SEED + 600000L + (row$design_id - 1) * 10000L,
      cache_file = cache_file,
      n_cores    = N_CORES,
      cell_label = cell_label,
      verbose    = VERBOSE
    )

    # Aggregate per method
    methods <- unique(raw_results$method[!is.na(raw_results$method)])
    for (m in methods) {
      m_reps <- raw_results[!is.na(raw_results$method) &
                            raw_results$method == m, ]
      agg <- aggregate_cell(m_reps, truth = row$effect_size)

      # Add timing summary
      agg$mean_time_s <- mean(m_reps$time_s, na.rm = TRUE)

      # Cell identifiers
      agg$method                <- m
      agg$censoring_target      <- row$censoring_target
      agg$n_subjects            <- row$n_subjects
      agg$effect_size           <- row$effect_size
      agg$lod_ratio             <- row$lod_ratio
      agg$cell_id               <- row$cell_id
      agg$informative_censoring <- TRUE

      # Per-group realized censoring (average across reps)
      if ("realized_censoring_control" %in% names(m_reps)) {
        agg$mean_realized_censoring_control   <-
          mean(m_reps$realized_censoring_control, na.rm = TRUE)
        agg$mean_realized_censoring_treatment <-
          mean(m_reps$realized_censoring_treatment, na.rm = TRUE)
      }

      phaseB_summaries[[length(phaseB_summaries) + 1]] <- agg
    }
  }

  phaseB_summary <- bind_rows(phaseB_summaries)

  # Add acceptance annotations
  phaseB_null <- phaseB_summary %>%
    filter(effect_size == 0) %>%
    mutate(
      pass_bias     = NA,
      pass_coverage = NA,
      pass_type1    = passes_threshold(rejection_rate, 0.025, 0.075)
    )

  phaseB_alt <- phaseB_summary %>%
    filter(effect_size != 0) %>%
    mutate(
      pass_bias     = passes_threshold(abs(bias_ratio), upper = 0.1),
      pass_coverage = passes_threshold(coverage, 0.93, 0.97),
      pass_type1    = NA
    )

  phaseB_annotated <- bind_rows(phaseB_null, phaseB_alt)

  # Save into the main summary file, tagged with informative_censoring = TRUE.
  # Existing main-grid rows get informative_censoring = FALSE if not already set.
  summary_file <- file.path(CACHE_DIR, paste0(BENCHMARK, "_v", CACHE_VERSION,
                                               "_summary.rds"))
  if (file.exists(summary_file)) {
    existing <- readRDS(summary_file)
    if (!"informative_censoring" %in% names(existing)) {
      existing$informative_censoring <- FALSE
    }
  }
  merge_with_existing_summary(phaseB_annotated, summary_file, verbose = VERBOSE)

  if (VERBOSE) {
    message("\n=== Phase B2: Informative Censoring Summary ===\n")

    method_order <- c("Oracle", "Tobit", "AFT", "Gamma+LOD/2",
                      "Gaussian+LOD/2", "Naive LOD/2")
    phaseB_annotated$method <- factor(phaseB_annotated$method,
                                       levels = method_order)

    phaseB_formatted <- format_summary_table(phaseB_annotated)
    print(
      phaseB_formatted %>%
        filter(effect_size %in% c(0, 0.5),
               n_subjects %in% c(40, 200)) %>%
        select(censoring_target, n_subjects, effect_size, lod_ratio,
               method, bias_fmt, coverage_fmt, rejection_fmt,
               starts_with("pass_")) %>%
        as.data.frame(),
      right = FALSE
    )

    message("\n=== Phase B2 complete ===")
    message("Summary merged -> ", summary_file)
  }

  # Exit after Phase B — skip main grid and other sub-analyses
  if (VERBOSE) message("Phase B mode: skipping main grid execution.")
  quit(save = "no", status = 0)
}


# ---- Execute ----------------------------------------------------------------

if (VERBOSE) message("\nStarting benchmark execution...")

summary_df <- run_benchmark_censoring(
  dgp_grid  = dgp_grid,
  n_reps    = N_REPS,
  base_seed = BASE_SEED,
  cache_dir = CACHE_DIR,
  n_cores   = N_CORES,
  verbose   = VERBOSE
)

# Add acceptance annotations
summary_null <- summary_df %>%
  filter(effect_size == 0) %>%
  mutate(
    pass_bias     = NA,
    pass_coverage = NA,
    pass_type1    = passes_threshold(rejection_rate, 0.025, 0.075)
  )

summary_alt <- summary_df %>%
  filter(effect_size != 0) %>%
  mutate(
    pass_bias     = passes_threshold(abs(bias_ratio), upper = 0.1),
    pass_coverage = passes_threshold(coverage, 0.93, 0.97),
    pass_type1    = NA
  )

summary_annotated <- bind_rows(summary_null, summary_alt)

# Tag main-grid rows as non-informative censoring for consistent schema
summary_annotated$informative_censoring <- FALSE

# Save final summary (merge with existing when running extra sizes)
summary_file <- file.path(CACHE_DIR, paste0(BENCHMARK, "_v", CACHE_VERSION,
                                             "_summary.rds"))
if (!is.null(EXTRA_SIZES)) {
  merge_with_existing_summary(summary_annotated, summary_file, verbose = VERBOSE)
} else {
  save_main_summary_preserving_supplements(summary_annotated, summary_file,
                                           verbose = VERBOSE)
}



# Method ordering for Formatted Summary Table
method_order <- c("Oracle", "Tobit", "AFT", "Gamma+LOD/2",
                  "Gaussian+LOD/2", "Naive LOD/2")
summary_annotated$method <- factor(summary_annotated$method,
                                    levels = method_order)

# Figures are generated by the standalone plot_B1.R script.
# Run: Rscript inst/simulations/plot_B1.R [--cache_dir PATH]


# ---- Sub-analysis: Auto-Family Model Selection -----------------------------

if (VERBOSE) message("\n--- Sub-analysis: Auto-Family Selection ---")

auto_grid <- make_dgp_grid(
  scenario      = c("gaussian_no_cens", "lognormal_light", "lognormal_heavy"),
  n_subjects    = if (!is.null(EXTRA_SIZES)) EXTRA_SIZES else c(50L, 100L),
  effect_size   = c(0, 0.5)
)
if (!is.null(EXTRA_SIZES)) auto_grid$cell_id <- auto_grid$cell_id + EXTRA_OFFSET

run_one_rep_auto <- function(seed, scenario, n_subjects, effect_size) {
  # Map scenario to DGP parameters
  dgp <- switch(scenario,
    gaussian_no_cens = list(censoring = 0,    residual_sd = 1.0),
    lognormal_light  = list(censoring = 0.15, residual_sd = 1.0),
    lognormal_heavy  = list(censoring = 0.50, residual_sd = 1.5)
  )

  sim <- simulate_immunoassay(
    n_subjects      = n_subjects,
    n_timepoints    = 1L,
    n_analytes      = 1L,
    design          = "cross_sectional",
    group_levels    = c("control", "treatment"),
    group_effects   = effect_size,
    signal_analytes = if (effect_size != 0) 1L else 0L,
    effect_direction = "up",
    re_intercept_sd = 0,
    residual_sd     = dgp$residual_sd,
    lod_quantile    = if (dgp$censoring > 0) dgp$censoring else NULL,
    lod_values      = if (dgp$censoring == 0) c("IL-1b" = 1e-10) else NULL,
    seed            = seed
  )

  d <- as.data.frame(sim$data[!is.na(sim$data$value), ])

  # Auto-family via fit_one
  auto_fit <- tryCatch(
    fit_one(d, family = "auto", fixed = "group", random = NULL),
    error = function(e) list(converged = FALSE, model = NULL, family = NA)
  )

  selected_family <- auto_fit$family

  # Extract coefficient based on selected family
  if (!auto_fit$converged || is.null(auto_fit$model)) {
    return(data.frame(
      method = "Auto", selected_family = selected_family,
      estimate = NA_real_, se = NA_real_,
      ci_lo = NA_real_, ci_hi = NA_real_,
      p_value = NA_real_, converged = FALSE,
      true_effect = effect_size, stringsAsFactors = FALSE
    ))
  }

  r <- if (selected_family %in% c("tobit", "aft")) {
    .extract_survreg(auto_fit, "Auto")
  } else {
    .extract_glmmtmb(auto_fit, "Auto")
  }

  r$selected_family <- selected_family
  r$true_effect     <- effect_size
  r
}

auto_summaries <- list()
for (i in seq_len(nrow(auto_grid))) {
  row <- auto_grid[i, ]
  cell_label <- sprintf("auto_%d [%s/n=%d/eff=%.1f]",
                        row$cell_id, row$scenario,
                        row$n_subjects, row$effect_size)

  rep_fn <- function(seed) {
    run_one_rep_auto(seed, row$scenario, row$n_subjects, row$effect_size)
  }

  cache_file <- file.path(CACHE_DIR,
                          paste0(BENCHMARK, "_v", CACHE_VERSION,
                                 "_auto_cell", row$cell_id, ".rds"))

  raw <- run_benchmark_cell(
    rep_fn     = rep_fn,
    n_reps     = N_REPS,
    base_seed  = BASE_SEED + 900000L + (i - 1) * 10000L,
    cache_file = cache_file,
    n_cores    = N_CORES,
    cell_label = cell_label,
    verbose    = VERBOSE
  )

  # Aggregate
  ok <- raw[!raw$error, ]
  agg <- aggregate_cell(ok, truth = row$effect_size)
  agg$scenario   <- row$scenario
  agg$n_subjects <- row$n_subjects
  agg$effect_size <- row$effect_size
  agg$cell_id    <- row$cell_id

  # Family selection frequencies
  if ("selected_family" %in% names(ok)) {
    freq <- table(ok$selected_family) / nrow(ok)
    agg$family_freq <- list(as.list(freq))
  }

  auto_summaries[[i]] <- agg
}

auto_summary <- bind_rows(auto_summaries)
auto_file <- file.path(CACHE_DIR, paste0(BENCHMARK, "_auto_summary.rds"))
if (!is.null(EXTRA_SIZES)) {
  merge_with_existing_summary(auto_summary, auto_file, verbose = VERBOSE)
} else {
  saveRDS(auto_summary, auto_file)
}
if (VERBOSE) message("Auto-family summary saved -> ", auto_file)


# ---- Sub-analysis: Noise Sensitivity (residual_sd = 1.5) -------------------
#
# Report §6c identified "single residual SD" as a design gap. This sub-analysis
# tests residual_sd = 1.5 (worse signal-to-noise) at a subset of conditions to
# check whether methods' relative performance changes under higher variability.

if (VERBOSE) message("\n--- Sub-analysis: Noise Sensitivity (residual_sd = 1.5) ---")

noise_grid <- make_dgp_grid(
  censoring_target = c(0, 0.30, 0.70),
  n_subjects       = if (!is.null(EXTRA_SIZES)) EXTRA_SIZES else c(50L, 200L),
  effect_size      = c(0, 0.5)
)

noise_grid <- noise_grid %>%
  group_by(n_subjects, effect_size) %>%
  mutate(design_id = cur_group_id()) %>%
  ungroup()
if (!is.null(EXTRA_SIZES)) noise_grid$cell_id <- noise_grid$cell_id + EXTRA_OFFSET

noise_summaries <- list()
for (i in seq_len(nrow(noise_grid))) {
  row <- noise_grid[i, ]
  cell_label <- sprintf("noise_%d [cens=%.0f%%/n=%d/eff=%.1f/sd=1.5]",
                        row$cell_id, row$censoring_target * 100,
                        row$n_subjects, row$effect_size)

  rep_fn <- function(seed) {
    run_one_rep(
      seed             = seed,
      n_subjects       = row$n_subjects,
      effect_size      = row$effect_size,
      censoring_target = row$censoring_target,
      residual_sd      = 1.5
    )
  }

  cache_file <- file.path(CACHE_DIR,
                          paste0(BENCHMARK, "_v", CACHE_VERSION,
                                 "_noise_cell", row$cell_id, ".rds"))

  raw_results <- run_benchmark_cell(
    rep_fn     = rep_fn,
    n_reps     = N_REPS,
    base_seed  = BASE_SEED + 800000L + (row$design_id - 1) * 10000L,
    cache_file = cache_file,
    n_cores    = N_CORES,
    cell_label = cell_label,
    verbose    = VERBOSE
  )

  methods <- unique(raw_results$method[!is.na(raw_results$method)])
  for (m in methods) {
    m_reps <- raw_results[!is.na(raw_results$method) &
                          raw_results$method == m, ]
    agg <- aggregate_cell(m_reps, truth = row$effect_size)
    agg$mean_time_s       <- mean(m_reps$time_s, na.rm = TRUE)
    agg$method            <- m
    agg$censoring_target  <- row$censoring_target
    agg$n_subjects        <- row$n_subjects
    agg$effect_size       <- row$effect_size
    agg$residual_sd       <- 1.5
    agg$cell_id           <- row$cell_id
    noise_summaries[[length(noise_summaries) + 1]] <- agg
  }
}

noise_summary <- bind_rows(noise_summaries)
noise_file <- file.path(CACHE_DIR, paste0(BENCHMARK, "_noise_summary.rds"))
if (!is.null(EXTRA_SIZES)) {
  merge_with_existing_summary(noise_summary, noise_file, verbose = VERBOSE)
} else {
  saveRDS(noise_summary, noise_file)
}
if (VERBOSE) message("Noise sensitivity summary saved -> ", noise_file)


# ---- Oracle Sanity Diagnostics -----------------------------------------------
#
# Post-hoc checks on Oracle performance. With the paired seed design, the
# Oracle should produce *identical* results across censoring levels for each
# (n_subjects, effect_size) cell. Any deviation indicates a DGP or seed issue.

if (VERBOSE) message("\n--- Oracle Sanity Diagnostics ---")

oracle_diag <- summary_annotated %>%
  filter(method == "Oracle")

# Check 1: Oracle bias should be near zero under the null
oracle_null <- oracle_diag %>% filter(effect_size == 0)
if (nrow(oracle_null) > 0) {
  worst_null_bias <- oracle_null %>%
    mutate(z_bias = abs(bias) / bias_mcse) %>%
    slice_max(z_bias, n = 1, with_ties = FALSE)
  if (VERBOSE) {
    message("  Oracle null bias (worst cell): ",
            sprintf("%.4f (MCSE %.4f, z = %.1f) at cens=%.0f%%, n=%d",
                    worst_null_bias$bias, worst_null_bias$bias_mcse,
                    worst_null_bias$z_bias,
                    worst_null_bias$censoring_target * 100,
                    worst_null_bias$n_subjects))
    if (worst_null_bias$z_bias > 3)
      message("  ** WARNING: Oracle null bias exceeds 3 MCSE **")
  }
}

# Check 2: Oracle Type I error should be near 0.05 across censoring levels
oracle_t1 <- oracle_diag %>% filter(effect_size == 0)
if (nrow(oracle_t1) > 0) {
  worst_t1 <- oracle_t1 %>%
    mutate(z_t1 = abs(rejection_rate - 0.05) / rejection_mcse) %>%
    slice_max(z_t1, n = 1, with_ties = FALSE)
  if (VERBOSE) {
    message("  Oracle Type I error (worst cell): ",
            sprintf("%.3f (MCSE %.3f, z = %.1f) at cens=%.0f%%, n=%d",
                    worst_t1$rejection_rate, worst_t1$rejection_mcse,
                    worst_t1$z_t1,
                    worst_t1$censoring_target * 100,
                    worst_t1$n_subjects))
    if (worst_t1$z_t1 > 3)
      message("  ** WARNING: Oracle Type I error deviates > 3 MCSE from 0.05 **")
  }
}

# Check 3: With paired design, Oracle should be constant across censoring
# (compare within each n_subjects x effect_size group)
oracle_paired <- oracle_diag %>%
  group_by(n_subjects, effect_size) %>%
  summarise(
    bias_range = max(bias, na.rm = TRUE) - min(bias, na.rm = TRUE),
    rmse_range = max(rmse, na.rm = TRUE) - min(rmse, na.rm = TRUE),
    .groups = "drop"
  )
max_bias_range <- max(oracle_paired$bias_range, na.rm = TRUE)
max_rmse_range <- max(oracle_paired$rmse_range, na.rm = TRUE)
if (VERBOSE) {
  message("  Oracle bias range across censoring (max): ",
          sprintf("%.6f", max_bias_range),
          if (max_bias_range < 1e-6) " [PASS: paired design confirmed]"
          else " [NOTE: non-zero range suggests seed pairing issue]")
  message("  Oracle RMSE range across censoring (max): ",
          sprintf("%.6f", max_rmse_range))
}


# ---- Formatted Summary Table ------------------------------------------------

summary_formatted <- format_summary_table(summary_annotated)

if (VERBOSE) {
  message("\n=== Benchmark 1: Summary (effect_size = 0.5, key cells) ===\n")
  print(
    summary_formatted %>%
      filter(effect_size %in% c(0, 0.5),
             n_subjects %in% c(40, 200)) %>%
      select(censoring_target, n_subjects, effect_size,
             method, bias_fmt, coverage_fmt, rejection_fmt,
             starts_with("pass_")) %>%
      as.data.frame(),
    right = FALSE
  )
}

if (VERBOSE) {
  message("\n=== Benchmark 1: Single-analyte analysis complete ===")
  if (exists("fig_dir")) message("Figures -> ", fig_dir)
  message("Summary -> ", summary_file)
}


# ---- Multi-Analyte Panel Extension ------------------------------------------
#
# Separate sub-grid for panel-level FDR analysis.
# n_analytes = 10, n_signal = 3 (30% signal prevalence).
# Per-analyte models are fit independently; BH FDR is applied across analytes
# within each replication x method by aggregate_panel_cell().

N_PANEL_ANALYTES <- 10L
N_PANEL_SIGNAL   <- 3L

panel_grid <- make_dgp_grid(
  censoring_target = c(0, 0.30, 0.70),
  n_subjects       = if (!is.null(EXTRA_SIZES)) EXTRA_SIZES else c(40L, 100L, 200L),
  effect_size      = c(0, 0.5, 1.0)
)

# Paired seed design for panel cells
panel_grid <- panel_grid %>%
  group_by(n_subjects, effect_size) %>%
  mutate(design_id = cur_group_id()) %>%
  ungroup()
if (!is.null(EXTRA_SIZES)) panel_grid$cell_id <- panel_grid$cell_id + EXTRA_OFFSET

if (VERBOSE) {
  message("\n=== Panel Extension: Multi-Analyte FDR Analysis ===")
  message("Panel grid: ", nrow(panel_grid), " cells x ", N_REPS, " reps = ",
          nrow(panel_grid) * N_REPS, " total runs")
  message("Analytes per sim: ", N_PANEL_ANALYTES,
          " (", N_PANEL_SIGNAL, " signal)")
}

panel_raw_all <- list()

for (i in seq_len(nrow(panel_grid))) {
  row <- panel_grid[i, ]
  cell_label <- sprintf("panel_%d [cens=%.0f%%/n=%d/eff=%.1f]",
                         row$cell_id, row$censoring_target * 100,
                         row$n_subjects, row$effect_size)

  # Use local() for safe closure binding in loop
  rep_fn <- local({
    ns <- row$n_subjects
    es <- row$effect_size
    ct <- row$censoring_target
    na <- N_PANEL_ANALYTES
    nsig <- N_PANEL_SIGNAL
    function(seed) {
      run_one_rep(seed = seed, n_subjects = ns, effect_size = es,
                  censoring_target = ct, n_analytes = na,
                  n_signal = nsig)
    }
  })

  cache_file <- file.path(CACHE_DIR,
                           paste0(BENCHMARK, "_v", CACHE_VERSION,
                                  "_panel_cell", row$cell_id, ".rds"))

  # Paired design: seed depends on design_id (n_subjects x effect_size),
  # NOT on censoring_target. Offset by 500000 to avoid main-grid collision.
  raw_results <- run_benchmark_cell(
    rep_fn     = rep_fn,
    n_reps     = N_REPS,
    base_seed  = BASE_SEED + 500000L + (row$design_id - 1) * 10000L,
    cache_file = cache_file,
    n_cores    = N_CORES,
    cell_label = cell_label,
    verbose    = VERBOSE
  )

  # Tag with DGP parameters (use same names as panel_summary for consistency)
  raw_results$censoring_target  <- row$censoring_target
  raw_results$n_subjects_dgp    <- row$n_subjects
  raw_results$effect_size_dgp   <- row$effect_size
  raw_results$panel_cell_id     <- row$cell_id

  panel_raw_all[[i]] <- raw_results
}


# ---- Panel Aggregation ------------------------------------------------------

if (VERBOSE) message("\nAggregating panel-level metrics...")

panel_summaries <- list()
for (i in seq_len(nrow(panel_grid))) {
  row <- panel_grid[i, ]
  raw <- panel_raw_all[[i]]

  panel_agg <- aggregate_panel_cell(raw, alpha = 0.05, fdr_method = "BH")
  panel_agg$censoring_target <- row$censoring_target
  panel_agg$n_subjects       <- row$n_subjects
  panel_agg$effect_size      <- row$effect_size
  panel_agg$cell_id          <- row$cell_id
  panel_summaries[[i]] <- panel_agg
}

panel_summary <- bind_rows(panel_summaries)

# Save panel summary
panel_file <- file.path(CACHE_DIR,
                         paste0(BENCHMARK, "_v", CACHE_VERSION,
                                "_panel_summary.rds"))
if (!is.null(EXTRA_SIZES)) {
  merge_with_existing_summary(panel_summary, panel_file, verbose = VERBOSE)
} else {
  saveRDS(panel_summary, panel_file)
}
if (VERBOSE) message("Panel summary saved -> ", panel_file)

# Save per-analyte averages for standalone plotting
panel_raw_combined <- bind_rows(panel_raw_all)
panel_analyte_avg <- panel_raw_combined %>%
  filter(!error, !is.na(method), !is.na(analyte)) %>%
  group_by(censoring_target, n_subjects_dgp, effect_size_dgp,
           method, analyte, is_signal) %>%
  summarise(
    mean_bias          = mean(estimate - true_effect, na.rm = TRUE),
    mean_realized_cens = mean(realized_censoring, na.rm = TRUE),
    .groups = "drop"
  )

panel_analyte_file <- file.path(CACHE_DIR,
                                paste0(BENCHMARK, "_v", CACHE_VERSION,
                                       "_panel_analyte_summary.rds"))
if (!is.null(EXTRA_SIZES)) {
  merge_with_existing_summary(panel_analyte_avg, panel_analyte_file, verbose = VERBOSE)
} else {
  saveRDS(panel_analyte_avg, panel_analyte_file)
}
if (VERBOSE) message("Panel analyte summary saved -> ", panel_analyte_file)



# Panel figures are generated by the standalone plot_B1.R script.


if (VERBOSE) {
  message("\n=== Benchmark 1 complete (single-analyte + panel) ===")
  if (exists("fig_dir")) message("Figures -> ", fig_dir)
  message("Summary -> ", summary_file)
  message("Panel summary -> ", panel_file)
}
