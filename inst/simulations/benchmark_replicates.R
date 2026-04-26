#!/usr/bin/env Rscript
# benchmark_replicates.R — Benchmark 6: Replicate Handling Strategy (v4)
#
# ADEMP Structure:
#   Aims:    Compare replicate handling strategies — naive averaging, QC-filtered
#            averaging, random-effect modeling, and ignoring replicate structure —
#            to quantify the impact on bias, coverage, and type I error.
#   DGP:     simulate_immunoassay() with n_replicates = 2 and varying replicate_sd.
#            Post-hoc outlier injection with varying magnitude (2 or 4 SD).
#            Single analyte, cross-sectional, 2 groups.  Sample size
#            varied (40/100/200).
#            Censoring varied (0% or 15%).  RE intercept SD varied (0 or 0.5).
#   Estimand: Marginal group effect beta_group on log-concentration.
#   Methods: (1) Naive average, (2) QC drop_both + average,
#            (3) QC drop_farther + average, (4) QC winsorize + average,
#            (5) Random effects (1|subject_id), (6) Ignore structure,
#            (7) Oracle average, (8) QC + Tobit, (9) Avg + Tobit.
#   Performance: bias, coverage, RMSE, power, type I error, relative efficiency,
#                QC detection rate.
#
# v2 changes (from B6 report recommendations):
#   - Added outlier_magnitude {2, 4} to test QC sensitivity near threshold
#   - Added lod_quantile {0, 0.15} to isolate censoring from replicate effects
#   - Added re_intercept_sd {0, 0.5} to test RE efficiency with real clustering
#   - Added --clear_cache flag to purge stale results
#   - Figures generated per (lod_quantile, re_intercept_sd) slice
#   - New cross-slice comparison figure (Fig 6)
#
# v3 changes (backbone harmonization):
#   - Added n_subjects = 200 (backbone alignment: 40/100/200)
#   - Grid: 360 -> 540 cells
#   - Cache prefix changed to "replicates_v3"
#
# v4 changes (benchmark extensions Phase A):
#   - Added Method 8: QC + Tobit (QC drop_both, average, then fit_one tobit)
#   - Added Method 9: Avg + Tobit (naive average, then fit_one tobit)
#   - Cache prefix changed to "replicates_v4"
#
# Phase B8 (weaker outlier magnitudes):
#   - When --phase_b is passed, runs ONLY a sub-grid with
#     outlier_magnitude = c(0.5, 1.0, 1.5) instead of c(2, 4).
#   - Tests whether QC methods can detect subtle outliers near the noise floor.
#   - Cell IDs offset by PHASE_B_OFFSET (20000).
#   - Cache prefix: "replicates_v4_phaseB"
#   - 648 cells (3 rep_sd × 2 out_rate × 3 n × 3 eff × 3 mag × 2 lod × 2 re)
#
# Usage:
#   Rscript benchmark_replicates.R [--n_reps N] [--n_cores N] [--cache_dir DIR]
#                                  [--clear_cache] [--phase_b]

# ---- Setup ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(glmmTMB)
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
has_cli_flag <- function(flag) flag %in% cli_args

N_REPS      <- as.integer(parse_cli_arg("--n_reps",    "200"))
N_CORES     <- as.integer(parse_cli_arg("--n_cores",   "1"))
CACHE_DIR   <- parse_cli_arg("--cache_dir", file.path(script_dir, "cache"))
EXTRA_N     <- parse_cli_arg("--extra_n", NULL)
PHASE_B     <- has_cli_flag("--phase_b")
CLEAR_CACHE <- has_cli_flag("--clear_cache")
BASE_SEED   <- 20250306L
BENCHMARK   <- if (PHASE_B) "replicates_v4_phaseB" else "replicates_v4"
VERBOSE     <- TRUE

# ---- Replicate-aggregation modes (scaffold) --------------------------------
#
# `--mode legacy` (default) runs B6 v4 exactly as before. `--mode mcnemar`
# and `--mode regression` swap in a different DGP that emits data with an
# explicit replicate axis -- the bivariate-normal latent from
# benchmark_mcnemar.R's run_one_rep(), then
#   value_ijkt_r = Z_{i,k,t}
#                + beta_t * signal_mask_k         (added on T2, signal analytes)
#                + epsilon_r ~ N(0, sigma_rep^2)  (rep-level noise)
#                + residual_{i,k,t,r} ~ N(0, sigma_res^2)
#
# mcnemar mode thresholds each rep value -> booleans, then compares
# replicate_agg in {majority_vote, any_detected, all_detected} via
# mcnemar_detection() (which aggregates rows within subject x cytokine x
# timepoint internally).
#
# regression mode feeds the continuous values to fit_models() both with
# rep_col = "rep_id" (RE append, glmmTMB adds (1|subject_id:rep_id)) and via
# aggregate_replicates(rule = "mean") + fit_models() pre-agg. Gaussian only
# in the scaffold; cens_lod / cens_ulod all FALSE, sigma_plate = 0 (plate
# deferred to a follow-up chat).
MODE        <- parse_cli_arg("--mode",       "legacy")
N_ANALYTES  <- as.integer(parse_cli_arg("--n_analytes", "10"))
N_SIGNAL    <- as.integer(parse_cli_arg("--n_signal",   "3"))

valid_modes <- c("legacy", "mcnemar", "regression")
if (!MODE %in% valid_modes) {
  stop("Invalid --mode '", MODE, "'; must be one of: ",
       paste(valid_modes, collapse = ", "), call. = FALSE)
}

# Replicate-agg modes use a distinct cache prefix so they never collide with
# replicates_v4_cell*.rds. BASE_SEED is also bumped so any cached legacy rep
# seeds don't accidentally overlap.
if (MODE != "legacy") {
  BENCHMARK <- paste0("replagg_", MODE, "_v1")
  BASE_SEED <- 20260423L
}

# Extra-n supplement mode
EXTRA_SIZES  <- if (!is.null(EXTRA_N)) as.integer(strsplit(EXTRA_N, ",")[[1]]) else NULL
EXTRA_OFFSET <- 10000L

# Phase B8: weaker outlier magnitudes
PHASE_B_OFFSET <- 20000L

if (!dir.exists(CACHE_DIR)) dir.create(CACHE_DIR, recursive = TRUE)

# Purge stale cache when requested (prevents reuse of results from
# earlier code versions — see B6 report Anomaly 1)
if (CLEAR_CACHE) {
  old_files <- list.files(CACHE_DIR, pattern = paste0("^", BENCHMARK, "_cell"),
                          full.names = TRUE)
  if (length(old_files) > 0) {
    file.remove(old_files)
    if (VERBOSE) message("Cleared ", length(old_files), " cached cell files.")
  }
  # Also remove stale summary
  old_summary <- file.path(CACHE_DIR, paste0(BENCHMARK, "_summary.rds"))
  if (file.exists(old_summary)) file.remove(old_summary)
}


# ---- DGP Grid ---------------------------------------------------------------
#
# v2 grid additions (from B6 report recommendations):
#   - outlier_magnitude {2, 4}: tests QC sensitivity near the CV threshold
#     (2 SD outliers may escape the 25% CV criterion)
#   - lod_quantile {0, 0.15}: isolates replicate-handling effects from
#     censoring-induced bias
#   - re_intercept_sd {0, 0.5}: tests whether the RE method gains meaningful
#     efficiency when true between-subject variance exists
#
# When outlier_rep_rate = 0, outlier_magnitude is irrelevant.  We keep a
# single magnitude value (4) for those cells and filter the duplicate.

if (PHASE_B) {
  # ---- Phase B8: Weaker outlier magnitudes (0.5, 1.0, 1.5 SD) ---------------
  # Only non-zero outlier_rep_rate — magnitude is irrelevant when rate = 0.
  # 3 × 2 × 3 × 3 × 3 × 2 × 2 = 648 cells
  dgp_grid <- make_dgp_grid(
    replicate_sd      = c(0.05, 0.12, 0.20),
    outlier_rep_rate  = c(0.05, 0.15),
    n_subjects        = c(40L, 100L, 200L),
    effect_size       = c(0, 0.5, 1.0),
    outlier_magnitude = c(0.5, 1.0, 1.5),
    lod_quantile      = c(0, 0.15),
    re_intercept_sd   = c(0, 0.5)
  )

  # Cell IDs offset by PHASE_B_OFFSET
  dgp_grid$cell_id <- seq_len(nrow(dgp_grid)) + PHASE_B_OFFSET

} else {
  # ---- Main grid (540 cells) -------------------------------------------------
  dgp_grid <- make_dgp_grid(
    replicate_sd      = c(0.05, 0.12, 0.20),
    outlier_rep_rate  = c(0, 0.05, 0.15),
    n_subjects        = if (!is.null(EXTRA_SIZES)) EXTRA_SIZES else c(40L, 100L, 200L),
    effect_size       = c(0, 0.5, 1.0),
    outlier_magnitude = c(2, 4),
    lod_quantile      = c(0, 0.15),
    re_intercept_sd   = c(0, 0.5)
  )

  # Drop redundant rows where magnitude varies but outlier_rate = 0
  dgp_grid <- dgp_grid %>%
    filter(!(outlier_rep_rate == 0 & outlier_magnitude != 4))

  # Re-number cell IDs after filtering
  dgp_grid$cell_id <- seq_len(nrow(dgp_grid))
  if (!is.null(EXTRA_SIZES)) dgp_grid$cell_id <- dgp_grid$cell_id + EXTRA_OFFSET
}

if (VERBOSE) {
  phase_label <- if (PHASE_B) " [Phase B8: weak outliers]" else ""
  message("=== Benchmark 6: Replicate Handling Strategy (v4)", phase_label, " ===")
  message("DGP grid: ", nrow(dgp_grid), " cells x ", N_REPS, " reps = ",
          nrow(dgp_grid) * N_REPS, " total runs")
  message("Cache dir: ", CACHE_DIR)
  if (CLEAR_CACHE) message("Cache cleared (--clear_cache)")
}


# ---- Outlier Injection Helper -----------------------------------------------

#' Inject outlier replicates into simulated data (post-hoc)
#'
#' Following B3's pattern: after simulation, with probability outlier_rep_rate,
#' one replicate per pair is shifted by +/- outlier_magnitude * SD(value).
#'
#' @param d                data.frame with replicate-level data.
#' @param outlier_rep_rate Fraction of subject-analyte pairs to corrupt.
#' @param outlier_magnitude Number of SDs to shift the corrupted replicate.
#' @param seed             RNG seed for reproducibility.
#' @return Modified data.frame with outlier replicates injected.
inject_outlier_reps <- function(d, outlier_rep_rate, outlier_magnitude = 4,
                                seed) {
  set.seed(seed)
  # Identify replicate pairs
  pairs <- unique(d[, c("subject_id", "cytokine")])
  n_pairs <- nrow(pairs)
  n_outliers <- round(n_pairs * outlier_rep_rate)

  if (n_outliers == 0) return(list(data = d, outlier_pairs = character(0)))

  # Select pairs to corrupt
  corrupt_idx <- sample(n_pairs, n_outliers)
  global_sd <- sd(d$value, na.rm = TRUE)

  outlier_pair_ids <- character(n_outliers)

  for (i in seq_along(corrupt_idx)) {
    ci <- corrupt_idx[i]
    mask <- d$subject_id == pairs$subject_id[ci] &
            d$cytokine == pairs$cytokine[ci]
    reps <- which(mask)
    # Pick one replicate to corrupt, shift by +/- magnitude * SD
    target_rep <- sample(reps, 1)
    direction <- sample(c(-1, 1), 1)
    d$value[target_rep] <- d$value[target_rep] +
                           direction * outlier_magnitude * global_sd
    d$value_raw[target_rep] <- exp(d$value[target_rep])

    outlier_pair_ids[i] <- paste0(pairs$subject_id[ci], "_",
                                  pairs$cytokine[ci])
  }

  list(data = d, outlier_pairs = outlier_pair_ids)
}


# ---- Coefficient Extraction Helpers ------------------------------------------

.extract_glmmtmb <- function(fit_obj, method_name) {
  fail_row <- data.frame(
    method = method_name, estimate = NA_real_, se = NA_real_,
    ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
    converged = FALSE, time_s = NA_real_, stringsAsFactors = FALSE
  )

  converged <- FALSE
  model <- NULL

  if (is.list(fit_obj) && "converged" %in% names(fit_obj)) {
    converged <- fit_obj$converged
    model <- fit_obj$model
  } else if (inherits(fit_obj, "glmmTMB")) {
    converged <- TRUE
    model <- fit_obj
  }

  if (!converged || is.null(model)) return(fail_row)

  coef_tab <- tryCatch(
    summary(model)$coefficients$cond,
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


# ---- Tobit Extraction Helper ------------------------------------------------

.extract_tobit <- function(fit_obj, method_name) {
  fail_row <- data.frame(
    method = method_name, estimate = NA_real_, se = NA_real_,
    ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
    converged = FALSE, time_s = NA_real_, stringsAsFactors = FALSE
  )

  if (!is.list(fit_obj) || !isTRUE(fit_obj$converged)) return(fail_row)

  model <- fit_obj$model
  if (!inherits(model, "survreg")) return(fail_row)

  coef_tab <- tryCatch(summary(model)$table, error = function(e) NULL)
  if (is.null(coef_tab)) return(fail_row)

  group_rows <- grep("^group", rownames(coef_tab))
  if (length(group_rows) == 0) return(fail_row)

  est  <- coef_tab[group_rows[1], "Value"]
  se   <- coef_tab[group_rows[1], "Std. Error"]
  pval <- coef_tab[group_rows[1], "p"]

  data.frame(
    method = method_name, estimate = est, se = se,
    ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se,
    p_value = pval, converged = TRUE, time_s = NA_real_,
    stringsAsFactors = FALSE
  )
}


# ---- QC + Average Helper ----------------------------------------------------

#' Apply flag_replicate_outliers() then average remaining replicates
#'
#' @param d           Replicate-level data.frame with value on raw scale in
#'                    value_raw and log scale in value.
#' @param action_type One of "drop_both", "drop_farther", "winsorize".
#' @param method_name Label for the method.
#' @return One-row data.frame with standard benchmark columns.
fit_qc_average <- function(d, action_type, method_name) {
  t0 <- proc.time()["elapsed"]

  # flag_replicate_outliers expects raw-scale values
  qc_result <- tryCatch(
    flag_replicate_outliers(
      data           = d,
      sample_id_col  = "subject_id",
      cytokine_col   = "cytokine",
      value_col      = "value_raw",
      cv_threshold   = 25,
      flag_rule      = "cv",
      action         = action_type,
      verbose        = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(qc_result)) {
    elapsed <- proc.time()["elapsed"] - t0
    return(data.frame(
      method = method_name, estimate = NA_real_, se = NA_real_,
      ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
      converged = FALSE, time_s = elapsed, stringsAsFactors = FALSE
    ))
  }

  # After QC: recalculate log values from (potentially modified) raw values
  qc_result$value <- log(qc_result$value_raw)

  # Filter dropped rows and average
  if ("replicate_action" %in% names(qc_result)) {
    qc_result <- qc_result[qc_result$replicate_action != "dropped", ]
  }

  d_agg <- qc_result %>%
    group_by(subject_id, cytokine, group) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

  # Fit Gaussian (no RE) on aggregated data
  fit <- tryCatch(
    stats::lm(value ~ group, data = d_agg),
    error = function(e) NULL
  )
  elapsed <- proc.time()["elapsed"] - t0

  if (is.null(fit)) {
    return(data.frame(
      method = method_name, estimate = NA_real_, se = NA_real_,
      ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
      converged = FALSE, time_s = elapsed, stringsAsFactors = FALSE
    ))
  }

  r <- .extract_lm(fit, method_name)
  r$time_s <- elapsed
  r
}


# ---- QC Detection Rate Helper -----------------------------------------------

#' Compute QC detection rate for injected outlier replicates
#'
#' @param d              Replicate-level data.frame.
#' @param outlier_pairs  Character vector of "subject_id_cytokine" pairs that
#'                       were corrupted.
#' @return Fraction of corrupted pairs that were flagged by QC.
compute_qc_detection <- function(d, outlier_pairs) {
  if (length(outlier_pairs) == 0) return(NA_real_)

  qc_result <- tryCatch(
    flag_replicate_outliers(
      data           = d,
      sample_id_col  = "subject_id",
      cytokine_col   = "cytokine",
      value_col      = "value_raw",
      cv_threshold   = 25,
      flag_rule      = "cv",
      action         = "flag",
      verbose        = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(qc_result)) return(NA_real_)

  # For each corrupted pair, check if it was flagged
  flagged_pairs <- qc_result %>%
    filter(replicate_discordant) %>%
    mutate(pair_id = paste0(subject_id, "_", cytokine)) %>%
    pull(pair_id) %>%
    unique()

  sum(outlier_pairs %in% flagged_pairs) / length(outlier_pairs)
}


# ---- Single-Replication Function --------------------------------------------

#' Run one replication for a single DGP cell
#'
#' Generates data with replicates, injects outliers, applies all 7 methods,
#' returns one row per method.
#'
#' @param seed              RNG seed.
#' @param n_subjects        Sample size.
#' @param effect_size       True group effect (log-scale); 0 = null.
#' @param replicate_sd      Replicate-to-replicate SD (log-scale).
#' @param outlier_rep_rate  Fraction of pairs to corrupt with outliers.
#' @param outlier_magnitude Number of SDs to shift corrupted replicates.
#' @param lod_quantile      LOD quantile (0 = no censoring).
#' @param re_intercept_sd   Between-subject random intercept SD.
#' @return data.frame with 9 rows (one per method).
run_one_rep <- function(seed, n_subjects, effect_size, replicate_sd,
                        outlier_rep_rate, outlier_magnitude = 4,
                        lod_quantile = 0.15, re_intercept_sd = 0) {

  # --- Simulate data with 2 replicates ---
  sim <- simulate_immunoassay(
    n_subjects       = n_subjects,
    n_timepoints     = 1L,
    n_analytes       = 1L,
    design           = "cross_sectional",
    group_levels     = c("control", "treatment"),
    group_effects    = effect_size,
    signal_analytes  = if (effect_size != 0) 1L else 0L,
    effect_direction = "up",
    re_intercept_sd  = re_intercept_sd,
    residual_sd      = 1.0,
    n_replicates     = 2L,
    replicate_sd     = replicate_sd,
    well_failure_rate = 0,
    lod_quantile     = lod_quantile,
    seed             = seed
  )

  d  <- as.data.frame(sim$data[!is.na(sim$data$value), ])
  tr <- as.data.frame(sim$truth[!is.na(sim$data$value), ])

  # --- Post-hoc outlier injection ---
  injection <- inject_outlier_reps(
    d, outlier_rep_rate = outlier_rep_rate,
    outlier_magnitude = outlier_magnitude, seed = seed + 999999L
  )
  d <- injection$data
  outlier_pairs <- injection$outlier_pairs

  # --- QC detection rate (computed once, reported with all methods) ---
  qc_detection <- compute_qc_detection(d, outlier_pairs)

  results_list <- vector("list", 9)

  # --- Method 1: Naive average (no QC) ---
  t0 <- proc.time()["elapsed"]
  d_naive <- d %>%
    group_by(subject_id, cytokine, group) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
  fit_naive <- tryCatch(
    stats::lm(value ~ group, data = d_naive),
    error = function(e) NULL
  )
  time_naive <- proc.time()["elapsed"] - t0
  if (!is.null(fit_naive)) {
    r <- .extract_lm(fit_naive, "Naive average")
    r$time_s <- time_naive
  } else {
    r <- data.frame(
      method = "Naive average", estimate = NA_real_, se = NA_real_,
      ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
      converged = FALSE, time_s = time_naive, stringsAsFactors = FALSE
    )
  }
  results_list[[1]] <- r

  # --- Methods 2-4: QC-filtered averaging ---
  results_list[[2]] <- fit_qc_average(d, "drop_both",    "QC drop_both")
  results_list[[3]] <- fit_qc_average(d, "drop_farther", "QC drop_farther")
  results_list[[4]] <- fit_qc_average(d, "winsorize",    "QC winsorize")

  # --- Method 5: Random effects (no QC) ---
  # Use glmmTMB directly with (1|subject_id), NOT via fit_one()
  t0 <- proc.time()["elapsed"]
  re_fit <- tryCatch({
    fit <- glmmTMB::glmmTMB(
      value ~ group + (1 | subject_id),
      data   = d,
      family = gaussian()
    )
    fit
  }, error = function(e) NULL)
  time_re <- proc.time()["elapsed"] - t0

  if (!is.null(re_fit)) {
    r <- .extract_glmmtmb(re_fit, "Random effects")
    r$time_s <- time_re
  } else {
    r <- data.frame(
      method = "Random effects", estimate = NA_real_, se = NA_real_,
      ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
      converged = FALSE, time_s = time_re, stringsAsFactors = FALSE
    )
  }
  results_list[[5]] <- r

  # --- Method 6: Ignore structure (treat replicates as independent) ---
  t0 <- proc.time()["elapsed"]
  fit_ignore <- tryCatch(
    stats::lm(value ~ group, data = d),
    error = function(e) NULL
  )
  time_ignore <- proc.time()["elapsed"] - t0
  if (!is.null(fit_ignore)) {
    r <- .extract_lm(fit_ignore, "Ignore structure")
    r$time_s <- time_ignore
  } else {
    r <- data.frame(
      method = "Ignore structure", estimate = NA_real_, se = NA_real_,
      ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
      converged = FALSE, time_s = time_ignore, stringsAsFactors = FALSE
    )
  }
  results_list[[6]] <- r

  # --- Method 7: Oracle average (mean of true latent, no noise/outliers) ---
  t0 <- proc.time()["elapsed"]
  d_oracle <- data.frame(
    value      = tr$true_latent,
    subject_id = d$subject_id,
    group      = d$group
  )
  d_oracle_agg <- d_oracle %>%
    group_by(subject_id, group) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
  fit_oracle <- tryCatch(
    stats::lm(value ~ group, data = d_oracle_agg),
    error = function(e) NULL
  )
  time_oracle <- proc.time()["elapsed"] - t0
  if (!is.null(fit_oracle)) {
    r <- .extract_lm(fit_oracle, "Oracle average")
    r$time_s <- time_oracle
  } else {
    r <- data.frame(
      method = "Oracle average", estimate = NA_real_, se = NA_real_,
      ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
      converged = FALSE, time_s = time_oracle, stringsAsFactors = FALSE
    )
  }
  results_list[[7]] <- r

  # --- Method 8: QC + Tobit ---
  # Apply QC drop_both to replicates, average surviving replicates per subject,
  # then fit via fit_one(family="tobit") instead of lm().
  t0 <- proc.time()["elapsed"]
  qc_tobit_result <- tryCatch({
    qc_out <- flag_replicate_outliers(
      data           = d,
      sample_id_col  = "subject_id",
      cytokine_col   = "cytokine",
      value_col      = "value_raw",
      cv_threshold   = 25,
      flag_rule      = "cv",
      action         = "drop_both",
      verbose        = FALSE
    )
    qc_out$value <- log(qc_out$value_raw)
    if ("replicate_action" %in% names(qc_out)) {
      qc_out <- qc_out[qc_out$replicate_action != "dropped", ]
    }
    # Average replicates per subject, preserving censoring columns
    d_qc_agg <- qc_out %>%
      group_by(subject_id, cytokine, group) %>%
      summarise(
        value    = mean(value, na.rm = TRUE),
        cens_lod  = any(cens_lod),
        cens_ulod = any(cens_ulod),
        lod       = first(lod),
        ulod      = first(ulod),
        .groups   = "drop"
      )
    fit_one(
      dat    = as.data.frame(d_qc_agg),
      family = "tobit",
      fixed  = "group",
      random = NULL
    )
  }, error = function(e) NULL)
  time_qc_tobit <- proc.time()["elapsed"] - t0

  if (!is.null(qc_tobit_result)) {
    r <- .extract_tobit(qc_tobit_result, "QC + Tobit")
    r$time_s <- time_qc_tobit
  } else {
    r <- data.frame(
      method = "QC + Tobit", estimate = NA_real_, se = NA_real_,
      ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
      converged = FALSE, time_s = time_qc_tobit, stringsAsFactors = FALSE
    )
  }
  results_list[[8]] <- r

  # --- Method 9: Avg + Tobit ---
  # Average replicates per subject (like naive average), then fit via
  # fit_one(family="tobit") on subject-level data.
  t0 <- proc.time()["elapsed"]
  avg_tobit_result <- tryCatch({
    d_avg_agg <- d %>%
      group_by(subject_id, cytokine, group) %>%
      summarise(
        value     = mean(value, na.rm = TRUE),
        cens_lod  = any(cens_lod),
        cens_ulod = any(cens_ulod),
        lod       = first(lod),
        ulod      = first(ulod),
        .groups   = "drop"
      )
    fit_one(
      dat    = as.data.frame(d_avg_agg),
      family = "tobit",
      fixed  = "group",
      random = NULL
    )
  }, error = function(e) NULL)
  time_avg_tobit <- proc.time()["elapsed"] - t0

  if (!is.null(avg_tobit_result)) {
    r <- .extract_tobit(avg_tobit_result, "Avg + Tobit")
    r$time_s <- time_avg_tobit
  } else {
    r <- data.frame(
      method = "Avg + Tobit", estimate = NA_real_, se = NA_real_,
      ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
      converged = FALSE, time_s = time_avg_tobit, stringsAsFactors = FALSE
    )
  }
  results_list[[9]] <- r

  # --- Combine ---
  results <- do.call(rbind, results_list)
  results$true_effect   <- effect_size
  results$qc_detection  <- qc_detection
  results
}


# ---- Run Benchmark ----------------------------------------------------------

run_benchmark_replicates <- function(dgp_grid, n_reps, base_seed,
                                     cache_dir, n_cores, verbose) {
  all_summaries <- list()

  for (i in seq_len(nrow(dgp_grid))) {
    row <- dgp_grid[i, ]
    cell_label <- sprintf(
      "cell_%d [rep_sd=%.2f/out_rate=%.2f/out_mag=%g/n=%d/eff=%.1f/lod=%.2f/re=%.1f]",
      row$cell_id, row$replicate_sd, row$outlier_rep_rate,
      row$outlier_magnitude, row$n_subjects, row$effect_size,
      row$lod_quantile, row$re_intercept_sd
    )

    rep_fn <- function(seed) {
      run_one_rep(
        seed              = seed,
        n_subjects        = row$n_subjects,
        effect_size       = row$effect_size,
        replicate_sd      = row$replicate_sd,
        outlier_rep_rate  = row$outlier_rep_rate,
        outlier_magnitude = row$outlier_magnitude,
        lod_quantile      = row$lod_quantile,
        re_intercept_sd   = row$re_intercept_sd
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

    # Aggregate per method
    methods <- unique(raw_results$method[!is.na(raw_results$method)])
    for (m in methods) {
      m_reps <- raw_results[!is.na(raw_results$method) &
                            raw_results$method == m, ]
      agg <- aggregate_cell(m_reps, truth = row$effect_size)

      # Timing summary
      agg$mean_time_s <- mean(m_reps$time_s, na.rm = TRUE)

      # CI width for relative efficiency analysis
      if (all(c("ci_lo", "ci_hi") %in% names(m_reps))) {
        ciw <- calc_ci_width(m_reps$ci_lo, m_reps$ci_hi)
        agg$ci_width      <- ciw$estimate
        agg$ci_width_mcse <- ciw$mcse
      }

      # QC detection rate (same for all methods in a cell, but tracked here)
      if ("qc_detection" %in% names(m_reps)) {
        qc_vals <- m_reps$qc_detection[is.finite(m_reps$qc_detection)]
        agg$qc_detection      <- if (length(qc_vals) > 0) mean(qc_vals) else NA_real_
        agg$qc_detection_mcse <- if (length(qc_vals) > 1) {
          sd(qc_vals) / sqrt(length(qc_vals))
        } else NA_real_
      }

      # Cell identifiers
      agg$method            <- m
      agg$replicate_sd      <- row$replicate_sd
      agg$outlier_rep_rate  <- row$outlier_rep_rate
      agg$outlier_magnitude <- row$outlier_magnitude
      agg$n_subjects        <- row$n_subjects
      agg$effect_size       <- row$effect_size
      agg$lod_quantile      <- row$lod_quantile
      agg$re_intercept_sd   <- row$re_intercept_sd
      agg$cell_id           <- row$cell_id
      all_summaries[[length(all_summaries) + 1]] <- agg
    }
  }

  bind_rows(all_summaries)
}


# =============================================================================
# ---- Replicate-Aggregation Scaffold (MODE = "mcnemar" | "regression") -------
# =============================================================================
#
# This block short-circuits out before the legacy (B6 v4) Execute section
# below. Scope of this chat: DGP + harness only -- no local run, no HPC
# submit, no report. Report + figure generation live in a follow-up chat.

if (MODE %in% c("mcnemar", "regression")) {

  # ---- DGP grid -----------------------------------------------------------
  if (MODE == "mcnemar") {
    replagg_grid <- make_dgp_grid(
      n_subjects         = c(30L, 60L, 120L),
      n_replicates       = c(2L, 3L),
      sigma_rep          = c(0.20, 0.50),
      baseline_detection = c(0.30, 0.50),
      delta_rate         = c(0, 0.15),
      analyte_rho        = 0.30,
      rho_t              = c(0.30, 0.50)
    )
  } else {
    replagg_grid <- make_dgp_grid(
      n_subjects   = c(30L, 60L, 120L),
      n_replicates = c(2L, 3L),
      sigma_rep    = c(0.20, 0.50),
      beta_t       = c(0, 0.50),
      analyte_rho  = 0.30,
      rho_t        = c(0.30, 0.50),
      residual_sd  = 1.0
    )
  }

  if (VERBOSE) {
    message("=== Replicate-aggregation benchmark (mode=", MODE, ") ===")
    message("DGP grid: ", nrow(replagg_grid), " cells x ", N_REPS, " reps = ",
            nrow(replagg_grid) * N_REPS, " total runs")
    message("Analytes per sim: ", N_ANALYTES, " (", N_SIGNAL, " signal)")
    message("Cache prefix: ", BENCHMARK, "_cell*.rds")
  }

  # ---- Shared latent generator (mirrors benchmark_mcnemar.R::run_one_rep) -
  # Kept inline rather than hoisted into helpers_benchmark.R to avoid a
  # cross-cutting change; matches B2 block-for-block so future maintainers
  # who fix one can grep and find the other.
  .replagg_make_latent <- function(n_subjects, n_analytes, analyte_rho, rho_t) {
    stopifnot(abs(rho_t) < 1,
              analyte_rho > -1 / (n_analytes - 1),
              analyte_rho < 1)
    cor_mat <- diag(n_analytes)
    if (analyte_rho != 0) {
      cor_mat[lower.tri(cor_mat)] <- analyte_rho
      cor_mat[upper.tri(cor_mat)] <- analyte_rho
    }
    chol_mat <- chol(cor_mat)
    Z1_raw   <- matrix(stats::rnorm(n_subjects * n_analytes),
                       n_subjects, n_analytes)
    Z1       <- Z1_raw %*% chol_mat
    Z2_innov <- matrix(stats::rnorm(n_subjects * n_analytes),
                       n_subjects, n_analytes)
    Z2_innov <- Z2_innov %*% chol_mat
    Z2       <- rho_t * Z1 + sqrt(1 - rho_t^2) * Z2_innov
    list(Z1 = Z1, Z2 = Z2)
  }

  # Expand a [subjects, analytes] latent matrix into a
  # [subjects, analytes, reps] array with iid N(0, sigma_rep^2) noise
  # added at the rep level (epsilon_r in the DGP equation).
  .replagg_expand_reps <- function(Z, n_replicates, sigma_rep) {
    n_s <- nrow(Z); n_a <- ncol(Z)
    arr <- array(0.0, dim = c(n_s, n_a, n_replicates))
    for (r in seq_len(n_replicates)) {
      noise <- matrix(stats::rnorm(n_s * n_a, sd = sigma_rep), n_s, n_a)
      arr[, , r] <- Z + noise
    }
    arr
  }

  .replagg_analyte_names <- function(n_analytes) {
    paste0("A", sprintf("%02d", seq_len(n_analytes)))
  }

  # ---- run_one_rep: mcnemar mode ------------------------------------------
  run_one_rep_replagg_mcnemar <- function(seed, n_subjects, baseline_detection,
                                          delta_rate, analyte_rho, rho_t,
                                          n_replicates, sigma_rep) {
    tryCatch({
      set.seed(seed)

      analyte_names <- .replagg_analyte_names(N_ANALYTES)
      signal_mask   <- rep(FALSE, N_ANALYTES)
      signal_mask[sort(sample.int(N_ANALYTES, N_SIGNAL))] <- TRUE

      t1       <- stats::qnorm(1 - baseline_detection)
      p2       <- min(baseline_detection + delta_rate, 0.99)
      t2_sig   <- stats::qnorm(1 - p2)
      true_dp  <- (p2 - baseline_detection) * 100

      Z   <- .replagg_make_latent(n_subjects, N_ANALYTES, analyte_rho, rho_t)
      Z1r <- .replagg_expand_reps(Z$Z1, n_replicates, sigma_rep)
      Z2r <- .replagg_expand_reps(Z$Z2, n_replicates, sigma_rep)

      det_T1 <- Z1r > t1
      det_T2 <- array(FALSE, dim = dim(Z2r))
      for (k in seq_len(N_ANALYTES)) {
        thr_k <- if (signal_mask[k]) t2_sig else t1
        det_T2[, k, ] <- Z2r[, k, ] > thr_k
      }

      n_rows_per_tp <- n_subjects * N_ANALYTES * n_replicates
      long_data <- data.frame(
        subject_id = rep(seq_len(n_subjects),
                         times = N_ANALYTES * n_replicates * 2L),
        cytokine   = rep(rep(analyte_names, each = n_subjects),
                         times = n_replicates * 2L),
        rep_id     = rep(rep(seq_len(n_replicates),
                             each = n_subjects * N_ANALYTES), times = 2L),
        timepoint  = rep(c("T1", "T2"), each = n_rows_per_tp),
        cens_lod   = c(as.vector(!det_T1), as.vector(!det_T2)),
        stringsAsFactors = FALSE
      )

      rules <- c("majority_vote", "any_detected", "all_detected")
      per_rule <- lapply(rules, function(rule) {
        res <- tryCatch(
          mcnemar_detection(
            data          = long_data,
            subject_col   = "subject_id",
            cytokine_col  = "cytokine",
            timepoint_col = "timepoint",
            censoring_col = "cens_lod",
            baseline      = "T1",
            comparison    = "T2",
            replicate_agg = rule,
            method        = "midp",
            fdr_method    = "BH",
            quiet         = TRUE
          ),
          error = function(e) NULL
        )
        if (is.null(res)) {
          return(data.frame(
            method       = paste0("replicate_agg: ", rule),
            analyte      = analyte_names,
            is_signal    = signal_mask,
            p_value      = NA_real_, q_value = NA_real_,
            delta        = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
            n_discordant = NA_integer_,
            true_delta   = ifelse(signal_mask, true_dp, 0),
            stringsAsFactors = FALSE
          ))
        }
        ki <- match(as.character(res$cytokine), analyte_names)
        q_col <- if ("q_mcnemar" %in% names(res)) res$q_mcnemar
                 else p.adjust(res$p_mcnemar, "BH")
        data.frame(
          method       = paste0("replicate_agg: ", rule),
          analyte      = as.character(res$cytokine),
          is_signal    = signal_mask[ki],
          p_value      = res$p_mcnemar,
          q_value      = q_col,
          delta        = res$delta_detection,
          ci_lo        = res$delta_ci_lo,
          ci_hi        = res$delta_ci_hi,
          n_discordant = as.integer(res$n_discordant),
          true_delta   = ifelse(signal_mask[ki], true_dp, 0),
          stringsAsFactors = FALSE
        )
      })
      do.call(rbind, per_rule)

    }, error = function(e) {
      data.frame(
        method = NA_character_, analyte = NA_character_,
        is_signal = NA, p_value = NA_real_, q_value = NA_real_,
        delta = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
        n_discordant = NA_integer_, true_delta = NA_real_,
        stringsAsFactors = FALSE
      )
    })
  }

  # ---- run_one_rep: regression mode ---------------------------------------
  # Small helper: extract the timepoint coefficient from an
  # immuno_model_set returned by fit_models(families = "gaussian", ...).
  .replagg_extract_tp <- function(fit_result, method_name, analyte) {
    fail <- data.frame(
      method = method_name, analyte = analyte,
      estimate = NA_real_, se = NA_real_,
      ci_lo = NA_real_, ci_hi = NA_real_,
      p_value = NA_real_, converged = FALSE,
      stringsAsFactors = FALSE
    )
    if (is.null(fit_result)) return(fail)
    fit_obj <- if (inherits(fit_result, "immuno_model_set"))
                 fit_result$models[["gaussian"]] else fit_result
    if (is.null(fit_obj) || !isTRUE(fit_obj$converged)) return(fail)
    model <- fit_obj$model
    coefs <- tryCatch(
      if (inherits(model, "glmmTMB")) summary(model)$coefficients$cond
      else summary(model)$coefficients,
      error = function(e) NULL
    )
    if (is.null(coefs)) return(fail)
    tp_rows <- grep("^timepoint", rownames(coefs))
    if (length(tp_rows) == 0) return(fail)
    est  <- coefs[tp_rows[1], "Estimate"]
    se   <- coefs[tp_rows[1], "Std. Error"]
    pval <- coefs[tp_rows[1], ncol(coefs)]
    data.frame(
      method = method_name, analyte = analyte,
      estimate = est, se = se,
      ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se,
      p_value = pval, converged = TRUE,
      stringsAsFactors = FALSE
    )
  }

  run_one_rep_replagg_regression <- function(seed, n_subjects, beta_t,
                                             analyte_rho, rho_t,
                                             n_replicates, sigma_rep,
                                             residual_sd) {
    tryCatch({
      set.seed(seed)

      analyte_names <- .replagg_analyte_names(N_ANALYTES)
      signal_mask   <- rep(FALSE, N_ANALYTES)
      signal_mask[sort(sample.int(N_ANALYTES, N_SIGNAL))] <- TRUE

      Z <- .replagg_make_latent(n_subjects, N_ANALYTES, analyte_rho, rho_t)
      Z2_shifted <- Z$Z2
      for (k in seq_len(N_ANALYTES)) {
        if (signal_mask[k]) Z2_shifted[, k] <- Z2_shifted[, k] + beta_t
      }

      Z1r <- .replagg_expand_reps(Z$Z1,       n_replicates, sigma_rep)
      Z2r <- .replagg_expand_reps(Z2_shifted, n_replicates, sigma_rep)

      # Additional iid residual on top of rep-level noise (separates the
      # subject-analyte-timepoint "true" level from the single-measurement
      # noise floor).
      Z1r <- Z1r + array(stats::rnorm(length(Z1r), sd = residual_sd),
                         dim = dim(Z1r))
      Z2r <- Z2r + array(stats::rnorm(length(Z2r), sd = residual_sd),
                         dim = dim(Z2r))

      n_rows_per_tp <- n_subjects * N_ANALYTES * n_replicates
      d <- data.frame(
        subject_id = rep(seq_len(n_subjects),
                         times = N_ANALYTES * n_replicates * 2L),
        cytokine   = rep(rep(analyte_names, each = n_subjects),
                         times = n_replicates * 2L),
        rep_id     = rep(rep(seq_len(n_replicates),
                             each = n_subjects * N_ANALYTES), times = 2L),
        timepoint  = factor(rep(c("T1", "T2"), each = n_rows_per_tp),
                            levels = c("T1", "T2")),
        value      = c(as.vector(Z1r), as.vector(Z2r)),
        cens_lod   = FALSE,
        cens_ulod  = FALSE,
        lod        = -Inf,
        ulod       = Inf,
        stringsAsFactors = FALSE
      )

      per_analyte <- vector("list", N_ANALYTES * 2L)
      idx <- 0L
      for (k in seq_len(N_ANALYTES)) {
        ak  <- analyte_names[k]
        d_k <- d[d$cytokine == ak, , drop = FALSE]

        # Method A: rep-level fit with rep_col RE append
        fit_rep <- tryCatch(
          fit_models(
            dat      = d_k,
            families = "gaussian",
            fixed    = "timepoint",
            random   = "(1|subject_id)",
            rep_col  = "rep_id",
            quiet    = TRUE
          ),
          error = function(e) NULL
        )
        idx <- idx + 1L
        per_analyte[[idx]] <- .replagg_extract_tp(fit_rep,
                                                  "fit_models + rep RE", ak)

        # Method B: pre-aggregate replicates via mean, then fit without rep_col
        d_k_agg <- tryCatch(
          aggregate_replicates(
            data          = d_k,
            rep_col       = "rep_id",
            subject_col   = "subject_id",
            cytokine_col  = "cytokine",
            timepoint_col = "timepoint",
            value_col     = "value",
            rule          = "mean"
          ),
          error = function(e) NULL
        )
        fit_agg <- NULL
        if (!is.null(d_k_agg)) {
          # aggregate_replicates may strip cens_* columns; fit_models needs
          # them. Re-add neutral defaults (scaffold defers censoring).
          for (col in c("cens_lod", "cens_ulod")) {
            if (!col %in% names(d_k_agg)) d_k_agg[[col]] <- FALSE
          }
          if (!"lod"  %in% names(d_k_agg)) d_k_agg$lod  <- -Inf
          if (!"ulod" %in% names(d_k_agg)) d_k_agg$ulod <- Inf
          fit_agg <- tryCatch(
            fit_models(
              dat      = d_k_agg,
              families = "gaussian",
              fixed    = "timepoint",
              random   = "(1|subject_id)",
              quiet    = TRUE
            ),
            error = function(e) NULL
          )
        }
        idx <- idx + 1L
        per_analyte[[idx]] <- .replagg_extract_tp(
          fit_agg, "aggregate_replicates + fit_models", ak
        )
      }

      res <- do.call(rbind, per_analyte)
      ki <- match(res$analyte, analyte_names)
      res$is_signal   <- signal_mask[ki]
      res$true_effect <- ifelse(res$is_signal, beta_t, 0)
      res

    }, error = function(e) {
      data.frame(
        method = NA_character_, analyte = NA_character_,
        estimate = NA_real_, se = NA_real_,
        ci_lo = NA_real_, ci_hi = NA_real_,
        p_value = NA_real_, converged = FALSE,
        is_signal = NA, true_effect = NA_real_,
        stringsAsFactors = FALSE
      )
    })
  }

  # ---- Cell-level runner ---------------------------------------------------
  run_benchmark_replagg <- function(mode, dgp_grid, n_reps, base_seed,
                                    cache_dir, n_cores, verbose) {
    all_raw <- list()
    for (i in seq_len(nrow(dgp_grid))) {
      row <- dgp_grid[i, ]
      cell_label <- if (mode == "mcnemar") {
        sprintf(paste0("replagg_mcnemar_cell_%d ",
                       "[n=%d/reps=%d/sigma=%.2f/base=%.2f/delta=%.2f/rho_t=%.2f]"),
                row$cell_id, row$n_subjects, row$n_replicates, row$sigma_rep,
                row$baseline_detection, row$delta_rate, row$rho_t)
      } else {
        sprintf(paste0("replagg_regression_cell_%d ",
                       "[n=%d/reps=%d/sigma=%.2f/beta_t=%.2f/rho_t=%.2f]"),
                row$cell_id, row$n_subjects, row$n_replicates, row$sigma_rep,
                row$beta_t, row$rho_t)
      }

      rep_fn <- local({
        r <- row
        force(r)
        if (mode == "mcnemar") {
          function(seed) {
            run_one_rep_replagg_mcnemar(
              seed               = seed,
              n_subjects         = r$n_subjects,
              baseline_detection = r$baseline_detection,
              delta_rate         = r$delta_rate,
              analyte_rho        = r$analyte_rho,
              rho_t              = r$rho_t,
              n_replicates       = r$n_replicates,
              sigma_rep          = r$sigma_rep
            )
          }
        } else {
          function(seed) {
            run_one_rep_replagg_regression(
              seed         = seed,
              n_subjects   = r$n_subjects,
              beta_t       = r$beta_t,
              analyte_rho  = r$analyte_rho,
              rho_t        = r$rho_t,
              n_replicates = r$n_replicates,
              sigma_rep    = r$sigma_rep,
              residual_sd  = r$residual_sd
            )
          }
        }
      })

      cache_file <- file.path(cache_dir,
                              paste0(BENCHMARK, "_cell", row$cell_id, ".rds"))

      raw <- run_benchmark_cell(
        rep_fn     = rep_fn,
        n_reps     = n_reps,
        base_seed  = base_seed + (i - 1L) * 10000L,
        cache_file = cache_file,
        n_cores    = n_cores,
        cell_label = cell_label,
        verbose    = verbose
      )

      for (nm in names(row)) raw[[nm]] <- row[[nm]]
      all_raw[[i]] <- raw
    }
    bind_rows(all_raw)
  }

  # ---- Execute + minimal aggregation (report is a follow-up chat) ---------
  if (VERBOSE) message("\nStarting ", MODE, " execution...")

  all_results <- run_benchmark_replagg(
    mode      = MODE,
    dgp_grid  = replagg_grid,
    n_reps    = N_REPS,
    base_seed = BASE_SEED,
    cache_dir = CACHE_DIR,
    n_cores   = N_CORES,
    verbose   = VERBOSE
  )

  if (VERBOSE) message("\nAggregating ", MODE, " results...")

  ok <- all_results %>% dplyr::filter(!error, !is.na(method))

  if (MODE == "mcnemar") {
    cell_summary <- ok %>%
      dplyr::group_by(n_subjects, n_replicates, sigma_rep, baseline_detection,
                      delta_rate, rho_t, method, cell_id) %>%
      dplyr::summarise(
        n_reps     = length(unique(rep_id)),
        type1      = mean((p_value < 0.05)[!is_signal], na.rm = TRUE),
        type1_mcse = sd((p_value < 0.05)[!is_signal], na.rm = TRUE) /
                     sqrt(max(sum(!is.na(p_value) & !is_signal), 1L)),
        power      = mean((p_value < 0.05)[is_signal], na.rm = TRUE),
        power_mcse = sd((p_value < 0.05)[is_signal], na.rm = TRUE) /
                     sqrt(max(sum(!is.na(p_value) & is_signal), 1L)),
        .groups = "drop"
      )
  } else {
    cell_summary <- ok %>%
      dplyr::group_by(n_subjects, n_replicates, sigma_rep, beta_t, rho_t,
                      residual_sd, method, cell_id) %>%
      dplyr::summarise(
        n_reps        = length(unique(rep_id)),
        bias          = mean((estimate - true_effect)[is_signal], na.rm = TRUE),
        bias_mcse     = sd((estimate - true_effect)[is_signal], na.rm = TRUE) /
                        sqrt(max(sum(is.finite(estimate) & is_signal), 1L)),
        coverage      = mean(((ci_lo <= true_effect) &
                              (ci_hi >= true_effect))[is_signal], na.rm = TRUE),
        coverage_mcse = sd(((ci_lo <= true_effect) &
                            (ci_hi >= true_effect))[is_signal], na.rm = TRUE) /
                        sqrt(max(sum(is.finite(ci_lo) & is_signal), 1L)),
        type1         = mean((p_value < 0.05)[!is_signal], na.rm = TRUE),
        type1_mcse    = sd((p_value < 0.05)[!is_signal], na.rm = TRUE) /
                        sqrt(max(sum(!is.na(p_value) & !is_signal), 1L)),
        .groups = "drop"
      )
  }

  summary_file <- file.path(CACHE_DIR, paste0(BENCHMARK, "_summary.rds"))
  saveRDS(cell_summary, summary_file)
  if (VERBOSE) {
    message("\n=== Replicate-agg scaffold (mode=", MODE, ") complete ===")
    message("Summary -> ", summary_file)
    message("(Figures + report are produced by a separate follow-up chat.)")
  }

  quit(save = "no", status = 0)
}


# ---- Execute ----------------------------------------------------------------

if (VERBOSE) message("\nStarting benchmark execution...")

summary_df <- run_benchmark_replicates(
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

# Save final summary (merge with existing when running extra sizes or Phase B)
# Phase B merges into the MAIN summary file (replicates_v4_summary.rds) so
# plot scripts see all data together
summary_file <- file.path(CACHE_DIR,
                          paste0(if (PHASE_B) "replicates_v4" else BENCHMARK,
                                 "_summary.rds"))
if (!is.null(EXTRA_SIZES) || PHASE_B) {
  merge_with_existing_summary(summary_annotated, summary_file, verbose = VERBOSE)
} else {
  save_main_summary_preserving_supplements(summary_annotated, summary_file,
                                           verbose = VERBOSE)
}



# Method ordering for Formatted Summary Table
method_order <- c("Oracle average", "QC drop_farther", "QC drop_both",
                  "QC winsorize", "Naive average", "Random effects",
                  "Ignore structure", "QC + Tobit", "Avg + Tobit")
summary_annotated$method <- factor(summary_annotated$method,
                                    levels = method_order)

# Figures are generated by the standalone plot_B6.R script.
# Run: Rscript inst/simulations/plot_B6.R [--cache_dir PATH]


# ---- Formatted Summary Table ------------------------------------------------

summary_formatted <- format_summary_table(summary_annotated)

if (VERBOSE) {
  message("\n=== Benchmark 6: Summary (selected cells) ===\n")
  print(
    summary_formatted %>%
      filter(effect_size %in% c(0, 0.5),
             n_subjects == 100,
             outlier_magnitude == 4,
             lod_quantile == 0.15,
             re_intercept_sd == 0) %>%
      select(replicate_sd, outlier_rep_rate, n_subjects, effect_size,
             method, bias_fmt, coverage_fmt, rejection_fmt,
             starts_with("pass_")) %>%
      as.data.frame(),
    right = FALSE
  )
}

if (VERBOSE) {
  phase_label <- if (PHASE_B) " [Phase B8]" else ""
  message("\n=== Benchmark 6 complete (v4)", phase_label, " ===")
  if (exists("fig_dir")) message("Figures -> ", fig_dir)
  message("Summary -> ", summary_file)
}
