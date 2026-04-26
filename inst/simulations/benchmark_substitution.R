#!/usr/bin/env Rscript
# benchmark_substitution.R — Benchmark 5: LOD Substitution Method Comparison (v3)
#
# ADEMP Structure:
#   Aims:    Quantify how the choice of LOD substitution method affects bias,
#            coverage, and power for group-effect estimation across censoring
#            rates and residual noise levels. Compare substitution-based
#            approaches against Tobit. Extension: evaluate panel-level FDR
#            control and per-analyte bias under heterogeneous LODs across a
#            10-analyte panel.
#   DGP:     simulate_immunoassay() with clean assay (no plate/lot/replicate
#            effects). Single or multi-analyte, cross-sectional, 2 groups.
#            Grid varies censoring rate, sample size, effect size, and
#            residual SD. Panel mode uses 10 analytes with heterogeneous
#            LODs from the built-in cytokine library.
#   Estimand: Marginal group effect beta_group on log-concentration.
#   Methods: (1) LOD/2+Gaussian, (2) LOD/sqrt2+Gaussian, (3) halfmin+Gaussian,
#            (4) LOD+Gaussian, (5) zero+Gaussian, (6) LOD/2+Gamma,
#            (7) Tobit, (8) Oracle.
#   Performance: Single-analyte (90 cells): bias, coverage, RMSE, power,
#                type I error, convergence, pd-Hessian rate.
#                Panel (18 cells): FDR, FWER, sensitivity, per-analyte bias
#                vs. realized censoring.
#                Phase B7 (81 cells): mis-specified DGP robustness — t(df=5)
#                errors, LOD ±10% measurement error, or both.
#
# Usage:
#   Rscript benchmark_substitution.R [--n_reps N] [--n_cores N] [--cache_dir DIR]
#   Rscript benchmark_substitution.R --phase_b B7 [--n_reps N] [--n_cores N]

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

N_REPS     <- as.integer(parse_cli_arg("--n_reps",    "500"))
N_CORES    <- as.integer(parse_cli_arg("--n_cores",   "1"))
CACHE_DIR  <- parse_cli_arg("--cache_dir", file.path(script_dir, "cache"))
EXTRA_N    <- parse_cli_arg("--extra_n", NULL)
PHASE_B    <- parse_cli_arg("--phase_b", NULL)
BASE_SEED  <- 20250305L
BENCHMARK  <- "substitution_v3"
VERBOSE    <- TRUE

# Extra-n supplement mode
EXTRA_SIZES  <- if (!is.null(EXTRA_N)) as.integer(strsplit(EXTRA_N, ",")[[1]]) else NULL
EXTRA_OFFSET <- 10000L

# Phase B offset for mis-specified DGP sub-grids
PHASE_B_OFFSET <- 20000L

if (!dir.exists(CACHE_DIR)) dir.create(CACHE_DIR, recursive = TRUE)


# ---- DGP Grid & Main Execution (skipped when --phase_b) --------------------

if (is.null(PHASE_B)) {

dgp_grid <- make_dgp_grid(
  censoring_target = c(0.05, 0.15, 0.30, 0.50, 0.70),
  n_subjects       = if (!is.null(EXTRA_SIZES)) EXTRA_SIZES else c(40L, 100L, 200L),
  effect_size      = c(0, 0.5, 1.0),
  residual_sd      = c(0.5, 1.0)
)
if (!is.null(EXTRA_SIZES)) dgp_grid$cell_id <- dgp_grid$cell_id + EXTRA_OFFSET

if (VERBOSE) {
  message("=== Benchmark 5: LOD Substitution Method Comparison (v3) ===")
  message("DGP grid: ", nrow(dgp_grid), " cells (5 cens x 3 n x 3 eff x 2 sd = 90)")
  message("Reps: ", N_REPS, " | Total runs: ", nrow(dgp_grid) * N_REPS)
  message("Cache dir: ", CACHE_DIR)
}

}  # end if (is.null(PHASE_B)) — grid-setup guard; function defs below are shared


# ---- LOD Substitution Helper ------------------------------------------------

#' Apply LOD substitution to simulated data
#' @param d  data.frame from simulate_immunoassay with value, value_raw, lod, cens_lod.
#' @param method One of "half", "sqrt2", "halfmin", "lod", "zero".
#' @return Modified data.frame with substituted values.
substitute_lod <- function(d, method) {
  cens <- d$cens_lod
  lod  <- d$lod[cens]
  raw  <- d$value_raw

  if (method == "half")    raw[cens] <- lod / 2
  if (method == "sqrt2")   raw[cens] <- lod / sqrt(2)
  if (method == "halfmin") {
    obs_vals <- raw[!cens]
    hm <- if (length(obs_vals) > 0) min(obs_vals) / 2 else lod / 2
    raw[cens] <- hm
  }
  if (method == "lod")     raw[cens] <- lod
  if (method == "zero")    { d$value[cens] <- 0; return(d) }

  d$value_raw <- raw
  d$value     <- log(raw)
  d
}


# ---- Coefficient Extraction Helpers ------------------------------------------

.extract_survreg <- function(fit_obj, method_name) {
  fail_row <- data.frame(
    method = method_name, estimate = NA_real_, se = NA_real_,
    ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
    converged = FALSE, pd_hess = NA, time_s = NA_real_,
    stringsAsFactors = FALSE
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
    p_value = pval, converged = TRUE, pd_hess = NA,
    time_s = NA_real_, stringsAsFactors = FALSE
  )
}

.extract_glmmtmb <- function(fit_obj, method_name) {
  fail_row <- data.frame(
    method = method_name, estimate = NA_real_, se = NA_real_,
    ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
    converged = FALSE, pd_hess = NA, time_s = NA_real_,
    stringsAsFactors = FALSE
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
  pd   <- tryCatch(isTRUE(fit_obj$model$sdr$pdHess), error = function(e) NA)

  data.frame(
    method = method_name, estimate = est, se = se,
    ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se,
    p_value = pval, converged = TRUE, pd_hess = pd,
    time_s = NA_real_, stringsAsFactors = FALSE
  )
}

.extract_lm <- function(fit, method_name) {
  fail_row <- data.frame(
    method = method_name, estimate = NA_real_, se = NA_real_,
    ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
    converged = FALSE, pd_hess = NA, time_s = NA_real_,
    stringsAsFactors = FALSE
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
    p_value = pval, converged = TRUE, pd_hess = NA,
    time_s = NA_real_, stringsAsFactors = FALSE
  )
}


# ---- Per-Analyte Method Application -----------------------------------------

#' Apply all 8 substitution/modeling methods to a single-analyte data subset
#'
#' Extracted to support both single-analyte and multi-analyte (panel) modes.
#'
#' @param d_k  Per-analyte data.frame (value, value_raw, group, cens_lod, lod, etc.).
#' @param tr_k Per-analyte truth data.frame (true_latent column aligned with d_k).
#' @return data.frame with 8 rows (one per method): method, estimate, se,
#'   ci_lo, ci_hi, p_value, converged, pd_hess, time_s.
.fit_all_substitution_methods <- function(d_k, tr_k) {

  results_list <- vector("list", 8)

  # Helper: fit Gaussian on substituted data via fit_one
  fit_substituted_gaussian <- function(d_orig, sub_method, method_label) {
    t0 <- proc.time()["elapsed"]
    d_sub <- substitute_lod(d_orig, sub_method)
    d_sub$cens_lod  <- FALSE
    d_sub$cens_ulod <- FALSE
    fit <- tryCatch(
      fit_one(d_sub, family = "gaussian", fixed = "group", random = NULL),
      error = function(e) list(converged = FALSE, model = NULL)
    )
    elapsed <- proc.time()["elapsed"] - t0
    r <- .extract_glmmtmb(fit, method_label)
    r$time_s <- elapsed
    r
  }

  # --- Methods 1-5: Substitution + Gaussian ---
  results_list[[1]] <- fit_substituted_gaussian(d_k, "half",    "LOD/2 + Gaussian")
  results_list[[2]] <- fit_substituted_gaussian(d_k, "sqrt2",   "LOD/sqrt2 + Gaussian")
  results_list[[3]] <- fit_substituted_gaussian(d_k, "halfmin", "halfmin + Gaussian")
  results_list[[4]] <- fit_substituted_gaussian(d_k, "lod",     "LOD + Gaussian")
  results_list[[5]] <- fit_substituted_gaussian(d_k, "zero",    "zero + Gaussian")

  # --- Method 6: LOD/2 + Gamma (raw-scale) ---
  t0 <- proc.time()["elapsed"]
  d_gamma <- d_k
  d_gamma$value <- d_k$value_raw
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
  r <- .extract_glmmtmb(gamma_fit, "LOD/2 + Gamma")
  r$time_s <- time_gamma
  results_list[[6]] <- r

  # --- Method 7: Tobit (censoring-aware) ---
  t0 <- proc.time()["elapsed"]
  tobit_fit <- tryCatch(
    fit_one(d_k, family = "tobit", fixed = "group", random = NULL),
    error = function(e) list(converged = FALSE, model = NULL)
  )
  time_tobit <- proc.time()["elapsed"] - t0
  r <- .extract_survreg(tobit_fit, "Tobit")
  r$time_s <- time_tobit
  results_list[[7]] <- r

  # --- Method 8: Oracle (uncensored truth) ---
  t0 <- proc.time()["elapsed"]
  d_oracle <- data.frame(value = tr_k$true_latent, group = d_k$group)
  oracle_fit <- tryCatch(
    stats::lm(value ~ group, data = d_oracle),
    error = function(e) NULL
  )
  time_oracle <- proc.time()["elapsed"] - t0
  if (!is.null(oracle_fit)) {
    r <- .extract_lm(oracle_fit, "Oracle")
    r$time_s <- time_oracle
  } else {
    r <- data.frame(
      method = "Oracle", estimate = NA_real_, se = NA_real_,
      ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_,
      converged = FALSE, pd_hess = NA, time_s = time_oracle,
      stringsAsFactors = FALSE
    )
  }
  results_list[[8]] <- r

  do.call(rbind, results_list)
}


# ---- Single-Replication Function --------------------------------------------

#' Run one replication for a single DGP cell
#'
#' Generates data, applies all 8 methods per analyte, returns one row per
#' method (single-analyte) or one row per method x analyte (panel).
#'
#' @param seed          RNG seed.
#' @param n_subjects    Sample size.
#' @param effect_size   True group effect (log-scale); 0 = null.
#' @param censoring_target Target censoring fraction (0.05--0.70).
#' @param residual_sd   Within-subject log-scale SD (default 1.0).
#' @param n_analytes    Number of analytes (1 = single-analyte, >1 = panel).
#' @param n_signal      Number of signal analytes (default: 1 if effect != 0).
#' @return data.frame with 8 rows (single-analyte) or 8 * n_analytes rows
#'   (panel mode).
run_one_rep <- function(seed, n_subjects, effect_size, censoring_target,
                        residual_sd = 1.0, n_analytes = 1L,
                        n_signal = NULL) {

  # Determine signal count
  if (is.null(n_signal)) {
    n_signal_use <- if (effect_size != 0) min(1L, n_analytes) else 0L
  } else {
    n_signal_use <- if (effect_size != 0) n_signal else 0L
  }

  # --- Simulate data (clean assay: no plate/lot/replicate effects) ---
  sim <- simulate_immunoassay(
    n_subjects      = n_subjects,
    n_timepoints    = 1L,
    n_analytes      = n_analytes,
    design          = "cross_sectional",
    group_levels    = c("control", "treatment"),
    group_effects   = effect_size,
    signal_analytes = n_signal_use,
    effect_direction = "up",
    re_intercept_sd = 0,
    residual_sd     = residual_sd,
    lod_quantile    = censoring_target,
    seed            = seed
  )

  d  <- as.data.frame(sim$data[!is.na(sim$data$value), ])
  tr <- as.data.frame(sim$truth[!is.na(sim$data$value), ])

  # --- Single-analyte path ---
  if (n_analytes == 1L) {
    realized_cens <- mean(d$cens_lod, na.rm = TRUE)
    results <- .fit_all_substitution_methods(d, tr)
    results$true_effect <- effect_size
    results$realized_censoring <- realized_cens
    return(results)
  }

  # --- Multi-analyte panel path ---
  # Shared LOD for heterogeneous per-analyte censoring:
  # lod_quantile calibrates each analyte's LOD independently, giving every
  # analyte ~the same censoring rate.  Instead, take the median per-analyte
  # LOD and apply it to all analytes.  Because analytes have different grand
  # means, the single threshold censors low-abundance analytes heavily and
  # high-abundance analytes lightly — the realistic heterogeneous-LOD
  # scenario that drives the B5 multi-analyte story.
  shared_lod <- median(unique(d$lod))
  true_raw   <- exp(tr$true_latent)
  new_cens   <- true_raw < shared_lod
  d$lod      <- shared_lod
  d$cens_lod <- new_cens
  d$value_raw <- ifelse(new_cens, shared_lod, true_raw)
  d$value     <- ifelse(new_cens, log(shared_lod), tr$true_latent)

  # Column is 'cytokine' in simulate_immunoassay output; rename to 'analyte'
  # for compatibility with aggregate_panel_cell()
  analyte_names <- unique(d$cytokine)
  all_results <- vector("list", length(analyte_names))

  for (k in seq_along(analyte_names)) {
    a_name <- analyte_names[k]
    d_k  <- d[d$cytokine == a_name, ]
    tr_k <- tr[d$cytokine == a_name, ]

    is_signal_k     <- tr_k$signal_analyte[1]
    true_effect_k   <- if (is_signal_k) effect_size else 0
    realized_cens_k <- mean(d_k$cens_lod, na.rm = TRUE)

    method_results <- .fit_all_substitution_methods(d_k, tr_k)
    method_results$analyte            <- a_name
    method_results$is_signal          <- is_signal_k
    method_results$true_effect        <- true_effect_k
    method_results$realized_censoring <- realized_cens_k

    all_results[[k]] <- method_results
  }

  do.call(rbind, all_results)
}


# ---- Run Benchmark ----------------------------------------------------------

run_benchmark_substitution <- function(dgp_grid, n_reps, base_seed,
                                       cache_dir, n_cores, verbose) {
  all_summaries <- list()

  for (i in seq_len(nrow(dgp_grid))) {
    row <- dgp_grid[i, ]
    cell_label <- sprintf("cell_%d [cens=%.0f%%/n=%d/eff=%.1f/sd=%.1f]",
                          row$cell_id, row$censoring_target * 100,
                          row$n_subjects, row$effect_size, row$residual_sd)

    rep_fn <- local({
      ns <- row$n_subjects
      es <- row$effect_size
      ct <- row$censoring_target
      sd <- row$residual_sd
      function(seed) {
        run_one_rep(
          seed             = seed,
          n_subjects       = ns,
          effect_size      = es,
          censoring_target = ct,
          residual_sd      = sd
        )
      }
    })

    # Descriptive cache filenames to avoid collisions when grid changes
    cache_file <- file.path(cache_dir,
                            sprintf("%s_cens%s_n%d_eff%s_sd%s.rds",
                                    BENCHMARK,
                                    format(row$censoring_target, nsmall = 2),
                                    row$n_subjects,
                                    format(row$effect_size, nsmall = 1),
                                    format(row$residual_sd, nsmall = 1)))

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

      # Add timing and CI width summaries
      agg$mean_time_s <- mean(m_reps$time_s, na.rm = TRUE)
      if (all(c("ci_lo", "ci_hi") %in% names(m_reps))) {
        ciw <- calc_ci_width(m_reps$ci_lo, m_reps$ci_hi)
        agg$ci_width      <- ciw$estimate
        agg$ci_width_mcse <- ciw$mcse
      }

      # Hessian diagnostic rate (Gamma/Gaussian via glmmTMB)
      if ("pd_hess" %in% names(m_reps)) {
        pd_vals <- m_reps$pd_hess[!is.na(m_reps$pd_hess)]
        agg$pd_hess_rate <- if (length(pd_vals) > 0) mean(pd_vals) else NA_real_
      }

      # Cell identifiers
      agg$method            <- m
      agg$censoring_target  <- row$censoring_target
      agg$n_subjects        <- row$n_subjects
      agg$effect_size       <- row$effect_size
      agg$residual_sd       <- row$residual_sd
      agg$cell_id           <- row$cell_id
      all_summaries[[length(all_summaries) + 1]] <- agg
    }
  }

  bind_rows(all_summaries)
}


# ---- Execute ----------------------------------------------------------------

if (is.null(PHASE_B)) {

if (VERBOSE) message("\nStarting single-analyte benchmark execution...")

summary_df <- run_benchmark_substitution(
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

# Save final summary (merge with existing when running extra sizes)
summary_file <- file.path(CACHE_DIR, paste0(BENCHMARK, "_summary.rds"))
if (!is.null(EXTRA_SIZES)) {
  merge_with_existing_summary(summary_annotated, summary_file, verbose = VERBOSE)
} else {
  save_main_summary_preserving_supplements(summary_annotated, summary_file,
                                           verbose = VERBOSE)
}


# ---- Figures ----------------------------------------------------------------

# Readable labellers (shared by main + panel figures)
cens_labeller <- function(x) paste0(as.numeric(x) * 100, "% cens.")
n_labeller    <- function(x) paste0("n = ", x)
eff_labeller  <- function(x) paste0("Effect = ", x)
sd_labeller   <- function(x) paste0("SD = ", x)

# Method ordering for consistent display (shared by main + panel figures)
method_order <- c("Oracle", "Tobit",
                  "LOD/2 + Gaussian", "LOD/sqrt2 + Gaussian",
                  "halfmin + Gaussian", "LOD + Gaussian",
                  "zero + Gaussian", "LOD/2 + Gamma")
summary_annotated$method <- factor(summary_annotated$method,
                                    levels = method_order)

fig_dir <- file.path(CACHE_DIR, "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)


# Figures are generated by the standalone plot_B5.R script.
# Run: Rscript inst/simulations/plot_B5.R [--cache_dir PATH]


# ===========================================================================
# PART 2: Multi-Analyte Panel Extension
# ===========================================================================
#
# 10-analyte panel with heterogeneous per-analyte LODs from the built-in
# cytokine library. 3 of 10 analytes are signal (receive group effect).
# Each analyte gets its own LOD, so realized censoring varies across
# analytes within the same dataset. The substitution method choice
# (LOD/2 vs Tobit etc.) affects each analyte differently — the key
# scientific insight for this benchmark.
#
# Panel sub-grid (18 cells):
#   n_subjects = {40, 100, 200}
#   effect_size = {0, 0.5, 1.0}
#   censoring_target = {0.15, 0.50}
#   residual_sd = 1.0 (fixed)
# ===========================================================================

if (VERBOSE) message("\n=== Panel Extension: 10-analyte heterogeneous LOD ===")

N_PANEL_ANALYTES  <- 10L
N_PANEL_SIGNAL    <- 3L
PANEL_RESIDUAL_SD <- 1.0

panel_grid <- make_dgp_grid(
  n_subjects       = if (!is.null(EXTRA_SIZES)) EXTRA_SIZES else BACKBONE_PARAMS$n_subjects,
  effect_size      = BACKBONE_PARAMS$effect_size,
  censoring_target = c(0.15, 0.50)
)
if (!is.null(EXTRA_SIZES)) panel_grid$cell_id <- panel_grid$cell_id + EXTRA_OFFSET

if (VERBOSE) {
  message("Panel grid: ", nrow(panel_grid), " cells (3 n x 3 eff x 2 cens = 18)")
  message("Analytes per sim: ", N_PANEL_ANALYTES,
          " (", N_PANEL_SIGNAL, " signal)")
  message("Reps: ", N_REPS, " | Total runs: ", nrow(panel_grid) * N_REPS)
}


# ---- Panel Execution --------------------------------------------------------

panel_summaries         <- vector("list", nrow(panel_grid))
panel_analyte_summaries <- vector("list", nrow(panel_grid))

for (i in seq_len(nrow(panel_grid))) {
  row <- panel_grid[i, ]
  cell_label <- sprintf("panel_%d [cens=%.0f%%/n=%d/eff=%.1f]",
                        row$cell_id, row$censoring_target * 100,
                        row$n_subjects, row$effect_size)

  # Use local() for safe closure binding in the loop
  rep_fn <- local({
    ns <- row$n_subjects
    es <- row$effect_size
    ct <- row$censoring_target
    function(seed) {
      run_one_rep(seed = seed, n_subjects = ns, effect_size = es,
                  censoring_target = ct, residual_sd = PANEL_RESIDUAL_SD,
                  n_analytes = N_PANEL_ANALYTES,
                  n_signal   = N_PANEL_SIGNAL)
    }
  })

  cache_file <- file.path(CACHE_DIR,
                          sprintf("%s_panel_cell%d.rds",
                                  BENCHMARK, row$cell_id))

  raw_results <- run_benchmark_cell(
    rep_fn     = rep_fn,
    n_reps     = N_REPS,
    base_seed  = BASE_SEED + 500000L + (i - 1) * 10000L,
    cache_file = cache_file,
    n_cores    = N_CORES,
    cell_label = cell_label,
    verbose    = VERBOSE
  )

  # --- Panel-level aggregation (FDR, FWER, sensitivity) ---
  panel_agg <- aggregate_panel_cell(raw_results)
  panel_agg$censoring_target <- row$censoring_target
  panel_agg$n_subjects       <- row$n_subjects
  panel_agg$effect_size      <- row$effect_size
  panel_agg$cell_id          <- row$cell_id
  panel_summaries[[i]] <- panel_agg

  # --- Per-analyte aggregation (for bias vs. realized censoring scatter) ---
  ok <- raw_results[!raw_results$error & !is.na(raw_results$method), ]
  analyte_agg <- ok %>%
    group_by(method, analyte) %>%
    summarise(
      n_reps             = n(),
      mean_bias          = mean(estimate - true_effect, na.rm = TRUE),
      mean_realized_cens = mean(realized_censoring, na.rm = TRUE),
      is_signal          = first(is_signal),
      .groups = "drop"
    )
  analyte_agg$censoring_target <- row$censoring_target
  analyte_agg$n_subjects       <- row$n_subjects
  analyte_agg$effect_size      <- row$effect_size
  panel_analyte_summaries[[i]] <- analyte_agg
}

panel_summary         <- bind_rows(panel_summaries)
panel_analyte_summary <- bind_rows(panel_analyte_summaries)

# Save panel summaries (merge with existing when running extra sizes)
panel_summary_file <- file.path(CACHE_DIR,
                                paste0(BENCHMARK, "_panel_summary.rds"))
panel_analyte_file <- file.path(CACHE_DIR,
                                paste0(BENCHMARK, "_panel_analyte_summary.rds"))
if (!is.null(EXTRA_SIZES)) {
  merge_with_existing_summary(panel_summary, panel_summary_file, verbose = VERBOSE)
  merge_with_existing_summary(panel_analyte_summary, panel_analyte_file, verbose = VERBOSE)
} else {
  saveRDS(panel_summary, panel_summary_file)
  saveRDS(panel_analyte_summary, panel_analyte_file)
}

if (VERBOSE) {
  message("Panel summary saved -> ", panel_summary_file)
  message("Panel analyte summary saved -> ", panel_analyte_file)
}



# Panel figures are generated by the standalone plot_B5.R script.


# ---- Panel Summary Table ----------------------------------------------------

if (VERBOSE) {
  message("\n=== Panel Summary (effect = 0, null scenario) ===\n")
  panel_null <- panel_summary %>%
    filter(effect_size == 0) %>%
    mutate(fdr_fmt = fmt_mcse_vec(panel_fdr, panel_fdr_mcse))
  print(
    panel_null %>%
      select(censoring_target, n_subjects, method,
             fdr_fmt, panel_fwer) %>%
      as.data.frame(),
    right = FALSE
  )

  message("\n=== Panel Summary (effect != 0, sensitivity) ===\n")
  panel_alt <- panel_summary %>%
    filter(effect_size != 0) %>%
    mutate(
      fdr_fmt  = fmt_mcse_vec(panel_fdr, panel_fdr_mcse),
      sens_fmt = fmt_mcse_vec(panel_sensitivity, panel_sensitivity_mcse)
    )
  print(
    panel_alt %>%
      select(censoring_target, n_subjects, effect_size, method,
             fdr_fmt, sens_fmt) %>%
      as.data.frame(),
    right = FALSE
  )
}


}  # end if (is.null(PHASE_B))


# ===========================================================================
# PART 3: Phase B7 — Mis-Specified DGP Sub-Grid
# ===========================================================================
#
# Tests robustness when DGP assumptions are violated:
#   (a) t_errors:  residuals follow t(df=5) instead of Gaussian
#   (b) lod_error: LOD values have ±10% measurement error
#   (c) both:      both violations simultaneously
#
# Phase B7 sub-grid (81 cells):
#   censoring_target = {0.15, 0.30, 0.50}
#   n_subjects       = {40, 100, 200}
#   effect_size      = {0, 0.5, 1.0}
#   residual_sd      = 1.0
#   dgp_violation    = {"t_errors", "lod_error", "both"}
# ===========================================================================

if (!is.null(PHASE_B)) {

if (VERBOSE) message("\n=== Phase B7: Mis-Specified DGP Sub-Grid ===")


# ---- Mis-Specification Helper ------------------------------------------------

#' Apply DGP mis-specification to simulated data
#'
#' Transforms normally-simulated data to introduce DGP violations:
#' - t_errors:  replace Gaussian residuals with t(df=5) via PIT
#' - lod_error: perturb each observation's LOD by U(-10%, +10%)
#' - both:      apply both transformations
#'
#' @param d     data.frame from simulate_immunoassay (value, value_raw, lod,
#'              cens_lod columns).
#' @param tr    truth data.frame (true_latent column aligned with d).
#' @param violation One of "t_errors", "lod_error", "both".
#' @param residual_sd Residual SD used in the original simulation.
#' @return list with modified \code{d} and \code{tr} (tr unchanged for lod_error).
apply_dgp_misspec <- function(d, tr, violation, residual_sd = 1.0) {

  n <- nrow(d)

  # --- t-distributed errors via probability integral transform ---
  if (violation %in% c("t_errors", "both")) {
    # Extract Gaussian residuals from the simulation
    gauss_resid <- d$value - (tr$true_latent - d$value + d$value)
    # Actually: true_latent already includes the Gaussian residual.
    # The structural part (without residual) is: true_latent - epsilon
    # But we don't have epsilon separately. Instead, use PIT:
    # 1. Convert Gaussian residuals to uniform via pnorm
    # 2. Convert uniform to t(df=5)
    # residuals = value - (true_latent - residual). Since true_latent
    # already includes the residual, we can't separate it directly.
    # Instead: generate fresh t-errors and replace the values.
    t_df <- 5
    # Scale factor so that Var(t_error) = residual_sd^2
    # Var(t(df)) = df/(df-2), so sd = sqrt(df/(df-2))
    # We want errors with sd = residual_sd, so scale by residual_sd / sqrt(df/(df-2))
    t_scale <- residual_sd * sqrt((t_df - 2) / t_df)
    t_errors <- stats::rt(n, df = t_df) * t_scale

    # Reconstruct: structural part = true_latent - original_residual
    # Since true_latent = structural + gaussian_error, and we generated the
    # sim with known residual_sd, the structural part is:
    #   structural = true_group_effect + intercept + ...
    # We can get structural = true_latent - (value - true_latent)... no.
    # Actually true_latent IS the full latent value (structural + residual).
    # The simulation stores: true_latent = mu + effects + epsilon
    # So to get the structural part without residual we'd need the raw epsilon.
    #
    # Simpler approach: compute residuals from the group mean structure,
    # then replace them. The Oracle model uses true_latent directly, so we
    # need to modify both d$value and tr$true_latent consistently.
    #
    # Approach: fit the known group structure to get structural means,
    # then add t-errors.
    group_means <- tapply(tr$true_latent, d$group, mean)
    structural <- group_means[as.character(d$group)]
    new_latent <- as.numeric(structural) + t_errors

    # Update truth and observed
    tr$true_latent <- new_latent
    d$value     <- new_latent
    d$value_raw <- exp(new_latent)

    # Re-apply LOD censoring
    below_lod <- d$value_raw < exp(log(d$lod[1]))  # lod is on raw scale
    d$cens_lod  <- below_lod
    d$value[below_lod]     <- log(d$lod[below_lod])
    d$value_raw[below_lod] <- d$lod[below_lod]
  }

  # --- LOD measurement error ---
  if (violation %in% c("lod_error", "both")) {
    lod_perturbation <- 1 + stats::runif(n, -0.10, 0.10)
    effective_lod <- d$lod * lod_perturbation

    # Re-apply censoring with perturbed LOD
    true_raw <- exp(tr$true_latent)
    below_lod <- true_raw < effective_lod
    d$lod       <- effective_lod
    d$cens_lod  <- below_lod
    d$value[below_lod]     <- log(effective_lod[below_lod])
    d$value_raw[below_lod] <- effective_lod[below_lod]
    # Uncensored values remain at true_latent
    d$value[!below_lod]     <- tr$true_latent[!below_lod]
    d$value_raw[!below_lod] <- true_raw[!below_lod]
  }

  list(d = d, tr = tr)
}


# ---- Mis-Specified Replication Function --------------------------------------

#' Run one replication with DGP mis-specification
#'
#' Generates data normally via simulate_immunoassay, then applies the specified
#' DGP violation before fitting all 8 methods.
#'
#' @inheritParams run_one_rep
#' @param dgp_violation One of "t_errors", "lod_error", "both".
#' @return data.frame with 8 rows (one per method), plus dgp_violation column.
run_one_rep_misspec <- function(seed, n_subjects, effect_size, censoring_target,
                                residual_sd = 1.0, dgp_violation = "both") {

  # Simulate data normally (clean assay, Gaussian errors)
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
    lod_quantile    = censoring_target,
    seed            = seed
  )

  d  <- as.data.frame(sim$data[!is.na(sim$data$value), ])
  tr <- as.data.frame(sim$truth[!is.na(sim$data$value), ])

  # Apply DGP mis-specification
  set.seed(seed + 999999L)
  misspec <- apply_dgp_misspec(d, tr, dgp_violation, residual_sd)
  d  <- misspec$d
  tr <- misspec$tr

  realized_cens <- mean(d$cens_lod, na.rm = TRUE)
  results <- .fit_all_substitution_methods(d, tr)
  results$true_effect        <- effect_size
  results$realized_censoring <- realized_cens
  results$dgp_violation      <- dgp_violation
  results
}


# ---- Phase B7 Grid -----------------------------------------------------------

# Build the base grid (without dgp_violation, which we loop over separately)
phaseB7_base <- make_dgp_grid(
  censoring_target = c(0.15, 0.30, 0.50),
  n_subjects       = c(40L, 100L, 200L),
  effect_size      = c(0, 0.5, 1.0),
  residual_sd      = 1.0
)

# Expand with dgp_violation
dgp_violations <- c("t_errors", "lod_error", "both")
phaseB7_grid <- do.call(rbind, lapply(seq_along(dgp_violations), function(v) {
  g <- phaseB7_base
  g$dgp_violation <- dgp_violations[v]
  g$cell_id <- g$cell_id + (v - 1) * nrow(phaseB7_base) + PHASE_B_OFFSET

  g
}))

if (VERBOSE) {
  message("Phase B7 grid: ", nrow(phaseB7_grid), " cells ",
          "(3 cens x 3 n x 3 eff x 1 sd x 3 dgp = 81)")
  message("Reps: ", N_REPS, " | Total runs: ", nrow(phaseB7_grid) * N_REPS)
}


# ---- Phase B7 Execution -----------------------------------------------------

phaseB7_summaries <- list()

for (i in seq_len(nrow(phaseB7_grid))) {
  row <- phaseB7_grid[i, ]
  cell_label <- sprintf("phaseB7_%d [cens=%.0f%%/n=%d/eff=%.1f/dgp=%s]",
                        row$cell_id, row$censoring_target * 100,
                        row$n_subjects, row$effect_size, row$dgp_violation)

  rep_fn <- local({
    ns  <- row$n_subjects
    es  <- row$effect_size
    ct  <- row$censoring_target
    sd  <- row$residual_sd
    dgp <- row$dgp_violation
    function(seed) {
      run_one_rep_misspec(
        seed             = seed,
        n_subjects       = ns,
        effect_size      = es,
        censoring_target = ct,
        residual_sd      = sd,
        dgp_violation    = dgp
      )
    }
  })

  cache_file <- file.path(CACHE_DIR,
                          sprintf("%s_phaseB_cell%d.rds",
                                  BENCHMARK, row$cell_id))

  raw_results <- run_benchmark_cell(
    rep_fn     = rep_fn,
    n_reps     = N_REPS,
    base_seed  = BASE_SEED + 700000L + (i - 1) * 10000L,
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

    # Add timing and CI width summaries
    agg$mean_time_s <- mean(m_reps$time_s, na.rm = TRUE)
    if (all(c("ci_lo", "ci_hi") %in% names(m_reps))) {
      ciw <- calc_ci_width(m_reps$ci_lo, m_reps$ci_hi)
      agg$ci_width      <- ciw$estimate
      agg$ci_width_mcse <- ciw$mcse
    }

    # Hessian diagnostic rate
    if ("pd_hess" %in% names(m_reps)) {
      pd_vals <- m_reps$pd_hess[!is.na(m_reps$pd_hess)]
      agg$pd_hess_rate <- if (length(pd_vals) > 0) mean(pd_vals) else NA_real_
    }

    # Cell identifiers
    agg$method            <- m
    agg$censoring_target  <- row$censoring_target
    agg$n_subjects        <- row$n_subjects
    agg$effect_size       <- row$effect_size
    agg$residual_sd       <- row$residual_sd
    agg$cell_id           <- row$cell_id
    agg$dgp_violation     <- row$dgp_violation
    phaseB7_summaries[[length(phaseB7_summaries) + 1]] <- agg
  }
}

phaseB7_summary <- bind_rows(phaseB7_summaries)

# Add acceptance annotations
phaseB7_null <- phaseB7_summary %>%
  filter(effect_size == 0) %>%
  mutate(
    pass_bias     = NA,
    pass_coverage = NA,
    pass_type1    = passes_threshold(rejection_rate, 0.025, 0.075)
  )

phaseB7_alt <- phaseB7_summary %>%
  filter(effect_size != 0) %>%
  mutate(
    pass_bias     = passes_threshold(abs(bias_ratio), upper = 0.1),
    pass_coverage = passes_threshold(coverage, 0.93, 0.97),
    pass_type1    = NA
  )

phaseB7_annotated <- bind_rows(phaseB7_null, phaseB7_alt)

# Save Phase B7 summary — merge into the MAIN summary file so plot scripts
# see all data together (Phase B cells have offset IDs, no collision)
phaseB7_summary_file <- file.path(CACHE_DIR, paste0(BENCHMARK, "_summary.rds"))
merge_with_existing_summary(phaseB7_annotated, phaseB7_summary_file,
                            verbose = VERBOSE)

if (VERBOSE) {
  message("Phase B7 summary merged -> ", phaseB7_summary_file)
  message("Phase B7: ", nrow(phaseB7_grid), " cells, ",
          nrow(phaseB7_annotated), " summary rows")
}

}  # end if (!is.null(PHASE_B))


# ---- Done --------------------------------------------------------------------

if (VERBOSE && is.null(PHASE_B)) {
  message("\n=== Benchmark 5 (v3) complete ===")
  message("Single-analyte: ", nrow(dgp_grid), " cells, ",
          nrow(summary_annotated), " summary rows")
  message("Panel: ", nrow(panel_grid), " cells, ",
          nrow(panel_summary), " summary rows")
  message("Figures -> ", fig_dir)
  message("Summary -> ", summary_file)
  message("Panel   -> ", panel_summary_file)
}
