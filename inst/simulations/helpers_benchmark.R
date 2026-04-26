# helpers_benchmark.R — Shared infrastructure for immunoPlex Monte Carlo benchmarks
#
# Provides:
#   - Performance metric functions (bias, RMSE, coverage, power, type I error)
#   - Monte Carlo standard error (MCSE) calculation
#   - Parallel execution wrapper with progress & caching
#   - Result aggregation and formatting
#   - Publication-quality ggplot2 theme and colour palette

# ---- Dependencies -----------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

# ---- Performance Metrics ----------------------------------------------------

#' Bias of an estimator relative to the true value
#' @param estimates Numeric vector of point estimates across replications.
#' @param truth     Scalar true value.
#' @return Named list: estimate, mcse.
calc_bias <- function(estimates, truth) {
  estimates <- estimates[is.finite(estimates)]
  n <- length(estimates)
  if (n == 0) return(list(estimate = NA_real_, mcse = NA_real_))
  bias <- mean(estimates) - truth
  mcse <- sd(estimates) / sqrt(n)
  list(estimate = bias, mcse = mcse)
}

#' Root mean squared error
#' @inheritParams calc_bias
#' @return Named list: estimate, mcse.
calc_rmse <- function(estimates, truth) {
  estimates <- estimates[is.finite(estimates)]
  n <- length(estimates)
  if (n == 0) return(list(estimate = NA_real_, mcse = NA_real_))
  sq_err <- (estimates - truth)^2
  rmse <- sqrt(mean(sq_err))
  # Delta-method MCSE for RMSE: sd(squared_errors) / (2 * RMSE * sqrt(n))
  mcse <- if (rmse > 0) sd(sq_err) / (2 * rmse * sqrt(n)) else NA_real_
  list(estimate = rmse, mcse = mcse)
}

#' Empirical coverage of confidence intervals
#' @param ci_lo   Numeric vector of lower CI bounds.
#' @param ci_hi   Numeric vector of upper CI bounds.
#' @param truth   Scalar true value.
#' @return Named list: estimate, mcse.
calc_coverage <- function(ci_lo, ci_hi, truth) {
  valid <- is.finite(ci_lo) & is.finite(ci_hi)
  ci_lo <- ci_lo[valid]; ci_hi <- ci_hi[valid]
  n <- length(ci_lo)
  if (n == 0) return(list(estimate = NA_real_, mcse = NA_real_))
  covers <- (ci_lo <= truth) & (ci_hi >= truth)
  p <- mean(covers)
  mcse <- sqrt(p * (1 - p) / n)
  list(estimate = p, mcse = mcse)
}

#' Empirical rejection rate (power or type I error)
#' @param p_values Numeric vector of p-values across replications.
#' @param alpha    Significance threshold (default 0.05).
#' @return Named list: estimate, mcse.
calc_rejection_rate <- function(p_values, alpha = 0.05) {
  p_values <- p_values[is.finite(p_values)]
  n <- length(p_values)
  if (n == 0) return(list(estimate = NA_real_, mcse = NA_real_))
  p <- mean(p_values < alpha)
  mcse <- sqrt(p * (1 - p) / n)
  list(estimate = p, mcse = mcse)
}

#' Convergence rate
#' @param converged Logical vector indicating convergence per replication.
#' @return Named list: estimate, mcse.
calc_convergence_rate <- function(converged) {
  converged <- converged[!is.na(converged)]
  n <- length(converged)
  if (n == 0) return(list(estimate = NA_real_, mcse = NA_real_))
  p <- mean(converged)
  mcse <- sqrt(p * (1 - p) / n)
  list(estimate = p, mcse = mcse)
}

#' Relative bias: bias / true value (meaningful when truth != 0)
#' @inheritParams calc_bias
#' @return Named list: estimate, mcse.
calc_relative_bias <- function(estimates, truth) {
  if (is.na(truth) || truth == 0) return(list(estimate = NA_real_, mcse = NA_real_))
  b <- calc_bias(estimates, truth)
  list(estimate = b$estimate / abs(truth), mcse = b$mcse / abs(truth))
}

#' Bias ratio: bias / empirical SE of the estimator
#' @inheritParams calc_bias
#' @return Named list: estimate, mcse.
calc_bias_ratio <- function(estimates, truth) {
  estimates <- estimates[is.finite(estimates)]
  n <- length(estimates)
  if (n < 2) return(list(estimate = NA_real_, mcse = NA_real_))
  bias <- mean(estimates) - truth
  emp_se <- sd(estimates)
  if (emp_se == 0) return(list(estimate = NA_real_, mcse = NA_real_))
  ratio <- bias / emp_se
  # Approximate MCSE via delta method
  mcse <- sqrt((1 + ratio^2 / 2) / n)
  list(estimate = ratio, mcse = mcse)
}

#' Mean CI width
#' @param ci_lo Numeric vector of lower CI bounds.
#' @param ci_hi Numeric vector of upper CI bounds.
#' @return Named list: estimate, mcse.
calc_ci_width <- function(ci_lo, ci_hi) {
  widths <- (ci_hi - ci_lo)[is.finite(ci_hi) & is.finite(ci_lo)]
  n <- length(widths)
  if (n == 0) return(list(estimate = NA_real_, mcse = NA_real_))
  list(estimate = mean(widths), mcse = sd(widths) / sqrt(n))
}


# ---- MCSE Utilities ---------------------------------------------------------

#' Required replications to achieve a target MCSE for a proportion
#' @param p     Anticipated proportion (e.g. 0.05 for type I error).
#' @param target_mcse Target MCSE (e.g. 0.005).
#' @return Integer number of replications needed.
reps_for_proportion_mcse <- function(p, target_mcse) {
  ceiling(p * (1 - p) / target_mcse^2)
}

#' Required replications to achieve a target MCSE for bias ratio
#' @param target_mcse Target MCSE for the bias ratio (e.g. 0.05).
#' @param bias_ratio  Anticipated bias ratio (default 0, worst case).
#' @return Integer number of replications needed.
reps_for_bias_mcse <- function(target_mcse, bias_ratio = 0) {
  ceiling((1 + bias_ratio^2 / 2) / target_mcse^2)
}

#' Format a metric ± MCSE for display
#' @param estimate Point estimate.
#' @param mcse     Monte Carlo standard error.
#' @param digits   Rounding digits (default 3).
#' @return Character string like "0.053 (0.005)".
fmt_mcse <- function(estimate, mcse, digits = 3) {
  if (is.na(estimate)) return("NA")
  paste0(round(estimate, digits), " (", round(mcse, digits), ")")
}

#' Vectorised version of fmt_mcse for use in dplyr mutate
#' @inheritParams fmt_mcse
#' @return Character vector.
fmt_mcse_vec <- function(estimate, mcse, digits = 3) {
  mapply(fmt_mcse, estimate, mcse, digits, USE.NAMES = FALSE)
}


# ---- Parallel Execution Wrapper ---------------------------------------------

#' Run a single replication of a benchmark
#'
#' Internal helper — executes \code{rep_fn} inside tryCatch so that a single
#' failing replication does not abort the whole benchmark.
#'
#' @param rep_id   Integer replication number (used as seed offset).
#' @param rep_fn   Function(seed) -> one-row data.frame of per-replication results.
#' @param base_seed Integer base seed; actual seed = base_seed + rep_id.
#' @return One-row data.frame on success, or one-row data.frame with
#'   \code{error = TRUE} and \code{error_msg} on failure. Warnings emitted
#'   during \code{rep_fn} are captured in \code{warning_msg}. If every row
#'   returned by \code{rep_fn} has \code{converged = FALSE} AND warnings were
#'   captured, the rep is treated as a *silent* failure: \code{error} is
#'   flipped to TRUE and the concatenated warning text is written to
#'   \code{error_msg}. This catches the pattern where a benchmark's inner
#'   tryCatch swallows an error, issues a warning, and returns fail_result
#'   with converged=FALSE — which previously looked indistinguishable from a
#'   legitimate non-convergence in the cache.
safe_replicate <- function(rep_id, rep_fn, base_seed) {
  seed <- base_seed + rep_id
  warnings_caught <- character(0)
  tryCatch(
    {
      result <- withCallingHandlers(
        rep_fn(seed),
        warning = function(w) {
          warnings_caught <<- c(warnings_caught, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      result$rep_id <- rep_id
      result$seed   <- seed
      warning_text <- if (length(warnings_caught) > 0) {
        paste(unique(warnings_caught), collapse = " | ")
      } else {
        NA_character_
      }
      result$warning_msg <- warning_text

      # Silent-failure detection: fail_result rows have converged=FALSE and
      # no real metrics. If *every* returned row is non-converged AND rep_fn
      # warned, treat as failure rather than legitimate non-convergence.
      silent_fail <- length(warnings_caught) > 0 &&
                     "converged" %in% names(result) &&
                     !any(result$converged %in% TRUE)
      result$error     <- silent_fail
      result$error_msg <- if (silent_fail) warning_text else NA_character_
      result
    },
    error = function(e) {
      data.frame(
        rep_id      = rep_id,
        seed        = seed,
        error       = TRUE,
        error_msg   = conditionMessage(e),
        warning_msg = if (length(warnings_caught) > 0) {
          paste(unique(warnings_caught), collapse = " | ")
        } else {
          NA_character_
        },
        stringsAsFactors = FALSE
      )
    }
  )
}

#' Execute a benchmark cell (one DGP configuration) across replications
#'
#' Runs \code{rep_fn} for \code{n_reps} replications, optionally in parallel.
#' Supports caching: if a cache file exists with >= \code{n_reps} rows, it is
#' returned immediately. If the cache has fewer rows, only the missing
#' replications are run and appended (top-up).
#'
#' @param rep_fn     Function(seed) -> one-row data.frame.
#' @param n_reps     Target number of replications.
#' @param base_seed  Integer base seed for reproducibility.
#' @param cache_file Optional path to an .rds cache file.
#' @param n_cores    Number of parallel cores (1 = sequential).
#' @param cell_label Character label for progress messages.
#' @param verbose    Logical; print progress messages.
#' @return data.frame of concatenated per-replication results.
run_benchmark_cell <- function(rep_fn,
                               n_reps,
                               base_seed,
                               cache_file = NULL,
                               n_cores    = 1L,
                               cell_label = "",
                               verbose    = TRUE) {

  existing <- NULL
  start_rep <- 1L


# --- Cache check ---
  if (!is.null(cache_file) && file.exists(cache_file)) {
    existing <- readRDS(cache_file)
    # For multi-row rep_fn (e.g. one row per method), count unique rep_ids
    n_existing <- if ("rep_id" %in% names(existing)) {
      length(unique(existing$rep_id))
    } else {
      nrow(existing)
    }
    if (n_existing >= n_reps) {
      if (verbose) message("[", cell_label, "] Cache hit: ", n_existing,
                           " reps (>= target ", n_reps, "). Skipping.")
      if ("rep_id" %in% names(existing)) {
        target_ids <- sort(unique(existing$rep_id))[seq_len(n_reps)]
        return(existing[existing$rep_id %in% target_ids, , drop = FALSE])
      }
      return(existing[seq_len(n_reps), , drop = FALSE])
    }
    start_rep <- n_existing + 1L
    if (verbose) message("[", cell_label, "] Cache has ", n_existing,
                         " reps. Topping up ", n_reps - n_existing, " more.")
  }

  rep_ids <- seq.int(start_rep, n_reps)

  if (verbose) {
    message("[", cell_label, "] Running ", length(rep_ids), " replications",
            if (n_cores > 1) paste0(" on ", n_cores, " cores") else " sequentially",
            " ...")
  }

  t0 <- proc.time()

  if (n_cores > 1L && requireNamespace("parallel", quietly = TRUE)) {
    results <- parallel::mclapply(
      rep_ids, safe_replicate,
      rep_fn = rep_fn, base_seed = base_seed,
      mc.cores = n_cores
    )
  } else {
    results <- lapply(rep_ids, safe_replicate,
                      rep_fn = rep_fn, base_seed = base_seed)
  }

  # Bind rows — handle heterogeneous columns from errors
  new_results <- do.call(rbind, lapply(results, function(x) {
    as.data.frame(x, stringsAsFactors = FALSE)
  }))

  elapsed <- (proc.time() - t0)["elapsed"]
  n_errors <- sum(new_results$error, na.rm = TRUE)
  if (verbose) {
    message("[", cell_label, "] Done in ", round(elapsed, 1), "s. ",
            nrow(new_results), " reps completed",
            if (n_errors > 0) paste0(" (", n_errors, " errors)") else "",
            ".")
  }

  # Combine with cached results
  if (!is.null(existing)) {
    # Align columns
    all_cols <- union(names(existing), names(new_results))
    for (col in setdiff(all_cols, names(existing)))
      existing[[col]] <- NA
    for (col in setdiff(all_cols, names(new_results)))
      new_results[[col]] <- NA
    combined <- rbind(existing[, all_cols, drop = FALSE],
                      new_results[, all_cols, drop = FALSE])
  } else {
    combined <- new_results
  }

  # --- Cache write ---
  if (!is.null(cache_file)) {
    cache_dir <- dirname(cache_file)
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
    saveRDS(combined, cache_file)
    if (verbose) message("[", cell_label, "] Cached ", nrow(combined),
                         " reps -> ", cache_file)
  }

  combined
}


# ---- Result Aggregation -----------------------------------------------------

#' Aggregate per-replication results into summary statistics
#'
#' Given a data.frame of per-replication results (from \code{run_benchmark_cell}),
#' computes standard performance metrics for a single method × DGP cell.
#'
#' @param reps         data.frame of per-replication results. Must contain
#'   \code{estimate} (point estimate) and optionally \code{ci_lo}, \code{ci_hi},
#'   \code{p_value}, \code{converged}.
#' @param truth        Scalar true parameter value.
#' @param alpha        Significance level for rejection rate (default 0.05).
#' @return One-row tibble with columns: n_reps, n_errors, convergence_rate (+mcse),
#'   bias (+mcse), rmse (+mcse), coverage (+mcse), rejection_rate (+mcse),
#'   bias_ratio (+mcse), mean_estimate, emp_se.
aggregate_cell <- function(reps, truth, alpha = 0.05) {

  # Exclude error rows
  ok <- reps[!reps$error, , drop = FALSE]

  n_total  <- nrow(reps)
  n_errors <- sum(reps$error, na.rm = TRUE)

  # Convergence
  conv <- if ("converged" %in% names(ok)) {
    calc_convergence_rate(ok$converged)
  } else {
    list(estimate = NA_real_, mcse = NA_real_)
  }

  # Bias, RMSE, bias ratio
  est_vec    <- if ("estimate" %in% names(ok)) ok$estimate else numeric(0)
  .bias      <- calc_bias(est_vec, truth)
  .rmse      <- calc_rmse(est_vec, truth)
  .bias_rat  <- calc_bias_ratio(est_vec, truth)

  # Coverage
  .cov <- if (all(c("ci_lo", "ci_hi") %in% names(ok))) {
    calc_coverage(ok$ci_lo, ok$ci_hi, truth)
  } else {
    list(estimate = NA_real_, mcse = NA_real_)
  }

  # Rejection rate (power or type I error depending on truth)
  .rej <- if ("p_value" %in% names(ok)) {
    calc_rejection_rate(ok$p_value, alpha)
  } else {
    list(estimate = NA_real_, mcse = NA_real_)
  }

  .est_finite <- est_vec[is.finite(est_vec)]

  tibble(
    n_reps           = n_total,
    n_errors         = n_errors,
    convergence_rate = conv$estimate,
    convergence_mcse = conv$mcse,
    bias             = .bias$estimate,
    bias_mcse        = .bias$mcse,
    rmse             = .rmse$estimate,
    rmse_mcse        = .rmse$mcse,
    coverage         = .cov$estimate,
    coverage_mcse    = .cov$mcse,
    rejection_rate   = .rej$estimate,
    rejection_mcse   = .rej$mcse,
    bias_ratio       = .bias_rat$estimate,
    bias_ratio_mcse  = .bias_rat$mcse,
    mean_estimate    = if (length(.est_finite) > 0) mean(.est_finite) else NA_real_,
    emp_se           = if (length(.est_finite) > 1) sd(.est_finite) else NA_real_
  )
}


# ---- Plotting Theme ---------------------------------------------------------

#' Publication-quality ggplot2 theme for benchmark figures
#'
#' Clean theme with minimal gridlines, legible text, and no default title.
#' @param base_size Base font size (default 11).
#' @return A ggplot2 theme object.
theme_benchmark <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      strip.background  = element_rect(fill = "grey95", colour = "grey70"),
      strip.text        = element_text(face = "bold", size = rel(0.9)),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(colour = "grey90", linewidth = 0.3),
      legend.position   = "bottom",
      legend.box        = "horizontal",
      legend.title      = element_text(face = "bold", size = rel(0.85)),
      legend.text       = element_text(size = rel(0.8)),
      plot.title        = element_text(face = "bold", size = rel(1.1), hjust = 0),
      plot.subtitle     = element_text(colour = "grey40", size = rel(0.9)),
      axis.title        = element_text(size = rel(0.95)),
      axis.text         = element_text(size = rel(0.85))
    )
}

#' Method colour palette for benchmark plots
#'
#' Returns a named character vector of colours, one per method label.
#' Designed for colour-blind accessibility (adapted from Okabe-Ito).
#' @return Named character vector.
method_colours <- function() {
  c(
    "Oracle"       = "#000000",
    "Tobit"        = "#E69F00",
    "AFT"          = "#56B4E9",
    "Gamma+LOD/2"  = "#009E73",
    "Gaussian+LOD/2" = "#F0E442",
    "Naive LOD/2"  = "#CC79A7",
    "Raw"          = "#CC79A7",
    "Preprocess (half)" = "#F0E442",
    "Preprocess (halfmin)" = "#56B4E9",
    "Preprocess (sqrt2)" = "#D55E00",
    "Plate-corrected" = "#009E73",
    "Plate RE" = "#0072B2",
    "McNemar exact"   = "#E69F00",
    "McNemar CC"      = "#56B4E9",
    "Chi-squared"     = "#CC79A7",
    "McNemar mid-p"   = "#009E73",
    "McNemar noCC"    = "#D55E00",
    "ANCOVA + outlier" = "#E69F00",
    "ANCOVA raw"       = "#56B4E9",
    "t-test/ANOVA"     = "#CC79A7",
    "Parametric ANCOVA" = "#009E73",
    "ANCOVA + tight fence" = "#D55E00",
    "Robust ANCOVA"    = "#0072B2",
    "PLS-DA nested CV"  = "#E69F00",
    "PLS-DA single CV"  = "#56B4E9",
    "Permutation null"  = "#CC79A7",
    # B5: Substitution methods
    "LOD/2 + Gaussian"     = "#F0E442",
    "LOD/sqrt2 + Gaussian" = "#D55E00",
    "halfmin + Gaussian"   = "#0072B2",
    "LOD + Gaussian"       = "#999999",
    "zero + Gaussian"      = "#882255",
    "LOD/2 + Gamma"        = "#009E73",
    # B6: Replicate methods
    "Naive average"        = "#E69F00",
    "QC drop_both"         = "#56B4E9",
    "QC drop_farther"      = "#009E73",
    "QC winsorize"         = "#0072B2",
    "Random effects"       = "#D55E00",
    "Ignore structure"     = "#CC79A7",
    "Oracle average"       = "#000000",
    # B0 v4: Tobit + plate FE
    "Tobit + plate FE"     = "#882255",
    # B6 v4: Tobit-based replicate methods
    "QC + Tobit"           = "#44AA99",
    "Avg + Tobit"          = "#999933"
  )
}

#' Scale colour/fill using method palette
#' @param ... Passed to ggplot2::scale_*_manual.
#' @return ggplot2 scale object.
scale_colour_method <- function(...) {
  scale_colour_manual(values = method_colours(), ...)
}

scale_fill_method <- function(...) {
  scale_fill_manual(values = method_colours(), ...)
}


# ---- Acceptance Threshold Helpers -------------------------------------------

#' Check whether a metric passes its acceptance criterion
#' @param value    Observed metric value.
#' @param lower    Lower acceptable bound (use -Inf for one-sided).
#' @param upper    Upper acceptable bound (use Inf for one-sided).
#' @return Logical.
passes_threshold <- function(value, lower = -Inf, upper = Inf) {
  !is.na(value) & value >= lower & value <= upper
}

#' Annotate a summary table with pass/fail flags for standard acceptance criteria
#' @param summary_df Aggregated summary tibble (output of \code{aggregate_cell}).
#' @param truth      True parameter value (0 = null scenario).
#' @param bias_ratio_threshold Max acceptable |bias/SE| (default 0.1).
#' @param coverage_bounds      Two-element vector c(lo, hi) for coverage (default c(0.93, 0.97)).
#' @param type1_bounds         Two-element vector for type I error under null (default c(0.025, 0.075)).
#' @return summary_df with additional logical columns: pass_bias, pass_coverage, pass_type1.
annotate_acceptance <- function(summary_df,
                                truth,
                                bias_ratio_threshold = 0.1,
                                coverage_bounds = c(0.93, 0.97),
                                type1_bounds    = c(0.025, 0.075)) {
  summary_df %>%
    mutate(
      pass_bias     = passes_threshold(abs(bias_ratio), upper = bias_ratio_threshold),
      pass_coverage = passes_threshold(coverage, coverage_bounds[1], coverage_bounds[2]),
      pass_type1    = if (truth == 0) {
        passes_threshold(rejection_rate, type1_bounds[1], type1_bounds[2])
      } else {
        NA
      }
    )
}


# ---- Grid Helpers -----------------------------------------------------------

#' Build a factorial DGP grid from named lists of parameter levels
#'
#' @param ... Named arguments, each a vector of levels. E.g.
#'   \code{n = c(20, 50), censoring = c(0.1, 0.3)}.
#' @return tibble with one row per combination and a \code{cell_id} column.
make_dgp_grid <- function(...) {
  grid <- expand.grid(..., stringsAsFactors = FALSE)
  grid$cell_id <- seq_len(nrow(grid))
  as_tibble(grid)
}

#' Generate a unique cache filename for a benchmark cell
#' @param benchmark_name Character benchmark identifier (e.g. "preprocessing").
#' @param cell_id        Integer cell ID within the grid.
#' @param method         Character method name.
#' @param cache_dir      Directory for cache files.
#' @return Character path.
cache_path <- function(benchmark_name, cell_id, method, cache_dir) {
  file.path(cache_dir, paste0(benchmark_name, "_cell", cell_id,
                               "_", gsub("[^a-zA-Z0-9]", "_", method),
                               ".rds"))
}


# ---- Summary Table Formatting -----------------------------------------------

#' Format an aggregated summary table for display/export
#'
#' Adds MCSE-annotated columns for key metrics.
#'
#' @param summary_df Aggregated summary tibble.
#' @param digits     Rounding digits (default 3).
#' @return tibble with additional character columns: bias_fmt, rmse_fmt,
#'   coverage_fmt, rejection_fmt.
format_summary_table <- function(summary_df, digits = 3) {
  summary_df %>%
    mutate(
      bias_fmt      = fmt_mcse_vec(bias, bias_mcse, digits),
      rmse_fmt      = fmt_mcse_vec(rmse, rmse_mcse, digits),
      coverage_fmt  = fmt_mcse_vec(coverage, coverage_mcse, digits),
      rejection_fmt = fmt_mcse_vec(rejection_rate, rejection_mcse, digits)
    )
}


# ---- Backbone Parameters ----------------------------------------------------
#
# Common parameter values shared across benchmarks (v3 harmonization).
# Each benchmark uses the backbone as its base and extends with
# benchmark-specific parameters. This enables direct cross-benchmark
# comparisons at shared parameter values.

BACKBONE_PARAMS <- list(
  n_subjects  = c(40L, 100L, 200L),
  effect_size = c(0, 0.5, 1.0),
  residual_sd = 1.0,
  n_reps_base = 500L
)


# ---- Panel-Level Metrics (Multi-Analyte) ------------------------------------

#' Panel false discovery rate from a single replication
#'
#' Given rejection decisions and the true signal mask for all analytes,
#' computes the false discovery proportion (FDP). Averaged across
#' replications, this estimates the panel FDR.
#'
#' @param rejected     Logical vector: was analyte k rejected (after FDR
#'   correction)?
#' @param is_signal    Logical vector: is analyte k a true signal analyte?
#' @return Scalar FDP in [0, 1]. Returns 0 if no discoveries were made
#'   (convention: FDP = 0 when R = 0).
calc_panel_fdp <- function(rejected, is_signal) {
  n_rejected <- sum(rejected, na.rm = TRUE)
  if (n_rejected == 0) return(0)
  n_false <- sum(rejected & !is_signal, na.rm = TRUE)
  n_false / n_rejected
}

#' Panel family-wise error rate indicator from a single replication
#'
#' Returns 1 if at least one null analyte was falsely rejected, 0 otherwise.
#' Averaged across replications, this estimates the FWER.
#'
#' @inheritParams calc_panel_fdp
#' @return 0 or 1.
calc_panel_fwer_indicator <- function(rejected, is_signal) {
  as.integer(any(rejected & !is_signal, na.rm = TRUE))
}

#' Panel sensitivity (true positive rate) from a single replication
#'
#' Fraction of true signal analytes that were correctly rejected.
#'
#' @inheritParams calc_panel_fdp
#' @return Scalar in [0, 1]. Returns NA if there are no signal analytes.
calc_panel_sensitivity <- function(rejected, is_signal) {
  n_signal <- sum(is_signal, na.rm = TRUE)
  if (n_signal == 0) return(NA_real_)
  sum(rejected & is_signal, na.rm = TRUE) / n_signal
}

#' Aggregate panel-level results across replications
#'
#' Takes the multi-row-per-replication output (one row per method x analyte),
#' applies BH FDR correction within each replication x method, computes
#' panel-level metrics (FDR, FWER, sensitivity), and averages across
#' replications with MCSEs.
#'
#' @param reps         data.frame with columns: rep_id, method, analyte,
#'   p_value, is_signal, estimate, se, ci_lo, ci_hi, true_effect, converged,
#'   error.
#' @param alpha        Significance threshold for BH FDR (default 0.05).
#' @param fdr_method   p.adjust method (default "BH").
#' @return tibble with one row per method, containing panel_fdr (+ mcse),
#'   panel_fwer (+ mcse), panel_sensitivity (+ mcse), plus per-analyte
#'   mean bias and coverage averaged across signal analytes.
aggregate_panel_cell <- function(reps, alpha = 0.05, fdr_method = "BH") {

  # Exclude error rows
  ok <- reps[!reps$error, , drop = FALSE]
  methods <- unique(ok$method[!is.na(ok$method)])

  panel_summaries <- list()

  for (m in methods) {
    m_data <- ok[ok$method == m, , drop = FALSE]
    rep_ids <- unique(m_data$rep_id)

    fdp_vec  <- numeric(length(rep_ids))
    fwer_vec <- numeric(length(rep_ids))
    sens_vec <- numeric(length(rep_ids))
    signal_bias_vec <- numeric(length(rep_ids))
    signal_cov_vec  <- numeric(length(rep_ids))

    for (ri in seq_along(rep_ids)) {
      r <- m_data[m_data$rep_id == rep_ids[ri], , drop = FALSE]

      # Apply BH FDR correction across analytes within this rep
      r$q_value <- p.adjust(r$p_value, method = fdr_method)
      r$rejected <- !is.na(r$q_value) & r$q_value < alpha

      fdp_vec[ri]  <- calc_panel_fdp(r$rejected, r$is_signal)
      fwer_vec[ri] <- calc_panel_fwer_indicator(r$rejected, r$is_signal)
      sens_vec[ri] <- calc_panel_sensitivity(r$rejected, r$is_signal)

      # Per-analyte bias and coverage for signal analytes only
      sig_rows <- r[r$is_signal & r$converged, , drop = FALSE]
      if (nrow(sig_rows) > 0) {
        signal_bias_vec[ri] <- mean(sig_rows$estimate - sig_rows$true_effect,
                                     na.rm = TRUE)
        covers <- (sig_rows$ci_lo <= sig_rows$true_effect) &
                  (sig_rows$ci_hi >= sig_rows$true_effect)
        signal_cov_vec[ri] <- mean(covers, na.rm = TRUE)
      } else {
        signal_bias_vec[ri] <- NA_real_
        signal_cov_vec[ri]  <- NA_real_
      }
    }

    n <- length(rep_ids)
    panel_summaries[[length(panel_summaries) + 1]] <- tibble(
      method           = m,
      n_reps           = n,
      panel_fdr        = mean(fdp_vec, na.rm = TRUE),
      panel_fdr_mcse   = sd(fdp_vec, na.rm = TRUE) / sqrt(n),
      panel_fwer       = mean(fwer_vec, na.rm = TRUE),
      panel_fwer_mcse  = sqrt(mean(fwer_vec) * (1 - mean(fwer_vec)) / n),
      panel_sensitivity      = mean(sens_vec, na.rm = TRUE),
      panel_sensitivity_mcse = sd(sens_vec, na.rm = TRUE) / sqrt(n),
      signal_mean_bias       = mean(signal_bias_vec, na.rm = TRUE),
      signal_mean_bias_mcse  = sd(signal_bias_vec, na.rm = TRUE) / sqrt(n),
      signal_mean_coverage       = mean(signal_cov_vec, na.rm = TRUE),
      signal_mean_coverage_mcse  = sd(signal_cov_vec, na.rm = TRUE) / sqrt(n)
    )
  }

  bind_rows(panel_summaries)
}


# ---- Extra-N (Supplement) Helpers ------------------------------------------

#' Read benchmark_config.yaml
#'
#' @param script_dir Directory containing benchmark_config.yaml.
#' @return Parsed list with \code{sample_sizes} and \code{extra_cell_id_offset}.
read_benchmark_config <- function(script_dir) {
  cfg_path <- file.path(script_dir, "benchmark_config.yaml")
  if (!file.exists(cfg_path)) {
    stop("Config not found: ", cfg_path)
  }
  yaml::read_yaml(cfg_path)
}

#' Get extra sample sizes for a benchmark from the config
#'
#' @param config     List returned by \code{read_benchmark_config()}.
#' @param benchmark  Config key, e.g. "B3_ancova".
#' @return Integer vector of extra sample sizes, or NULL if none defined.
get_extra_sizes <- function(config, benchmark) {
  sizes <- config$sample_sizes[[benchmark]]$extra
  if (is.null(sizes) || length(sizes) == 0) return(NULL)
  as.integer(sizes)
}

#' Filter a DGP grid to extra sample sizes and offset cell IDs
#'
#' Keeps only rows where the sample-size column matches \code{extra_sizes},
#' then offsets cell_id by \code{offset} so filenames don't collide with
#' the original v3 cache.
#'
#' @param grid         tibble from \code{make_dgp_grid()}.
#' @param size_col     Unquoted column name (e.g. \code{n_subjects}).
#' @param extra_sizes  Integer vector of sample sizes to keep.
#' @param offset       Integer offset added to cell_id (default 10000).
#' @return Filtered tibble with offset cell_ids.
filter_grid_extra <- function(grid, size_col, extra_sizes, offset = 10000L) {
  size_col_chr <- deparse(substitute(size_col))
  if (!size_col_chr %in% names(grid)) {
    stop("Column '", size_col_chr, "' not found in grid")
  }
  out <- grid[grid[[size_col_chr]] %in% extra_sizes, , drop = FALSE]
  out$cell_id <- out$cell_id + as.integer(offset)
  out
}

#' Merge new summary rows with an existing summary RDS file
#'
#' Loads the existing summary, appends new rows (deduplicating by cell_id
#' and method), and saves the result back to the same file.
#'
#' @param new_summary  tibble of summary rows from the supplement run.
#' @param summary_file Path to the existing summary .rds file.
#' @param verbose      Logical; print merge diagnostics.
#' @return Invisibly, the merged tibble.
merge_with_existing_summary <- function(new_summary, summary_file,
                                        verbose = TRUE) {
  get_usable_key <- function(df, key = "cell_id") {
    key_values <- df[[key]]
    if (is.null(key_values)) {
      return(NULL)
    }
    if (length(key_values) != nrow(df)) {
      return(NULL)
    }
    key_values
  }

  if (file.exists(summary_file)) {
    existing <- readRDS(summary_file)

    existing_ids <- get_usable_key(existing, "cell_id")
    new_ids <- get_usable_key(new_summary, "cell_id")

    if (!is.null(existing_ids) && !is.null(new_ids)) {
      existing <- existing[!(existing_ids %in% unique(new_ids)), , drop = FALSE]
    } else if (verbose) {
      missing_from <- if (is.null(existing_ids) && is.null(new_ids)) {
        "existing and new"
      } else if (is.null(existing_ids)) {
        "existing"
      } else {
        "new"
      }
      message("merge_with_existing_summary: usable cell_id unavailable from ",
              missing_from, " summary; appending without dedup")
    }
    merged <- dplyr::bind_rows(existing, new_summary)
    if (verbose) {
      message("Merged summary: ", nrow(existing), " existing + ",
              nrow(new_summary), " new = ", nrow(merged), " total rows")
    }
  } else {
    merged <- new_summary
    if (verbose) {
      message("No existing summary found; saving ", nrow(merged), " rows")
    }
  }
  saveRDS(merged, summary_file)
  invisible(merged)
}


#' Save a main-grid summary while preserving high-offset supplement rows
#'
#' Main benchmark reruns typically regenerate only the base grid (cell_ids
#' 1..N). Supplemental runs such as `--extra_n`, Phase B, or Phase C usually
#' write rows with offset cell_ids (10000+, 20000+, 30000+). This helper keeps
#' those offset rows when the main summary is refreshed, preventing later main
#' reruns from silently dropping already-completed supplements.
#'
#' @param new_summary  tibble of freshly generated main-grid summary rows.
#' @param summary_file Path to the summary .rds file.
#' @param verbose      Logical; print save diagnostics.
#' @return Invisibly, the saved tibble.
save_main_summary_preserving_supplements <- function(new_summary, summary_file,
                                                     verbose = TRUE) {
  get_usable_key <- function(df, key = "cell_id") {
    key_values <- df[[key]]
    if (is.null(key_values)) {
      return(NULL)
    }
    if (length(key_values) != nrow(df)) {
      return(NULL)
    }
    key_values
  }

  saved <- new_summary

  if (file.exists(summary_file)) {
    existing <- readRDS(summary_file)
    existing_ids <- get_usable_key(existing, "cell_id")
    new_ids <- get_usable_key(new_summary, "cell_id")

    if (!is.null(existing_ids) && !is.null(new_ids) && length(new_ids) > 0) {
      main_max <- suppressWarnings(max(new_ids, na.rm = TRUE))
      if (is.finite(main_max)) {
        supplement_rows <- existing[existing_ids > main_max, , drop = FALSE]
        saved <- dplyr::bind_rows(new_summary, supplement_rows)
        if (verbose) {
          message("Preserved ", nrow(supplement_rows),
                  " supplemental rows while saving main summary")
        }
      } else if (verbose) {
        message("save_main_summary_preserving_supplements: main summary has no finite cell_id; saving new rows only")
      }
    } else if (verbose) {
      message("save_main_summary_preserving_supplements: usable cell_id unavailable; saving new rows only")
    }
  }

  saved <- carry_forward_small_n_from_prev_version(saved, summary_file)

  saveRDS(saved, summary_file)
  if (verbose) {
    message("Summary saved -> ", summary_file, " (", nrow(saved), " rows)")
  }
  invisible(saved)
}


# =============================================================================
# Small-n supplement registry + validation
# =============================================================================
#
# Small-n runs (submit_benchmark_*_small_n.sh) add cells with cell_id offset
# +10000 on top of the main grid. They are the single most common source of
# "silent regression" when CACHE_VERSION is bumped: the main grid gets rerun
# under the new prefix, but the small-n jobs are not, and the new summary
# ends up missing n=10,20 rows that older downstream figures showed.
#
# The registry below is the source of truth for which benchmarks should have
# small-n rows and at which sample sizes. Two helpers use it:
#
#   1. validate_small_n_supplements() — called from plot_*.R. Warns (or
#      stops, with `strict = TRUE`) if the summary lacks expected small-n
#      rows. Keeps the regression loud instead of silent.
#   2. carry_forward_small_n_from_prev_version() — called from
#      save_main_summary_preserving_supplements(). When the current-version
#      summary has zero offset-10k rows but a predecessor-version summary
#      on disk has them, carries them forward with a LOUD warning so the
#      operator knows rerunning the small-n job is still the right thing.

#' Registry of expected small-n supplements per summary prefix.
#'
#' `prefix` is the cache basename without the trailing `_summary.rds` — e.g.
#' the active preprocessing summary is `preprocessing_v4_summary.rds`, so
#' the prefix is `preprocessing_v4`.
#'
#' `n_col` is the sample-size column name inside that summary.
#' `extra_n` is the set of sample sizes the small-n submit script runs.
#' `prev_prefix` is the previous-version summary to carry forward from when
#' no fresh supplement exists (NULL = no predecessor).
small_n_registry <- function() {
  list(
    preprocessing_v4 = list(n_col = "n_subjects",  extra_n = c(10L, 20L),
                            prev_prefix = "preprocessing_v3"),
    censoring_v3     = list(n_col = "n_subjects",  extra_n = c(10L),
                            prev_prefix = NULL),
    ancova_v3        = list(n_col = "n_per_group", extra_n = c(10L),
                            prev_prefix = NULL),
    mcnemar_v2       = list(n_col = "n_subjects",  extra_n = c(10L, 20L),
                            prev_prefix = NULL),
    plsda_v3         = list(n_col = "n_per_group", extra_n = c(10L),
                            prev_prefix = NULL),
    replicates_v4    = list(n_col = "n_subjects",  extra_n = c(10L, 20L),
                            prev_prefix = "replicates_v3"),
    substitution_v3  = list(n_col = "n_subjects",  extra_n = c(10L, 20L),
                            prev_prefix = NULL)
  )
}

#' Check that a summary has the small-n rows the registry expects.
#'
#' Looks up `prefix` in `small_n_registry()`. If the registry entry is missing
#' the function is a no-op. Otherwise it warns (or stops, with
#' `strict = TRUE`) if every expected small-n value is absent from the summary.
#'
#' @param summary_df  Aggregated summary tibble.
#' @param prefix      Character. Summary prefix (e.g. "preprocessing_v4").
#' @param strict      If TRUE, raise an error instead of warning.
#' @return invisibly, a list describing the check outcome.
validate_small_n_supplements <- function(summary_df, prefix, strict = FALSE) {
  reg <- small_n_registry()
  if (!prefix %in% names(reg)) {
    return(invisible(list(status = "no_registry_entry", prefix = prefix)))
  }
  entry <- reg[[prefix]]
  n_col <- entry$n_col
  if (!n_col %in% names(summary_df)) {
    msg <- sprintf("validate_small_n_supplements: column '%s' not in %s summary",
                   n_col, prefix)
    if (strict) stop(msg) else warning(msg)
    return(invisible(list(status = "missing_column", prefix = prefix)))
  }
  present <- unique(summary_df[[n_col]])
  missing <- setdiff(entry$extra_n, present)
  if (length(missing) == 0) {
    return(invisible(list(status = "ok", prefix = prefix,
                          present_small_n = intersect(entry$extra_n, present))))
  }
  msg <- sprintf(
    "%s summary is missing small-n rows: expected %s in column '%s', present %s.\n  Re-run submit_benchmark_%s_small_n.sh, OR confirm the registry entry in helpers_benchmark.R is still correct.",
    prefix,
    paste(entry$extra_n, collapse = ", "),
    n_col,
    paste(sort(present), collapse = ", "),
    sub("_v[0-9]+$", "", prefix)
  )
  if (strict) stop(msg) else warning(msg, call. = FALSE)
  invisible(list(status = "missing_small_n", prefix = prefix,
                 missing = missing, present = sort(present)))
}

#' Carry forward offset-10k rows from a predecessor-version summary.
#'
#' Called from save_main_summary_preserving_supplements() when the current
#' summary has no offset-10k (small-n) rows but the registry says a
#' predecessor exists and has them. Only rows whose `method` values are
#' present in the current summary are carried forward — this prevents
#' grafting methods that were removed between versions. A loud warning is
#' emitted so operators know to rerun the small-n job at the new version.
#'
#' @param current_df  Current-version summary tibble (after preserving
#'                    in-version supplements).
#' @param summary_file Full path to the current summary file. Used to find
#'                    the predecessor file in the same directory and to
#'                    derive the prefix via basename.
#' @return current_df with predecessor offset-10k rows appended where
#'   methods intersect.
carry_forward_small_n_from_prev_version <- function(current_df, summary_file) {
  prefix <- sub("_summary\\.rds$", "", basename(summary_file))
  reg <- small_n_registry()
  if (!prefix %in% names(reg)) return(current_df)
  entry <- reg[[prefix]]
  if (is.null(entry$prev_prefix)) return(current_df)

  cur_ids <- current_df[["cell_id"]]
  if (is.null(cur_ids)) return(current_df)
  has_supp_now <- any(cur_ids >= 10000 & cur_ids < 20000, na.rm = TRUE)
  if (has_supp_now) return(current_df)

  prev_file <- file.path(dirname(summary_file),
                         paste0(entry$prev_prefix, "_summary.rds"))
  if (!file.exists(prev_file)) return(current_df)
  prev_df <- tryCatch(readRDS(prev_file), error = function(e) NULL)
  if (is.null(prev_df) || !"cell_id" %in% names(prev_df)) return(current_df)

  supp_rows <- prev_df[prev_df$cell_id >= 10000 & prev_df$cell_id < 20000, , drop = FALSE]
  if (nrow(supp_rows) == 0) return(current_df)

  if ("method" %in% names(supp_rows) && "method" %in% names(current_df)) {
    keep <- supp_rows$method %in% unique(current_df$method)
    supp_rows <- supp_rows[keep, , drop = FALSE]
  }
  if (nrow(supp_rows) == 0) return(current_df)

  keep_cols <- intersect(names(current_df), names(supp_rows))
  supp_rows <- supp_rows[, keep_cols, drop = FALSE]

  warning(sprintf(
    paste0("CARRY-FORWARD: %s has no small-n supplement; grafting %d rows ",
           "from %s. Rerun submit_benchmark_%s_small_n.sh at the current ",
           "version to replace these rows with fresh fits."),
    prefix, nrow(supp_rows), basename(prev_file),
    sub("_v[0-9]+$", "", prefix)
  ), call. = FALSE)

  dplyr::bind_rows(current_df, supp_rows)
}
