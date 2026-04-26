#!/usr/bin/env Rscript
# benchmark_mcnemar.R — Benchmark 2: McNemar Calibration (v2)
#
# ADEMP Structure:
#   Aims:    Evaluate type I error calibration, power, and FDR control of
#            five testing variants under realistic paired detection data:
#            exact McNemar, McNemar with/without continuity correction,
#            mid-p McNemar, and marginal chi-squared. Assess robustness
#            across a range of within-subject temporal correlations.
#   DGP:     Latent bivariate normal threshold model with asymmetric detection
#            shift for signal analytes. Factorial over sample size, baseline
#            detection rate, detection shift (delta), analyte correlation,
#            and temporal correlation (rho_t).
#   Estimand: Change in detection probability (delta) between paired timepoints.
#   Methods: (1) Exact McNemar via exact2x2, (2) McNemar with continuity
#            correction, (3) Marginal chi-squared (ignoring pairing),
#            (4) Mid-p McNemar, (5) McNemar without continuity correction.
#   Performance: type I error, power, size-adjusted power, CI coverage for
#                delta, MPOR bias (Haldane-corrected), FDR control at q = 0.05.
#
# Revisions from v1 (per B2 report section 6 recommendations):
#   - Added mid-p McNemar (Method 4) to address exact-test conservatism
#   - Added asymptotic McNemar without CC (Method 5) for completeness
#   - Varied rho_T over {0, 0.3, 0.5, 0.7} to stress-test marginal chi-squared
#   - Increased default reps from 200 to 1000 for sharper MCSE
#   - Haldane MPOR correction (cc+0.5)/(b+0.5) avoids undefined ratios
#   - Size-adjusted power panel (calibrated to exact 5% type I error)
#   - MPOR computability diagnostic tracks fraction with b > 0
#   - N_ANALYTES and N_SIGNAL configurable via CLI
#
# Usage:
#   Rscript benchmark_mcnemar.R [--n_reps N] [--n_cores N] [--cache_dir DIR]
#                                [--n_analytes N] [--n_signal N]
#                                [--extra_n N1,N2,...] [--phase_b TRUE]
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
  library(exact2x2)
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

N_REPS     <- as.integer(parse_cli_arg("--n_reps",     "1000"))
N_CORES    <- as.integer(parse_cli_arg("--n_cores",    "1"))
CACHE_DIR  <- parse_cli_arg("--cache_dir", file.path(script_dir, "cache"))
EXTRA_N    <- parse_cli_arg("--extra_n", NULL)
PHASE_B    <- parse_cli_arg("--phase_b", NULL)
S1_MODE    <- parse_cli_arg("--s1", NULL)
N_ANALYTES <- as.integer(parse_cli_arg("--n_analytes", "10"))
N_SIGNAL   <- as.integer(parse_cli_arg("--n_signal",   "3"))
BASE_SEED  <- 20240201L
BENCHMARK  <- "mcnemar_v2"
VERBOSE    <- TRUE

# Extra-n supplement mode
EXTRA_SIZES  <- if (!is.null(EXTRA_N)) as.integer(strsplit(EXTRA_N, ",")[[1]]) else NULL
EXTRA_OFFSET <- 10000L

# Phase B supplement mode (negative temporal correlation)
PHASE_B_OFFSET <- 20000L

# S1 flow-through validation mode (review plan chat 2, section 6c)
S1_OFFSET <- 30000L

# Fixed design parameters
P2_CAP     <- 0.99    # maximum detection probability at T2

if (!dir.exists(CACHE_DIR)) dir.create(CACHE_DIR, recursive = TRUE)


# ---- DGP Grid ---------------------------------------------------------------

dgp_grid <- make_dgp_grid(
  n_subjects         = if (!is.null(EXTRA_SIZES)) EXTRA_SIZES else c(15L, 30L, 50L, 100L, 200L),
  baseline_detection = c(0.30, 0.50, 0.80),
  delta_rate         = c(0, 0.05, 0.15, 0.25),
  analyte_rho        = c(0, 0.3, 0.6),
  rho_t              = c(0, 0.3, 0.5, 0.7)
)
if (!is.null(EXTRA_SIZES)) dgp_grid$cell_id <- dgp_grid$cell_id + EXTRA_OFFSET

if (VERBOSE) {
  message("=== Benchmark 2: McNemar Calibration (v2) ===")
  message("DGP grid: ", nrow(dgp_grid), " cells x ", N_REPS, " reps = ",
          nrow(dgp_grid) * N_REPS, " total runs")
  message("Analytes per sim: ", N_ANALYTES, " (", N_SIGNAL, " signal)")
  message("Temporal correlation grid: 0, 0.3, 0.5, 0.7")
  message("Methods: 5 (exact, CC, chi-sq, mid-p, noCC)")
  message("Cache dir: ", CACHE_DIR)
}


# ---- Bivariate Normal CDF Helper -------------------------------------------

# P(X < a, Y < b | rho) for standard bivariate normal
.bvnorm_cdf_local <- function(a, b, rho) {
  if (abs(rho) < 1e-10) return(stats::pnorm(a) * stats::pnorm(b))

  if (requireNamespace("mvtnorm", quietly = TRUE)) {
    sigma <- matrix(c(1, rho, rho, 1), 2, 2)
    return(mvtnorm::pmvnorm(upper = c(a, b), sigma = sigma)[1])
  }

  # Fallback: numerical integration via conditional distribution
  sq <- sqrt(1 - rho^2)
  f <- function(x) stats::pnorm((b - rho * x) / sq) * stats::dnorm(x)
  stats::integrate(f, lower = -8, upper = a, subdivisions = 200,
                   rel.tol = 1e-8)$value
}


# ---- True MPOR Computation --------------------------------------------------

#' Compute theoretical MPOR from the bivariate normal threshold model.
#' Under H0 (delta=0), true MPOR = 1. Under H1 (delta>0), MPOR > 1.
#'
#' @param baseline Detection probability at T1.
#' @param delta_rate Detection shift for signal analytes (0-1 scale).
#' @param rho_t Within-subject temporal correlation.
#' @return True matched-pair odds ratio (gain/loss).
compute_true_mpor <- function(baseline, delta_rate, rho_t) {
  if (delta_rate == 0) return(1.0)

  p1 <- baseline
  p2 <- min(baseline + delta_rate, P2_CAP)
  t1 <- stats::qnorm(1 - p1)
  t2 <- stats::qnorm(1 - p2)

  phi2 <- .bvnorm_cdf_local(t1, t2, rho_t)

  p_gain <- stats::pnorm(t1) - phi2   # P(Z1 < t1, Z2 > t2)
  p_loss <- stats::pnorm(t2) - phi2   # P(Z1 > t1, Z2 < t2)

  if (p_loss < 1e-10) return(NA_real_)
  p_gain / p_loss
}


# ---- Per-Analyte Method Helpers ---------------------------------------------

#' Wald CI for paired proportion difference delta = (c - b) / n
paired_wald_ci <- function(b, cc, n, conf_level = 0.95) {
  delta <- (cc - b) / n
  var_di <- (b + cc) / n - delta^2
  se <- sqrt(max(var_di / n, 0))
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  list(
    delta = delta * 100,
    ci_lo = (delta - z * se) * 100,
    ci_hi = (delta + z * se) * 100
  )
}


#' Exact McNemar test for one analyte (Method 1)
exact_mcnemar_analyte <- function(a, b, cc, d, n, conf_level = 0.95) {
  n_disc <- b + cc
  delta  <- (cc - b) / n * 100

  if (n_disc == 0) {
    return(list(p_value = 1.0, delta = 0, ci_lo = 0, ci_hi = 0,
                mpor = 1.0, mpor_defined = TRUE, n_discordant = 0L))
  }

  tab <- matrix(c(a, b, cc, d), nrow = 2, byrow = TRUE)

  # Exact McNemar p-value
  ex <- tryCatch(
    exact2x2::mcnemar.exact(tab, conf.level = conf_level),
    error = function(e) NULL
  )
  if (!is.null(ex)) {
    p_value <- ex$p.value
  } else {
    # Fallback: binomial test on discordant pairs
    bt <- stats::binom.test(cc, n_disc, p = 0.5)
    p_value <- bt$p.value
  }

  # Exact CI for delta via mcnemarExactDP
  dp <- tryCatch(
    exact2x2::mcnemarExactDP(n = n, m = n_disc, x = cc,
                              conf.level = conf_level),
    error = function(e) NULL
  )
  if (!is.null(dp)) {
    ci_lo <- dp$conf.int[1] * 100
    ci_hi <- dp$conf.int[2] * 100
  } else {
    wald <- paired_wald_ci(b, cc, n, conf_level)
    ci_lo <- wald$ci_lo
    ci_hi <- wald$ci_hi
  }

  # Haldane-corrected MPOR (always finite)
  mpor <- (cc + 0.5) / (b + 0.5)
  mpor_defined <- b > 0

  list(p_value = p_value, delta = delta, ci_lo = ci_lo, ci_hi = ci_hi,
       mpor = mpor, mpor_defined = mpor_defined,
       n_discordant = as.integer(n_disc))
}


#' McNemar with continuity correction for one analyte (Method 2)
mcnemar_cc_analyte <- function(a, b, cc, d, n, conf_level = 0.95) {
  n_disc <- b + cc
  delta  <- (cc - b) / n * 100

  if (n_disc == 0) {
    return(list(p_value = 1.0, delta = 0, ci_lo = 0, ci_hi = 0,
                mpor = 1.0, mpor_defined = TRUE, n_discordant = 0L))
  }

  chi2 <- max((abs(b - cc) - 1)^2 / (b + cc), 0)
  p_value <- stats::pchisq(chi2, df = 1, lower.tail = FALSE)

  wald <- paired_wald_ci(b, cc, n, conf_level)
  mpor <- (cc + 0.5) / (b + 0.5)
  mpor_defined <- b > 0

  list(p_value = p_value, delta = delta,
       ci_lo = wald$ci_lo, ci_hi = wald$ci_hi,
       mpor = mpor, mpor_defined = mpor_defined,
       n_discordant = as.integer(n_disc))
}


#' Marginal chi-squared test for one analyte (Method 3)
#' Treats T1 and T2 as independent samples (ignoring pairing).
marginal_chisq_analyte <- function(a, b, cc, d, n, conf_level = 0.95) {
  delta <- (cc - b) / n * 100
  n_disc <- b + cc

  n_det_T1 <- a + b
  n_det_T2 <- a + cc

  pt <- tryCatch(
    stats::prop.test(c(n_det_T2, n_det_T1), c(n, n),
                      conf.level = conf_level, correct = FALSE),
    error = function(e) NULL
  )

  if (!is.null(pt)) {
    p_value <- pt$p.value
    ci_lo   <- pt$conf.int[1] * 100
    ci_hi   <- pt$conf.int[2] * 100
  } else {
    p_value <- NA_real_
    ci_lo   <- NA_real_
    ci_hi   <- NA_real_
  }

  mpor <- (cc + 0.5) / (b + 0.5)
  mpor_defined <- b > 0

  list(p_value = p_value, delta = delta, ci_lo = ci_lo, ci_hi = ci_hi,
       mpor = mpor, mpor_defined = mpor_defined,
       n_discordant = as.integer(n_disc))
}


#' Mid-p McNemar test for one analyte (Method 4)
#'
#' Halves the excess conservatism of the exact test by splitting the
#' probability mass at the observed value between rejection and non-rejection.
#' Under H0, cc ~ Binom(n_disc, 0.5); mid-p subtracts the point probability
#' from the exact two-sided p-value.
midp_mcnemar_analyte <- function(a, b, cc, d, n, conf_level = 0.95) {
  n_disc <- b + cc
  delta  <- (cc - b) / n * 100

  if (n_disc == 0) {
    return(list(p_value = 1.0, delta = 0, ci_lo = 0, ci_hi = 0,
                mpor = 1.0, mpor_defined = TRUE, n_discordant = 0L))
  }

  # Exact two-sided p-value from binomial test
  p_exact <- stats::binom.test(cc, n_disc, p = 0.5)$p.value
  # Mid-p correction: subtract point mass at observed value
  p_value <- max(p_exact - stats::dbinom(cc, n_disc, 0.5), 0)

  wald <- paired_wald_ci(b, cc, n, conf_level)
  mpor <- (cc + 0.5) / (b + 0.5)
  mpor_defined <- b > 0

  list(p_value = p_value, delta = delta,
       ci_lo = wald$ci_lo, ci_hi = wald$ci_hi,
       mpor = mpor, mpor_defined = mpor_defined,
       n_discordant = as.integer(n_disc))
}


#' McNemar without continuity correction for one analyte (Method 5)
#'
#' Standard asymptotic McNemar chi-squared: chi2 = (b - c)^2 / (b + c).
#' Less conservative than the CC version; widely used in practice.
mcnemar_nocc_analyte <- function(a, b, cc, d, n, conf_level = 0.95) {
  n_disc <- b + cc
  delta  <- (cc - b) / n * 100

  if (n_disc == 0) {
    return(list(p_value = 1.0, delta = 0, ci_lo = 0, ci_hi = 0,
                mpor = 1.0, mpor_defined = TRUE, n_discordant = 0L))
  }

  chi2 <- (b - cc)^2 / (b + cc)
  p_value <- stats::pchisq(chi2, df = 1, lower.tail = FALSE)

  wald <- paired_wald_ci(b, cc, n, conf_level)
  mpor <- (cc + 0.5) / (b + 0.5)
  mpor_defined <- b > 0

  list(p_value = p_value, delta = delta,
       ci_lo = wald$ci_lo, ci_hi = wald$ci_hi,
       mpor = mpor, mpor_defined = mpor_defined,
       n_discordant = as.integer(n_disc))
}


# ---- S1: Wrapper Flow-Through Helpers ---------------------------------------
#
# S1 checks that B2's calibration guarantees on the 2x2 statistical kernel
# carry over to the exported wrapper mcnemar_detection(). The standalone
# helpers above operate on raw cell counts (a, b, cc, d); the wrapper
# additionally runs preprocessing — long-to-wide pivot, replicate
# aggregation via majority vote, NA filtering, timepoint defaulting. S1
# pipes the same DGP through both paths at matched seeds so per-rep
# parity and cell-aggregate drift can both be measured.

#' Call mcnemar_detection() and return per-analyte rows in the same schema
#' used by exact_mcnemar_analyte() / midp_mcnemar_analyte() / mcnemar_nocc_analyte(),
#' but as a single data.frame (one row per analyte) rather than a list per call.
#'
#' The S1 DGP produces exactly one observation per (subject, cytokine, timepoint),
#' so the wrapper's majority-vote replicate aggregation is an identity transform
#' in this harness — per-rep parity with the standalone helpers is expected to
#' hold to machine precision for the kernel calls (mcnemar.exact, binom.test
#' mid-p, pchisq noCC).
wrapper_per_analyte <- function(long_data, baseline, comparison, method,
                                conf_level = 0.95) {
  res <- mcnemar_detection(
    data           = long_data,
    subject_col    = "subject_id",
    cytokine_col   = "cytokine",
    timepoint_col  = "timepoint",
    censoring_col  = "cens_lod",
    baseline       = baseline,
    comparison     = comparison,
    conf_level     = conf_level,
    fdr_method     = "BH",
    method         = method,
    quiet          = TRUE
  )

  data.frame(
    analyte      = as.character(res$cytokine),
    p_value      = res$p_mcnemar,
    delta        = res$delta_detection,
    ci_lo        = res$delta_ci_lo,
    ci_hi        = res$delta_ci_hi,
    mpor         = res$mpor,
    mpor_defined = !is.na(res$mpor) & res$loss > 0,
    n_discordant = as.integer(res$n_discordant),
    stringsAsFactors = FALSE
  )
}


#' Run one replication through the wrapper path.
#'
#' Regenerates the same latent bivariate normal DGP as run_one_rep() (so
#' identical set.seed(seed) produces identical 2x2 counts per analyte),
#' assembles the detection booleans into long-format data with columns
#' (subject_id, cytokine, timepoint, cens_lod = !detected), and pipes
#' through wrapper_per_analyte() for each of method in {"exact","midp","noCC"}.
#'
#' @return data.frame with 3 * N_ANALYTES rows, schema aligned to run_one_rep()
#'   for direct rbind after tagging a `path` column upstream.
run_one_rep_via_wrapper <- function(seed, n_subjects, baseline_detection,
                                    delta_rate, analyte_rho, rho_t) {

  tryCatch({
    set.seed(seed)

    # --- Analyte setup (identical to run_one_rep) ---
    analyte_names <- paste0("A", sprintf("%02d", seq_len(N_ANALYTES)))
    signal_mask   <- rep(FALSE, N_ANALYTES)
    signal_idx    <- sort(sample.int(N_ANALYTES, N_SIGNAL))
    signal_mask[signal_idx] <- TRUE

    # --- Thresholds ---
    t1 <- stats::qnorm(1 - baseline_detection)
    p2_signal <- min(baseline_detection + delta_rate, P2_CAP)
    t2_signal <- stats::qnorm(1 - p2_signal)
    actual_delta_pct <- (p2_signal - baseline_detection) * 100
    true_mpor_signal <- compute_true_mpor(baseline_detection, delta_rate, rho_t)

    # --- Generate latent bivariate normal ---
    stopifnot(abs(rho_t) < 1,
              analyte_rho > -1 / (N_ANALYTES - 1),
              analyte_rho < 1)
    cor_mat <- diag(N_ANALYTES)
    if (analyte_rho != 0) {
      cor_mat[lower.tri(cor_mat)] <- analyte_rho
      cor_mat[upper.tri(cor_mat)] <- analyte_rho
    }
    chol_mat <- chol(cor_mat)

    Z1_raw   <- matrix(stats::rnorm(n_subjects * N_ANALYTES),
                       n_subjects, N_ANALYTES)
    Z1       <- Z1_raw %*% chol_mat
    Z2_innov <- matrix(stats::rnorm(n_subjects * N_ANALYTES),
                       n_subjects, N_ANALYTES)
    Z2_innov <- Z2_innov %*% chol_mat
    Z2       <- rho_t * Z1 + sqrt(1 - rho_t^2) * Z2_innov

    # --- Threshold to detection ---
    det_T1 <- Z1 > t1
    det_T2 <- matrix(FALSE, n_subjects, N_ANALYTES)
    for (k in seq_len(N_ANALYTES)) {
      threshold_k <- if (signal_mask[k]) t2_signal else t1
      det_T2[, k] <- Z2[, k] > threshold_k
    }

    # --- Assemble long-format data for the wrapper ---
    #   Rows: n_subjects * N_ANALYTES * 2 timepoints
    #   Timepoint labels "T1"/"T2" passed explicitly to mcnemar_detection() so
    #   the alphabetical-default path is not taken (it would still pick T1/T2
    #   correctly, but we suppress the message(...) hint by being explicit).
    long_data <- data.frame(
      subject_id = rep(seq_len(n_subjects), times = N_ANALYTES * 2L),
      cytokine   = rep(rep(analyte_names, each = n_subjects), times = 2L),
      timepoint  = rep(c("T1", "T2"), each = n_subjects * N_ANALYTES),
      cens_lod   = c(as.vector(!det_T1), as.vector(!det_T2)),
      stringsAsFactors = FALSE
    )

    # --- Call wrapper for each method ---
    methods <- c("exact", "midp", "noCC")
    method_labels <- c(exact = "McNemar exact",
                       midp  = "McNemar mid-p",
                       noCC  = "McNemar noCC")

    per_method <- lapply(methods, function(m) {
      w <- wrapper_per_analyte(long_data, baseline = "T1", comparison = "T2",
                               method = m, conf_level = 0.95)
      # Align by analyte name to recover is_signal and truths
      k_idx <- match(w$analyte, analyte_names)
      data.frame(
        method       = method_labels[[m]],
        analyte      = w$analyte,
        is_signal    = signal_mask[k_idx],
        p_value      = w$p_value,
        q_value      = NA_real_,
        delta        = w$delta,
        delta_ci_lo  = w$ci_lo,
        delta_ci_hi  = w$ci_hi,
        mpor         = w$mpor,
        mpor_defined = w$mpor_defined,
        n_discordant = w$n_discordant,
        true_delta   = ifelse(signal_mask[k_idx], actual_delta_pct, 0),
        true_mpor    = ifelse(signal_mask[k_idx], true_mpor_signal, 1.0),
        stringsAsFactors = FALSE
      )
    })

    out <- do.call(rbind, per_method)

    # Apply BH FDR correction per method (match run_one_rep() semantics)
    for (lbl in unname(method_labels)) {
      sel <- out$method == lbl
      out$q_value[sel] <- p.adjust(out$p_value[sel], "BH")
    }
    out

  }, error = function(e) {
    data.frame(
      method = NA_character_, analyte = NA_character_,
      is_signal = NA, p_value = NA_real_, q_value = NA_real_,
      delta = NA_real_, delta_ci_lo = NA_real_, delta_ci_hi = NA_real_,
      mpor = NA_real_, mpor_defined = NA,
      n_discordant = NA_integer_,
      true_delta = NA_real_, true_mpor = NA_real_,
      stringsAsFactors = FALSE
    )
  })
}


# ---- Single-Replication Function --------------------------------------------

#' Run one replication for a single DGP cell.
#'
#' Generates paired detection data from a latent bivariate normal model,
#' computes 2x2 tables per analyte, runs all 5 methods, and returns
#' a data.frame with one row per (method x analyte).
#'
#' @param seed             RNG seed.
#' @param n_subjects       Sample size.
#' @param baseline_detection Marginal detection probability at T1.
#' @param delta_rate       Net detection shift for signal analytes (0-1 scale).
#' @param analyte_rho      Between-analyte correlation (shared latent factor).
#' @param rho_t            Within-subject temporal correlation.
#' @return data.frame with 5 * N_ANALYTES rows.
run_one_rep <- function(seed, n_subjects, baseline_detection, delta_rate,
                        analyte_rho, rho_t) {

  # Error-safe wrapper: always return a valid data.frame
  tryCatch({
    set.seed(seed)

    # --- Analyte setup ---
    analyte_names <- paste0("A", sprintf("%02d", seq_len(N_ANALYTES)))
    signal_mask   <- rep(FALSE, N_ANALYTES)
    signal_idx    <- sort(sample.int(N_ANALYTES, N_SIGNAL))
    signal_mask[signal_idx] <- TRUE

    # --- Thresholds ---
    t1 <- stats::qnorm(1 - baseline_detection)
    p2_signal <- min(baseline_detection + delta_rate, P2_CAP)
    t2_signal <- stats::qnorm(1 - p2_signal)
    actual_delta_pct <- (p2_signal - baseline_detection) * 100

    # True MPOR for signal analytes
    true_mpor_signal <- compute_true_mpor(baseline_detection, delta_rate, rho_t)

    # --- Generate latent bivariate normal ---
    # Joint (Z1, Z2) covariance = A ⊗ [[1, rho_t],[rho_t, 1]], PD iff A is PD
    # and |rho_t| < 1. A (equicorrelation) is PD iff -1/(N-1) < analyte_rho < 1.
    stopifnot(abs(rho_t) < 1,
              analyte_rho > -1 / (N_ANALYTES - 1),
              analyte_rho < 1)
    cor_mat <- diag(N_ANALYTES)
    if (analyte_rho != 0) {
      cor_mat[lower.tri(cor_mat)] <- analyte_rho
      cor_mat[upper.tri(cor_mat)] <- analyte_rho
    }
    chol_mat <- chol(cor_mat)

    Z1_raw   <- matrix(stats::rnorm(n_subjects * N_ANALYTES),
                        n_subjects, N_ANALYTES)
    Z1       <- Z1_raw %*% chol_mat

    Z2_innov <- matrix(stats::rnorm(n_subjects * N_ANALYTES),
                        n_subjects, N_ANALYTES)
    Z2_innov <- Z2_innov %*% chol_mat
    Z2       <- rho_t * Z1 + sqrt(1 - rho_t^2) * Z2_innov

    # --- Threshold to detection ---
    det_T1 <- Z1 > t1
    det_T2 <- matrix(FALSE, n_subjects, N_ANALYTES)
    for (k in seq_len(N_ANALYTES)) {
      threshold_k <- if (signal_mask[k]) t2_signal else t1
      det_T2[, k] <- Z2[, k] > threshold_k
    }

    # --- Compute 2x2 tables and run methods per analyte ---
    n_methods <- 5L
    results_list <- vector("list", N_ANALYTES * n_methods)
    p_values_m1  <- numeric(N_ANALYTES)
    p_values_m2  <- numeric(N_ANALYTES)
    p_values_m3  <- numeric(N_ANALYTES)
    p_values_m4  <- numeric(N_ANALYTES)
    p_values_m5  <- numeric(N_ANALYTES)
    idx <- 0L

    for (k in seq_len(N_ANALYTES)) {
      a  <- sum( det_T1[, k] &  det_T2[, k])
      b  <- sum( det_T1[, k] & !det_T2[, k])   # loss
      cc <- sum(!det_T1[, k] &  det_T2[, k])    # gain
      d  <- sum(!det_T1[, k] & !det_T2[, k])
      n  <- n_subjects

      true_delta_k <- if (signal_mask[k]) actual_delta_pct else 0
      true_mpor_k  <- if (signal_mask[k]) true_mpor_signal else 1.0

      make_row <- function(method_name, res) {
        data.frame(
          method = method_name, analyte = analyte_names[k],
          is_signal = signal_mask[k],
          p_value = res$p_value, q_value = NA_real_,
          delta = res$delta, delta_ci_lo = res$ci_lo, delta_ci_hi = res$ci_hi,
          mpor = res$mpor, mpor_defined = res$mpor_defined,
          n_discordant = res$n_discordant,
          true_delta = true_delta_k, true_mpor = true_mpor_k,
          stringsAsFactors = FALSE
        )
      }

      # Method 1: Exact McNemar
      m1 <- exact_mcnemar_analyte(a, b, cc, d, n)
      p_values_m1[k] <- m1$p_value
      idx <- idx + 1L
      results_list[[idx]] <- make_row("McNemar exact", m1)

      # Method 2: McNemar CC
      m2 <- mcnemar_cc_analyte(a, b, cc, d, n)
      p_values_m2[k] <- m2$p_value
      idx <- idx + 1L
      results_list[[idx]] <- make_row("McNemar CC", m2)

      # Method 3: Marginal chi-squared
      m3 <- marginal_chisq_analyte(a, b, cc, d, n)
      p_values_m3[k] <- m3$p_value
      idx <- idx + 1L
      results_list[[idx]] <- make_row("Chi-squared", m3)

      # Method 4: Mid-p McNemar
      m4 <- midp_mcnemar_analyte(a, b, cc, d, n)
      p_values_m4[k] <- m4$p_value
      idx <- idx + 1L
      results_list[[idx]] <- make_row("McNemar mid-p", m4)

      # Method 5: McNemar noCC
      m5 <- mcnemar_nocc_analyte(a, b, cc, d, n)
      p_values_m5[k] <- m5$p_value
      idx <- idx + 1L
      results_list[[idx]] <- make_row("McNemar noCC", m5)
    }

    results <- do.call(rbind, results_list)

    # Apply BH FDR correction per method (across analytes within this rep)
    results$q_value[results$method == "McNemar exact"]  <- p.adjust(p_values_m1, "BH")
    results$q_value[results$method == "McNemar CC"]     <- p.adjust(p_values_m2, "BH")
    results$q_value[results$method == "Chi-squared"]    <- p.adjust(p_values_m3, "BH")
    results$q_value[results$method == "McNemar mid-p"]  <- p.adjust(p_values_m4, "BH")
    results$q_value[results$method == "McNemar noCC"]   <- p.adjust(p_values_m5, "BH")

    results

  }, error = function(e) {
    # Return single-row failure indicator with all expected columns
    data.frame(
      method = NA_character_, analyte = NA_character_,
      is_signal = NA, p_value = NA_real_, q_value = NA_real_,
      delta = NA_real_, delta_ci_lo = NA_real_, delta_ci_hi = NA_real_,
      mpor = NA_real_, mpor_defined = NA,
      n_discordant = NA_integer_,
      true_delta = NA_real_, true_mpor = NA_real_,
      stringsAsFactors = FALSE
    )
  })
}


# ---- Benchmark Runner -------------------------------------------------------

run_benchmark_mcnemar <- function(dgp_grid, n_reps, base_seed,
                                   cache_dir, n_cores, verbose) {

  all_raw <- list()

  for (i in seq_len(nrow(dgp_grid))) {
    row <- dgp_grid[i, ]
    cell_label <- paste0("cell_", row$cell_id,
                          " [n=", row$n_subjects,
                          "/base=", row$baseline_detection,
                          "/delta=", row$delta_rate,
                          "/rho=", row$analyte_rho,
                          "/rho_t=", row$rho_t, "]")

    # Closure over DGP parameters (use local() for safe loop binding)
    rep_fn <- local({
      ns <- row$n_subjects
      bd <- row$baseline_detection
      dr <- row$delta_rate
      ar <- row$analyte_rho
      rt <- row$rho_t
      function(seed) {
        run_one_rep(seed = seed, n_subjects = ns, baseline_detection = bd,
                    delta_rate = dr, analyte_rho = ar, rho_t = rt)
      }
    })

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

    # Tag with DGP parameters
    raw_results$n_subjects         <- row$n_subjects
    raw_results$baseline_detection <- row$baseline_detection
    raw_results$delta_rate         <- row$delta_rate
    raw_results$analyte_rho        <- row$analyte_rho
    raw_results$rho_t              <- row$rho_t
    raw_results$cell_id            <- row$cell_id

    all_raw[[i]] <- raw_results
  }

  bind_rows(all_raw)
}


# ---- Aggregation Helpers ----------------------------------------------------

#' Safe mean of log-ratio (for MPOR bias)
safe_log_bias <- function(est, truth) {
  valid <- is.finite(est) & est > 0 & is.finite(truth) & truth > 0
  if (sum(valid) == 0) return(NA_real_)
  mean(log(est[valid]) - log(truth[valid]))
}


#' Aggregate raw results into cell-level summaries.
#'
#' Two-step aggregation:
#'   1. Per-replication summary (one row per DGP cell x method x rep)
#'   2. Across-replication means and MCSEs (one row per DGP cell x method)
aggregate_mcnemar <- function(all_results) {

  ok <- all_results %>%
    filter(!error, !is.na(method))

  # Step 1: Per-replication summaries
  rep_summary <- ok %>%
    group_by(n_subjects, baseline_detection, delta_rate, analyte_rho, rho_t,
             method, rep_id) %>%
    summarise(
      # Type I error: rejection rate among null analytes
      null_reject = mean((p_value < 0.05)[!is_signal], na.rm = TRUE),
      # Power: rejection rate among signal analytes
      sig_reject  = mean((p_value < 0.05)[is_signal], na.rm = TRUE),
      # CI coverage (null): true delta = 0
      null_cover  = mean(
        ((delta_ci_lo <= 0) & (delta_ci_hi >= 0))[!is_signal],
        na.rm = TRUE
      ),
      # CI coverage (signal): true delta = actual_delta_pct
      sig_cover   = mean(
        ((delta_ci_lo <= true_delta) & (delta_ci_hi >= true_delta))[is_signal],
        na.rm = TRUE
      ),
      # Delta bias (signal analytes, percentage points)
      sig_delta_bias = mean((delta - true_delta)[is_signal], na.rm = TRUE),
      # MPOR log-bias (signal)
      mpor_sig_bias = safe_log_bias(mpor[is_signal], true_mpor[is_signal]),
      # MPOR log-bias (null; true MPOR = 1, so log(true) = 0)
      mpor_null_bias = safe_log_bias(mpor[!is_signal], true_mpor[!is_signal]),
      # MPOR computability: fraction of discordant analytes with b > 0
      mpor_frac_defined = ifelse(
        sum(n_discordant > 0) > 0,
        mean(mpor_defined[n_discordant > 0], na.rm = TRUE),
        NA_real_
      ),
      # FDR components
      n_reject_q = sum(q_value < 0.05, na.rm = TRUE),
      n_false_q  = sum(q_value < 0.05 & !is_signal, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      fdr = ifelse(n_reject_q > 0, n_false_q / n_reject_q, 0)
    )

  # Step 2: Aggregate across replications
  cell_summary <- rep_summary %>%
    group_by(n_subjects, baseline_detection, delta_rate, analyte_rho, rho_t,
             method) %>%
    summarise(
      n_reps         = n(),
      type1          = mean(null_reject, na.rm = TRUE),
      type1_mcse     = sd(null_reject, na.rm = TRUE) / sqrt(n()),
      power          = mean(sig_reject, na.rm = TRUE),
      power_mcse     = sd(sig_reject, na.rm = TRUE) / sqrt(n()),
      cov_null       = mean(null_cover, na.rm = TRUE),
      cov_null_mcse  = sd(null_cover, na.rm = TRUE) / sqrt(n()),
      cov_sig        = mean(sig_cover, na.rm = TRUE),
      cov_sig_mcse   = sd(sig_cover, na.rm = TRUE) /
        sqrt(max(sum(is.finite(sig_cover)), 1)),
      delta_bias     = mean(sig_delta_bias, na.rm = TRUE),
      delta_bias_mcse = sd(sig_delta_bias, na.rm = TRUE) /
        sqrt(max(sum(is.finite(sig_delta_bias)), 1)),
      mpor_bias_sig  = mean(mpor_sig_bias, na.rm = TRUE),
      mpor_bias_sig_mcse = sd(mpor_sig_bias, na.rm = TRUE) /
        sqrt(max(sum(is.finite(mpor_sig_bias)), 1)),
      mpor_bias_null = mean(mpor_null_bias, na.rm = TRUE),
      mpor_bias_null_mcse = sd(mpor_null_bias, na.rm = TRUE) /
        sqrt(max(sum(is.finite(mpor_null_bias)), 1)),
      mpor_computable = mean(mpor_frac_defined, na.rm = TRUE),
      mean_fdr       = mean(fdr, na.rm = TRUE),
      fdr_mcse       = sd(fdr, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )

  list(rep_summary = rep_summary, cell_summary = cell_summary)
}


# ---- Execute ----------------------------------------------------------------

if (!is.null(S1_MODE)) {

# ---- S1: Flow-through validation of mcnemar_detection() wrapper -------------
#
# Review plan (`.claude/projects/.../plans/20260422_mcnemar_detection_review_plan.md`
# §6c, "Proposed targeted sims — S1"). B2 does not call the exported function
# mcnemar_detection(); it runs the kernel on raw 2x2 counts. S1 re-runs a
# small subset of B2 cells through both paths — standalone `exact/midp/noCC`
# analyte helpers AND the full wrapper — at matched seeds, so that per-rep
# parity and cell-aggregate drift can both be measured. If they agree within
# MCSE, B2's guarantees flow through to the exported API.
#
# Reduced S1 grid: 3 (n) x 1 (baseline) x 2 (delta) x 1 (rho) x 2 (rho_t) = 12 cells.

s1_grid <- make_dgp_grid(
  n_subjects         = c(30L, 50L, 100L),
  baseline_detection = 0.50,
  delta_rate         = c(0, 0.15),
  analyte_rho        = 0.3,
  rho_t              = c(0, 0.5)
)
s1_grid$cell_id <- s1_grid$cell_id + S1_OFFSET

if (VERBOSE) {
  message("=== Benchmark 2: S1 flow-through validation ===")
  message("DGP grid: ", nrow(s1_grid), " cells x ", N_REPS, " reps")
  message("Methods (standalone + wrapper): exact, mid-p, noCC")
  message("Cache prefix: mcnemar_S1_")
}

# Rep function that runs both paths at the same seed and rbinds their results.
# The wrapper DGP calls set.seed(seed) internally; the standalone DGP (via
# run_one_rep()) also calls set.seed(seed) first. Identical RNG state ->
# identical 2x2 counts per analyte, so per-rep p-value differences isolate
# the wrapper's preprocessing / backend differences, nothing else.
s1_kept_methods <- c("McNemar exact", "McNemar mid-p", "McNemar noCC")

make_s1_rep_fn <- function(ns, bd, dr, ar, rt) {
  force(ns); force(bd); force(dr); force(ar); force(rt)
  function(seed) {
    standalone <- run_one_rep(seed = seed, n_subjects = ns,
                              baseline_detection = bd, delta_rate = dr,
                              analyte_rho = ar, rho_t = rt)
    standalone <- standalone[standalone$method %in% s1_kept_methods, ,
                              drop = FALSE]
    standalone$path <- "standalone"

    wrapper <- run_one_rep_via_wrapper(seed = seed, n_subjects = ns,
                                       baseline_detection = bd, delta_rate = dr,
                                       analyte_rho = ar, rho_t = rt)
    wrapper$path <- "wrapper"

    # Align columns before rbind
    cols_union <- union(names(standalone), names(wrapper))
    for (col in setdiff(cols_union, names(standalone))) standalone[[col]] <- NA
    for (col in setdiff(cols_union, names(wrapper)))    wrapper[[col]]    <- NA
    rbind(standalone[, cols_union, drop = FALSE],
          wrapper[, cols_union, drop = FALSE])
  }
}

run_benchmark_s1 <- function(dgp_grid, n_reps, base_seed, cache_dir,
                              n_cores, verbose) {
  all_raw <- list()
  for (i in seq_len(nrow(dgp_grid))) {
    row <- dgp_grid[i, ]
    cell_label <- paste0("S1_cell_", row$cell_id,
                          " [n=", row$n_subjects,
                          "/delta=", row$delta_rate,
                          "/rho_t=", row$rho_t, "]")
    rep_fn <- make_s1_rep_fn(row$n_subjects, row$baseline_detection,
                             row$delta_rate, row$analyte_rho, row$rho_t)
    cache_file <- file.path(cache_dir,
                             paste0("mcnemar_S1_cell", row$cell_id, ".rds"))
    raw <- run_benchmark_cell(
      rep_fn     = rep_fn,
      n_reps     = n_reps,
      base_seed  = base_seed + (i - 1L) * 10000L,
      cache_file = cache_file,
      n_cores    = n_cores,
      cell_label = cell_label,
      verbose    = verbose
    )
    raw$n_subjects         <- row$n_subjects
    raw$baseline_detection <- row$baseline_detection
    raw$delta_rate         <- row$delta_rate
    raw$analyte_rho        <- row$analyte_rho
    raw$rho_t              <- row$rho_t
    raw$cell_id            <- row$cell_id
    all_raw[[i]] <- raw
  }
  bind_rows(all_raw)
}

if (VERBOSE) message("\nStarting S1 execution...")

all_results <- run_benchmark_s1(
  dgp_grid  = s1_grid,
  n_reps    = N_REPS,
  base_seed = BASE_SEED,
  cache_dir = CACHE_DIR,
  n_cores   = N_CORES,
  verbose   = VERBOSE
)

if (VERBOSE) message("\nAggregating S1 results...")

# --- Cell-level aggregation: per (path, method, cell) ---
ok <- all_results %>%
  filter(!error, !is.na(method), method %in% s1_kept_methods)

rep_summary_s1 <- ok %>%
  group_by(n_subjects, delta_rate, rho_t, path, method, rep_id) %>%
  summarise(
    null_reject = mean((p_value < 0.05)[!is_signal], na.rm = TRUE),
    sig_reject  = mean((p_value < 0.05)[is_signal],  na.rm = TRUE),
    .groups = "drop"
  )

cell_summary_s1 <- rep_summary_s1 %>%
  group_by(n_subjects, delta_rate, rho_t, path, method) %>%
  summarise(
    n_reps     = n(),
    type1      = mean(null_reject, na.rm = TRUE),
    type1_mcse = sd(null_reject, na.rm = TRUE)  / sqrt(n()),
    power      = mean(sig_reject, na.rm = TRUE),
    power_mcse = sd(sig_reject, na.rm = TRUE)   / sqrt(n()),
    .groups = "drop"
  )

s1_summary_file <- file.path(CACHE_DIR, "mcnemar_S1_summary.rds")
saveRDS(cell_summary_s1, s1_summary_file)
if (VERBOSE) message("S1 summary saved -> ", s1_summary_file)

# --- Per-rep parity diagnostics (wrapper vs standalone, matched seed) ---
parity <- ok %>%
  select(rep_id, n_subjects, delta_rate, rho_t, method, analyte, path,
         p_value, n_discordant) %>%
  tidyr::pivot_wider(
    names_from  = path,
    values_from = c(p_value, n_discordant),
    names_sep   = "_"
  ) %>%
  mutate(
    p_abs_diff    = abs(p_value_wrapper - p_value_standalone),
    discord_match = n_discordant_wrapper == n_discordant_standalone
  )

# na.rm removed from max()/sum() below: an NA in either path means a backend
# failed on that rep/analyte, which is exactly what we want to surface rather
# than silently drop. n_p_na and n_discord_na count the dropped rows explicitly.
parity_summary <- parity %>%
  group_by(method, n_subjects, delta_rate, rho_t) %>%
  summarise(
    n_rows             = n(),
    n_p_na             = sum(is.na(p_abs_diff)),
    n_discord_na       = sum(is.na(discord_match)),
    max_p_abs_diff     = if (all(is.na(p_abs_diff))) NA_real_
                         else max(p_abs_diff, na.rm = TRUE),
    med_p_abs_diff     = if (all(is.na(p_abs_diff))) NA_real_
                         else median(p_abs_diff, na.rm = TRUE),
    n_discord_mismatch = sum(!discord_match & !is.na(discord_match)),
    .groups = "drop"
  )

parity_file <- file.path(CACHE_DIR, "mcnemar_S1_parity.rds")
saveRDS(parity_summary, parity_file)
if (VERBOSE) {
  message("S1 parity diagnostics saved -> ", parity_file)
  message("\nPer-rep parity summary (max |p_wrap - p_stand|):")
  print(as.data.frame(parity_summary), right = FALSE)
}

# --- Cell-level wrapper-vs-standalone comparison (paired-diff MCSE) ---
#
# Wrapper and standalone share rep seeds by design: every (cell, method, rep)
# scores both paths on the same DGP draw. The correct MCSE for the drift band
# is therefore the per-rep *paired* difference, not the independent-samples
# sqrt(mcse_std^2 + mcse_wrap^2) formula. Under bit-for-bit path identity the
# per-rep diff vector is all zeros, sd = 0, and 2 * MCSE collapses to 0 -
# making pass_2mcse a true equality test that will flag *any* non-zero
# disagreement in a future run. The old independent-samples formula sat at
# ~0.008-0.012 regardless, so would not have caught a preprocessing drift
# that disagreed on roughly half of reps.

rep_paired_diff <- rep_summary_s1 %>%
  select(n_subjects, delta_rate, rho_t, method, rep_id, path,
         null_reject, sig_reject) %>%
  tidyr::pivot_wider(
    names_from  = path,
    values_from = c(null_reject, sig_reject),
    names_sep   = "_"
  ) %>%
  mutate(
    null_reject_diff = null_reject_wrapper - null_reject_standalone,
    sig_reject_diff  = sig_reject_wrapper  - sig_reject_standalone
  )

paired_mcse <- rep_paired_diff %>%
  group_by(n_subjects, delta_rate, rho_t, method) %>%
  summarise(
    n_reps                 = n(),
    type1_paired_mean_diff = mean(null_reject_diff, na.rm = TRUE),
    type1_paired_mcse      = sd(null_reject_diff,   na.rm = TRUE) /
      sqrt(sum(is.finite(null_reject_diff))),
    power_paired_mean_diff = mean(sig_reject_diff,  na.rm = TRUE),
    power_paired_mcse      = sd(sig_reject_diff,    na.rm = TRUE) /
      sqrt(max(sum(is.finite(sig_reject_diff)), 1)),
    .groups = "drop"
  )

type1_wide <- cell_summary_s1 %>%
  filter(delta_rate == 0) %>%
  select(method, n_subjects, rho_t, path, type1, type1_mcse) %>%
  tidyr::pivot_wider(
    names_from  = path,
    values_from = c(type1, type1_mcse),
    names_sep   = "_"
  ) %>%
  left_join(paired_mcse %>% filter(delta_rate == 0) %>%
              select(method, n_subjects, rho_t,
                     type1_paired_mcse, type1_paired_mean_diff),
            by = c("method", "n_subjects", "rho_t")) %>%
  mutate(
    abs_diff   = abs(type1_wrapper - type1_standalone),
    mcse_2     = 2 * type1_paired_mcse,
    pass_2mcse = abs_diff <= mcse_2
  )

type1_compare_file <- file.path(CACHE_DIR, "mcnemar_S1_type1_compare.rds")
saveRDS(type1_wide, type1_compare_file)
if (VERBOSE) {
  message("\nCell-level Type I comparison (standalone vs wrapper, paired MCSE):")
  print(as.data.frame(type1_wide), right = FALSE)
}

power_wide <- cell_summary_s1 %>%
  filter(delta_rate > 0) %>%
  select(method, n_subjects, rho_t, path, power, power_mcse) %>%
  tidyr::pivot_wider(
    names_from  = path,
    values_from = c(power, power_mcse),
    names_sep   = "_"
  ) %>%
  left_join(paired_mcse %>% filter(delta_rate > 0) %>%
              select(method, n_subjects, rho_t,
                     power_paired_mcse, power_paired_mean_diff),
            by = c("method", "n_subjects", "rho_t")) %>%
  mutate(
    abs_diff   = abs(power_wrapper - power_standalone),
    mcse_2     = 2 * power_paired_mcse,
    pass_2mcse = abs_diff <= mcse_2
  )

power_compare_file <- file.path(CACHE_DIR, "mcnemar_S1_power_compare.rds")
saveRDS(power_wide, power_compare_file)
if (VERBOSE) {
  message("\nCell-level Power comparison (standalone vs wrapper):")
  print(as.data.frame(power_wide), right = FALSE)
}

# --- Figure: Type I curves, wrapper vs standalone overlay ---
fig_dir <- file.path(script_dir, "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

type1_plot_data <- cell_summary_s1 %>%
  filter(delta_rate == 0) %>%
  mutate(rho_t_lbl = paste0("rho_t = ", rho_t))

p_s1 <- ggplot(type1_plot_data,
               aes(x = n_subjects, y = type1,
                   colour = method, linetype = path,
                   shape = path, group = interaction(method, path))) +
  geom_hline(yintercept = 0.05, colour = "grey40", linewidth = 0.3) +
  geom_errorbar(aes(ymin = type1 - 2 * type1_mcse,
                    ymax = type1 + 2 * type1_mcse),
                width = 3, alpha = 0.55) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2.2) +
  facet_wrap(~ rho_t_lbl, nrow = 1) +
  scale_colour_method() +
  scale_linetype_manual(values = c(standalone = "solid", wrapper = "22")) +
  scale_shape_manual(values = c(standalone = 16, wrapper = 1)) +
  scale_x_continuous(breaks = c(30, 50, 100)) +
  labs(
    title    = "B2 S1: Wrapper flow-through — Type I error",
    subtitle = "Standalone 2x2 kernel vs exported mcnemar_detection() wrapper (matched seed, 500 reps)",
    x = "Sample size (n_subjects)",
    y = "Empirical Type I error",
    colour = "Method", linetype = "Path", shape = "Path"
  ) +
  theme_benchmark()

fig_path <- file.path(fig_dir, "B2_S1_flowthrough.png")
ggsave(fig_path, p_s1, width = 8, height = 4, dpi = 140)
if (VERBOSE) message("\nS1 figure saved -> ", fig_path)

if (VERBOSE) message("\n=== Benchmark 2: S1 flow-through validation complete ===")
quit(save = "no", status = 0)

} else if (!is.null(PHASE_B)) {

# ---- Phase B4: Negative within-subject temporal correlation -----------------
#
# Tests McNemar calibration under rho_t = -0.3 (negative within-subject
# pairing), e.g. compensatory immune responses where high T1 detection
# predicts low T2 detection. The latent bivariate normal model handles
# negative rho_t correctly: Z2 = rho_t * Z1 + sqrt(1 - rho_t^2) * Z2_innov.
# With rho_t = -0.3, sqrt(1 - 0.09) = sqrt(0.91) is valid.

dgp_grid_b <- make_dgp_grid(
  n_subjects         = c(15L, 30L, 50L, 100L, 200L),
  baseline_detection = c(0.30, 0.50, 0.80),
  delta_rate         = c(0, 0.05, 0.15, 0.25),
  analyte_rho        = c(0, 0.3, 0.6),
  rho_t              = -0.3
)
dgp_grid_b$cell_id <- dgp_grid_b$cell_id + PHASE_B_OFFSET

if (VERBOSE) {
  message("=== Benchmark 2: Phase B4 — Negative temporal correlation ===")
  message("DGP grid: ", nrow(dgp_grid_b), " cells x ", N_REPS, " reps = ",
          nrow(dgp_grid_b) * N_REPS, " total runs")
  message("Analytes per sim: ", N_ANALYTES, " (", N_SIGNAL, " signal)")
  message("Temporal correlation grid: -0.3 (Phase B4)")
  message("Methods: 5 (exact, CC, chi-sq, mid-p, noCC)")
  message("Cache dir: ", CACHE_DIR)
}

if (VERBOSE) message("\nStarting Phase B4 execution...")

# Swap BENCHMARK for Phase B cache filenames (mcnemar_v2_phaseB_cell*.rds)
BENCHMARK_MAIN <- BENCHMARK
BENCHMARK <- paste0(BENCHMARK_MAIN, "_phaseB")

all_results <- run_benchmark_mcnemar(
  dgp_grid  = dgp_grid_b,
  n_reps    = N_REPS,
  base_seed = BASE_SEED,
  cache_dir = CACHE_DIR,
  n_cores   = N_CORES,
  verbose   = VERBOSE
)

# Restore original BENCHMARK name
BENCHMARK <- BENCHMARK_MAIN

if (VERBOSE) message("\nAggregating Phase B4 results...")

agg <- aggregate_mcnemar(all_results)
cell_summary <- agg$cell_summary

# Acceptance annotations
cell_annotated <- cell_summary %>%
  mutate(
    pass_type1 = ifelse(
      delta_rate == 0,
      type1 >= 0.025 & type1 <= 0.075,
      NA
    ),
    pass_cov_null = cov_null >= 0.93 & cov_null <= 0.97,
    pass_cov_sig  = ifelse(
      delta_rate > 0,
      cov_sig >= 0.93 & cov_sig <= 0.97,
      NA
    ),
    pass_fdr = mean_fdr <= 0.06  # slight tolerance for Monte Carlo noise
  )

# Save and merge Phase B4 into existing summary
summary_file <- file.path(CACHE_DIR, paste0(BENCHMARK, "_summary.rds"))
merge_with_existing_summary(cell_annotated, summary_file, verbose = VERBOSE)
if (VERBOSE) message("Phase B4 summary merged -> ", summary_file)

if (VERBOSE) message("\nSkipping figures in --phase_b mode (partial data)")

} else {
# ---- Main grid execution ----------------------------------------------------

if (VERBOSE) message("\nStarting benchmark execution...")

all_results <- run_benchmark_mcnemar(
  dgp_grid  = dgp_grid,
  n_reps    = N_REPS,
  base_seed = BASE_SEED,
  cache_dir = CACHE_DIR,
  n_cores   = N_CORES,
  verbose   = VERBOSE
)

if (VERBOSE) message("\nAggregating results...")

agg <- aggregate_mcnemar(all_results)
cell_summary <- agg$cell_summary

# Acceptance annotations
cell_annotated <- cell_summary %>%
  mutate(
    pass_type1 = ifelse(
      delta_rate == 0,
      type1 >= 0.025 & type1 <= 0.075,
      NA
    ),
    pass_cov_null = cov_null >= 0.93 & cov_null <= 0.97,
    pass_cov_sig  = ifelse(
      delta_rate > 0,
      cov_sig >= 0.93 & cov_sig <= 0.97,
      NA
    ),
    pass_fdr = mean_fdr <= 0.06  # slight tolerance for Monte Carlo noise
  )

# Save final summary (merge with existing when running extra sizes)
summary_file <- file.path(CACHE_DIR, paste0(BENCHMARK, "_summary.rds"))
if (!is.null(EXTRA_SIZES)) {
  merge_with_existing_summary(cell_annotated, summary_file, verbose = VERBOSE)
} else {
  save_main_summary_preserving_supplements(cell_annotated, summary_file,
                                           verbose = VERBOSE)
}


# ---- Size-Adjusted Power + Figures ------------------------------------------

if (!is.null(EXTRA_SIZES)) {
  if (VERBOSE) message("\nSkipping figures in --extra_n mode (partial data)")
} else {

if (VERBOSE) message("\nComputing size-adjusted power...")

# Step 1: Calibrated thresholds from null (delta = 0) data.
# Under delta = 0, all analytes are truly null regardless of signal_mask.
null_pvals <- all_results %>%
  filter(!error, !is.na(method), delta_rate == 0)

alpha_star <- null_pvals %>%
  group_by(n_subjects, baseline_detection, analyte_rho, rho_t, method) %>%
  summarise(
    threshold = quantile(p_value, probs = 0.05, na.rm = TRUE),
    .groups = "drop"
  )

# Step 2: Apply calibrated thresholds to signal analytes under alternative
sig_alt <- all_results %>%
  filter(!error, !is.na(method), delta_rate > 0, is_signal) %>%
  inner_join(alpha_star, by = c("n_subjects", "baseline_detection",
                                 "analyte_rho", "rho_t", "method"))

size_adj_power <- sig_alt %>%
  group_by(n_subjects, baseline_detection, delta_rate, analyte_rho, rho_t,
           method) %>%
  summarise(
    adj_power = mean(p_value < threshold, na.rm = TRUE),
    adj_power_mcse = sqrt(adj_power * (1 - adj_power) / n()),
    .groups = "drop"
  )

# Save size-adjusted power
adj_power_file <- file.path(CACHE_DIR, paste0(BENCHMARK, "_adj_power.rds"))
saveRDS(size_adj_power, adj_power_file)
if (VERBOSE) message("Size-adjusted power saved -> ", adj_power_file)



# Figures are generated by the standalone plot_B2.R script.
# Run: Rscript inst/simulations/plot_B2.R [--cache_dir PATH]

} # end if (!is.null(EXTRA_SIZES)) — skip figures

} # end if (!is.null(PHASE_B)) — main vs Phase B

# ---- Formatted Summary Table ------------------------------------------------

summary_formatted <- cell_annotated %>%
  mutate(
    type1_fmt     = fmt_mcse_vec(type1, type1_mcse),
    power_fmt     = fmt_mcse_vec(power, power_mcse),
    cov_null_fmt  = fmt_mcse_vec(cov_null, cov_null_mcse),
    cov_sig_fmt   = fmt_mcse_vec(cov_sig, cov_sig_mcse),
    fdr_fmt       = fmt_mcse_vec(mean_fdr, fdr_mcse),
    mpor_bias_fmt = fmt_mcse_vec(mpor_bias_sig, mpor_bias_sig_mcse),
    mpor_comp_pct = round(mpor_computable * 100, 1)
  )

# Print concise table (subset: rho=0, baseline=0.50, rho_t=0.5 for brevity)
if (VERBOSE) {
  message("\n=== Benchmark 2: Summary (baseline=50%, rho=0, rho_t=0.5) ===\n")
  print(
    summary_formatted %>%
      filter(baseline_detection == 0.50, analyte_rho == 0, rho_t == 0.5) %>%
      select(n_subjects, delta_rate, method,
             type1_fmt, power_fmt, cov_null_fmt, cov_sig_fmt,
             fdr_fmt, mpor_comp_pct, starts_with("pass_")) %>%
      as.data.frame(),
    right = FALSE
  )
}

if (VERBOSE) {
  message("\n=== Benchmark 2 (v2) complete ===")
  if (exists("fig_dir")) message("Figures -> ", fig_dir)
  message("Summary -> ", summary_file)
}
