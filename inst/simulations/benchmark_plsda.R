#!/usr/bin/env Rscript
# benchmark_plsda.R — Benchmark 4: PLS-DA Discrimination Recovery
#
# ADEMP Structure:
#   Aims:    Evaluate whether PLS-DA correctly identifies discriminatory analytes
#            and achieves valid classification, with proper cross-validation and
#            permutation testing.
#   DGP:     simulate_plsda_data() with factorial grid over n_per_group,
#            n_discriminatory, effect_size, correlation, censoring, class_imbalance.
#   Estimand: Classification of group membership; identification of truly
#             discriminatory analytes.
#   Methods: (1) PLS-DA nested 5-fold CV, (2) PLS-DA single 7-fold CV,
#            (3) Permutation null test (via ropls built-in permI).
#   Performance: balanced accuracy (nested CV), VIP AUROC, Q2, feature stability
#                (Jaccard), Q2 optimism gap, permutation FPR under null.
#
# Usage:
#   Rscript benchmark_plsda.R [--n_reps N] [--n_cores N] [--cache_dir DIR]
#                              [--n_perm N] [--phase_b B5] [--dev]
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

N_REPS     <- as.integer(parse_cli_arg("--n_reps",  "200"))
N_CORES    <- as.integer(parse_cli_arg("--n_cores", "1"))
CACHE_DIR  <- parse_cli_arg("--cache_dir", file.path(script_dir, "cache"))
EXTRA_N    <- parse_cli_arg("--extra_n", NULL)
PHASE_B    <- parse_cli_arg("--phase_b", NULL)
PHASE_C    <- parse_cli_arg("--phase_c", NULL)
DEV_MODE   <- "--dev" %in% cli_args
BASE_SEED  <- 20240004L
BENCHMARK  <- "plsda_v3"

# Phase C defaults: c3 = 1000 permutations; c4 = alternative transforms
if (!is.null(PHASE_C) && identical(PHASE_C, "c3")) {
  # Default N_PERM = 1000 in c3 mode unless user overrode it
  N_PERM <- as.integer(parse_cli_arg("--n_perm", "1000"))
} else {
  N_PERM <- as.integer(parse_cli_arg("--n_perm",  "100"))
}

# Extra-n supplement mode
EXTRA_SIZES  <- if (!is.null(EXTRA_N)) as.integer(strsplit(EXTRA_N, ",")[[1]]) else NULL
EXTRA_OFFSET <- 10000L
PHASE_B_OFFSET <- 20000L
PHASE_C_OFFSET <- 30000L
N_ANALYTES <- 20L
MAX_COMP   <- 3L   # Max components tested in nested CV inner loop
K_OUTER    <- 5L   # Outer CV folds for nested CV
K_INNER    <- 5L   # Inner CV folds for component selection
VERBOSE    <- TRUE

if (!dir.exists(CACHE_DIR)) dir.create(CACHE_DIR, recursive = TRUE)


# ---- PLS-DA Helper Functions ------------------------------------------------

#' Fit ropls::opls silently (suppress all console and file output)
#' @param ... Arguments forwarded to ropls::opls.
#' @return The fitted opls model object.
fit_ropls_silent <- function(...) {
  model <- NULL
  invisible(capture.output(
    model <- suppressMessages(suppressWarnings(
      ropls::opls(...)
    )),
    type = "output"
  ))
  model
}

# Predict from an opls model (works around S3/S4 dispatch issue where
# UseMethod("predict") fails to find the S4 method registered by ropls).
# Uses manual coefficient multiplication as a robust fallback.
predict_opls <- function(model, newdata) {
  # Try S4 dispatch first
  pred <- tryCatch(
    selectMethod("predict", "opls")(model, newdata),
    error = function(e) NULL
  )
  if (!is.null(pred)) return(as.character(pred))

  # Manual fallback: coefficients × scaled newdata + y intercept
  pred_y <- scale(newdata, center = model@xMeanVn, scale = model@xSdVn) %*%
    model@coefficientMN
  pred_y <- sweep(pred_y, 2, model@yMeanVn, "+")

  resp_levels <- levels(model@suppLs$y)
  if (length(resp_levels) >= 2) {
    # Binary DA: threshold at 0.5 (y coded as 0/1)
    ifelse(pred_y[, 1] >= 0.5, resp_levels[2], resp_levels[1])
  } else {
    rep(NA_character_, nrow(newdata))
  }
}


#' Area under the ROC curve via Wilcoxon-Mann-Whitney statistic
#'
#' @param scores Numeric vector (higher = more likely positive).
#' @param labels Logical vector (TRUE = positive class).
#' @return Scalar AUROC in [0, 1], or NA if degenerate.
calc_auroc <- function(scores, labels) {
  ok <- is.finite(scores) & !is.na(labels)
  scores <- scores[ok]; labels <- labels[ok]
  n_pos <- sum(labels)
  n_neg <- sum(!labels)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  ranks <- rank(scores)
  (sum(ranks[labels]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}


#' Mean pairwise Jaccard index of feature sets
#'
#' @param set_list List of character vectors (one per CV fold).
#' @return Scalar mean Jaccard index, or NA if < 2 valid sets.
calc_jaccard <- function(set_list) {
  set_list <- set_list[vapply(set_list, function(s) length(s) > 0, logical(1))]
  n <- length(set_list)
  if (n < 2) return(NA_real_)
  jaccards <- numeric(choose(n, 2))
  idx <- 1L
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      isect <- length(intersect(set_list[[i]], set_list[[j]]))
      union_n <- length(union(set_list[[i]], set_list[[j]]))
      jaccards[idx] <- if (union_n > 0) isect / union_n else 0
      idx <- idx + 1L
    }
  }
  mean(jaccards)
}


#' Create stratified k-fold assignments
#'
#' @param labels Factor or character vector of class labels.
#' @param k      Integer number of folds.
#' @return Integer vector of fold assignments (values in 1:k).
create_stratified_folds <- function(labels, k) {
  n <- length(labels)
  folds <- integer(n)
  for (level in unique(labels)) {
    idx <- which(labels == level)
    folds[idx] <- sample(rep_len(seq_len(k), length(idx)))
  }
  folds
}


#' Balanced accuracy: mean per-class accuracy
calc_balanced_accuracy <- function(predicted, actual) {
  valid <- !is.na(predicted) & !is.na(actual)
  if (sum(valid) == 0) return(NA_real_)
  pred <- predicted[valid]; true <- actual[valid]
  per_class <- tapply(pred == true, true, mean, na.rm = TRUE)
  mean(per_class, na.rm = TRUE)
}


#' Custom aggregation for PLS-DA benchmark cells
#'
#' Computes mean +/- MCSE for each PLS-DA metric, split by method.
#'
#' @param raw_results data.frame from run_benchmark_cell.
#' @return tibble with one row per method.
aggregate_plsda_cell <- function(raw_results) {
  ok <- raw_results[!raw_results$error, , drop = FALSE]
  if (nrow(ok) == 0) return(NULL)

  summ <- function(x) {
    x <- x[is.finite(x)]
    n <- length(x)
    if (n == 0) return(c(mean = NA_real_, mcse = NA_real_))
    c(mean = mean(x), mcse = sd(x) / sqrt(n))
  }

  methods <- unique(ok$method)
  out <- list()

  for (m in methods) {
    md <- ok[ok$method == m, , drop = FALSE]
    acc  <- summ(md$balanced_accuracy)
    vip  <- summ(md$vip_auroc)
    q2v  <- summ(md$q2)
    stab <- summ(md$feature_stability)
    gap  <- summ(md$q2_optimism_gap)

    # Permutation rejection rate (power or FPR)
    pvals <- md$perm_pvalue[is.finite(md$perm_pvalue)]
    if (length(pvals) > 0) {
      rej <- mean(pvals < 0.05)
      rej_mcse <- sqrt(rej * (1 - rej) / length(pvals))
    } else {
      rej <- NA_real_; rej_mcse <- NA_real_
    }

    conv <- calc_convergence_rate(md$converged)

    out[[length(out) + 1]] <- tibble(
      method                 = m,
      n_reps                 = nrow(md),
      convergence_rate       = conv$estimate,
      convergence_mcse       = conv$mcse,
      balanced_accuracy      = unname(acc["mean"]),
      balanced_accuracy_mcse = unname(acc["mcse"]),
      vip_auroc              = unname(vip["mean"]),
      vip_auroc_mcse         = unname(vip["mcse"]),
      q2                     = unname(q2v["mean"]),
      q2_mcse                = unname(q2v["mcse"]),
      feature_stability      = unname(stab["mean"]),
      feature_stability_mcse = unname(stab["mcse"]),
      q2_optimism_gap        = unname(gap["mean"]),
      q2_optimism_gap_mcse   = unname(gap["mcse"]),
      perm_rejection_rate    = rej,
      perm_rejection_mcse    = rej_mcse
    )
  }
  bind_rows(out)
}


# ---- DGP Grid ---------------------------------------------------------------

if (DEV_MODE) {
  # Minimal grid for pipeline validation
  dgp_grid <- make_dgp_grid(
    n_per_group      = if (!is.null(EXTRA_SIZES)) EXTRA_SIZES else c(25L, 50L),
    n_discriminatory = 5L,
    effect_size      = c(0, 0.7, 1.2),
    correlation      = c("none", "block_0.5"),
    censoring        = 0.01,
    imbalance        = 1L
  )
  if (VERBOSE) message("DEV MODE: reduced grid (", nrow(dgp_grid), " cells)")
} else {
  dgp_grid <- make_dgp_grid(
    n_per_group      = if (!is.null(EXTRA_SIZES)) EXTRA_SIZES else c(20L, 50L, 100L),
    n_discriminatory = c(3L, 5L, 10L),
    effect_size      = c(0, 0.3, 0.7, 1.2),
    correlation      = c("none", "block_0.5", "block_0.8"),
    censoring        = c(0.01, 0.20),
    imbalance        = c(1L, 2L, 3L)
  )
}
if (!is.null(EXTRA_SIZES)) dgp_grid$cell_id <- dgp_grid$cell_id + EXTRA_OFFSET

if (VERBOSE) {
  message("=== Benchmark 4: PLS-DA Discrimination Recovery ===")
  message("DGP grid: ", nrow(dgp_grid), " cells x ", N_REPS, " reps = ",
          nrow(dgp_grid) * N_REPS, " total runs")
  message("Permutations per single-CV fit: ", N_PERM)
  message("Cache dir: ", CACHE_DIR)
}


# ---- Per-Replication Function ------------------------------------------------

#' Run one replication for a single DGP cell
#'
#' Simulates data, preprocesses for PLS-DA, runs nested CV (Method 1) and
#' single CV with permutation test (Method 2), returns per-method metrics.
#'
#' @param seed            Integer RNG seed.
#' @param n_per_group     Samples per group (minority group for imbalanced).
#' @param n_discriminatory Number of truly discriminatory analytes.
#' @param effect_size     Log-scale group effect (0 = null).
#' @param correlation     "none", "block_0.5", or "block_0.8".
#' @param censoring       LOD quantile (0.01 ~ 0%, 0.20 = 20%).
#' @param imbalance       Class imbalance ratio (1 = balanced, 2 = 2:1, 3 = 3:1).
#' @param n_perm          Number of permutations for single-CV model.
#' @return data.frame with one row per method.
run_one_rep <- function(seed, n_per_group, n_discriminatory, effect_size,
                        correlation, censoring, imbalance, n_perm,
                        n_analytes = N_ANALYTES,
                        transform = NULL) {

  # Failure result returned on any unrecoverable error
  fail_result <- data.frame(
    method            = c("PLS-DA nested CV", "PLS-DA single CV"),
    balanced_accuracy = NA_real_,
    vip_auroc         = NA_real_,
    q2                = NA_real_,
    feature_stability = NA_real_,
    perm_pvalue       = NA_real_,
    q2_optimism_gap   = NA_real_,
    converged         = FALSE,
    stringsAsFactors  = FALSE
  )

  tryCatch({

    # --- Parse DGP parameters ---
    cor_type  <- if (correlation == "none") NULL else "block"
    block_rho <- if (correlation == "none") 0 else as.numeric(sub("block_", "", correlation))

    # Total n: balanced → 2*n, 2:1 → 3*n, 3:1 → 4*n
    total_n     <- n_per_group * (1L + as.integer(imbalance))
    class_imbal <- if (imbalance == 1L) NULL else as.numeric(imbalance)

    # --- Simulate ---
    sim <- simulate_plsda_data(
      n_subjects       = total_n,
      n_analytes       = n_analytes,
      n_discriminatory = n_discriminatory,
      group_effects    = effect_size,
      analyte_correlation = cor_type,
      block_rho        = block_rho,
      class_imbalance  = class_imbal,
      lod_quantile     = censoring,
      seed             = seed
    )

    signal_names  <- sim$meta$dgp_params$signal_analyte_names
    analyte_names <- sim$meta$dgp_params$analyte_names

    # --- Prepare wide format ---
    d <- sim$data[!is.na(sim$data$value_raw), ]

    wide <- tidyr::pivot_wider(
      d[, c("subject_id", "group", "cytokine", "value_raw")],
      names_from  = "cytokine",
      values_from = "value_raw",
      values_fn   = mean      # handle any duplicate rows gracefully
    )

    lod_lookup <- data.frame(
      cytokine = analyte_names,
      lod = vapply(analyte_names, function(nm) {
        vals <- d$lod[d$cytokine == nm]
        vals[!is.na(vals)][1]
      }, numeric(1)),
      stringsAsFactors = FALSE
    )

    # --- Preprocess (LOD sub + optional transform/scaling) ---
    pp <- plsda_preprocess(
      data          = as.data.frame(wide),
      cytokine_cols = analyte_names,
      metadata_cols = c("subject_id", "group"),
      lod_lookup    = lod_lookup,
      lod_method    = "half",
      log_transform = TRUE,
      scale_data    = TRUE,
      transform     = transform,
      verbose       = FALSE
    )

    expr_mat <- pp$expression
    response <- as.factor(pp$metadata$group)
    n_total  <- nrow(expr_mat)
    analyte_names_pp <- colnames(expr_mat)
    is_signal_pp     <- analyte_names_pp %in% signal_names
    resp_levels      <- levels(response)
    has_signal       <- any(is_signal_pp)

    # ===================================================================
    # METHOD 1: Nested (double) cross-validation
    # ===================================================================
    outer_folds    <- create_stratified_folds(response, K_OUTER)
    outer_preds    <- rep(NA_character_, n_total)
    outer_vip_sets <- vector("list", K_OUTER)
    outer_q2s      <- rep(NA_real_, K_OUTER)
    vip_matrix     <- matrix(NA_real_, nrow = K_OUTER, ncol = length(analyte_names_pp),
                             dimnames = list(NULL, analyte_names_pp))

    for (k in seq_len(K_OUTER)) {
      test_idx  <- which(outer_folds == k)
      train_idx <- which(outer_folds != k)
      train_expr <- expr_mat[train_idx, , drop = FALSE]
      train_resp <- response[train_idx]
      test_expr  <- expr_mat[test_idx, , drop = FALSE]

      # Guard: minimum samples per group in training
      min_train_grp <- min(table(train_resp))
      if (min_train_grp < 3) next

      # Inner CV: select best n_components by Q2
      max_comp_k <- min(MAX_COMP, n_discriminatory, floor(min_train_grp / 2))
      max_comp_k <- max(1L, max_comp_k)
      inner_cv_k <- max(2L, min(K_INNER, min_train_grp))

      best_q2    <- -Inf
      best_ncomp <- 1L

      for (nc in seq_len(max_comp_k)) {
        inner_model <- tryCatch(
          fit_ropls_silent(
            x = train_expr, y = train_resp,
            predI = nc, orthoI = 0,
            crossvalI = inner_cv_k,
            permI = 0, scaleC = "none"
          ),
          error = function(e) NULL
        )
        if (!is.null(inner_model)) {
          q2 <- tryCatch(inner_model@summaryDF$`Q2(cum)`[1], error = function(e) NA_real_)
          if (!is.na(q2) && q2 > best_q2) {
            best_q2    <- q2
            best_ncomp <- nc
          }
        }
      }
      outer_q2s[k] <- best_q2

      # Fit on full training with best ncomp
      outer_model <- tryCatch(
        fit_ropls_silent(
          x = train_expr, y = train_resp,
          predI = best_ncomp, orthoI = 0,
          crossvalI = inner_cv_k,
          permI = 0, scaleC = "none"
        ),
        error = function(e) NULL
      )

      if (!is.null(outer_model)) {
        # Out-of-sample predictions
        preds <- tryCatch(
          predict_opls(outer_model, test_expr),
          error = function(e) rep(NA_character_, length(test_idx))
        )
        outer_preds[test_idx] <- preds

        # VIP scores for this fold
        vip <- tryCatch(ropls::getVipVn(outer_model), error = function(e) NULL)
        if (!is.null(vip)) {
          vip_matrix[k, names(vip)] <- vip
          top_k_n <- min(n_discriminatory, length(vip))
          outer_vip_sets[[k]] <- names(sort(vip, decreasing = TRUE))[seq_len(top_k_n)]
        }
      }
    } # end outer CV

    # --- Nested CV metrics ---
    nested_balanced_acc <- calc_balanced_accuracy(outer_preds, as.character(response))

    mean_vip <- colMeans(vip_matrix, na.rm = TRUE)
    nested_vip_auroc <- if (has_signal) calc_auroc(mean_vip, is_signal_pp) else NA_real_

    nested_jaccard <- calc_jaccard(outer_vip_sets)
    nested_q2 <- mean(outer_q2s[is.finite(outer_q2s)], na.rm = TRUE)
    if (is.nan(nested_q2)) nested_q2 <- NA_real_

    # ===================================================================
    # METHOD 2: Single 7-fold CV with permutation test
    # ===================================================================
    single_cv_k  <- min(7L, min(table(response)))
    single_cv_k  <- max(2L, single_cv_k)
    single_ncomp <- min(2L, MAX_COMP, n_discriminatory)

    single_model <- tryCatch(
      fit_ropls_silent(
        x = expr_mat, y = response,
        predI = single_ncomp, orthoI = 0,
        crossvalI = single_cv_k,
        permI = n_perm, scaleC = "none"
      ),
      error = function(e) NULL
    )

    if (!is.null(single_model)) {
      single_q2 <- tryCatch(single_model@summaryDF$`Q2(cum)`[1], error = function(e) NA_real_)
      sv <- tryCatch(ropls::getVipVn(single_model), error = function(e) NULL)
      single_vip_auroc <- if (!is.null(sv) && has_signal) {
        calc_auroc(sv[analyte_names_pp], is_signal_pp)
      } else {
        NA_real_
      }
      perm_pval <- tryCatch(single_model@summaryDF$pQ2[1], error = function(e) NA_real_)
      single_converged <- TRUE
    } else {
      single_q2 <- NA_real_
      single_vip_auroc <- NA_real_
      perm_pval <- NA_real_
      single_converged <- FALSE
    }

    # Q2 optimism gap
    q2_gap <- if (is.finite(single_q2) && is.finite(nested_q2)) {
      single_q2 - nested_q2
    } else {
      NA_real_
    }

    # --- Return combined results ---
    data.frame(
      method            = c("PLS-DA nested CV", "PLS-DA single CV"),
      balanced_accuracy = c(nested_balanced_acc, NA_real_),
      vip_auroc         = c(nested_vip_auroc, single_vip_auroc),
      q2                = c(nested_q2, single_q2),
      feature_stability = c(nested_jaccard, NA_real_),
      perm_pvalue       = c(perm_pval, perm_pval),
      q2_optimism_gap   = c(q2_gap, q2_gap),
      converged         = c(!is.na(nested_balanced_acc), single_converged),
      stringsAsFactors  = FALSE
    )

  }, error = function(e) {
    warning("Replication failed (seed=", seed, "): ", conditionMessage(e))
    fail_result
  })
}


# ---- Run Benchmark -----------------------------------------------------------

run_benchmark_plsda <- function(dgp_grid, n_reps, base_seed,
                                cache_dir, n_cores, n_perm, verbose,
                                n_analytes = N_ANALYTES,
                                cache_prefix = BENCHMARK,
                                transform = NULL) {

  all_summaries <- list()

  for (i in seq_len(nrow(dgp_grid))) {
    row <- dgp_grid[i, ]
    cell_label <- paste0(
      "cell_", row$cell_id,
      " [n=", row$n_per_group,
      "/disc=", row$n_discriminatory,
      "/eff=", row$effect_size,
      "/cor=", row$correlation,
      "/cens=", row$censoring,
      "/imb=", row$imbalance, "]"
    )

    # Closure over DGP parameters
    cell_transform <- if ("transform" %in% names(row)) row$transform else transform
    rep_fn <- local({
      rr <- row; np <- n_perm; na <- n_analytes; tr <- cell_transform
      function(seed) {
        run_one_rep(
          seed             = seed,
          n_per_group      = rr$n_per_group,
          n_discriminatory = rr$n_discriminatory,
          effect_size      = rr$effect_size,
          correlation      = rr$correlation,
          censoring        = rr$censoring,
          imbalance        = rr$imbalance,
          n_perm           = np,
          n_analytes       = na,
          transform        = tr
        )
      }
    })

    cache_file <- file.path(cache_dir, paste0(cache_prefix, "_cell", row$cell_id, ".rds"))

    # Prefer explicit seed_index when provided (used by Phase C to pair seeds
    # with the main grid / across transforms); fall back to loop position.
    seed_idx <- if ("seed_index" %in% names(row)) row$seed_index else i

    raw_results <- run_benchmark_cell(
      rep_fn     = rep_fn,
      n_reps     = n_reps,
      base_seed  = base_seed + (seed_idx - 1) * 10000L,
      cache_file = cache_file,
      n_cores    = n_cores,
      cell_label = cell_label,
      verbose    = verbose
    )

    # Aggregate per method
    agg <- aggregate_plsda_cell(raw_results)
    if (!is.null(agg)) {
      agg$n_per_group      <- row$n_per_group
      agg$n_discriminatory <- row$n_discriminatory
      agg$effect_size      <- row$effect_size
      agg$correlation      <- row$correlation
      agg$censoring        <- row$censoring
      agg$imbalance        <- row$imbalance
      agg$cell_id          <- row$cell_id
      agg$n_analytes       <- n_analytes
      if (!is.null(cell_transform)) agg$transform <- cell_transform
      all_summaries[[length(all_summaries) + 1]] <- agg
    }
  }

  bind_rows(all_summaries)
}


# ---- Execute -----------------------------------------------------------------

#' Annotate summary with acceptance criteria
annotate_acceptance <- function(df) {
  df %>%
    mutate(
      # VIP AUROC > 0.80 when effect >= 0.7 and n >= 20
      pass_vip_auroc = ifelse(
        effect_size >= 0.7 & n_per_group >= 20 & method == "PLS-DA nested CV",
        vip_auroc > 0.80,
        NA
      ),
      # Feature stability Jaccard > 0.6 when effect >= 0.7 and n >= 20
      pass_stability = ifelse(
        effect_size >= 0.7 & n_per_group >= 20 & method == "PLS-DA nested CV",
        feature_stability > 0.6,
        NA
      ),
      # Permutation FPR <= 0.075 under null (allow slight slack for MCSE)
      pass_perm_fpr = ifelse(
        effect_size == 0,
        perm_rejection_rate <= 0.075,
        NA
      )
    )
}

if (!is.null(PHASE_C)) {

  # ---- Phase C: parameter tweaks (c3 = 1000 permutations; c4 = transforms) ---

  if (identical(PHASE_C, "c3")) {

    if (VERBOSE) {
      message("\n=== Phase C3: 1000-permutation subset ===")
      message("N_PERM = ", N_PERM)
    }

    # Subset: filter the MAIN dgp_grid to effect_size in {0, 0.7} and carry
    # the original row index as seed_index so C3 replicates the main-grid's
    # exact seed stream. This lets C3 (1000 perms) be compared to main-grid
    # (100 perms) on identical replicates — otherwise differences conflate
    # permutation count with Monte-Carlo draw variance.
    c3_idx  <- which(dgp_grid$effect_size %in% c(0, 0.7))
    c3_grid <- dgp_grid[c3_idx, , drop = FALSE]
    c3_grid$seed_index <- c3_idx
    c3_grid$cell_id    <- c3_grid$cell_id + PHASE_C_OFFSET

    if (VERBOSE) {
      message("C3 grid: ", nrow(c3_grid), " cells x ", N_REPS, " reps x ",
              N_PERM, " perms")
    }

    c3_prefix <- "plsda_v4_c3"
    c3_summary <- run_benchmark_plsda(
      dgp_grid     = c3_grid,
      n_reps       = N_REPS,
      base_seed    = BASE_SEED,
      cache_dir    = CACHE_DIR,
      n_cores      = N_CORES,
      n_perm       = N_PERM,
      verbose      = VERBOSE,
      n_analytes   = N_ANALYTES,
      cache_prefix = c3_prefix
    )

    c3_annotated <- annotate_acceptance(c3_summary)
    c3_file <- file.path(CACHE_DIR, paste0(c3_prefix, "_summary.rds"))
    saveRDS(c3_annotated, c3_file)
    if (VERBOSE) message("Phase C3 summary saved -> ", c3_file)

    quit(save = "no", status = 0)
  }

  if (identical(PHASE_C, "c4")) {

    if (VERBOSE) message("\n=== Phase C4: alternative PLS-DA transforms ===")

    transforms <- c("log2p1_zscore", "cuberoot_zscore", "pareto")

    # Reduced grid to keep runtime manageable (drop block_0.8, imbalance=1 only)
    c4_base <- make_dgp_grid(
      n_per_group      = if (!is.null(EXTRA_SIZES)) EXTRA_SIZES else c(20L, 50L, 100L),
      n_discriminatory = c(3L, 5L, 10L),
      effect_size      = c(0, 0.7),
      correlation      = c("none", "block_0.5"),
      censoring        = c(0.01, 0.20),
      imbalance        = 1L
    )

    # Cross transforms with the base grid. Each transform subset shares the
    # same seed_index per DGP-factor tuple (base cell_id), so differences
    # between transforms are attributable to the transform rather than to
    # independent Monte-Carlo draws.
    c4_base$seed_index <- c4_base$cell_id
    c4_grid_list <- lapply(transforms, function(tr) {
      g <- c4_base
      g$transform <- tr
      g
    })
    c4_grid <- bind_rows(c4_grid_list)
    c4_grid$cell_id <- seq_len(nrow(c4_grid)) + PHASE_C_OFFSET + 5000L

    if (VERBOSE) {
      message("C4 grid: ", nrow(c4_grid), " cells (",
              nrow(c4_base), " base x ", length(transforms), " transforms)")
    }

    c4_prefix <- "plsda_v4_c4"
    c4_summary <- run_benchmark_plsda(
      dgp_grid     = c4_grid,
      n_reps       = N_REPS,
      base_seed    = BASE_SEED,
      cache_dir    = CACHE_DIR,
      n_cores      = N_CORES,
      n_perm       = N_PERM,
      verbose      = VERBOSE,
      n_analytes   = N_ANALYTES,
      cache_prefix = c4_prefix
    )

    c4_annotated <- annotate_acceptance(c4_summary)
    c4_file <- file.path(CACHE_DIR, paste0(c4_prefix, "_summary.rds"))
    saveRDS(c4_annotated, c4_file)
    if (VERBOSE) message("Phase C4 summary saved -> ", c4_file)

    quit(save = "no", status = 0)
  }

  stop("Unknown PHASE_C mode: ", PHASE_C, " (expected 'c3' or 'c4')")
}

if (!is.null(PHASE_B)) {
  # ---- Phase B: n_analytes sub-grids ------------------------------------------
  phase_b_analytes <- c(10L, 40L)
  phase_b_offsets  <- setNames(
    PHASE_B_OFFSET + (seq_along(phase_b_analytes) - 1L) * 5000L,
    as.character(phase_b_analytes)
  )  # 10 -> 20000, 40 -> 25000

  if (VERBOSE) message("\n=== Phase B5: n_analytes sub-grids ===")

  phase_b_summaries <- list()

  for (na_val in phase_b_analytes) {
    if (VERBOSE) message("\n--- n_analytes = ", na_val, " ---")

    # Build grid (same as main but clamp n_discriminatory to n_analytes)
    if (DEV_MODE) {
      pb_grid <- make_dgp_grid(
        n_per_group      = if (!is.null(EXTRA_SIZES)) EXTRA_SIZES else c(25L, 50L),
        n_discriminatory = 5L,
        effect_size      = c(0, 0.7, 1.2),
        correlation      = c("none", "block_0.5"),
        censoring        = 0.01,
        imbalance        = 1L
      )
    } else {
      # Clamp n_discriminatory values to those that fit within n_analytes
      disc_vals <- c(3L, 5L, 10L)
      disc_vals <- disc_vals[disc_vals <= na_val]
      pb_grid <- make_dgp_grid(
        n_per_group      = if (!is.null(EXTRA_SIZES)) EXTRA_SIZES else c(20L, 50L, 100L),
        n_discriminatory = disc_vals,
        effect_size      = c(0, 0.3, 0.7, 1.2),
        correlation      = c("none", "block_0.5", "block_0.8"),
        censoring        = c(0.01, 0.20),
        imbalance        = c(1L, 2L, 3L)
      )
    }

    # Apply offset for Phase B
    pb_grid$cell_id <- pb_grid$cell_id + phase_b_offsets[as.character(na_val)]

    if (VERBOSE) {
      message("DGP grid: ", nrow(pb_grid), " cells x ", N_REPS, " reps = ",
              nrow(pb_grid) * N_REPS, " total runs")
    }

    pb_cache_prefix <- paste0(BENCHMARK, "_phaseB_na", na_val)

    pb_summary <- run_benchmark_plsda(
      dgp_grid     = pb_grid,
      n_reps       = N_REPS,
      base_seed    = BASE_SEED,
      cache_dir    = CACHE_DIR,
      n_cores      = N_CORES,
      n_perm       = N_PERM,
      verbose      = VERBOSE,
      n_analytes   = na_val,
      cache_prefix = pb_cache_prefix
    )

    pb_annotated <- annotate_acceptance(pb_summary)

    # Save per-n_analytes summary
    pb_summary_file <- file.path(CACHE_DIR, paste0(pb_cache_prefix, "_summary.rds"))
    saveRDS(pb_annotated, pb_summary_file)
    if (VERBOSE) message("Phase B5 summary saved -> ", pb_summary_file)

    phase_b_summaries[[length(phase_b_summaries) + 1]] <- pb_annotated
  }

  # Merge Phase B summaries with main summary (add n_analytes=20 to existing)
  summary_file <- file.path(CACHE_DIR, paste0(BENCHMARK, "_summary.rds"))
  if (file.exists(summary_file)) {
    main_summary <- readRDS(summary_file)
    if (!"n_analytes" %in% names(main_summary)) {
      main_summary$n_analytes <- 20L
    }
    combined <- bind_rows(main_summary, bind_rows(phase_b_summaries))
    saveRDS(combined, summary_file)
    if (VERBOSE) message("Merged Phase B into main summary -> ", summary_file)
  } else {
    if (VERBOSE) message("Main summary not found; Phase B summaries saved separately.")
  }

  summary_annotated <- bind_rows(phase_b_summaries)

} else {
  # ---- Main grid execution -----------------------------------------------------
  if (VERBOSE) message("\nStarting benchmark execution...")

  summary_df <- run_benchmark_plsda(
    dgp_grid  = dgp_grid,
    n_reps    = N_REPS,
    base_seed = BASE_SEED,
    cache_dir = CACHE_DIR,
    n_cores   = N_CORES,
    n_perm    = N_PERM,
    verbose   = VERBOSE
  )

  summary_annotated <- annotate_acceptance(summary_df)

  # Save final summary (merge with existing when running extra sizes)
  summary_file <- file.path(CACHE_DIR, paste0(BENCHMARK, "_summary.rds"))
  if (!is.null(EXTRA_SIZES)) {
    merge_with_existing_summary(summary_annotated, summary_file, verbose = VERBOSE)
  } else {
    if (!"n_analytes" %in% names(summary_annotated)) {
      summary_annotated$n_analytes <- 20L
    }
    save_main_summary_preserving_supplements(summary_annotated, summary_file,
                                             verbose = VERBOSE)
  }
}



# Figures are generated by the standalone plot_B4.R script.
# Run: Rscript inst/simulations/plot_B4.R [--cache_dir PATH]


# ---- Formatted Summary Table -------------------------------------------------

format_plsda_summary <- function(df, digits = 3) {
  df %>%
    mutate(
      accuracy_fmt  = fmt_mcse_vec(balanced_accuracy, balanced_accuracy_mcse, digits),
      vip_auroc_fmt = fmt_mcse_vec(vip_auroc, vip_auroc_mcse, digits),
      q2_fmt        = fmt_mcse_vec(q2, q2_mcse, digits),
      stability_fmt = fmt_mcse_vec(feature_stability, feature_stability_mcse, digits),
      perm_rej_fmt  = fmt_mcse_vec(perm_rejection_rate, perm_rejection_mcse, digits)
    )
}

summary_formatted <- format_plsda_summary(summary_annotated)

if (VERBOSE) {
  message("\n=== Benchmark 4: Summary ===\n")
  print(
    summary_formatted %>%
      select(any_of("n_analytes"),
             n_per_group, n_discriminatory, effect_size, correlation,
             censoring, imbalance, method, n_reps,
             accuracy_fmt, vip_auroc_fmt, q2_fmt, stability_fmt, perm_rej_fmt,
             starts_with("pass_")) %>%
      as.data.frame(),
    right = FALSE
  )
}

if (VERBOSE) {
  message("\n=== Benchmark 4 complete ===")
  if (exists("fig_dir")) message("Figures -> ", fig_dir)
  message("Summary -> ", summary_file)
}
