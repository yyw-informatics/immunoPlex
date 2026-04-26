#!/usr/bin/env Rscript
# benchmark_preprocessing.R — Benchmark 0: Preprocessing Recovery (v4)
#
# ADEMP Structure:
#   Aims:    Does preprocessing improve downstream group-effect estimation when
#            data contain known technical artifacts (plate/lot/background)?
#   DGP:     simulate_immunoassay() with factorial grid over artifact severity
#            (none/moderate/heavy), censoring rate (10/25/40%), sample size
#            (40/100/200), and effect size (0/0.5/1.0).
#   Estimand: Marginal group effect on log-concentration for a single analyte.
#   Methods: (1) Raw, (2) LOD/2, (3) LOD/sqrt(2), (4) plate-corrected,
#            (5) Tobit (censoring-aware), (6) plate random-effect + LOD/2,
#            (7) Tobit + plate fixed-effect.
#   Performance: bias, coverage, RMSE, power, convergence, FPR under null.
#
# Usage:
#   Rscript benchmark_preprocessing.R [--n_reps N] [--n_cores N] [--cache_dir DIR]
#
# Results are cached as .rds files; re-running tops up to the target replication
# count without re-running existing replications.
#
# v2 changes (from report B0_preprocessing_recovery.md):
#   - Replaced halfmin with LOD/sqrt(2) (report anomaly #2: methods identical)
#   - Fixed plate-corrected ordering: LOD/2 first, then centering (anomaly #3)
#   - Added Tobit as censoring-aware reference method (missing methods)
#   - Added plate RE method for fair plate-effect comparison (method fairness)
#   - Added 25% censoring level and effect size 1.0 (missing DGP conditions)
#   - Default reps 200 -> 500 for tighter MCSEs (Monte Carlo precision)
#   - Added power and RMSE figures (missing figures)
#   - Cache prefix changed to "preprocessing_v2" (grid/method incompatibility)
#
# v3 changes (backbone harmonization):
#   - Added n_subjects = 200 (backbone alignment: 40/100/200)
#   - Grid: 54 -> 81 cells
#   - Cache prefix changed to "preprocessing_v3"
#
# v4 changes (benchmark extensions Phase A):
#   - Added Method 7: Tobit + plate FE (survreg with factor(plate_id))
#   - Cache prefix changed to "preprocessing_v4"

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

# Parse command-line arguments
cli_args <- commandArgs(trailingOnly = TRUE)
parse_cli_arg <- function(flag, default) {
  idx <- which(cli_args == flag)
  if (length(idx) > 0 && idx < length(cli_args)) {
    return(cli_args[idx + 1])
  }
  default
}

N_REPS     <- as.integer(parse_cli_arg("--n_reps",    "500"))
N_CORES    <- as.integer(parse_cli_arg("--n_cores",   "1"))
CACHE_DIR  <- parse_cli_arg("--cache_dir", file.path(script_dir, "cache"))
EXTRA_N    <- parse_cli_arg("--extra_n", NULL)
PHASE_B    <- parse_cli_arg("--phase_b", NULL)
PHASE_C    <- parse_cli_arg("--phase_c", NULL)
BASE_SEED  <- 20240001L
BENCHMARK  <- "preprocessing_v4"
VERBOSE    <- TRUE

# Extra-n supplement mode: run only the specified sample sizes with offset cell_ids
EXTRA_SIZES  <- if (!is.null(EXTRA_N)) as.integer(strsplit(EXTRA_N, ",")[[1]]) else NULL
EXTRA_OFFSET <- 10000L

# Phase B extension mode: run only Phase B sub-grids (skip main grid)
PHASE_B_OFFSET <- 20000L
# Phase C extension mode: parameter tweaks (c1 = null topup, c5 = Wald vs LRT)
PHASE_C_OFFSET <- 30000L

if (!dir.exists(CACHE_DIR)) dir.create(CACHE_DIR, recursive = TRUE)


# ---- DGP Grid ---------------------------------------------------------------

# Artifact severity presets
artifact_presets <- list(
  none = list(
    n_plates = 1, plate_sd = 0,
    n_lots = 1, lot_sd = 0,
    bg_shape = 0, bg_rate = 1,
    n_replicates = 1, replicate_sd = 0
  ),
  moderate = list(
    n_plates = 4, plate_sd = 0.10,
    n_lots = 2, lot_sd = 0.08,
    bg_shape = 1.5, bg_rate = 1,
    n_replicates = 2, replicate_sd = 0.08
  ),
  heavy = list(
    n_plates = 4, plate_sd = 0.15,
    n_lots = 2, lot_sd = 0.10,
    bg_shape = 2, bg_rate = 1,
    n_replicates = 2, replicate_sd = 0.10
  )
)

dgp_grid <- make_dgp_grid(
  artifact_severity = c("none", "moderate", "heavy"),
  censoring_target  = c(0.10, 0.25, 0.40),
  n_subjects        = if (!is.null(EXTRA_SIZES)) EXTRA_SIZES else c(40L, 100L, 200L),
  effect_size       = c(0, 0.5, 1.0)
)
if (!is.null(EXTRA_SIZES)) dgp_grid$cell_id <- dgp_grid$cell_id + EXTRA_OFFSET

if (VERBOSE) {
  message("=== Benchmark 0: Preprocessing Recovery (v4) ===")
  message("DGP grid: ", nrow(dgp_grid), " cells x ", N_REPS, " reps = ",
          nrow(dgp_grid) * N_REPS, " total runs")
  message("Methods:  7 (Raw, LOD/2, LOD/sqrt2, Plate-corrected, Tobit, Plate RE, Tobit+Plate FE)")
  message("Cache dir: ", CACHE_DIR)
}


# ---- Single-Analyte Replication Function ------------------------------------

#' Run one replication for a single DGP cell
#'
#' Generates data, applies all 7 methods, fits the model for each,
#' and returns a one-row-per-method data.frame.
#'
#' @param seed         RNG seed.
#' @param n_subjects   Sample size.
#' @param effect_size  True group effect (log-scale); 0 = null.
#' @param censoring_target Target censoring fraction.
#' @param artifact_preset  Named list of assay artifact parameters.
#' @return data.frame with columns: method, estimate, se, ci_lo, ci_hi,
#'   p_value, converged, true_effect.
run_one_rep <- function(seed, n_subjects, effect_size, censoring_target,
                        artifact_preset) {

  # --- Simulate data ---
  sim <- simulate_immunoassay(
    n_subjects     = n_subjects,
    n_timepoints   = 1L,
    n_analytes     = 1L,
    design         = "cross_sectional",
    group_levels   = c("control", "treatment"),
    group_effects  = effect_size,
    signal_analytes = if (effect_size != 0) 1L else 0L,
    effect_direction = "up",
    re_intercept_sd = 0,
    residual_sd     = 1.0,
    lod_quantile    = censoring_target,
    # Assay parameters from preset
    n_plates       = artifact_preset$n_plates,
    plate_sd       = artifact_preset$plate_sd,
    n_lots         = artifact_preset$n_lots,
    lot_sd         = artifact_preset$lot_sd,
    bg_shape       = artifact_preset$bg_shape,
    bg_rate        = artifact_preset$bg_rate,
    n_replicates   = artifact_preset$n_replicates,
    replicate_sd   = artifact_preset$replicate_sd,
    seed           = seed
  )

  d <- as.data.frame(sim$data)

  # True group effect (from the simulation)
  true_effect <- effect_size

  # --- Aggregate technical replicates (mean of log-scale values) ---
  if (artifact_preset$n_replicates > 1) {
    d_agg <- d %>%
      filter(!is.na(value), qc_flag == "pass") %>%
      group_by(subject_id, cytokine, group, timepoint) %>%
      summarise(
        value     = mean(value, na.rm = TRUE),
        value_raw = mean(value_raw, na.rm = TRUE),
        cens_lod  = any(cens_lod),
        cens_ulod = any(cens_ulod),
        lod       = first(lod),
        ulod      = first(ulod),
        plate_id  = first(plate_id),
        .groups   = "drop"
      )
    d_agg <- as.data.frame(d_agg)
  } else {
    d_agg <- d %>%
      filter(!is.na(value), qc_flag == "pass")
    d_agg <- as.data.frame(d_agg)
  }

  # Common formula (cross-sectional, no random effects)
  fit_formula <- "group"

  # --- Method 1: Raw (no preprocessing, use observed values as-is) ---
  m1 <- fit_method(d_agg, fit_formula, "Raw")

  # --- Method 2: LOD/2 substitution ---
  d_half <- d_agg
  below_lod <- !is.na(d_half$cens_lod) & d_half$cens_lod
  if (any(below_lod)) {
    d_half$value[below_lod] <- log(d_half$lod[below_lod] / 2)
  }
  m2 <- fit_method(d_half, fit_formula, "Preprocess (half)")

  # --- Method 3: LOD/sqrt(2) substitution ---
  # Replaces halfmin (report anomaly #2: LOD/2 and halfmin were nearly identical
  # because min-detected ~ LOD; LOD/sqrt(2) provides a meaningfully different
  # substitution point between LOD and LOD/2)
  d_sqrt2 <- d_agg
  below_lod_s2 <- !is.na(d_sqrt2$cens_lod) & d_sqrt2$cens_lod
  if (any(below_lod_s2)) {
    d_sqrt2$value[below_lod_s2] <- log(d_sqrt2$lod[below_lod_s2] / sqrt(2))
  }
  m3 <- fit_method(d_sqrt2, fit_formula, "Preprocess (sqrt2)")

  # --- Method 4: Plate-corrected + LOD/2 ---
  # v2 fix (report anomaly #3): apply LOD/2 FIRST, then plate centering.
  # In v1 the order was reversed — LOD/2 after centering overwrote the
  # plate-corrected values for censored observations, neutralising the
  # plate correction for exactly the observations that needed it most.
  d_plate <- d_agg
  below_lod_p <- !is.na(d_plate$cens_lod) & d_plate$cens_lod
  if (any(below_lod_p)) {
    d_plate$value[below_lod_p] <- log(d_plate$lod[below_lod_p] / 2)
  }
  if (length(unique(d_plate$plate_id)) > 1) {
    plate_means <- tapply(d_plate$value, d_plate$plate_id, mean, na.rm = TRUE)
    grand_mean  <- mean(d_plate$value, na.rm = TRUE)
    d_plate$value <- d_plate$value -
      plate_means[as.character(d_plate$plate_id)] + grand_mean
  }
  m4 <- fit_method(d_plate, fit_formula, "Plate-corrected")

  # --- Method 5: Tobit (censoring-aware, no substitution) ---
  # Uses the censoring flags directly via survival::survreg; does not
  # replace censored values with a constant. This is the reference method
  # the report flagged as missing from B0.
  m5 <- fit_method(d_agg, fit_formula, "Tobit", family = "tobit")

  # --- Method 6: Plate random effect + LOD/2 ---
  # Models the plate effect inside the likelihood rather than pre-subtracting.
  # Report flagged that pre-subtraction may be unfair because the downstream
  # model does not account for uncertainty in the plate-effect estimate.
  d_plate_re <- d_agg
  below_lod_re <- !is.na(d_plate_re$cens_lod) & d_plate_re$cens_lod
  if (any(below_lod_re)) {
    d_plate_re$value[below_lod_re] <- log(d_plate_re$lod[below_lod_re] / 2)
  }
  m6 <- if (length(unique(d_plate_re$plate_id)) > 1) {
    fit_method(d_plate_re, fit_formula, "Plate RE",
               random = "(1|plate_id)")
  } else {
    # Single plate (artifact = "none"): identical to LOD/2
    fit_method(d_plate_re, fit_formula, "Plate RE")
  }

  # --- Method 7: Tobit + plate fixed effect ---
  # Combines censoring-aware Tobit likelihood with plate adjustment via
  # fixed-effect plate indicators (factor(plate_id)) inside survreg.
  # Falls back to plain Tobit when only 1 plate (artifact = "none").
  m7 <- if (length(unique(d_agg$plate_id)) > 1) {
    d_tobit_fe <- d_agg
    d_tobit_fe$plate_id <- factor(d_tobit_fe$plate_id)
    fit_method(d_tobit_fe, "group + plate_id", "Tobit + plate FE",
               family = "tobit")
  } else {
    # Single plate: identical to plain Tobit
    fit_method(d_agg, fit_formula, "Tobit + plate FE", family = "tobit")
  }

  # --- Combine ---
  results <- rbind(m1, m2, m3, m4, m5, m6, m7)
  results$true_effect <- true_effect
  results
}


#' Fit a single method and extract coefficient summary
#'
#' Calls fit_one() with the specified family and extracts the group
#' coefficient estimate, SE, CI, and p-value.
#'
#' @param d_method    data.frame ready for fit_one().
#' @param formula_str Fixed-effects formula string.
#' @param method_name Character label for this method.
#' @param family      Model family passed to fit_one(): "gaussian" (default)
#'                    or "tobit".
#' @param random      Random-effects formula string, or NULL.
#' @return One-row data.frame.
fit_method <- function(d_method, formula_str, method_name,
                       family = "gaussian", random = NULL,
                       compute_lrt = getOption("immunoplex.compute_lrt", FALSE)) {

  # Default return on failure
  fail_row <- data.frame(
    method       = method_name,
    estimate     = NA_real_,
    se           = NA_real_,
    ci_lo        = NA_real_,
    ci_hi        = NA_real_,
    p_value      = NA_real_,
    p_value_lrt  = NA_real_,
    converged    = FALSE,
    stringsAsFactors = FALSE
  )

  tryCatch({
    # Temporarily lower random-effects thresholds for plate RE
    # (only 4 plates, well below the default min_subjects = 30)
    old_opts <- NULL
    if (!is.null(random) && nzchar(random)) {
      old_opts <- options(fit_one.min_subjects = 2L, fit_one.min_reps = 1L)
      on.exit(options(old_opts), add = TRUE)
    }

    fit <- fit_one(
      dat    = d_method,
      family = family,
      fixed  = formula_str,
      random = random
    )

    if (!fit$converged) return(fail_row)

    pval_lrt <- NA_real_

    # Extract group coefficient — dispatch on model class
    if (inherits(fit$model, "survreg")) {
      # Tobit / AFT: coefficient table from survival::survreg
      coef_tab <- summary(fit$model)$table
      group_rows <- grep("^group", rownames(coef_tab))
      if (length(group_rows) == 0) return(fail_row)
      est  <- coef_tab[group_rows[1], "Value"]
      se   <- coef_tab[group_rows[1], "Std. Error"]
      pval <- coef_tab[group_rows[1], "p"]

      # Optional LRT: fit a null model that drops group terms and compare
      if (isTRUE(compute_lrt)) {
        pval_lrt <- tryCatch({
          # Build null formula by removing every term involving `group` via
          # the formula-terms machinery. Handles interactions, transformed
          # terms, and whitespace quirks correctly.
          full_fml  <- stats::as.formula(paste("~", formula_str))
          labs      <- attr(stats::terms(full_fml), "term.labels")
          keep_labs <- labs[vapply(labs, function(lab) {
            !("group" %in% all.vars(stats::reformulate(lab)))
          }, logical(1))]
          null_rhs  <- if (length(keep_labs) == 0) {
            "1"
          } else {
            paste(keep_labs, collapse = " + ")
          }
          null_fit <- fit_one(
            dat    = d_method,
            family = family,
            fixed  = null_rhs,
            random = random
          )
          if (!null_fit$converged) {
            NA_real_
          } else {
            # anova(null, full) for survreg returns a table with Pr(>Chi)
            an <- anova(null_fit$model, fit$model)
            # Find the p-value column
            p_col <- grep("^Pr\\(>", colnames(an), value = TRUE)
            if (length(p_col) == 0) {
              NA_real_
            } else {
              # Use the last row (full model vs null)
              as.numeric(an[nrow(an), p_col[1]])
            }
          }
        }, error = function(e) NA_real_)
      }

    } else if (inherits(fit$model, "glmmTMB")) {
      # Gaussian glmmTMB (with or without random effects)
      coef_tab <- summary(fit$model)$coefficients$cond
      group_rows <- grep("^group", rownames(coef_tab))
      if (length(group_rows) == 0) return(fail_row)
      est  <- coef_tab[group_rows[1], "Estimate"]
      se   <- coef_tab[group_rows[1], "Std. Error"]
      pval <- coef_tab[group_rows[1], "Pr(>|z|)"]

    } else {
      return(fail_row)
    }

    data.frame(
      method       = method_name,
      estimate     = est,
      se           = se,
      ci_lo        = est - 1.96 * se,
      ci_hi        = est + 1.96 * se,
      p_value      = pval,
      p_value_lrt  = pval_lrt,
      converged    = TRUE,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    fail_row
  })
}


# ---- Run Benchmark ----------------------------------------------------------

run_benchmark_preprocessing <- function(dgp_grid, n_reps, base_seed,
                                        cache_dir, n_cores, verbose) {

  all_summaries <- list()

  for (i in seq_len(nrow(dgp_grid))) {
    row <- dgp_grid[i, ]
    art_preset <- artifact_presets[[row$artifact_severity]]
    cell_label <- paste0("cell_", row$cell_id, " [",
                         row$artifact_severity, "/cens=",
                         row$censoring_target, "/n=",
                         row$n_subjects, "/eff=",
                         row$effect_size, "]")

    # Build the per-replication function (closure over DGP parameters)
    rep_fn <- function(seed) {
      run_one_rep(
        seed             = seed,
        n_subjects       = row$n_subjects,
        effect_size      = row$effect_size,
        censoring_target = row$censoring_target,
        artifact_preset  = art_preset
      )
    }

    cache_file <- file.path(cache_dir, paste0(BENCHMARK, "_cell", row$cell_id, ".rds"))

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
      m_reps <- raw_results[!is.na(raw_results$method) & raw_results$method == m, ]
      agg <- aggregate_cell(m_reps, truth = row$effect_size)
      agg$method             <- m
      agg$artifact_severity  <- row$artifact_severity
      agg$censoring_target   <- row$censoring_target
      agg$n_subjects         <- row$n_subjects
      agg$effect_size        <- row$effect_size
      agg$cell_id            <- row$cell_id
      all_summaries[[length(all_summaries) + 1]] <- agg
    }
  }

  bind_rows(all_summaries)
}


# ---- Execute ----------------------------------------------------------------

if (!is.null(PHASE_C)) {

  # ---- Phase C Extensions (parameter tweaks) --------------------------------

  if (VERBOSE) message("\n=== Phase C Extensions (mode: ", PHASE_C, ") ===")

  if (identical(PHASE_C, "c1")) {

    # ---- Phase C1: 1000 reps for null cells (top-up existing v4 caches) -----
    # Filter full grid to effect_size == 0 (27 cells) and target 1000 reps.
    # IMPORTANT: we preserve the *original* row index in the full dgp_grid so
    # the seed stream `base_seed + (i-1)*10000L` matches the earlier 500-rep
    # run and `run_benchmark_cell()` continues the sequence rather than
    # duplicating seeds.
    #
    # Cache-column skew: the first 500 reps in each `preprocessing_v4_cellN.rds`
    # were generated before `fit_method` grew the `p_value_lrt` column. The new
    # 500 reps added by this top-up will contain `p_value_lrt = NA` unless
    # `options(immunoplex.compute_lrt = TRUE)` is set (C1 leaves it off to
    # avoid needless compute — LRT is a C5 concern). `run_benchmark_cell`
    # NA-fills missing columns on merge, so the combined 1000-rep cache is
    # consistent. Downstream aggregation uses `is.finite()` filtering, so the
    # denominator for rejection_rate remains 1000 while rejection_rate_lrt (if
    # ever computed) would have a 500-rep denominator — flag for anyone
    # tempted to read the raw cache directly.

    if (VERBOSE) message("\n--- Phase C1: Null-cell top-up to 1000 reps ---")

    null_idx <- which(dgp_grid$effect_size == 0)
    null_grid <- dgp_grid[null_idx, , drop = FALSE]
    null_grid$orig_index <- null_idx

    if (VERBOSE) {
      message("Null cells: ", nrow(null_grid),
              " (target ", 1000L, " reps each, topping up from ~500)")
    }

    all_summaries_c1 <- list()
    for (j in seq_len(nrow(null_grid))) {
      row <- null_grid[j, ]
      art_preset <- artifact_presets[[row$artifact_severity]]
      cell_label <- paste0("C1_cell_", row$cell_id, " [",
                           row$artifact_severity, "/cens=",
                           row$censoring_target, "/n=",
                           row$n_subjects, "/eff=",
                           row$effect_size, "]")

      rep_fn <- local({
        ns <- row$n_subjects; es <- row$effect_size
        ct <- row$censoring_target; ap <- art_preset
        function(seed) run_one_rep(seed, ns, es, ct, ap)
      })

      # Reuse the existing v4 cache file so previous 500 reps are preserved.
      cache_file <- file.path(CACHE_DIR,
                              paste0(BENCHMARK, "_cell", row$cell_id, ".rds"))

      raw_results <- run_benchmark_cell(
        rep_fn     = rep_fn,
        n_reps     = 1000L,
        base_seed  = BASE_SEED + (row$orig_index - 1) * 10000L,
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
        agg$method             <- m
        agg$artifact_severity  <- row$artifact_severity
        agg$censoring_target   <- row$censoring_target
        agg$n_subjects         <- row$n_subjects
        agg$effect_size        <- row$effect_size
        agg$cell_id            <- row$cell_id
        all_summaries_c1[[length(all_summaries_c1) + 1]] <- agg
      }
    }

    summary_c1 <- bind_rows(all_summaries_c1)
    summary_c1_annotated <- summary_c1 %>%
      annotate_acceptance(truth = 0, type1_bounds = c(0.025, 0.075))

    c1_file <- file.path(CACHE_DIR, paste0(BENCHMARK, "_c1_summary.rds"))
    saveRDS(summary_c1_annotated, c1_file)
    if (VERBOSE) message("Phase C1 summary saved -> ", c1_file)

    # Also merge into main summary (top-up overwrites the null rows)
    summary_file <- file.path(CACHE_DIR, paste0(BENCHMARK, "_summary.rds"))
    merge_with_existing_summary(summary_c1_annotated, summary_file,
                                verbose = VERBOSE)
    if (VERBOSE) message("Phase C1 merged -> ", summary_file)

    quit(save = "no", status = 0)
  }

  if (identical(PHASE_C, "c5")) {

    # ---- Phase C5: Wald vs LRT for Tobit (targeted sub-grid) ----------------
    # Target cells: n_subjects == 200 & censoring_target == 0.40.
    # Cache prefix "preprocessing_v4_c5" with offset cell_ids so v4 caches
    # (which lack the new p_value_lrt column) are not overwritten.

    if (VERBOSE) message("\n--- Phase C5: Wald vs LRT for Tobit ---")

    # Enable LRT computation in fit_method for C5 only (opt-in).
    # Script exits via quit() at end of block, so no restoration needed.
    options(immunoplex.compute_lrt = TRUE)

    c5_idx <- which(dgp_grid$n_subjects == 200L &
                    dgp_grid$censoring_target == 0.40)
    c5_grid <- dgp_grid[c5_idx, , drop = FALSE]
    c5_grid$orig_index <- c5_idx
    c5_grid$cell_id    <- c5_grid$cell_id + PHASE_C_OFFSET

    if (VERBOSE) {
      message("C5 cells: ", nrow(c5_grid),
              " (n=200, cens=0.40) x ", N_REPS, " reps")
    }

    all_summaries_c5 <- list()
    for (j in seq_len(nrow(c5_grid))) {
      row <- c5_grid[j, ]
      art_preset <- artifact_presets[[row$artifact_severity]]
      cell_label <- paste0("C5_cell_", row$cell_id, " [",
                           row$artifact_severity, "/cens=",
                           row$censoring_target, "/n=",
                           row$n_subjects, "/eff=",
                           row$effect_size, "]")

      rep_fn <- local({
        ns <- row$n_subjects; es <- row$effect_size
        ct <- row$censoring_target; ap <- art_preset
        function(seed) run_one_rep(seed, ns, es, ct, ap)
      })

      cache_file <- file.path(CACHE_DIR,
                              paste0(BENCHMARK, "_c5_cell", row$cell_id, ".rds"))

      raw_results <- run_benchmark_cell(
        rep_fn     = rep_fn,
        n_reps     = N_REPS,
        base_seed  = BASE_SEED + (row$orig_index - 1) * 10000L,
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

        # Add LRT rejection rate (Tobit methods only)
        if ("p_value_lrt" %in% names(m_reps)) {
          pvals_lrt <- m_reps$p_value_lrt[is.finite(m_reps$p_value_lrt)]
          if (length(pvals_lrt) > 0) {
            rej_lrt <- mean(pvals_lrt < 0.05)
            agg$rejection_rate_lrt <- rej_lrt
            agg$rejection_mcse_lrt <- sqrt(rej_lrt * (1 - rej_lrt) /
                                           length(pvals_lrt))
          } else {
            agg$rejection_rate_lrt <- NA_real_
            agg$rejection_mcse_lrt <- NA_real_
          }
        }

        agg$method             <- m
        agg$artifact_severity  <- row$artifact_severity
        agg$censoring_target   <- row$censoring_target
        agg$n_subjects         <- row$n_subjects
        agg$effect_size        <- row$effect_size
        agg$cell_id            <- row$cell_id
        all_summaries_c5[[length(all_summaries_c5) + 1]] <- agg
      }
    }

    summary_c5 <- bind_rows(all_summaries_c5)

    # Apply the same null/alt acceptance annotation used by the main grid so
    # downstream reporting can read pass_type1 / pass_bias / pass_coverage
    # columns consistently.
    c5_null <- summary_c5 %>%
      filter(effect_size == 0) %>%
      annotate_acceptance(truth = 0, type1_bounds = c(0.025, 0.075))
    c5_alt  <- summary_c5 %>% filter(effect_size != 0)
    c5_alt_list <- lapply(unique(c5_alt$effect_size), function(es) {
      c5_alt %>%
        filter(effect_size == es) %>%
        annotate_acceptance(truth = es, bias_ratio_threshold = 0.1,
                            coverage_bounds = c(0.93, 0.97))
    })
    summary_c5_annotated <- bind_rows(c5_null, bind_rows(c5_alt_list))

    c5_file <- file.path(CACHE_DIR, paste0(BENCHMARK, "_c5_summary.rds"))
    saveRDS(summary_c5_annotated, c5_file)
    if (VERBOSE) message("Phase C5 summary saved -> ", c5_file)

    quit(save = "no", status = 0)
  }

  stop("Unknown PHASE_C mode: ", PHASE_C, " (expected 'c1' or 'c5')")
}

if (is.null(PHASE_B)) {

  # ---- Main Grid Execution (Phase A) ----------------------------------------

  if (VERBOSE) message("\nStarting benchmark execution...")

  summary_df <- run_benchmark_preprocessing(
    dgp_grid  = dgp_grid,
    n_reps    = N_REPS,
    base_seed = BASE_SEED,
    cache_dir = CACHE_DIR,
    n_cores   = N_CORES,
    verbose   = VERBOSE
  )

  # Add acceptance annotations — per effect-size group for non-null cells
  summary_null <- summary_df %>%
    filter(effect_size == 0) %>%
    annotate_acceptance(truth = 0, type1_bounds = c(0.025, 0.075))

  summary_alt <- summary_df %>%
    filter(effect_size != 0)

  summary_alt_list <- lapply(unique(summary_alt$effect_size), function(es) {
    summary_alt %>%
      filter(effect_size == es) %>%
      annotate_acceptance(truth = es, bias_ratio_threshold = 0.1,
                          coverage_bounds = c(0.93, 0.97))
  })
  summary_alt <- bind_rows(summary_alt_list)

  summary_annotated <- bind_rows(summary_null, summary_alt)

  # Save final summary (merge with existing when running extra sizes)
  summary_file <- file.path(CACHE_DIR, paste0(BENCHMARK, "_summary.rds"))
  if (!is.null(EXTRA_SIZES)) {
    merge_with_existing_summary(summary_annotated, summary_file, verbose = VERBOSE)
  } else {
    save_main_summary_preserving_supplements(summary_annotated, summary_file,
                                             verbose = VERBOSE)
  }


  # Figures are generated by the standalone plot_B0.R script.
  # Run: Rscript inst/simulations/plot_B0.R [--cache_dir PATH]

  # ---- Formatted Summary Table ------------------------------------------------

  summary_formatted <- format_summary_table(summary_annotated)

  # Print concise table
  if (VERBOSE) {
    message("\n=== Benchmark 0: Summary ===\n")
    print(
      summary_formatted %>%
        select(artifact_severity, censoring_target, n_subjects, effect_size,
               method, n_reps, bias_fmt, coverage_fmt, rejection_fmt,
               starts_with("pass_")) %>%
        as.data.frame(),
      right = FALSE
    )
  }

} else {

  # ---- Phase B Extensions ----------------------------------------------------

  if (VERBOSE) message("\n=== Phase B Extensions ===")


  # ---- Phase B1: Artifact-only (0% censoring) --------------------------------
  #
  # Isolates artifact effects from censoring effects by setting censoring to

  # near-zero. Uses moderate and heavy artifact presets only (none has no
  # artifacts to isolate). Sub-grid: 2 artifact x 1 censoring(0) x 3 n x 3 eff
  # = 18 cells.

  dgp_grid_b1 <- make_dgp_grid(
    artifact_severity = c("moderate", "heavy"),
    censoring_target  = 0,
    n_subjects        = c(40L, 100L, 200L),
    effect_size       = c(0, 0.5, 1.0)
  )
  dgp_grid_b1$cell_id <- dgp_grid_b1$cell_id + PHASE_B_OFFSET

  if (VERBOSE) {
    message("\n--- Phase B1: Artifact-only (0% censoring) ---")
    message("DGP grid: ", nrow(dgp_grid_b1), " cells x ", N_REPS, " reps = ",
            nrow(dgp_grid_b1) * N_REPS, " total runs")
  }

  all_summaries_b1 <- list()

  for (i in seq_len(nrow(dgp_grid_b1))) {
    row <- dgp_grid_b1[i, ]
    art_preset <- artifact_presets[[row$artifact_severity]]
    cell_label <- paste0("B1_cell_", row$cell_id, " [",
                         row$artifact_severity, "/cens=0/n=",
                         row$n_subjects, "/eff=",
                         row$effect_size, "]")

    # Build the per-replication function
    # Use lod_quantile = 0.001 to get near-zero censoring
    rep_fn <- function(seed) {
      run_one_rep(
        seed             = seed,
        n_subjects       = row$n_subjects,
        effect_size      = row$effect_size,
        censoring_target = 0.001,
        artifact_preset  = art_preset
      )
    }

    cache_file <- file.path(CACHE_DIR,
                            paste0(BENCHMARK, "_phase_b1_cell", row$cell_id, ".rds"))

    raw_results <- run_benchmark_cell(
      rep_fn     = rep_fn,
      n_reps     = N_REPS,
      base_seed  = BASE_SEED + (i - 1) * 10000L,
      cache_file = cache_file,
      n_cores    = N_CORES,
      cell_label = cell_label,
      verbose    = VERBOSE
    )

    # Aggregate per method
    methods <- unique(raw_results$method[!is.na(raw_results$method)])
    for (m in methods) {
      m_reps <- raw_results[!is.na(raw_results$method) & raw_results$method == m, ]
      agg <- aggregate_cell(m_reps, truth = row$effect_size)
      agg$method             <- m
      agg$artifact_severity  <- row$artifact_severity
      agg$censoring_target   <- row$censoring_target
      agg$n_subjects         <- row$n_subjects
      agg$effect_size        <- row$effect_size
      agg$cell_id            <- row$cell_id
      all_summaries_b1[[length(all_summaries_b1) + 1]] <- agg
    }
  }

  summary_b1 <- bind_rows(all_summaries_b1)

  # Annotate B1
  summary_b1_null <- summary_b1 %>%
    filter(effect_size == 0) %>%
    annotate_acceptance(truth = 0, type1_bounds = c(0.025, 0.075))

  summary_b1_alt <- summary_b1 %>%
    filter(effect_size != 0)
  summary_b1_alt_list <- lapply(unique(summary_b1_alt$effect_size), function(es) {
    summary_b1_alt %>%
      filter(effect_size == es) %>%
      annotate_acceptance(truth = es, bias_ratio_threshold = 0.1,
                          coverage_bounds = c(0.93, 0.97))
  })
  summary_b1_alt <- bind_rows(summary_b1_alt_list)
  summary_b1_annotated <- bind_rows(summary_b1_null, summary_b1_alt)

  # Merge B1 into existing summary
  summary_file <- file.path(CACHE_DIR, paste0(BENCHMARK, "_summary.rds"))
  merge_with_existing_summary(summary_b1_annotated, summary_file, verbose = VERBOSE)
  if (VERBOSE) message("Phase B1 summary merged -> ", summary_file)


  # ---- Phase B3: Multi-analyte plate correction (10 analytes) ----------------
  #
  # Tests whether plate-correction methods maintain FDR control when applied
  # across a panel of 10 analytes with shared plate effects. Uses BH FDR
  # correction within each replication. Sub-grid: 2 artifact x 2 censoring x
  # 3 n x 2 effect = 24 cells.

  dgp_grid_b3 <- make_dgp_grid(
    artifact_severity = c("moderate", "heavy"),
    censoring_target  = c(0.10, 0.40),
    n_subjects        = c(40L, 100L, 200L),
    effect_size       = c(0, 0.5)
  )
  dgp_grid_b3$cell_id <- dgp_grid_b3$cell_id + PHASE_B_OFFSET + 1000L

  if (VERBOSE) {
    message("\n--- Phase B3: Multi-analyte plate correction (10 analytes) ---")
    message("DGP grid: ", nrow(dgp_grid_b3), " cells x ", N_REPS, " reps = ",
            nrow(dgp_grid_b3) * N_REPS, " total runs")
  }

  #' Run one replication for the multi-analyte panel (Phase B3)
  #'
  #' Generates data with 10 analytes and shared plate effects, fits all 7
  #' preprocessing methods for each analyte independently, and returns one
  #' row per (method x analyte).
  #'
  #' @param seed              RNG seed.
  #' @param n_subjects        Sample size.
  #' @param effect_size       True group effect (log-scale); 0 = null.
  #' @param censoring_target  Target censoring fraction.
  #' @param artifact_preset   Named list of assay artifact parameters.
  #' @param n_signal_use      Number of signal analytes (default: computed from
  #'   effect_size).
  #' @return data.frame with columns: method, analyte, estimate, se, ci_lo,
  #'   ci_hi, p_value, converged, is_signal, true_effect.
  run_one_rep_panel <- function(seed, n_subjects, effect_size, censoring_target,
                                artifact_preset, n_signal_use = NULL) {

    # Determine number of signal analytes: half the panel when effect > 0
    if (is.null(n_signal_use)) {
      n_signal_use <- if (effect_size != 0) 5L else 0L
    }

    sim <- simulate_immunoassay(
      n_subjects      = n_subjects,
      n_timepoints    = 1L,
      n_analytes      = 10L,
      design          = "cross_sectional",
      group_levels    = c("control", "treatment"),
      group_effects   = effect_size,
      signal_analytes = n_signal_use,
      effect_direction = "up",
      re_intercept_sd = 0,
      residual_sd     = 1.0,
      lod_quantile    = censoring_target,
      # Assay parameters from preset (shared plate effects across analytes)
      n_plates        = artifact_preset$n_plates,
      plate_sd        = artifact_preset$plate_sd,
      n_lots          = artifact_preset$n_lots,
      lot_sd          = artifact_preset$lot_sd,
      bg_shape        = artifact_preset$bg_shape,
      bg_rate         = artifact_preset$bg_rate,
      n_replicates    = artifact_preset$n_replicates,
      replicate_sd    = artifact_preset$replicate_sd,
      seed            = seed
    )

    d  <- as.data.frame(sim$data)
    tr <- as.data.frame(sim$truth)
    analytes <- unique(d$cytokine)

    # Determine which analytes are signal from truth table
    signal_analytes_vec <- if ("signal_analyte" %in% names(tr)) {
      unique(d$cytokine[tr$signal_analyte])
    } else {
      character(0)
    }

    all_results <- list()

    for (analyte_name in analytes) {
      d_a <- d[d$cytokine == analyte_name, , drop = FALSE]
      is_signal <- analyte_name %in% signal_analytes_vec
      true_eff  <- if (is_signal) effect_size else 0

      # Aggregate technical replicates
      if (artifact_preset$n_replicates > 1) {
        d_agg <- d_a %>%
          filter(!is.na(value), qc_flag == "pass") %>%
          group_by(subject_id, cytokine, group, timepoint) %>%
          summarise(
            value     = mean(value, na.rm = TRUE),
            value_raw = mean(value_raw, na.rm = TRUE),
            cens_lod  = any(cens_lod),
            cens_ulod = any(cens_ulod),
            lod       = first(lod),
            ulod      = first(ulod),
            plate_id  = first(plate_id),
            .groups   = "drop"
          )
        d_agg <- as.data.frame(d_agg)
      } else {
        d_agg <- d_a %>%
          filter(!is.na(value), qc_flag == "pass")
        d_agg <- as.data.frame(d_agg)
      }

      fit_formula <- "group"

      # Method 1: Raw
      m1 <- fit_method(d_agg, fit_formula, "Raw")

      # Method 2: LOD/2 substitution
      d_half <- d_agg
      below_lod <- !is.na(d_half$cens_lod) & d_half$cens_lod
      if (any(below_lod)) {
        d_half$value[below_lod] <- log(d_half$lod[below_lod] / 2)
      }
      m2 <- fit_method(d_half, fit_formula, "Preprocess (half)")

      # Method 3: LOD/sqrt(2) substitution
      d_sqrt2 <- d_agg
      below_lod_s2 <- !is.na(d_sqrt2$cens_lod) & d_sqrt2$cens_lod
      if (any(below_lod_s2)) {
        d_sqrt2$value[below_lod_s2] <- log(d_sqrt2$lod[below_lod_s2] / sqrt(2))
      }
      m3 <- fit_method(d_sqrt2, fit_formula, "Preprocess (sqrt2)")

      # Method 4: Plate-corrected + LOD/2
      d_plate <- d_agg
      below_lod_p <- !is.na(d_plate$cens_lod) & d_plate$cens_lod
      if (any(below_lod_p)) {
        d_plate$value[below_lod_p] <- log(d_plate$lod[below_lod_p] / 2)
      }
      if (length(unique(d_plate$plate_id)) > 1) {
        plate_means <- tapply(d_plate$value, d_plate$plate_id, mean, na.rm = TRUE)
        grand_mean  <- mean(d_plate$value, na.rm = TRUE)
        d_plate$value <- d_plate$value -
          plate_means[as.character(d_plate$plate_id)] + grand_mean
      }
      m4 <- fit_method(d_plate, fit_formula, "Plate-corrected")

      # Method 5: Tobit
      m5 <- fit_method(d_agg, fit_formula, "Tobit", family = "tobit")

      # Method 6: Plate RE + LOD/2
      d_plate_re <- d_agg
      below_lod_re <- !is.na(d_plate_re$cens_lod) & d_plate_re$cens_lod
      if (any(below_lod_re)) {
        d_plate_re$value[below_lod_re] <- log(d_plate_re$lod[below_lod_re] / 2)
      }
      m6 <- if (length(unique(d_plate_re$plate_id)) > 1) {
        fit_method(d_plate_re, fit_formula, "Plate RE",
                   random = "(1|plate_id)")
      } else {
        fit_method(d_plate_re, fit_formula, "Plate RE")
      }

      # Method 7: Tobit + plate FE
      m7 <- if (length(unique(d_agg$plate_id)) > 1) {
        d_tobit_fe <- d_agg
        d_tobit_fe$plate_id <- factor(d_tobit_fe$plate_id)
        fit_method(d_tobit_fe, "group + plate_id", "Tobit + plate FE",
                   family = "tobit")
      } else {
        fit_method(d_agg, fit_formula, "Tobit + plate FE", family = "tobit")
      }

      analyte_results <- rbind(m1, m2, m3, m4, m5, m6, m7)
      analyte_results$analyte     <- analyte_name
      analyte_results$is_signal   <- is_signal
      analyte_results$true_effect <- true_eff
      all_results[[length(all_results) + 1]] <- analyte_results
    }

    do.call(rbind, all_results)
  }


  all_summaries_b3_panel <- list()
  all_summaries_b3_analyte <- list()

  for (i in seq_len(nrow(dgp_grid_b3))) {
    row <- dgp_grid_b3[i, ]
    art_preset <- artifact_presets[[row$artifact_severity]]
    cell_label <- paste0("B3_cell_", row$cell_id, " [",
                         row$artifact_severity, "/cens=",
                         row$censoring_target, "/n=",
                         row$n_subjects, "/eff=",
                         row$effect_size, "]")

    rep_fn <- function(seed) {
      run_one_rep_panel(
        seed             = seed,
        n_subjects       = row$n_subjects,
        effect_size      = row$effect_size,
        censoring_target = row$censoring_target,
        artifact_preset  = art_preset
      )
    }

    cache_file <- file.path(CACHE_DIR,
                            paste0(BENCHMARK, "_phase_b3_cell", row$cell_id, ".rds"))

    raw_results <- run_benchmark_cell(
      rep_fn     = rep_fn,
      n_reps     = N_REPS,
      base_seed  = BASE_SEED + (i - 1) * 10000L,
      cache_file = cache_file,
      n_cores    = N_CORES,
      cell_label = cell_label,
      verbose    = VERBOSE
    )

    # Panel-level aggregation (FDR, FWER, sensitivity via aggregate_panel_cell)
    panel_agg <- aggregate_panel_cell(raw_results)
    panel_agg$artifact_severity <- row$artifact_severity
    panel_agg$censoring_target  <- row$censoring_target
    panel_agg$n_subjects        <- row$n_subjects
    panel_agg$effect_size       <- row$effect_size
    panel_agg$cell_id           <- row$cell_id
    all_summaries_b3_panel[[length(all_summaries_b3_panel) + 1]] <- panel_agg

    # Per-analyte aggregation (bias, coverage, etc. per method x analyte)
    methods <- unique(raw_results$method[!is.na(raw_results$method)])
    analytes <- unique(raw_results$analyte[!is.na(raw_results$analyte)])
    for (m in methods) {
      for (a in analytes) {
        ma_reps <- raw_results[!is.na(raw_results$method) &
                               raw_results$method == m &
                               !is.na(raw_results$analyte) &
                               raw_results$analyte == a, , drop = FALSE]
        if (nrow(ma_reps) == 0) next
        truth_a <- if (nrow(ma_reps) > 0 && "true_effect" %in% names(ma_reps)) {
          ma_reps$true_effect[1]
        } else {
          0
        }
        agg <- aggregate_cell(ma_reps, truth = truth_a)
        agg$method             <- m
        agg$analyte            <- a
        agg$artifact_severity  <- row$artifact_severity
        agg$censoring_target   <- row$censoring_target
        agg$n_subjects         <- row$n_subjects
        agg$effect_size        <- row$effect_size
        agg$cell_id            <- row$cell_id
        all_summaries_b3_analyte[[length(all_summaries_b3_analyte) + 1]] <- agg
      }
    }
  }

  summary_b3_panel   <- bind_rows(all_summaries_b3_panel)
  summary_b3_analyte <- bind_rows(all_summaries_b3_analyte)

  # Save B3 panel summary
  panel_summary_file <- file.path(CACHE_DIR,
                                  paste0(BENCHMARK, "_phase_b3_panel_summary.rds"))
  saveRDS(summary_b3_panel, panel_summary_file)
  if (VERBOSE) message("Phase B3 panel summary saved -> ", panel_summary_file)

  # Save B3 per-analyte summary
  analyte_summary_file <- file.path(CACHE_DIR,
                                    paste0(BENCHMARK, "_phase_b3_analyte_summary.rds"))
  saveRDS(summary_b3_analyte, analyte_summary_file)
  if (VERBOSE) message("Phase B3 per-analyte summary saved -> ", analyte_summary_file)

  # Merge B3 per-analyte summary into main summary (for cross-phase comparison)
  summary_file <- file.path(CACHE_DIR, paste0(BENCHMARK, "_summary.rds"))
  # Annotate B3 analyte summary for merging
  summary_b3_null <- summary_b3_analyte %>%
    filter(effect_size == 0) %>%
    annotate_acceptance(truth = 0, type1_bounds = c(0.025, 0.075))
  summary_b3_alt <- summary_b3_analyte %>%
    filter(effect_size != 0)
  summary_b3_alt_list <- lapply(unique(summary_b3_alt$effect_size), function(es) {
    summary_b3_alt %>%
      filter(effect_size == es) %>%
      annotate_acceptance(truth = es, bias_ratio_threshold = 0.1,
                          coverage_bounds = c(0.93, 0.97))
  })
  summary_b3_alt <- bind_rows(summary_b3_alt_list)
  summary_b3_annotated <- bind_rows(summary_b3_null, summary_b3_alt)
  merge_with_existing_summary(summary_b3_annotated, summary_file, verbose = VERBOSE)
  if (VERBOSE) message("Phase B3 analyte summary merged -> ", summary_file)

  # Print B3 panel summary
  if (VERBOSE) {
    message("\n=== Phase B3: Panel Summary ===\n")
    print(as.data.frame(summary_b3_panel), right = FALSE)
  }

}  # end Phase B vs main grid conditional

if (VERBOSE) {
  summary_file <- file.path(CACHE_DIR, paste0(BENCHMARK, "_summary.rds"))
  message("\n=== Benchmark 0 complete ===")
  if (exists("fig_dir")) message("Figures -> ", fig_dir)
  message("Summary -> ", summary_file)
}
