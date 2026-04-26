# Package startup hooks for immunoPlex
# Sets sensible defaults on package load and provides informative startup messages

# Declare global variables used in NSE contexts to avoid R CMD check NOTEs
utils::globalVariables(c(
  # Variables used in dplyr/tidyr NSE contexts
  "sample_id", "concentration", "lod", "timepoint", "subject_id", "cytokine", 
  "detected", "p_mcnemar", "q_mcnemar", "censored", "pct_censored",
  # Variables used in ggplot2 aes() contexts  
  "fitted", "residuals", "censor_status", "model", "aic", "approach_type",
  "family", "delta_aic", "immunoplex_example",
  # Variables used in ANCOVA contexts
  "outcome_rank", "covariate_rank", "omega_sq_partial",
  "omega_magnitude", "q_value", "q_omnibus", "fill_color",
  "analyte", "is_significant",
  # Variables used in winsorize_by_group / LOO sensitivity
  "q1", "q3", "iqr", "upper_bound", "lower_bound",
  "excluded_subject", "influential", "pct_change",
  "full_estimate", "loo_estimate", "difference",
  # Variables used in replicate QC
  "replicate_cv", "replicate_discordant", "replicate_action", "replicate_n",
  "pair_cv", "pair_abs_diff", "pair_mean", "pair_n", "pair_sd",
  "concentration_raw", "population_median", "pop_median_diff", "pop_mad_diff",
  "mad_cutoff", "flag_cv", "flag_mad", "dist_to_median", "is_farther",
  # Variables used in well-level detection (replicate QC)
  "well_failure", "well_failure_type", "well_score", "well_sign_bias",
  "well_lod_fraction", "frac_flagged", "n_flagged_analytes",
  "n_total_analytes", "population_mad", "population_sd", "good_well",
  "sign_bias", "lod_fraction", "mean_sign_z", "score_ratio",
  ".scale", ".z", ".is_bad_well", ".is_rwf", "n_wells", "min_score", "max_score"
))

.onLoad <- function(libname, pkgname) {
  op <- options()
  op_pkg <- list(
    fit_one.min_subjects = 30L,   # keep random effects if ≥ 30 clusters
    fit_one.min_reps     = 3L     # and ≥ 3 rows per cluster
  )
  toset <- !(names(op_pkg) %in% names(op))
  if (any(toset)) options(op_pkg[toset])
  invisible()
}

.onAttach <- function(libname, pkgname) {
  v <- utils::packageDescription(pkgname, fields = "Version")
  msg <- sprintf(
    "immunoPlex %s loaded | Random-effects thresholds: %d subjects, %d reps/subject\n",
    v,
    getOption("fit_one.min_subjects"),
    getOption("fit_one.min_reps")
  )
  packageStartupMessage(msg)
}