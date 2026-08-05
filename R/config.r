# Configuration management for immunoPlex package
# Centralized parameter control for model thresholds, plotting defaults, and diagnostic settings

#' immunoPlex Package Configuration
#'
#' @description
#' Central configuration file for immunoPlex package parameters.
#' This allows users to modify key thresholds and defaults in one place.
#'
#' @name immunoplex_config
NULL

#' Get immunoPlex configuration
#'
#' @description
#' Returns a list of current configuration parameters for the immunoPlex package.
#' These parameters control various thresholds and defaults throughout the package.
#'
#' @return Named list of configuration parameters:
#' \describe{
#'   \item{random_effects}{List of thresholds for random effects modeling}
#'   \item{auto_family}{Thresholds for automatic family selection}
#'   \item{plotting}{Default plotting parameters}
#'   \item{dharma}{DHARMa simulation defaults}
#'   \item{ancova}{Rank-based ANCOVA defaults (outlier removal, FDR, prevalence)}
#'   \item{testing}{Parameters used in test scenarios}
#' }
#'
#' @examples
#' # Get current configuration
#' config <- get_immunoplex_config()
#' 
#' # View random effects thresholds
#' config$random_effects
#' 
#' # Modify configuration (affects subsequent function calls)
#' options(fit_one.min_subjects = 25)  # Adjust minimum subject threshold
#' 
#' @export
get_immunoplex_config <- function() {
  list(
    # Random effects modeling thresholds
    random_effects = list(
      min_subjects = getOption("fit_one.min_subjects", 30L),
      min_reps = getOption("fit_one.min_reps", 3L)
    ),
    
    # Automatic family selection thresholds
    auto_family = list(
      high_censoring_threshold = getOption("immunoplex.high_censoring_pct", 70),
      high_skewness_threshold = getOption("immunoplex.high_skewness", 2)
    ),
    
    # Plotting defaults
    plotting = list(
      default_point_size = getOption("immunoplex.point_size", 1.2),
      default_alpha_cens = getOption("immunoplex.alpha_cens", 0.7),
      outlier_threshold = getOption("immunoplex.outlier_threshold", 3),
      loess_color = getOption("immunoplex.loess_color", "gray50"),
      text_size = getOption("immunoplex.text_size", 3)
    ),
    
    # DHARMa simulation defaults
    dharma = list(
      default_nsim = as.integer(getOption("immunoplex.dharma_nsim", 1000L)),
      test_nsim = as.integer(getOption("immunoplex.test_dharma_nsim", 100L))
    ),
    
    # Rank-based ANCOVA defaults
    ancova = list(
      iqr_multiplier = getOption("immunoplex.ancova_iqr_multiplier", 3),
      min_n = as.integer(getOption("immunoplex.ancova_min_n", 12L)),
      fdr_method = getOption("immunoplex.ancova_fdr_method", "BH"),
      fdr_threshold = getOption("immunoplex.ancova_fdr_threshold", 0.05),
      prevalence_threshold = getOption("immunoplex.ancova_prevalence_threshold", 0.03)
    ),

    # Gaussian glmmTMB defaults
    glmmTMB = list(
      gradient_warn = getOption("immunoplex.glmmTMB_gradient_warn", 0.001),
      gradient_fail = getOption("immunoplex.glmmTMB_gradient_fail", 0.01),
      loo_influence_threshold = getOption("immunoplex.glmmTMB_loo_threshold", 20),
      loo_subject_col = getOption("immunoplex.glmmTMB_loo_subject_col", "subject_id")
    ),

    # Replicate QC defaults
    replicate_qc = list(
      cv_threshold = getOption("immunoplex.replicate_cv_threshold", 25),
      mad_threshold = getOption("immunoplex.replicate_mad_threshold", 3),
      flag_rule = getOption("immunoplex.replicate_flag_rule", "cv"),
      action = getOption("immunoplex.replicate_action", "flag"),
      well_detection = getOption("immunoplex.replicate_well_detection", "analyte"),
      well_params = getOption("immunoplex.replicate_well_params", list())
    ),

    # Testing scenarios (not for end users, but for test consistency)
    testing = list(
      small_dataset_size = getOption("immunoplex.test_small_n", 10),
      medium_dataset_size = getOption("immunoplex.test_medium_n", 50),
      large_dataset_size = getOption("immunoplex.test_large_n", 500),
      default_lod = getOption("immunoplex.test_lod", 1.0),
      default_ulod = getOption("immunoplex.test_ulod", 10.0),
      performance_dataset_size = getOption("immunoplex.test_perf_n", 200)
    )
  )
}

#' Set immunoPlex configuration
#'
#' @description
#' Convenience function to set multiple immunoPlex configuration parameters at once.
#' Individual parameters can also be set using base R `options()`.
#'
#' @param random_effects_min_subjects Minimum number of subjects required for random effects (default: 30)
#' @param random_effects_min_reps Minimum replicates per subject for random effects (default: 3)
#' @param high_censoring_threshold Percentage above which to prefer Tobit models (default: 70)
#' @param high_skewness_threshold Skewness above which to prefer Gamma models (default: 2)
#' @param default_dharma_nsim Default number of DHARMa simulations (default: 1000)
#' @param plotting_point_size Default point size for plots (default: 1.2)
#' @param plotting_alpha_cens Default alpha for censored points (default: 0.7)
#' @param ancova_iqr_multiplier IQR multiplier for ANCOVA outlier removal (default: 3)
#' @param ancova_min_n Minimum observations per analyte for ANCOVA (default: 12)
#' @param ancova_fdr_method FDR correction method for ANCOVA (default: "BH")
#' @param ancova_fdr_threshold FDR significance threshold for ANCOVA (default: 0.05)
#' @param ancova_prevalence_threshold Minimum prevalence for conditional covariates (default: 0.03)
#' @param glmmTMB_gradient_warn Gradient threshold for convergence warning (default: 0.001)
#' @param glmmTMB_gradient_fail Gradient threshold for convergence failure (default: 0.01)
#' @param glmmTMB_loo_influence_threshold Percent change threshold for LOO influential subject detection (default: 20)
#' @param replicate_cv_threshold CV threshold (%) for replicate QC flagging (default: 25)
#' @param replicate_mad_threshold MAD multiplier for replicate QC flagging (default: 3)
#' @param replicate_flag_rule Flag rule for replicate QC: "cv", "mad", "either", "both" (default: "cv")
#' @param replicate_action Action for flagged replicates: "flag", "drop_both", "drop_farther", "winsorize" (default: "flag")
#' @param replicate_well_detection Well detection mode: "analyte" (per-analyte) or "sample" (sample-level well failure detection) (default: "analyte")
#' @param replicate_well_params Named list of well detection tuning parameters (default: list())
#'
#' @examples
#' # Adjust thresholds for small sample sizes
#' set_immunoplex_config(
#'   random_effects_min_subjects = 20,
#'   random_effects_min_reps = 2
#' )
#' 
#' # Modify censoring threshold for model selection
#' set_immunoplex_config(high_censoring_threshold = 50)
#' 
#' # Faster DHARMa for interactive use
#' set_immunoplex_config(default_dharma_nsim = 250)
#' 
#' @export
set_immunoplex_config <- function(random_effects_min_subjects = NULL,
                                  random_effects_min_reps = NULL,
                                  high_censoring_threshold = NULL,
                                  high_skewness_threshold = NULL,
                                  default_dharma_nsim = NULL,
                                  plotting_point_size = NULL,
                                  plotting_alpha_cens = NULL,
                                  ancova_iqr_multiplier = NULL,
                                  ancova_min_n = NULL,
                                  ancova_fdr_method = NULL,
                                  ancova_fdr_threshold = NULL,
                                  ancova_prevalence_threshold = NULL,
                                  glmmTMB_gradient_warn = NULL,
                                  glmmTMB_gradient_fail = NULL,
                                  glmmTMB_loo_influence_threshold = NULL,
                                  replicate_cv_threshold = NULL,
                                  replicate_mad_threshold = NULL,
                                  replicate_flag_rule = NULL,
                                  replicate_action = NULL,
                                  replicate_well_detection = NULL,
                                  replicate_well_params = NULL) {
  
  if (!is.null(random_effects_min_subjects)) {
    options(fit_one.min_subjects = as.integer(random_effects_min_subjects))
  }
  
  if (!is.null(random_effects_min_reps)) {
    options(fit_one.min_reps = as.integer(random_effects_min_reps))
  }
  
  if (!is.null(high_censoring_threshold)) {
    options(immunoplex.high_censoring_pct = high_censoring_threshold)
  }
  
  if (!is.null(high_skewness_threshold)) {
    options(immunoplex.high_skewness = high_skewness_threshold)
  }
  
  if (!is.null(default_dharma_nsim)) {
    options(immunoplex.dharma_nsim = as.integer(default_dharma_nsim))
  }
  
  if (!is.null(plotting_point_size)) {
    options(immunoplex.point_size = plotting_point_size)
  }
  
  if (!is.null(plotting_alpha_cens)) {
    options(immunoplex.alpha_cens = plotting_alpha_cens)
  }

  if (!is.null(ancova_iqr_multiplier)) {
    options(immunoplex.ancova_iqr_multiplier = ancova_iqr_multiplier)
  }
  if (!is.null(ancova_min_n)) {
    options(immunoplex.ancova_min_n = as.integer(ancova_min_n))
  }
  if (!is.null(ancova_fdr_method)) {
    options(immunoplex.ancova_fdr_method = ancova_fdr_method)
  }
  if (!is.null(ancova_fdr_threshold)) {
    options(immunoplex.ancova_fdr_threshold = ancova_fdr_threshold)
  }
  if (!is.null(ancova_prevalence_threshold)) {
    options(immunoplex.ancova_prevalence_threshold = ancova_prevalence_threshold)
  }

  if (!is.null(glmmTMB_gradient_warn)) {
    options(immunoplex.glmmTMB_gradient_warn = glmmTMB_gradient_warn)
  }
  if (!is.null(glmmTMB_gradient_fail)) {
    options(immunoplex.glmmTMB_gradient_fail = glmmTMB_gradient_fail)
  }
  if (!is.null(glmmTMB_loo_influence_threshold)) {
    options(immunoplex.glmmTMB_loo_threshold = as.numeric(glmmTMB_loo_influence_threshold))
  }

  if (!is.null(replicate_cv_threshold)) {
    options(immunoplex.replicate_cv_threshold = replicate_cv_threshold)
  }
  if (!is.null(replicate_mad_threshold)) {
    options(immunoplex.replicate_mad_threshold = replicate_mad_threshold)
  }
  if (!is.null(replicate_flag_rule)) {
    options(immunoplex.replicate_flag_rule = replicate_flag_rule)
  }
  if (!is.null(replicate_action)) {
    options(immunoplex.replicate_action = replicate_action)
  }
  if (!is.null(replicate_well_detection)) {
    options(immunoplex.replicate_well_detection = replicate_well_detection)
  }
  if (!is.null(replicate_well_params)) {
    options(immunoplex.replicate_well_params = replicate_well_params)
  }

  invisible(get_immunoplex_config())
}

#' Reset immunoPlex configuration to defaults
#'
#' @description
#' Resets all immunoPlex configuration parameters to their default values.
#'
#' @examples
#' # After experimenting with settings
#' reset_immunoplex_config()
#' 
#' @export
reset_immunoplex_config <- function() {
  # Remove all immunoplex options to restore defaults
  option_names <- c(
    "fit_one.min_subjects", "fit_one.min_reps",
    "immunoplex.high_censoring_pct", "immunoplex.high_skewness",
    "immunoplex.dharma_nsim", "immunoplex.test_dharma_nsim",
    "immunoplex.point_size", "immunoplex.alpha_cens",
    "immunoplex.outlier_threshold", "immunoplex.loess_color",
    "immunoplex.text_size", "immunoplex.test_small_n",
    "immunoplex.test_medium_n", "immunoplex.test_large_n",
    "immunoplex.test_lod", "immunoplex.test_ulod",
    "immunoplex.test_perf_n",
    "immunoplex.ancova_iqr_multiplier", "immunoplex.ancova_min_n",
    "immunoplex.ancova_fdr_method", "immunoplex.ancova_fdr_threshold",
    "immunoplex.ancova_prevalence_threshold",
    "immunoplex.glmmTMB_gradient_warn", "immunoplex.glmmTMB_gradient_fail",
    "immunoplex.glmmTMB_loo_threshold", "immunoplex.glmmTMB_loo_subject_col",
    "immunoplex.replicate_cv_threshold", "immunoplex.replicate_mad_threshold",
    "immunoplex.replicate_flag_rule", "immunoplex.replicate_action",
    "immunoplex.replicate_well_detection", "immunoplex.replicate_well_params"
  )
  
  # Set all to NULL to restore defaults
  opts <- as.list(rep(list(NULL), length(option_names)))
  names(opts) <- option_names
  do.call(options, opts)
  
  invisible(get_immunoplex_config())
}

# Internal helper for getting test configuration values
# (Not exported - for internal package use)
.get_test_config <- function() {
  list(
    small_n = getOption("immunoplex.test_small_n", 10),
    medium_n = getOption("immunoplex.test_medium_n", 50), 
    large_n = getOption("immunoplex.test_large_n", 500),
    perf_n = getOption("immunoplex.test_perf_n", 200),
    default_lod = getOption("immunoplex.test_lod", 1.0),
    default_ulod = getOption("immunoplex.test_ulod", 10.0),
    test_dharma_nsim = getOption("immunoplex.test_dharma_nsim", 100)
  )
}