# Extended diagnostics for glmmTMB models
# DHARMa simulation-based diagnostics, convergence checks, and random effects validation

# Internal helper: extract glmmTMB model from immuno_fit or raw model object
.extract_glmmTMB_model <- function(x) {
  if (inherits(x, "glmmTMB")) return(x)
  if (inherits(x, "immuno_fit") && inherits(x$model, "glmmTMB")) return(x$model)
  stop("Expected a glmmTMB model or immuno_fit object, got: ", class(x)[1])
}

#' Comprehensive DHARMa diagnostics for glmmTMB models
#'
#' Runs a full DHARMa simulation-based diagnostic suite on a fitted glmmTMB
#' model: uniformity (KS test), dispersion, outlier detection, and quantile
#' deviation. Extends the basic DHARMa checks in
#' \code{\link{compute_residual_diagnostics}} with configurable thresholds
#' and richer output.
#'
#' @param model A fitted \code{glmmTMB} model object or an \code{immuno_fit}
#'   object (the model is extracted automatically).
#' @param label Character label for the analyte (used in warning messages).
#' @param n_sim Integer number of DHARMa simulations. Default uses
#'   \code{getOption("immunoplex.dharma_nsim", 1000L)}.
#' @param alpha_uniformity Numeric significance threshold for the KS
#'   uniformity test. Default 0.01.
#' @param alpha_dispersion Numeric significance threshold for the dispersion
#'   test. Default 0.01.
#' @param alpha_outliers Numeric significance threshold for the outlier
#'   test. Default 0.01.
#' @param alpha_quantiles Numeric significance threshold for the quantile
#'   deviation test. Default 0.05.
#'
#' @return A list with components:
#' \describe{
#'   \item{tests}{Named list of test results, each with \code{test},
#'     \code{statistic} or \code{ratio}, \code{p_value}, and \code{passed}.}
#'   \item{residuals}{The DHARMa residual object (for plotting).}
#'   \item{passed_all}{Logical: did all tests pass?}
#'   \item{warnings}{Character vector of warning messages.}
#' }
#'
#' @export
dharma_diagnostics <- function(model,
                               label = "",
                               n_sim = NULL,
                               alpha_uniformity = 0.01,
                               alpha_dispersion = 0.01,
                               alpha_outliers = 0.01,
                               alpha_quantiles = 0.05) {

  if (!requireNamespace("DHARMa", quietly = TRUE)) {
    stop("DHARMa package required for dharma_diagnostics()")
  }

  model <- .extract_glmmTMB_model(model)

  if (is.null(n_sim)) {
    n_sim <- as.integer(getOption("immunoplex.dharma_nsim", 1000L))
  }

  diagnostics <- list(
    label = label,
    passed_all = TRUE,
    tests = list(),
    residuals = NULL,
    warnings = character()
  )

  # Create scaled residuals
  sim_res <- tryCatch(
    DHARMa::simulateResiduals(model, n = as.integer(n_sim),
                              plot = FALSE, refit = FALSE),
    error = function(e) {
      diagnostics$warnings <<- c(diagnostics$warnings,
                                  paste0("Failed to create DHARMa residuals: ", e$message))
      diagnostics$passed_all <<- FALSE
      NULL
    }
  )

  if (is.null(sim_res)) return(diagnostics)
  diagnostics$residuals <- sim_res

  # 1. Uniformity test (KS test)
  ks_test <- tryCatch(DHARMa::testUniformity(sim_res, plot = FALSE),
                       error = function(e) NULL)
  if (!is.null(ks_test)) {
    diagnostics$tests$uniformity <- list(
      test = "Kolmogorov-Smirnov uniformity",
      statistic = ks_test$statistic,
      p_value = ks_test$p.value,
      passed = ks_test$p.value >= alpha_uniformity
    )
    if (!diagnostics$tests$uniformity$passed) {
      diagnostics$warnings <- c(diagnostics$warnings,
                                 sprintf("Uniformity violated (KS p=%.3f)", ks_test$p.value))
    }
  }

  # 2. Dispersion test
  disp_test <- tryCatch(DHARMa::testDispersion(sim_res, plot = FALSE),
                          error = function(e) NULL)
  if (!is.null(disp_test)) {
    diagnostics$tests$dispersion <- list(
      test = "Dispersion",
      ratio = as.numeric(disp_test$statistic),
      p_value = disp_test$p.value,
      passed = disp_test$p.value >= alpha_dispersion
    )
    if (!diagnostics$tests$dispersion$passed) {
      diagnostics$warnings <- c(diagnostics$warnings,
                                 sprintf("Dispersion issue (ratio=%.2f, p=%.3f)",
                                         disp_test$statistic, disp_test$p.value))
    }
  }

  # 3. Outlier test
  outlier_test <- tryCatch(DHARMa::testOutliers(sim_res, plot = FALSE),
                            error = function(e) NULL)
  if (!is.null(outlier_test)) {
    diagnostics$tests$outliers <- list(
      test = "Outliers (bootstrap)",
      statistic = outlier_test$statistic,
      p_value = outlier_test$p.value,
      passed = outlier_test$p.value >= alpha_outliers
    )
    if (!diagnostics$tests$outliers$passed) {
      diagnostics$warnings <- c(diagnostics$warnings, "Significant outliers detected")
    }
  }

  # 4. Quantile deviation test (catches heteroscedasticity)
  quant_test <- tryCatch(DHARMa::testQuantiles(sim_res, plot = FALSE),
                          error = function(e) NULL)
  if (!is.null(quant_test) && !is.null(quant_test$p.value) && !is.na(quant_test$p.value)) {
    diagnostics$tests$quantiles <- list(
      test = "Quantile deviation (combined)",
      p_value = quant_test$p.value,
      passed = quant_test$p.value >= alpha_quantiles
    )
    if (!diagnostics$tests$quantiles$passed) {
      diagnostics$warnings <- c(diagnostics$warnings,
                                 sprintf("Quantile deviations detected (p=%.3f)",
                                         quant_test$p.value))
    }
  }

  # Summary: passed all if no test failed
  diagnostics$passed_all <- all(vapply(diagnostics$tests,
                                        function(x) x$passed, logical(1)))

  diagnostics
}


#' Check convergence and numerical stability of a glmmTMB model
#'
#' Examines Hessian positive-definiteness, gradient magnitude, and NA
#' standard errors to assess whether a glmmTMB model converged properly.
#'
#' @param model A fitted \code{glmmTMB} model object or an \code{immuno_fit}
#'   object.
#' @param gradient_warn Numeric gradient threshold for issuing a warning.
#'   Default 0.001.
#' @param gradient_fail Numeric gradient threshold for declaring convergence
#'   failure. Default 0.01.
#'
#' @return A list with components:
#' \describe{
#'   \item{converged}{Logical: did the model converge?}
#'   \item{hessian_pd}{Logical: is the Hessian positive-definite?}
#'   \item{max_gradient}{Numeric: maximum absolute gradient component.}
#'   \item{na_se_count}{Integer: number of NA standard errors.}
#'   \item{warnings}{Character vector of diagnostic warnings.}
#' }
#'
#' @export
check_convergence <- function(model,
                              gradient_warn = NULL,
                              gradient_fail = NULL) {

  model <- .extract_glmmTMB_model(model)

  if (is.null(gradient_warn)) {
    gradient_warn <- getOption("immunoplex.glmmTMB_gradient_warn", 0.001)
  }
  if (is.null(gradient_fail)) {
    gradient_fail <- getOption("immunoplex.glmmTMB_gradient_fail", 0.01)
  }

  diag <- list(
    converged = TRUE,
    hessian_pd = NA,
    max_gradient = NA_real_,
    na_se_count = 0L,
    warnings = character()
  )

  # Check Hessian positive-definiteness
  if (!is.null(model$sdr)) {
    diag$hessian_pd <- isTRUE(model$sdr$pdHess)
    if (!diag$hessian_pd) {
      diag$converged <- FALSE
      diag$warnings <- c(diag$warnings, "Hessian not positive-definite")
    }
  }

  # Check for boundary/singular fit
  if (isTRUE(model$sdr$boundary)) {
    diag$converged <- FALSE
    diag$warnings <- c(diag$warnings, "Boundary/singular fit detected")
  }

  # Check gradients
  grad <- tryCatch(model$fit$gradient, error = function(e) NULL)
  if (!is.null(grad)) {
    diag$max_gradient <- max(abs(grad), na.rm = TRUE)
    if (diag$max_gradient > gradient_warn) {
      diag$warnings <- c(diag$warnings,
                          sprintf("Large gradient: %.4f", diag$max_gradient))
    }
    if (diag$max_gradient > gradient_fail) {
      diag$converged <- FALSE
    }
  }

  # Check for NA standard errors
  summ <- tryCatch(summary(model), error = function(e) NULL)
  if (!is.null(summ)) {
    coef_tab <- summ$coefficients$cond
    if (!is.null(coef_tab)) {
      diag$na_se_count <- sum(is.na(coef_tab[, "Std. Error"]))
      if (diag$na_se_count > 0L) {
        diag$converged <- FALSE
        diag$warnings <- c(diag$warnings,
                            sprintf("%d coefficient(s) with NA standard errors",
                                    diag$na_se_count))
      }
    }
  }

  diag
}


#' Check random effects diagnostics for a glmmTMB model
#'
#' Extracts variance components, computes ICC, and tests normality of
#' random effects via Shapiro-Wilk. Auto-detects grouping factor(s) from
#' the model's \code{VarCorr()} rather than hardcoding a column name.
#'
#' @param model A fitted \code{glmmTMB} model object or an \code{immuno_fit}
#'   object.
#' @param zero_var_threshold Numeric variance below which a "near zero"
#'   warning is issued. Default 1e-6.
#' @param min_re_for_normality Integer minimum number of random effects
#'   required to run the Shapiro-Wilk test. Default 8.
#'
#' @return A list with components:
#' \describe{
#'   \item{has_random_effects}{Logical.}
#'   \item{re_variance}{Numeric: random effect variance.}
#'   \item{re_sd}{Numeric: random effect standard deviation.}
#'   \item{icc}{Numeric: intraclass correlation coefficient.}
#'   \item{re_normality_p}{Numeric: Shapiro-Wilk p-value (or NA).}
#'   \item{warnings}{Character vector.}
#' }
#'
#' @export
check_random_effects <- function(model,
                                 zero_var_threshold = 1e-6,
                                 min_re_for_normality = 8L) {

  model <- .extract_glmmTMB_model(model)

  diag <- list(
    has_random_effects = FALSE,
    re_variance = NA_real_,
    re_sd = NA_real_,
    icc = NA_real_,
    re_normality_p = NA_real_,
    warnings = character()
  )

  # Extract variance components (VarCorr is from nlme, re-exported by glmmTMB)
  var_comp <- tryCatch(glmmTMB::VarCorr(model), error = function(e) NULL)
  if (is.null(var_comp) || is.null(var_comp$cond) ||
      length(var_comp$cond) == 0L) {
    return(diag)
  }

  diag$has_random_effects <- TRUE

  # Auto-detect first grouping factor
  grp_name <- names(var_comp$cond)[1]
  re_sd <- attr(var_comp$cond[[grp_name]], "stddev")
  if (is.null(re_sd) || length(re_sd) == 0L) return(diag)

  # Use the first variance component (intercept)
  re_sd <- re_sd[1]
  re_var <- re_sd^2
  diag$re_variance <- re_var
  diag$re_sd <- re_sd

  # Near-zero variance warning

  if (re_var < zero_var_threshold) {
    diag$warnings <- c(diag$warnings,
                        "Random effect variance near zero - consider removing")
  }

  # Extract random effects for normality check (ranef is from nlme, re-exported by glmmTMB)
  re <- tryCatch(glmmTMB::ranef(model)$cond[[grp_name]], error = function(e) NULL)
  if (!is.null(re) && nrow(re) >= min_re_for_normality) {
    sw_test <- tryCatch(stats::shapiro.test(re[, 1]), error = function(e) NULL)
    if (!is.null(sw_test)) {
      diag$re_normality_p <- sw_test$p.value
      if (sw_test$p.value < 0.01) {
        diag$warnings <- c(diag$warnings,
                            sprintf("Random effects may not be normal (SW p=%.3f)",
                                    sw_test$p.value))
      }
    }
  }

  # Compute ICC
  resid_var <- tryCatch(stats::sigma(model)^2, error = function(e) NA_real_)
  if (!is.na(resid_var) && !is.na(re_var) && (re_var + resid_var) > 0) {
    diag$icc <- re_var / (re_var + resid_var)
  }

  diag
}
