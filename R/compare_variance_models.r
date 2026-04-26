# Variance model comparison for Gaussian glmmTMB fits
# Compares homogeneous vs heterogeneous dispersion models using LRT, AIC, BIC

#' Compare homogeneous vs heterogeneous variance models
#'
#' Fits both a homogeneous (\code{dispformula = ~1}) and a heterogeneous
#' dispersion model, then compares them using a likelihood ratio test (LRT),
#' AIC, and BIC. Uses \code{\link{fit_one}} with \code{family = "gaussian"}
#' internally, inheriting the optimizer cascade.
#'
#' @param data A data frame containing the response and predictor columns.
#' @param fixed Character string of fixed effects (e.g. \code{"timepoint + age"}).
#' @param random Character string of random effects (e.g. \code{"(1|subject_id)"}).
#' @param dispformula_hetero Formula for the heterogeneous dispersion model.
#'   Default \code{~timepoint}.
#' @param label Character label for the analyte (used in output).
#' @param lrt_alpha Numeric significance threshold for the LRT. Default 0.05.
#'
#' @return A list with components:
#' \describe{
#'   \item{model_homo}{\code{immuno_fit} object with homogeneous dispersion.}
#'   \item{model_hetero}{\code{immuno_fit} object with heterogeneous dispersion.}
#'   \item{comparison}{A tibble with AIC, BIC, logLik, LRT statistic, df,
#'     and p-value for both models.}
#'   \item{preferred_model}{Character: \code{"homogeneous"} or
#'     \code{"heterogeneous"}.}
#'   \item{preference_reason}{Character explanation of the selection.}
#' }
#'
#' @export
compare_variance_models <- function(data,
                                    fixed,
                                    random = "(1|subject_id)",
                                    dispformula_hetero = ~timepoint,
                                    label = "",
                                    lrt_alpha = 0.05) {

  if (!is.data.frame(data)) stop("'data' must be a data frame")

  # Fit homogeneous variance model
  fit_homo <- tryCatch(
    fit_one(data, family = "gaussian", fixed = fixed, random = random,
            dispformula = ~1),
    error = function(e) NULL
  )

  # Fit heterogeneous variance model
  fit_hetero <- tryCatch(
    fit_one(data, family = "gaussian", fixed = fixed, random = random,
            dispformula = dispformula_hetero),
    error = function(e) NULL
  )

  homo_ok <- !is.null(fit_homo) && isTRUE(fit_homo$converged)
  hetero_ok <- !is.null(fit_hetero) && isTRUE(fit_hetero$converged)

  result <- list(
    label = label,
    model_homo = fit_homo,
    model_hetero = fit_hetero,
    comparison = NULL,
    preferred_model = NA_character_,
    preference_reason = NA_character_
  )

  if (homo_ok && hetero_ok) {
    aic_homo <- fit_homo$aic
    aic_hetero <- fit_hetero$aic
    bic_homo <- fit_homo$bic
    bic_hetero <- fit_hetero$bic
    ll_homo <- fit_homo$logLik
    ll_hetero <- fit_hetero$logLik

    # LRT: hetero is the more complex model
    lrt_stat <- 2 * (ll_hetero - ll_homo)
    # df = number of extra dispersion parameters
    disp_vars <- all.vars(dispformula_hetero)
    if (length(disp_vars) > 0 && disp_vars[1] %in% names(data)) {
      df_diff <- length(unique(data[[disp_vars[1]]])) - 1L
    } else {
      # Fallback: count parameters from model objects
      df_diff <- max(1L, length(stats::coef(fit_hetero$model)$disp) -
                          length(stats::coef(fit_homo$model)$disp))
    }
    lrt_p <- stats::pchisq(lrt_stat, df = df_diff, lower.tail = FALSE)

    result$comparison <- tibble::tibble(
      label = label,
      aic_homo = aic_homo, aic_hetero = aic_hetero,
      aic_diff = aic_hetero - aic_homo,
      bic_homo = bic_homo, bic_hetero = bic_hetero,
      bic_diff = bic_hetero - bic_homo,
      loglik_homo = ll_homo, loglik_hetero = ll_hetero,
      lrt_statistic = lrt_stat, lrt_df = df_diff, lrt_pvalue = lrt_p
    )

    # Selection criteria
    if (is.na(lrt_p) || is.na(aic_hetero) || is.na(aic_homo)) {
      result$preferred_model <- "homogeneous"
      result$preference_reason <- "comparison statistics unavailable, prefer simpler model"
    } else if (lrt_p < lrt_alpha && aic_hetero < aic_homo) {
      result$preferred_model <- "heterogeneous"
      result$preference_reason <- sprintf(
        "LRT significant (p=%.3f) and AIC favors hetero (delta=%.1f)",
        lrt_p, aic_hetero - aic_homo)
    } else {
      result$preferred_model <- "homogeneous"
      if (lrt_p >= lrt_alpha) {
        result$preference_reason <- sprintf(
          "LRT not significant (p=%.3f), prefer simpler model", lrt_p)
      } else {
        result$preference_reason <- sprintf(
          "BIC favors homogeneous (delta=%.1f)", bic_hetero - bic_homo)
      }
    }

  } else if (homo_ok && !hetero_ok) {
    result$preferred_model <- "homogeneous"
    result$preference_reason <- "heterogeneous model failed to converge"
  } else if (!homo_ok && hetero_ok) {
    result$preferred_model <- "heterogeneous"
    result$preference_reason <- "homogeneous model failed to converge"
  } else {
    result$preferred_model <- NA_character_
    result$preference_reason <- "both models failed to converge"
  }

  result
}
