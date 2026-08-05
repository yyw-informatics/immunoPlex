# Summary method for immuno_fit objects
# Comprehensive model diagnostics and coefficient reporting for fitted immunoPlex models

#' Summarize an immuno_fit object
#'
#' @description
#' Provides a comprehensive summary of a fitted immunoPlex model including
#' model diagnostics, convergence information, estimand type, and coefficient estimates.
#'
#' @param object An `immuno_fit` object from `fit_one()`.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns a list with summary information, but primarily
#'   called for its side effect of printing a formatted summary.
#'
#' @examples
#' \dontrun{
#' fit <- fit_one(data, family = "gamma")
#' summary(fit)
#' }
#'
#' @export
summary.immuno_fit <- function(object, ...) {
  
  if (!inherits(object, "immuno_fit")) {
    stop("Object must be of class 'immuno_fit'")
  }
  
  cat("immunoPlex Model Summary\n")
  cat("========================\n\n")
  
  # Basic model information
  cat("Family:", toupper(object$family), "\n")
  if (!is.null(object$estimand)) {
    cat("Estimand:", object$estimand, "\n")
  }
  cat("Converged:", if (object$converged) "Yes" else "No", "\n")
  
  # Handle case where data_used might be NULL (small datasets)
  n_obs <- if (!is.null(object$data_used)) nrow(object$data_used) else object$n_obs
  if (is.null(n_obs)) n_obs <- object$n_cens_lod + object$n_cens_ulod  # fallback
  
  cat("Observations:", n_obs, "\n")
  cat("Left-censored:", object$n_cens_lod, 
      sprintf("(%.1f%%)", 100 * object$n_cens_lod / n_obs), "\n")
  if (object$n_cens_ulod > 0) {
    cat("Right-censored:", object$n_cens_ulod, 
        sprintf("(%.1f%%)", 100 * object$n_cens_ulod / n_obs), "\n")
  }
  cat("\n")
  
  # Model fit statistics
  if (object$converged) {
    cat("Model Fit:\n")
    cat("----------\n")
    cat("AIC:", sprintf("%.2f", object$aic), "\n")
    cat("BIC:", sprintf("%.2f", object$bic), "\n")
    cat("Log-likelihood:", sprintf("%.2f", object$logLik), "\n")
    cat("\n")
    
    # Model-specific details
    if (object$family == "gamma" && inherits(object$model, "glmmTMB")) {
      cat("Model Details (Gamma GLMM):\n")
      cat("---------------------------\n")
      m <- object$model
      cat("Convergence code:", m$fit$convergence, "\n")
      cat("Hessian positive definite:", isTRUE(m$sdr$pdHess), "\n")
      
      # Random effects info if present
      theta <- m$fit$par[names(m$fit$par) == "theta"]
      if (length(theta) > 0) {
        cat("Random effect variance:", sprintf("%.6f", theta^2), "\n")
        cat("Singular fit:", any(abs(theta) < 1e-4), "\n")
      }
      
      # Coefficient summary
      cat("\nFixed Effects:\n")
      print(summary(m)$coefficients$cond)
      
    } else if (object$family %in% c("tobit", "aft") && inherits(object$model, "survreg")) {
      cat("Model Details (Survival Regression):\n")
      cat("------------------------------------\n")
      m <- object$model
      cat("Iterations:", m$iter, "\n")
      cat("Distribution:", m$dist, "\n")
      
      # Coefficient summary  
      cat("\nCoefficients:\n")
      print(summary(m)$table)
      
    } else if (object$family == "tobit_censreg" && inherits(object$model, "censReg")) {
      cat("Model Details (censReg):\n")
      cat("------------------------\n")
      m <- object$model
      
      # Basic estimates
      cat("\nEstimates:\n")
      print(m$estimate)
      
    } else {
      cat("Model-specific details not available for", class(object$model)[1], "\n")
    }
    
  } else {
    cat("Model did not converge - no fit statistics available\n")
  }
  
  # Return summary invisibly
  invisible(list(
    family = object$family,
    estimand = object$estimand,
    converged = object$converged,
    n_obs = n_obs,
    n_cens_lod = object$n_cens_lod,
    n_cens_ulod = object$n_cens_ulod,
    aic = object$aic,
    bic = object$bic,
    logLik = object$logLik
  ))
}