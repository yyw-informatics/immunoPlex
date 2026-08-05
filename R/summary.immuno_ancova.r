# S3 summary methods for immuno_ancova and immuno_ancova_set objects

#' Summarize a single rank-based ANCOVA result
#'
#' @description
#' Full coefficient table, effect sizes with CIs, formula, and sample sizes.
#'
#' @param object An \code{immuno_ancova} object from \code{ancova_one()}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns a list with summary information.
#' @export
summary.immuno_ancova <- function(object, ...) {
  label <- if (!is.null(object$analyte)) object$analyte else "unnamed"

  cat("Rank-Based ANCOVA Summary\n")
  cat("=========================\n\n")
  cat("Analyte:", label, "\n")
  cat("Group variable:", object$group_var, sprintf("(%d levels: %s)\n",
      object$n_groups, paste(object$group_levels, collapse = ", ")))
  cat("Formula:", object$formula_used, "\n")
  cat("Observations:", object$n_obs, sprintf("(%d outliers removed)\n", object$n_removed))

  if (length(object$covariates_included) > 0) {
    cat("Covariates included:", paste(object$covariates_included, collapse = ", "), "\n")
  }
  if (length(object$covariates_dropped) > 0) {
    cat("Covariates dropped:", paste(object$covariates_dropped, collapse = ", "), "\n")
  }
  cat("\n")

  # Model fit
  if (!is.null(object$model)) {
    cat("Model Fit:\n")
    cat("----------\n")
    cat(sprintf("R-squared: %.4f\n", object$r_squared))
    omega_display <- max(object$omega_sq_partial, 0, na.rm = TRUE)
    cat(sprintf("Partial omega-squared: %.4f [%.4f, %.4f] (%s)\n",
                omega_display, max(object$omega_sq_ci[1], 0, na.rm = TRUE),
                object$omega_sq_ci[2], object$omega_magnitude))
    cat("\n")

    if (object$n_groups == 2) {
      cat("Group Effect (binary):\n")
      cat("----------------------\n")
      cat(sprintf("  Estimate: %.4f\n", object$group_effect))
      cat(sprintf("  SE:       %.4f\n", object$group_se))
      cat(sprintf("  p-value:  %.6f\n", object$group_p))
      cat(sprintf("  CI:       [%.4f, %.4f]\n", object$group_ci[1], object$group_ci[2]))
    } else {
      cat(sprintf("Omnibus F: %.3f (p = %.6f)\n\n", object$omnibus_f, object$omnibus_p))
      if (!is.null(object$contrasts)) {
        cat("Pairwise Contrasts:\n")
        cat("-------------------\n")
        print(object$contrasts, row.names = FALSE)
      }
    }

    if (!is.null(object$log_effect)) {
      cat("\nLog-Scale Validation:\n")
      cat("---------------------\n")
      cat(sprintf("  Log effect:   %.4f\n", object$log_effect))
      cat(sprintf("  Fold change:  %.3f\n", object$log_fold_change))
      cat(sprintf("  Log p-value:  %.6f\n", object$log_p))
    }

    cat("\nCoefficient Table:\n")
    print(summary(object$model)$coefficients)
  } else {
    cat("Model did not converge - no fit statistics available\n")
  }

  if (length(object$warnings) > 0) {
    cat("\nWarnings:\n")
    for (w in object$warnings) cat("  -", w, "\n")
  }

  invisible(list(
    analyte = label,
    n_obs = object$n_obs,
    n_removed = object$n_removed,
    r_squared = object$r_squared,
    omega_sq_partial = object$omega_sq_partial,
    omega_magnitude = object$omega_magnitude,
    group_effect = object$group_effect,
    group_p = object$group_p
  ))
}


#' Summarize a batch ANCOVA result set
#'
#' @description
#' Returns the full results tibble with all statistics.
#'
#' @param object An \code{immuno_ancova_set} object from \code{ancova_fit()}.
#' @param ... Additional arguments (currently unused).
#'
#' @return The results tibble (invisibly).
#' @export
summary.immuno_ancova_set <- function(object, ...) {
  cat("Rank-Based ANCOVA Batch Summary\n")
  cat("===============================\n\n")
  cat(sprintf("Analytes: %d | Significant: %d (FDR %s < %.2f)\n\n",
              object$n_analytes, object$n_significant,
              object$fdr_method, object$fdr_threshold))

  print(object$results)

  invisible(object$results)
}
