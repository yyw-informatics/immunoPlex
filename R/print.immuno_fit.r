# S3 print method for immuno_fit objects
# Displays concise model information including convergence status and fit statistics

#' @title Print method for immuno_fit
#' @description Shows key information about the fitted model.
#' @param x   An `immuno_fit` object.
#' @param ... Unused.
#' @return `x` invisibly.
#' @export
print.immuno_fit <- function(x, ...) {
  cat("immunoPlex Model (immuno_fit)\n")
  cat("Family:", x$family, "\n")
  if (!is.null(x$estimand)) {
    cat("Estimand:", x$estimand, "\n")
  }
  cat("Converged:", x$converged, "\n")
  
  if (x$converged) {
    cat("AIC:", round(x$aic, 2), "\n")
    cat("Censored obs:", x$n_cens_lod, "left,", x$n_cens_ulod, "right\n")
  } else {
    cat("Model failed to converge\n")
  }
  
  cat("Use summary(x) for details\n")
  invisible(x)
}
