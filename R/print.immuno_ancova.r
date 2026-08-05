# S3 print methods for immuno_ancova and immuno_ancova_set objects

#' @title Print method for immuno_ancova
#' @description One-line summary of a single rank-based ANCOVA result.
#' @param x An \code{immuno_ancova} object.
#' @param ... Unused.
#' @return \code{x} invisibly.
#' @export
print.immuno_ancova <- function(x, ...) {
  label <- if (!is.null(x$analyte)) x$analyte else "unnamed"

  if (x$n_groups == 2 && !is.na(x$group_p)) {
    sig <- if (x$group_p < 0.001) "***" else if (x$group_p < 0.01) "**" else if (x$group_p < 0.05) "*" else ""
    cat(sprintf("ANCOVA: %s | effect = %.2f (p = %.4f%s) | omega_sq = %.3f (%s)\n",
                label, x$group_effect, x$group_p, sig,
                max(x$omega_sq_partial, 0, na.rm = TRUE), x$omega_magnitude))
  } else if (x$n_groups > 2 && !is.na(x$omnibus_p)) {
    sig <- if (x$omnibus_p < 0.001) "***" else if (x$omnibus_p < 0.01) "**" else if (x$omnibus_p < 0.05) "*" else ""
    cat(sprintf("ANCOVA: %s | F = %.2f (p = %.4f%s) | omega_sq = %.3f (%s) | %d contrasts\n",
                label, x$omnibus_f, x$omnibus_p, sig,
                max(x$omega_sq_partial, 0, na.rm = TRUE), x$omega_magnitude,
                if (!is.null(x$contrasts)) nrow(x$contrasts) else 0L))
  } else {
    cat(sprintf("ANCOVA: %s | model failed (n = %d)\n", label, x$n_obs))
  }
  invisible(x)
}


#' @title Print method for immuno_ancova_set
#' @description Summary counts for a batch ANCOVA result set.
#' @param x An \code{immuno_ancova_set} object.
#' @param ... Unused.
#' @return \code{x} invisibly.
#' @export
print.immuno_ancova_set <- function(x, ...) {
  cat(sprintf("Rank-based ANCOVA: %d analytes, %d significant (FDR %s < %.2f)\n",
              x$n_analytes, x$n_significant, x$fdr_method, x$fdr_threshold))

  # Effect size distribution
  omega_vals <- x$results$omega_sq_partial
  omega_vals <- pmax(omega_vals, 0, na.rm = TRUE)
  n_large <- sum(omega_vals >= 0.14, na.rm = TRUE)
  n_medium <- sum(omega_vals >= 0.06 & omega_vals < 0.14, na.rm = TRUE)
  n_small <- sum(omega_vals >= 0.01 & omega_vals < 0.06, na.rm = TRUE)
  n_negl <- sum(omega_vals < 0.01, na.rm = TRUE)
  cat(sprintf("  Effect sizes: %d Large, %d Medium, %d Small, %d Negligible\n",
              n_large, n_medium, n_small, n_negl))

  invisible(x)
}
