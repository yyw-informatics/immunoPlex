# Batch ANCOVA wrapper for immunoPlex
# Iterates ancova_one() over analytes with FDR correction

#' Fit Rank-Based ANCOVA Across Multiple Analytes
#'
#' Iterates \code{\link{ancova_one}} over analytes in a long-format data frame,
#' applies FDR correction, and returns a structured result set.
#'
#' @importFrom stats p.adjust
#'
#' @param data A data frame in long format with one row per subject-analyte
#'   combination, containing outcome, group, covariate, and analyte columns.
#' @param outcome Character. Column name for the outcome variable.
#' @param group Character. Column name for the group factor.
#' @param covariate Character. Column name for the primary covariate.
#' @param analyte_col Character. Column name identifying the analyte/cytokine.
#'   Default \code{"cytokine"}.
#' @param covariates Character vector. Additional fixed covariates. Default \code{NULL}.
#' @param conditional_covariates Character vector. Binary covariates included
#'   conditionally based on prevalence. Default \code{NULL}.
#' @param prevalence_threshold Numeric. Minimum prevalence for conditional
#'   covariates. Default from config (0.03).
#' @param outlier_removal Logical. Apply IQR outlier removal. Default \code{TRUE}.
#' @param iqr_multiplier Numeric. IQR multiplier for outlier fencing. Default
#'   from config (3).
#' @param min_n Integer. Minimum observations per analyte. Default from config (12).
#' @param conf_level Numeric. Confidence level. Default 0.95.
#' @param fdr_method Character. FDR correction method passed to \code{p.adjust()}.
#'   Default from config (\code{"BH"}).
#' @param fdr_threshold Numeric. Significance threshold for FDR-corrected
#'   p-values. Default from config (0.05).
#' @param log_validate Logical. Also fit log-scale validation models.
#'   Default \code{FALSE}.
#' @param rep_col Optional character. Column name identifying technical
#'   replicates within a subject. Forwarded to \code{\link{ancova_one}};
#'   when supplied, each per-analyte fit switches to \code{glmmTMB::glmmTMB()}
#'   with a nested random intercept \code{(1 | subject_id:rep_col)}.
#'   Default \code{NULL}. See \code{\link{aggregate_replicates}} for a
#'   pre-aggregation alternative.
#' @param plate_col Optional character. Column name identifying plates.
#'   Forwarded to \code{\link{ancova_one}}; when supplied, each per-analyte
#'   fit appends a flat random intercept \code{(1 | plate_col)}.
#'   Default \code{NULL}.
#' @param quiet Logical. Suppress progress messages. Default \code{FALSE}.
#'
#' @return An S3 object of class \code{"immuno_ancova_set"}, a list containing:
#'   \describe{
#'     \item{results}{A \code{tibble} with one row per analyte and columns for
#'       all statistics, plus \code{q_value} (FDR-corrected p-value) columns}
#'     \item{models}{Named list of \code{immuno_ancova} objects, one per analyte}
#'     \item{n_analytes}{Integer: total analytes analyzed}
#'     \item{n_significant}{Integer: count with q-value below \code{fdr_threshold}}
#'     \item{fdr_method}{Character: FDR method used}
#'     \item{fdr_threshold}{Numeric: significance threshold}
#'     \item{call}{The matched call}
#'     \item{config}{List: config snapshot for reproducibility}
#'   }
#'
#' @examples
#' \dontrun{
#' results <- ancova_fit(
#'   data = my_long_data,
#'   outcome = "cord_value",
#'   group = "covid_status",
#'   covariate = "maternal_value",
#'   analyte_col = "cytokine",
#'   covariates = c("age", "gest_delivery"),
#'   conditional_covariates = c("HIV", "malaria"),
#'   log_validate = TRUE
#' )
#' print(results)
#' summary(results)
#' }
#'
#' @export
ancova_fit <- function(data,
                       outcome,
                       group,
                       covariate,
                       analyte_col = "cytokine",
                       covariates = NULL,
                       conditional_covariates = NULL,
                       prevalence_threshold = getOption("immunoplex.ancova_prevalence_threshold", 0.03),
                       outlier_removal = TRUE,
                       iqr_multiplier = getOption("immunoplex.ancova_iqr_multiplier", 3),
                       min_n = getOption("immunoplex.ancova_min_n", 12L),
                       conf_level = 0.95,
                       fdr_method = getOption("immunoplex.ancova_fdr_method", "BH"),
                       fdr_threshold = getOption("immunoplex.ancova_fdr_threshold", 0.05),
                       log_validate = FALSE,
                       rep_col = NULL,
                       plate_col = NULL,
                       quiet = FALSE) {

  matched_call <- match.call()

  # --- Input validation ---
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame", call. = FALSE)
  }
  if (!analyte_col %in% names(data)) {
    stop("Analyte column '", analyte_col, "' not found in data", call. = FALSE)
  }

  analytes <- unique(data[[analyte_col]])
  n_analytes <- length(analytes)

  if (!quiet) {
    message(sprintf("ancova_fit: analyzing %d analytes", n_analytes))
  }

  # --- Iterate over analytes ---
  models <- setNames(vector("list", n_analytes), analytes)
  results_list <- vector("list", n_analytes)

  for (i in seq_along(analytes)) {
    a <- analytes[i]
    if (!quiet) {
      message(sprintf("  [%d/%d] %s", i, n_analytes, a))
    }

    # Subset to this analyte
    d <- data[data[[analyte_col]] == a, , drop = FALSE]

    # Fit single-analyte ANCOVA
    fit <- ancova_one(
      data = d,
      outcome = outcome,
      group = group,
      covariate = covariate,
      covariates = covariates,
      conditional_covariates = conditional_covariates,
      prevalence_threshold = prevalence_threshold,
      outlier_removal = outlier_removal,
      iqr_multiplier = iqr_multiplier,
      min_n = min_n,
      conf_level = conf_level,
      log_validate = log_validate,
      rep_col = rep_col,
      plate_col = plate_col
    )

    # Tag with analyte name
    fit$analyte <- a
    models[[a]] <- fit

    # Build summary row
    row <- data.frame(
      analyte = a,
      n_obs = fit$n_obs,
      n_removed = fit$n_removed,
      r_squared = fit$r_squared,
      omega_sq_partial = fit$omega_sq_partial,
      omega_sq_ci_low = fit$omega_sq_ci[1],
      omega_sq_ci_high = fit$omega_sq_ci[2],
      omega_magnitude = fit$omega_magnitude,
      stringsAsFactors = FALSE
    )

    if (fit$n_groups == 2) {
      row$group_effect <- fit$group_effect
      row$group_se <- fit$group_se
      row$group_p <- fit$group_p
      row$group_ci_lower <- fit$group_ci[1]
      row$group_ci_upper <- fit$group_ci[2]
      if (log_validate) {
        row$log_effect <- if (!is.null(fit$log_effect)) fit$log_effect else NA_real_
        row$log_fold_change <- if (!is.null(fit$log_fold_change)) fit$log_fold_change else NA_real_
        row$log_p <- if (!is.null(fit$log_p)) fit$log_p else NA_real_
      }
    } else {
      row$omnibus_f <- fit$omnibus_f
      row$omnibus_p <- fit$omnibus_p
    }

    row$n_warnings <- length(fit$warnings)
    results_list[[i]] <- row
  }

  # --- Combine results ---
  results <- do.call(rbind, results_list)

  # --- FDR correction ---
  if (n_analytes > 0) {
    first_fit <- models[[1]]

    if (first_fit$n_groups == 2) {
      results$q_value <- p.adjust(results$group_p, method = fdr_method)
      if (log_validate && "log_p" %in% names(results)) {
        results$q_log <- p.adjust(results$log_p, method = fdr_method)
      }
      n_significant <- sum(results$q_value < fdr_threshold, na.rm = TRUE)
    } else {
      # k-level: FDR on omnibus + each pairwise contrast
      results$q_omnibus <- p.adjust(results$omnibus_p, method = fdr_method)

      # Also add FDR-corrected pairwise contrasts
      contrast_terms <- if (!is.null(first_fit$contrasts)) first_fit$contrasts$term else character(0)
      for (ct in contrast_terms) {
        col_p <- paste0("p_", gsub(" ", "_", ct))
        col_q <- paste0("q_", gsub(" ", "_", ct))
        # Extract p-values for this contrast across all analytes
        p_vals <- vapply(models, function(m) {
          if (!is.null(m$contrasts)) {
            idx <- which(m$contrasts$term == ct)
            if (length(idx) > 0) return(m$contrasts$p[idx])
          }
          NA_real_
        }, numeric(1))
        results[[col_p]] <- p_vals
        results[[col_q]] <- p.adjust(p_vals, method = fdr_method)
      }

      n_significant <- sum(results$q_omnibus < fdr_threshold, na.rm = TRUE)
    }
  } else {
    n_significant <- 0L
  }

  # Convert to tibble if available
  if (requireNamespace("tibble", quietly = TRUE)) {
    results <- tibble::as_tibble(results)
  }

  if (!quiet) {
    message(sprintf("ancova_fit: %d/%d significant (FDR %s < %.2f)",
                    n_significant, n_analytes, fdr_method, fdr_threshold))
  }

  # --- Construct result ---
  structure(list(
    results = results,
    models = models,
    n_analytes = n_analytes,
    n_significant = n_significant,
    fdr_method = fdr_method,
    fdr_threshold = fdr_threshold,
    call = matched_call,
    config = list(
      outcome = outcome,
      group = group,
      covariate = covariate,
      analyte_col = analyte_col,
      covariates = covariates,
      conditional_covariates = conditional_covariates,
      prevalence_threshold = prevalence_threshold,
      outlier_removal = outlier_removal,
      iqr_multiplier = iqr_multiplier,
      min_n = min_n,
      conf_level = conf_level,
      log_validate = log_validate
    )
  ), class = "immuno_ancova_set")
}
