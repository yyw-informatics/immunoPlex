# Leave-one-out sensitivity analysis for mixed models
# Identifies influential subjects whose removal substantially changes estimates

#' Leave-one-out sensitivity analysis
#'
#' For each unique subject in the data, refits the model with that subject
#' excluded and records the change in coefficient estimates for the target
#' terms. Subjects whose removal changes an estimate by more than
#' \code{influence_threshold} percent are flagged as influential.
#'
#' Uses \code{\link{fit_one}} internally, inheriting the optimizer cascade
#' and convergence checking.
#'
#' @param data A data frame containing all observations.
#' @param fixed Character string of fixed effects formula.
#' @param target_terms Character vector of coefficient names to assess
#'   (e.g. \code{c("timepointAcute", "timepointConvalescent")}).
#' @param random Character string of random effects. Default
#'   \code{"(1|subject_id)"}.
#' @param family Character family for \code{fit_one()}. Default
#'   \code{"gaussian"}.
#' @param dispformula Dispersion formula passed to \code{fit_one()}.
#'   Default \code{~1}.
#' @param subject_col Character name of the column defining the LOO
#'   grouping (i.e. which rows to exclude together). Default
#'   \code{"subject_id"}.
#' @param influence_threshold Numeric percent change threshold to flag
#'   a subject as influential. Default 20 (i.e. 20%).
#' @param progress Logical; if \code{TRUE} (default), prints progress
#'   messages.
#'
#' @return A tibble with columns: \code{excluded_subject}, \code{term},
#'   \code{full_estimate}, \code{loo_estimate}, \code{difference},
#'   \code{pct_change}, \code{influential}.
#'
#' @export
leave_one_out_sensitivity <- function(data,
                            fixed,
                            target_terms,
                            random = "(1|subject_id)",
                            family = "gaussian",
                            dispformula = ~1,
                            subject_col = NULL,
                            influence_threshold = NULL,
                            progress = TRUE) {

  if (!is.data.frame(data)) stop("'data' must be a data frame")
  stopifnot(is.character(fixed), is.character(target_terms))

  if (is.null(subject_col)) {
    subject_col <- getOption("immunoplex.glmmTMB_loo_subject_col", "subject_id")
  }
  if (is.null(influence_threshold)) {
    influence_threshold <- getOption("immunoplex.glmmTMB_loo_threshold", 20)
  }

  if (!subject_col %in% names(data)) {
    stop("Column '", subject_col, "' not found in data")
  }

  subjects <- unique(data[[subject_col]])

  # Fit full model
  full_fit <- tryCatch(
    fit_one(data, family = family, fixed = fixed, random = random,
            dispformula = dispformula),
    error = function(e) NULL
  )

  if (is.null(full_fit) || !isTRUE(full_fit$converged) ||
      !inherits(full_fit$model, "glmmTMB")) {
    warning("Full model did not converge; cannot perform LOO sensitivity")
    return(NULL)
  }

  full_coef <- glmmTMB::fixef(full_fit$model)$cond

  # Check that target_terms exist in model coefficients
  missing_terms <- setdiff(target_terms, names(full_coef))
  if (length(missing_terms) > 0) {
    warning("Terms not found in model coefficients: ",
            paste(missing_terms, collapse = ", "))
    target_terms <- intersect(target_terms, names(full_coef))
    if (length(target_terms) == 0L) return(NULL)
  }

  # LOO for each subject
  loo_results <- vector("list", length(subjects) * length(target_terms))
  idx <- 0L
  n_subj <- length(subjects)

  for (i in seq_along(subjects)) {
    subj <- subjects[i]
    if (progress && i %% 10 == 0) {
      message(sprintf("  LOO: %d/%d subjects", i, n_subj))
    }

    loo_data <- data[data[[subject_col]] != subj, , drop = FALSE]

    loo_fit <- tryCatch(
      suppressMessages(
        fit_one(loo_data, family = family, fixed = fixed, random = random,
                dispformula = dispformula)
      ),
      error = function(e) NULL
    )

    if (is.null(loo_fit) || !isTRUE(loo_fit$converged) ||
        !inherits(loo_fit$model, "glmmTMB")) {
      next
    }

    loo_coef <- glmmTMB::fixef(loo_fit$model)$cond

    for (term in target_terms) {
      if (term %in% names(loo_coef)) {
        full_est <- full_coef[term]
        loo_est <- loo_coef[term]
        pct_change <- if (abs(full_est) > 0.001) {
          100 * (loo_est - full_est) / abs(full_est)
        } else {
          NA_real_
        }

        idx <- idx + 1L
        loo_results[[idx]] <- tibble::tibble(
          excluded_subject = as.character(subj),
          term = term,
          full_estimate = full_est,
          loo_estimate = loo_est,
          difference = loo_est - full_est,
          pct_change = pct_change,
          influential = !is.na(pct_change) & abs(pct_change) > influence_threshold
        )
      }
    }
  }

  # Remove NULL entries and bind
  loo_results <- loo_results[seq_len(idx)]
  if (length(loo_results) == 0L) return(NULL)

  dplyr::bind_rows(loo_results)
}
