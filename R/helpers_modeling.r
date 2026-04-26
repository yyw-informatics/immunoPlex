# Preprocessing helpers for immunoPlex modeling
# Covariate selection and outlier winsorization utilities

#' Conditionally include covariates based on prevalence
#'
#' Checks each candidate covariate's prevalence (proportion coded as 1) in
#' the data and includes it only if it meets the threshold. Useful for
#' conditionally adding rare binary covariates (e.g., HIV, malaria) to
#' a model formula.
#'
#' @param data A data frame containing the candidate covariate columns.
#' @param base_covars Character vector of covariate names to always include.
#' @param candidates Character vector of covariate names to conditionally
#'   include based on prevalence. Each must be a column in \code{data} with
#'   binary (0/1) coding. Defaults to \code{NULL} (no conditional covariates).
#' @param rare_threshold Numeric prevalence threshold (proportion of 1s).
#'   Candidates with \code{mean(x == 1) >= rare_threshold} are included.
#'   Default is 0.03 (3%).
#'
#' @return Character vector of covariates to include in the model formula.
#'
#' @examples
#' \dontrun{
#' covars <- choose_covariates(
#'   data = my_data,
#'   base_covars = c("age", "gest_acute"),
#'   candidates = c("HIV", "malaria"),
#'   rare_threshold = 0.03
#' )
#' # Returns e.g. c("age", "gest_acute", "HIV") if HIV >= 3% but malaria < 3%
#' }
#'
#' @export
choose_covariates <- function(data,
                              base_covars,
                              candidates = NULL,
                              rare_threshold = 0.03) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  stopifnot(is.character(base_covars))
  stopifnot(is.numeric(rare_threshold), length(rare_threshold) == 1L,
            rare_threshold >= 0, rare_threshold <= 1)

  covars <- base_covars

  if (!is.null(candidates)) {
    stopifnot(is.character(candidates))
    for (cand in candidates) {
      if (cand %in% names(data) &&
          mean(data[[cand]] == 1, na.rm = TRUE) >= rare_threshold) {
        covars <- c(covars, cand)
      }
    }
  }

  covars
}


#' Winsorize outliers within groups using IQR bounds
#'
#' For each combination of grouping columns, computes the IQR and caps
#' values beyond \code{iqr_mult * IQR} from Q1/Q3. This is a robust
#' approach to handling extreme values in repeated-measures cytokine data
#' where outliers may be group-specific.
#'
#' @param data A data frame containing the value and grouping columns.
#' @param value_col Character name of the column to winsorize. Default
#'   \code{"value"}.
#' @param group_cols Character vector of column names defining groups.
#'   Default \code{c("cytokine", "timepoint")}.
#' @param iqr_mult Numeric multiplier for IQR to set bounds. Default 3.
#' @param verbose Logical; if \code{TRUE} (default), prints outlier count.
#'
#' @return The data frame with outlier values in \code{value_col} replaced
#'   by group-specific IQR bounds.
#'
#' @examples
#' \dontrun{
#' cleaned <- winsorize_by_group(
#'   data = my_data,
#'   value_col = "value",
#'   group_cols = c("cytokine", "timepoint"),
#'   iqr_mult = 3
#' )
#' }
#'
#' @export
winsorize_by_group <- function(data,
                               value_col = "value",
                               group_cols = c("cytokine", "timepoint"),
                               iqr_mult = 3,
                               verbose = TRUE) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  stopifnot(value_col %in% names(data))
  stopifnot(all(group_cols %in% names(data)))
  stopifnot(is.numeric(iqr_mult), length(iqr_mult) == 1L, iqr_mult > 0)

  # Compute group-level IQR bounds
  outlier_summary <- data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      q1 = stats::quantile(.data[[value_col]], 0.25, na.rm = TRUE),
      q3 = stats::quantile(.data[[value_col]], 0.75, na.rm = TRUE),
      iqr = q3 - q1,
      lower_bound = q1 - iqr_mult * iqr,
      upper_bound = q3 + iqr_mult * iqr,
      .groups = "drop"
    )

  # Count outliers
  n_outliers <- data %>%
    dplyr::left_join(outlier_summary, by = group_cols) %>%
    dplyr::summarise(
      n = sum(.data[[value_col]] < lower_bound |
                .data[[value_col]] > upper_bound,
              na.rm = TRUE)
    ) %>%
    dplyr::pull(n)

  if (verbose) {
    if (n_outliers > 0) {
      message(sprintf("Winsorizing %d outlier values (%.1fx IQR rule)",
                      n_outliers, iqr_mult))
    } else {
      message("No outliers detected for winsorization")
    }
  }

  # Apply winsorization
  data %>%
    dplyr::left_join(outlier_summary, by = group_cols) %>%
    dplyr::mutate(
      "{value_col}" := dplyr::case_when(
        .data[[value_col]] > upper_bound ~ upper_bound,
        .data[[value_col]] < lower_bound ~ lower_bound,
        TRUE ~ .data[[value_col]]
      )
    ) %>%
    dplyr::select(-q1, -q3, -iqr, -lower_bound, -upper_bound)
}
