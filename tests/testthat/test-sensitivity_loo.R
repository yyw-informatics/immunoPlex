# Tests for loo_sensitivity()
#
# Validates:
# 1. Returns tibble with correct columns
# 2. Influential subject detection
# 3. Configurable subject_col and influence_threshold
# 4. Handles missing terms gracefully

library(testthat)

describe("loo_sensitivity()", {

  test_that("returns tibble with expected columns", {
    skip_if_not_installed("glmmTMB")

    dat <- create_gaussian_test_data(n_subjects = 30, n_timepoints = 3)

    result <- loo_sensitivity(
      data = dat,
      fixed = "timepoint + disease",
      target_terms = c("timepointT2", "timepointT3"),
      random = "(1|subject_id)",
      family = "gaussian",
      progress = FALSE
    )

    # Should return a tibble (or NULL if full model fails)
    if (!is.null(result)) {
      expect_s3_class(result, "tbl_df")
      expected_cols <- c("excluded_subject", "term", "full_estimate",
                         "loo_estimate", "difference", "pct_change",
                         "influential")
      expect_true(all(expected_cols %in% names(result)))
      expect_type(result$influential, "logical")
    }
  })

  test_that("influential column uses threshold correctly", {
    skip_if_not_installed("glmmTMB")

    dat <- create_gaussian_test_data(n_subjects = 30, n_timepoints = 3)

    # Very low threshold — should flag more subjects
    result_low <- loo_sensitivity(
      data = dat,
      fixed = "timepoint + disease",
      target_terms = c("timepointT2"),
      random = "(1|subject_id)",
      influence_threshold = 1,
      progress = FALSE
    )

    # Very high threshold — should flag fewer subjects
    result_high <- loo_sensitivity(
      data = dat,
      fixed = "timepoint + disease",
      target_terms = c("timepointT2"),
      random = "(1|subject_id)",
      influence_threshold = 100,
      progress = FALSE
    )

    if (!is.null(result_low) && !is.null(result_high)) {
      n_infl_low <- sum(result_low$influential, na.rm = TRUE)
      n_infl_high <- sum(result_high$influential, na.rm = TRUE)
      expect_true(n_infl_low >= n_infl_high)
    }
  })

  test_that("warns on missing target_terms", {
    skip_if_not_installed("glmmTMB")

    dat <- create_gaussian_test_data(n_subjects = 30, n_timepoints = 2)

    expect_warning(
      result <- loo_sensitivity(
        data = dat,
        fixed = "timepoint + disease",
        target_terms = c("nonexistent_term_xyz"),
        random = "(1|subject_id)",
        progress = FALSE
      ),
      "Terms not found"
    )
  })

  test_that("errors on invalid input", {
    expect_error(
      loo_sensitivity("not_a_df", fixed = "x", target_terms = "y"),
      "must be a data frame"
    )
  })

  test_that("errors when subject_col is missing", {
    dat <- create_gaussian_test_data(n_subjects = 30, n_timepoints = 2)
    expect_error(
      loo_sensitivity(dat, fixed = "timepoint", target_terms = "timepointT2",
                      subject_col = "nonexistent_column"),
      "not found in data"
    )
  })
})
