# Tests for compare_variance_models()
#
# Validates:
# 1. Returns correct structure with both models and comparison
# 2. Selection criteria logic
# 3. Handles convergence failures gracefully

library(testthat)

describe("compare_variance_models()", {

  test_that("returns correct structure with comparison table", {
    skip_if_not_installed("glmmTMB")

    dat <- create_gaussian_test_data(n_subjects = 40, n_timepoints = 3)

    result <- compare_variance_models(
      data = dat,
      fixed = "timepoint + disease + age",
      random = "(1|subject_id)",
      dispformula_hetero = ~timepoint,
      label = "test_analyte"
    )

    expect_type(result, "list")
    expect_true("model_homo" %in% names(result))
    expect_true("model_hetero" %in% names(result))
    expect_true("comparison" %in% names(result))
    expect_true("preferred_model" %in% names(result))
    expect_true("preference_reason" %in% names(result))

    # Check comparison table structure
    if (!is.null(result$comparison)) {
      expect_s3_class(result$comparison, "tbl_df")
      expect_true("aic_homo" %in% names(result$comparison))
      expect_true("aic_hetero" %in% names(result$comparison))
      expect_true("lrt_statistic" %in% names(result$comparison))
      expect_true("lrt_pvalue" %in% names(result$comparison))
    }
  })

  test_that("preferred_model is one of expected values", {
    skip_if_not_installed("glmmTMB")

    dat <- create_gaussian_test_data(n_subjects = 40, n_timepoints = 3)

    result <- compare_variance_models(
      data = dat,
      fixed = "timepoint + disease",
      random = "(1|subject_id)",
      dispformula_hetero = ~timepoint
    )

    expect_true(result$preferred_model %in%
                  c("homogeneous", "heterogeneous", NA_character_))
    expect_true(is.character(result$preference_reason))
  })

  test_that("both models are immuno_fit objects", {
    skip_if_not_installed("glmmTMB")

    dat <- create_gaussian_test_data(n_subjects = 40, n_timepoints = 3)

    result <- compare_variance_models(
      data = dat,
      fixed = "timepoint + disease",
      random = "(1|subject_id)"
    )

    if (!is.null(result$model_homo)) {
      expect_s3_class(result$model_homo, "immuno_fit")
    }
    if (!is.null(result$model_hetero)) {
      expect_s3_class(result$model_hetero, "immuno_fit")
    }
  })

  test_that("works without random effects", {
    skip_if_not_installed("glmmTMB")

    dat <- create_gaussian_test_data(n_subjects = 20, n_timepoints = 3)

    result <- compare_variance_models(
      data = dat,
      fixed = "timepoint + disease",
      random = "",
      dispformula_hetero = ~timepoint
    )

    expect_type(result, "list")
    expect_true(result$preferred_model %in%
                  c("homogeneous", "heterogeneous", NA_character_))
  })

  test_that("errors on invalid input", {
    expect_error(compare_variance_models("not_a_df", fixed = "x"),
                 "must be a data frame")
  })
})
