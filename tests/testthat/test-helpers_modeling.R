# Tests for choose_covariates() and winsorize_by_group()
#
# Validates:
# 1. choose_covariates() includes/excludes based on prevalence
# 2. winsorize_by_group() caps outliers correctly
# 3. Edge cases and error handling

library(testthat)

describe("choose_covariates()", {

  test_that("returns base_covars when no candidates", {
    dat <- data.frame(age = 1:100, x = 1:100)
    result <- choose_covariates(dat, base_covars = c("age"), candidates = NULL)
    expect_equal(result, c("age"))
  })

  test_that("includes candidates above threshold", {
    dat <- data.frame(
      age = 1:100,
      HIV = c(rep(1, 10), rep(0, 90)),      # 10% prevalence
      malaria = c(rep(1, 2), rep(0, 98))     # 2% prevalence
    )

    result <- choose_covariates(
      dat,
      base_covars = c("age"),
      candidates = c("HIV", "malaria"),
      rare_threshold = 0.03
    )

    expect_true("HIV" %in% result)       # 10% >= 3%
    expect_false("malaria" %in% result)  # 2% < 3%
    expect_true("age" %in% result)
  })

  test_that("handles missing candidates gracefully", {
    dat <- data.frame(age = 1:100, HIV = rep(1, 100))

    result <- choose_covariates(
      dat,
      base_covars = c("age"),
      candidates = c("HIV", "nonexistent_var"),
      rare_threshold = 0.03
    )

    expect_true("HIV" %in% result)
    expect_false("nonexistent_var" %in% result)
  })

  test_that("threshold = 0 includes all candidates", {
    dat <- data.frame(
      age = 1:100,
      HIV = c(1, rep(0, 99))  # 1% prevalence
    )

    result <- choose_covariates(
      dat,
      base_covars = c("age"),
      candidates = c("HIV"),
      rare_threshold = 0
    )

    expect_true("HIV" %in% result)
  })

  test_that("errors on invalid data", {
    expect_error(choose_covariates("not_a_df", base_covars = "x"),
                 "must be a data frame")
  })
})


describe("winsorize_by_group()", {

  test_that("caps outliers to IQR bounds", {
    dat <- data.frame(
      cytokine = rep("A", 20),
      timepoint = rep("T1", 20),
      value = c(rnorm(18, 5, 1), 50, -40)  # Two extreme outliers
    )

    result <- suppressMessages(
      winsorize_by_group(dat, value_col = "value",
                         group_cols = c("cytokine", "timepoint"),
                         iqr_mult = 3)
    )

    expect_true(max(result$value) < 50)
    expect_true(min(result$value) > -40)
  })

  test_that("no change when no outliers", {
    set.seed(42)
    dat <- data.frame(
      cytokine = rep("A", 50),
      timepoint = rep("T1", 50),
      value = rnorm(50, 5, 0.5)
    )

    result <- suppressMessages(
      winsorize_by_group(dat, iqr_mult = 3)
    )

    expect_equal(unname(result$value), unname(dat$value))
  })

  test_that("respects group_cols", {
    set.seed(42)
    dat <- data.frame(
      cytokine = rep(c("A", "B"), each = 20),
      timepoint = rep("T1", 40),
      value = c(rnorm(19, 5, 1), 100,    # outlier in A
                rnorm(20, 50, 1))          # no outlier in B
    )

    result <- suppressMessages(
      winsorize_by_group(dat, iqr_mult = 3)
    )

    # Group B values should be unchanged
    b_original <- dat$value[dat$cytokine == "B"]
    b_result <- result$value[result$cytokine == "B"]
    expect_equal(unname(b_original), unname(b_result))

    # Group A's outlier should be capped
    expect_true(max(result$value[result$cytokine == "A"]) < 100)
  })

  test_that("verbose controls message output", {
    dat <- data.frame(
      cytokine = rep("A", 20),
      timepoint = rep("T1", 20),
      value = rnorm(20, 5, 1)
    )

    expect_message(
      winsorize_by_group(dat, verbose = TRUE),
      "outlier"
    )

    expect_silent(
      winsorize_by_group(dat, verbose = FALSE)
    )
  })

  test_that("errors on missing columns", {
    dat <- data.frame(x = 1:10)
    expect_error(
      winsorize_by_group(dat, value_col = "value"),
      "value_col"
    )
  })

  test_that("errors on invalid data", {
    expect_error(winsorize_by_group("not_a_df"), "must be a data frame")
  })
})
