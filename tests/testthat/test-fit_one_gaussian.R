# Tests for Gaussian family in fit_one()
#
# Validates:
# 1. Basic gaussian fit returns immuno_fit with correct family/estimand
# 2. Optimizer cascade convergence
# 3. dispformula parameter (homo vs hetero)
# 4. Random effects dropping for small samples
# 5. Auto family selection picks gaussian when no censoring
# 6. Integration with fit_models()
# 7. Error handling

library(testthat)

describe("fit_one() gaussian family", {

  test_that("basic gaussian fit returns correct immuno_fit structure", {
    skip_if_not_installed("glmmTMB")

    dat <- create_gaussian_test_data(n_subjects = 30, n_timepoints = 3)

    fit <- fit_one(dat, family = "gaussian",
                   fixed = "timepoint + disease + age",
                   random = "(1|subject_id)")

    expect_s3_class(fit, "immuno_fit")
    expect_equal(fit$family, "gaussian")
    expect_equal(fit$estimand, "mean_difference")
    expect_true(fit$converged)
    expect_true(is.numeric(fit$aic))
    expect_true(is.numeric(fit$bic))
    expect_true(is.numeric(fit$logLik))
    expect_true(inherits(fit$model, "glmmTMB"))
  })

  test_that("gaussian converges with default optimizer (nlminb)", {
    skip_if_not_installed("glmmTMB")

    dat <- create_gaussian_test_data(n_subjects = 40, n_timepoints = 3)

    fit <- fit_one(dat, family = "gaussian",
                   fixed = "timepoint + disease",
                   random = "(1|subject_id)")

    expect_true(fit$converged)
    expect_true(isTRUE(fit$model$sdr$pdHess))
  })

  test_that("dispformula = ~timepoint produces different fit than ~1", {
    skip_if_not_installed("glmmTMB")

    dat <- create_gaussian_test_data(n_subjects = 40, n_timepoints = 3)

    fit_homo <- fit_one(dat, family = "gaussian",
                        fixed = "timepoint + disease",
                        random = "(1|subject_id)",
                        dispformula = ~1)

    fit_hetero <- fit_one(dat, family = "gaussian",
                          fixed = "timepoint + disease",
                          random = "(1|subject_id)",
                          dispformula = ~timepoint)

    expect_true(fit_homo$converged)
    expect_true(fit_hetero$converged)
    # Models should differ (different dispersion structure)
    expect_false(identical(fit_homo$aic, fit_hetero$aic))
  })

  test_that("random effects dropped for small samples", {
    skip_if_not_installed("glmmTMB")

    dat <- create_small_gaussian_data(n_subjects = 10)

    # With default min_subjects=30, should drop RE
    expect_message(
      fit <- fit_one(dat, family = "gaussian",
                     fixed = "timepoint + disease",
                     random = "(1|subject_id)"),
      "Dropping random term"
    )

    expect_s3_class(fit, "immuno_fit")
    expect_equal(fit$family, "gaussian")
  })

  test_that("gaussian without random effects works", {
    skip_if_not_installed("glmmTMB")

    dat <- create_gaussian_test_data(n_subjects = 20, n_timepoints = 2)

    fit <- fit_one(dat, family = "gaussian",
                   fixed = "timepoint + disease",
                   random = "")

    expect_s3_class(fit, "immuno_fit")
    expect_equal(fit$family, "gaussian")
    expect_true(fit$converged)
  })

  test_that("auto family selects gaussian when no censoring", {
    skip_if_not_installed("glmmTMB")
    skip_if_not_installed("moments")

    dat <- create_gaussian_test_data(n_subjects = 30, n_timepoints = 3)
    # All cens_lod and cens_ulod are FALSE

    fit <- fit_one(dat, family = "auto",
                   fixed = "timepoint + disease",
                   random = "(1|subject_id)")

    expect_equal(fit$family, "gaussian")
    expect_equal(fit$estimand, "mean_difference")
  })

  test_that("gaussian fit works with fit_models()", {
    skip_if_not_installed("glmmTMB")

    dat <- create_gaussian_test_data(n_subjects = 30, n_timepoints = 3)

    models <- fit_models(dat,
                         families = c("gaussian", "gamma"),
                         fixed = "timepoint + disease + age",
                         random = "(1|subject_id)",
                         quiet = TRUE)

    expect_s3_class(models, "immuno_model_set")
    expect_true("gaussian" %in% models$comparison$family)
    # Check estimand column
    gauss_row <- models$comparison[models$comparison$family == "gaussian", ]
    expect_equal(gauss_row$estimand, "mean_difference")
  })

  test_that("gaussian fit works with compute_residual_diagnostics()", {
    skip_if_not_installed("glmmTMB")
    skip_if_not_installed("DHARMa")

    dat <- create_gaussian_test_data(n_subjects = 30, n_timepoints = 3)

    fit <- fit_one(dat, family = "gaussian",
                   fixed = "timepoint + disease",
                   random = "(1|subject_id)")

    diag <- compute_residual_diagnostics(fit, nsim = 100)
    expect_type(diag, "list")
    expect_true("p_uniform" %in% names(diag))
    expect_true("dispersion_proxy" %in% names(diag))
  })

  test_that("invalid family still errors", {
    dat <- create_gaussian_test_data()
    expect_error(fit_one(dat, family = "invalid"), "should be one of")
  })

  test_that("gaussian estimand is mean_difference", {
    skip_if_not_installed("glmmTMB")

    dat <- create_gaussian_test_data(n_subjects = 30, n_timepoints = 2)

    fit <- fit_one(dat, family = "gaussian",
                   fixed = "timepoint + disease",
                   random = "(1|subject_id)")

    expect_equal(fit$estimand, "mean_difference")
  })
})
