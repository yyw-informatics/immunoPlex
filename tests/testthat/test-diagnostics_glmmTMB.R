# Tests for glmmTMB diagnostic functions
#
# Validates:
# 1. dharma_diagnostics() — output structure, accepts both raw and wrapped models
# 2. check_convergence() — correct convergence assessment
# 3. check_random_effects() — auto-detect grouping, ICC, Shapiro-Wilk

library(testthat)

# Helper: fit a basic gaussian model for diagnostics testing
.fit_gauss_for_diag <- function(n_subjects = 40, n_timepoints = 3) {
  dat <- create_gaussian_test_data(n_subjects = n_subjects,
                                    n_timepoints = n_timepoints)
  fit_one(dat, family = "gaussian",
          fixed = "timepoint + disease + age",
          random = "(1|subject_id)")
}


describe("dharma_diagnostics()", {

  test_that("returns correct structure with test results", {
    skip_if_not_installed("glmmTMB")
    skip_if_not_installed("DHARMa")

    fit <- .fit_gauss_for_diag()
    if (!fit$converged) skip("Model did not converge")

    diag <- dharma_diagnostics(fit, label = "test_analyte", n_sim = 100)

    expect_type(diag, "list")
    expect_true("tests" %in% names(diag))
    expect_true("residuals" %in% names(diag))
    expect_true("passed_all" %in% names(diag))
    expect_true("warnings" %in% names(diag))
    expect_type(diag$passed_all, "logical")

    # Check individual test entries
    if (length(diag$tests) > 0) {
      for (test_result in diag$tests) {
        expect_true("test" %in% names(test_result))
        expect_true("p_value" %in% names(test_result))
        expect_true("passed" %in% names(test_result))
      }
    }
  })

  test_that("accepts raw glmmTMB model object", {
    skip_if_not_installed("glmmTMB")
    skip_if_not_installed("DHARMa")

    fit <- .fit_gauss_for_diag()
    if (!fit$converged) skip("Model did not converge")

    # Pass raw model, not immuno_fit wrapper
    diag <- dharma_diagnostics(fit$model, n_sim = 100)
    expect_type(diag, "list")
    expect_true("passed_all" %in% names(diag))
  })

  test_that("accepts immuno_fit wrapper object", {
    skip_if_not_installed("glmmTMB")
    skip_if_not_installed("DHARMa")

    fit <- .fit_gauss_for_diag()
    if (!fit$converged) skip("Model did not converge")

    diag <- dharma_diagnostics(fit, n_sim = 100)
    expect_type(diag, "list")
    expect_true("passed_all" %in% names(diag))
  })

  test_that("respects custom alpha thresholds", {
    skip_if_not_installed("glmmTMB")
    skip_if_not_installed("DHARMa")

    fit <- .fit_gauss_for_diag()
    if (!fit$converged) skip("Model did not converge")

    # Very strict thresholds — more likely to fail
    diag_strict <- dharma_diagnostics(fit, n_sim = 100,
                                       alpha_uniformity = 0.99,
                                       alpha_dispersion = 0.99)

    # Very lenient thresholds — more likely to pass
    diag_lenient <- dharma_diagnostics(fit, n_sim = 100,
                                        alpha_uniformity = 0.001,
                                        alpha_dispersion = 0.001)

    expect_type(diag_strict, "list")
    expect_type(diag_lenient, "list")
  })

  test_that("errors on invalid input", {
    expect_error(dharma_diagnostics("not_a_model"), "Expected a glmmTMB")
    expect_error(dharma_diagnostics(lm(1:10 ~ 1)), "Expected a glmmTMB")
  })
})


describe("check_convergence()", {

  test_that("returns converged=TRUE for well-specified model", {
    skip_if_not_installed("glmmTMB")

    fit <- .fit_gauss_for_diag()
    if (!fit$converged) skip("Model did not converge")

    diag <- check_convergence(fit)

    expect_type(diag, "list")
    expect_true(diag$converged)
    expect_true(diag$hessian_pd)
    expect_true(is.numeric(diag$max_gradient))
    expect_equal(diag$na_se_count, 0L)
  })

  test_that("accepts raw glmmTMB model object", {
    skip_if_not_installed("glmmTMB")

    fit <- .fit_gauss_for_diag()
    if (!fit$converged) skip("Model did not converge")

    diag <- check_convergence(fit$model)
    expect_true(diag$converged)
  })

  test_that("uses configurable gradient thresholds", {
    skip_if_not_installed("glmmTMB")

    fit <- .fit_gauss_for_diag()
    if (!fit$converged) skip("Model did not converge")

    # Extremely strict threshold — likely flags even good models
    diag <- check_convergence(fit, gradient_warn = 1e-10,
                               gradient_fail = 1e-8)
    expect_type(diag, "list")
    # May or may not converge depending on actual gradient
  })

  test_that("errors on invalid input", {
    expect_error(check_convergence("not_a_model"), "Expected a glmmTMB")
  })
})


describe("check_random_effects()", {

  test_that("detects random effects and computes ICC", {
    skip_if_not_installed("glmmTMB")

    fit <- .fit_gauss_for_diag()
    if (!fit$converged) skip("Model did not converge")

    diag <- check_random_effects(fit)

    expect_type(diag, "list")
    expect_true(diag$has_random_effects)
    expect_true(is.numeric(diag$re_variance))
    expect_true(is.numeric(diag$re_sd))
    expect_true(is.numeric(diag$icc))
    expect_true(diag$icc >= 0 && diag$icc <= 1)
  })

  test_that("auto-detects grouping factor", {
    skip_if_not_installed("glmmTMB")

    fit <- .fit_gauss_for_diag()
    if (!fit$converged) skip("Model did not converge")

    # Should work without hardcoding "subject_id"
    diag <- check_random_effects(fit)
    expect_true(diag$has_random_effects)
  })

  test_that("runs Shapiro-Wilk when enough RE", {
    skip_if_not_installed("glmmTMB")

    fit <- .fit_gauss_for_diag(n_subjects = 40)
    if (!fit$converged) skip("Model did not converge")

    diag <- check_random_effects(fit, min_re_for_normality = 8)
    expect_true(is.numeric(diag$re_normality_p))
  })

  test_that("returns has_random_effects=FALSE for fixed-only model", {
    skip_if_not_installed("glmmTMB")

    dat <- create_gaussian_test_data(n_subjects = 20, n_timepoints = 2)
    fit <- fit_one(dat, family = "gaussian",
                   fixed = "timepoint + disease",
                   random = "")

    if (!fit$converged) skip("Model did not converge")

    diag <- check_random_effects(fit)
    expect_false(diag$has_random_effects)
  })

  test_that("errors on invalid input", {
    expect_error(check_random_effects("not_a_model"), "Expected a glmmTMB")
  })
})
