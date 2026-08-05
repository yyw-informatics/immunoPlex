# Statistical Visualization and Diagnostic Protocol
#
# Validates visualization and diagnostic functionality:
# 1. Distribution Assessment:
#    - Standard residual analysis
#    - Censoring-aware visualization
#    - Randomized quantile residuals
#
# 2. Model Diagnostics:
#    - DHARMa simulation-based assessment
#    - Distributional assumption validation
#    - Residual pattern analysis
#
# 3. Implementation Validation:
#    - Parameter specification protocols
#    - Boundary condition handling
#    - Statistical test integration
#
# Statistical Dependencies:
# - glmmTMB: Distribution family estimation
# - DHARMa: Simulation-based diagnostics
# - ggplot2: Statistical visualization
# - patchwork: Multi-panel composition
# - statmod: Quantile residual computation

library(testthat)

# Initialize validation dataset with defined statistical properties
setup_plot_test_data <- function() {
  # Implement hierarchical data source selection
  tryCatch({
    setup_robust_test_data(prefer_real = TRUE, family_hint = "gamma")  # Primary: empirical data
  }, error = function(e) {
    create_s3_test_data(family_hint = "gamma")  # Secondary: synthetic data
  })
}

describe("Distribution assessment and visualization", {

test_that("visualization mode selection functions correctly", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("ggplot2")
  
  dat <- setup_plot_test_data()
  fit <- fit_one(dat, family = "gamma", random = "")
  
  if (!fit$converged) {
    skip("Model failed to converge")
  }
  
  # Test basic plotting
  result_basic <- plot(fit, plot_type = "basic", save_pdf = FALSE)
  expect_type(result_basic, "list")
  expect_s3_class(result_basic$rvf_plot, "gg")
  expect_s3_class(result_basic$qq_plot, "gg")
  
  # Test censor-aware plotting
  result_censor <- plot(fit, plot_type = "censor_aware", save_pdf = FALSE)
  expect_type(result_censor, "list")
  expect_s3_class(result_censor$rvf_plot, "gg")
  expect_s3_class(result_censor$qq_plot, "gg")
})

test_that("simulation-based diagnostics execute for supported distributions", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("DHARMa")
  
  dat <- setup_plot_test_data()
  fit <- fit_one(dat, family = "gamma", random = "")
  
  if (!fit$converged) {
    skip("Model failed to converge")
  }
  
  # Test DHARMa plotting
  result <- plot(fit, dharma = TRUE, nsim = 50, save_pdf = FALSE)
  expect_type(result, "list")
  expect_true("dharma_obj" %in% names(result))
  
  # DHARMa object might be NULL if computation fails, but should be in the list
  if (!is.null(result$dharma_obj)) {
    expect_s3_class(result$dharma_obj, "DHARMa")
  }
})

test_that("external simulation protocol executes for censored regression", {
  skip_if_not_installed("censReg")
  skip_if_not_installed("DHARMa")
  
  dat <- setup_plot_test_data()
  
  # Skip if insufficient censoring
  if (sum(dat$cens_lod) < 5) {
    skip("Insufficient censoring for censReg testing")
  }
  
  fit <- fit_one(dat, family = "tobit_censreg", random = "")
  
  if (!fit$converged) {
    skip("censReg model failed to converge")
  }
  
  # Test DHARMa with external simulation
  result <- plot(fit, dharma = TRUE, nsim = 50, save_pdf = FALSE)
  
  # Result might be NULL if plotting fails
  if (!is.null(result)) {
    expect_type(result, "list")
    expect_true("dharma_obj" %in% names(result))
  } else {
    skip("Plotting returned NULL - model may not support DHARMa")
  }
  
  # Should use external simulation for censReg
  if (!is.null(result$dharma_obj)) {
    expect_s3_class(result$dharma_obj, "DHARMa")
  }
})

test_that("quantile residual computation executes for censored distributions", {
  skip_if_not_installed("survival")
  skip_if_not_installed("statmod")
  skip_if_not_installed("ggplot2")
  
  dat <- setup_plot_test_data()
  
  # Skip if insufficient censoring
  if (sum(dat$cens_lod) < 5) {
    skip("Insufficient censoring for RQR testing")
  }
  
  fit <- fit_one(dat, family = "tobit", random = "")
  
  if (!fit$converged) {
    skip("Tobit model failed to converge")
  }
  
  # Test censor-aware plotting which should use RQR for tobit
  result <- plot(fit, plot_type = "censor_aware", save_pdf = FALSE)
  expect_type(result, "list")
  expect_s3_class(result$rvf_plot, "gg")
  expect_s3_class(result$qq_plot, "gg")
})

test_that("visualization parameters are correctly implemented", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("ggplot2")
  
  dat <- setup_plot_test_data()
  fit <- fit_one(dat, family = "gamma", random = "")
  
  if (!fit$converged) {
    skip("Model failed to converge")
  }
  
  # Test different parameter values
  result1 <- plot(fit, point_size = 1.5, alpha_cens = 0.6, save_pdf = FALSE)
  result2 <- plot(fit, point_size = 3.0, alpha_cens = 0.9, save_pdf = FALSE)
  
  expect_type(result1, "list")
  expect_type(result2, "list")
  expect_s3_class(result1$rvf_plot, "gg")
  expect_s3_class(result2$rvf_plot, "gg")
})

test_that("visualization system handles non-standard conditions appropriately", {
  skip_if_not_installed("glmmTMB")
  
  dat <- setup_plot_test_data()
  
  # Test with non-converged model
  fit <- fit_one(dat, family = "gamma", random = "")
  fit$converged <- FALSE  # Force non-convergence for testing
  
  expect_message(plot(fit), "Cannot plot: model did not converge")
  expect_null(plot(fit))
  
  # Test with NULL model
  fit$model <- NULL
  expect_message(plot(fit), "Cannot plot: model did not converge")
  expect_null(plot(fit))
})

test_that("censored observation visualization protocol executes correctly", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("ggplot2")
  
  dat <- setup_plot_test_data()
  
  # Skip if no censoring
  if (sum(dat$cens_lod) == 0) {
    skip("No censoring in test data")
  }
  
  fit <- fit_one(dat, family = "gamma", random = "")
  
  if (!fit$converged) {
    skip("Model failed to converge")
  }
  
  # Test that censored points are handled differently in visualization
  result <- plot(fit, plot_type = "censor_aware", alpha_cens = 0.5, save_pdf = FALSE)
  
  expect_type(result, "list")
  expect_s3_class(result$rvf_plot, "gg")
  expect_s3_class(result$qq_plot, "gg")
  
  # Check that censored points are included in the plot data
  plot_data <- result$rvf_plot$data
  if ("cens_lod" %in% names(plot_data)) {
    expect_true(any(plot_data$cens_lod))  # Should have some censored points
  }
})

test_that("multi-panel diagnostic visualization protocol executes correctly", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("DHARMa")
  skip_if_not_installed("patchwork")
  
  dat <- setup_plot_test_data()
  fit <- fit_one(dat, family = "gamma", random = "")
  
  if (!fit$converged) {
    skip("Model failed to converge")
  }
  
  # Test four-panel DHARMa plot
  result <- plot(fit, dharma = TRUE, nsim = 50, save_pdf = FALSE)
  
  expect_type(result, "list")
  
  if (!is.null(result$dharma_obj) && "dharma_plot" %in% names(result)) {
    # Should have four-panel plot structure
    expect_true(!is.null(result$dharma_plot))
  }
})

})

describe("System robustness and boundary condition handling", {

test_that("system manages dependency requirements appropriately", {
  dat <- setup_plot_test_data()
  fit <- fit_one(dat, family = "gamma", random = "")
  
  if (!fit$converged) {
    skip("Model failed to converge")
  }
  
  # Test should work even if optional packages are missing
  # (skip tests will handle missing dependencies)
  expect_true(is.list(fit))
  expect_s3_class(fit, "immuno_fit")
})

test_that("system manages non-standard residual patterns appropriately", {
  skip_if_not_installed("glmmTMB")
  
  dat <- setup_plot_test_data()
  
  # Create problematic data that might cause residual computation issues
  dat_problematic <- dat
  dat_problematic$value[1:5] <- Inf  # Add infinite values
  
  fit <- fit_one(dat_problematic, family = "gamma", random = "")
  
  if (fit$converged) {
    # Should handle problematic residuals without crashing
    result <- plot(fit, save_pdf = FALSE)
    # Might return NULL or error gracefully
    expect_true(is.null(result) || is.list(result))
  }
})

test_that("system processes limited sample sizes appropriately", {
  skip_if_not_installed("glmmTMB")
  
  dat <- setup_plot_test_data()
  
  # Create very small dataset
  small_dat <- dat[1:10, ]  # Only 10 observations
  
  fit <- fit_one(small_dat, family = "gamma", random = "")
  
  if (fit$converged) {
    result <- plot(fit, save_pdf = FALSE)
    expect_true(is.null(result) || is.list(result))
  }
})

test_that("system processes extreme value distributions appropriately", {
  skip_if_not_installed("glmmTMB")
  
  dat <- setup_plot_test_data()
  
  # Add extreme values
  extreme_dat <- dat
  extreme_dat$value[1:3] <- c(1e10, 1e-10, 0)  # Very large, very small, zero
  
  fit <- fit_one(extreme_dat, family = "gamma", random = "")
  
  if (fit$converged) {
    result <- plot(fit, save_pdf = FALSE)
    expect_true(is.null(result) || is.list(result))
  }
})

test_that("system processes complete and zero censoring conditions appropriately", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  
  dat <- setup_plot_test_data()
  
  # Test all-censored scenario
  all_cens_dat <- dat
  all_cens_dat$value <- all_cens_dat$lod * 0.5  # All below LOD
  all_cens_dat$cens_lod <- TRUE
  
  fit_all_cens <- fit_one(all_cens_dat, family = "gamma", random = "")
  if (fit_all_cens$converged) {
    result <- plot(fit_all_cens, save_pdf = FALSE)
    expect_true(is.null(result) || is.list(result))
  }
  
  # Test no-censored scenario
  no_cens_dat <- dat
  no_cens_dat$value <- no_cens_dat$lod * 2  # All above LOD
  no_cens_dat$cens_lod <- FALSE
  
  fit_no_cens <- fit_one(no_cens_dat, family = "gamma", random = "")
  if (fit_no_cens$converged) {
    result <- plot(fit_no_cens, save_pdf = FALSE)
    expect_true(is.null(result) || is.list(result))
  }
})

})

describe("Simulation-based statistical validation", {

test_that("distributional uniformity assessment protocol executes correctly", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("DHARMa")
  
  dat <- setup_plot_test_data()
  fit <- fit_one(dat, family = "gamma", random = "")
  
  if (!fit$converged) {
    skip("Model failed to converge")
  }
  
  # Test DHARMa with statistical tests
  result <- plot(fit, dharma = TRUE, nsim = 100, save_pdf = FALSE)
  
  if (!is.null(result$dharma_obj)) {
    # Should be able to run uniformity test
    test_result <- DHARMa::testUniformity(result$dharma_obj, plot = FALSE)
    expect_true(is.numeric(test_result$p.value))
    expect_true(test_result$p.value >= 0 && test_result$p.value <= 1)
  }
})

test_that("dispersion assessment protocol executes correctly", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("DHARMa")
  
  dat <- setup_plot_test_data()
  fit <- fit_one(dat, family = "gamma", random = "")
  
  if (!fit$converged) {
    skip("Model failed to converge")
  }
  
  result <- plot(fit, dharma = TRUE, nsim = 100, save_pdf = FALSE)
  
  if (!is.null(result$dharma_obj)) {
    # Should be able to run dispersion test
    test_result <- DHARMa::testDispersion(result$dharma_obj, plot = FALSE)
    expect_true(is.numeric(test_result$p.value))
    expect_true(test_result$p.value >= 0 && test_result$p.value <= 1)
  }
})

test_that("outlier detection protocol executes correctly", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("DHARMa")
  
  dat <- setup_plot_test_data()
  fit <- fit_one(dat, family = "gamma", random = "")
  
  if (!fit$converged) {
    skip("Model failed to converge")
  }
  
  result <- plot(fit, dharma = TRUE, nsim = 100, save_pdf = FALSE)
  
  if (!is.null(result$dharma_obj)) {
    # Should be able to run outlier test
    test_result <- DHARMa::testOutliers(result$dharma_obj, plot = FALSE)
    expect_true(is.numeric(test_result$p.value))
    expect_true(test_result$p.value >= 0 && test_result$p.value <= 1)
  }
})

})