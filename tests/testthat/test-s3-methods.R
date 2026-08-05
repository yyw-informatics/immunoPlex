# Object-Oriented Interface Validation Protocol
#
# Validates method implementation and inheritance structure:
# 1. Object Representation:
#    - String representation
#    - Attribute accessibility
#    - State validation
#
# 2. Object Analysis:
#    - Statistical summaries
#    - Diagnostic information
#    - Metadata extraction
#
# 3. Object Visualization:
#    - Diagnostic plots
#    - Comparative analysis
#    - Interactive display
#
# 4. Interface Stability:
#    - Version compatibility
#    - Method deprecation
#    - Error handling

library(testthat)

# Initialize validation dataset with defined statistical properties
setup_s3_test_data <- function() {
  setup_robust_test_data(prefer_real = TRUE, family_hint = "gamma")  # Primary: empirical data
}

describe("Single model object interface implementation", {

test_that("string representation protocol executes correctly", {
  skip_if_not_installed("glmmTMB")
  
  dat <- setup_s3_test_data()
  fit <- fit_one(dat, family = "gamma", random = "")
  
  # Test print output (including new estimand field)
  expect_output(print(fit), "immunoPlex Model")
  expect_output(print(fit), "Family:")
  expect_output(print(fit), "Estimand:")
  expect_output(print(fit), "Converged:")
  expect_output(print(fit), "AIC:")
})

test_that("statistical summary protocol generates comprehensive output", {
  skip_if_not_installed("glmmTMB")
  
  dat <- setup_s3_test_data()
  fit <- fit_one(dat, family = "gamma", random = "")
  
  if (!fit$converged) {
    skip("Model failed to converge")
  }
  
  # Test summary output (including new estimand field)
  expect_output(summary(fit), "immunoPlex Model Summary")
  expect_output(summary(fit), "Family:")
  expect_output(summary(fit), "Estimand:")
  expect_output(summary(fit), "Converged:")
  expect_output(summary(fit), "Model Fit:")
  expect_output(summary(fit), "AIC:")
  expect_output(summary(fit), "Fixed Effects:")
  
  # Test return value (including new estimand field)
  result <- summary(fit)
  expect_type(result, "list")
  expect_true(all(c("family", "estimand", "converged", "n_obs", "aic") %in% names(result)))
  expect_equal(result$estimand, "ratio_of_means")  # Gamma should be ratio_of_means
})

test_that("visualization protocol supports multiple diagnostic modes", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  
  dat <- setup_s3_test_data()
  fit <- fit_one(dat, family = "gamma", random = "")
  
  if (!fit$converged) {
    skip("Model failed to converge")
  }
  
  # Test basic plotting
  result_basic <- plot(fit, plot_type = "basic", save_pdf = FALSE)
  expect_type(result_basic, "list")
  expect_s3_class(result_basic$rvf_plot, "gg")
  expect_s3_class(result_basic$qq_plot, "gg")
  expect_s3_class(result_basic$combined_plot, "patchwork")
  
  # Test censor-aware plotting
  result_censor <- plot(fit, plot_type = "censor_aware", save_pdf = FALSE)
  expect_type(result_censor, "list")
  expect_s3_class(result_censor$rvf_plot, "gg")
  expect_s3_class(result_censor$qq_plot, "gg")
})

test_that("visualization system manages non-convergent models appropriately", {
  skip_if_not_installed("glmmTMB")
  
  dat <- setup_s3_test_data()
  fit <- fit_one(dat, family = "gamma", random = "")
  fit$converged <- FALSE  # Force non-convergence for testing
  
  expect_message(plot(fit), "Cannot plot: model did not converge")
  expect_null(plot(fit))
})

})

describe("Model set object interface implementation", {

test_that("model set string representation protocol executes correctly", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  
  dat <- setup_s3_test_data()
  models <- fit_models(dat, families = c("gamma", "tobit"), random = "")
  
  expect_output(print(models), "immunoPlex Model Set")
  expect_output(print(models), "Models fitted:")
  expect_output(print(models), "Best model:")
  expect_output(print(models), "Model Comparison:")
  
  # Test that estimand information appears in best model display
  output <- capture.output(print(models))
  best_model_line <- output[grepl("Best model:", output)]
  expect_true(any(grepl("ratio_of_", best_model_line)))
})

test_that("model set summary protocol generates comprehensive output", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  
  dat <- setup_s3_test_data()
  models <- fit_models(dat, families = c("gamma", "tobit"), random = "")
  
  expect_output(summary(models), "immuno_model_set Summary")
  expect_output(summary(models), "Number of models")
  expect_output(summary(models), "Best Model Details")
})

test_that("model set visualization protocol supports multiple diagnostic modes", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  skip_if_not_installed("ggplot2")
  
  dat <- setup_s3_test_data()
  models <- fit_models(dat, families = c("gamma", "tobit"), random = "")
  
  # Test comparison plot
  comp_plot <- plot(models, plot_type = "comparison", save_pdf = FALSE)
  expect_s3_class(comp_plot, "gg")
  
  # Test best model residuals
  expect_message(plot(models, plot_type = "best_residuals", save_pdf = FALSE),
                "Plotting residuals for best model")
  
  # Test all residuals
  expect_message(plot(models, plot_type = "all_residuals", save_pdf = FALSE),
                "Plotting residuals for")
})

})

describe("Interface version compatibility validation", {

test_that("legacy interface implementation maintains compatibility", {
  skip_if_not_installed("glmmTMB")
  
  dat <- setup_s3_test_data()
  fit <- fit_one(dat, family = "gamma", random = "")
  
  if (!fit$converged) {
    skip("Model failed to converge")
  }
  
  # Should work with legacy method name
  expect_output(summary.immuno_model(fit), "immunoPlex Model Summary")
})

test_that("legacy interface validates object class requirements", {
  expect_error(summary.immuno_model(list()), 
               "Object must be of class 'immuno_fit' or 'immuno_model'")
})

})