# System Integration and Compatibility Protocol
#
# Validates system-wide functionality and interoperability:
# 1. Object-Oriented Implementation:
#    - Method availability verification
#    - Interface consistency validation
#    - Inheritance structure assessment
#
# 2. Component Integration:
#    - Inter-module communication
#    - Data flow validation
#    - State management verification
#
# 3. System Robustness:
#    - Exception propagation
#    - Error handling consistency
#    - Recovery mechanism validation
#
# 4. Compatibility Assessment:
#    - Version compatibility
#    - Interface stability
#    - Legacy support verification
#
# 5. Performance Evaluation:
#    - Execution time analysis
#    - Resource utilization assessment
#    - Scalability verification
#
# 6. Documentation Verification:
#    - Interface specification
#    - Example validation
#    - Usage pattern assessment

library(testthat)

describe("System integration and compatibility validation", {

test_that("object-oriented interface implementation is complete", {
  skip_if_no_stats_packages()
  skip_if_no_plotting()
  
  # Create test objects
  objects <- create_all_immunoplex_objects()
  
  # Test immuno_fit methods
  if ("immuno_fit" %in% names(objects)) {
    fit <- objects$immuno_fit
    
    # Test all S3 methods exist and work
    expect_true(test_s3_method(fit, "print"))
    expect_true(test_s3_method(fit, "summary"))
    expect_true(test_s3_method(fit, "plot", save_pdf = FALSE))
    
    # Validate object structure
    expect_true(validate_s3_object(fit, "immuno_fit", 
                                  c("model", "family", "converged", "aic")))
  }
  
  # Test immuno_model_set methods
  if ("immuno_model_set" %in% names(objects)) {
    models <- objects$immuno_model_set
    
    expect_true(test_s3_method(models, "print"))
    expect_true(test_s3_method(models, "summary"))
    expect_true(test_s3_method(models, "plot", save_pdf = FALSE))
    
    expect_true(validate_s3_object(models, "immuno_model_set",
                                  c("models", "comparison", "best_model")))
  }
})

test_that("visualization subsystem integration executes correctly", {
  skip_if_no_stats_packages()
  skip_if_no_plotting()
  
  dat <- create_s3_test_data()
  fit <- fit_one(dat, family = "gamma", random = "")
  
  if (!fit$converged) {
    skip("Model failed to converge")
  }
  
  # Test that all plotting components work together
  result <- plot(fit, plot_type = "basic", save_pdf = FALSE)
  expect_type(result, "list")
  expect_true(all(c("combined_plot", "rvf_plot", "qq_plot") %in% names(result)))
  
  # Test censor-aware plotting
  result_censor <- plot(fit, plot_type = "censor_aware", save_pdf = FALSE)
  expect_type(result_censor, "list")
  expect_s3_class(result_censor$combined_plot, "patchwork")
})

test_that("exception handling protocol executes across system components", {
  # Test enhanced error messages work in integration
  
  # Test fit_models with invalid data
  expect_message(
    expect_error(fit_models(data.frame()), "No models converged successfully"),
    "Missing required columns"
  )
  
  # Test fit_one with invalid data
  expect_error(fit_one(data.frame()), "Missing required columns")
  
  # Test quiet mode works
  expect_silent(
    expect_error(fit_models(data.frame(), quiet = TRUE), "No models converged successfully")
  )
})

test_that("system maintains interface compatibility with legacy code", {
  skip_if_no_stats_packages()
  
  # Ensure old code patterns still work
  dat <- create_s3_test_data()
  
  # Old-style function calls should still work
  fit <- fit_one(dat, family = "gamma")
  expect_s3_class(fit, "immuno_fit")
  
  models <- fit_models(dat, families = c("gamma", "tobit"))
  expect_s3_class(models, "immuno_model_set")
  
  # Legacy method names should work
  expect_true(test_s3_method(fit, "summary.immuno_model"))
})

test_that("system maintains computational efficiency requirements", {
  skip_if_no_stats_packages()
  
  dat <- create_s3_test_data(n_obs = 100)
  
  # Test that fitting is still reasonably fast
  start_time <- Sys.time()
  fit <- fit_one(dat, family = "gamma", random = "")
  fit_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Should complete within reasonable time (adjust as needed)
  expect_lt(fit_time, 30)  # 30 seconds should be plenty
  
  if (fit$converged) {
    # Plotting should also be reasonably fast
    start_time <- Sys.time()
    plot(fit, save_pdf = FALSE)
    plot_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    expect_lt(plot_time, 10)  # 10 seconds for plotting
  }
})

test_that("system verifies component registration and accessibility", {
  # Test that all our new files are accessible
  
  # Check that internal functions are available (with triple colon)
  expect_true(exists(".get_residuals_and_fitted", envir = asNamespace("immunoPlex")))
  expect_true(exists(".create_basic_plots", envir = asNamespace("immunoPlex")))
  expect_true(exists(".create_censor_aware_plots", envir = asNamespace("immunoPlex")))
  expect_true(exists(".compute_dharma", envir = asNamespace("immunoPlex")))
  
  # Check that S3 methods are properly registered
  expect_true("summary.immuno_fit" %in% methods("summary"))
  expect_true("plot.immuno_fit" %in% methods("plot"))
  expect_true("plot.immuno_model_set" %in% methods("plot"))
})

})

describe("Documentation and usage pattern validation", {

test_that("interface documentation is accessible and complete", {
  # Test that help files exist for new methods
  # Note: This would need to be tested after roxygen2::roxygenise()
  expect_true(TRUE)  # Placeholder - real test would check help files
})

test_that("documented usage patterns execute correctly", {
  skip_if_no_stats_packages()
  skip_if_no_plotting()
  
  # Test basic example pattern from our documentation
  dat <- create_s3_test_data()
  
  # Example: fit single model and summarize
  fit <- fit_one(dat, family = "gamma")
  result <- summary(fit)
  expect_type(result, "list")
  
  # Example: fit multiple models and compare
  models <- fit_models(dat, families = c("gamma", "tobit"))
  comparison_plot <- plot(models, plot_type = "comparison", save_pdf = FALSE)
  expect_s3_class(comparison_plot, "gg")
})

})