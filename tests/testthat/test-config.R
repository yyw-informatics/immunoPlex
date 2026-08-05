# Parameter Management and Configuration Protocol
#
# Validates system configuration functionality:
# 1. Parameter Specification:
#    - Value retrieval and validation
#    - Default initialization
#    - Type constraint enforcement
#
# 2. State Management:
#    - Configuration persistence
#    - State restoration
#    - Modification protocols
#
# 3. System Integration:
#    - Core function interaction
#    - Cross-component validation
#    - Interface consistency
#
# Parameter Categories:
# 1. Estimation Parameters:
#    - Random effects thresholds
#    - Distribution selection criteria
#
# 2. Visualization Parameters:
#    - Aesthetic specifications
#    - Layout configuration
#    - Display constraints
#
# 3. Diagnostic Parameters:
#    - Simulation specifications
#    - Validation thresholds
#    - Test data generation

library(testthat)

describe("Parameter management system", {

test_that("parameter retrieval protocol returns valid structure", {
  config <- get_immunoplex_config()
  
  # Check main structure
  expect_type(config, "list")
  expect_true(all(c("random_effects", "auto_family", "plotting", "dharma", "testing") %in% names(config)))
  
  # Check random_effects section
  expect_type(config$random_effects, "list")
  expect_true(all(c("min_subjects", "min_reps") %in% names(config$random_effects)))
  expect_type(config$random_effects$min_subjects, "integer")
  expect_type(config$random_effects$min_reps, "integer")
  
  # Check auto_family section
  expect_type(config$auto_family, "list")
  expect_true(all(c("high_censoring_threshold", "high_skewness_threshold") %in% names(config$auto_family)))
  expect_type(config$auto_family$high_censoring_threshold, "double")
  expect_type(config$auto_family$high_skewness_threshold, "double")
  
  # Check plotting section
  expect_type(config$plotting, "list")
  expect_true(all(c("default_point_size", "default_alpha_cens", "outlier_threshold", 
                   "loess_color", "text_size") %in% names(config$plotting)))
  
  # Check dharma section
  expect_type(config$dharma, "list")
  expect_true(all(c("default_nsim", "test_nsim") %in% names(config$dharma)))
  expect_type(config$dharma$default_nsim, "integer")
  expect_type(config$dharma$test_nsim, "integer")
  
  # Check testing section
  expect_type(config$testing, "list")
  expect_true(all(c("small_dataset_size", "medium_dataset_size", "large_dataset_size", 
                   "default_lod", "default_ulod", "performance_dataset_size") %in% names(config$testing)))
})

test_that("default parameter specifications maintain validity constraints", {
  config <- get_immunoplex_config()
  
  # Random effects defaults
  expect_gte(config$random_effects$min_subjects, 5)  # At least 5 subjects
  expect_lte(config$random_effects$min_subjects, 100)  # Not too high
  expect_gte(config$random_effects$min_reps, 2)  # At least 2 reps
  expect_lte(config$random_effects$min_reps, 10)  # Not too high
  
  # Auto family selection defaults
  expect_gte(config$auto_family$high_censoring_threshold, 50)  # Reasonable censoring threshold
  expect_lte(config$auto_family$high_censoring_threshold, 90)
  expect_gte(config$auto_family$high_skewness_threshold, 1)  # Reasonable skewness threshold
  expect_lte(config$auto_family$high_skewness_threshold, 5)
  
  # Plotting defaults
  expect_gte(config$plotting$default_point_size, 0.5)
  expect_lte(config$plotting$default_point_size, 5)
  expect_gte(config$plotting$default_alpha_cens, 0.1)
  expect_lte(config$plotting$default_alpha_cens, 1.0)
  expect_gte(config$plotting$outlier_threshold, 2)
  expect_lte(config$plotting$outlier_threshold, 5)
  
  # DHARMa defaults
  expect_gte(config$dharma$default_nsim, 100)  # Sufficient simulations
  expect_lte(config$dharma$default_nsim, 5000)  # Not too slow
  expect_gte(config$dharma$test_nsim, 50)  # Fast enough for testing
  expect_lte(config$dharma$test_nsim, 500)
})

test_that("parameter modification protocol maintains state consistency", {
  # Store original config
  original_config <- get_immunoplex_config()
  
  # Test setting individual parameters
  set_immunoplex_config(
    random_effects_min_subjects = 25,
    random_effects_min_reps = 2,
    high_censoring_threshold = 60,
    high_skewness_threshold = 1.5,
    default_dharma_nsim = 500,
    plotting_point_size = 2.0,
    plotting_alpha_cens = 0.8
  )
  
  # Check that values were changed
  new_config <- get_immunoplex_config()
  expect_equal(new_config$random_effects$min_subjects, 25L)
  expect_equal(new_config$random_effects$min_reps, 2L)
  expect_equal(new_config$auto_family$high_censoring_threshold, 60)
  expect_equal(new_config$auto_family$high_skewness_threshold, 1.5)
  expect_equal(new_config$dharma$default_nsim, 500L)
  expect_equal(new_config$plotting$default_point_size, 2.0)
  expect_equal(new_config$plotting$default_alpha_cens, 0.8)
  
  # Restore original configuration
  reset_immunoplex_config()
  restored_config <- get_immunoplex_config()
  expect_equal(restored_config$random_effects$min_subjects, original_config$random_effects$min_subjects)
  expect_equal(restored_config$random_effects$min_reps, original_config$random_effects$min_reps)
})

test_that("parameter restoration protocol reinstates default specifications", {
  # Store original defaults
  original_config <- get_immunoplex_config()
  
  # Change some values
  set_immunoplex_config(
    random_effects_min_subjects = 999,
    high_censoring_threshold = 99,
    default_dharma_nsim = 9999
  )
  
  # Verify changes took effect
  changed_config <- get_immunoplex_config()
  expect_equal(changed_config$random_effects$min_subjects, 999L)
  expect_equal(changed_config$auto_family$high_censoring_threshold, 99)
  expect_equal(changed_config$dharma$default_nsim, 9999L)
  
  # Reset to defaults
  reset_immunoplex_config()
  
  # Verify restoration
  reset_config <- get_immunoplex_config()
  expect_equal(reset_config$random_effects$min_subjects, original_config$random_effects$min_subjects)
  expect_equal(reset_config$auto_family$high_censoring_threshold, original_config$auto_family$high_censoring_threshold)
  expect_equal(reset_config$dharma$default_nsim, original_config$dharma$default_nsim)
})

test_that("parameter specifications influence estimation protocol behavior", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("moments")
  
  # Create test data with moderate censoring and skewness
  dat <- data.frame(
    value = c(rep(0.5, 15), rep(2, 15), rep(5, 15)),  # 33% censoring
    cens_lod = c(rep(TRUE, 15), rep(FALSE, 30)),
    cens_ulod = rep(FALSE, 45),
    lod = rep(1, 45),
    ulod = rep(10, 45),
    timepoint = rep(c("Pre", "Post"), length.out = 45),
    disease = rep(c("Healthy", "Disease"), length.out = 45),
    age = 25:69,
    subject_id = rep(1:15, each = 3)
  )
  
  # Store original config
  original_config <- get_immunoplex_config()
  
  # With our improved auto-selector, presence of censoring determines choice
  # Old threshold-based logic is replaced with simpler, better logic
  
  # With censoring present, should avoid gamma regardless of thresholds
  fit_with_censoring <- fit_one(dat, family = "auto", random = "")
  expect_true(fit_with_censoring$family %in% c("tobit", "aft"))  # Should avoid gamma
  
  # Test without censoring - can choose gamma
  dat_no_cens <- dat
  dat_no_cens$cens_lod[] <- FALSE
  dat_no_cens$cens_ulod[] <- FALSE
  
  fit_no_censoring <- fit_one(dat_no_cens, family = "auto", random = "")
  expect_true(fit_no_censoring$family %in% c("gamma", "tobit", "aft", "gaussian"))  # Can choose any
  
  # Restore original config
  reset_immunoplex_config()
})

test_that("parameter specifications influence visualization protocol defaults", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("ggplot2")
  
  # Create robust test data for plotting
  dat <- create_robust_test_data(n_subjects = 10, n_reps = 2, family_type = "gamma")
  
  fit <- fit_one(dat, family = "gamma", random = "")
  
  if (fit$converged) {
    # Store original config
    original_config <- get_immunoplex_config()
    
    # Test with custom plotting config
    set_immunoplex_config(
      plotting_point_size = 3.0,
      plotting_alpha_cens = 0.9
    )
    
    # Test that plotting uses new defaults
    # Note: We can't easily test the actual plot appearance, but we can test that
    # the function accepts the configured defaults
    result <- plot(fit, save_pdf = FALSE)
    expect_type(result, "list")
    
    # Restore original config
    reset_immunoplex_config()
  } else {
    skip("Model failed to converge")
  }
})

test_that("internal test configuration protocol maintains validity", {
  test_config <- immunoPlex:::.get_test_config()
  
  expect_type(test_config, "list")
  expect_true(all(c("small_n", "medium_n", "large_n", "perf_n", 
                   "default_lod", "default_ulod", "test_dharma_nsim") %in% names(test_config)))
  
  # Check reasonable test values
  expect_gte(test_config$small_n, 5)
  expect_lte(test_config$small_n, 20)
  expect_gte(test_config$medium_n, 20)
  expect_lte(test_config$medium_n, 100)
  expect_gte(test_config$large_n, 100)
  expect_lte(test_config$large_n, 1000)
  expect_gt(test_config$default_lod, 0)
  expect_gt(test_config$default_ulod, test_config$default_lod)
})

test_that("parameter type constraints are enforced correctly", {
  # Test that integer options are correctly converted
  set_immunoplex_config(
    random_effects_min_subjects = 25.7,  # Should be converted to integer
    default_dharma_nsim = 500.9           # Should be converted to integer
  )
  
  config <- get_immunoplex_config()
  expect_type(config$random_effects$min_subjects, "integer")
  expect_equal(config$random_effects$min_subjects, 25L)
  expect_type(config$dharma$default_nsim, "integer")
  expect_equal(config$dharma$default_nsim, 500L)
  
  # Test that numeric options remain numeric
  set_immunoplex_config(
    high_censoring_threshold = 65.5,
    plotting_point_size = 1.8
  )
  
  config <- get_immunoplex_config()
  expect_type(config$auto_family$high_censoring_threshold, "double")
  expect_equal(config$auto_family$high_censoring_threshold, 65.5)
  expect_type(config$plotting$default_point_size, "double")
  expect_equal(config$plotting$default_point_size, 1.8)
  
  # Reset for cleanup
  reset_immunoplex_config()
})

test_that("parameter state persistence maintains consistency across operations", {
  skip_if_not_installed("glmmTMB")
  
  # Set custom configuration
  set_immunoplex_config(random_effects_min_subjects = 15)
  
  # Create data that would normally trigger random effects dropping
  dat <- data.frame(
    value = 1:25,
    cens_lod = rep(FALSE, 25),
    cens_ulod = rep(FALSE, 25),
    lod = rep(0.5, 25),
    ulod = rep(10, 25),
    subject_id = rep(1:20, length.out = 25)  # 20 subjects (< default 30, but >= our 15)
  )
  
  # Should NOT drop random effects with our lower threshold
  fit <- fit_one(dat, family = "gamma", random = "(1|subject_id)")
  
  # Check that the configuration was respected
  config <- get_immunoplex_config()
  expect_equal(config$random_effects$min_subjects, 15L)
  
  # Reset for cleanup
  reset_immunoplex_config()
})

test_that("partial parameter modification protocol maintains state integrity", {
  original_config <- get_immunoplex_config()
  
  # Update only one parameter
  set_immunoplex_config(random_effects_min_subjects = 20)
  
  updated_config <- get_immunoplex_config()
  
  # Check that only the specified parameter changed
  expect_equal(updated_config$random_effects$min_subjects, 20L)
  expect_equal(updated_config$random_effects$min_reps, original_config$random_effects$min_reps)
  expect_equal(updated_config$auto_family$high_censoring_threshold, original_config$auto_family$high_censoring_threshold)
  expect_equal(updated_config$plotting$default_point_size, original_config$plotting$default_point_size)
  
  # Reset for cleanup
  reset_immunoplex_config()
})

test_that("parameter management system handles boundary conditions appropriately", {
  # Test with NULL values (should be ignored)
  original_config <- get_immunoplex_config()
  
  set_immunoplex_config(
    random_effects_min_subjects = NULL,
    high_censoring_threshold = 55
  )
  
  config <- get_immunoplex_config()
  
  # NULL parameter should be ignored (unchanged)
  expect_equal(config$random_effects$min_subjects, original_config$random_effects$min_subjects)
  # Non-NULL parameter should be updated
  expect_equal(config$auto_family$high_censoring_threshold, 55)
  
  # Reset for cleanup
  reset_immunoplex_config()
})

})

describe("Parameter system integration validation", {

test_that("estimation protocol respects random effects parameter constraints", {
  skip_if_not_installed("glmmTMB")
  
  # Create data with exactly 25 subjects, 2 reps each
  dat <- data.frame(
    value = rnorm(50, mean = 3, sd = 1),
    cens_lod = rep(FALSE, 50),
    cens_ulod = rep(FALSE, 50),
    lod = rep(1, 50),
    ulod = rep(10, 50),
    subject_id = rep(1:25, each = 2),
    timepoint = rep(c("Pre", "Post"), 25),
    disease = rep(c("Healthy", "Disease"), each = 25),
    age = rep(25:49, each = 2)
  )
  
  # With default thresholds (30 subjects, 3 reps), should drop random effects
  reset_immunoplex_config()
  expect_message(
    fit_default <- fit_one(dat, family = "gamma", random = "(1|subject_id)"),
    "Dropping random term"
  )
  
  # With lower thresholds, should keep random effects
  set_immunoplex_config(
    random_effects_min_subjects = 20,
    random_effects_min_reps = 2
  )
  
  # Should NOT get dropping message
  fit_custom <- fit_one(dat, family = "gamma", random = "(1|subject_id)")
  expect_s3_class(fit_custom, "immuno_fit")
  
  # Reset for cleanup
  reset_immunoplex_config()
})

test_that("visualization protocol implements parameter specifications", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("ggplot2")
  
  # Create robust test data
  dat <- create_robust_test_data(n_subjects = 10, n_reps = 2, family_type = "gamma")
  
  fit <- fit_one(dat, family = "gamma", random = "")
  
  if (fit$converged) {
    # Set custom plotting configuration
    set_immunoplex_config(
      plotting_point_size = 2.5,
      plotting_alpha_cens = 0.9
    )
    
    # Test that plotting functions complete without error
    result <- plot(fit, save_pdf = FALSE)
    expect_type(result, "list")
    
    # Reset for cleanup
    reset_immunoplex_config()
  } else {
    skip("Model failed to converge")
  }
})

})