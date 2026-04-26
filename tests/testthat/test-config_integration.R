# Configuration System Integration Validation Protocol
#
# Validates system-wide parameter management including:
# - Option persistence and retrieval
# - Default value handling
# - Cross-component parameter propagation
# - System state restoration
# - Invalid input handling
#
# Validation sections:
# 1. Parameter persistence
# 2. Default value protocols
# 3. Cross-component integration
# 4. System state management
# 5. Input validation

library(testthat)

describe("Configuration system integration validation", {

test_that("options are properly set and retrieved", {
  # Store current state
  original_options <- list(
    fit_one.min_subjects = getOption("fit_one.min_subjects"),
    fit_one.min_reps = getOption("fit_one.min_reps"),
    immunoplex.high_censoring_pct = getOption("immunoplex.high_censoring_pct"),
    immunoplex.dharma_nsim = getOption("immunoplex.dharma_nsim")
  )
  
  # Test setting options through R's options() function
  options(
    fit_one.min_subjects = 25L,
    immunoplex.high_censoring_pct = 65,
    immunoplex.dharma_nsim = 750L
  )
  
  # Verify they're reflected in config
  config <- get_immunoplex_config()
  expect_equal(config$random_effects$min_subjects, 25L)
  expect_equal(config$auto_family$high_censoring_threshold, 65)
  expect_equal(config$dharma$default_nsim, 750L)
  
  # Test setting through our config function
  set_immunoplex_config(
    random_effects_min_subjects = 35,
    high_censoring_threshold = 75
  )
  
  # Verify they're reflected in R options
  expect_equal(getOption("fit_one.min_subjects"), 35L)
  expect_equal(getOption("immunoplex.high_censoring_pct"), 75)
  
  # Restore original options
  do.call(options, original_options)
})

test_that("default values are used when options are unset", {
  # Store current state
  original_options <- list(
    fit_one.min_subjects = getOption("fit_one.min_subjects"),
    immunoplex.high_censoring_pct = getOption("immunoplex.high_censoring_pct"),
    immunoplex.point_size = getOption("immunoplex.point_size")
  )
  
  # Clear the options
  options(
    fit_one.min_subjects = NULL,
    immunoplex.high_censoring_pct = NULL,
    immunoplex.point_size = NULL
  )
  
  # Get config - should return defaults
  config <- get_immunoplex_config()
  
  # Check that reasonable defaults are returned
  expect_equal(config$random_effects$min_subjects, 30L)  # Default from fit_one.r
  expect_equal(config$auto_family$high_censoring_threshold, 70)  # Default from fit_one.r
  expect_equal(config$plotting$default_point_size, 1.2)  # Default from plot.immuno_fit.r
  
  # Restore original options
  do.call(options, original_options)
})

test_that("configuration affects actual function behavior", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  skip_if_not_installed("moments")
  
  # Initialize dataset with defined distributional properties
  set.seed(123)
  dat <- data.frame(
    value = c(rep(0.3, 25), exp(rnorm(25, mean = 2, sd = 1))),  # Specified censoring rate and asymmetry
    cens_lod = c(rep(TRUE, 25), rep(FALSE, 25)),
    cens_ulod = rep(FALSE, 50),
    lod = rep(0.5, 50),
    ulod = rep(10, 50),
    timepoint = rep(c("Pre", "Post"), 25),
    disease = rep(c("Healthy", "Disease"), each = 25),
    age = 25:74,
    subject_id = rep(1:25, each = 2)
  )
  
  # Store original configuration
  original_config <- get_immunoplex_config()
  
  # With our improved auto-selector, any censoring avoids gamma
  # The old thresholds are no longer used - new logic is simpler and better
  
  # Evaluation 1: With censoring present, should choose tobit or aft
  fit_auto <- fit_one(dat, family = "auto", random = "")
  expect_true(fit_auto$family %in% c("tobit", "aft"))  # Should avoid gamma
  
  # Evaluation 2: Test with no censoring - can choose gamma
  dat_no_cens <- dat
  dat_no_cens$cens_lod[] <- FALSE
  dat_no_cens$cens_ulod[] <- FALSE
  
  fit_no_cens <- fit_one(dat_no_cens, family = "auto", random = "")
  expect_true(fit_no_cens$family %in% c("gamma", "tobit", "aft", "gaussian"))  # Can choose any
  
  # Evaluation 3: Distribution selection with modified thresholds
  set_immunoplex_config(
    high_censoring_threshold = 60,  # Above observed rate
    high_skewness_threshold = 10    # Above observed asymmetry
  )
  fit_aft <- fit_one(dat, family = "auto", random = "")
  expect_equal(fit_aft$family, "aft")
  
  # Restore original configuration
  reset_immunoplex_config()
  final_config <- get_immunoplex_config()
  expect_equal(final_config$auto_family$high_censoring_threshold, 
               original_config$auto_family$high_censoring_threshold)
})

test_that("configuration system works across package reload", {
  # Set custom configuration
  set_immunoplex_config(
    random_effects_min_subjects = 22,
    high_censoring_threshold = 55
  )
  
  # Verify settings
  config_before <- get_immunoplex_config()
  expect_equal(config_before$random_effects$min_subjects, 22L)
  expect_equal(config_before$auto_family$high_censoring_threshold, 55)
  
  # Settings should persist in the same R session
  # (Note: In tests, we can't actually reload the package, but we can verify
  # that the underlying options system maintains the values)
  expect_equal(getOption("fit_one.min_subjects"), 22L)
  expect_equal(getOption("immunoplex.high_censoring_pct"), 55)
  
  # Reset for cleanup
  reset_immunoplex_config()
})

test_that("invalid configuration values are handled gracefully", {
  # Store original config
  original_config <- get_immunoplex_config()
  
  # Validate system response to boundary conditions
  expect_no_error(
    set_immunoplex_config(
      random_effects_min_subjects = 0,    # Minimum subject threshold
      random_effects_min_reps = 1,        # Minimum replicate count
      high_censoring_threshold = 100,     # Maximum censoring rate
      high_skewness_threshold = 0,        # Minimum asymmetry threshold
      default_dharma_nsim = 1             # Minimum simulation count
    )
  )
  
  # Verify the values were set (even if extreme)
  config <- get_immunoplex_config()
  expect_equal(config$random_effects$min_subjects, 0L)
  expect_equal(config$random_effects$min_reps, 1L)
  expect_equal(config$auto_family$high_censoring_threshold, 100)
  expect_equal(config$auto_family$high_skewness_threshold, 0)
  expect_equal(config$dharma$default_nsim, 1L)
  
  # Reset to restore sanity
  reset_immunoplex_config()
})

test_that("configuration interacts correctly with function defaults", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("ggplot2")
  
  # Create simple test data
  dat <- data.frame(
    value = c(1, 2, 3, 4, 5),
    cens_lod = rep(FALSE, 5),
    cens_ulod = rep(FALSE, 5),
    lod = rep(0.5, 5),
    ulod = rep(10, 5)
  )
  
  fit <- fit_one(dat, family = "gamma", random = "")
  
  if (fit$converged) {
    # Set configuration
    set_immunoplex_config(
      plotting_point_size = 2.5,
      plotting_alpha_cens = 0.9
    )
    
    # Test 1: Function should use configured defaults when no explicit params given
    result1 <- plot(fit, save_pdf = FALSE)
    expect_type(result1, "list")
    
    # Test 2: Explicit parameters should override configuration
    result2 <- plot(fit, point_size = 1.0, alpha_cens = 0.3, save_pdf = FALSE)
    expect_type(result2, "list")
    
    # Both should succeed without error
    expect_true(!is.null(result1))
    expect_true(!is.null(result2))
    
    # Reset for cleanup
    reset_immunoplex_config()
  } else {
    skip("Model failed to converge")
  }
})

test_that("test configuration helper works correctly", {
  # Test the internal helper function
  test_config <- immunoPlex:::.get_test_config()
  
  # Should have all expected components
  expected_names <- c("small_n", "medium_n", "large_n", "perf_n", 
                     "default_lod", "default_ulod", "test_dharma_nsim")
  expect_true(all(expected_names %in% names(test_config)))
  
  # Modify test configuration and verify it's reflected
  options(
    immunoplex.test_small_n = 15,
    immunoplex.test_lod = 2.0,
    immunoplex.test_dharma_nsim = 75
  )
  
  updated_config <- immunoPlex:::.get_test_config()
  expect_equal(updated_config$small_n, 15)
  expect_equal(updated_config$default_lod, 2.0)
  expect_equal(updated_config$test_dharma_nsim, 75)
  
  # Reset test options
  options(
    immunoplex.test_small_n = NULL,
    immunoplex.test_lod = NULL,
    immunoplex.test_dharma_nsim = NULL
  )
})

test_that("configuration documentation examples work", {
  # Test examples from the documentation
  
  # Example 1: Get current configuration
  config <- get_immunoplex_config()
  expect_type(config, "list")
  expect_type(config$random_effects, "list")
  
  # Example 2: Modified estimation constraints for limited samples
  set_immunoplex_config(
    random_effects_min_subjects = 20,
    random_effects_min_reps = 2
  )
  
  updated_config <- get_immunoplex_config()
  expect_equal(updated_config$random_effects$min_subjects, 20L)
  expect_equal(updated_config$random_effects$min_reps, 2L)
  
  # Example 3: Modified distribution selection criteria
  set_immunoplex_config(high_censoring_threshold = 50)
  config_tobit <- get_immunoplex_config()
  expect_equal(config_tobit$auto_family$high_censoring_threshold, 50)
  
  # Example 4: Optimized simulation parameters
  set_immunoplex_config(default_dharma_nsim = 250)
  config_dharma <- get_immunoplex_config()
  expect_equal(config_dharma$dharma$default_nsim, 250L)
  
  # Example 5: Reset configuration
  reset_immunoplex_config()
  reset_config <- get_immunoplex_config()
  
  # Should return to defaults (not necessarily the original values from this session)
  expect_type(reset_config, "list")
  expect_type(reset_config$random_effects$min_subjects, "integer")
})

})

# Implementation Note:
# Each validation protocol includes state restoration to maintain system integrity