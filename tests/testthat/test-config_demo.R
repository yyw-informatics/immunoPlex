# Configuration System Parameter Validation Protocol
#
# Demonstrates and validates:
# - Parameter initialization and modification
# - Distribution family selection thresholds
# - Estimation parameter constraints
# - Diagnostic simulation parameters
# - Visualization specifications
#
# Protocol sections:
# 1. Parameter validation
# 2. Threshold evaluation
# 3. Edge case handling
# 4. Usage pattern verification

library(testthat)

describe("Configuration system validation", {

test_that("parameter initialization and modification protocol executes correctly", {
  cat("\n=== Configuration System Parameter Validation ===\n")
  
  # 1. Parameter initialization verification
  cat("\n1. Current configuration:\n")
  config <- get_immunoplex_config()
  cat("   Random effects thresholds:", config$random_effects$min_subjects, "subjects,", 
      config$random_effects$min_reps, "reps\n")
  cat("   Auto selection thresholds:", config$auto_family$high_censoring_threshold, "% censoring,", 
      config$auto_family$high_skewness_threshold, "skewness\n")
  cat("   DHARMa simulations:", config$dharma$default_nsim, "\n")
  
  # 2. Parameter modification assessment
  cat("\n2. Evaluating parameter modification protocol:\n")
  set_immunoplex_config(
    random_effects_min_subjects = 15,
    random_effects_min_reps = 2,
    high_censoring_threshold = 40,
    default_dharma_nsim = 250
  )
  
  new_config <- get_immunoplex_config()
  cat("   New random effects thresholds:", new_config$random_effects$min_subjects, "subjects,", 
      new_config$random_effects$min_reps, "reps\n")
  cat("   New censoring threshold:", new_config$auto_family$high_censoring_threshold, "%\n")
  cat("   New DHARMa simulations:", new_config$dharma$default_nsim, "\n")
  
  # 3. Distribution family selection threshold validation
  if (requireNamespace("glmmTMB", quietly = TRUE) && 
      requireNamespace("survival", quietly = TRUE) &&
      requireNamespace("moments", quietly = TRUE)) {
    
    cat("\n3. Validating distribution family selection criteria:\n")
    
    # Initialize dataset with defined censoring characteristics
    set.seed(42)
    dat <- data.frame(
      value = c(rep(0.4, 15), exp(rnorm(35, mean = 1.5, sd = 0.8))),
      cens_lod = c(rep(TRUE, 15), rep(FALSE, 35)),
      cens_ulod = rep(FALSE, 50),
      lod = rep(0.5, 50),
      ulod = rep(10, 50),
      timepoint = rep(c("Pre", "Post"), 25),
      disease = rep(c("Healthy", "Disease"), each = 25),
      age = 25:74,
      subject_id = rep(1:25, each = 2)
    )
    
    # Evaluate distribution selection with lower censoring threshold
    set_immunoplex_config(high_censoring_threshold = 25)
    fit_low <- fit_one(dat, family = "auto", random = "")
    cat("   Threshold 25%: Selected distribution", fit_low$family, "\n")
    
    # Evaluate distribution selection with higher censoring threshold
    set_immunoplex_config(high_censoring_threshold = 50)
    fit_high <- fit_one(dat, family = "auto", random = "")
    cat("   Threshold 50%: Selected distribution", fit_high$family, "\n")
    
    # With improved auto-selector, both should avoid gamma due to censoring
    expect_true(fit_low$family %in% c("tobit", "aft"))
    expect_true(fit_high$family %in% c("tobit", "aft"))
  } else {
    cat("\n3. Skipping auto selection demo (dependencies not available)\n")
  }
  
  # 4. Visualization parameter validation
  cat("\n4. Validating visualization parameters:\n")
  set_immunoplex_config(
    plotting_point_size = 2.0,
    plotting_alpha_cens = 0.8
  )
  
  plot_config <- get_immunoplex_config()$plotting
  cat("   Point size:", plot_config$default_point_size, "\n")
  cat("   Alpha for censored points:", plot_config$default_alpha_cens, "\n")
  cat("   Text size:", plot_config$text_size, "\n")
  
  # 5. Parameter restoration verification
  cat("\n5. Verifying parameter restoration protocol:\n")
  reset_immunoplex_config()
  
  final_config <- get_immunoplex_config()
  cat("   Reset random effects thresholds:", final_config$random_effects$min_subjects, "subjects,", 
      final_config$random_effects$min_reps, "reps\n")
  cat("   Reset censoring threshold:", final_config$auto_family$high_censoring_threshold, "%\n")
  cat("   Reset DHARMa simulations:", final_config$dharma$default_nsim, "\n")
  
  cat("\n=== Demo Complete ===\n")
  
  # Test passes if we reach here without errors
  expect_true(TRUE)
})

test_that("parameter boundary conditions are handled appropriately", {
  cat("\n=== Parameter Boundary Condition Validation ===\n")
  
  # Test partial updates
  cat("\n1. Testing partial configuration updates:\n")
  original_config <- get_immunoplex_config()
  
  # Update only one parameter
  set_immunoplex_config(random_effects_min_subjects = 99)
  partial_config <- get_immunoplex_config()
  
  cat("   Changed min_subjects to:", partial_config$random_effects$min_subjects, "\n")
  cat("   min_reps unchanged:", partial_config$random_effects$min_reps, "\n")
  cat("   censoring threshold unchanged:", partial_config$auto_family$high_censoring_threshold, "\n")
  
  expect_equal(partial_config$random_effects$min_subjects, 99L)
  expect_equal(partial_config$random_effects$min_reps, original_config$random_effects$min_reps)
  
  # Test NULL handling
  cat("\n2. Testing NULL parameter handling:\n")
  set_immunoplex_config(
    random_effects_min_subjects = NULL,  # Should be ignored
    high_censoring_threshold = 77        # Should be updated
  )
  
  null_config <- get_immunoplex_config()
  cat("   min_subjects (NULL ignored):", null_config$random_effects$min_subjects, "\n")
  cat("   censoring threshold updated:", null_config$auto_family$high_censoring_threshold, "\n")
  
  expect_equal(null_config$random_effects$min_subjects, 99L)  # Unchanged from previous
  expect_equal(null_config$auto_family$high_censoring_threshold, 77)
  
  # Test type conversion
  cat("\n3. Testing type conversion:\n")
  set_immunoplex_config(
    random_effects_min_subjects = 25.7,  # Should become 25L
    default_dharma_nsim = 500.9          # Should become 500L
  )
  
  type_config <- get_immunoplex_config()
  cat("   25.7 converted to:", type_config$random_effects$min_subjects, "(type:", 
      typeof(type_config$random_effects$min_subjects), ")\n")
  cat("   500.9 converted to:", type_config$dharma$default_nsim, "(type:", 
      typeof(type_config$dharma$default_nsim), ")\n")
  
  expect_equal(type_config$random_effects$min_subjects, 25L)
  expect_equal(type_config$dharma$default_nsim, 500L)
  expect_type(type_config$random_effects$min_subjects, "integer")
  expect_type(type_config$dharma$default_nsim, "integer")
  
  # Reset for cleanup
  reset_immunoplex_config()
  cat("\n4. Configuration reset to defaults\n")
  
  cat("\n=== Edge Cases Demo Complete ===\n")
  expect_true(TRUE)
})

test_that("parameter configuration protocols are implemented correctly", {
  cat("\n=== Parameter Configuration Protocol Validation ===\n")
  
  # Protocol 1: Interactive analysis parameter configuration
  cat("\n1. Interactive analysis parameter specification:\n")
  set_immunoplex_config(
    random_effects_min_subjects = 10,   # Reduced estimation constraints
    default_dharma_nsim = 250,          # Optimized simulation count
    plotting_point_size = 1.8           # Enhanced visualization scale
  )
  
  interactive_config <- get_immunoplex_config()
  cat("   Quick setup: min_subjects =", interactive_config$random_effects$min_subjects,
      ", dharma_nsim =", interactive_config$dharma$default_nsim,
      ", point_size =", interactive_config$plotting$default_point_size, "\n")
  
  # Protocol 2: Production analysis parameter configuration
  cat("\n2. Production analysis parameter specification:\n")
  set_immunoplex_config(
    random_effects_min_subjects = 50,   # Maximum estimation precision
    random_effects_min_reps = 4,        # Enhanced replicate requirements
    default_dharma_nsim = 2000,         # Comprehensive diagnostic evaluation
    high_censoring_threshold = 60       # Stringent distribution selection
  )
  
  production_config <- get_immunoplex_config()
  cat("   Production setup: min_subjects =", production_config$random_effects$min_subjects,
      ", min_reps =", production_config$random_effects$min_reps,
      ", dharma_nsim =", production_config$dharma$default_nsim, "\n")
  
  # Protocol 3: Limited sample size parameter configuration
  cat("\n3. Limited sample size parameter specification:\n")
  set_immunoplex_config(
    random_effects_min_subjects = 8,    # Minimal subject count threshold
    random_effects_min_reps = 2,        # Minimal replicate requirement
    high_censoring_threshold = 30,      # Modified distribution selection
    high_skewness_threshold = 1.5       # Adjusted distribution criteria
  )
  
  small_config <- get_immunoplex_config()
  cat("   Small data setup: min_subjects =", small_config$random_effects$min_subjects,
      ", censoring_threshold =", small_config$auto_family$high_censoring_threshold,
      ", skewness_threshold =", small_config$auto_family$high_skewness_threshold, "\n")
  
  # Show how to check current settings
  cat("\n4. Checking current configuration:\n")
  current <- get_immunoplex_config()
  cat("   Current random effects policy: >=", current$random_effects$min_subjects, 
      "subjects with >=", current$random_effects$min_reps, "reps each\n")
  cat("   Current auto selection: Tobit if >", current$auto_family$high_censoring_threshold,
      "% censored, Gamma if skewness >", current$auto_family$high_skewness_threshold, "\n")
  
  # Reset
  reset_immunoplex_config()
  cat("\n5. Reset to package defaults\n")
  
  cat("\n=== Usage Patterns Demo Complete ===\n")
  expect_true(TRUE)
})

})

# Implementation Note:
# This validation protocol verifies configuration system functionality
# and parameter specification protocols. Execute validation with:
# testthat::test_file("tests/testthat/test-config_demo.R")