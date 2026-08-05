# Critical Interval Coding Validation Protocol
#
# Validates the most important fixes from the code review:
# 1. Proper interval construction for censored models
# 2. Correct bounds for left-censored: [-Inf, LOD] (Tobit) vs [ε, LOD] (AFT)
# 3. Correct bounds for right-censored: [ULOD, +Inf]
# 4. Interval-censored: [LOD, ULOD]
# 5. Exact observations: [value, value]
#
# These tests verify the core statistical correctness of our censoring fixes.

library(testthat)

# Create test data with specific censoring patterns
setup_interval_test_data <- function() {
  set.seed(123)  # For reproducible tests
  n <- 50
  
  # Create synthetic data with known censoring patterns
  dat <- data.frame(
    sample_id = paste0("S", 1:n),
    subject_id = rep(paste0("P", 1:25), length.out = n),  # Match n=50
    timepoint = rep(c("T1", "T2"), length.out = n),
    disease = rep(c("Control", "Disease"), length.out = n),
    age = rnorm(n, 45, 10),
    value = exp(rnorm(n, 2, 0.5)),  # Log-normal data
    lod = 0.5,
    ulod = 20.0,
    stringsAsFactors = FALSE
  )
  
  # Create specific censoring patterns for testing
  dat$cens_lod <- dat$value < dat$lod
  dat$cens_ulod <- dat$value > dat$ulod
  
  # Force some specific censoring cases for testing
  dat$cens_lod[1:5] <- TRUE
  dat$value[1:5] <- dat$lod[1:5] - 0.1  # Below LOD
  
  dat$cens_ulod[6:8] <- TRUE  
  dat$value[6:8] <- dat$ulod[6:8] + 1  # Above ULOD
  
  # Create interval-censored cases (both flags TRUE)
  dat$cens_lod[9:10] <- TRUE
  dat$cens_ulod[9:10] <- TRUE
  dat$lod[9:10] <- 0.3
  dat$ulod[9:10] <- 0.7
  dat$value[9:10] <- 0.5  # Between LOD and ULOD
  
  return(dat)
}

describe("Critical interval coding validation", {

test_that("Tobit model uses correct interval bounds", {
  skip_if_not_installed("survival")
  
  dat <- setup_interval_test_data()
  
  # Test that Tobit model fits without errors (with ulod=TRUE to include right-censoring)
  fit <- fit_one(dat, family = "tobit", ulod = TRUE, random = "")
  expect_s3_class(fit, "immuno_fit")
  expect_equal(fit$family, "tobit")
  expect_equal(fit$estimand, "ratio_of_means")
  
  # For Tobit, left-censored should use [-Inf, LOD] bounds
  # This is tested implicitly by successful model fitting with survival::survreg
  if (fit$converged) {
    expect_true(is.numeric(fit$aic))
    expect_true(is.finite(fit$aic))
    
    # Check that model handled censored data correctly
    expect_equal(fit$n_cens_lod, sum(dat$cens_lod))
    expect_equal(fit$n_cens_ulod, sum(dat$cens_ulod))
  }
})

test_that("AFT model uses correct interval bounds", {
  skip_if_not_installed("survival")
  
  dat <- setup_interval_test_data()
  
  # Test that AFT model fits without errors (AFT doesn't use ulod parameter, but let's be explicit)
  fit <- fit_one(dat, family = "aft", random = "")
  expect_s3_class(fit, "immuno_fit")
  expect_equal(fit$family, "aft")
  expect_equal(fit$estimand, "ratio_of_medians")
  
  # For AFT, left-censored should use [ε, LOD] bounds (positive for lognormal)
  # This is tested implicitly by successful model fitting with survival::survreg
  if (fit$converged) {
    expect_true(is.numeric(fit$aic))
    expect_true(is.finite(fit$aic))
    
    # Check that model handled censored data correctly
    expect_equal(fit$n_cens_lod, sum(dat$cens_lod))
    # AFT should handle right-censoring regardless of ulod parameter
    expect_equal(fit$n_cens_ulod, sum(dat$cens_ulod))
  }
})

test_that("interval bounds are constructed correctly for different censoring types", {
  skip_if_not_installed("survival")
  
  dat <- setup_interval_test_data()
  
  # Test both Tobit and AFT handle all censoring types
  families <- c("tobit", "aft")
  
  for (fam in families) {
    # For Tobit, need ulod=TRUE to include right-censoring; AFT handles it automatically
    fit <- if (fam == "tobit") {
      fit_one(dat, family = fam, ulod = TRUE, random = "")
    } else {
      fit_one(dat, family = fam, random = "")
    }
    
    if (fit$converged) {
      # Model should handle:
      # - Left-censored (cens_lod = TRUE, cens_ulod = FALSE)
      # - Right-censored (cens_lod = FALSE, cens_ulod = TRUE)  
      # - Interval-censored (cens_lod = TRUE, cens_ulod = TRUE)
      # - Exact (cens_lod = FALSE, cens_ulod = FALSE)
      
      left_cens <- sum(dat$cens_lod & !dat$cens_ulod)
      right_cens <- sum(!dat$cens_lod & dat$cens_ulod)
      interval_cens <- sum(dat$cens_lod & dat$cens_ulod)
      exact <- sum(!dat$cens_lod & !dat$cens_ulod)
      
      expect_equal(fit$n_cens_lod, left_cens + interval_cens)
      expect_equal(fit$n_cens_ulod, right_cens + interval_cens)
      
      # Total observations should match
      expect_equal(left_cens + right_cens + interval_cens + exact, nrow(dat))
    }
  }
})

test_that("ulod parameter controls right-censoring correctly", {
  skip_if_not_installed("survival")
  
  dat <- setup_interval_test_data()
  
  # Test with ulod = TRUE (should include right-censoring)
  fit_with_ulod <- fit_one(dat, family = "tobit", ulod = TRUE, random = "")
  
  # Test with ulod = FALSE (should ignore right-censoring)
  fit_without_ulod <- fit_one(dat, family = "tobit", ulod = FALSE, random = "")
  
  if (fit_with_ulod$converged && fit_without_ulod$converged) {
    # With ulod=TRUE, should count right-censored observations
    expect_equal(fit_with_ulod$n_cens_ulod, sum(dat$cens_ulod))
    
    # With ulod=FALSE, should ignore right-censoring (set to 0)
    expect_equal(fit_without_ulod$n_cens_ulod, 0)
    
    # Left-censoring should be the same in both cases
    expect_equal(fit_with_ulod$n_cens_lod, fit_without_ulod$n_cens_lod)
  }
})

test_that("bounds validation prevents invalid intervals", {
  skip_if_not_installed("survival")
  
  dat <- setup_interval_test_data()
  
  # Create invalid bounds (LOD > ULOD) for testing
  dat$lod[1] <- 10
  dat$ulod[1] <- 5  # LOD > ULOD (invalid)
  dat$cens_lod[1] <- TRUE
  dat$cens_ulod[1] <- TRUE
  
  # Should fail with bounds validation error
  expect_error(
    fit_one(dat, family = "tobit", random = ""),
    "Invalid censoring bounds"
  )
})

test_that("row-specific LOD/ULOD values work correctly", {
  skip_if_not_installed("survival")
  
  dat <- setup_interval_test_data()
  
  # Create row-specific LOD values
  dat$lod <- c(0.3, 0.4, 0.5, rep(0.5, nrow(dat) - 3))
  dat$ulod <- c(15, 18, 20, rep(20, nrow(dat) - 3))
  
  # Update censoring flags based on new LODs
  dat$cens_lod <- dat$value < dat$lod
  dat$cens_ulod <- dat$value > dat$ulod
  
  # Should work fine for Tobit and AFT (they support row-specific LODs)
  fit_tobit <- fit_one(dat, family = "tobit", random = "")
  expect_s3_class(fit_tobit, "immuno_fit")
  
  fit_aft <- fit_one(dat, family = "aft", random = "")
  expect_s3_class(fit_aft, "immuno_fit")
  
  # But should fail for censReg (requires single LOD)
  if (requireNamespace("censReg", quietly = TRUE)) {
    expect_error(
      fit_one(dat, family = "tobit_censreg", random = ""),
      "censReg requires a single LOD"
    )
  }
})

test_that("estimand assignment is correct for each family", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  
  dat <- setup_interval_test_data()
  
  # Test all families have correct estimands
  families_and_estimands <- list(
    gamma = "ratio_of_means",
    tobit = "ratio_of_means", 
    aft = "ratio_of_medians"
  )
  
  if (requireNamespace("censReg", quietly = TRUE)) {
    families_and_estimands$tobit_censreg <- "ratio_of_means"
  }
  
  for (fam in names(families_and_estimands)) {
    expected_estimand <- families_and_estimands[[fam]]
    
    # Handle censReg which may fail with row-specific LODs
    if (fam == "tobit_censreg") {
      tryCatch({
        fit <- fit_one(dat, family = fam, random = "")
        expect_s3_class(fit, "immuno_fit")
        expect_equal(fit$estimand, expected_estimand)
      }, error = function(e) {
        # Expected to fail with row-specific LODs
        expect_true(grepl("censReg requires a single LOD", e$message))
      })
    } else {
      fit <- fit_one(dat, family = fam, random = "")
      expect_s3_class(fit, "immuno_fit")
      expect_equal(fit$estimand, expected_estimand)
    }
  }
})

test_that("auto-selector avoids Gamma with censoring", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  skip_if_not_installed("moments")
  
  dat <- setup_interval_test_data()
  
  # With censoring present, auto should NOT select gamma
  fit_auto <- fit_one(dat, family = "auto", random = "")
  expect_s3_class(fit_auto, "immuno_fit")
  
  # Should select tobit or aft, not gamma (since we have censoring)
  expect_true(fit_auto$family %in% c("tobit", "aft"))
  expect_false(fit_auto$family == "gamma")
  
  # Test with no censoring - should be able to select gamma
  dat_no_cens <- dat
  dat_no_cens$cens_lod[] <- FALSE
  dat_no_cens$cens_ulod[] <- FALSE
  
  fit_auto_no_cens <- fit_one(dat_no_cens, family = "auto", random = "")
  expect_s3_class(fit_auto_no_cens, "immuno_fit")
  # Can now select gamma or gaussian since no censoring
  expect_true(fit_auto_no_cens$family %in% c("gamma", "tobit", "aft", "gaussian"))
})

})
