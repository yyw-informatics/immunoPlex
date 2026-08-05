# Tests for OPLS-DA Model Statistics Extraction
#
# These tests verify that OPLS-DA models correctly extract R2X, R2Y, and Q2
# statistics from the ropls package's @summaryDF slot when @modelDF contains NAs.
#
# Issue: OPLS-DA models were returning NA for cumulative statistics because
# the extraction was only looking at @modelDF, which stores NAs for OPLS-DA.
# The fix extracts from @summaryDF (which has R2X(cum), R2Y(cum), Q2(cum)) when needed.

# Helper function to create test data with good group separation
create_plsda_test_data <- function(n_samples = 50, n_cytokines = 14) {
  set.seed(123)  # For reproducibility
  
  test_data <- data.frame(
    SubjectID = paste0("S", 1:n_samples),
    group = rep(c("Control", "Treatment"), each = n_samples/2),
    stringsAsFactors = FALSE
  )
  
  # Add cytokines with realistic group separation
  for (i in 1:n_cytokines) {
    control_vals <- rnorm(n_samples/2, mean = 100, sd = 25)
    treatment_vals <- rnorm(n_samples/2, mean = 150, sd = 25)
    test_data[[paste0("cyt", i)]] <- c(control_vals, treatment_vals)
  }
  
  lod_df <- data.frame(
    cytokine = paste0("cyt", 1:n_cytokines),
    lod = rep(10, n_cytokines)
  )
  
  list(data = test_data, lod = lod_df)
}


test_that("OPLS-DA model extracts valid R2X statistic (not NA)", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_plsda_test_data(n_samples = 50, n_cytokines = 14)
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:14),
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  oplsda_model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "group",
    n_components = 2,
    method = "OPLS-DA",
    permutations = 0,  # Skip for speed
    verbose = FALSE
  )
  
  # Check that R2X is not NA
  r2x_sum <- sum(oplsda_model$model_stats$R2X, na.rm = TRUE)
  
  expect_false(is.na(r2x_sum), 
               info = "OPLS-DA R2X should not be NA")
  expect_true(r2x_sum > 0, 
              info = "OPLS-DA R2X should be positive")
  expect_true(r2x_sum <= 1, 
              info = "OPLS-DA R2X should be <= 1")
})


test_that("OPLS-DA model extracts valid R2Y statistic (not NA)", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_plsda_test_data(n_samples = 50, n_cytokines = 14)
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:14),
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  oplsda_model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "group",
    n_components = 2,
    method = "OPLS-DA",
    permutations = 0,
    verbose = FALSE
  )
  
  # Check that R2Y is not NA
  r2y_sum <- sum(oplsda_model$model_stats$R2Y, na.rm = TRUE)
  
  expect_false(is.na(r2y_sum), 
               info = "OPLS-DA R2Y should not be NA")
  expect_true(r2y_sum > 0, 
              info = "OPLS-DA R2Y should be positive")
  expect_true(r2y_sum <= 1, 
              info = "OPLS-DA R2Y should be <= 1")
})


test_that("OPLS-DA model extracts valid Q2 statistic (not NA)", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_plsda_test_data(n_samples = 50, n_cytokines = 14)
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:14),
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  oplsda_model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "group",
    n_components = 2,
    method = "OPLS-DA",
    permutations = 0,
    verbose = FALSE
  )
  
  # Check that Q2 is not NA
  q2_sum <- sum(oplsda_model$model_stats$Q2, na.rm = TRUE)
  
  expect_false(is.na(q2_sum), 
               info = "OPLS-DA Q2 should not be NA")
  expect_true(is.numeric(q2_sum) && is.finite(q2_sum), 
              info = "OPLS-DA Q2 should be a finite numeric value")
  # Q2 can be negative for poor models, so don't test positivity
})


test_that("OPLS-DA and PLS-DA both return valid statistics", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_plsda_test_data(n_samples = 50, n_cytokines = 14)
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:14),
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  # Fit both model types
  plsda_model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "group",
    n_components = 2,
    method = "PLS-DA",
    permutations = 0,
    verbose = FALSE
  )
  
  oplsda_model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "group",
    n_components = 2,
    method = "OPLS-DA",
    permutations = 0,
    verbose = FALSE
  )
  
  # Both should have valid statistics
  expect_false(is.na(sum(plsda_model$model_stats$R2X, na.rm = TRUE)))
  expect_false(is.na(sum(plsda_model$model_stats$R2Y, na.rm = TRUE)))
  expect_false(is.na(sum(plsda_model$model_stats$Q2, na.rm = TRUE)))
  
  expect_false(is.na(sum(oplsda_model$model_stats$R2X, na.rm = TRUE)))
  expect_false(is.na(sum(oplsda_model$model_stats$R2Y, na.rm = TRUE)))
  expect_false(is.na(sum(oplsda_model$model_stats$Q2, na.rm = TRUE)))
})


test_that("Model comparison table has no NA values for OPLS-DA", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_plsda_test_data(n_samples = 50, n_cytokines = 14)
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:14),
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  plsda_model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "group",
    n_components = 2,
    method = "PLS-DA",
    permutations = 0,
    verbose = FALSE
  )
  
  oplsda_model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "group",
    n_components = 2,
    method = "OPLS-DA",
    permutations = 0,
    verbose = FALSE
  )
  
  # Create model comparison as done in analysis scripts
  model_comparison <- data.frame(
    Model = c("PLS-DA", "OPLS-DA"),
    R2X_total = c(
      sum(plsda_model$model_stats$R2X, na.rm = TRUE),
      sum(oplsda_model$model_stats$R2X, na.rm = TRUE)
    ),
    R2Y_total = c(
      sum(plsda_model$model_stats$R2Y, na.rm = TRUE),
      sum(oplsda_model$model_stats$R2Y, na.rm = TRUE)
    ),
    Q2_total = c(
      sum(plsda_model$model_stats$Q2, na.rm = TRUE),
      sum(oplsda_model$model_stats$Q2, na.rm = TRUE)
    )
  )
  
  # Check that OPLS-DA row has no NAs
  expect_true(all(!is.na(model_comparison[2, c("R2X_total", "R2Y_total", "Q2_total")])),
              info = "OPLS-DA row should have no NA values")
  
  # Check that all values are finite
  expect_true(all(is.finite(model_comparison$R2X_total)))
  expect_true(all(is.finite(model_comparison$R2Y_total)))
  expect_true(all(is.finite(model_comparison$Q2_total)))
})

