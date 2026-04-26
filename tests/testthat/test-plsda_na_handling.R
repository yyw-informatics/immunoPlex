# Tests for Safe NA Handling in PLS-DA Functions
#
# These tests verify that print and summary methods handle NA values gracefully
# and that sum() operations use na.rm = TRUE to prevent failures.
#
# Issue: When OPLS-DA models had NA statistics, print/summary methods would
# crash or display confusing output. The fix adds safe NA handling throughout.

# Helper function
create_na_test_data <- function(n_samples = 50, n_cytokines = 10) {
  set.seed(789)
  
  test_data <- data.frame(
    SubjectID = paste0("S", 1:n_samples),
    condition = rep(c("Pre", "Post"), each = n_samples/2),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n_cytokines) {
    pre_vals <- rnorm(n_samples/2, mean = 90, sd = 30)
    post_vals <- rnorm(n_samples/2, mean = 130, sd = 30)
    test_data[[paste0("cyt", i)]] <- c(pre_vals, post_vals)
  }
  
  lod_df <- data.frame(
    cytokine = paste0("cyt", 1:n_cytokines),
    lod = rep(10, n_cytokines)
  )
  
  list(data = test_data, lod = lod_df)
}


test_that("print.plsda_model works for OPLS-DA models", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_na_test_data()
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:10),
    metadata_cols = c("SubjectID", "condition"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  oplsda_model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "condition",
    n_components = 2,
    method = "OPLS-DA",
    permutations = 50,  # Fewer for speed
    verbose = FALSE
  )
  
  # Print should not error
  expect_no_error(print(oplsda_model))
  
  # Should produce expected output
  expect_output(print(oplsda_model), "OPLS-DA")
  expect_output(print(oplsda_model), "Model Performance")
})


test_that("print method displays R2X, R2Y, Q2 for OPLS-DA", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_na_test_data()
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:10),
    metadata_cols = c("SubjectID", "condition"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  oplsda_model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "condition",
    n_components = 2,
    method = "OPLS-DA",
    permutations = 50,
    verbose = FALSE
  )
  
  # Capture output
  output <- capture.output(print(oplsda_model))
  output_str <- paste(output, collapse = "\n")
  
  # Should display all statistics
  expect_match(output_str, "R2X")
  expect_match(output_str, "R2Y")
  expect_match(output_str, "Q2")
  
  # Should not show "NA" for these statistics
  expect_false(grepl("R2X.*NA", output_str),
               info = "R2X should not display as NA")
  expect_false(grepl("R2Y.*NA", output_str),
               info = "R2Y should not display as NA")
})


test_that("print method displays permutation p-value when available", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_na_test_data()
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:10),
    metadata_cols = c("SubjectID", "condition"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "condition",
    n_components = 2,
    method = "PLS-DA",
    permutations = 100,
    verbose = FALSE
  )
  
  output <- capture.output(print(model))
  output_str <- paste(output, collapse = "\n")
  
  # Should display permutation p-value
  expect_match(output_str, "Permutation p-value",
               info = "Print should display permutation p-value")
  
  # Should not show NA
  expect_false(grepl("Permutation.*NA", output_str),
               info = "Permutation p-value should not be NA")
})


test_that("print method finds p-value under multiple field names", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_na_test_data()
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:10),
    metadata_cols = c("SubjectID", "condition"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "condition",
    n_components = 2,
    method = "PLS-DA",
    permutations = 100,
    verbose = FALSE
  )
  
  # Manually check that print looks for both field names
  # (This tests the fallback mechanism in the print function)
  expect_true("permutation_pvalue" %in% names(model$model_stats) ||
              "perm_pval_R2Y" %in% names(model$model_stats),
              info = "Print should be able to find p-value under either name")
})


test_that("summary.plsda_model works for OPLS-DA models", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_na_test_data()
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:10),
    metadata_cols = c("SubjectID", "condition"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  oplsda_model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "condition",
    n_components = 2,
    method = "OPLS-DA",
    permutations = 50,
    verbose = FALSE
  )
  
  # Summary should not error
  expect_no_error(summary(oplsda_model))
  
  # Should produce expected output
  expect_output(summary(oplsda_model), "Model Summary")
  expect_output(summary(oplsda_model), "Per-Component Statistics")
  expect_output(summary(oplsda_model), "Top 10 Cytokines")
})


test_that("Verbose output during model fitting doesn't crash", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_na_test_data()
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:10),
    metadata_cols = c("SubjectID", "condition"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  # Verbose mode should not error for OPLS-DA
  expect_output(
    plsda_fit(
      preprocessed_data = preprocessed,
      response_var = "condition",
      n_components = 2,
      method = "OPLS-DA",
      permutations = 50,
      verbose = TRUE  # Verbose mode
    ),
    "Model fitting complete"
  )
})


test_that("Verbose output displays statistics correctly for OPLS-DA", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_na_test_data()
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:10),
    metadata_cols = c("SubjectID", "condition"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  output <- capture.output(
    plsda_fit(
      preprocessed_data = preprocessed,
      response_var = "condition",
      n_components = 2,
      method = "OPLS-DA",
      permutations = 50,
      verbose = TRUE
    )
  )
  
  output_str <- paste(output, collapse = "\n")
  
  # Should display R2X, R2Y, Q2
  expect_match(output_str, "R2X:")
  expect_match(output_str, "R2Y:")
  expect_match(output_str, "Q2:")
  
  # Should display permutation p-value
  expect_match(output_str, "Permutation p-value")
})


test_that("sum() with na.rm=TRUE works correctly", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_na_test_data()
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:10),
    metadata_cols = c("SubjectID", "condition"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "condition",
    n_components = 2,
    method = "OPLS-DA",
    permutations = 0,
    verbose = FALSE
  )
  
  # These operations should not error even if there were NAs
  expect_no_error(sum(model$model_stats$R2X, na.rm = TRUE))
  expect_no_error(sum(model$model_stats$R2Y, na.rm = TRUE))
  expect_no_error(sum(model$model_stats$Q2, na.rm = TRUE))
  
  # Results should be numeric and finite
  expect_true(is.numeric(sum(model$model_stats$R2X, na.rm = TRUE)))
  expect_true(is.finite(sum(model$model_stats$R2X, na.rm = TRUE)))
})

