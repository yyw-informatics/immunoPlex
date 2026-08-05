# Tests for Permutation P-value Field Naming
#
# These tests verify that permutation p-values are stored with user-friendly
# field names that users can easily access in their analysis scripts.
#
# Issue: The function stored p-values as 'perm_pval_R2Y', but user scripts
# expected 'permutation_pvalue'. The fix stores both names for convenience
# and backwards compatibility.

# Helper function to create test data
create_perm_test_data <- function(n_samples = 50, n_cytokines = 14) {
  set.seed(456)  # Different seed from other tests
  
  test_data <- data.frame(
    SubjectID = paste0("S", 1:n_samples),
    group = rep(c("Baseline", "Followup"), each = n_samples/2),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n_cytokines) {
    baseline_vals <- rnorm(n_samples/2, mean = 110, sd = 20)
    followup_vals <- rnorm(n_samples/2, mean = 140, sd = 20)
    test_data[[paste0("cyt", i)]] <- c(baseline_vals, followup_vals)
  }
  
  lod_df <- data.frame(
    cytokine = paste0("cyt", 1:n_cytokines),
    lod = rep(10, n_cytokines)
  )
  
  list(data = test_data, lod = lod_df)
}


test_that("Permutation p-value stored with user-friendly field name", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_perm_test_data()
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:14),
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "group",
    n_components = 2,
    method = "PLS-DA",
    permutations = 100,
    verbose = FALSE
  )
  
  # Check that user-friendly field exists
  expect_true("permutation_pvalue" %in% names(model$model_stats),
              info = "permutation_pvalue field should exist for easy user access")
})


test_that("Permutation p-value for Q2 also stored with user-friendly name", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_perm_test_data()
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:14),
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "group",
    n_components = 2,
    method = "PLS-DA",
    permutations = 100,
    verbose = FALSE
  )
  
  # Check that Q2 p-value field exists
  expect_true("permutation_pvalue_Q2" %in% names(model$model_stats),
              info = "permutation_pvalue_Q2 field should exist")
})


test_that("Technical field names maintained for backwards compatibility", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_perm_test_data()
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:14),
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "group",
    n_components = 2,
    method = "PLS-DA",
    permutations = 100,
    verbose = FALSE
  )
  
  # Check that technical names still exist
  expect_true("perm_pval_R2Y" %in% names(model$model_stats),
              info = "perm_pval_R2Y (technical name) should exist for backwards compatibility")
  expect_true("perm_pval_Q2" %in% names(model$model_stats),
              info = "perm_pval_Q2 (technical name) should exist for backwards compatibility")
})


test_that("User-friendly and technical names point to same values", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_perm_test_data()
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:14),
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "group",
    n_components = 2,
    method = "PLS-DA",
    permutations = 100,
    verbose = FALSE
  )
  
  # Both names should be identical
  expect_equal(model$model_stats$permutation_pvalue[1],
               model$model_stats$perm_pval_R2Y[1],
               info = "permutation_pvalue and perm_pval_R2Y should be identical")
  
  expect_equal(model$model_stats$permutation_pvalue_Q2[1],
               model$model_stats$perm_pval_Q2[1],
               info = "permutation_pvalue_Q2 and perm_pval_Q2 should be identical")
})


test_that("Permutation p-values are in valid range [0, 1]", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_perm_test_data()
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:14),
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "group",
    n_components = 2,
    method = "PLS-DA",
    permutations = 100,
    verbose = FALSE
  )
  
  # Check R2Y p-value range
  pval_r2y <- model$model_stats$permutation_pvalue[1]
  expect_true(pval_r2y >= 0 && pval_r2y <= 1,
              info = "Permutation p-value for R2Y should be between 0 and 1")
  
  # Check Q2 p-value range
  pval_q2 <- model$model_stats$permutation_pvalue_Q2[1]
  expect_true(pval_q2 >= 0 && pval_q2 <= 1,
              info = "Permutation p-value for Q2 should be between 0 and 1")
})


test_that("OPLS-DA models also have permutation p-value fields", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_perm_test_data()
  
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
    permutations = 100,
    verbose = FALSE
  )
  
  # OPLS-DA should also have user-friendly field names
  expect_true("permutation_pvalue" %in% names(oplsda_model$model_stats),
              info = "OPLS-DA should have permutation_pvalue field")
  
  pval <- oplsda_model$model_stats$permutation_pvalue[1]
  expect_false(is.na(pval),
               info = "OPLS-DA permutation p-value should not be NA")
  expect_true(pval >= 0 && pval <= 1,
              info = "OPLS-DA permutation p-value should be in valid range")
})


test_that("Models without permutations don't have p-value fields", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_perm_test_data()
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:14),
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "group",
    n_components = 2,
    method = "PLS-DA",
    permutations = 0,  # No permutations
    verbose = FALSE
  )
  
  # Should not have permutation fields when permutations = 0
  expect_false("permutation_pvalue" %in% names(model$model_stats),
               info = "permutation_pvalue should not exist when permutations = 0")
  expect_false("perm_pval_R2Y" %in% names(model$model_stats),
               info = "perm_pval_R2Y should not exist when permutations = 0")
})


test_that("Analysis script pattern for accessing p-value works", {
  skip_if_not_installed("ropls")
  
  test_obj <- create_perm_test_data()
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:14),
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "group",
    n_components = 2,
    method = "PLS-DA",
    permutations = 100,
    verbose = FALSE
  )
  
  # Test the exact pattern used in analysis scripts
  pval_result <- ifelse(
    "permutation_pvalue" %in% names(model$model_stats) && 
    !is.null(model$model_stats$permutation_pvalue) &&
    !is.na(model$model_stats$permutation_pvalue[1]),
    sprintf("%.4f", model$model_stats$permutation_pvalue[1]),
    "Not computed"
  )
  
  # Should successfully return a p-value, not "Not computed"
  expect_false(pval_result == "Not computed",
               info = "Analysis script should successfully access permutation p-value")
  
  # Parse and validate
  pval_numeric <- as.numeric(pval_result)
  expect_true(!is.na(pval_numeric))
  expect_true(pval_numeric >= 0 && pval_numeric <= 1)
})

