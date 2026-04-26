# Tests for plsda_preprocess function

test_that("plsda_preprocess validates inputs correctly", {
  # Create minimal test data
  test_data <- data.frame(
    SubjectID = paste0("S", 1:10),
    group = rep(c("A", "B"), each = 5),
    cyt1 = rnorm(10, 100, 20),
    cyt2 = rnorm(10, 50, 10),
    stringsAsFactors = FALSE
  )
  
  lod_df <- data.frame(
    cytokine = c("cyt1", "cyt2"),
    lod = c(10, 5)
  )
  
  # Should work with valid inputs
  expect_no_error(
    plsda_preprocess(
      data = test_data,
      cytokine_cols = c("cyt1", "cyt2"),
      metadata_cols = c("SubjectID", "group"),
      lod_lookup = lod_df,
      verbose = FALSE
    )
  )
  
  # Should error with missing cytokine column
  expect_error(
    plsda_preprocess(
      data = test_data,
      cytokine_cols = c("cyt1", "cyt_missing"),
      metadata_cols = c("SubjectID", "group"),
      lod_lookup = lod_df,
      verbose = FALSE
    ),
    "Cytokine columns not found"
  )
  
  # Should error with missing metadata column
  expect_error(
    plsda_preprocess(
      data = test_data,
      cytokine_cols = c("cyt1", "cyt2"),
      metadata_cols = c("SubjectID", "missing_meta"),
      lod_lookup = lod_df,
      verbose = FALSE
    ),
    "Metadata columns not found"
  )
  
  # Should warn (not error) with missing LOD values - function now allows processing
  # without LOD values for some cytokines
  lod_incomplete <- data.frame(
    cytokine = c("cyt1"),
    lod = c(10)
  )
  
  # This should work now, just without LOD substitution for cyt2
  result <- plsda_preprocess(
    data = test_data,
    cytokine_cols = c("cyt1", "cyt2"),
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = lod_incomplete,
    verbose = FALSE
  )
  expect_s3_class(result, "plsda_preprocessed")
  
  # Should error with invalid LOD method
  expect_error(
    plsda_preprocess(
      data = test_data,
      cytokine_cols = c("cyt1", "cyt2"),
      metadata_cols = c("SubjectID", "group"),
      lod_lookup = lod_df,
      lod_method = "invalid_method",
      verbose = FALSE
    ),
    "lod_method must be one of"
  )
})


test_that("plsda_preprocess handles string cytokine values correctly", {
  # Create test data with string values (like "< 12.8")
  test_data <- data.frame(
    SubjectID = paste0("S", 1:10),
    group = rep(c("A", "B"), each = 5),
    cyt1 = c("< 10", "15", "20", "< 10", "25", "30", "< 10", "35", "40", "45"),
    cyt2 = c("< 5", "< 5", "8", "10", "12", "< 5", "15", "18", "20", "22"),
    stringsAsFactors = FALSE
  )
  
  lod_df <- data.frame(
    cytokine = c("cyt1", "cyt2"),
    lod = c(10, 5)
  )
  
  result <- plsda_preprocess(
    data = test_data,
    cytokine_cols = c("cyt1", "cyt2"),
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = lod_df,
    lod_method = "half",
    log_transform = FALSE,
    scale_data = FALSE,
    verbose = FALSE
  )
  
  # Check that string values were converted
  expect_true(is.numeric(result$expression[, "cyt1"]))
  expect_true(is.numeric(result$expression[, "cyt2"]))
  
  # Check that censored values (originally "< X") were substituted with LOD/2
  # Sample 1 had "< 10" for cyt1, should now be 5
  expect_equal(result$expression[1, "cyt1"], 5)
  
  # Check that non-censored values are preserved
  # Sample 2 had "15" for cyt1, should still be 15
  expect_equal(result$expression[2, "cyt1"], 15)
})


test_that("plsda_preprocess applies LOD methods correctly", {
  test_data <- data.frame(
    SubjectID = paste0("S", 1:5),
    group = c("A", "A", "B", "B", "B"),
    cyt1 = c(5, 15, 25, 35, 45),
    stringsAsFactors = FALSE
  )
  
  lod_df <- data.frame(
    cytokine = "cyt1",
    lod = 10
  )
  
  # Test "half" method
  result_half <- plsda_preprocess(
    data = test_data,
    cytokine_cols = "cyt1",
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = lod_df,
    lod_method = "half",
    log_transform = FALSE,
    scale_data = FALSE,
    verbose = FALSE
  )
  expect_equal(result_half$expression[1, "cyt1"], 5)  # LOD/2 = 5
  
  # Test "zero" method
  result_zero <- plsda_preprocess(
    data = test_data,
    cytokine_cols = "cyt1",
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = lod_df,
    lod_method = "zero",
    log_transform = FALSE,
    scale_data = FALSE,
    verbose = FALSE
  )
  expect_equal(result_zero$expression[1, "cyt1"], 0)
  
  # Test "lod" method
  result_lod <- plsda_preprocess(
    data = test_data,
    cytokine_cols = "cyt1",
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = lod_df,
    lod_method = "lod",
    log_transform = FALSE,
    scale_data = FALSE,
    verbose = FALSE
  )
  expect_equal(result_lod$expression[1, "cyt1"], 10)
  
  # Test "sqrt" method
  result_sqrt <- plsda_preprocess(
    data = test_data,
    cytokine_cols = "cyt1",
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = lod_df,
    lod_method = "sqrt",
    log_transform = FALSE,
    scale_data = FALSE,
    verbose = FALSE
  )
  expect_equal(result_sqrt$expression[1, "cyt1"], sqrt(10))
})


test_that("plsda_preprocess applies log transformation correctly", {
  test_data <- data.frame(
    SubjectID = paste0("S", 1:5),
    group = c("A", "A", "B", "B", "B"),
    cyt1 = c(0, 1, 3, 7, 15),
    stringsAsFactors = FALSE
  )
  
  lod_df <- data.frame(
    cytokine = "cyt1",
    lod = 0.5
  )
  
  result <- plsda_preprocess(
    data = test_data,
    cytokine_cols = "cyt1",
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = lod_df,
    lod_method = "half",
    log_transform = TRUE,
    scale_data = FALSE,
    verbose = FALSE
  )
  
  # Check log2(x + 1) transformation
  # Note: value 0 is < LOD (0.5), so it gets substituted with LOD/2 = 0.25
  expect_equal(result$expression[1, "cyt1"], log2(0.25 + 1))
  expect_equal(result$expression[2, "cyt1"], log2(1 + 1))
  expect_equal(result$expression[5, "cyt1"], log2(15 + 1))
})


test_that("plsda_preprocess applies z-score scaling correctly", {
  test_data <- data.frame(
    SubjectID = paste0("S", 1:10),
    group = rep(c("A", "B"), each = 5),
    cyt1 = seq(10, 100, by = 10),
    cyt2 = seq(5, 50, by = 5),
    stringsAsFactors = FALSE
  )
  
  lod_df <- data.frame(
    cytokine = c("cyt1", "cyt2"),
    lod = c(5, 2)
  )
  
  result <- plsda_preprocess(
    data = test_data,
    cytokine_cols = c("cyt1", "cyt2"),
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = lod_df,
    lod_method = "half",
    log_transform = FALSE,
    scale_data = TRUE,
    verbose = FALSE
  )
  
  # Check that mean is approximately 0 and sd is approximately 1 for each cytokine
  expect_equal(mean(result$expression[, "cyt1"]), 0, tolerance = 1e-10)
  expect_equal(sd(result$expression[, "cyt1"]), 1, tolerance = 1e-10)
  expect_equal(mean(result$expression[, "cyt2"]), 0, tolerance = 1e-10)
  expect_equal(sd(result$expression[, "cyt2"]), 1, tolerance = 1e-10)
})


test_that("plsda_preprocess returns correct structure", {
  test_data <- data.frame(
    SubjectID = paste0("S", 1:5),
    group = c("A", "A", "B", "B", "B"),
    cyt1 = c(5, 15, 25, 35, 45),
    stringsAsFactors = FALSE
  )
  
  lod_df <- data.frame(
    cytokine = "cyt1",
    lod = 10
  )
  
  result <- plsda_preprocess(
    data = test_data,
    cytokine_cols = "cyt1",
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = lod_df,
    verbose = FALSE
  )
  
  # Check class
  expect_s3_class(result, "plsda_preprocessed")
  expect_type(result, "list")
  
  # Check components
  expect_true("expression" %in% names(result))
  expect_true("metadata" %in% names(result))
  expect_true("preprocessing_info" %in% names(result))
  expect_true("sample_ids" %in% names(result))
  
  # Check expression is a matrix
  expect_true(is.matrix(result$expression))
  expect_equal(nrow(result$expression), 5)
  expect_equal(ncol(result$expression), 1)
  
  # Check metadata is a data frame
  expect_true(is.data.frame(result$metadata))
  expect_equal(nrow(result$metadata), 5)
  
  # Check preprocessing info
  expect_type(result$preprocessing_info, "list")
  expect_equal(result$preprocessing_info$lod_method, "half")
})


test_that("plsda_preprocess averages replicates correctly", {
  test_data <- data.frame(
    SubjectID = c("S1", "S1", "S2", "S2", "S3"),
    group = c("A", "A", "B", "B", "B"),
    cyt1 = c(10, 12, 20, 22, 30),
    cyt2 = c(5, 7, 15, 17, 25),
    stringsAsFactors = FALSE
  )
  
  lod_df <- data.frame(
    cytokine = c("cyt1", "cyt2"),
    lod = c(5, 2)
  )
  
  result <- plsda_preprocess(
    data = test_data,
    cytokine_cols = c("cyt1", "cyt2"),
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = lod_df,
    lod_method = "half",
    log_transform = FALSE,
    scale_data = FALSE,
    replicate_col = "SubjectID",
    verbose = FALSE
  )
  
  # Should have 3 samples after averaging (S1, S2, S3)
  expect_equal(nrow(result$expression), 3)
  
  # Check that S1's replicates were averaged
  # S1 had values 10 and 12 for cyt1, mean = 11
  expect_equal(result$expression["S1", "cyt1"], 11)
  
  # Check that S2's replicates were averaged
  # S2 had values 20 and 22 for cyt1, mean = 21
  expect_equal(result$expression["S2", "cyt1"], 21)
})


test_that("plsda_preprocess handles missing LOD values with fallback", {
  # Create test data where cyt2 has no LOD value but has NA values
  test_data <- data.frame(
    SubjectID = paste0("S", 1:10),
    group = rep(c("A", "B"), each = 5),
    cyt1 = c(5, 15, 25, 35, 45, 55, 65, 75, 85, 95),
    cyt2 = c(NA, 100, 200, NA, 300, 400, 500, NA, 600, 700),  # Has NAs and no LOD
    stringsAsFactors = FALSE
  )
  
  # Only provide LOD for cyt1, not cyt2
  lod_df <- data.frame(
    cytokine = "cyt1",
    lod = 10
  )
  
  # Should work and use minimum value as empirical LOD for cyt2
  result <- plsda_preprocess(
    data = test_data,
    cytokine_cols = c("cyt1", "cyt2"),
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = lod_df,
    lod_method = "half",
    log_transform = TRUE,
    scale_data = TRUE,
    verbose = FALSE
  )
  
  # Check that result is valid
  expect_s3_class(result, "plsda_preprocessed")
  
  # Check that no non-finite values remain (the key test!)
  expect_equal(sum(!is.finite(result$expression)), 0,
               info = "Should have no non-finite values after fallback LOD handling")
  
  # Check that no NA values remain
  expect_equal(sum(is.na(result$expression)), 0,
               info = "Should have no NA values after fallback LOD handling")
  
  # Check that all samples are present
  expect_equal(nrow(result$expression), 10)
  
  # Check that both cytokines are present
  expect_equal(ncol(result$expression), 2)
  expect_true("cyt1" %in% colnames(result$expression))
  expect_true("cyt2" %in% colnames(result$expression))
})


test_that("plsda_preprocess fallback uses minimum value correctly", {
  # Create test data with known minimum value
  test_data <- data.frame(
    SubjectID = paste0("S", 1:8),
    group = rep(c("A", "B"), each = 4),
    cyt1 = c(NA, 50, 100, 150, NA, 200, 250, 300),  # min = 50, has NAs
    stringsAsFactors = FALSE
  )
  
  # No LOD provided for cyt1
  lod_df <- data.frame(
    cytokine = character(0),
    lod = numeric(0)
  )
  
  result <- plsda_preprocess(
    data = test_data,
    cytokine_cols = "cyt1",
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = lod_df,
    lod_method = "half",
    log_transform = FALSE,
    scale_data = FALSE,
    verbose = FALSE
  )
  
  # The minimum positive value is 50, so empirical LOD = 50
  # NAs are converted to 0, which is < 50, so they get substituted with LOD/2 = 25
  expect_equal(result$expression[1, "cyt1"], 25,
               info = "NA should be substituted with half of minimum value (50/2 = 25)")
  expect_equal(result$expression[5, "cyt1"], 25,
               info = "NA should be substituted with half of minimum value (50/2 = 25)")
  
  # Non-NA values should be preserved
  expect_equal(result$expression[2, "cyt1"], 50)
  expect_equal(result$expression[3, "cyt1"], 100)
})


test_that("plsda_preprocess handles character LOD values", {
  # Test that LOD values stored as characters in CSV are handled
  test_data <- data.frame(
    SubjectID = paste0("S", 1:5),
    group = c("A", "A", "B", "B", "B"),
    cyt1 = c(5, 15, 25, 35, 45),
    stringsAsFactors = FALSE
  )
  
  # LOD values as character (as they come from CSV)
  lod_df <- data.frame(
    cytokine = "cyt1",
    lod = "10",  # Character, not numeric
    stringsAsFactors = FALSE
  )
  
  # Should work and convert character LOD to numeric
  result <- plsda_preprocess(
    data = test_data,
    cytokine_cols = "cyt1",
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = lod_df,
    lod_method = "half",
    log_transform = FALSE,
    scale_data = FALSE,
    verbose = FALSE
  )
  
  # Should substitute value 5 (< LOD of 10) with 5 (LOD/2)
  expect_equal(result$expression[1, "cyt1"], 5)
  
  # Should preserve non-censored values
  expect_equal(result$expression[2, "cyt1"], 15)
})


test_that("plsda_preprocess handles NA in character LOD values", {
  # Test handling when LOD is "N/A" or similar in CSV
  test_data <- data.frame(
    SubjectID = paste0("S", 1:6),
    group = rep(c("A", "B"), each = 3),
    cyt1 = c(10, 20, 30, 40, 50, 60),
    cyt2 = c(NA, 100, 200, NA, 300, 400),
    stringsAsFactors = FALSE
  )
  
  # cyt1 has numeric LOD, cyt2 has "N/A" (as in your actual data)
  lod_df <- data.frame(
    cytokine = c("cyt1", "cyt2"),
    lod = c("10", "N/A"),
    stringsAsFactors = FALSE
  )
  
  # Should work with fallback for cyt2
  result <- plsda_preprocess(
    data = test_data,
    cytokine_cols = c("cyt1", "cyt2"),
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = lod_df,
    lod_method = "half",
    log_transform = TRUE,
    scale_data = FALSE,
    verbose = FALSE
  )
  
  # Should have no non-finite values
  expect_equal(sum(!is.finite(result$expression)), 0)
  expect_equal(sum(is.na(result$expression)), 0)
  
  # cyt1 value 10 should be right at LOD, preserved as is
  expect_equal(result$expression[1, "cyt1"], log2(10 + 1))
  
  # cyt2 NA should be substituted with half of minimum (100/2 = 50)
  expect_equal(result$expression[1, "cyt2"], log2(50 + 1))
})


test_that("plsda_preprocess print and summary methods work", {
  test_data <- data.frame(
    SubjectID = paste0("S", 1:5),
    group = c("A", "A", "B", "B", "B"),
    cyt1 = c(5, 15, 25, 35, 45),
    stringsAsFactors = FALSE
  )
  
  lod_df <- data.frame(
    cytokine = "cyt1",
    lod = 10
  )
  
  result <- plsda_preprocess(
    data = test_data,
    cytokine_cols = "cyt1",
    metadata_cols = c("SubjectID", "group"),
    lod_lookup = lod_df,
    verbose = FALSE
  )
  
  # Print should not error
  expect_output(print(result), "PLS-DA Preprocessed Data")
  expect_output(print(result), "Samples: 5")
  
  # Summary should not error
  expect_output(summary(result), "PLS-DA Preprocessed Data Summary")
  expect_output(summary(result), "Expression matrix summary")
})
