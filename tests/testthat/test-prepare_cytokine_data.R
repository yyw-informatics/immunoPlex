# Analyte Data Processing and Validation Protocol
#
# Validates data transformation and integrity:
# 1. Data Structure Processing:
#    - Matrix transformation
#    - Feature extraction
#    - Metadata integration
#
# 2. Analyte-Specific Validation:
#    - Detection limit verification
#    - Measurement consistency
#    - Cross-reference validation
#
# 3. System Integration:
#    - Statistical model compatibility
#    - Data format verification
#    - Pipeline integration
#
# 4. Quality Assessment:
#    - Data completeness
#    - Reference integrity
#    - Constraint validation
#
# 5. Implementation Verification:
#    - Protocol documentation
#    - Workflow validation
#    - Error handling assessment

library(testthat)

describe("Analyte data processing functionality", {

test_that("core data transformation protocol executes correctly", {
  data("immunoplex_example", package = "immunoPlex")
  
  # Get the first available cytokine with LOD
  lod_lookup <- immunoplex_example$lod_lookup
  available_cytokines <- lod_lookup$cytokine[!is.na(lod_lookup$lod)]
  
  if (length(available_cytokines) == 0) {
    skip("No cytokines with LOD information found in test data")
  }
  
  target_cyto <- available_cytokines[1]
  
  dat <- prepare_cytokine_data(target_cyto)
  
  # Check basic structure
  expect_true(is.data.frame(dat))
  expect_true(all(c("sample_id", "subject_id", "value", "timepoint", "disease", 
                   "age", "lod", "ulod", "cens_lod", "cens_ulod") %in% names(dat)))
  
  # Check data types
  expect_type(dat$value, "double")
  expect_type(dat$lod, "double")
  expect_type(dat$cens_lod, "logical")
  expect_type(dat$cens_ulod, "logical")
  
  # Check censoring logic
  expect_true(all(dat$cens_lod == (dat$value < dat$lod)))
  if (!all(is.na(dat$ulod))) {
    expect_true(all(dat$cens_ulod == (!is.na(dat$ulod) & dat$value > dat$ulod)))
  }
  
  # Check attributes
  expect_equal(attr(dat, "cytokine"), target_cyto)
  expect_type(attr(dat, "n_censored"), "integer")
  expect_type(attr(dat, "pct_censored"), "double")
})

test_that("system processes multiple analytes appropriately", {
  data("immunoplex_example", package = "immunoPlex")
  
  lod_lookup <- immunoplex_example$lod_lookup
  available_cytokines <- lod_lookup$cytokine[!is.na(lod_lookup$lod)]
  
  if (length(available_cytokines) < 2) {
    skip("Need at least 2 cytokines with LOD for this test")
  }
  
  # Test first two cytokines
  dat1 <- prepare_cytokine_data(available_cytokines[1])
  dat2 <- prepare_cytokine_data(available_cytokines[2])
  
  # Should have different cytokine attributes
  expect_equal(attr(dat1, "cytokine"), available_cytokines[1])
  expect_equal(attr(dat2, "cytokine"), available_cytokines[2])
  
  # Should have different LOD values (if they differ in lookup)
  lod1 <- lod_lookup$lod[lod_lookup$cytokine == available_cytokines[1]]
  lod2 <- lod_lookup$lod[lod_lookup$cytokine == available_cytokines[2]]
  
  expect_equal(unique(dat1$lod), lod1)
  expect_equal(unique(dat2$lod), lod2)
  
  # May have different censoring percentages
  expect_true(is.numeric(attr(dat1, "pct_censored")))
  expect_true(is.numeric(attr(dat2, "pct_censored")))
})

test_that("system manages invalid input specifications appropriately", {
  data("immunoplex_example", package = "immunoPlex")
  
  # Non-existent cytokine
  expect_error(
    prepare_cytokine_data("NonExistentCytokine"),
    "not found"
  )
  
  # Cytokine without LOD information
  all_cytokines <- names(immunoplex_example$expression)
  lod_cytokines <- immunoplex_example$lod_lookup$cytokine[!is.na(immunoplex_example$lod_lookup$lod)]
  no_lod_cytokines <- setdiff(all_cytokines, lod_cytokines)
  
  if (length(no_lod_cytokines) > 0) {
    # Test should expect error OR success (since some cytokines might have NA LOD vs missing)
    result <- tryCatch({
      prepare_cytokine_data(no_lod_cytokines[1])
      "success"
    }, error = function(e) {
      e$message
    })
    
    # Either it works (cytokine found with NA LOD) or gives appropriate error
    expect_true(result == "success" || grepl("No LOD information found", result))
  } else {
    # All cytokines have LOD info in test data - skip this test
    skip("All cytokines in test data have LOD information")
  }
})

test_that("custom data source works", {
  data("immunoplex_example", package = "immunoPlex")
  
  # Use custom data source
  custom_data <- immunoplex_example
  
  available_cytokines <- custom_data$lod_lookup$cytokine[!is.na(custom_data$lod_lookup$lod)]
  if (length(available_cytokines) == 0) {
    skip("No cytokines with LOD information found")
  }
  
  dat <- prepare_cytokine_data(available_cytokines[1], data_source = custom_data)
  
  expect_true(is.data.frame(dat))
  expect_equal(attr(dat, "cytokine"), available_cytokines[1])
})

test_that("system maintains data integrity across transformations", {
  data("immunoplex_example", package = "immunoPlex")
  
  lod_lookup <- immunoplex_example$lod_lookup
  available_cytokines <- lod_lookup$cytokine[!is.na(lod_lookup$lod)]
  
  if (length(available_cytokines) == 0) {
    skip("No cytokines with LOD information found")
  }
  
  target_cyto <- available_cytokines[1]
  dat <- prepare_cytokine_data(target_cyto)
  
  # Check that sample IDs match between expression and metadata
  expect_true(all(dat$sample_id %in% rownames(immunoplex_example$expression)))
  expect_true(all(dat$sample_id %in% immunoplex_example$metadata$sample_id))
  
  # Check that values match original expression data
  original_values <- immunoplex_example$expression[dat$sample_id, target_cyto]
  expect_equal(dat$value, as.numeric(original_values))
  
  # Check that LOD matches lookup table
  expected_lod <- lod_lookup$lod[lod_lookup$cytokine == target_cyto]
  expect_equal(unique(dat$lod), expected_lod)
})

test_that("system processes incomplete data matrices appropriately", {
  data("immunoplex_example", package = "immunoPlex")
  
  lod_lookup <- immunoplex_example$lod_lookup
  available_cytokines <- lod_lookup$cytokine[!is.na(lod_lookup$lod)]
  
  if (length(available_cytokines) == 0) {
    skip("No cytokines with LOD information found")
  }
  
  dat <- prepare_cytokine_data(available_cytokines[1])
  
  # Should handle NA values in metadata gracefully
  na_timepoint <- is.na(dat$timepoint)
  na_disease <- is.na(dat$disease)
  na_age <- is.na(dat$age)
  
  # Test that the function completes without error even with NAs
  expect_true(is.data.frame(dat))
  expect_equal(nrow(dat), nrow(immunoplex_example$expression))
})

})

describe("Analyte enumeration and characterization", {

test_that("core analyte listing protocol executes correctly", {
  data("immunoplex_example", package = "immunoPlex")
  
  cytokine_info <- list_cytokines()
  
  # Check structure
  expect_true(is.data.frame(cytokine_info))
  expect_true(all(c("cytokine", "lod", "n_obs", "n_censored", "pct_censored") %in% names(cytokine_info)))
  
  # Check data types
  expect_type(cytokine_info$cytokine, "character")
  expect_type(cytokine_info$lod, "double")
  expect_type(cytokine_info$n_obs, "integer")
  expect_type(cytokine_info$n_censored, "integer")
  expect_type(cytokine_info$pct_censored, "double")
  
  # Should include all cytokines from expression data
  all_cytokines <- names(immunoplex_example$expression)
  expect_true(all(all_cytokines %in% cytokine_info$cytokine))
})

test_that("censoring metric computation protocol executes correctly", {
  data("immunoplex_example", package = "immunoPlex")
  
  cytokine_info <- list_cytokines()
  
  # Pick a cytokine with LOD information for manual verification
  lod_lookup <- immunoplex_example$lod_lookup
  test_cytokines <- lod_lookup$cytokine[!is.na(lod_lookup$lod)]
  
  if (length(test_cytokines) > 0) {
    test_cyto <- test_cytokines[1]
    info_row <- cytokine_info[cytokine_info$cytokine == test_cyto, ]
    
    # Manual calculation
    expr_values <- immunoplex_example$expression[, test_cyto]
    lod_value <- lod_lookup$lod[lod_lookup$cytokine == test_cyto]
    
    expected_n_obs <- sum(!is.na(expr_values))
    expected_n_censored <- sum(expr_values < lod_value, na.rm = TRUE)
    expected_pct_censored <- round(100 * expected_n_censored / expected_n_obs, 1)
    
    expect_equal(info_row$n_obs, expected_n_obs)
    expect_equal(info_row$n_censored, expected_n_censored)
    expect_equal(info_row$pct_censored, expected_pct_censored)
  }
})

test_that("censoring-based ordering protocol executes correctly", {
  data("immunoplex_example", package = "immunoPlex")
  
  cytokine_info <- list_cytokines()
  
  # Should be sorted by pct_censored in descending order
  pct_values <- cytokine_info$pct_censored[!is.na(cytokine_info$pct_censored)]
  if (length(pct_values) > 1) {
    expect_true(all(diff(pct_values) <= 0))
  }
})

test_that("system processes analytes without detection limits appropriately", {
  data("immunoplex_example", package = "immunoPlex")
  
  cytokine_info <- list_cytokines()
  
  # Cytokines without LOD should have NA for LOD and n_censored
  no_lod_rows <- cytokine_info[is.na(cytokine_info$lod), ]
  
  if (nrow(no_lod_rows) > 0) {
    expect_true(all(is.na(no_lod_rows$n_censored)))
    expect_true(all(is.na(no_lod_rows$pct_censored)))
  }
})

test_that("custom data source works", {
  data("immunoplex_example", package = "immunoPlex")
  
  # Use custom data source
  custom_data <- immunoplex_example
  cytokine_info <- list_cytokines(data_source = custom_data)
  
  expect_true(is.data.frame(cytokine_info))
  expect_true(all(names(custom_data$expression) %in% cytokine_info$cytokine))
})

test_that("system provides comprehensive analyte characteristics for estimation", {
  data("immunoplex_example", package = "immunoPlex")
  
  cytokine_info <- list_cytokines()
  
  # Should help identify good candidates for censored modeling
  high_censoring <- cytokine_info[cytokine_info$pct_censored > 10 & !is.na(cytokine_info$pct_censored), ]
  moderate_censoring <- cytokine_info[cytokine_info$pct_censored > 5 & cytokine_info$pct_censored <= 10 & 
                                     !is.na(cytokine_info$pct_censored), ]
  
  # Should have some cytokines in different categories
  total_with_lod <- sum(!is.na(cytokine_info$lod))
  expect_true(total_with_lod > 0)  # At least some should have LOD info
  
  # Information should be sufficient for choosing test cases
  expect_true(all(cytokine_info$n_obs >= 0))
  expect_true(all(cytokine_info$n_censored[!is.na(cytokine_info$n_censored)] >= 0))
  expect_true(all(cytokine_info$pct_censored[!is.na(cytokine_info$pct_censored)] >= 0))
  expect_true(all(cytokine_info$pct_censored[!is.na(cytokine_info$pct_censored)] <= 100))
})

})

describe("Statistical modeling integration validation", {

test_that("processed data structure supports single distribution estimation", {
  skip_if_not_installed("glmmTMB")
  
  data("immunoplex_example", package = "immunoPlex")
  
  lod_lookup <- immunoplex_example$lod_lookup
  available_cytokines <- lod_lookup$cytokine[!is.na(lod_lookup$lod)]
  
  if (length(available_cytokines) == 0) {
    skip("No cytokines with LOD information found")
  }
  
  dat <- prepare_cytokine_data(available_cytokines[1])
  
  # Should work with fit_one without errors
  fit <- fit_one(dat, family = "gamma", random = "")
  
  expect_s3_class(fit, "immuno_fit")
  expect_equal(fit$family, "gamma")
})

test_that("processed data structure supports multi-distribution comparison", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  
  data("immunoplex_example", package = "immunoPlex")
  
  lod_lookup <- immunoplex_example$lod_lookup
  available_cytokines <- lod_lookup$cytokine[!is.na(lod_lookup$lod)]
  
  if (length(available_cytokines) == 0) {
    skip("No cytokines with LOD information found")
  }
  
  dat <- prepare_cytokine_data(available_cytokines[1])
  
  # Should work with fit_models
  models <- fit_models(dat, families = c("gamma", "tobit"), random = "")
  
  expect_s3_class(models, "immuno_model_set")
  expect_true(length(models$models) >= 1)
})

test_that("processed data structure supports detection limit methodology comparison", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  
  data("immunoplex_example", package = "immunoPlex")
  
  lod_lookup <- immunoplex_example$lod_lookup
  available_cytokines <- lod_lookup$cytokine[!is.na(lod_lookup$lod)]
  
  if (length(available_cytokines) == 0) {
    skip("No cytokines with LOD information found")
  }
  
  dat <- prepare_cytokine_data(available_cytokines[1])
  
  # Skip if insufficient censoring
  if (attr(dat, "n_censored") < 3) {
    skip("Insufficient censoring for LOD comparison")
  }
  
  # Should work with compare_lod_models
  comparison <- compare_lod_models(
    dat,
    families = c("gamma", "tobit"),
    lod_methods = c("half"),
    include_log_transform = FALSE,
    random = ""
  )
  
  expect_s3_class(comparison, "immuno_lod_comparison")
  expect_true(length(comparison$models) >= 1)
})

test_that("documented analysis protocol executes correctly", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  skip_if_not_installed("ggplot2")
  
  data("immunoplex_example", package = "immunoPlex")
  
  # List cytokines to find good test case
  cyto_info <- list_cytokines()
  good_cytokines <- cyto_info[cyto_info$pct_censored >= 5 & cyto_info$pct_censored <= 50 & 
                             !is.na(cyto_info$pct_censored), ]
  
  if (nrow(good_cytokines) == 0) {
    skip("No cytokines with suitable censoring levels found")
  }
  
  # Use the first suitable cytokine
  target_cyto <- good_cytokines$cytokine[1]
  
  # Prepare data
  dat <- prepare_cytokine_data(target_cyto)
  
  # Fit models
  models <- fit_models(dat, families = c("gamma", "tobit"), random = "")
  
  # Plot comparison
  p <- plot(models, plot_best = FALSE)
  
  # Should all work without errors
  expect_s3_class(models, "immuno_model_set")
  expect_s3_class(p, "gg")
  expect_true(length(models$models) >= 1)
})

})