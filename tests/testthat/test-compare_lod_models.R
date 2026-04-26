# Limit of Detection (LOD) Model Comparison Validation Protocol
#
# Validates statistical model performance under censoring conditions:
# - Data preprocessing protocols:
#   * LOD fraction methods (half, zero)
#   * Minimum-based methods (halfmin)
#   * Direct substitution (lod)
# - Censoring-aware estimation:
#   * Tobit regression (left-censoring)
#   * Accelerated Failure Time models
# - Distribution selection criteria
# - Model evaluation metrics
# - Object-oriented implementation validation
#
# Protocol sections:
# 1. Core estimation functionality
# 2. Preprocessing method validation
# 3. Object-oriented interface verification
#
# Statistical dependencies:
# - glmmTMB: Gamma distribution estimation
# - survival: Accelerated Failure Time models
# - censReg: Tobit regression implementation
# - ggplot2: Statistical visualization

library(testthat)

# Initialize validation dataset with defined censoring characteristics
setup_lod_test_data <- function() {
  data("immunoplex_example", package = "immunoPlex")
  
  raw_expr <- immunoplex_example$expression
  meta <- immunoplex_example$metadata
  lod_lookup <- immunoplex_example$lod_lookup
  
  # Identify analytes with defined detection limits
  lod_cytokines <- lod_lookup$cytokine[!is.na(lod_lookup$lod)]
  
  # Select analyte with optimal censoring proportion
  best_cyto <- NULL
  best_pct <- 0
  for (cyto in lod_cytokines) {
    lod_val <- lod_lookup$lod[lod_lookup$cytokine == cyto]
    pct_cens <- mean(raw_expr[, cyto] < lod_val, na.rm = TRUE)
    if (pct_cens > best_pct && pct_cens < 0.9) {  # Target: moderate censoring rate
      best_cyto <- cyto
      best_pct <- pct_cens
    }
  }
  
  if (is.null(best_cyto)) {
    best_cyto <- lod_cytokines[1]  # Fallback
  }
  
  target_lod <- lod_lookup$lod[lod_lookup$cytokine == best_cyto]
  
  # Create test data
  dat <- data.frame(
    sample_id = rownames(raw_expr),
    subject_id = meta$subject_id[match(rownames(raw_expr), meta$sample_id)],
    value = raw_expr[, best_cyto],
    timepoint = meta$timepoint[match(rownames(raw_expr), meta$sample_id)],
    disease = meta$disease[match(rownames(raw_expr), meta$sample_id)],
    age = meta$age[match(rownames(raw_expr), meta$sample_id)],
    lod = target_lod,
    ulod = lod_lookup$ulod[lod_lookup$cytokine == best_cyto],
    stringsAsFactors = FALSE
  )
  
  dat$cens_lod <- dat$value < dat$lod
  dat$cens_ulod <- !is.na(dat$ulod) & dat$value > dat$ulod
  
  return(dat)
}

describe("Model comparison under censoring conditions", {

test_that("core estimation functionality executes correctly", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  
  dat <- setup_lod_test_data()
  
  # Validate minimum censoring requirements
  if (sum(dat$cens_lod) < 5) {
    skip("Insufficient censored observations for statistical validation")
  }
  
  # Execute model comparison protocol
  comparison <- compare_lod_models(
    dat,
    families = c("gamma", "tobit"),      # Distribution specifications
    lod_methods = c("half", "halfmin"),  # Preprocessing methods
    include_log_transform = FALSE,       # Direct scale estimation
    random = ""                          # Fixed effects only
  )
  
  expect_s3_class(comparison, "immuno_lod_comparison")
  expect_true("models" %in% names(comparison))
  expect_true("comparison" %in% names(comparison))
  expect_true("best_model" %in% names(comparison))
  expect_true("best_name" %in% names(comparison))
  
  # Check comparison table structure
  expect_true(is.data.frame(comparison$comparison))
  expect_true(all(c("model", "family", "lod_method", "aic", "delta_aic") %in% names(comparison$comparison)))
  
  # Should have both preprocessing and censoring approaches
  expect_true(any(grepl("gamma_", comparison$comparison$model)))
  expect_true(any(comparison$comparison$lod_method == "native_censoring"))
})

test_that("log transformation option works", {
  skip_if_not_installed("glmmTMB")
  
  dat <- setup_lod_test_data()
  
  # Skip if insufficient censoring
  if (sum(dat$cens_lod) < 3) {
    skip("Insufficient censoring in test data")
  }
  
  comparison <- compare_lod_models(
    dat,
    families = c("gamma"),
    lod_methods = c("half"),
    include_log_transform = TRUE,
    random = ""
  )
  
  # Should have both regular and log-transformed gamma models
  model_names <- names(comparison$models)
  expect_true(any(grepl("gamma_half$", model_names)))
  expect_true(any(grepl("gamma_log_half", model_names)))
})

test_that("multiple LOD methods are tested", {
  skip_if_not_installed("glmmTMB")
  
  dat <- setup_lod_test_data()
  
  # Skip if insufficient censoring
  if (sum(dat$cens_lod) < 3) {
    skip("Insufficient censoring in test data")
  }
  
  lod_methods <- c("half", "halfmin", "zero")
  comparison <- compare_lod_models(
    dat,
    families = c("gamma"),
    lod_methods = lod_methods,
    include_log_transform = FALSE,
    random = ""
  )
  
  # Should have models for each LOD method
  tested_methods <- unique(comparison$comparison$lod_method[grepl("gamma_", comparison$comparison$model)])
  expect_true(length(tested_methods) >= 2)  # At least some should work
  expect_true(all(tested_methods %in% lod_methods))
})

test_that("censoring-aware models are included", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  
  dat <- setup_lod_test_data()
  
  comparison <- compare_lod_models(
    dat,
    families = c("gamma", "tobit", "aft"),
    lod_methods = c("half"),
    include_log_transform = FALSE,
    random = ""
  )
  
  # Should have native censoring models
  native_models <- comparison$comparison$model[comparison$comparison$lod_method == "native_censoring"]
  expect_true(length(native_models) >= 1)
  expect_true(any(c("tobit", "aft") %in% native_models))
})

test_that("censReg integration works in comparison", {
  skip_if_not_installed("censReg")
  
  dat <- setup_lod_test_data()
  
  # Skip if insufficient censoring
  if (sum(dat$cens_lod) < 5) {
    skip("Insufficient censoring in test data")
  }
  
  comparison <- compare_lod_models(
    dat,
    families = c("gamma", "tobit_censreg"),
    lod_methods = c("half"),
    include_log_transform = FALSE,
    random = ""
  )
  
  # Check if tobit_censreg was fitted successfully
  if ("tobit_censreg" %in% names(comparison$models)) {
    censreg_model <- comparison$models$tobit_censreg
    expect_equal(censreg_model$family, "tobit_censreg")
    expect_s3_class(censreg_model$model, "censReg")
  }
})

test_that("input validation works", {
  dat <- setup_lod_test_data()
  
  # Missing required columns
  dat_bad <- dat[, !names(dat) %in% "cens_lod"]
  expect_error(
    compare_lod_models(dat_bad),
    "cens_lod"
  )
  
  # Invalid LOD methods
  expect_error(
    compare_lod_models(dat, lod_methods = "invalid_method"),
    "Unsupported LOD method"
  )
})

test_that("comparison table is sorted by AIC", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  
  dat <- setup_lod_test_data()
  
  # Skip if insufficient censoring
  if (sum(dat$cens_lod) < 3) {
    skip("Insufficient censoring in test data")
  }
  
  comparison <- compare_lod_models(
    dat,
    families = c("gamma", "tobit"),
    lod_methods = c("half", "halfmin"),
    include_log_transform = FALSE,
    random = ""
  )
  
  # Comparison table should be sorted by AIC
  aic_values <- comparison$comparison$aic
  expect_true(all(diff(aic_values) >= 0))
  
  # Best model should have lowest AIC
  expect_equal(comparison$best_model$aic, min(aic_values))
})

test_that("attributes are preserved correctly", {
  skip_if_not_installed("glmmTMB")
  
  dat <- setup_lod_test_data()
  
  lod_methods <- c("half", "halfmin")
  families <- c("gamma", "tobit")
  
  comparison <- compare_lod_models(
    dat,
    families = families,
    lod_methods = lod_methods,
    include_log_transform = TRUE,
    random = ""
  )
  
  expect_equal(comparison$lod_methods_tested, lod_methods)
  expect_equal(comparison$families_tested, families)
  expect_true(comparison$include_log_transform)
})

test_that("handles model failures gracefully", {
  skip_if_not_installed("glmmTMB")
  
  dat <- setup_lod_test_data()
  
  # Create problematic data
  dat_problematic <- dat
  dat_problematic$value[1:80] <- NA  # Add many missing values
  
  # Should still work for some models (no specific message expected)
  comparison <- compare_lod_models(
    dat_problematic,
    families = c("gamma", "tobit"),
    lod_methods = c("half"),
    include_log_transform = FALSE,
    random = ""
  )
  
  # At least one model should succeed
  expect_true(length(comparison$models) >= 1)
  expect_s3_class(comparison, "immuno_lod_comparison")
})

})

describe("Preprocessing method implementation", {

test_that("LOD fraction method implementation is correct", {
  dat <- data.frame(
    value = c(1, 2, 0.5, 3, 0.3),
    lod = rep(1, 5),
    cens_lod = c(FALSE, FALSE, TRUE, FALSE, TRUE)
  )
  
  result <- immunoPlex:::.apply_lod_method(dat, "half")
  
  # Censored values should be replaced with LOD/2 = 0.5
  expect_equal(result$value[dat$cens_lod], c(0.5, 0.5))
  # Non-censored values should be unchanged
  expect_equal(result$value[!dat$cens_lod], c(1, 2, 3))
})

test_that("LOD method 'zero' works correctly", {
  dat <- data.frame(
    value = c(1, 2, 0.5, 3, 0.3),
    lod = rep(1, 5),
    cens_lod = c(FALSE, FALSE, TRUE, FALSE, TRUE)
  )
  
  result <- immunoPlex:::.apply_lod_method(dat, "zero")
  
  # Censored values should be replaced with 0
  expect_equal(result$value[dat$cens_lod], c(0, 0))
})

test_that("LOD method 'lod' works correctly", {
  dat <- data.frame(
    value = c(1, 2, 0.5, 3, 0.3),
    lod = rep(1, 5),
    cens_lod = c(FALSE, FALSE, TRUE, FALSE, TRUE)
  )
  
  result <- immunoPlex:::.apply_lod_method(dat, "lod")
  
  # Censored values should be replaced with LOD value
  expect_equal(result$value[dat$cens_lod], c(1, 1))
})

test_that("LOD method 'halfmin' works correctly", {
  dat <- data.frame(
    value = c(4, 2, 0.5, 3, 0.3),  # min positive uncensored = 2
    lod = rep(1, 5),
    cens_lod = c(FALSE, FALSE, TRUE, FALSE, TRUE)
  )
  
  result <- immunoPlex:::.apply_lod_method(dat, "halfmin")
  
  # Censored values should be replaced with min_positive/2 = 1
  expect_equal(result$value[dat$cens_lod], c(1, 1))
})

test_that("invalid LOD method throws error", {
  dat <- data.frame(
    value = c(1, 2, 0.5),
    lod = rep(1, 3),
    cens_lod = c(FALSE, FALSE, TRUE)
  )
  
  expect_error(
    immunoPlex:::.apply_lod_method(dat, "invalid"),
    "Unsupported LOD method"
  )
})

})

describe("Object-oriented interface implementation", {

test_that("print method generates appropriate output", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  
  dat <- setup_lod_test_data()
  
  # Skip if insufficient censoring
  if (sum(dat$cens_lod) < 3) {
    skip("Insufficient censoring in test data")
  }
  
  comparison <- compare_lod_models(
    dat,
    families = c("gamma", "tobit"),
    lod_methods = c("half"),
    include_log_transform = FALSE,
    random = ""
  )
  
  # Should print without error
  expect_output(print(comparison), "immuno_lod_comparison")
  expect_output(print(comparison), "LOD methods tested")
  expect_output(print(comparison), "Best approach")
})

test_that("summary method provides comprehensive model assessment", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  
  dat <- setup_lod_test_data()
  
  # Skip if insufficient censoring
  if (sum(dat$cens_lod) < 3) {
    skip("Insufficient censoring in test data")
  }
  
  comparison <- compare_lod_models(
    dat,
    families = c("gamma", "tobit"),
    lod_methods = c("half"),
    include_log_transform = FALSE,
    random = ""
  )
  
  # Should summarize without error
  expect_output(summary(comparison), "immuno_lod_comparison Summary")
  expect_output(summary(comparison), "Total models fitted")
  expect_output(summary(comparison), "Best Model Details")
})

test_that("visualization method produces valid statistical graphics", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  skip_if_not_installed("ggplot2")
  
  dat <- setup_lod_test_data()
  
  # Skip if insufficient censoring
  if (sum(dat$cens_lod) < 3) {
    skip("Insufficient censoring in test data")
  }
  
  comparison <- compare_lod_models(
    dat,
    families = c("gamma", "tobit"),
    lod_methods = c("half"),
    include_log_transform = FALSE,
    random = ""
  )
  
  # Should plot without error and return ggplot object
  p <- plot(comparison, plot_best = FALSE)
  expect_s3_class(p, "gg")
  
  # Test with plot_best = TRUE
  expect_output(
    plot(comparison, plot_best = TRUE, save_pdf = FALSE),
    "Diagnostic plots for best approach"
  )
})

})