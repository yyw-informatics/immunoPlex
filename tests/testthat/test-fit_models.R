# Statistical Model Estimation and Selection Protocol
#
# Validates model estimation and selection procedures:
# 1. Distribution Family Assessment:
#    - Gamma regression
#    - Tobit regression
#    - Accelerated Failure Time models
#
# 2. Model Selection Criteria:
#    - Information criterion evaluation
#    - Likelihood-based comparison
#    - Cross-validation metrics
#
# 3. Estimation Robustness:
#    - Convergence diagnostics
#    - Numerical stability
#    - Parameter constraint validation
#
# 4. Implementation Validation:
#    - Object-oriented interface
#    - Diagnostic visualization
#    - External package integration

library(testthat)

# Initialize validation dataset with defined detection limits
setup_test_data <- function() {
  data("immunoplex_example", package = "immunoPlex")
  
  raw_expr <- immunoplex_example$expression
  meta <- immunoplex_example$metadata
  lod_lookup <- immunoplex_example$lod_lookup
  
  # Pin the analyte deterministically so row-reordering in
  # immunoplex_example$lod_lookup cannot silently change which cytokine drives
  # the fixture. Selected by sorted name rather than by a literal, so renaming
  # the example analytes cannot break the fixture either — a hard-coded
  # "Cytokine_01" is what broke it when the dataset moved to real cytokine names.
  complete_lod <- sort(lod_lookup$cytokine[!is.na(lod_lookup$lod)])
  stopifnot(length(complete_lod) > 0L)
  target_cyto <- complete_lod[1]
  stopifnot(target_cyto %in% colnames(raw_expr))
  target_lod <- lod_lookup$lod[lod_lookup$cytokine == target_cyto]
  
  # Create test data
  dat <- data.frame(
    sample_id = rownames(raw_expr),
    subject_id = meta$subject_id[match(rownames(raw_expr), meta$sample_id)],
    value = raw_expr[, target_cyto],
    timepoint = meta$timepoint[match(rownames(raw_expr), meta$sample_id)],
    disease = meta$disease[match(rownames(raw_expr), meta$sample_id)],
    age = meta$age[match(rownames(raw_expr), meta$sample_id)],
    lod = target_lod,
    ulod = lod_lookup$ulod[lod_lookup$cytokine == target_cyto],
    stringsAsFactors = FALSE
  )
  
  dat$cens_lod <- dat$value < dat$lod
  dat$cens_ulod <- !is.na(dat$ulod) & dat$value > dat$ulod
  
  return(dat)
}

describe("Model estimation and selection functionality", {

test_that("core model comparison protocol executes correctly", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  
  dat <- setup_test_data()
  
  # Test with basic families
  models <- fit_models(dat, families = c("gamma", "tobit"), random = "")
  
  expect_s3_class(models, "immuno_model_set")
  expect_true("models" %in% names(models))
  expect_true("comparison" %in% names(models))
  expect_true("best_model" %in% names(models))
  expect_true("best_family" %in% names(models))
  
  # Check comparison table structure (including new estimand column)
  expect_true(is.data.frame(models$comparison))
  expect_true(all(c("family", "estimand", "aic", "bic", "loglik", "delta_aic") %in% names(models$comparison)))
  expect_equal(nrow(models$comparison), length(models$models))
  
  # Test estimand values are correct
  gamma_row <- models$comparison$family == "gamma"
  tobit_row <- models$comparison$family == "tobit"
  if (any(gamma_row)) {
    expect_equal(models$comparison$estimand[gamma_row], "ratio_of_means")
  }
  if (any(tobit_row)) {
    expect_equal(models$comparison$estimand[tobit_row], "ratio_of_means")
  }
  
  # Best model should be the one with lowest AIC
  expect_equal(models$best_family, models$comparison$family[1])
  expect_equal(min(models$comparison$aic), models$best_model$aic)
})

test_that("optimal model selection parameter functions correctly", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  
  dat <- setup_test_data()
  
  # Test return_best = TRUE
  best_model <- fit_models(dat, families = c("gamma", "tobit"), 
                          random = "", return_best = TRUE)
  
  expect_s3_class(best_model, "immuno_fit")
  expect_true(best_model$converged)
  
  # Test return_best = FALSE (default)
  all_models <- fit_models(dat, families = c("gamma", "tobit"), 
                          random = "", return_best = FALSE)
  
  expect_s3_class(all_models, "immuno_model_set")
  expect_equal(best_model$family, all_models$best_family)
})

test_that("option_state snapshot captures fit_one thresholds at call time", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")

  dat <- setup_test_data()

  # Snapshot current options; restore at exit.
  orig_subjects <- getOption("fit_one.min_subjects")
  orig_reps     <- getOption("fit_one.min_reps")
  on.exit({
    options(fit_one.min_subjects = orig_subjects,
            fit_one.min_reps     = orig_reps)
  }, add = TRUE)

  set_immunoplex_config(random_effects_min_subjects = 25L,
                        random_effects_min_reps     = 4L)

  models <- fit_models(dat, families = c("gamma", "tobit"),
                      random = "", quiet = TRUE)

  expect_true("option_state" %in% names(models))
  expect_equal(models$option_state$min_subjects, 25L)
  expect_equal(models$option_state$min_reps, 4L)
})

test_that("return_best = TRUE attaches comparison table as attribute", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")

  dat <- setup_test_data()

  best_model <- fit_models(dat, families = c("gamma", "tobit"),
                           random = "", return_best = TRUE)

  comparison <- attr(best_model, "comparison")
  expect_true(is.data.frame(comparison))
  expect_true(all(c("family", "estimand", "aic", "bic", "loglik",
                    "n_cens_lod", "n_cens_ulod", "delta_aic")
                  %in% names(comparison)))
  expect_equal(comparison$family[1], best_model$family)
  expect_equal(comparison$delta_aic[1], 0)
})

test_that("system manages estimation non-convergence appropriately", {
  skip_if_not_installed("glmmTMB")
  
  dat <- setup_test_data()
  
  # Create problematic data that might cause some models to fail
  dat_problematic <- dat
  dat_problematic$value[1:50] <- NA  # Add many missing values
  
  # Should still work for at least some models (no specific message expected)
  models <- fit_models(dat_problematic, families = c("gamma", "tobit", "aft"), 
                      random = "")
  
  # At least one model should succeed
  expect_true(length(models$models) >= 1)
  expect_s3_class(models, "immuno_model_set")
})

test_that("censored regression integration protocol executes correctly", {
  skip_if_not_installed("censReg")
  
  dat <- setup_test_data()
  
  # Skip if no censored values
  if (sum(dat$cens_lod) == 0) {
    skip("No left-censored values found in test data")
  }
  
  models <- fit_models(dat, families = c("gamma", "tobit_censreg"), random = "")
  
  expect_s3_class(models, "immuno_model_set")
  
  # Check if tobit_censreg model was fitted
  if ("tobit_censreg" %in% names(models$models)) {
    censreg_model <- models$models$tobit_censreg
    expect_equal(censreg_model$family, "tobit_censreg")
    expect_s3_class(censreg_model$model, "censReg")
  }
})

test_that("system validates input parameter specifications", {
  dat <- setup_test_data()
  
  # Invalid families
  expect_error(
    fit_models(dat, families = "invalid_family"),
    "should be one of"
  )
  
  # Test error handling: Empty data frame should fail gracefully
  # This test intentionally passes invalid data to verify error handling
  cat("Testing error handling with empty data frame (expected to fail)...\n")
  expect_error(
    fit_models(data.frame(), quiet = TRUE),
    "No models converged successfully"
  )
})

test_that("model ranking protocol implements information criterion ordering", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  
  dat <- setup_test_data()
  
  models <- fit_models(dat, families = c("gamma", "tobit", "aft"), random = "")
  
  # Comparison table should be sorted by AIC (ascending)
  aic_values <- models$comparison$aic
  expect_true(all(diff(aic_values) >= 0))
  
  # Delta AIC should start from 0
  expect_equal(models$comparison$delta_aic[1], 0)
  expect_true(all(models$comparison$delta_aic >= 0))
})

test_that("system validates estimation completion count", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  
  dat <- setup_test_data()
  
  families <- c("gamma", "tobit")
  models <- fit_models(dat, families = families, random = "")
  
  expect_equal(models$n_models, length(models$models))
  expect_true(models$n_models <= length(families))  # Some models might fail
  expect_true(models$n_models >= 1)  # At least one should succeed
})

test_that("system processes all distribution family combinations correctly", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  
  dat <- setup_test_data()
  
  # Test different family combinations
  family_combinations <- list(
    c("gamma"),
    c("tobit"),
    c("aft"),
    c("gamma", "tobit"),
    c("gamma", "aft"),
    c("tobit", "aft"),
    c("gamma", "tobit", "aft")
  )
  
  for (families in family_combinations) {
    models <- fit_models(dat, families = families, random = "")
    expect_s3_class(models, "immuno_model_set")
    expect_true(length(models$models) >= 1)
    expect_true(all(models$comparison$family %in% families))
  }
})

test_that("ulod parameter works correctly in fit_models", {
  skip_if_not_installed("survival")
  
  dat <- setup_test_data()
  
  # Add some right-censored values for testing
  dat$cens_ulod[1:2] <- TRUE
  dat$ulod[1:2] <- max(dat$value, na.rm = TRUE) * 1.5
  
  # Test with ulod = TRUE
  models_with_ulod <- fit_models(dat, families = c("tobit", "aft"), 
                                ulod = TRUE, random = "")
  expect_s3_class(models_with_ulod, "immuno_model_set")
  
  # Test with ulod = FALSE (should ignore right censoring)
  models_without_ulod <- fit_models(dat, families = c("tobit", "aft"), 
                                   ulod = FALSE, random = "")
  expect_s3_class(models_without_ulod, "immuno_model_set")
})

test_that("estimand tracking works across all families", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  
  dat <- setup_test_data()
  
  models <- fit_models(dat, families = c("gamma", "tobit", "aft"), random = "")
  
  # Check that each family has the correct estimand
  for (i in seq_len(nrow(models$comparison))) {
    family <- models$comparison$family[i]
    estimand <- models$comparison$estimand[i]
    
    if (family %in% c("gamma", "tobit", "tobit_censreg")) {
      expect_equal(estimand, "ratio_of_means")
    } else if (family == "aft") {
      expect_equal(estimand, "ratio_of_medians")
    }
  }
})

test_that("mixed estimand comparison emits a warning (unless quiet)", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")

  dat <- setup_test_data()

  # Collect only warnings that mention estimand heterogeneity, so unrelated
  # warnings emitted by fit_one (e.g. gamma-with-censoring advisory) don't
  # interfere with assertions.
  collect_estimand_warnings <- function(expr) {
    msgs <- character()
    withCallingHandlers(
      expr,
      warning = function(w) {
        if (grepl("mixes estimands", conditionMessage(w))) {
          msgs <<- c(msgs, conditionMessage(w))
        }
        invokeRestart("muffleWarning")
      }
    )
    msgs
  }

  # gamma (ratio_of_means) + aft (ratio_of_medians) → mixed, should warn
  w_mixed <- collect_estimand_warnings(
    fit_models(dat, families = c("gamma", "aft"), random = "")
  )
  expect_true(length(w_mixed) >= 1)
  expect_match(w_mixed[1], "ratio_of_means")
  expect_match(w_mixed[1], "ratio_of_medians")

  # quiet = TRUE suppresses the estimand warning
  w_quiet <- collect_estimand_warnings(
    fit_models(dat, families = c("gamma", "aft"), random = "", quiet = TRUE)
  )
  expect_length(w_quiet, 0)

  # gamma + tobit are both ratio_of_means → no estimand warning
  w_same <- collect_estimand_warnings(
    fit_models(dat, families = c("gamma", "tobit"), random = "")
  )
  expect_length(w_same, 0)
})

test_that("fit_models accepts gaussian family alone and in combination", {
  skip_if_not_installed("glmmTMB")

  dat <- setup_test_data()
  # Gaussian expects log-scale response (see fit_one gaussian branch); apply
  # log1p so values stay strictly positive and gamma remains fittable in
  # the combined candidate set.
  dat$value <- log1p(dat$value)
  dat$lod   <- log1p(dat$lod)
  dat$ulod  <- log1p(dat$ulod)
  dat$cens_lod  <- dat$value < dat$lod
  dat$cens_ulod <- !is.na(dat$ulod) & dat$value > dat$ulod

  # Gaussian alone
  m1 <- fit_models(dat, families = "gaussian", random = "", quiet = TRUE)
  expect_s3_class(m1, "immuno_model_set")
  expect_true("gaussian" %in% names(m1$models))
  expect_equal(m1$best_family, "gaussian")
  gauss_row <- m1$comparison[m1$comparison$family == "gaussian", ]
  expect_equal(gauss_row$estimand, "mean_difference")

  # Gaussian + gamma combination
  m2 <- fit_models(dat, families = c("gamma", "gaussian"),
                   random = "", quiet = TRUE)
  expect_s3_class(m2, "immuno_model_set")
  expect_true("gaussian" %in% m2$comparison$family)
})

test_that("enhanced error handling works correctly", {
  dat <- setup_test_data()

  # Test with row-specific LODs for censReg (should produce informative error)
  if (requireNamespace("censReg", quietly = TRUE)) {
    # Create data with varying LODs
    dat$lod[1:5] <- dat$lod[1:5] * 0.5  # Different LOD values
    dat$cens_lod[1:5] <- TRUE

    # Should get informative error message for censReg
    models <- fit_models(dat, families = c("gamma", "tobit_censreg"),
                        quiet = TRUE, random = "")

    # censReg should fail, but gamma should succeed
    expect_true("gamma" %in% names(models$models))
    expect_false("tobit_censreg" %in% names(models$models))
  }
})

})

describe("Classed-condition dispatch between fit_one and fit_models", {

  # Minimal valid frame — enough rows & columns that fit_one reaches its
  # family-specific branches without tripping the generic column/dimension
  # guards. Each per-branch smoke test mutates one field to trigger the
  # corresponding class.
  make_minimal_dat <- function() {
    data.frame(
      value      = log1p(seq(1, 15)),
      cens_lod   = FALSE,
      cens_ulod  = FALSE,
      lod        = log1p(0.5),
      ulod       = log1p(100),
      timepoint  = rep("T1", 15),
      disease    = rep("Control", 15),
      age        = rep(30, 15),
      subject_id = seq_len(15),
      stringsAsFactors = FALSE
    )
  }

  test_that("fit_one_missing_cols class actually fires from fit_one", {
    # Gap 1 anchor test: prove the classed condition is reachable via tryCatch
    # on the tag (not just via a grepl() of the message text).
    got <- tryCatch(
      fit_one(data.frame(value = 1:5)),
      fit_one_missing_cols = function(e) TRUE,
      error = function(e) FALSE
    )
    expect_true(got)
  })

  test_that("fit_one_not_dataframe fires on non-data.frame input", {
    got <- tryCatch(
      fit_one(list(value = 1:5)),
      fit_one_not_dataframe = function(e) TRUE,
      error = function(e) FALSE
    )
    expect_true(got)
  })

  test_that("fit_one_formula_vars fires on bogus placeholder variable", {
    dat <- make_minimal_dat()
    got <- tryCatch(
      fit_one(dat, fixed = "nonexistent_var", random = ""),
      fit_one_formula_vars = function(e) TRUE,
      error = function(e) FALSE
    )
    expect_true(got)
  })

  test_that("fit_one_missing_lod fires when tobit sees cens_lod=TRUE but lod=NA", {
    skip_if_not_installed("survival")
    dat <- make_minimal_dat()
    dat$cens_lod <- TRUE
    dat$lod      <- NA_real_
    got <- tryCatch(
      fit_one(dat, family = "tobit", random = ""),
      fit_one_missing_lod = function(e) TRUE,
      error = function(e) FALSE
    )
    expect_true(got)
  })

  test_that("fit_one_missing_ulod fires when tobit sees cens_ulod=TRUE but ulod=NA", {
    skip_if_not_installed("survival")
    dat <- make_minimal_dat()
    dat$cens_ulod <- TRUE
    dat$ulod      <- NA_real_
    got <- tryCatch(
      fit_one(dat, family = "tobit", random = "", ulod = TRUE),
      fit_one_missing_ulod = function(e) TRUE,
      error = function(e) FALSE
    )
    expect_true(got)
  })

  test_that("fit_one_censreg_row_lod fires on row-specific LODs", {
    skip_if_not_installed("censReg")
    dat <- make_minimal_dat()
    dat$cens_lod <- TRUE
    dat$lod      <- rep(c(log1p(0.3), log1p(0.6)), length.out = nrow(dat))
    got <- tryCatch(
      fit_one(dat, family = "tobit_censreg", random = ""),
      fit_one_censreg_row_lod = function(e) TRUE,
      error = function(e) FALSE
    )
    expect_true(got)
  })

  test_that("fit_one_invalid_bounds fires on inverted interval censoring", {
    dat <- make_minimal_dat()
    dat$cens_lod  <- TRUE
    dat$cens_ulod <- TRUE
    dat$lod       <- 100    # inverted: lod > ulod
    dat$ulod      <- 0.5
    got <- tryCatch(
      fit_one(dat, family = "tobit", random = ""),
      fit_one_invalid_bounds = function(e) TRUE,
      error = function(e) FALSE
    )
    expect_true(got)
  })

  test_that("fit_one_error parent class catches any input-validation failure", {
    # Unified parent class lets callers handle all fit_one input failures
    # with a single handler when finer-grained routing is unnecessary.
    caught_class <- tryCatch(
      fit_one(list(value = 1:5)),  # triggers fit_one_not_dataframe
      fit_one_error = function(e) class(e)[[1]]
    )
    expect_equal(caught_class, "fit_one_not_dataframe")

    caught_class <- tryCatch(
      fit_one(data.frame(value = 1:5)),  # triggers fit_one_missing_cols
      fit_one_error = function(e) class(e)[[1]]
    )
    expect_equal(caught_class, "fit_one_missing_cols")
  })

  test_that("fit_one_missing_package class exists and is reachable via class tag", {
    # Exercising the glmmTMB/censReg absence directly would require teardown
    # of the loaded namespace. Instead verify the condition object honors the
    # class when constructed via .fit_one_stop(), matching the dispatch
    # contract fit_models relies on.
    cond <- tryCatch(
      immunoPlex:::.fit_one_stop("fit_one_missing_package", "test"),
      fit_one_missing_package = function(e) e
    )
    expect_s3_class(cond, "fit_one_missing_package")
    expect_s3_class(cond, "fit_one_error")
    expect_s3_class(cond, "error")
  })

  test_that("partial failure populates $failures on the returned model set", {
    # Gap 2 anchor test: pair a family that must succeed (gamma) with one that
    # must fail (censReg on row-specific LODs), and confirm the post-mortem
    # reason for the failing family survives on the return object.
    skip_if_not_installed("glmmTMB")
    skip_if_not_installed("censReg")

    dat <- setup_test_data()
    dat$lod[1:5]      <- dat$lod[1:5] * 0.5
    dat$cens_lod[1:5] <- TRUE

    models <- fit_models(dat,
                         families = c("gamma", "tobit_censreg"),
                         quiet    = TRUE,
                         random   = "")

    expect_true("failures" %in% names(models))
    expect_true(is.character(models$failures))
    expect_true(length(models$failures) >= 1)
    expect_true("tobit_censreg" %in% names(models$failures))
    expect_match(models$failures[["tobit_censreg"]],
                 "single LOD|row-specific LODs")
  })

  test_that("total failure throws fit_models_all_failed carrying failures", {
    # Gap 2 anchor test: the total-failure path must preserve per-family
    # diagnostics on the condition so callers can recover what went wrong.
    caught <- tryCatch(
      fit_models(data.frame(), quiet = TRUE),
      fit_models_all_failed = function(e) e
    )

    expect_s3_class(caught, "fit_models_all_failed")
    expect_s3_class(caught, "error")
    expect_true(is.character(caught$failures))
    expect_true(length(caught$failures) >= 1)
    # At least one family must have recorded the "missing required columns"
    # reason — that is the class that fires on data.frame().
    expect_true(any(grepl("Missing required columns", caught$failures)))
  })

  test_that("fit_models_all_failed preserves the legacy error message", {
    # Backwards-compat guard: pre-existing test suites match on
    # "No models converged successfully" — the class upgrade must not
    # change the user-visible text.
    expect_error(
      fit_models(data.frame(), quiet = TRUE),
      "No models converged successfully"
    )
  })

})

describe("Object-oriented interface implementation", {

test_that("print method generates appropriate output", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  
  dat <- setup_test_data()
  models <- fit_models(dat, families = c("gamma", "tobit"), random = "")
  
  # Should print without error
  expect_output(print(models), "immuno_model_set")
  expect_output(print(models), "Models fitted")
  expect_output(print(models), "Best model")
  expect_output(print(models), "Model Comparison")
})

test_that("summary method provides comprehensive model assessment", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  
  dat <- setup_test_data()
  models <- fit_models(dat, families = c("gamma", "tobit"), random = "")
  
  # Should summarize without error
  expect_output(summary(models), "immuno_model_set Summary")
  expect_output(summary(models), "Number of models")
  expect_output(summary(models), "Best Model Details")
})

test_that("visualization method produces valid statistical graphics", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  skip_if_not_installed("ggplot2")
  
  dat <- setup_test_data()
  models <- fit_models(dat, families = c("gamma", "tobit"), random = "")
  
  # Test comparison plot (default)
  p <- plot(models, plot_type = "comparison", save_pdf = FALSE)
  expect_s3_class(p, "gg")
  
  # Test best model residuals
  expect_message(
    plot(models, plot_type = "best_residuals", save_pdf = FALSE),
    "Plotting residuals for best model"
  )
  
  # Test all residuals
  expect_message(
    plot(models, plot_type = "all_residuals", save_pdf = FALSE),
    "Plotting residuals for"
  )
})

})