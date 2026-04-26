# Statistical Model Estimation Protocol
#
# Validates estimation procedures and convergence criteria:
# 1. Distribution Family Implementation:
#    - Gamma regression
#    - Tobit regression
#    - Accelerated Failure Time models
#
# 2. Detection Limit Processing:
#    - Threshold-based imputation
#    - Censoring-aware estimation
#    - Boundary value handling
#
# 3. Model Specification:
#    - Fixed effects estimation
#    - Random effects validation
#    - Parameter constraint verification
#
# 4. Diagnostic Evaluation:
#    - Residual analysis
#    - Convergence assessment
#    - Numerical stability verification

library(testthat)

# Initialize validation dataset with defined statistical properties
setup_fit_one_data <- function() {
  # Implement hierarchical data source selection
  return(setup_robust_test_data(prefer_real = TRUE, family_hint = "gamma"))  # Primary: empirical data
}

describe("Single distribution estimation functionality", {

test_that("core estimation protocol executes correctly", {
  skip_if_not_installed("glmmTMB")
  
  dat <- setup_fit_one_data()
  
  # Test gamma family
  fit <- fit_one(dat, family = "gamma", random = "")
  expect_s3_class(fit, "immuno_fit")
  expect_equal(fit$family, "gamma")
  expect_true(is.numeric(fit$n_cens_lod))
  expect_true(is.numeric(fit$n_cens_ulod))
  
  # Test new estimand field
  expect_true("estimand" %in% names(fit))
  expect_equal(fit$estimand, "ratio_of_means")
})

test_that("detection limit imputation protocol functions correctly", {
  skip_if_not_installed("glmmTMB")
  
  dat <- setup_fit_one_data()
  
  # Skip if no censored values in this dataset
  if (sum(dat$cens_lod) == 0) {
    skip("No left-censored values found in test data")
  }
  
  # Without imputation
  fit_no_imp <- fit_one(dat, family = "gamma", random = "", impute_lod = FALSE)
  
  # With imputation
  fit_imp <- fit_one(dat, family = "gamma", random = "", impute_lod = TRUE)
  
  if (fit_no_imp$converged && fit_imp$converged) {
    # Check that imputation occurred for censored values
    censored_indices <- which(dat$cens_lod)
    expect_true(all(fit_imp$data_used$value[censored_indices] == dat$lod[censored_indices]))
  }
})

test_that("distribution family estimation protocols achieve convergence", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  
  dat <- setup_fit_one_data()
  families <- c("gamma", "tobit", "aft")
  expected_estimands <- c("ratio_of_means", "ratio_of_means", "ratio_of_medians")
  
  for (i in seq_along(families)) {
    fam <- families[i]
    expected_estimand <- expected_estimands[i]
    
    fit <- fit_one(dat, family = fam, random = "")
    expect_s3_class(fit, "immuno_fit")
    expect_equal(fit$family, fam)
    expect_equal(fit$estimand, expected_estimand)
    
    if (fit$converged) {
      expect_true(is.numeric(fit$aic))
      expect_true(is.numeric(fit$bic))
      expect_true(is.numeric(fit$logLik))
    }
  }
})

test_that("censored regression estimation executes with required dependencies", {
  skip_if_not_installed("censReg")
  
  dat <- setup_fit_one_data()
  
  # Skip if no censored values
  if (sum(dat$cens_lod) == 0) {
    skip("No left-censored values found in test data")
  }
  
  fit <- fit_one(dat, family = "tobit_censreg", random = "")
  expect_s3_class(fit, "immuno_fit")
  expect_equal(fit$family, "tobit_censreg")
  
  if (fit$converged) {
    expect_true(is.numeric(fit$aic))
    expect_true(is.numeric(fit$bic))
    expect_true(is.numeric(fit$logLik))
    expect_s3_class(fit$model, "censReg")
  }
})

test_that("system manages missing dependency requirements appropriately", {
  # Skip test if censReg is available (can't easily mock missing package in R)
  skip_if(requireNamespace("censReg", quietly = TRUE), 
          "censReg package is available - cannot test missing package scenario")
  
  dat <- setup_fit_one_data()
  expect_error(
    fit_one(dat, family = "tobit_censreg", random = ""),
    "censReg package required"
  )
})

test_that("automatic distribution selection protocol executes correctly", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
  skip_if_not_installed("moments")
  
  dat <- setup_fit_one_data()
  
  fit <- fit_one(dat, family = "auto", random = "")
  expect_s3_class(fit, "immuno_fit")
  expect_true(fit$family %in% c("gamma", "tobit", "aft", "gaussian"))

  # Test improved auto-selector logic
  if (any(dat$cens_lod | dat$cens_ulod)) {
    # With censoring, should not select gaussian
    expect_true(fit$family %in% c("tobit", "aft"))
  } else {
    # Without censoring, gaussian is preferred
    expect_true(fit$family %in% c("gamma", "tobit", "aft", "gaussian"))
  }
})

test_that("random effects parameter validation protocol executes correctly", {
  skip_if_not_installed("glmmTMB")
  
  dat <- setup_fit_one_data()
  
  # Should drop random effects due to small sample size (< 30 subjects threshold)
  expect_message(
    fit_one(dat, family = "gamma", random = "(1|subject_id)"),
    "Dropping random term"
  )
})

test_that("AFT model manages non-positive observations appropriately", {
  skip_if_not_installed("survival")
  
  dat <- setup_fit_one_data()
  dat$value[1:3] <- 0  # Add some zeros
  
  fit <- fit_one(dat, family = "aft", random = "")
  expect_s3_class(fit, "immuno_fit")
  
  # AFT should now use proper bounds instead of zero rescue
  if (fit$converged) {
    expect_true(is.numeric(fit$aic))
  }
})

test_that("interval coding works correctly for censored models", {
  skip_if_not_installed("survival")
  
  dat <- setup_fit_one_data()
  
  # Skip if no censored values
  if (sum(dat$cens_lod) == 0) {
    skip("No left-censored values found in test data")
  }
  
  # Test Tobit with proper interval coding
  fit_tobit <- fit_one(dat, family = "tobit", random = "")
  expect_s3_class(fit_tobit, "immuno_fit")
  expect_equal(fit_tobit$estimand, "ratio_of_means")
  
  # Test AFT with proper interval coding  
  fit_aft <- fit_one(dat, family = "aft", random = "")
  expect_s3_class(fit_aft, "immuno_fit")
  expect_equal(fit_aft$estimand, "ratio_of_medians")
})

test_that("censReg validation works correctly", {
  skip_if_not_installed("censReg")
  
  dat <- setup_fit_one_data()
  
  # Skip if no censored values
  if (sum(dat$cens_lod) == 0) {
    skip("No left-censored values found in test data")
  }
  
  # Create data with row-specific LODs (should fail for censReg)
  if (length(unique(dat$lod[dat$cens_lod])) > 1) {
    expect_error(
      fit_one(dat, family = "tobit_censreg", random = ""),
      "censReg requires a single LOD"
    )
  } else {
    # If LODs are already consistent, test should pass
    fit <- fit_one(dat, family = "tobit_censreg", random = "")
    expect_s3_class(fit, "immuno_fit")
    expect_equal(fit$estimand, "ratio_of_means")
  }
})

test_that("ulod parameter works correctly", {
  skip_if_not_installed("survival")
  
  dat <- setup_fit_one_data()
  
  # Add some right-censored values for testing
  dat$cens_ulod[1:2] <- TRUE
  dat$ulod[1:2] <- max(dat$value, na.rm = TRUE) * 1.5
  
  # Test with ulod = TRUE
  fit_with_ulod <- fit_one(dat, family = "tobit", ulod = TRUE, random = "")
  expect_s3_class(fit_with_ulod, "immuno_fit")
  
  # Test with ulod = FALSE (should ignore right censoring)
  fit_without_ulod <- fit_one(dat, family = "tobit", ulod = FALSE, random = "")
  expect_s3_class(fit_without_ulod, "immuno_fit")
})

test_that("fit_one handles random=NULL without error", {
  skip_if_not_installed("glmmTMB")

  dat <- setup_fit_one_data()

  # random=NULL should be treated as "no random effects" (same as random="")
  fit <- fit_one(dat, family = "gamma", random = NULL)
  expect_s3_class(fit, "immuno_fit")
  expect_equal(fit$family, "gamma")
})

test_that("system validates input parameter specifications", {
  dat <- setup_fit_one_data()
  
  # Missing required columns
  dat_bad <- dat[, !names(dat) %in% "cens_lod"]
  expect_error(fit_one(dat_bad), "cens_lod")
  
  # Invalid family
  expect_error(fit_one(dat, family = "invalid"), "should be one of")
})

test_that("convergence indicator functions correctly", {
  skip_if_not_installed("glmmTMB")
  
  dat <- setup_fit_one_data()
  
  fit <- fit_one(dat, family = "gamma", random = "")
  expect_type(fit$converged, "logical")
  
  if (fit$converged) {
    expect_false(inherits(fit$model, "error"))
    expect_true(is.numeric(fit$aic))
  }
})

test_that("core visualization protocol executes correctly", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("ggplot2")
  
  dat <- setup_fit_one_data()
  fit <- fit_one(dat, family = "gamma", random = "")
  
  if (fit$converged) {
    # Test basic plotting
    result <- plot(fit, save_pdf = FALSE)
    expect_type(result, "list")
    expect_true("combined_plot" %in% names(result))
    expect_true("rvf_plot" %in% names(result))
    expect_true("qq_plot" %in% names(result))
    expect_s3_class(result$rvf_plot, "gg")
    expect_s3_class(result$qq_plot, "gg")
  }
})

test_that("censoring-aware visualization protocol executes correctly", {
  skip_if_not_installed("survival")
  skip_if_not_installed("ggplot2")
  
  dat <- setup_fit_one_data()
  
  # Skip if no censored values
  if (sum(dat$cens_lod) == 0) {
    skip("No left-censored values found in test data")
  }
  
  fit <- fit_one(dat, family = "tobit", random = "")
  
  if (fit$converged) {
    result <- plot(fit, plot_type = "censor_aware", save_pdf = FALSE)
    expect_type(result, "list")
    expect_s3_class(result$rvf_plot, "gg")
    expect_s3_class(result$qq_plot, "gg")
  }
})

test_that("simulation-based diagnostic protocol executes correctly", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("DHARMa")
  
  dat <- setup_fit_one_data()
  fit <- fit_one(dat, family = "gamma", random = "")
  
  if (fit$converged) {
    result <- plot(fit, dharma = TRUE, nsim = 100, save_pdf = FALSE)
    expect_type(result, "list")
    expect_true("dharma_obj" %in% names(result))
    # DHARMa object might be NULL if computation fails, so we don't require it
  }
})

test_that("visualization parameter specifications function correctly", {
  skip_if_not_installed("glmmTMB")
  
  dat <- setup_fit_one_data()
  fit <- fit_one(dat, family = "gamma", random = "")
  
  if (fit$converged) {
    # Test different parameter combinations
    result1 <- plot(fit, point_size = 2.0, alpha_cens = 0.8, save_pdf = FALSE)
    result2 <- plot(fit, plot_type = "basic", save_pdf = FALSE)
    
    expect_type(result1, "list")
    expect_type(result2, "list")
    expect_s3_class(result1$rvf_plot, "gg")
    expect_s3_class(result2$rvf_plot, "gg")
  }
})

test_that("visualization system manages non-convergent models appropriately", {
  skip_if_not_installed("glmmTMB")

  dat <- setup_fit_one_data()

  # Create a fit object with converged = FALSE for testing
  fit <- fit_one(dat, family = "gamma", random = "")
  fit$converged <- FALSE

  expect_message(plot(fit), "Cannot plot: model did not converge")
  expect_null(plot(fit))
})

})


# ---------------------------------------------------------------------------
# rep_col / plate_col: replicate- and plate-level random effects
# ---------------------------------------------------------------------------

# Inline replicate-level DGP kept local to this test block so the new tests do
# not depend on any shared helper. Layout: n_subj subjects x n_tp timepoints x
# n_rep technical replicates, with subject, nested (subject:rep) and plate
# random effects baked in.
.make_repdat <- function(n_subj = 40, n_tp = 2, n_rep = 3, n_plates = 4,
                         seed = 1234,
                         sigma_subj = 0.4, sigma_rep = 0.35,
                         sigma_plate = 0.3, sigma_err = 0.2) {
  set.seed(seed)
  n_obs <- n_subj * n_tp * n_rep
  subj_id   <- rep(seq_len(n_subj), each = n_tp * n_rep)
  timepoint <- rep(rep(paste0("T", seq_len(n_tp)), each = n_rep), times = n_subj)
  rep_id    <- rep(paste0("R", seq_len(n_rep)), times = n_subj * n_tp)
  disease   <- rep(rep(c("Case", "Control"), length.out = n_subj),
                   each = n_tp * n_rep)
  age       <- rep(round(rnorm(n_subj, 30, 5)), each = n_tp * n_rep)
  # plate assigned per (subject, timepoint) cell — plates cross subjects
  plate     <- rep(paste0("P", sample(seq_len(n_plates), n_subj * n_tp, TRUE)),
                   each = n_rep)

  subj_re  <- rnorm(n_subj, 0, sigma_subj)
  rep_key  <- paste0(subj_id, ":", rep_id)
  rep_re   <- setNames(rnorm(length(unique(rep_key)), 0, sigma_rep),
                       unique(rep_key))
  plate_re <- setNames(rnorm(length(unique(plate)), 0, sigma_plate),
                       unique(plate))

  mu <- 5 + (timepoint == "T2") * 0.3 + (disease == "Case") * 0.25 +
        0.02 * (age - 30) + subj_re[subj_id] +
        rep_re[rep_key] + plate_re[plate]
  value <- mu + rnorm(n_obs, 0, sigma_err)

  data.frame(
    subject_id = factor(subj_id),
    timepoint  = factor(timepoint),
    rep_id     = factor(rep_id),
    disease    = factor(disease),
    age        = age,
    plate      = factor(plate),
    value      = value,
    cens_lod   = FALSE,
    cens_ulod  = FALSE,
    lod        = 0.1,
    ulod       = 100,
    stringsAsFactors = FALSE
  )
}

test_that("fit_one: rep_col = NULL reproduces existing behavior exactly", {
  skip_if_not_installed("glmmTMB")
  dat <- .make_repdat()

  base <- fit_one(dat, family = "gaussian",
                  fixed = "timepoint + disease + age",
                  random = "(1|subject_id)")
  new  <- fit_one(dat, family = "gaussian",
                  fixed = "timepoint + disease + age",
                  random = "(1|subject_id)",
                  rep_col = NULL, plate_col = NULL)

  expect_true(base$converged)
  expect_true(new$converged)
  expect_equal(glmmTMB::fixef(base$model)$cond,
               glmmTMB::fixef(new$model)$cond)
  expect_equal(base$logLik, new$logLik)
})

test_that("fit_one: rep_col appends nested (subject_id:rep_col) random intercept", {
  skip_if_not_installed("glmmTMB")
  dat <- .make_repdat()
  fit <- fit_one(dat, family = "gaussian",
                 fixed = "timepoint + disease + age",
                 random = "(1|subject_id)",
                 rep_col = "rep_id")

  expect_true(fit$converged)

  vc <- glmmTMB::VarCorr(fit$model)$cond
  expect_true("subject_id:rep_id" %in% names(vc))
  sigma_rep <- attr(vc[["subject_id:rep_id"]], "stddev")
  expect_true(is.numeric(sigma_rep) && sigma_rep > 0)
})

test_that("fit_one: plate_col appends flat (plate) random intercept without altering fixed effects much", {
  skip_if_not_installed("glmmTMB")
  dat <- .make_repdat()

  base <- fit_one(dat, family = "gaussian",
                  fixed = "timepoint + disease + age",
                  random = "(1|subject_id)")
  fit  <- fit_one(dat, family = "gaussian",
                  fixed = "timepoint + disease + age",
                  random = "(1|subject_id)",
                  plate_col = "plate")

  expect_true(fit$converged)

  vc <- glmmTMB::VarCorr(fit$model)$cond
  expect_true("plate" %in% names(vc))
  sigma_plate <- attr(vc[["plate"]], "stddev")
  expect_true(is.numeric(sigma_plate) && sigma_plate > 0)

  # Fixed-effect point estimates should not move materially
  b_base <- glmmTMB::fixef(base$model)$cond
  b_fit  <- glmmTMB::fixef(fit$model)$cond
  # Compare on shared terms; both models have identical fixed part
  shared <- intersect(names(b_base), names(b_fit))
  expect_true(max(abs(b_base[shared] - b_fit[shared])) < 0.1)
})

test_that("fit_one: rep_col = nonexistent column errors clearly", {
  skip_if_not_installed("glmmTMB")
  dat <- .make_repdat()
  expect_error(
    fit_one(dat, family = "gaussian", rep_col = "no_such_col"),
    "`rep_col` not found in data"
  )
})

test_that("fit_one: plate_col = nonexistent column errors clearly", {
  skip_if_not_installed("glmmTMB")
  dat <- .make_repdat()
  expect_error(
    fit_one(dat, family = "gaussian", plate_col = "no_such_col"),
    "`plate_col` not found in data"
  )
})

test_that("fit_one: non-character rep_col / plate_col errors clearly", {
  skip_if_not_installed("glmmTMB")
  dat <- .make_repdat()
  expect_error(fit_one(dat, family = "gaussian", rep_col = 42),
               "must be a single column name")
  expect_error(fit_one(dat, family = "gaussian", plate_col = c("a", "b")),
               "must be a single column name")
})

test_that("fit_models: forwards rep_col / plate_col through to fit_one", {
  skip_if_not_installed("glmmTMB")
  dat <- .make_repdat()

  res <- fit_models(dat, families = "gaussian",
                    fixed = "timepoint + disease + age",
                    random = "(1|subject_id)",
                    rep_col = "rep_id", plate_col = "plate",
                    quiet = TRUE)
  expect_s3_class(res, "immuno_model_set")
  best <- res$best_model
  expect_true(best$converged)

  vc <- glmmTMB::VarCorr(best$model)$cond
  expect_true("subject_id:rep_id" %in% names(vc))
  expect_true("plate" %in% names(vc))
})

test_that("fit_models: rep_col = nonexistent column surfaces the error via fit_one", {
  skip_if_not_installed("glmmTMB")
  dat <- .make_repdat()
  # fit_models swallows fit_one errors with tryCatch, but the per-family fit
  # fails to converge, so no models converge → fit_models itself errors.
  expect_error(
    fit_models(dat, families = "gaussian",
               fixed = "timepoint + disease + age",
               random = "(1|subject_id)",
               rep_col = "no_such_col", quiet = TRUE),
    "No models converged"
  )
})