# Integration tests: verify that analysis functions recover known ground-truth
# effects from simulated data at large n.
#
# Each test generates data with a known effect, runs the analysis, and checks
# that the estimate is within 2 SE of the true value (or meets an equivalent
# acceptance criterion).


# ---- Test 1: simulate_immunoassay() -> fit_one() recovers group effect --------

test_that("fit_one() recovers known group effect from simulated data", {
  true_effect <- 0.8

  sim <- simulate_immunoassay(
    n_subjects       = 500,
    n_timepoints     = 1,
    n_analytes       = 1,
    design           = "cross_sectional",
    group_levels     = c("control", "treatment"),
    group_effects    = true_effect,
    signal_analytes  = 1,
    effect_direction = "up",
    re_intercept_sd  = 0.3,
    residual_sd      = 0.8,
    lod_values       = c("IL-1b" = 0.0001),
    seed             = 12345
  )

  # Verify the DGP encoded the expected effect
  expect_equal(sim$meta$dgp_params$group_effects[1], true_effect)

  d <- as.data.frame(sim$data[!is.na(sim$data$value), ])

  # Gaussian model, no random effects (1 obs per subject)
  result <- fit_one(d, family = "gaussian", fixed = "group", random = "")

  expect_s3_class(result, "immuno_fit")
  expect_true(result$converged)

  # Extract estimated group effect and SE
  coef_tbl <- summary(result$model)$coefficients$cond
  est <- coef_tbl["grouptreatment", "Estimate"]
  se  <- coef_tbl["grouptreatment", "Std. Error"]

  # Ground truth recovery: estimated effect within 2 SE of true effect
  expect_true(
    abs(est - true_effect) < 2 * se,
    info = sprintf("est=%.4f, true=%.4f, 2*SE=%.4f, gap=%.4f",
                   est, true_effect, 2 * se, abs(est - true_effect))
  )
})


# ---- Test 2: simulate_ancova_data() -> ancova_one() recovers group effect ----

test_that("ancova_one() recovers group effect and omega-squared from simulated data", {
  true_effect <- 0.8


  sim <- simulate_ancova_data(
    n_subjects          = 200,
    n_analytes          = 1,
    signal_analytes     = 1,
    group_effects       = true_effect,
    covariate_group_cor = 0.0,
    outlier_rate        = 0,
    seed                = 23456
  )

  d <- as.data.frame(sim$data[!is.na(sim$data$value), ])

  result <- ancova_one(
    data             = d,
    outcome          = "value",
    group            = "group",
    covariate        = "age",
    covariates       = "bmi",
    outlier_removal  = FALSE,
    log_validate     = TRUE
  )

  expect_s3_class(result, "immuno_ancova")

  # Rank-based checks: group effect detected, correct direction, omega > 0

  expect_true(result$group_p < 0.05,
    info = sprintf("group_p=%.4f, expected < 0.05", result$group_p))
  expect_true(result$group_effect > 0,
    info = "Rank-based group effect should be positive (effect_direction = up)")
  expect_true(result$omega_sq_partial > 0,
    info = sprintf("omega_sq_partial=%.4f, expected > 0", result$omega_sq_partial))

  # Log-scale validation: log_effect within 2 SE of true effect
  expect_false(is.null(result$log_model),
    info = "log_model should exist when log_validate = TRUE")

  log_coefs <- summary(result$log_model)$coefficients
  term_name <- paste0("group", result$group_levels[2])
  expect_true(term_name %in% rownames(log_coefs),
    info = sprintf("Expected term '%s' in log model coefficients", term_name))

  log_est <- log_coefs[term_name, "Estimate"]
  log_se  <- log_coefs[term_name, "Std. Error"]

  expect_true(
    abs(log_est - true_effect) < 2 * log_se,
    info = sprintf("log_est=%.4f, true=%.4f, 2*SE=%.4f, gap=%.4f",
                   log_est, true_effect, 2 * log_se, abs(log_est - true_effect))
  )
})


# ---- Test 3: simulate_mcnemar_data() -> mcnemar_detection() recovers delta ---

test_that("mcnemar_detection() recovers detection delta from simulated data", {
  target_disc <- 0.20

  sim <- simulate_mcnemar_data(
    n_subjects         = 500,
    n_analytes         = 5,
    baseline_detection = 0.5,
    discordant_rate    = target_disc,
    signal_analytes    = 2,
    seed               = 34567
  )

  result <- mcnemar_detection(
    data          = sim$data,
    subject_col   = "subject_id",
    cytokine_col  = "cytokine",
    timepoint_col = "timepoint",
    censoring_col = "cens_lod",
    baseline      = "T1",
    comparison    = "T2",
    quiet         = TRUE
  )

  signal_names <- sim$meta$dgp_params$signal_analyte_names
  null_names   <- setdiff(unique(sim$data$cytokine), signal_names)

  # Signal analytes: observed discordant fraction close to target
  for (nm in signal_names) {
    row <- result[result$cytokine == nm, ]
    obs_disc <- row$n_discordant / row$n_pairs
    se_disc  <- sqrt(obs_disc * (1 - obs_disc) / row$n_pairs)
    expect_true(
      abs(obs_disc - target_disc) < 2 * se_disc + 0.02,
      info = sprintf("%s: obs disc=%.3f, target=%.3f, 2*SE=%.3f",
                     nm, obs_disc, target_disc, 2 * se_disc)
    )
  }

  # Symmetric discordance model => true delta = 0 => CI should contain 0
  for (nm in signal_names) {
    row <- result[result$cytokine == nm, ]
    expect_true(
      row$delta_ci_lo <= 0 && row$delta_ci_hi >= 0,
      info = sprintf("%s: delta=%.1f, CI=[%.1f, %.1f] should contain 0",
                     nm, row$delta_detection, row$delta_ci_lo, row$delta_ci_hi)
    )
  }

  # Null analytes: discordant rate should be near 0
  for (nm in null_names) {
    row <- result[result$cytokine == nm, ]
    obs_disc <- row$n_discordant / row$n_pairs
    expect_true(
      obs_disc < 0.10,
      info = sprintf("%s (null): disc rate=%.3f, expected < 0.10", nm, obs_disc)
    )
  }
})


# ---- Test 3b: independent-DGP recovery (review plan gap, not blocking) -------
#
# The previous test uses simulate_mcnemar_data() from the same package as
# mcnemar_detection(). That is co-authored validation: a bug shared by DGP
# and fitter would cancel out. This test uses an independent, self-contained
# bivariate-normal-threshold generator defined inline (no simulate_mcnemar_data,
# no simulate_immunoassay) so the DGP and the fitter cannot share an
# implementation bug.

test_that("mcnemar_detection() recovers a known paired-proportion difference from an independent DGP", {
  skip_if_not_installed("exact2x2")

  # Independent DGP: bivariate normal threshold model.
  #   (Z1, Z2) ~ MVN(0, [[1, rho_t], [rho_t, 1]])
  #   detected_T1_k = Z1[, k] > t1_k
  #   detected_T2_k = Z2[, k] > t2_k    (t2_k < t1_k for signal analytes)
  #
  # Theoretical true Delta = P(Z > t2) - P(Z > t1) on a standard normal.
  set.seed(20260423L)
  n_subjects <- 800L
  analyte_names <- c("S1", "S2", "N1", "N2")   # 2 signal, 2 null
  is_signal <- c(TRUE, TRUE, FALSE, FALSE)
  baseline_detect <- 0.50                       # p1 -> t1 = qnorm(1 - p1)
  delta_true      <- 0.15                       # p2 - p1 for signal
  rho_t           <- 0.4

  t1 <- stats::qnorm(1 - baseline_detect)
  t2_sig  <- stats::qnorm(1 - (baseline_detect + delta_true))
  t2_null <- t1

  n_ana <- length(analyte_names)
  Z1 <- matrix(stats::rnorm(n_subjects * n_ana), n_subjects, n_ana)
  Z2_innov <- matrix(stats::rnorm(n_subjects * n_ana), n_subjects, n_ana)
  Z2 <- rho_t * Z1 + sqrt(1 - rho_t^2) * Z2_innov

  det_T1 <- Z1 > t1
  det_T2 <- matrix(FALSE, n_subjects, n_ana)
  for (k in seq_len(n_ana)) {
    thr_k <- if (is_signal[k]) t2_sig else t2_null
    det_T2[, k] <- Z2[, k] > thr_k
  }

  # Assemble long-format data for mcnemar_detection()
  long_dat <- data.frame(
    subject_id = rep(seq_len(n_subjects), times = n_ana * 2L),
    cytokine   = rep(rep(analyte_names, each = n_subjects), times = 2L),
    timepoint  = rep(c("T1", "T2"), each = n_subjects * n_ana),
    cens_lod   = c(as.vector(!det_T1), as.vector(!det_T2)),
    stringsAsFactors = FALSE
  )

  result <- mcnemar_detection(
    data        = long_dat,
    baseline    = "T1",
    comparison  = "T2",
    quiet       = TRUE
  )

  # Signal analytes: Delta point estimate near true 15pp (within ~2*SE on a
  # proportion at n=800 ~= 2 * 100*sqrt(0.5*0.5/800) ~= 3.5pp + margin).
  for (nm in analyte_names[is_signal]) {
    row <- result[result$cytokine == nm, ]
    obs_delta_pp <- row$delta_detection   # already in percentage points
    expect_true(
      abs(obs_delta_pp - 100 * delta_true) < 5,
      info = sprintf("%s (signal): delta=%.2fpp, truth=%.2fpp",
                     nm, obs_delta_pp, 100 * delta_true)
    )
    # Exact Delta CI should cover the true 15pp
    expect_true(
      row$delta_ci_lo <= 100 * delta_true && row$delta_ci_hi >= 100 * delta_true,
      info = sprintf("%s: CI=[%.1f, %.1f] should cover truth %.1f",
                     nm, row$delta_ci_lo, row$delta_ci_hi, 100 * delta_true)
    )
    # Power at n=800, delta=0.15 is essentially 1. Should reject.
    expect_true(row$p_mcnemar < 0.01,
                info = sprintf("%s: p=%.4g should be small", nm, row$p_mcnemar))
  }

  # Null analytes: true Delta = 0. CI should cover 0; p should not reject.
  for (nm in analyte_names[!is_signal]) {
    row <- result[result$cytokine == nm, ]
    expect_true(
      row$delta_ci_lo <= 0 && row$delta_ci_hi >= 0,
      info = sprintf("%s (null): CI=[%.1f, %.1f] should cover 0",
                     nm, row$delta_ci_lo, row$delta_ci_hi)
    )
    expect_gt(row$p_mcnemar, 0.01)
  }
})


# ---- Test 4: simulate_plsda_data() -> PLS-DA VIP ranks signal > noise --------

test_that("PLS-DA VIP correctly ranks signal vs. noise analytes (AUROC > 0.8)", {
  skip_if_not_installed("ropls")

  sim <- simulate_plsda_data(
    n_subjects      = 200,
    n_analytes      = 15,
    n_discriminatory = 5,
    group_effects   = 1.0,
    analyte_correlation = NULL,
    seed            = 45678
  )

  signal_names  <- sim$meta$dgp_params$signal_analyte_names
  analyte_names <- sim$meta$dgp_params$analyte_names

  # Pivot to wide format for plsda_preprocess
  d <- sim$data[!is.na(sim$data$value), ]
  wide <- tidyr::pivot_wider(
    d[, c("subject_id", "group", "cytokine", "value")],
    names_from  = "cytokine",
    values_from = "value"
  )

  lod_lookup <- data.frame(
    cytokine = analyte_names,
    lod = vapply(analyte_names, function(nm) {
      d$lod[d$cytokine == nm][1]
    }, numeric(1)),
    stringsAsFactors = FALSE
  )

  pp <- plsda_preprocess(
    data          = as.data.frame(wide),
    cytokine_cols = analyte_names,
    metadata_cols = c("subject_id", "group"),
    lod_lookup    = lod_lookup,
    lod_method    = "half",
    log_transform = TRUE,
    scale_data    = TRUE,
    verbose       = FALSE
  )

  fit <- plsda_fit(
    preprocessed_data = pp,
    response_var      = "group",
    n_components      = 2,
    method            = "PLS-DA",
    permutations      = 0,
    verbose           = FALSE
  )

  expect_s3_class(fit, "plsda_model")

  vip <- fit$vip_scores
  vip_signal <- vip$vip_score[vip$cytokine %in% signal_names]
  vip_noise  <- vip$vip_score[!vip$cytokine %in% signal_names]

  # AUROC: probability that a signal analyte has higher VIP than a noise analyte
  auroc <- mean(vapply(vip_signal, function(s) {
    mean(s > vip_noise) + 0.5 * mean(s == vip_noise)
  }, numeric(1)))

  expect_true(
    auroc > 0.8,
    info = sprintf("VIP AUROC=%.3f, expected > 0.8", auroc)
  )
})


# ---- Test 5: fit_one(family="aft") recovers group effect with censoring ------

test_that("fit_one(family='aft') recovers known group effect from censored data", {
  skip_if_not_installed("survival")

  true_effect <- 0.8

  sim <- simulate_immunoassay(
    n_subjects       = 500,
    n_timepoints     = 1,
    n_analytes       = 1,
    design           = "cross_sectional",
    group_levels     = c("control", "treatment"),
    group_effects    = true_effect,
    signal_analytes  = 1,
    effect_direction = "up",
    re_intercept_sd  = 0,
    residual_sd      = 0.8,
    lod_quantile     = 0.20,
    seed             = 54321
  )

  d <- as.data.frame(sim$data[!is.na(sim$data$value), ])
  # Confirm value_raw is present (simulate_immunoassay always provides it)
  expect_true("value_raw" %in% names(d))

  result <- fit_one(d, family = "aft", fixed = "group", random = "")

  expect_s3_class(result, "immuno_fit")
  expect_true(result$converged)
  expect_equal(result$family, "aft")
  expect_equal(result$estimand, "ratio_of_medians")

  # Extract survreg coefficient for group effect
  coef_tab <- summary(result$model)$table
  group_row <- grep("^group", rownames(coef_tab))
  expect_true(length(group_row) > 0,
    info = "Expected 'group*' row in survreg coefficient table")

  est <- coef_tab[group_row[1], "Value"]
  se  <- coef_tab[group_row[1], "Std. Error"]

  # Ground truth: estimate within 2 SE of true effect
  expect_true(
    abs(est - true_effect) < 2 * se,
    info = sprintf("AFT: est=%.4f, true=%.4f, 2*SE=%.4f, gap=%.4f",
                   est, true_effect, 2 * se, abs(est - true_effect))
  )
})
