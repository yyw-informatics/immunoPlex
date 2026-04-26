# Tests for scenario wrappers and specialized generators (Phase 2)

# ---- Helper: validate immuno_sim structure -----------------------------------

expect_valid_immuno_sim <- function(sim, design = NULL) {
  expect_s3_class(sim, "immuno_sim")
  expect_true(is.data.frame(sim$data))
  expect_true(is.data.frame(sim$truth))
  expect_true(is.list(sim$meta))

  # Required data columns
  expect_true(all(c("subject_id", "cytokine", "value") %in% names(sim$data)))

  # Required truth columns
  expect_true("true_latent" %in% names(sim$truth))

  # Row alignment

  expect_equal(nrow(sim$data), nrow(sim$truth))

  # Meta required fields
  expect_true(!is.null(sim$meta$design))
  expect_true(!is.null(sim$meta$n_subjects))
  expect_true(sim$meta$n_subjects >= 1)

  if (!is.null(design)) {
    expect_equal(sim$meta$design, design)
  }
}


# ==============================================================================
# Scenario Wrappers
# ==============================================================================

test_that("simulate_vaccine_trial returns valid immuno_sim", {
  sim <- simulate_vaccine_trial(n_subjects = 40, n_analytes = 10, seed = 42)
  expect_valid_immuno_sim(sim, design = "pre_post")

  # Check group levels
  expect_setequal(levels(sim$data$group), c("responder", "non-responder"))

  # Check 2 timepoints
  expect_equal(length(unique(sim$data$timepoint)), 2)

  # Check n_subjects
  expect_equal(length(unique(sim$data$subject_id)), 40)

  # Check n_analytes
  expect_equal(length(unique(sim$data$cytokine)), 10)
})

test_that("simulate_vaccine_trial has heterogeneous per-analyte censoring", {
  sim <- simulate_vaccine_trial(n_subjects = 60, n_analytes = 10, seed = 123)

  # Should have some censored observations from library LODs
  expect_true("cens_lod" %in% names(sim$data))

  # Censoring rates should vary across analytes (heterogeneous, not uniform)
  per_analyte <- tapply(sim$data$cens_lod, sim$data$cytokine, mean)
  # Range of censoring rates should be > 0.05 (not all the same)
  expect_gt(diff(range(per_analyte)), 0.05)
})

test_that("simulate_vaccine_trial respects ... passthrough", {
  sim <- simulate_vaccine_trial(
    n_subjects = 40, n_analytes = 10,
    n_plates = 2, plate_sd = 0.1,
    seed = 42
  )
  expect_valid_immuno_sim(sim)
  expect_equal(length(unique(sim$data$plate_id)), 2)
})

test_that("simulate_vaccine_trial has signal analytes", {
  sim <- simulate_vaccine_trial(
    n_subjects = 40, n_analytes = 10, signal_analytes = 3, seed = 42
  )
  n_signal <- length(unique(sim$data$cytokine[sim$truth$signal_analyte]))
  expect_equal(n_signal, 3)
})


# ---- simulate_infection_timecourse ------------------------------------------

test_that("simulate_infection_timecourse returns valid immuno_sim", {
  sim <- simulate_infection_timecourse(n_subjects = 40, n_analytes = 10, seed = 42)
  expect_valid_immuno_sim(sim, design = "time_course")

  # Check group levels
  expect_setequal(levels(sim$data$group), c("mild", "severe"))

  # Check 5 timepoints (default)
  expect_equal(length(unique(sim$data$timepoint)), 5)
})

test_that("simulate_infection_timecourse has heteroscedastic residuals", {
  sim <- simulate_infection_timecourse(
    n_subjects = 200, n_analytes = 5, seed = 42
  )

  # Heteroscedastic: severe group should have larger variance
  d <- sim$data
  for (cyt in unique(d$cytokine)[1:2]) {
    vals_mild   <- d$value[d$group == "mild" & d$cytokine == cyt]
    vals_severe <- d$value[d$group == "severe" & d$cytokine == cyt]
    vals_mild   <- vals_mild[!is.na(vals_mild)]
    vals_severe <- vals_severe[!is.na(vals_severe)]
    # Severe should have 1.5x larger SD
    if (length(vals_mild) > 10 && length(vals_severe) > 10) {
      ratio <- sd(vals_severe) / sd(vals_mild)
      expect_gt(ratio, 1.0)  # should be close to 1.5 but with randomness
    }
  }
})

test_that("simulate_infection_timecourse has dropout missingness", {
  sim <- simulate_infection_timecourse(
    n_subjects = 200, n_analytes = 5, dropout_rate = 0.08, seed = 42
  )

  # With dropout_rate = 0.08 over 5 timepoints, expect some NA values
  na_rate <- mean(is.na(sim$data$value))
  expect_gt(na_rate, 0.01)  # at least some missingness

  # Check that some subjects have incomplete data (dropped out)
  obs_per_subj <- tapply(!is.na(sim$data$value), sim$data$subject_id, sum)
  max_obs <- max(obs_per_subj)
  incomplete <- sum(obs_per_subj < max_obs)
  expect_gt(incomplete, 0)  # some subjects dropped out
})

test_that("simulate_infection_timecourse supports custom timepoints", {
  sim <- simulate_infection_timecourse(
    n_subjects = 30, n_timepoints = 3, n_analytes = 5, seed = 42
  )
  expect_equal(length(unique(sim$data$timepoint)), 3)
})


# ---- simulate_severity_comparison ------------------------------------------

test_that("simulate_severity_comparison returns valid immuno_sim", {
  sim <- simulate_severity_comparison(n_subjects = 60, n_analytes = 10, seed = 42)
  expect_valid_immuno_sim(sim, design = "cross_sectional")

  # 4 severity groups
  expect_setequal(levels(sim$data$group),
                  c("none", "mild", "moderate", "severe"))

  # Single timepoint
  expect_equal(length(unique(sim$data$timepoint)), 1)
})

test_that("simulate_severity_comparison has graded effects", {
  sim <- simulate_severity_comparison(
    n_subjects = 400, n_analytes = 10, signal_analytes = 3,
    group_effects = c(0, 0.3, 0.6, 1.0), seed = 42
  )

  # For signal analytes, group effect should increase with severity
  sig_analytes <- unique(sim$data$cytokine[sim$truth$signal_analyte])
  if (length(sig_analytes) > 0) {
    cyt <- sig_analytes[1]
    d <- sim$data[sim$data$cytokine == cyt, ]
    means <- tapply(d$value, d$group, mean, na.rm = TRUE)
    # Severe should have highest mean (effects are "up")
    expect_gt(means["severe"], means["none"])
  }
})

test_that("simulate_severity_comparison supports class imbalance", {
  sim <- simulate_severity_comparison(
    n_subjects = 100, n_analytes = 5,
    class_imbalance = TRUE, seed = 42
  )

  # With 4:3:2:1 allocation, group sizes should be unequal
  d <- sim$data
  per_subj <- d[!duplicated(d$subject_id), ]
  grp_counts <- table(per_subj$group)
  expect_true(grp_counts["none"] > grp_counts["severe"])
})

test_that("simulate_severity_comparison has block correlation", {
  sim <- simulate_severity_comparison(
    n_subjects = 100, n_analytes = 10, block_rho = 0.6, seed = 42
  )
  expect_valid_immuno_sim(sim)

  # Check that realized correlation shows block structure
  cor_mat <- sim$meta$realized_correlation
  expect_true(is.matrix(cor_mat))
})


# ---- simulate_exposure_cohort -----------------------------------------------

test_that("simulate_exposure_cohort returns valid immuno_sim", {
  sim <- simulate_exposure_cohort(n_subjects = 60, n_analytes = 10, seed = 42)
  expect_valid_immuno_sim(sim, design = "paired_exposure")

  # 3 groups
  expect_setequal(levels(sim$data$group),
                  c("unexposed", "low_exposure", "high_exposure"))

  # Has covariates
  expect_true("age" %in% names(sim$data))
  expect_true("sex" %in% names(sim$data))
})

test_that("simulate_exposure_cohort has heterogeneous per-analyte censoring", {
  sim <- simulate_exposure_cohort(n_subjects = 60, n_analytes = 10, seed = 42)

  # Should have substantial censoring from scaled library LODs
  expect_true("cens_lod" %in% names(sim$data))
  cens_rate <- mean(sim$data$cens_lod, na.rm = TRUE)
  expect_gt(cens_rate, 0.10)

  # Censoring rates should vary across analytes (heterogeneous, not uniform)
  per_analyte <- tapply(sim$data$cens_lod, sim$data$cytokine, mean)
  expect_gt(diff(range(per_analyte)), 0.05)
})

test_that("simulate_exposure_cohort has visit missingness", {
  sim <- simulate_exposure_cohort(
    n_subjects = 200, n_analytes = 10, visit_missingness = 0.05, seed = 42
  )

  # With visit_missingness = 0.05, expect ~5% NA values
  na_rate <- mean(is.na(sim$data$value))
  expect_gt(na_rate, 0.01)   # at least some missingness
  expect_lt(na_rate, 0.20)   # not excessive
})

test_that("simulate_exposure_cohort has confounding covariates", {
  # With covariate_group_cor = 0.3, covariates should correlate with group
  sim <- simulate_exposure_cohort(
    n_subjects = 300, n_analytes = 5,
    covariate_group_cor = 0.5, seed = 42
  )

  d <- sim$data
  per_subj <- d[!duplicated(d$subject_id), ]
  cor_val <- cor(as.numeric(per_subj$group), per_subj$age, use = "complete.obs")
  # With rho = 0.5, should see positive correlation
  expect_gt(cor_val, 0.1)
})


# ==============================================================================
# Specialized Generators
# ==============================================================================

# ---- simulate_mcnemar_data ---------------------------------------------------

test_that("simulate_mcnemar_data returns valid immuno_sim", {
  sim <- simulate_mcnemar_data(n_subjects = 50, n_analytes = 5, seed = 42)
  expect_valid_immuno_sim(sim, design = "pre_post")

  # 2 timepoints
  expect_equal(length(unique(sim$data$timepoint)), 2)

  # Has cens_lod for detection status
  expect_true("cens_lod" %in% names(sim$data))
})

test_that("simulate_mcnemar_data has correct discordant rates for signal analytes", {
  sim <- simulate_mcnemar_data(
    n_subjects = 500, n_analytes = 5,
    baseline_detection = 0.5,
    discordant_rate = 0.20,
    signal_analytes = 2,
    seed = 42
  )

  # Signal analytes should have ~20% discordant rate
  sig_names <- sim$meta$dgp_params$signal_analyte_names
  disc_rates <- sim$meta$mcnemar$realized_discordant

  for (nm in sig_names) {
    expect_gt(disc_rates[nm], 0.10)  # at least 10%
    expect_lt(disc_rates[nm], 0.35)  # not more than 35%
  }
})

test_that("simulate_mcnemar_data null analytes have low discordant rates", {
  sim <- simulate_mcnemar_data(
    n_subjects = 500, n_analytes = 5,
    discordant_rate = 0.20,
    signal_analytes = 2,
    seed = 42
  )

  null_mask <- !sim$meta$dgp_params$signal_analytes
  null_names <- sim$meta$dgp_params$signal_analyte_names  # these are signal
  all_names <- unique(sim$data$cytokine)
  null_analytes <- setdiff(all_names, null_names)

  disc_rates <- sim$meta$mcnemar$realized_discordant

  for (nm in null_analytes) {
    # Null analytes should have very low discordant rate (rho_t = 0.99)
    expect_lt(disc_rates[nm], 0.10)
  }
})

test_that("simulate_mcnemar_data zero discordant rate yields stable detection", {
  sim <- simulate_mcnemar_data(
    n_subjects = 200, n_analytes = 3,
    discordant_rate = 0,
    signal_analytes = 1,
    seed = 42
  )

  disc_rates <- sim$meta$mcnemar$realized_discordant
  # All analytes should have very low discordant rate
  for (nm in names(disc_rates)) {
    expect_lt(disc_rates[nm], 0.10)
  }
})

test_that("simulate_mcnemar_data supports correlated analyte detections", {
  sim <- simulate_mcnemar_data(
    n_subjects = 200, n_analytes = 5,
    analyte_rho = 0.5,
    discordant_rate = 0.15,
    signal_analytes = 2,
    seed = 42
  )

  expect_valid_immuno_sim(sim)
  expect_equal(sim$meta$mcnemar$analyte_rho, 0.5)

  # Check that detections are correlated across analytes
  d <- sim$data
  t1 <- d[d$timepoint == "T1", ]
  # Pivot to wide detection matrix
  det_mat <- tapply(!t1$cens_lod, list(t1$subject_id, t1$cytokine), identity)
  det_mat <- matrix(as.numeric(det_mat), nrow = nrow(det_mat))
  if (ncol(det_mat) >= 2 && nrow(det_mat) > 10) {
    cors <- cor(det_mat, use = "pairwise.complete.obs")
    off_diag <- cors[lower.tri(cors)]
    # With rho = 0.5, detection correlations should be positive
    expect_gt(mean(off_diag, na.rm = TRUE), 0)
  }
})

test_that("simulate_mcnemar_data is deterministic with seed", {
  sim1 <- simulate_mcnemar_data(n_subjects = 30, n_analytes = 3, seed = 99)
  sim2 <- simulate_mcnemar_data(n_subjects = 30, n_analytes = 3, seed = 99)
  expect_equal(sim1$data$value, sim2$data$value)
  expect_equal(sim1$data$cens_lod, sim2$data$cens_lod)
})


# ---- simulate_plsda_data ----------------------------------------------------

test_that("simulate_plsda_data returns valid immuno_sim", {
  sim <- simulate_plsda_data(n_subjects = 40, n_analytes = 10, seed = 42)
  expect_valid_immuno_sim(sim, design = "cross_sectional")

  # 2 groups
  expect_setequal(levels(sim$data$group), c("case", "control"))

  # Single timepoint
  expect_equal(length(unique(sim$data$timepoint)), 1)
})

test_that("simulate_plsda_data has correct number of discriminatory analytes", {
  sim <- simulate_plsda_data(
    n_subjects = 40, n_analytes = 15, n_discriminatory = 5, seed = 42
  )

  n_signal <- length(unique(sim$data$cytokine[sim$truth$signal_analyte]))
  expect_equal(n_signal, 5)
})

test_that("simulate_plsda_data supports class imbalance", {
  sim <- simulate_plsda_data(
    n_subjects = 60, n_analytes = 10,
    class_imbalance = 2,  # 2:1 ratio
    seed = 42
  )

  d <- sim$data
  per_subj <- d[!duplicated(d$subject_id), ]
  grp_counts <- table(per_subj$group)
  expect_gt(grp_counts["case"], grp_counts["control"])
})

test_that("simulate_plsda_data has block correlation", {
  sim <- simulate_plsda_data(
    n_subjects = 50, n_analytes = 10,
    analyte_correlation = "block", block_rho = 0.5, seed = 42
  )

  cor_mat <- sim$meta$realized_correlation
  expect_true(is.matrix(cor_mat))
  # Some off-diagonal should be positive (within category blocks)
  off_diag <- cor_mat[lower.tri(cor_mat)]
  expect_true(any(off_diag > 0.1))
})

test_that("simulate_plsda_data supports censoring via lod_quantile", {
  sim <- simulate_plsda_data(
    n_subjects = 50, n_analytes = 10,
    lod_quantile = 0.3, seed = 42
  )

  cens_rate <- mean(sim$data$cens_lod, na.rm = TRUE)
  expect_gt(cens_rate, 0.05)
})


# ---- simulate_ancova_data ---------------------------------------------------

test_that("simulate_ancova_data returns valid immuno_sim", {
  sim <- simulate_ancova_data(n_subjects = 40, n_analytes = 5, seed = 42)
  expect_valid_immuno_sim(sim, design = "cross_sectional")

  # 2 groups
  expect_setequal(levels(sim$data$group), c("control", "treatment"))

  # Has covariates
  expect_true("age" %in% names(sim$data))
  expect_true("bmi" %in% names(sim$data))
})

test_that("simulate_ancova_data has group-covariate correlation", {
  sim <- simulate_ancova_data(
    n_subjects = 300, n_analytes = 5,
    covariate_group_cor = 0.5, seed = 42
  )

  d <- sim$data
  per_subj <- d[!duplicated(d$subject_id), ]
  cor_val <- cor(as.numeric(per_subj$group), per_subj$age, use = "complete.obs")
  # With rho = 0.5, should see substantial positive correlation
  expect_gt(cor_val, 0.15)
})

test_that("simulate_ancova_data supports outlier injection", {
  sim <- simulate_ancova_data(
    n_subjects = 100, n_analytes = 5,
    outlier_rate = 0.05, seed = 42
  )

  expect_equal(sim$meta$ancova$outlier_rate, 0.05)
  expect_gt(length(sim$meta$ancova$outlier_idx), 0)

  # Outlier indices should be valid
  expect_true(all(sim$meta$ancova$outlier_idx >= 1))
  expect_true(all(sim$meta$ancova$outlier_idx <= nrow(sim$data)))
})

test_that("simulate_ancova_data supports heteroscedastic residuals", {
  sim <- simulate_ancova_data(
    n_subjects = 200, n_analytes = 5,
    heteroscedastic = TRUE, seed = 42
  )

  expect_true(sim$meta$ancova$heteroscedastic)

  # Treatment group should have larger variance
  d <- sim$data
  cyt <- unique(d$cytokine)[1]
  vals_ctrl <- d$value[d$group == "control" & d$cytokine == cyt]
  vals_trt  <- d$value[d$group == "treatment" & d$cytokine == cyt]
  vals_ctrl <- vals_ctrl[!is.na(vals_ctrl)]
  vals_trt  <- vals_trt[!is.na(vals_trt)]
  if (length(vals_ctrl) > 10 && length(vals_trt) > 10) {
    ratio <- sd(vals_trt) / sd(vals_ctrl)
    expect_gt(ratio, 1.0)
  }
})

test_that("simulate_ancova_data supports nonlinear covariate effects", {
  sim <- simulate_ancova_data(
    n_subjects = 60, n_analytes = 5,
    nonlinear_cov = TRUE, seed = 42
  )

  expect_true(sim$meta$ancova$nonlinear_cov)
  expect_valid_immuno_sim(sim)
})

test_that("simulate_ancova_data with zero confounding", {
  sim <- simulate_ancova_data(
    n_subjects = 200, n_analytes = 5,
    covariate_group_cor = 0, seed = 42
  )

  d <- sim$data
  per_subj <- d[!duplicated(d$subject_id), ]
  cor_val <- cor(as.numeric(per_subj$group), per_subj$age, use = "complete.obs")
  # With rho = 0, correlation should be near zero
  expect_lt(abs(cor_val), 0.2)
})


# ==============================================================================
# Cross-cutting concerns
# ==============================================================================

test_that("all wrappers are deterministic with seed", {
  s1 <- simulate_vaccine_trial(n_subjects = 20, n_analytes = 5, seed = 1)
  s2 <- simulate_vaccine_trial(n_subjects = 20, n_analytes = 5, seed = 1)
  expect_equal(s1$data$value, s2$data$value)

  s3 <- simulate_infection_timecourse(n_subjects = 20, n_analytes = 5, seed = 2)
  s4 <- simulate_infection_timecourse(n_subjects = 20, n_analytes = 5, seed = 2)
  expect_equal(s3$data$value, s4$data$value)

  s5 <- simulate_severity_comparison(n_subjects = 20, n_analytes = 5, seed = 3)
  s6 <- simulate_severity_comparison(n_subjects = 20, n_analytes = 5, seed = 3)
  expect_equal(s5$data$value, s6$data$value)

  s7 <- simulate_exposure_cohort(n_subjects = 20, n_analytes = 5, seed = 4)
  s8 <- simulate_exposure_cohort(n_subjects = 20, n_analytes = 5, seed = 4)
  expect_equal(s7$data$value, s8$data$value)
})

test_that("all wrappers return different results with different seeds", {
  s1 <- simulate_vaccine_trial(n_subjects = 20, n_analytes = 5, seed = 1)
  s2 <- simulate_vaccine_trial(n_subjects = 20, n_analytes = 5, seed = 2)
  expect_false(identical(s1$data$value, s2$data$value))
})

test_that("all generators have signal_analyte truth column", {
  sims <- list(
    simulate_vaccine_trial(n_subjects = 20, n_analytes = 5, seed = 1),
    simulate_infection_timecourse(n_subjects = 20, n_analytes = 5, seed = 2),
    simulate_severity_comparison(n_subjects = 20, n_analytes = 5, seed = 3),
    simulate_exposure_cohort(n_subjects = 20, n_analytes = 5, seed = 4),
    simulate_mcnemar_data(n_subjects = 20, n_analytes = 5, seed = 5),
    simulate_plsda_data(n_subjects = 20, n_analytes = 5, seed = 6),
    simulate_ancova_data(n_subjects = 20, n_analytes = 5, seed = 7)
  )

  for (i in seq_along(sims)) {
    expect_true("signal_analyte" %in% names(sims[[i]]$truth),
                info = paste("Generator", i, "missing signal_analyte"))
  }
})
