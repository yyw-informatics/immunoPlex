# Tests for simulate_latent_biology() — Phase 1b

# ---- Determinism -------------------------------------------------------------

test_that("same seed produces identical output", {
  a <- immunoPlex:::simulate_latent_biology(n_subjects = 20, n_analytes = 5,
                                            seed = 123)
  b <- immunoPlex:::simulate_latent_biology(n_subjects = 20, n_analytes = 5,
                                            seed = 123)
  expect_identical(a$latent_data, b$latent_data)
  expect_identical(a$truth_data, b$truth_data)
  # params contain realized summaries computed from the same RNG stream
  expect_equal(a$params$realized_icc, b$params$realized_icc)
})

test_that("different seeds produce different output", {
  a <- immunoPlex:::simulate_latent_biology(n_subjects = 20, n_analytes = 5,
                                            seed = 1)
  b <- immunoPlex:::simulate_latent_biology(n_subjects = 20, n_analytes = 5,
                                            seed = 2)
  expect_false(identical(a$latent_data$value, b$latent_data$value))
})


# ---- Output structure --------------------------------------------------------

test_that("return value has correct top-level structure", {
  res <- immunoPlex:::simulate_latent_biology(n_subjects = 10, n_analytes = 3,
                                              seed = 1)
  expect_type(res, "list")
  expect_named(res, c("latent_data", "truth_data", "params"))
  expect_s3_class(res$latent_data, "tbl_df")
  expect_s3_class(res$truth_data, "tbl_df")
  expect_type(res$params, "list")
})

test_that("latent_data has required columns", {
  res <- immunoPlex:::simulate_latent_biology(n_subjects = 10, n_analytes = 3,
                                              seed = 1)
  required <- c("subject_id", "timepoint", "cytokine", "group", "value",
                "missing_reason")
  expect_true(all(required %in% names(res$latent_data)))
})

test_that("truth_data has required columns", {
  res <- immunoPlex:::simulate_latent_biology(n_subjects = 10, n_analytes = 3,
                                              seed = 1)
  required <- c("true_latent", "signal_analyte", "true_group_effect",
                "true_time_effect", "true_interaction")
  expect_true(all(required %in% names(res$truth_data)))
})


# ---- Dimensions --------------------------------------------------------------

test_that("dimensions match n_subjects x n_timepoints x n_analytes", {
  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 30, n_timepoints = 3, n_analytes = 10, seed = 42
  )
  expected_rows <- 30 * 3 * 10
  expect_equal(nrow(res$latent_data), expected_rows)
  expect_equal(nrow(res$truth_data), expected_rows)
  expect_equal(length(unique(res$latent_data$subject_id)), 30)
  expect_equal(length(unique(res$latent_data$timepoint)), 3)
  expect_equal(length(unique(res$latent_data$cytokine)), 10)
})

test_that("cross-sectional design (1 timepoint) works", {
  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 20, n_timepoints = 1, n_analytes = 5,
    design = "cross_sectional", seed = 1
  )
  expect_equal(nrow(res$latent_data), 20 * 1 * 5)
  expect_equal(length(unique(res$latent_data$timepoint)), 1)
})

test_that("single analyte works", {
  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 10, n_analytes = 1, signal_analytes = 1, seed = 1
  )
  expect_equal(nrow(res$latent_data), 10 * 2 * 1)
  expect_equal(length(unique(res$latent_data$cytokine)), 1)
})

test_that("more analytes than library recycles correctly", {
  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 5, n_analytes = 25, signal_analytes = 2, seed = 1
  )
  expect_equal(length(unique(res$latent_data$cytokine)), 25)
  expect_true(all(grepl("^Analyte_", unique(res$latent_data$cytokine))))
})


# ---- Group allocation -------------------------------------------------------

test_that("balanced groups by default", {
  res <- immunoPlex:::simulate_latent_biology(n_subjects = 60, n_analytes = 2,
                                              seed = 1)
  grp_counts <- table(res$latent_data$group[
    res$latent_data$timepoint == "T1" & res$latent_data$cytokine == "IL-1b"
  ])
  expect_equal(as.integer(grp_counts), c(30, 30))
})

test_that("unbalanced allocation works", {
  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 60, n_analytes = 2,
    group_allocation = c(1, 2), seed = 1
  )
  grp_counts <- table(res$latent_data$group[
    res$latent_data$timepoint == "T1" & res$latent_data$cytokine == "IL-1b"
  ])
  expect_equal(as.integer(grp_counts), c(20, 40))
})

test_that("three groups work", {
  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 30, n_analytes = 2,
    group_levels = c("A", "B", "C"), seed = 1
  )
  grp_counts <- table(res$latent_data$group[
    res$latent_data$timepoint == "T1" & res$latent_data$cytokine == "IL-1b"
  ])
  expect_equal(length(grp_counts), 3)
  expect_equal(sum(grp_counts), 30)
})


# ---- Signal analytes ---------------------------------------------------------

test_that("signal_analytes = 0 produces no group effects", {
  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 20, n_analytes = 5, signal_analytes = 0, seed = 1
  )
  expect_true(all(res$truth_data$true_group_effect == 0))
  expect_true(all(!res$truth_data$signal_analyte))
})

test_that("signal_analytes as character vector selects correct analytes", {
  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 20, n_analytes = 5,
    signal_analytes = c("IL-1b", "IL-6"), seed = 1
  )
  sig_names <- unique(res$latent_data$cytokine[res$truth_data$signal_analyte])
  expect_setequal(sig_names, c("IL-1b", "IL-6"))
})

test_that("signal_analytes as logical vector works", {
  mask <- c(TRUE, FALSE, TRUE, FALSE, FALSE)
  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 20, n_analytes = 5, signal_analytes = mask, seed = 1
  )
  expect_equal(res$params$signal_analytes, mask)
})


# ---- Large-n group effect recovery ------------------------------------------

test_that("large-n run recovers true group effects within 2 SE", {
  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 1000,
    n_timepoints = 1,
    n_analytes = 5,
    group_effects = 0.8,
    signal_analytes = c("IL-1b", "IL-2"),
    effect_direction = "up",
    re_intercept_sd = 0.3,
    residual_sd = 0.8,
    design = "cross_sectional",
    seed = 202
  )

  d <- res$latent_data
  tr <- res$truth_data
  analytes <- unique(d$cytokine)

  for (cyt in analytes) {
    rows <- d$cytokine == cyt
    sub <- d[rows, ]
    true_eff <- unique(tr$true_group_effect[rows & d$group == "treatment"])

    if (true_eff[1] != 0) {
      fit <- stats::lm(value ~ group, data = sub)
      est <- stats::coef(fit)["grouptreatment"]
      se  <- summary(fit)$coefficients["grouptreatment", "Std. Error"]
      expect_lt(abs(est - true_eff[1]), 2 * se,
                label = paste(cyt, ": |est - true| < 2 SE"))
    } else {
      # Non-signal: effect should be near zero
      fit <- stats::lm(value ~ group, data = sub)
      est <- stats::coef(fit)["grouptreatment"]
      se  <- summary(fit)$coefficients["grouptreatment", "Std. Error"]
      expect_lt(abs(est), 2 * se,
                label = paste(cyt, ": non-signal effect near 0"))
    }
  }
})


# ---- ICC recovery ------------------------------------------------------------

test_that("achieved ICC is near target", {
  re_sd   <- 0.6
  res_sd  <- 0.8
  target_icc <- re_sd^2 / (re_sd^2 + res_sd^2)

  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 200,
    n_timepoints = 3,
    n_analytes = 5,
    re_intercept_sd = re_sd,
    residual_sd = res_sd,
    signal_analytes = 0,
    seed = 42
  )

  for (k in seq_len(5)) {
    expect_lt(abs(res$params$realized_icc[k] - target_icc), 0.1,
              label = paste("ICC for analyte", k))
  }
})


# ---- Analyte correlation recovery -------------------------------------------

test_that("block correlation is achieved near target", {
  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 300,
    n_timepoints = 1,
    n_analytes = 10,
    analyte_correlation = "block",
    block_rho = 0.5,
    signal_analytes = 0,
    re_intercept_sd = 0.5,
    residual_sd = 1.0,
    design = "cross_sectional",
    seed = 77
  )

  # Check that within-block pairs in the realized correlation matrix are near
  # target. The realized_correlation is based on random intercepts only, so
  # check the cor_matrix stored in params.
  lib <- immunoPlex:::analyte_library
  cats <- lib$category[1:10]
  cor_realized <- res$params$realized_correlation

  for (cat in unique(cats)) {
    idx <- which(cats == cat)
    if (length(idx) >= 2) {
      pairs <- cor_realized[idx, idx]
      off_diag <- pairs[upper.tri(pairs)]
      # Within-block correlations should be near block_rho
      expect_true(all(abs(off_diag - 0.5) < 0.15),
                  label = paste("Block correlation for", cat))
    }
  }
})

test_that("toeplitz correlation produces decaying pattern", {
  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 300,
    n_analytes = 5,
    n_timepoints = 1,
    analyte_correlation = "toeplitz",
    block_rho = 0.6,
    signal_analytes = 0,
    re_intercept_sd = 0.5,
    residual_sd = 1.0,
    design = "cross_sectional",
    seed = 55
  )
  cor_realized <- res$params$realized_correlation
  # Lag-1 should be greater than lag-2
  expect_gt(abs(cor_realized[1, 2]), abs(cor_realized[1, 4]))
})

test_that("custom correlation matrix is respected", {
  custom <- diag(3)
  custom[1, 2] <- custom[2, 1] <- 0.7
  custom[1, 3] <- custom[3, 1] <- 0.3
  custom[2, 3] <- custom[3, 2] <- 0.2

  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 500,
    n_analytes = 3,
    n_timepoints = 1,
    analyte_correlation = custom,
    signal_analytes = 0,
    re_intercept_sd = 0.5,
    residual_sd = 1.0,
    design = "cross_sectional",
    seed = 33
  )
  cor_r <- res$params$realized_correlation
  expect_lt(abs(cor_r[1, 2] - 0.7), 0.15)
  expect_lt(abs(cor_r[1, 3] - 0.3), 0.15)
})


# ---- Dropout -----------------------------------------------------------------

test_that("dropout rate is near target", {
  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 200,
    n_timepoints = 2,
    n_analytes = 5,
    dropout_rate = 0.3,
    signal_analytes = 0,
    seed = 10
  )

  # Fraction of subjects with T2 missing due to dropout
  d <- res$latent_data
  t2 <- d[d$timepoint == "T2", ]
  # Each subject has n_analytes rows at T2; check subject-level dropout
  subj_t2 <- split(t2$missing_reason, t2$subject_id)
  subj_dropped <- vapply(subj_t2, function(x) all(x %in% "dropout"), logical(1))
  achieved_rate <- mean(subj_dropped)

  # Should be within 10pp of target for n=200
  expect_lt(abs(achieved_rate - 0.3), 0.10,
            label = "subject-level dropout rate")
})

test_that("dropout is monotone", {
  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 100,
    n_timepoints = 4,
    n_analytes = 3,
    dropout_rate = 0.2,
    signal_analytes = 0,
    seed = 20
  )
  d <- res$latent_data
  # For every subject, if they are missing at timepoint t, they must be missing

  # at all subsequent timepoints
  for (sid in unique(d$subject_id)) {
    subj <- d[d$subject_id == sid & d$cytokine == unique(d$cytokine)[1], ]
    subj <- subj[order(subj$timepoint), ]
    vals <- !is.na(subj$value)
    # Once FALSE, must stay FALSE
    if (any(!vals)) {
      first_miss <- which(!vals)[1]
      expect_true(all(!vals[first_miss:length(vals)]),
                  label = paste("monotone dropout for", sid))
    }
  }
})

test_that("visit missingness is non-monotone", {
  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 200,
    n_timepoints = 4,
    n_analytes = 3,
    visit_missingness = 0.15,
    dropout_rate = 0,
    signal_analytes = 0,
    seed = 30
  )
  d <- res$latent_data
  n_missed <- sum(d$missing_reason %in% "missed_visit")
  expect_gt(n_missed, 0, label = "some visits missed")

  # Check non-monotone: at least one subject has a gap followed by observation
  found_gap <- FALSE
  for (sid in unique(d$subject_id)) {
    subj <- d[d$subject_id == sid & d$cytokine == unique(d$cytokine)[1], ]
    subj <- subj[order(subj$timepoint), ]
    vals <- !is.na(subj$value)
    # Gap pattern: FALSE followed by TRUE
    if (any(!vals[-length(vals)] & vals[-1])) {
      found_gap <- TRUE
      break
    }
  }
  expect_true(found_gap, label = "non-monotone pattern exists")
})


# ---- Covariates --------------------------------------------------------------

test_that("covariate columns appear in latent_data", {
  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 20, n_analytes = 3,
    covariates = list(age = list(type = "continuous", effect = 0.3),
                      sex = list(type = "binary", effect = -0.2)),
    seed = 1
  )
  expect_true("age" %in% names(res$latent_data))
  expect_true("sex" %in% names(res$latent_data))
  expect_true(all(res$latent_data$sex %in% c(0L, 1L)))
})

test_that("covariate-group correlation works", {
  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 500, n_analytes = 1,
    n_timepoints = 1,
    covariates = list(confounder = list(type = "continuous", effect = 0.5)),
    covariate_group_cor = 0.6,
    signal_analytes = 0,
    design = "cross_sectional",
    seed = 44
  )
  d <- res$latent_data
  # Correlation between covariate and group should be approximately 0.6
  group_num <- as.integer(d$group) - 1
  r <- stats::cor(d$confounder, group_num)
  expect_lt(abs(r - 0.6), 0.15)
})


# ---- Time / interaction effects ---------------------------------------------

test_that("time effects produce expected pattern", {
  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 300,
    n_timepoints = 3,
    n_analytes = 3,
    signal_analytes = c("IL-1b"),
    time_effects = 0.5,
    group_effects = 0,
    effect_direction = "up",
    re_intercept_sd = 0.2,
    residual_sd = 0.5,
    seed = 55
  )
  d <- res$latent_data
  sig <- d[d$cytokine == "IL-1b", ]
  means <- tapply(sig$value, sig$timepoint, mean)
  # T2 should be higher than T1, T3 higher than T2
  expect_gt(means["T2"], means["T1"])
  expect_gt(means["T3"], means["T2"])
})

test_that("interaction effects modify group difference over time", {
  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 400,
    n_timepoints = 2,
    n_analytes = 3,
    signal_analytes = c("IL-1b"),
    group_effects = 0.5,
    interaction_effects = 0.8,
    effect_direction = "up",
    re_intercept_sd = 0.2,
    residual_sd = 0.5,
    seed = 66
  )
  d <- res$latent_data
  sig <- d[d$cytokine == "IL-1b", ]
  # Group difference at T1: ~0.5 (group effect only)
  # Group difference at T2: ~0.5 + 0.8 = 1.3 (group + interaction)
  t1 <- sig[sig$timepoint == "T1", ]
  t2 <- sig[sig$timepoint == "T2", ]
  diff_t1 <- mean(t1$value[t1$group == "treatment"]) -
    mean(t1$value[t1$group == "control"])
  diff_t2 <- mean(t2$value[t2$group == "treatment"]) -
    mean(t2$value[t2$group == "control"])
  expect_gt(diff_t2, diff_t1 + 0.3)
})


# ---- AR(1) residuals --------------------------------------------------------

test_that("AR(1) produces positively autocorrelated residuals", {
  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 200,
    n_timepoints = 10,
    n_analytes = 2,
    residual_ar1 = 0.6,
    signal_analytes = 0,
    re_intercept_sd = 0.01,
    residual_sd = 1.0,
    seed = 77
  )
  d <- res$latent_data
  tr <- res$truth_data
  # Pool lag-1 residual pairs across all subjects for one analyte
  cyt <- unique(d$cytokine)[1]
  lag1_x <- c()
  lag1_y <- c()
  for (sid in unique(d$subject_id)) {
    sub <- d[d$subject_id == sid & d$cytokine == cyt, ]
    sub <- sub[order(sub$timepoint), ]
    v <- sub$value - mean(sub$value)  # demean per-subject
    if (length(v) >= 2 && all(!is.na(v))) {
      lag1_x <- c(lag1_x, v[-length(v)])
      lag1_y <- c(lag1_y, v[-1])
    }
  }
  pooled_ar <- stats::cor(lag1_x, lag1_y)
  # Should be substantially positive (true AR = 0.6)
  expect_gt(pooled_ar, 0.3)
})


# ---- Heteroscedastic --------------------------------------------------------

test_that("heteroscedastic = TRUE inflates variance for treatment group", {
  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 400,
    n_timepoints = 1,
    n_analytes = 3,
    heteroscedastic = TRUE,
    signal_analytes = 0,
    re_intercept_sd = 0,
    residual_sd = 1.0,
    design = "cross_sectional",
    seed = 88
  )
  d <- res$latent_data
  for (cyt in unique(d$cytokine)) {
    sub <- d[d$cytokine == cyt, ]
    var_ctrl <- stats::var(sub$value[sub$group == "control"])
    var_trt  <- stats::var(sub$value[sub$group == "treatment"])
    # Treatment group should have ~1.5^2 = 2.25x the variance
    expect_gt(var_trt / var_ctrl, 1.5)
  }
})


# ---- Random slopes -----------------------------------------------------------

test_that("random slopes add subject-level time variation", {
  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 100,
    n_timepoints = 3,
    n_analytes = 2,
    re_slope_sd = 0.5,
    signal_analytes = 0,
    re_intercept_sd = 0.3,
    residual_sd = 0.5,
    seed = 99
  )
  d <- res$latent_data
  # With random slopes, the variance at later timepoints should be higher
  # (because slope variance contributes (t-1)^2 * re_slope_sd^2)
  cyt1 <- d[d$cytokine == unique(d$cytokine)[1], ]
  var_t1 <- stats::var(cyt1$value[cyt1$timepoint == "T1"])
  var_t3 <- stats::var(cyt1$value[cyt1$timepoint == "T3"])
  expect_gt(var_t3, var_t1 * 1.2)
})


# ---- Params storage ---------------------------------------------------------

test_that("params stores all key parameters", {
  res <- immunoPlex:::simulate_latent_biology(
    n_subjects = 10, n_analytes = 3, seed = 1
  )
  p <- res$params
  expect_equal(p$design, "pre_post")
  expect_equal(p$n_subjects, 10)
  expect_equal(p$n_analytes, 3)
  expect_equal(p$seed, 1)
  expect_true(!is.null(p$realized_icc))
  expect_true(!is.null(p$realized_correlation))
  expect_true(is.matrix(p$cor_matrix))
})
