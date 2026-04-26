# Tests for apply_assay_process() — Phase 1c

# ---- Helper: generate latent biology for assay tests -------------------------

make_latent <- function(n_subjects = 60, n_analytes = 5, n_timepoints = 1,
                        seed = 42, signal_analytes = 0, ...) {
  immunoPlex:::simulate_latent_biology(
    n_subjects = n_subjects,
    n_analytes = n_analytes,
    n_timepoints = n_timepoints,
    design = "cross_sectional",
    signal_analytes = signal_analytes,
    seed = seed,
    ...
  )
}


# ---- Identity / passthrough test --------------------------------------------

test_that("zero assay effects preserve latent values", {
  lat <- make_latent(n_subjects = 30, n_analytes = 5, seed = 100)
  an <- lat$params$analyte_names

  result <- immunoPlex:::apply_assay_process(
    lat,
    lod_values = stats::setNames(rep(0, 5), an),
    n_plates = 1, plate_sd = 0,
    n_lots = 1, lot_sd = 0,
    bg_shape = 0,
    n_replicates = 1, replicate_sd = 0,
    well_failure_rate = 0,
    covariate_missingness = 0
  )

  # Log-scale values should round-trip through exp/log
  expect_equal(result$latent_data$value, lat$latent_data$value)
  # Raw-scale = exp(latent)
  non_na <- !is.na(lat$latent_data$value)
  expect_equal(result$latent_data$value_raw[non_na],
               exp(lat$latent_data$value[non_na]))
  # No censoring
  expect_true(all(!result$latent_data$cens_lod))
  expect_true(all(!result$latent_data$cens_ulod))
  # All pass QC
  expect_true(all(result$latent_data$qc_flag == "pass"))
  # Truth columns: zero plate/lot effects
  expect_true(all(result$truth_data$plate_effect == 0))
  expect_true(all(result$truth_data$lot_effect == 0))
  # true_raw matches
  expect_equal(result$truth_data$true_raw, exp(lat$truth_data$true_latent))
})


# ---- Output structure -------------------------------------------------------

test_that("output has correct structure and columns", {
  lat <- make_latent(n_subjects = 20, n_analytes = 3, seed = 1)
  result <- immunoPlex:::apply_assay_process(lat)

  expect_type(result, "list")
  expect_named(result, c("latent_data", "truth_data", "params"))
  expect_s3_class(result$latent_data, "tbl_df")
  expect_s3_class(result$truth_data, "tbl_df")

  # Data columns
  expected_cols <- c("subject_id", "timepoint", "cytokine", "group", "value",
                     "missing_reason", "plate_id", "lot_id", "replicate_id",
                     "value_raw", "cens_lod", "cens_ulod", "lod", "ulod",
                     "qc_flag")
  expect_true(all(expected_cols %in% names(result$latent_data)))

  # Truth columns
  truth_cols <- c("true_latent", "signal_analyte", "true_group_effect",
                  "true_time_effect", "true_interaction", "true_raw",
                  "plate_effect", "lot_effect")
  expect_true(all(truth_cols %in% names(result$truth_data)))

  # Rows aligned

  expect_equal(nrow(result$latent_data), nrow(result$truth_data))
})

test_that("params stores assay parameters and realized censoring", {
  lat <- make_latent(n_subjects = 20, n_analytes = 3, seed = 1)
  result <- immunoPlex:::apply_assay_process(lat, n_plates = 2, plate_sd = 0.1)

  expect_true(!is.null(result$params$assay))
  expect_equal(result$params$assay$n_plates, 2)
  expect_equal(result$params$assay$plate_sd, 0.1)
  expect_true(!is.null(result$params$realized_censoring))
  expect_equal(length(result$params$realized_censoring), 3)
})


# ---- Censoring within 5pp of target -----------------------------------------

test_that("censoring rate within 5pp of target via lod_quantile", {
  # Use large n, no assay effects, set LOD at 30th percentile
  lat <- make_latent(n_subjects = 500, n_analytes = 5, seed = 200)
  an <- lat$params$analyte_names

  result <- immunoPlex:::apply_assay_process(
    lat,
    lod_quantile = 0.30,
    n_plates = 1, plate_sd = 0,
    n_lots = 1, lot_sd = 0,
    bg_shape = 0,
    n_replicates = 1, replicate_sd = 0,
    well_failure_rate = 0
  )

  rc <- result$params$realized_censoring
  for (nm in an) {
    expect_lt(abs(rc[nm] - 0.30), 0.05,
              label = paste("censoring rate for", nm))
  }
})

test_that("explicit LOD values produce expected censoring", {
  lat <- make_latent(n_subjects = 200, n_analytes = 3, seed = 300)
  an <- lat$params$analyte_names

  # Set very high LOD for first analyte (should censor nearly all)
  # and very low for others (should censor nearly none)
  lods <- stats::setNames(c(1e6, 0.001, 0.001), an)
  result <- immunoPlex:::apply_assay_process(
    lat, lod_values = lods,
    n_plates = 1, plate_sd = 0, n_lots = 1, lot_sd = 0,
    bg_shape = 0, n_replicates = 1, replicate_sd = 0
  )

  rc <- result$params$realized_censoring
  expect_gt(rc[an[1]], 0.95)  # Nearly all censored
  expect_lt(rc[an[2]], 0.05)  # Nearly none censored
  expect_lt(rc[an[3]], 0.05)
})

test_that("ULOD censoring works", {
  lat <- make_latent(n_subjects = 100, n_analytes = 3, seed = 400)
  an <- lat$params$analyte_names

  # Low ULOD for first analyte
  ulods <- stats::setNames(c(1.0, Inf, Inf), an)
  result <- immunoPlex:::apply_assay_process(
    lat, lod_values = stats::setNames(rep(0, 3), an),
    ulod_values = ulods,
    n_plates = 1, plate_sd = 0, n_lots = 1, lot_sd = 0,
    bg_shape = 0, n_replicates = 1, replicate_sd = 0
  )

  # First analyte should have some ULOD censoring (mean ~ exp(1.5) = 4.5 pg/mL,
  # ULOD = 1.0 so most should be censored)
  ulod_frac <- mean(result$latent_data$cens_ulod[result$latent_data$cytokine == an[1]])
  expect_gt(ulod_frac, 0.3)
  # Others should have no ULOD censoring
  expect_equal(sum(result$latent_data$cens_ulod[result$latent_data$cytokine == an[2]]), 0)
})


# ---- QC flags at expected rate -----------------------------------------------

test_that("well failures appear at expected rate", {
  lat <- make_latent(n_subjects = 200, n_analytes = 5, seed = 500)
  an <- lat$params$analyte_names

  result <- immunoPlex:::apply_assay_process(
    lat,
    lod_values = stats::setNames(rep(0, 5), an),
    well_failure_rate = 0.10,
    n_plates = 1, plate_sd = 0, n_lots = 1, lot_sd = 0,
    bg_shape = 0, n_replicates = 1, replicate_sd = 0
  )

  n_total <- nrow(result$latent_data)
  n_failed <- sum(result$latent_data$qc_flag == "failed_well")
  achieved_rate <- n_failed / n_total
  expect_lt(abs(achieved_rate - 0.10), 0.02,
            label = "well failure rate")
  # Failed wells should have NA value
  expect_true(all(is.na(result$latent_data$value[
    result$latent_data$qc_flag == "failed_well"])))
})

test_that("high-CV flag works with replicates", {
  lat <- make_latent(n_subjects = 100, n_analytes = 3, seed = 600)
  an <- lat$params$analyte_names

  # High replicate noise should produce some high-CV flags
  result <- immunoPlex:::apply_assay_process(
    lat,
    lod_values = stats::setNames(rep(0, 3), an),
    n_replicates = 3, replicate_sd = 0.5,
    cv_threshold = 0.15,
    n_plates = 1, plate_sd = 0, n_lots = 1, lot_sd = 0,
    bg_shape = 0, well_failure_rate = 0
  )

  expect_gt(sum(result$latent_data$qc_flag == "high_cv"), 0,
            label = "some high-CV flags expected")

  # With no replicate noise, no high-CV flags
  result_clean <- immunoPlex:::apply_assay_process(
    lat,
    lod_values = stats::setNames(rep(0, 3), an),
    n_replicates = 3, replicate_sd = 0,
    cv_threshold = 0.15,
    n_plates = 1, plate_sd = 0, n_lots = 1, lot_sd = 0,
    bg_shape = 0, well_failure_rate = 0
  )
  expect_equal(sum(result_clean$latent_data$qc_flag == "high_cv"), 0)
})


# ---- Lot-specific LODs differ ------------------------------------------------

test_that("lot-specific LODs differ when lot_sd > 0", {
  lat <- make_latent(n_subjects = 40, n_analytes = 3, seed = 700)
  an <- lat$params$analyte_names

  result <- immunoPlex:::apply_assay_process(
    lat,
    n_lots = 2, lot_sd = 0.15,
    n_plates = 1, plate_sd = 0,
    bg_shape = 0, n_replicates = 1, replicate_sd = 0,
    well_failure_rate = 0
  )

  d <- result$latent_data
  # For each analyte, LOD should differ between lots
  for (nm in an) {
    sub <- d[d$cytokine == nm, ]
    lod_by_lot <- tapply(sub$lod, sub$lot_id, unique)
    # With 2 lots, should have 2 distinct LOD values
    expect_equal(length(lod_by_lot), 2)
    expect_false(lod_by_lot[[1]] == lod_by_lot[[2]],
                 label = paste("LODs differ for", nm))
  }

  # Lot effects stored in params
  expect_true(!is.null(result$params$assay$lot_effects))
  expect_equal(nrow(result$params$assay$lot_effects), 2)
})

test_that("lot LODs equal base LODs when lot_sd = 0", {
  lat <- make_latent(n_subjects = 20, n_analytes = 3, seed = 710)
  an <- lat$params$analyte_names

  result <- immunoPlex:::apply_assay_process(
    lat,
    n_lots = 2, lot_sd = 0,
    n_plates = 1, plate_sd = 0,
    bg_shape = 0, n_replicates = 1, replicate_sd = 0,
    well_failure_rate = 0
  )

  d <- result$latent_data
  # All LODs should be the same (base LOD) regardless of lot
  for (nm in an) {
    lods <- unique(d$lod[d$cytokine == nm])
    expect_equal(length(lods), 1)
  }
})


# ---- Replicate noise increases CV -------------------------------------------

test_that("replicate noise increases CV", {
  lat <- make_latent(n_subjects = 100, n_analytes = 3, seed = 800)
  an <- lat$params$analyte_names

  # Without noise: all replicates identical, CV = 0
  res_clean <- immunoPlex:::apply_assay_process(
    lat,
    lod_values = stats::setNames(rep(0, 3), an),
    n_replicates = 3, replicate_sd = 0,
    n_plates = 1, plate_sd = 0, n_lots = 1, lot_sd = 0,
    bg_shape = 0, well_failure_rate = 0
  )

  # With noise: replicates differ, CV > 0
  res_noisy <- immunoPlex:::apply_assay_process(
    lat,
    lod_values = stats::setNames(rep(0, 3), an),
    n_replicates = 3, replicate_sd = 0.15,
    n_plates = 1, plate_sd = 0, n_lots = 1, lot_sd = 0,
    bg_shape = 0, well_failure_rate = 0
  )

  # Compute mean CV across observations for each case
  compute_mean_cv <- function(d) {
    orig_n <- nrow(d) / 3
    orig_id <- rep(seq_len(orig_n), each = 3)
    cvs <- tapply(d$value_raw, orig_id, function(x) {
      xv <- x[!is.na(x)]
      if (length(xv) >= 2) stats::sd(xv) / mean(xv) else 0
    })
    mean(cvs, na.rm = TRUE)
  }

  cv_clean <- compute_mean_cv(res_clean$latent_data)
  cv_noisy <- compute_mean_cv(res_noisy$latent_data)

  expect_equal(cv_clean, 0)
  expect_gt(cv_noisy, 0.05)
})


# ---- Replicate expansion ----------------------------------------------------

test_that("replicate expansion produces correct row count", {
  lat <- make_latent(n_subjects = 20, n_analytes = 3, seed = 810)
  an <- lat$params$analyte_names

  result <- immunoPlex:::apply_assay_process(
    lat,
    lod_values = stats::setNames(rep(0, 3), an),
    n_replicates = 4,
    n_plates = 1, plate_sd = 0, n_lots = 1, lot_sd = 0,
    bg_shape = 0, replicate_sd = 0, well_failure_rate = 0
  )

  expect_equal(nrow(result$latent_data), nrow(lat$latent_data) * 4)
  expect_equal(nrow(result$truth_data), nrow(lat$truth_data) * 4)
  expect_true(all(result$latent_data$replicate_id %in% 1:4))
})

test_that("single replicate produces replicate_id = 1", {
  lat <- make_latent(n_subjects = 10, n_analytes = 2, seed = 820)
  result <- immunoPlex:::apply_assay_process(lat)
  expect_true(all(result$latent_data$replicate_id == 1L))
  expect_equal(nrow(result$latent_data), nrow(lat$latent_data))
})


# ---- Plate effects ----------------------------------------------------------

test_that("plate effects introduce between-plate variation", {
  lat <- make_latent(n_subjects = 200, n_analytes = 3, seed = 900,
                     re_intercept_sd = 0)
  an <- lat$params$analyte_names

  result <- immunoPlex:::apply_assay_process(
    lat,
    lod_values = stats::setNames(rep(0, 3), an),
    n_plates = 4, plate_sd = 0.20,
    n_lots = 1, lot_sd = 0,
    bg_shape = 0, n_replicates = 1, replicate_sd = 0,
    well_failure_rate = 0
  )

  d <- result$latent_data
  # Plate means should differ for same analyte
  for (nm in an) {
    sub <- d[d$cytokine == nm & !is.na(d$value), ]
    plate_means <- tapply(sub$value_raw, sub$plate_id, mean)
    # With plate_sd = 0.20, plate means should vary
    expect_gt(stats::sd(plate_means) / mean(plate_means), 0.05,
              label = paste("plate CV for", nm))
  }

  # Plate effects recorded in truth
  expect_true(any(result$truth_data$plate_effect != 0))
  # 4 unique plates assigned
  expect_equal(length(unique(d$plate_id)), 4)
})

test_that("plate_sd = 0 produces zero plate effects", {
  lat <- make_latent(n_subjects = 20, n_analytes = 2, seed = 910)
  result <- immunoPlex:::apply_assay_process(
    lat, n_plates = 3, plate_sd = 0
  )
  expect_true(all(result$truth_data$plate_effect == 0))
})


# ---- Background binding -----------------------------------------------------

test_that("background binding increases raw values", {
  lat <- make_latent(n_subjects = 100, n_analytes = 3, seed = 1000)
  an <- lat$params$analyte_names

  res_no_bg <- immunoPlex:::apply_assay_process(
    lat,
    lod_values = stats::setNames(rep(0, 3), an),
    bg_shape = 0,
    n_plates = 1, plate_sd = 0, n_lots = 1, lot_sd = 0,
    n_replicates = 1, replicate_sd = 0, well_failure_rate = 0
  )

  res_bg <- immunoPlex:::apply_assay_process(
    lat,
    lod_values = stats::setNames(rep(0, 3), an),
    bg_shape = 2, bg_rate = 1,
    n_plates = 1, plate_sd = 0, n_lots = 1, lot_sd = 0,
    n_replicates = 1, replicate_sd = 0, well_failure_rate = 0
  )

  # Background-added values should be >= no-background values
  non_na <- !is.na(res_no_bg$latent_data$value_raw)
  expect_true(all(res_bg$latent_data$value_raw[non_na] >=
                    res_no_bg$latent_data$value_raw[non_na] - 1e-10))
  # Mean should be higher
  expect_gt(mean(res_bg$latent_data$value_raw[non_na]),
            mean(res_no_bg$latent_data$value_raw[non_na]))
})


# ---- Covariate missingness --------------------------------------------------

test_that("covariate missingness introduces NAs", {
  lat <- make_latent(
    n_subjects = 100, n_analytes = 3, seed = 1100,
    covariates = list(age = list(type = "continuous", effect = 0.1))
  )
  an <- lat$params$analyte_names

  # Before assay: no NA in covariates
  expect_true(all(!is.na(lat$latent_data$age)))

  result <- immunoPlex:::apply_assay_process(
    lat,
    lod_values = stats::setNames(rep(0, 3), an),
    covariate_missingness = 0.20,
    n_plates = 1, plate_sd = 0, n_lots = 1, lot_sd = 0,
    bg_shape = 0, n_replicates = 1, replicate_sd = 0,
    well_failure_rate = 0
  )

  frac_missing <- mean(is.na(result$latent_data$age))
  expect_lt(abs(frac_missing - 0.20), 0.05)
})


# ---- Determinism -------------------------------------------------------------

test_that("same seed in latent biology produces deterministic assay process", {
  lat1 <- make_latent(n_subjects = 30, n_analytes = 3, seed = 1200)
  lat2 <- make_latent(n_subjects = 30, n_analytes = 3, seed = 1200)

  # apply_assay_process uses RNG for plate/lot effects, bg, replicate noise, etc.
  # With the same latent input and same RNG state, results should be identical.
  set.seed(99)
  res1 <- immunoPlex:::apply_assay_process(
    lat1, n_plates = 2, plate_sd = 0.1,
    n_lots = 2, lot_sd = 0.05,
    bg_shape = 1, n_replicates = 2, replicate_sd = 0.1,
    well_failure_rate = 0.05
  )
  set.seed(99)
  res2 <- immunoPlex:::apply_assay_process(
    lat2, n_plates = 2, plate_sd = 0.1,
    n_lots = 2, lot_sd = 0.05,
    bg_shape = 1, n_replicates = 2, replicate_sd = 0.1,
    well_failure_rate = 0.05
  )

  expect_identical(res1$latent_data$value, res2$latent_data$value)
  expect_identical(res1$truth_data$plate_effect, res2$truth_data$plate_effect)
  expect_identical(res1$truth_data$lot_effect, res2$truth_data$lot_effect)
})


# ---- Missingness preservation -----------------------------------------------

test_that("pre-existing missingness is preserved", {
  lat <- make_latent(
    n_subjects = 100, n_analytes = 3, n_timepoints = 2,
    seed = 1300, dropout_rate = 0.3
  )
  an <- lat$params$analyte_names

  # Count existing NAs
  n_na_before <- sum(is.na(lat$latent_data$value))
  expect_gt(n_na_before, 0)

  result <- immunoPlex:::apply_assay_process(
    lat,
    lod_values = stats::setNames(rep(0, 3), an),
    n_plates = 1, plate_sd = 0, n_lots = 1, lot_sd = 0,
    bg_shape = 0, n_replicates = 1, replicate_sd = 0,
    well_failure_rate = 0
  )

  # Originally missing values should still be NA
  was_na <- is.na(lat$latent_data$value)
  expect_true(all(is.na(result$latent_data$value[was_na])))
  # Missing reason preserved for dropout rows
  dropout_rows <- lat$latent_data$missing_reason %in% "dropout"
  expect_true(all(result$latent_data$missing_reason[dropout_rows] == "dropout"))
})


# ---- Combined effects -------------------------------------------------------

test_that("full pipeline with all effects produces valid output", {
  lat <- make_latent(
    n_subjects = 80, n_analytes = 5, seed = 1400,
    signal_analytes = 2, group_effects = 0.6,
    covariates = list(age = list(type = "continuous", effect = 0.1))
  )

  result <- immunoPlex:::apply_assay_process(
    lat,
    n_plates = 3, plate_sd = 0.15,
    n_lots = 2, lot_sd = 0.10,
    bg_shape = 2, bg_rate = 1,
    n_replicates = 2, replicate_sd = 0.10,
    well_failure_rate = 0.03,
    cv_threshold = 0.25,
    covariate_missingness = 0.10
  )

  d <- result$latent_data
  tr <- result$truth_data

  # Correct dimensions
  expect_equal(nrow(d), nrow(lat$latent_data) * 2)  # 2 replicates
  expect_equal(nrow(d), nrow(tr))

  # All required columns present
  expect_true("plate_id" %in% names(d))
  expect_true("lot_id" %in% names(d))
  expect_true("replicate_id" %in% names(d))
  expect_true("value_raw" %in% names(d))
  expect_true("cens_lod" %in% names(d))
  expect_true("qc_flag" %in% names(d))
  expect_true("plate_effect" %in% names(tr))
  expect_true("lot_effect" %in% names(tr))
  expect_true("true_raw" %in% names(tr))

  # QC flags are valid values
  expect_true(all(d$qc_flag %in% c("pass", "failed_well", "high_cv")))

  # Censoring flags are logical
  expect_type(d$cens_lod, "logical")
  expect_type(d$cens_ulod, "logical")

  # Censored values are clamped to LOD
  lod_rows <- d$cens_lod & !is.na(d$value_raw)
  if (sum(lod_rows) > 0) {
    expect_equal(d$value_raw[lod_rows], d$lod[lod_rows])
  }

  # Some well failures occurred
  expect_gt(sum(d$qc_flag == "failed_well"), 0)

  # Some covariate missingness
  expect_gt(sum(is.na(d$age)), 0)

  # Realized censoring is stored
  expect_equal(length(result$params$realized_censoring), 5)
})
