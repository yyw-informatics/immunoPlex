# Integration tests for simulate_immunoassay() — Phase 1d
#
# Tests that simulate_immunoassay() returns a valid immuno_sim object and
# that $data can be passed directly to fit_one(), ancova_one(),
# mcnemar_detection(), and plsda_preprocess() without error.

# ---- Unit tests for simulate_immunoassay() -----------------------------------

test_that("simulate_immunoassay returns immuno_sim object", {
  sim <- simulate_immunoassay(
    n_subjects = 30, n_timepoints = 2, n_analytes = 5,
    design = "pre_post", signal_analytes = 2, seed = 101
  )
  expect_s3_class(sim, "immuno_sim")
  expect_named(sim, c("data", "truth", "meta"))
  expect_s3_class(sim$data, "tbl_df")
  expect_s3_class(sim$truth, "tbl_df")
  expect_equal(nrow(sim$data), nrow(sim$truth))
})

test_that("simulate_immunoassay $data has required columns", {
  sim <- simulate_immunoassay(
    n_subjects = 20, n_timepoints = 1, n_analytes = 3,
    design = "cross_sectional", seed = 102
  )
  expected_cols <- c("subject_id", "timepoint", "cytokine", "group",
                     "value", "value_raw", "cens_lod", "cens_ulod",
                     "lod", "ulod", "plate_id", "lot_id",
                     "replicate_id", "qc_flag", "missing_reason")
  for (col in expected_cols) {
    expect_true(col %in% names(sim$data), info = paste("Missing column:", col))
  }
})

test_that("simulate_immunoassay $truth has required columns", {
  sim <- simulate_immunoassay(
    n_subjects = 20, n_timepoints = 1, n_analytes = 3,
    design = "cross_sectional", seed = 103
  )
  expected_truth <- c("true_latent", "true_raw", "signal_analyte",
                       "true_group_effect", "true_time_effect",
                       "true_interaction", "plate_effect", "lot_effect")
  for (col in expected_truth) {
    expect_true(col %in% names(sim$truth), info = paste("Missing column:", col))
  }
})

test_that("simulate_immunoassay $meta has required fields", {
  sim <- simulate_immunoassay(
    n_subjects = 20, n_timepoints = 1, n_analytes = 5,
    design = "cross_sectional", seed = 104
  )
  expect_equal(sim$meta$design, "cross_sectional")
  expect_equal(sim$meta$n_subjects, 20L)
  expect_equal(sim$meta$n_analytes, 5L)
  expect_true(!is.null(sim$meta$realized_censoring))
  expect_true(!is.null(sim$meta$realized_icc))
  expect_true(!is.null(sim$meta$realized_correlation))
  expect_true(!is.null(sim$meta$dgp_params))
  expect_true(!is.null(sim$meta$timestamp))
})

test_that("simulate_immunoassay column types are correct", {
  sim <- simulate_immunoassay(
    n_subjects = 20, n_timepoints = 2, n_analytes = 3,
    design = "pre_post", seed = 105
  )
  d <- sim$data
  expect_true(is.factor(d$subject_id))
  expect_true(is.factor(d$timepoint))
  expect_true(is.factor(d$group))
  expect_true(is.character(d$cytokine))
  expect_true(is.numeric(d$value))
  expect_true(is.numeric(d$value_raw))
  expect_true(is.logical(d$cens_lod))
  expect_true(is.logical(d$cens_ulod))
  expect_true(is.numeric(d$lod))
  expect_true(is.numeric(d$ulod))
  expect_true(is.factor(d$plate_id))
  expect_true(is.factor(d$lot_id))
  expect_true(is.integer(d$replicate_id))
  expect_true(is.character(d$qc_flag))
})

test_that("simulate_immunoassay is deterministic with seed", {
  sim1 <- simulate_immunoassay(
    n_subjects = 20, n_analytes = 3, n_timepoints = 1,
    design = "cross_sectional", seed = 999
  )
  sim2 <- simulate_immunoassay(
    n_subjects = 20, n_analytes = 3, n_timepoints = 1,
    design = "cross_sectional", seed = 999
  )
  expect_equal(sim1$data$value, sim2$data$value)
  expect_equal(sim1$truth$true_latent, sim2$truth$true_latent)
})

test_that("simulate_immunoassay dimensions match parameters", {
  sim <- simulate_immunoassay(
    n_subjects = 40, n_timepoints = 3, n_analytes = 8,
    design = "time_course", seed = 106
  )
  d <- sim$data
  expect_equal(length(unique(d$subject_id)), 40)
  expect_equal(length(unique(d$cytokine)), 8)
  expect_equal(length(unique(d$timepoint)), 3)
  expect_equal(nrow(d), 40 * 3 * 8)
})

test_that("simulate_immunoassay with assay effects works", {
  sim <- simulate_immunoassay(
    n_subjects = 30, n_timepoints = 1, n_analytes = 5,
    design = "cross_sectional", signal_analytes = 2,
    n_plates = 2, plate_sd = 0.1,
    n_lots = 2, lot_sd = 0.05,
    n_replicates = 2, replicate_sd = 0.05,
    seed = 107
  )
  expect_s3_class(sim, "immuno_sim")
  # With replicates, each obs is duplicated
  expect_equal(nrow(sim$data), 30 * 1 * 5 * 2)
  expect_true(any(sim$data$replicate_id == 2))
})

test_that("simulate_immunoassay with covariates works", {
  sim <- simulate_immunoassay(
    n_subjects = 30, n_timepoints = 1, n_analytes = 3,
    design = "cross_sectional", signal_analytes = 1,
    covariates = list(
      age = list(type = "continuous", effect = 0.3),
      sex = list(type = "binary", effect = 0.2)
    ),
    seed = 108
  )
  expect_true("age" %in% names(sim$data))
  expect_true("sex" %in% names(sim$data))
})

test_that("simulate_immunoassay print and summary work", {
  sim <- simulate_immunoassay(
    n_subjects = 20, n_timepoints = 1, n_analytes = 5,
    design = "cross_sectional", seed = 109
  )
  expect_output(print(sim), "immuno_sim:")
  expect_output(summary(sim), "Simulated Immunoassay Summary")
})


# ---- Integration: fit_one() -------------------------------------------------

test_that("simulated data passes to fit_one() without error", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")

  sim <- simulate_immunoassay(
    n_subjects = 60, n_timepoints = 2, n_analytes = 5,
    design = "pre_post", signal_analytes = 2, seed = 201
  )

  # Subset to one cytokine
  first_cyt <- unique(sim$data$cytokine)[1]
  d_one <- sim$data[sim$data$cytokine == first_cyt, ]
  d_one <- as.data.frame(d_one)

  # Remove rows with NA value (dropout/missingness)
  d_one <- d_one[!is.na(d_one$value), ]

  # fit_one with gaussian family (log-scale data, no censoring model)
  result <- fit_one(
    d_one,
    family = "gaussian",
    fixed = "timepoint * group",
    random = "(1|subject_id)"
  )
  expect_s3_class(result, "immuno_fit")
  expect_true(is.logical(result$converged))

  # fit_one with tobit family (censoring-aware)
  result_tobit <- fit_one(
    d_one,
    family = "tobit",
    fixed = "timepoint * group",
    random = ""
  )
  expect_s3_class(result_tobit, "immuno_fit")
})


# ---- Integration: ancova_one() ----------------------------------------------

test_that("simulated data passes to ancova_one() without error", {
  skip_if_not_installed("effectsize")

  # Cross-sectional design with a covariate
  sim <- simulate_immunoassay(
    n_subjects = 40, n_timepoints = 1, n_analytes = 5,
    design = "cross_sectional", signal_analytes = 2,
    covariates = list(age = list(type = "continuous", effect = 0.3)),
    seed = 202
  )

  # Subset to one cytokine — gives one row per subject
  first_cyt <- unique(sim$data$cytokine)[1]
  d_one <- sim$data[sim$data$cytokine == first_cyt, ]
  d_one <- as.data.frame(d_one)
  d_one <- d_one[!is.na(d_one$value), ]

  result <- ancova_one(
    data = d_one,
    outcome = "value",
    group = "group",
    covariate = "age",
    outlier_removal = FALSE
  )
  expect_s3_class(result, "immuno_ancova")
  expect_true(is.numeric(result$omega_sq_partial))
})


# ---- Integration: mcnemar_detection() ---------------------------------------

test_that("simulated data passes to mcnemar_detection() without error", {
  skip_if_not_installed("exact2x2")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")

  # Pre-post design: need 2 timepoints for paired analysis
  sim <- simulate_immunoassay(
    n_subjects = 40, n_timepoints = 2, n_analytes = 5,
    design = "pre_post", signal_analytes = 2, seed = 203
  )

  d <- as.data.frame(sim$data)
  d <- d[!is.na(d$value), ]

  tp_levels <- levels(d$timepoint)

  result <- mcnemar_detection(
    data = d,
    subject_col = "subject_id",
    cytokine_col = "cytokine",
    timepoint_col = "timepoint",
    censoring_col = "cens_lod",
    baseline = tp_levels[1],
    comparison = tp_levels[2],
    quiet = TRUE
  )
  expect_s3_class(result, "data.frame")
  expect_true("p_mcnemar" %in% names(result))
  expect_true("cytokine" %in% names(result))
  expect_equal(nrow(result), length(unique(d$cytokine)))
})


# ---- Integration: plsda_preprocess() ----------------------------------------

test_that("simulated data passes to plsda_preprocess() without error", {
  # Cross-sectional to get one row per subject per cytokine
  sim <- simulate_immunoassay(
    n_subjects = 30, n_timepoints = 1, n_analytes = 5,
    design = "cross_sectional", signal_analytes = 2, seed = 204
  )

  d <- sim$data
  d <- d[!is.na(d$value_raw), ]

  # Pivot to wide format: one row per subject, cytokines as columns
  analyte_names <- unique(d$cytokine)
  wide <- tidyr::pivot_wider(
    d,
    id_cols = c("subject_id", "group"),
    names_from = "cytokine",
    values_from = "value_raw"
  )
  wide <- as.data.frame(wide)

  # Build lod_lookup from simulated LOD values
  lod_lookup <- data.frame(
    cytokine = analyte_names,
    lod = vapply(analyte_names, function(nm) {
      vals <- d$lod[d$cytokine == nm]
      vals[!is.na(vals)][1]
    }, numeric(1)),
    stringsAsFactors = FALSE
  )

  result <- plsda_preprocess(
    data = wide,
    cytokine_cols = analyte_names,
    metadata_cols = c("subject_id", "group"),
    lod_lookup = lod_lookup,
    lod_method = "half",
    log_transform = TRUE,
    scale_data = TRUE,
    verbose = FALSE
  )
  expect_s3_class(result, "plsda_preprocessed")
  expect_equal(ncol(result$expression), length(analyte_names))
  expect_true(nrow(result$expression) > 0)
})
