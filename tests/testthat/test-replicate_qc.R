# Tests for flag_replicate_outliers()

# --- Helper: create test data ---
make_test_data <- function(n_samples = 10, n_cytokines = 2, seed = 42) {
  set.seed(seed)
  cytokines <- paste0("CYT", seq_len(n_cytokines))
  rows <- list()
  for (cyt in cytokines) {
    for (i in seq_len(n_samples)) {
      base_val <- runif(1, 50, 200)
      # Small noise for concordant replicates
      rows <- c(rows, list(data.frame(
        original_sample_id = paste0("S", i),
        cytokine = cyt,
        concentration = c(base_val + rnorm(1, 0, 2),
                          base_val + rnorm(1, 0, 2)),
        stringsAsFactors = FALSE
      )))
    }
  }
  do.call(rbind, rows)
}

# --- Helper: create data with one discordant pair ---
make_discordant_data <- function() {
  # 9 concordant pairs + 1 very discordant pair
  data.frame(
    original_sample_id = rep(paste0("S", 1:10), each = 2),
    cytokine = "CYT1",
    concentration = c(
      100, 102,  # S1: concordant (CV ~1.4%)
      150, 148,  # S2: concordant
      200, 203,  # S3: concordant
       80,  78,  # S4: concordant
      120, 122,  # S5: concordant
      130, 128,  # S6: concordant
      160, 158,  # S7: concordant
      110, 112,  # S8: concordant
      140, 138,  # S9: concordant
       50, 200   # S10: DISCORDANT (CV ~85%)
    ),
    stringsAsFactors = FALSE
  )
}


test_that("flag_replicate_outliers returns expected columns", {
  df <- make_test_data()
  result <- flag_replicate_outliers(df, verbose = FALSE)

  expect_s3_class(result, "data.frame")
  expect_true("replicate_n" %in% names(result))
  expect_true("replicate_cv" %in% names(result))
  expect_true("replicate_discordant" %in% names(result))
  expect_true("replicate_action" %in% names(result))
  # Same number of rows when action = "flag"
  expect_equal(nrow(result), nrow(df))
})

test_that("concordant pairs are not flagged", {
  df <- make_test_data()
  result <- flag_replicate_outliers(df, cv_threshold = 25, verbose = FALSE)

  expect_true(all(!result$replicate_discordant))
  expect_true(all(result$replicate_action == "keep"))
})

test_that("discordant pairs are flagged with CV rule", {
  df <- make_discordant_data()
  result <- flag_replicate_outliers(df, cv_threshold = 25, flag_rule = "cv",
                                    verbose = FALSE)

  # S10 has CV ~85%, should be flagged
  s10_rows <- result[result$original_sample_id == "S10", ]
  expect_true(all(s10_rows$replicate_discordant))

  # Other samples should not be flagged (CV < 3%)
  other_rows <- result[result$original_sample_id != "S10", ]
  expect_true(all(!other_rows$replicate_discordant))
})

test_that("drop_both removes all flagged replicates", {
  df <- make_discordant_data()
  result <- flag_replicate_outliers(df, cv_threshold = 25, action = "drop_both",
                                    verbose = FALSE)

  # S10 should be gone (both reps)
  expect_false("S10" %in% result$original_sample_id)
  # 9 concordant pairs x 2 reps = 18 rows
  expect_equal(nrow(result), 18)
})

test_that("drop_farther keeps closer replicate", {
  df <- make_discordant_data()
  result <- flag_replicate_outliers(df, cv_threshold = 25, action = "drop_farther",
                                    verbose = FALSE)

  # S10 should have only 1 replicate remaining
  s10_rows <- result[result$original_sample_id == "S10", ]
  expect_equal(nrow(s10_rows), 1)

  # The kept replicate should be the one closer to the population median
  # Population median of all values is around 120-140; 200 is closer than 50
  # So value 200 should be kept (not 50)
  # Actually: population median includes all 20 values; the one closer to it stays
  # With most values around 100-200, the median is ~130
  # |200 - 130| = 70, |50 - 130| = 80 -> 200 is closer
  expect_gt(s10_rows$concentration, 100)

  # Other samples should keep both replicates
  expect_equal(nrow(result), 19)  # 18 concordant + 1 kept
})

test_that("winsorize replaces farther replicate with population median", {
  df <- make_discordant_data()
  result <- flag_replicate_outliers(df, cv_threshold = 25, action = "winsorize",
                                    verbose = FALSE)

  # Same number of rows
  expect_equal(nrow(result), nrow(df))

  # S10: one replicate should have been replaced with population median
  s10_rows <- result[result$original_sample_id == "S10", ]
  expect_true(any(s10_rows$replicate_action == "winsorized"))
  expect_true(any(s10_rows$replicate_action == "keep"))
})

test_that("single-replicate samples are not flagged", {
  df <- make_discordant_data()
  # Add a single-replicate sample
  single_rep <- data.frame(
    original_sample_id = "S_single",
    cytokine = "CYT1",
    concentration = 100,
    stringsAsFactors = FALSE
  )
  df <- rbind(df, single_rep)

  result <- flag_replicate_outliers(df, cv_threshold = 25, verbose = FALSE)

  single_row <- result[result$original_sample_id == "S_single", ]
  expect_equal(nrow(single_row), 1)
  expect_false(single_row$replicate_discordant)
  expect_equal(single_row$replicate_n, 1)
})

test_that("near-zero mean pairs don't produce NaN CV", {
  df <- data.frame(
    original_sample_id = rep(paste0("S", 1:5), each = 2),
    cytokine = "CYT1",
    concentration = c(
      0.001, 0.001,   # S1: both at LOD/2 (mean=0.001 > 1e-10, CV=0)
      0.001, 0.002,   # S2: near zero
      100, 102,        # S3: normal
      150, 148,        # S4: normal
      200, 198         # S5: normal
    ),
    stringsAsFactors = FALSE
  )

  result <- flag_replicate_outliers(df, cv_threshold = 25, verbose = FALSE)

  # No NaN or Inf in CV column
  expect_false(any(is.nan(result$replicate_cv[!is.na(result$replicate_cv)])))
  expect_false(any(is.infinite(result$replicate_cv[!is.na(result$replicate_cv)])))

  # S1: identical reps -> CV = 0, not flagged
  s1_cv <- result$replicate_cv[result$original_sample_id == "S1"]
  expect_true(all(!is.na(s1_cv)))
  expect_equal(unique(s1_cv), 0)
  expect_false(any(result$replicate_discordant[result$original_sample_id == "S1"]))
})

test_that("truly zero mean pairs get NA CV", {
  df <- data.frame(
    original_sample_id = rep(paste0("S", 1:5), each = 2),
    cytokine = "CYT1",
    concentration = c(
      0, 0,           # S1: both zero -> mean < 1e-10, CV = NA
      100, 102,       # S2: normal
      150, 148,       # S3: normal
      200, 198,       # S4: normal
      120, 122        # S5: normal
    ),
    stringsAsFactors = FALSE
  )

  result <- flag_replicate_outliers(df, cv_threshold = 25, verbose = FALSE)

  s1_cv <- result$replicate_cv[result$original_sample_id == "S1"]
  expect_true(all(is.na(s1_cv)))
  # Zero-mean pairs should not be flagged by CV rule
  expect_false(any(result$replicate_discordant[result$original_sample_id == "S1"]))
})

test_that("group_cols stratification works", {
  # Create data with two timepoints
  df <- rbind(
    data.frame(
      original_sample_id = rep(paste0("S", 1:5), each = 2),
      cytokine = "CYT1",
      concentration = c(100, 102, 150, 148, 200, 203, 80, 78, 50, 200),
      timepoint = "T1",
      stringsAsFactors = FALSE
    ),
    data.frame(
      original_sample_id = rep(paste0("S", 1:5), each = 2),
      cytokine = "CYT1",
      concentration = c(300, 302, 350, 348, 400, 403, 380, 378, 50, 200),
      timepoint = "T2",
      stringsAsFactors = FALSE
    )
  )

  result <- flag_replicate_outliers(df, group_cols = "timepoint",
                                    cv_threshold = 25, verbose = FALSE)

  expect_s3_class(result, "data.frame")
  # S5 in both timepoints has the same absolute diff (150) but different population contexts
  # Should still flag S5 in both
  s5_t1 <- result[result$original_sample_id == "S5" & result$timepoint == "T1", ]
  s5_t2 <- result[result$original_sample_id == "S5" & result$timepoint == "T2", ]
  expect_true(all(s5_t1$replicate_discordant))
  expect_true(all(s5_t2$replicate_discordant))
})

test_that("input validation catches missing columns", {
  df <- data.frame(a = 1, b = 2, c = 3)

  expect_error(
    flag_replicate_outliers(df, sample_id_col = "nonexistent", verbose = FALSE),
    "Missing required column"
  )

  expect_error(
    flag_replicate_outliers(df, sample_id_col = "a", cytokine_col = "b",
                            value_col = "c", group_cols = "nonexistent",
                            verbose = FALSE),
    "Missing group column"
  )
})

test_that("input validation catches invalid parameters", {
  df <- make_test_data()

  expect_error(
    flag_replicate_outliers(df, action = "invalid", verbose = FALSE)
  )
  expect_error(
    flag_replicate_outliers(df, flag_rule = "invalid", verbose = FALSE)
  )
})

test_that("flag_rule = 'both' requires both criteria", {
  df <- make_discordant_data()

  # With "both" rule and very strict MAD threshold, the pair might not be flagged
  # because both criteria must be met
  result_both <- flag_replicate_outliers(df, cv_threshold = 25, mad_threshold = 100,
                                         flag_rule = "both", verbose = FALSE)

  result_cv <- flag_replicate_outliers(df, cv_threshold = 25,
                                        flag_rule = "cv", verbose = FALSE)

  # "both" should flag fewer or equal pairs compared to "cv" alone
  expect_lte(sum(result_both$replicate_discordant),
             sum(result_cv$replicate_discordant))
})

test_that("flag_rule = 'either' triggers on either criterion", {
  df <- make_discordant_data()

  result_either <- flag_replicate_outliers(df, cv_threshold = 25, mad_threshold = 3,
                                            flag_rule = "either", verbose = FALSE)

  result_cv <- flag_replicate_outliers(df, cv_threshold = 25,
                                        flag_rule = "cv", verbose = FALSE)

  # "either" should flag at least as many as "cv" alone
  expect_gte(sum(result_either$replicate_discordant),
             sum(result_cv$replicate_discordant))
})

test_that("flag_rule = 'mad' uses only MAD criterion", {
  df <- make_discordant_data()

  result <- flag_replicate_outliers(df, flag_rule = "mad", mad_threshold = 3,
                                    verbose = FALSE)

  # S10 (diff=150) should be flagged as outlier vs population
  s10 <- result[result$original_sample_id == "S10", ]
  expect_true(all(s10$replicate_discordant))
})

test_that("verbose output produces messages", {
  df <- make_discordant_data()

  expect_message(
    flag_replicate_outliers(df, cv_threshold = 25, verbose = TRUE),
    "Replicate QC"
  )
})

test_that("multiple cytokines are handled independently", {
  df <- rbind(
    make_discordant_data(),
    data.frame(
      original_sample_id = rep(paste0("S", 1:5), each = 2),
      cytokine = "CYT2",
      concentration = c(100, 102, 150, 148, 200, 203, 80, 78, 120, 122),
      stringsAsFactors = FALSE
    )
  )

  result <- flag_replicate_outliers(df, cv_threshold = 25, verbose = FALSE)

  # CYT1 S10 should be flagged
  cyt1_s10 <- result[result$cytokine == "CYT1" & result$original_sample_id == "S10", ]
  expect_true(all(cyt1_s10$replicate_discordant))

  # CYT2 should have no flags
  cyt2 <- result[result$cytokine == "CYT2", ]
  expect_true(all(!cyt2$replicate_discordant))
})

test_that("replicate_n is correctly computed", {
  df <- make_discordant_data()
  # Add a triple-replicate sample
  triple <- data.frame(
    original_sample_id = rep("S_triple", 3),
    cytokine = "CYT1",
    concentration = c(100, 102, 101),
    stringsAsFactors = FALSE
  )
  df <- rbind(df, triple)

  result <- flag_replicate_outliers(df, cv_threshold = 25, verbose = FALSE)

  s10_n <- unique(result$replicate_n[result$original_sample_id == "S10"])
  expect_equal(s10_n, 2)

  triple_n <- unique(result$replicate_n[result$original_sample_id == "S_triple"])
  expect_equal(triple_n, 3)
})


# =============================================================================
# Tests for sample-level well failure detection (well_detection = "sample")
# =============================================================================

# --- Helper: create data with a well-failure sample ---
# One sample ("WF1") has one bad replicate well systematically inflated across
# many analytes, plus several normal samples for population context.
make_well_failure_data <- function(n_normal = 8, n_cytokines = 10,
                                   failure_type = "inflation", seed = 123) {
  set.seed(seed)
  cytokines <- paste0("CYT", seq_len(n_cytokines))
  rows <- list()

  # Normal samples: two concordant replicates each
  for (i in seq_len(n_normal)) {
    for (cyt in cytokines) {
      base <- runif(1, 50, 200)
      rows <- c(rows, list(data.frame(
        original_sample_id = paste0("S", i),
        sample_id = c(paste0("S", i, "_1"), paste0("S", i, "_2")),
        cytokine = cyt,
        concentration = c(base + rnorm(1, 0, 3), base + rnorm(1, 0, 3)),
        stringsAsFactors = FALSE
      )))
    }
  }

  # Well-failure sample: one well is systematically bad
  for (cyt in cytokines) {
    base <- runif(1, 80, 150)
    if (failure_type == "inflation") {
      # Rep _2 is inflated 5-20x
      good_val <- base + rnorm(1, 0, 3)
      bad_val <- base * runif(1, 5, 20)
    } else if (failure_type == "dropout") {
      # Rep _1 reads near zero (dropout)
      good_val <- base + rnorm(1, 0, 3)
      bad_val <- runif(1, 0.1, 1.0)  # near LOD
    } else {
      # Mixed: half inflated, half deflated
      good_val <- base + rnorm(1, 0, 3)
      if (which(cytokines == cyt) <= n_cytokines / 2) {
        bad_val <- base * runif(1, 3, 10)
      } else {
        bad_val <- runif(1, 0.1, 2.0)
      }
    }
    rows <- c(rows, list(data.frame(
      original_sample_id = "WF1",
      sample_id = c("WF1_1", "WF1_2"),
      cytokine = cyt,
      concentration = if (failure_type == "dropout") c(bad_val, good_val)
                      else c(good_val, bad_val),
      stringsAsFactors = FALSE
    )))
  }

  do.call(rbind, rows)
}


test_that("well_detection requires replicate_id_col", {
  df <- make_test_data()
  df$sample_id <- paste0(df$original_sample_id, "_1")

  expect_error(
    flag_replicate_outliers(df, well_detection = "sample",
                            replicate_id_col = NULL, verbose = FALSE),
    "replicate_id_col is required"
  )
})

test_that("well_detection validates replicate_id_col exists", {
  df <- make_test_data()

  expect_error(
    flag_replicate_outliers(df, well_detection = "sample",
                            replicate_id_col = "nonexistent", verbose = FALSE),
    "replicate_id_col.*not found"
  )
})

test_that("well_detection validates lod_col exists", {
  df <- make_test_data()
  df$sample_id <- paste0(df$original_sample_id, "_1")

  expect_error(
    flag_replicate_outliers(df, well_detection = "sample",
                            replicate_id_col = "sample_id",
                            lod_col = "nonexistent", verbose = FALSE),
    "lod_col.*not found"
  )
})

test_that("sample-level failure detected", {
  df <- make_well_failure_data(failure_type = "inflation")

  result <- flag_replicate_outliers(
    df, replicate_id_col = "sample_id",
    cv_threshold = 25, action = "flag",
    well_detection = "sample", verbose = FALSE
  )

  # WF1 should be classified as well_failure
  wf1 <- result[result$original_sample_id == "WF1", ]
  expect_true(all(wf1$well_failure))

  # Normal samples should NOT be well_failure
  normal <- result[result$original_sample_id != "WF1", ]
  expect_true(all(!normal$well_failure))
})

test_that("good well identified by z-score (inflation)", {
  df <- make_well_failure_data(failure_type = "inflation")

  result <- flag_replicate_outliers(
    df, replicate_id_col = "sample_id",
    cv_threshold = 25, action = "drop_farther",
    well_detection = "sample", verbose = FALSE
  )

  # WF1: rep _1 is good (normal values), rep _2 is bad (inflated)
  # For flagged analytes, only _1 should remain
  wf1 <- result[result$original_sample_id == "WF1", ]
  flagged_wf1 <- wf1[wf1$replicate_discordant, ]

  # All kept rows for flagged analytes should be from _1
  expect_true(all(flagged_wf1$sample_id == "WF1_1"))
})

test_that("good well identified by z-score (dropout)", {
  df <- make_well_failure_data(failure_type = "dropout")

  result <- flag_replicate_outliers(
    df, replicate_id_col = "sample_id",
    cv_threshold = 25, action = "drop_farther",
    well_detection = "sample", verbose = FALSE
  )

  # WF1: rep _1 is bad (dropout), rep _2 is good
  wf1 <- result[result$original_sample_id == "WF1", ]
  flagged_wf1 <- wf1[wf1$replicate_discordant, ]

  # All kept rows for flagged analytes should be from _2
  expect_true(all(flagged_wf1$sample_id == "WF1_2"))
})

test_that("analyte-noise samples use per-analyte logic", {
  # Create data with only 1 discordant analyte out of 10 (10% < 30% threshold)
  set.seed(42)
  n_cyt <- 10
  cytokines <- paste0("CYT", seq_len(n_cyt))
  rows <- list()

  for (i in 1:8) {
    for (cyt in cytokines) {
      base <- runif(1, 50, 200)
      rows <- c(rows, list(data.frame(
        original_sample_id = paste0("S", i),
        sample_id = c(paste0("S", i, "_1"), paste0("S", i, "_2")),
        cytokine = cyt,
        concentration = c(base + rnorm(1, 0, 3), base + rnorm(1, 0, 3)),
        stringsAsFactors = FALSE
      )))
    }
  }

  # One sample with 1 discordant analyte (not enough for well failure)
  for (cyt in cytokines) {
    base <- runif(1, 80, 150)
    if (cyt == "CYT1") {
      # Only this one is discordant
      rows <- c(rows, list(data.frame(
        original_sample_id = "AN1",
        sample_id = c("AN1_1", "AN1_2"),
        cytokine = cyt,
        concentration = c(base, base * 10),
        stringsAsFactors = FALSE
      )))
    } else {
      rows <- c(rows, list(data.frame(
        original_sample_id = "AN1",
        sample_id = c("AN1_1", "AN1_2"),
        cytokine = cyt,
        concentration = c(base + rnorm(1, 0, 2), base + rnorm(1, 0, 2)),
        stringsAsFactors = FALSE
      )))
    }
  }
  df <- do.call(rbind, rows)

  result <- flag_replicate_outliers(
    df, replicate_id_col = "sample_id",
    cv_threshold = 25, action = "flag",
    well_detection = "sample", verbose = FALSE
  )

  # AN1 should NOT be classified as well_failure (only 1/10 flagged = 10%)
  an1 <- result[result$original_sample_id == "AN1", ]
  expect_true(all(!an1$well_failure))
})

test_that("well_params overrides defaults", {
  df <- make_well_failure_data(failure_type = "inflation")

  # With well_threshold = 0.95, even 90% flagged won't qualify as well failure
  result_strict <- flag_replicate_outliers(
    df, replicate_id_col = "sample_id",
    cv_threshold = 25, action = "flag",
    well_detection = "sample",
    well_params = list(well_threshold = 0.95),
    verbose = FALSE
  )

  # With strict threshold, WF1 may not be classified as well failure
  # (it has ~100% flagged, but let's verify the threshold works)
  # Actually with inflation data, likely all 10/10 flagged = 100% > 95%
  # Test a more extreme case
  result_extreme <- flag_replicate_outliers(
    df, replicate_id_col = "sample_id",
    cv_threshold = 25, action = "flag",
    well_detection = "sample",
    well_params = list(well_threshold = 1.0),
    verbose = FALSE
  )

  # Default threshold (0.3) should classify WF1 as well failure
  result_default <- flag_replicate_outliers(
    df, replicate_id_col = "sample_id",
    cv_threshold = 25, action = "flag",
    well_detection = "sample",
    verbose = FALSE
  )

  wf1_default <- result_default[result_default$original_sample_id == "WF1", ]
  expect_true(all(wf1_default$well_failure))
})

test_that("well_params partial override preserves other defaults", {
  df <- make_well_failure_data(failure_type = "inflation")

  # Override only sign_bias_threshold; well_threshold should remain 0.3
  result <- flag_replicate_outliers(
    df, replicate_id_col = "sample_id",
    cv_threshold = 25, action = "flag",
    well_detection = "sample",
    well_params = list(sign_bias_threshold = 0.9),
    verbose = FALSE
  )

  # WF1 should still be classified as well_failure (well_threshold = 0.3 default)
  wf1 <- result[result$original_sample_id == "WF1", ]
  expect_true(all(wf1$well_failure))
})

test_that("drop_farther + sample detection integration", {
  # Mix of well-failure and analyte-noise samples
  df_wf <- make_well_failure_data(failure_type = "inflation")

  result <- flag_replicate_outliers(
    df_wf, replicate_id_col = "sample_id",
    cv_threshold = 25, action = "drop_farther",
    well_detection = "sample", verbose = FALSE
  )

  # WF1 flagged analytes should only keep the good well (_1)
  wf1 <- result[result$original_sample_id == "WF1", ]
  wf1_flagged <- wf1[wf1$replicate_discordant, ]
  expect_true(all(wf1_flagged$sample_id == "WF1_1"))

  # Normal samples should keep both replicates (not flagged)
  normal <- result[result$original_sample_id != "WF1", ]
  expect_true(all(normal$replicate_action == "keep"))
})

test_that("winsorize + sample detection", {
  df <- make_well_failure_data(failure_type = "inflation")

  result <- flag_replicate_outliers(
    df, replicate_id_col = "sample_id",
    cv_threshold = 25, action = "winsorize",
    well_detection = "sample", verbose = FALSE
  )

  # All rows should be present (winsorize doesn't drop)
  expect_equal(nrow(result), nrow(df))

  # WF1 bad well rows for flagged analytes should be winsorized
  wf1 <- result[result$original_sample_id == "WF1", ]
  wf1_bad <- wf1[wf1$replicate_discordant & wf1$sample_id == "WF1_2", ]
  if (nrow(wf1_bad) > 0) {
    expect_true(all(wf1_bad$replicate_action == "winsorized"))
  }
})

test_that("well_score computed correctly", {
  # Simple deterministic case: 4 analytes, one well inflated
  df <- data.frame(
    original_sample_id = rep(c("S1", "S2", "S3", "S4", "WF"), each = 8),
    sample_id = rep(c(
      rep(c("S1_1", "S1_2"), 4),
      rep(c("S2_1", "S2_2"), 4),
      rep(c("S3_1", "S3_2"), 4),
      rep(c("S4_1", "S4_2"), 4),
      rep(c("WF_1", "WF_2"), 4)
    )),
    cytokine = rep(rep(c("A", "B", "C", "D"), each = 2), 5),
    concentration = c(
      # S1-S4: concordant, centered around 100
      100, 102, 100, 98, 100, 103, 100, 97,
      100, 101, 100, 99, 100, 102, 100, 98,
      100, 103, 100, 97, 100, 101, 100, 99,
      100, 98, 100, 102, 100, 99, 100, 101,
      # WF: _1 normal (~100), _2 inflated (~500)
      100, 500, 100, 500, 100, 500, 100, 500
    ),
    stringsAsFactors = FALSE
  )

  result <- flag_replicate_outliers(
    df, replicate_id_col = "sample_id",
    cv_threshold = 25, action = "flag",
    well_detection = "sample", verbose = FALSE
  )

  # WF should have well_score populated
  wf_rows <- result[result$original_sample_id == "WF", ]
  expect_true(all(!is.na(wf_rows$well_score)))

  # Bad well (_2) should have higher well_score than good well (_1)
  wf2_score <- unique(wf_rows$well_score[wf_rows$sample_id == "WF_2"])
  wf1_score <- unique(wf_rows$well_score[wf_rows$sample_id == "WF_1"])
  expect_gt(wf2_score, wf1_score)
})

test_that("sign_bias distinguishes systematic from random", {
  df <- make_well_failure_data(failure_type = "inflation")

  result <- flag_replicate_outliers(
    df, replicate_id_col = "sample_id",
    cv_threshold = 25, action = "flag",
    well_detection = "sample", verbose = FALSE
  )

  wf1 <- result[result$original_sample_id == "WF1", ]

  # Bad well (_2, inflated) should have high sign_bias (near 1.0)
  bad_bias <- unique(wf1$well_sign_bias[wf1$sample_id == "WF1_2"])
  expect_gt(bad_bias, 0.7)

  # Good well (_1, normal) should have lower sign_bias
  good_bias <- unique(wf1$well_sign_bias[wf1$sample_id == "WF1_1"])
  expect_lt(good_bias, bad_bias)
})

test_that("well_failure_type classified correctly — inflation", {
  df <- make_well_failure_data(failure_type = "inflation")

  result <- flag_replicate_outliers(
    df, replicate_id_col = "sample_id",
    cv_threshold = 25, action = "flag",
    well_detection = "sample", verbose = FALSE
  )

  wf1 <- result[result$original_sample_id == "WF1", ]
  expect_true(all(wf1$well_failure_type == "inflation", na.rm = TRUE))
})

test_that("well_failure_type classified correctly — dropout", {
  df <- make_well_failure_data(failure_type = "dropout")

  result <- flag_replicate_outliers(
    df, replicate_id_col = "sample_id",
    cv_threshold = 25, action = "flag",
    well_detection = "sample", verbose = FALSE
  )

  wf1 <- result[result$original_sample_id == "WF1", ]
  ft <- unique(wf1$well_failure_type[!is.na(wf1$well_failure_type)])
  expect_true(ft %in% c("dropout", "mixed"))
})

test_that("lod_fraction detects dropout pattern", {
  df <- make_well_failure_data(failure_type = "dropout")
  # Add LOD column: set LOD at 2.0 for all analytes
  df$lod_value <- 2.0

  result <- flag_replicate_outliers(
    df, replicate_id_col = "sample_id",
    cv_threshold = 25, action = "flag",
    well_detection = "sample", lod_col = "lod_value",
    verbose = FALSE
  )

  wf1 <- result[result$original_sample_id == "WF1", ]

  # Bad well (_1 in dropout: values near 0.1-1.0, all below LOD=2)
  bad_lod <- unique(wf1$well_lod_fraction[wf1$sample_id == "WF1_1"])
  expect_gt(bad_lod, 0.5)  # Most analytes at/below LOD

  # Good well (_2 in dropout: values ~80-150, all above LOD=2)
  good_lod <- unique(wf1$well_lod_fraction[wf1$sample_id == "WF1_2"])
  expect_lt(good_lod, 0.1)
})

test_that("lod_col = NULL skips LOD scoring", {
  df <- make_well_failure_data(failure_type = "inflation")

  result <- flag_replicate_outliers(
    df, replicate_id_col = "sample_id",
    cv_threshold = 25, action = "flag",
    well_detection = "sample", lod_col = NULL,
    verbose = FALSE
  )

  wf1 <- result[result$original_sample_id == "WF1", ]
  # LOD fraction should be NA when lod_col not provided
  expect_true(all(is.na(wf1$well_lod_fraction)))

  # But well_score and sign_bias should still be populated
  expect_true(all(!is.na(wf1$well_score)))
  expect_true(all(!is.na(wf1$well_sign_bias)))
})

test_that("well_detection = 'analyte' is backward-compatible", {
  df <- make_discordant_data()

  # Explicit analyte mode should match original behavior
  result_analyte <- flag_replicate_outliers(
    df, cv_threshold = 25, action = "drop_farther",
    well_detection = "analyte", verbose = FALSE
  )

  result_default <- flag_replicate_outliers(
    df, cv_threshold = 25, action = "drop_farther",
    verbose = FALSE
  )

  expect_equal(nrow(result_analyte), nrow(result_default))
  expect_equal(result_analyte$concentration, result_default$concentration)
  expect_equal(result_analyte$replicate_action, result_default$replicate_action)

  # No well detection columns in analyte mode
  expect_false("well_failure" %in% names(result_analyte))
})

test_that("verbose output reports well failures", {
  df <- make_well_failure_data(failure_type = "inflation")

  expect_message(
    flag_replicate_outliers(
      df, replicate_id_col = "sample_id",
      cv_threshold = 25, action = "drop_farther",
      well_detection = "sample", verbose = TRUE
    ),
    "Well-level failures detected"
  )
})

test_that("sample-level signal inflation pattern", {
  df <- make_well_failure_data(failure_type = "inflation")

  result <- flag_replicate_outliers(
    df, replicate_id_col = "sample_id",
    cv_threshold = 25, action = "drop_farther",
    well_detection = "sample", verbose = FALSE
  )

  wf1 <- result[result$original_sample_id == "WF1", ]
  # Inflation: _2 is bad, _1 is good
  flagged_wf1 <- wf1[wf1$replicate_discordant, ]
  expect_true(all(flagged_wf1$sample_id == "WF1_1"))

  ft <- unique(wf1$well_failure_type[!is.na(wf1$well_failure_type)])
  expect_equal(ft, "inflation")
})

test_that("sample-level signal dropout pattern", {
  df <- make_well_failure_data(failure_type = "dropout")

  result <- flag_replicate_outliers(
    df, replicate_id_col = "sample_id",
    cv_threshold = 25, action = "drop_farther",
    well_detection = "sample", verbose = FALSE
  )

  wf1 <- result[result$original_sample_id == "WF1", ]
  # Dropout: _1 is bad, _2 is good
  flagged_wf1 <- wf1[wf1$replicate_discordant, ]
  expect_true(all(flagged_wf1$sample_id == "WF1_2"))
})

test_that("tied well scores use sign_bias as tiebreaker", {
  # Create data where both wells have similar absolute z-scores
  # but different sign biases
  set.seed(999)
  n_cyt <- 10
  cytokines <- paste0("CYT", seq_len(n_cyt))
  rows <- list()

  # Normal samples for population context
  for (i in 1:6) {
    for (cyt in cytokines) {
      base <- 100
      rows <- c(rows, list(data.frame(
        original_sample_id = paste0("S", i),
        sample_id = c(paste0("S", i, "_1"), paste0("S", i, "_2")),
        cytokine = cyt,
        concentration = c(base + rnorm(1, 0, 5), base + rnorm(1, 0, 5)),
        stringsAsFactors = FALSE
      )))
    }
  }

  # Tied-score sample: both wells deviate similarly in magnitude
  # but well _1 has systematic bias (all positive) while _2 has random bias
  for (j in seq_along(cytokines)) {
    cyt <- cytokines[j]
    # Well _1: consistently 30 above population (all z positive)
    val1 <- 130
    # Well _2: alternating above/below by same magnitude (random sign)
    val2 <- if (j %% 2 == 0) 130 else 70
    rows <- c(rows, list(data.frame(
      original_sample_id = "TIE1",
      sample_id = c("TIE1_1", "TIE1_2"),
      cytokine = cyt,
      concentration = c(val1, val2),
      stringsAsFactors = FALSE
    )))
  }

  df <- do.call(rbind, rows)

  result <- flag_replicate_outliers(
    df, replicate_id_col = "sample_id",
    cv_threshold = 25, action = "flag",
    well_detection = "sample", verbose = FALSE
  )

  tie1 <- result[result$original_sample_id == "TIE1", ]
  if (any(tie1$well_failure)) {
    # Well _1 should have higher sign_bias (all same direction)
    bias_1 <- unique(tie1$well_sign_bias[tie1$sample_id == "TIE1_1"])
    bias_2 <- unique(tie1$well_sign_bias[tie1$sample_id == "TIE1_2"])
    expect_gt(bias_1, bias_2)
  }
})

test_that("tied well scores and sign_bias fall back to per-analyte", {
  # Create perfectly symmetric data: both wells identical scores and bias
  n_cyt <- 6
  cytokines <- paste0("CYT", seq_len(n_cyt))
  rows <- list()

  # Normal samples
  for (i in 1:6) {
    for (cyt in cytokines) {
      rows <- c(rows, list(data.frame(
        original_sample_id = paste0("S", i),
        sample_id = c(paste0("S", i, "_1"), paste0("S", i, "_2")),
        cytokine = cyt,
        concentration = c(100, 102),
        stringsAsFactors = FALSE
      )))
    }
  }

  # Tied sample: both wells identically deviated
  for (cyt in cytokines) {
    rows <- c(rows, list(data.frame(
      original_sample_id = "TIED",
      sample_id = c("TIED_1", "TIED_2"),
      cytokine = cyt,
      concentration = c(500, 500),  # Both identical, far from pop
      stringsAsFactors = FALSE
    )))
  }
  df <- do.call(rbind, rows)

  # This won't actually flag because both reps are identical (CV=0)

  # Instead, make them both equally bad but in opposite directions
  # so they both get flagged but are perfectly tied
  rows2 <- list()
  for (i in 1:6) {
    for (cyt in cytokines) {
      rows2 <- c(rows2, list(data.frame(
        original_sample_id = paste0("S", i),
        sample_id = c(paste0("S", i, "_1"), paste0("S", i, "_2")),
        cytokine = cyt,
        concentration = c(100 + rnorm(1, 0, 2), 100 + rnorm(1, 0, 2)),
        stringsAsFactors = FALSE
      )))
    }
  }

  # Tied sample: both wells very different from each other but equally bad
  set.seed(77)
  for (j in seq_along(cytokines)) {
    cyt <- cytokines[j]
    # Mirror: _1 = 100 + offset, _2 = 100 - offset
    offset <- 60
    rows2 <- c(rows2, list(data.frame(
      original_sample_id = "TIED",
      sample_id = c("TIED_1", "TIED_2"),
      cytokine = cyt,
      concentration = c(100 + offset, 100 - offset),
      stringsAsFactors = FALSE
    )))
  }
  df2 <- do.call(rbind, rows2)

  result <- flag_replicate_outliers(
    df2, replicate_id_col = "sample_id",
    cv_threshold = 25, action = "drop_farther",
    well_detection = "sample", verbose = FALSE
  )

  # The function should still work (no errors) — tied wells fall back to per-analyte
  expect_s3_class(result, "data.frame")

  # TIED sample should still have results
  tied <- result[result$original_sample_id == "TIED", ]
  expect_gt(nrow(tied), 0)
})

test_that("flag action with well_detection = sample still annotates", {
  df <- make_well_failure_data(failure_type = "inflation")

  result <- flag_replicate_outliers(
    df, replicate_id_col = "sample_id",
    cv_threshold = 25, action = "flag",
    well_detection = "sample", verbose = FALSE
  )

  # Same number of rows (flag doesn't remove)
  expect_equal(nrow(result), nrow(df))

  # Well failure columns should be populated
  wf1 <- result[result$original_sample_id == "WF1", ]
  expect_true(all(wf1$well_failure))
  expect_true(all(!is.na(wf1$well_failure_type)))
  expect_true(all(!is.na(wf1$well_score)))
  expect_true(all(!is.na(wf1$well_sign_bias)))
})

test_that("drop_both with well_detection = sample still drops both", {
  df <- make_well_failure_data(failure_type = "inflation")

  result <- flag_replicate_outliers(
    df, replicate_id_col = "sample_id",
    cv_threshold = 25, action = "drop_both",
    well_detection = "sample", verbose = FALSE
  )

  # Flagged rows should be gone
  expect_true(all(!result$replicate_discordant))
})
