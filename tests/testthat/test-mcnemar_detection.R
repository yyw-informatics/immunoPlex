# McNemar Detection Analysis Tests
#
# Validates paired detection frequency analysis using exact McNemar's test:
# 1. Core McNemar detection functionality
# 2. Input validation and error handling
# 3. Custom parameters (column names, confidence levels)
# 4. Edge cases (stable detection, all changes)
# 5. Output structure and statistical correctness

library(testthat)

# Helper function to create test data
create_mcnemar_test_data <- function(n_subjects = 20) {
  set.seed(42)
  
  expand.grid(
    subject_id = 1:n_subjects,
    cytokine = c("IL6", "TNFa", "IL10", "IFNg"),
    timepoint = c("Baseline", "Follow-up")
  ) %>%
    dplyr::mutate(
      # Simulate detection status with different patterns per cytokine
      cens_lod = dplyr::case_when(
        # IL6: Increased detection at follow-up
        cytokine == "IL6" & timepoint == "Baseline" ~ runif(dplyr::n()) > 0.4,
        cytokine == "IL6" & timepoint == "Follow-up" ~ runif(dplyr::n()) > 0.7,
        
        # TNFa: Decreased detection at follow-up
        cytokine == "TNFa" & timepoint == "Baseline" ~ runif(dplyr::n()) > 0.7,
        cytokine == "TNFa" & timepoint == "Follow-up" ~ runif(dplyr::n()) > 0.4,
        
        # IL10: Stable detection
        cytokine == "IL10" ~ runif(dplyr::n()) > 0.5,
        
        # IFNg: Stable non-detection
        cytokine == "IFNg" ~ runif(dplyr::n()) > 0.2,
        
        TRUE ~ TRUE
      )
    )
}

describe("McNemar detection analysis", {
  
  test_that("basic mcnemar_detection executes successfully", {
    skip_if_not_installed("exact2x2")
    skip_if_not_installed("dplyr")
    
    test_data <- create_mcnemar_test_data()
    
    results <- mcnemar_detection(
      data = test_data,
      baseline = "Baseline",
      comparison = "Follow-up",
      quiet = TRUE
    )
    
    # Check output structure
    expect_s3_class(results, "data.frame")
    expect_true(nrow(results) == 4)  # 4 cytokines
    
    # Check required columns
    required_cols <- c("cytokine", "n_pairs", "both_detect", "loss", "gain", 
                      "neither_detect", "prop_baseline", "prop_comparison", 
                      "delta_detection", "p_mcnemar", "q_mcnemar")
    expect_true(all(required_cols %in% names(results)))
    
    # Check data types
    expect_type(results$p_mcnemar, "double")
    expect_type(results$delta_detection, "double")
    expect_type(results$n_pairs, "integer")
  })
  
  test_that("custom column names work correctly", {
    skip_if_not_installed("exact2x2")
    
    test_data <- create_mcnemar_test_data() %>%
      dplyr::rename(
        patient_id = subject_id,
        analyte = cytokine,
        visit = timepoint,
        below_lod = cens_lod
      )
    
    results <- mcnemar_detection(
      data = test_data,
      subject_col = "patient_id",
      cytokine_col = "analyte",
      timepoint_col = "visit",
      censoring_col = "below_lod",
      baseline = "Baseline",
      comparison = "Follow-up",
      quiet = TRUE
    )
    
    expect_s3_class(results, "data.frame")
    expect_true(nrow(results) == 4)
  })
  
  test_that("pre-computed detection status works", {
    skip_if_not_installed("exact2x2")
    
    test_data <- create_mcnemar_test_data() %>%
      dplyr::mutate(is_detected = !cens_lod)
    
    results <- mcnemar_detection(
      data = test_data,
      detection_col = "is_detected",
      censoring_col = NULL,
      baseline = "Baseline",
      comparison = "Follow-up",
      quiet = TRUE
    )
    
    expect_s3_class(results, "data.frame")
    expect_true(all(c("gain", "loss", "both_detect") %in% names(results)))
  })
  
  test_that("custom confidence level and thresholds work", {
    skip_if_not_installed("exact2x2")
    
    test_data <- create_mcnemar_test_data()
    
    results_95 <- mcnemar_detection(
      data = test_data,
      baseline = "Baseline",
      comparison = "Follow-up",
      conf_level = 0.95,
      quiet = TRUE
    )
    
    results_99 <- mcnemar_detection(
      data = test_data,
      baseline = "Baseline",
      comparison = "Follow-up",
      conf_level = 0.99,
      quiet = TRUE
    )
    
    # 99% CIs should be wider than 95% CIs
    ci_width_95 <- mean(results_95$delta_ci_hi - results_95$delta_ci_lo, na.rm = TRUE)
    ci_width_99 <- mean(results_99$delta_ci_hi - results_99$delta_ci_lo, na.rm = TRUE)
    expect_true(ci_width_99 > ci_width_95)
    
    # P-values should be the same
    expect_equal(results_95$p_mcnemar, results_99$p_mcnemar)
  })
  
  test_that("edge case: all stable detection", {
    skip_if_not_installed("exact2x2")
    
    test_data <- create_mcnemar_test_data() %>%
      dplyr::mutate(cens_lod = FALSE)  # All detected
    
    results <- mcnemar_detection(
      data = test_data,
      baseline = "Baseline",
      comparison = "Follow-up",
      quiet = TRUE
    )
    
    # No changes expected
    expect_true(all(results$gain == 0))
    expect_true(all(results$loss == 0))
    expect_true(all(results$p_mcnemar == 1.0))
  })
  
  test_that("error handling: missing required column", {
    test_data <- create_mcnemar_test_data() %>%
      dplyr::select(-cytokine)
    
    expect_error(
      mcnemar_detection(
        data = test_data,
        baseline = "Baseline",
        comparison = "Follow-up",
        quiet = TRUE
      ),
      "Missing required columns"
    )
  })
  
  test_that("error handling: invalid timepoint", {
    test_data <- create_mcnemar_test_data()
    
    expect_error(
      mcnemar_detection(
        data = test_data,
        baseline = "InvalidTimepoint",
        comparison = "Follow-up",
        quiet = TRUE
      ),
      "not found in timepoint column"
    )
  })
  
  test_that("statistical correctness: contingency table consistency", {
    skip_if_not_installed("exact2x2")
    
    test_data <- create_mcnemar_test_data()
    
    results <- mcnemar_detection(
      data = test_data,
      baseline = "Baseline",
      comparison = "Follow-up",
      quiet = TRUE
    )
    
    # Check that all counts sum to n_pairs
    for (i in 1:nrow(results)) {
      row <- results[i, ]
      total <- row$both_detect + row$loss + row$gain + row$neither_detect
      expect_equal(total, row$n_pairs)
      
      # Check discordant pairs
      expect_equal(row$n_discordant, row$loss + row$gain)
      
      # Check baseline detection
      expect_equal(row$n_baseline_detect, row$both_detect + row$loss)
      
      # Check comparison detection
      expect_equal(row$n_comparison_detect, row$both_detect + row$gain)
    }
  })
  
  test_that("FDR correction methods work", {
    skip_if_not_installed("exact2x2")
    
    test_data <- create_mcnemar_test_data()
    
    results_bh <- mcnemar_detection(
      data = test_data,
      baseline = "Baseline",
      comparison = "Follow-up",
      fdr_method = "BH",
      quiet = TRUE
    )
    
    results_bonf <- mcnemar_detection(
      data = test_data,
      baseline = "Baseline",
      comparison = "Follow-up",
      fdr_method = "bonferroni",
      quiet = TRUE
    )
    
    # Bonferroni should be more conservative (higher or equal q-values, ignoring NAs)
    # Remove NA values for comparison
    valid_idx <- !is.na(results_bonf$q_mcnemar) & !is.na(results_bh$q_mcnemar)
    expect_true(all(results_bonf$q_mcnemar[valid_idx] >= results_bh$q_mcnemar[valid_idx] - 1e-10))
  })

  test_that("default baseline/comparison: character column picks alphabetically and messages", {
    skip_if_not_installed("exact2x2")

    test_data <- create_mcnemar_test_data() %>%
      dplyr::mutate(timepoint = ifelse(timepoint == "Baseline", "Pre", "Post"))
    expect_type(test_data$timepoint, "character")

    msgs <- testthat::capture_messages(
      mcnemar_detection(
        data = test_data,
        baseline = NULL,
        comparison = NULL,
        quiet = TRUE
      )
    )

    # sort(c("Pre", "Post")) = c("Post", "Pre"), so baseline -> "Post"
    expect_true(any(grepl("Using baseline timepoint: Post", msgs, fixed = TRUE)))
    expect_true(any(grepl("Using comparison timepoint: Pre", msgs, fixed = TRUE)))
    expect_true(any(grepl("alphabetical order", msgs)))
  })

  test_that(">2 timepoints with defaults emits a warning naming the dropped levels", {
    skip_if_not_installed("exact2x2")

    base <- create_mcnemar_test_data()
    extra <- base %>%
      dplyr::filter(timepoint == "Follow-up") %>%
      dplyr::mutate(timepoint = "Recovery")
    test_data <- dplyr::bind_rows(base, extra)

    expect_warning(
      suppressMessages(
        mcnemar_detection(
          data = test_data,
          baseline = NULL,
          comparison = NULL,
          quiet = TRUE
        )
      ),
      "Dropped: Recovery"
    )
  })

  test_that("default baseline/comparison: ordered factor uses level order without alpha warning", {
    skip_if_not_installed("exact2x2")

    test_data <- create_mcnemar_test_data() %>%
      dplyr::mutate(timepoint = factor(timepoint,
                                       levels = c("Baseline", "Follow-up"),
                                       ordered = TRUE))

    msgs <- testthat::capture_messages(
      results <- mcnemar_detection(
        data = test_data,
        baseline = NULL,
        comparison = NULL,
        quiet = TRUE
      )
    )

    expect_true(any(grepl("factor level order", msgs)))
    expect_false(any(grepl("alphabetical", msgs)))
    expect_s3_class(results, "mcnemar_detection")
  })

  test_that("return is class-tagged 'mcnemar_detection' and dispatches to print method", {
    skip_if_not_installed("exact2x2")

    test_data <- create_mcnemar_test_data()

    results <- mcnemar_detection(
      data = test_data,
      baseline = "Baseline",
      comparison = "Follow-up",
      quiet = TRUE
    )

    expect_true(inherits(results, "mcnemar_detection"))
    expect_true(inherits(results, "data.frame"))

    out <- utils::capture.output(print(results))
    expect_true(any(grepl("McNemar Detection Analysis Results", out)))
  })

  test_that("method switch: each variant runs and produces relative ordering on a prepared cell", {
    skip_if_not_installed("exact2x2")

    test_data <- create_mcnemar_test_data()

    res_exact <- mcnemar_detection(
      data = test_data, baseline = "Baseline", comparison = "Follow-up",
      method = "exact", quiet = TRUE
    )
    res_midp <- mcnemar_detection(
      data = test_data, baseline = "Baseline", comparison = "Follow-up",
      method = "midp", quiet = TRUE
    )
    res_nocc <- mcnemar_detection(
      data = test_data, baseline = "Baseline", comparison = "Follow-up",
      method = "noCC", quiet = TRUE
    )

    # All three return the same shape and contingency-table arithmetic.
    expect_s3_class(res_exact, "mcnemar_detection")
    expect_s3_class(res_midp, "mcnemar_detection")
    expect_s3_class(res_nocc, "mcnemar_detection")

    expect_equal(nrow(res_exact), nrow(res_midp))
    expect_equal(nrow(res_exact), nrow(res_nocc))

    # The 2x2 table itself is method-independent. Compare per-cytokine.
    by_cyt <- function(x) x[order(x$cytokine), ]
    re <- by_cyt(res_exact); rm <- by_cyt(res_midp); rn <- by_cyt(res_nocc)
    expect_equal(re$gain, rm$gain)
    expect_equal(re$loss, rn$loss)
    expect_equal(re$delta_detection, rm$delta_detection)
    expect_equal(re$delta_detection, rn$delta_detection)

    # Relative ordering: mid-p <= exact is mathematically guaranteed, since
    # mid-p is exact minus the point mass at the observed value (clamped at 0).
    # noCC vs exact has no strict guarantee, so just check both are valid p's.
    has_disc <- re$n_discordant > 0
    expect_true(any(has_disc))
    expect_true(all(rm$p_mcnemar[has_disc] <= re$p_mcnemar[has_disc] + 1e-10))
    expect_true(all(rn$p_mcnemar[has_disc] >= 0 & rn$p_mcnemar[has_disc] <= 1))
    expect_true(all(rm$p_mcnemar[has_disc] >= 0 & rm$p_mcnemar[has_disc] <= 1))

    # Non-exact methods report Haldane-corrected MPOR (always finite when
    # n_discordant > 0); they do not provide an exact MPOR CI -> NA.
    expect_true(all(!is.na(rm$mpor[has_disc])))
    expect_true(all(is.na(rm$mpor_ci_lo[has_disc])))
    expect_true(all(is.na(rn$mpor_ci_lo[has_disc])))
  })

  test_that("method = invalid is rejected by match.arg", {
    skip_if_not_installed("exact2x2")

    test_data <- create_mcnemar_test_data()
    expect_error(
      mcnemar_detection(
        data = test_data, baseline = "Baseline", comparison = "Follow-up",
        method = "wald", quiet = TRUE
      ),
      "should be one of"
    )
  })

  # ---- Degenerate-regime tests (review plan chat 1, gap 6) ------------------
  #
  # Exercise edge cases that operationally arise in immunoassay panels:
  # all-censored analyte, 1-cytokine panel, 0%/100% baseline detection.
  # Previously flagged as "no test probes" in the review plan.

  test_that("degenerate: all-censored analyte returns p=1, n_discordant=0", {
    skip_if_not_installed("exact2x2")

    test_data <- create_mcnemar_test_data()
    # Force one cytokine to be 100% censored at both timepoints
    test_data$cens_lod[test_data$cytokine == "IFNg"] <- TRUE

    results <- mcnemar_detection(
      data = test_data,
      baseline = "Baseline", comparison = "Follow-up",
      quiet = TRUE
    )

    row <- results[results$cytokine == "IFNg", ]
    expect_equal(row$both_detect, 0L)
    expect_equal(row$loss, 0L)
    expect_equal(row$gain, 0L)
    expect_equal(row$n_discordant, 0L)
    expect_equal(row$p_mcnemar, 1.0)
    # Cohen's kappa is defined by the zero-discordant branch (p_observed=1,
    # p_expected=1 => returns 1.0). Check we don't divide-by-zero.
    expect_false(is.nan(row$cohen_kappa))
  })

  test_that("degenerate: 1-cytokine panel -> q_mcnemar equals p_mcnemar", {
    skip_if_not_installed("exact2x2")

    test_data <- create_mcnemar_test_data() %>%
      dplyr::filter(cytokine == "IL6")

    results <- mcnemar_detection(
      data = test_data,
      baseline = "Baseline", comparison = "Follow-up",
      quiet = TRUE
    )

    expect_equal(nrow(results), 1L)
    # BH on a single p-value is the identity
    expect_equal(results$q_mcnemar, results$p_mcnemar)
  })

  test_that("degenerate: 0% baseline detection (everyone censored at T1)", {
    skip_if_not_installed("exact2x2")

    # Single cytokine, n=20 subjects: T1 always censored, T2 half detected
    set.seed(1)
    n <- 20L
    td <- data.frame(
      subject_id = rep(seq_len(n), times = 2L),
      cytokine   = "ONE",
      timepoint  = rep(c("T1", "T2"), each = n),
      cens_lod   = c(rep(TRUE, n),
                     sample(c(TRUE, FALSE), n, replace = TRUE)),
      stringsAsFactors = FALSE
    )

    results <- mcnemar_detection(
      data = td, baseline = "T1", comparison = "T2", quiet = TRUE
    )

    # n_baseline_detect = 0 -> prop_baseline = 0; loss (b) = 0 by construction
    expect_equal(results$n_baseline_detect, 0L)
    expect_equal(results$loss, 0L)
    expect_equal(results$prop_baseline, 0)
    # The only discordance is gain (c); p_mcnemar is the one-sided-ish exact
    # test on c/n_discordant with p=0.5 - should be in (0, 1].
    expect_true(results$p_mcnemar > 0 && results$p_mcnemar <= 1)
    # Clopper-Pearson CI on 0/n is [0, x]
    expect_equal(results$prop_baseline_ci_lo, 0)
  })

  test_that("degenerate: 100% baseline detection (everyone detected at T1)", {
    skip_if_not_installed("exact2x2")

    set.seed(2)
    n <- 20L
    td <- data.frame(
      subject_id = rep(seq_len(n), times = 2L),
      cytokine   = "ONE",
      timepoint  = rep(c("T1", "T2"), each = n),
      cens_lod   = c(rep(FALSE, n),
                     sample(c(TRUE, FALSE), n, replace = TRUE)),
      stringsAsFactors = FALSE
    )

    results <- mcnemar_detection(
      data = td, baseline = "T1", comparison = "T2", quiet = TRUE
    )

    expect_equal(results$n_baseline_detect, n)
    expect_equal(results$prop_baseline, 100)
    # gain (c) = 0 by construction (can't gain if already at 100%)
    expect_equal(results$gain, 0L)
    expect_true(results$p_mcnemar > 0 && results$p_mcnemar <= 1)
  })

  # ---- Guard tests (review plan chat 1, gaps: subject_id collision,
  #                  fdr_method validation) -----------------------------------

  test_that("internal columns do not clobber an unrelated user 'subject_id' column", {
    skip_if_not_installed("exact2x2")

    test_data <- create_mcnemar_test_data() %>%
      dplyr::rename(patient_id = subject_id) %>%
      # Introduce an unrelated 'subject_id' column that MUST NOT be touched
      dplyr::mutate(subject_id = "sentinel_value")

    # Before the fix this would either error or silently rename-over
    # the decoy column inside the pipeline. After the fix (internal .mcn_*
    # names) the call completes and the decoy is untouched in the input.
    results <- mcnemar_detection(
      data         = test_data,
      subject_col  = "patient_id",
      baseline     = "Baseline",
      comparison   = "Follow-up",
      quiet        = TRUE
    )

    expect_s3_class(results, "mcnemar_detection")
    expect_equal(nrow(results), 4L)
    # The decoy column in the INPUT must still be intact (dplyr doesn't mutate
    # by reference, but verify explicitly in case a future refactor changes
    # that).
    expect_true(all(test_data$subject_id == "sentinel_value"))
  })

  test_that("fdr_method: invalid string errors up front, partial match expands", {
    skip_if_not_installed("exact2x2")

    test_data <- create_mcnemar_test_data()

    expect_error(
      mcnemar_detection(
        data = test_data, baseline = "Baseline", comparison = "Follow-up",
        fdr_method = "notamethod", quiet = TRUE
      ),
      "should be one of"
    )

    # Partial match: "bonf" -> "bonferroni"
    res <- mcnemar_detection(
      data = test_data, baseline = "Baseline", comparison = "Follow-up",
      fdr_method = "bonf", quiet = TRUE
    )
    expect_equal(attr(res, "fdr_method"), "bonferroni")
  })

  test_that("fdr_method attribute survives on the return object", {
    skip_if_not_installed("exact2x2")

    test_data <- create_mcnemar_test_data()
    res_bh <- mcnemar_detection(
      data = test_data, baseline = "Baseline", comparison = "Follow-up",
      fdr_method = "BH", quiet = TRUE
    )
    res_bf <- mcnemar_detection(
      data = test_data, baseline = "Baseline", comparison = "Follow-up",
      fdr_method = "bonferroni", quiet = TRUE
    )

    expect_equal(attr(res_bh, "fdr_method"), "BH")
    expect_equal(attr(res_bf, "fdr_method"), "bonferroni")

    # print() header mentions the method (differentiates the two returns)
    out_bh <- utils::capture.output(print(res_bh))
    out_bf <- utils::capture.output(print(res_bf))
    expect_true(any(grepl("BH-adjusted", out_bh, fixed = TRUE)))
    expect_true(any(grepl("bonferroni-adjusted", out_bf, fixed = TRUE)))
  })

  # ---- replicate_agg argument (review plan chat 1, gap 4) -------------------

  test_that("replicate_agg: three rules diverge on asymmetric 2-rep input", {
    skip_if_not_installed("exact2x2")

    # 30 subjects, 1 cytokine, 2 replicates per (subject, timepoint).
    # T1: (detected, censored)  -> perfect replicate disagreement
    # T2: (censored, censored)  -> all-censored
    n <- 30L
    base <- expand.grid(subject_id = seq_len(n),
                        cytokine   = "A",
                        timepoint  = c("T1", "T2"),
                        rep_idx    = c(1L, 2L),
                        KEEP.OUT.ATTRS = FALSE,
                        stringsAsFactors = FALSE)
    base$cens_lod <- NA
    base$cens_lod[base$timepoint == "T1" & base$rep_idx == 1L] <- FALSE
    base$cens_lod[base$timepoint == "T1" & base$rep_idx == 2L] <- TRUE
    base$cens_lod[base$timepoint == "T2"]                       <- TRUE

    r_maj <- mcnemar_detection(
      base, baseline = "T1", comparison = "T2",
      replicate_agg = "majority_vote", quiet = TRUE
    )
    r_any <- mcnemar_detection(
      base, baseline = "T1", comparison = "T2",
      replicate_agg = "any_detected", quiet = TRUE
    )
    r_all <- mcnemar_detection(
      base, baseline = "T1", comparison = "T2",
      replicate_agg = "all_detected", quiet = TRUE
    )

    # Majority vote: mean=0.5 -> round(0.5)=0 -> T1 FALSE for every subject.
    # T2 also FALSE. No discordance.
    expect_equal(r_maj$loss, 0L)
    expect_equal(r_maj$gain, 0L)
    expect_equal(r_maj$p_mcnemar, 1.0)

    # Any: T1 TRUE (1 of 2 reps detected), T2 FALSE -> every subject is a loss.
    expect_equal(r_any$loss, n)
    expect_equal(r_any$gain, 0L)
    expect_true(r_any$p_mcnemar < 0.01)

    # All: T1 requires both reps detected -> FALSE. T2 FALSE. No discordance.
    expect_equal(r_all$loss, 0L)
    expect_equal(r_all$gain, 0L)
    expect_equal(r_all$p_mcnemar, 1.0)
  })

  test_that("replicate_agg: fully-missing replicate group drops the subject for all three rules", {
    # Prior behavior: any_detected on an all-NA group returned FALSE via
    # any(x, na.rm = TRUE), while majority_vote and all_detected returned NA.
    # That asymmetry silently coded fully-missing subjects as below-LOD
    # under any_detected only. Now every rule returns NA for all-NA input,
    # and the pivot_wider + !is.na() filter drops the subject uniformly.
    skip_if_not_installed("exact2x2")

    n <- 10L
    base <- expand.grid(subject_id = seq_len(n),
                        cytokine   = "A",
                        timepoint  = c("T1", "T2"),
                        rep_idx    = c(1L, 2L),
                        KEEP.OUT.ATTRS = FALSE,
                        stringsAsFactors = FALSE)
    base$cens_lod <- FALSE
    # Subject 1 is fully missing at T1 (both reps NA on cens_lod).
    base$cens_lod[base$subject_id == 1L & base$timepoint == "T1"] <- NA

    for (rl in c("majority_vote", "any_detected", "all_detected")) {
      res <- mcnemar_detection(
        base, baseline = "T1", comparison = "T2",
        replicate_agg = rl, quiet = TRUE
      )
      # Subject 1's T1 cell is NA -> pair drops -> n_pairs = n - 1 for every
      # rule. Pre-fix, any_detected kept subject 1 at T1 = FALSE and yielded
      # n_pairs = n.
      expect_equal(res$n_pairs, n - 1L,
                   info = paste0("rule = ", rl))
    }
  })

  test_that("replicate_agg: identity on one-obs-per-triple input", {
    skip_if_not_installed("exact2x2")

    test_data <- create_mcnemar_test_data()  # one row per triple

    r_maj <- mcnemar_detection(
      test_data, baseline = "Baseline", comparison = "Follow-up",
      replicate_agg = "majority_vote", quiet = TRUE
    )
    r_any <- mcnemar_detection(
      test_data, baseline = "Baseline", comparison = "Follow-up",
      replicate_agg = "any_detected", quiet = TRUE
    )
    r_all <- mcnemar_detection(
      test_data, baseline = "Baseline", comparison = "Follow-up",
      replicate_agg = "all_detected", quiet = TRUE
    )

    # With one observation per group, all three rules should yield identical
    # detection calls and therefore identical p-values per cytokine.
    ord <- function(x) x[order(x$cytokine), ]
    expect_equal(ord(r_maj)$p_mcnemar, ord(r_any)$p_mcnemar)
    expect_equal(ord(r_maj)$p_mcnemar, ord(r_all)$p_mcnemar)
  })

  # ---- plate_col confounding check -----------------------------------------

  test_that("plate_col: warns when plate is perfectly confounded with timepoint", {
    skip_if_not_installed("exact2x2")

    # Every 'Baseline' row is on plate P1; every 'Follow-up' row is on P2.
    d <- create_mcnemar_test_data()
    d$plate <- ifelse(d$timepoint == "Baseline", "P1", "P2")

    expect_warning(
      mcnemar_detection(d, baseline = "Baseline", comparison = "Follow-up",
                        plate_col = "plate", quiet = TRUE),
      "perfectly confounded with timepoint"
    )
  })

  test_that("plate_col: does not warn when plate spans both timepoints", {
    skip_if_not_installed("exact2x2")

    d <- create_mcnemar_test_data()
    # Assign plates in a way that every plate sees both timepoints.
    d$plate <- ifelse(d$subject_id %% 2L == 0L, "P1", "P2")

    expect_no_warning(
      mcnemar_detection(d, baseline = "Baseline", comparison = "Follow-up",
                        plate_col = "plate", quiet = TRUE)
    )
  })

  test_that("plate_col: single-plate input does not warn (no confounding possible)", {
    skip_if_not_installed("exact2x2")

    d <- create_mcnemar_test_data()
    d$plate <- "P1"  # one plate for everything

    expect_no_warning(
      mcnemar_detection(d, baseline = "Baseline", comparison = "Follow-up",
                        plate_col = "plate", quiet = TRUE)
    )
  })

  test_that("plate_col: missing column is an error", {
    skip_if_not_installed("exact2x2")
    d <- create_mcnemar_test_data()

    expect_error(
      mcnemar_detection(d, baseline = "Baseline", comparison = "Follow-up",
                        plate_col = "nope", quiet = TRUE),
      "plate_col 'nope' not found"
    )
  })

  test_that("plate_col: NULL (default) preserves legacy behavior (no plate check)", {
    skip_if_not_installed("exact2x2")
    d <- create_mcnemar_test_data()
    d$plate <- ifelse(d$timepoint == "Baseline", "P1", "P2")  # would confound

    # plate_col not supplied -> we should get the standard default-pick
    # message(s) but no plate-confound warning.
    expect_no_warning(
      suppressMessages(
        mcnemar_detection(d, baseline = "Baseline", comparison = "Follow-up",
                          quiet = TRUE)
      )
    )
  })

})


