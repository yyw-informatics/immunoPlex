# Wald CI for the paired-proportion difference delta = (c - b) / n,
# expressed on the percentage scale to match the rest of the return tibble.
# Used by method = "midp" and method = "noCC". Mirrors paired_wald_ci() in the
# B2 benchmark script of the separate benchmark suite, for cross-comparison.
.paired_wald_delta_ci <- function(b, c, n, conf_level = 0.95) {
  delta <- (c - b) / n
  var_di <- (b + c) / n - delta^2
  se <- sqrt(max(var_di / n, 0))
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  list(
    ci_lo = (delta - z * se) * 100,
    ci_hi = (delta + z * se) * 100
  )
}


#' McNemar's Test for Paired Detection Frequency Analysis
#'
#' Performs exact McNemar's test using the exact2x2 package to analyze changes 
#' in detection frequency between two paired timepoints or conditions. Designed 
#' for censored cytokine data where detection status (above/below limit of 
#' detection) may change over time.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats binom.test p.adjust
#' @importFrom utils head
#' @importFrom exact2x2 mcnemar.exact mcnemarExactDP
#'
#' @param data A data frame in long format containing paired observations
#' @param subject_col Character string. Name of the column identifying subjects (default: "subject_id")
#' @param cytokine_col Character string. Name of the column identifying cytokines or analytes (default: "cytokine")
#' @param timepoint_col Character string. Name of the column identifying timepoints/conditions (default: "timepoint")
#' @param detection_col Character string. Name of the column indicating detection status.
#'   If NULL (default), will be derived from censoring_col (detected = !censored)
#' @param censoring_col Character string. Name of the column indicating censoring status
#'   (TRUE = censored/below LOD, FALSE = detected). Only used if detection_col is NULL (default: "cens_lod")
#' @param baseline Character string. Name of the baseline/first timepoint.
#'   If NULL (default), the first level is picked from \code{sort(unique(data[[timepoint_col]]))}
#'   and a \code{message()} reports the choice. For factor / ordered-factor columns
#'   this respects level order; for character columns the choice is alphabetical, which
#'   may not match the semantic baseline (e.g. for c("Pre", "Post") the default is "Post").
#'   Pass \code{baseline} explicitly or supply an ordered factor to make this deterministic.
#' @param comparison Character string. Name of the comparison/second timepoint.
#'   If NULL (default), picked from \code{sort(unique(data[[timepoint_col]]))} after
#'   removing \code{baseline}; same caveats as \code{baseline} apply.
#' @param conf_level Numeric. Confidence level for confidence intervals (default: 0.95)
#' @param detection_change_threshold Numeric. Minimum percentage point change in detection
#'   to classify as "Increased" or "Decreased" detection (default: 15)
#' @param plate_col Character string. Optional name of a column identifying
#'   assay plate (or any other batch/technical layer). When supplied,
#'   \code{mcnemar_detection} checks whether plate is perfectly confounded
#'   with timepoint - i.e., every plate value maps to a single timepoint
#'   across all subjects and cytokines - and emits a \code{warning()} if
#'   so. The function does \strong{not} attempt to adjust for plate: a 2x2
#'   paired test has no free parameter for plate drift, so confounding
#'   with timepoint cannot be separated from the biological effect.
#'   Residualize against plate upstream (e.g., via
#'   \code{plsda_preprocess}'s \code{covariates} argument, or a
#'   per-analyte linear model) if you need to account for it. Default:
#'   \code{NULL}.
#' @param replicate_agg Character string. How to collapse technical replicates
#'   sharing the same \code{(subject, cytokine, timepoint)} key to a single
#'   detection call. One of:
#'   \itemize{
#'     \item \code{"majority_vote"} (default, backward compatible):
#'       \code{as.logical(round(mean(detected, na.rm = TRUE)))}. Ties hit R's
#'       IEC-60559 banker's rounding (\code{round(0.5) == 0}) so an even-count
#'       replicate split 50/50 resolves to \code{FALSE}.
#'     \item \code{"any_detected"}: \code{any(detected)} on the non-NA
#'       replicates. A subject is called detected if at least one replicate
#'       was detected. Most sensitive; protective against false-censoring
#'       from a bad replicate.
#'     \item \code{"all_detected"}: \code{all(detected)} on the non-NA
#'       replicates. Requires every replicate to be detected. Most specific;
#'       conservative against false-detection from a spurious above-LOD
#'       replicate.
#'   }
#'   All three rules return \code{NA} when a \code{(subject, cytokine,
#'   timepoint)} group has zero non-NA replicates, rather than emitting
#'   an informative default from the empty set
#'   (\code{any(logical(0))} is \code{FALSE}, \code{all(logical(0))} is
#'   \code{TRUE}). Such a subject then falls out of \code{paired_detection}
#'   via the \code{!is.na(...)} filter instead of being silently coded
#'   below-LOD.
#'   The degenerate case of one observation per
#'   \code{(subject, cytokine, timepoint)} triple collapses all three rules
#'   to the identity. Calibration of the three rules under a realistic
#'   replicate-emitting DGP is a separate benchmark
#'   (review plan section 6c, S1 follow-up).
#' @param fdr_method Character string. Method for multiplicity correction,
#'   validated against \code{\link[stats]{p.adjust.methods}} via
#'   \code{\link{match.arg}} (partial matches accepted, e.g. "bonf"). Default
#'   "BH" (Benjamini-Hochberg). The chosen method is attached to the return
#'   object as \code{attr(., "fdr_method")} so that downstream code can tell
#'   BH-adjusted q-values apart from bonferroni-adjusted ones (both arrive in
#'   the same \code{q_mcnemar} column).
#' @param method Character string. McNemar test variant. One of:
#'   \itemize{
#'     \item \code{"exact"} (default): exact McNemar p-value via
#'       \code{exact2x2::mcnemar.exact} and exact CI for Delta via
#'       \code{exact2x2::mcnemarExactDP}. Conservative, especially at small n.
#'     \item \code{"midp"}: mid-p correction to the exact two-sided p-value on
#'       discordant pairs (\code{max(p_exact - dbinom(c, b+c, 0.5), 0)}).
#'       Delta CI reported as Wald on the paired-proportion difference.
#'     \item \code{"noCC"}: asymptotic McNemar chi-squared without continuity
#'       correction (\code{chi2 = (b - c)^2 / (b + c)}). Wald CI for Delta.
#'   }
#'   The B2 calibration benchmark finds that \code{"midp"} and
#'   \code{"noCC"} are better calibrated at small n while \code{"exact"}
#'   is conservative; we keep \code{"exact"} as the default since
#'   conservatism fails safely.
#'   For non-exact methods, MPOR is reported as the Haldane-corrected
#'   \code{(c + 0.5) / (b + 0.5)} to match the B2 benchmark and to stay
#'   finite when \code{b = 0}.
#' @param quiet Logical. If TRUE, suppress progress \code{cat()}s (defaults reporting,
#'   per-cytokine progress, and the closing summary). Note: \code{quiet} does \strong{not}
#'   suppress error/warning signals - backend failures from \code{exact2x2::*} continue
#'   to fire \code{warning()} regardless. (default: FALSE)
#'
#' @return A tibble with one row per cytokine containing:
#'   \item{cytokine}{Cytokine/analyte name}
#'   \item{n_pairs}{Number of complete subject pairs (N)}
#'   \item{both_detect}{Count a: detected at both timepoints}
#'   \item{loss}{Count b: detected at baseline only (lost detection)}
#'   \item{gain}{Count c: detected at comparison only (gained detection)}
#'   \item{neither_detect}{Count d: detected at neither timepoint}
#'   \item{n_discordant}{Number of discordant pairs (b + c)}
#'   \item{n_baseline_detect}{Number detected at baseline (a + b)}
#'   \item{n_comparison_detect}{Number detected at comparison (a + c)}
#'   \item{prop_baseline}{Proportion detected at baseline (percentage)}
#'   \item{prop_baseline_ci_lo}{Lower 95% exact CI for baseline proportion}
#'   \item{prop_baseline_ci_hi}{Upper 95% exact CI for baseline proportion}
#'   \item{prop_comparison}{Proportion detected at comparison (percentage)}
#'   \item{prop_comparison_ci_lo}{Lower 95% exact CI for comparison proportion}
#'   \item{prop_comparison_ci_hi}{Upper 95% exact CI for comparison proportion}
#'   \item{delta_detection}{Paired difference Delta = prop_comparison - prop_baseline (percentage points)}
#'   \item{delta_ci_lo}{Lower 95% exact CI for Delta (from mcnemarExactDP)}
#'   \item{delta_ci_hi}{Upper 95% exact CI for Delta (from mcnemarExactDP)}
#'   \item{rate_upward}{Upward transition rate c/N (percentage)}
#'   \item{rate_upward_ci_lo}{Lower 95% exact CI for upward rate}
#'   \item{rate_upward_ci_hi}{Upper 95% exact CI for upward rate}
#'   \item{rate_downward}{Downward transition rate b/N (percentage)}
#'   \item{rate_downward_ci_lo}{Lower 95% exact CI for downward rate}
#'   \item{rate_downward_ci_hi}{Upper 95% exact CI for downward rate}
#'   \item{mpor}{Matched-pair odds ratio c/b (exact, central method)}
#'   \item{mpor_ci_lo}{Lower 95% exact CI for matched-pair OR}
#'   \item{mpor_ci_hi}{Upper 95% exact CI for matched-pair OR}
#'   \item{p_mcnemar}{Exact McNemar p-value (two-sided, central)}
#'   \item{q_mcnemar}{FDR-corrected q-value (BH method)}
#'   \item{cohen_kappa}{Cohen's kappa (agreement measure)}
#'   \item{detection_pattern}{Categorical pattern based on delta_detection threshold}
#'   \item{mcnemar_significance}{Significance annotation: "***", "**", "*", or ""}
#'
#' @details
#' This function implements a comprehensive paired detection analysis using exact methods
#' from the \code{exact2x2} package, following best practices for McNemar's test.
#'
#' \strong{Statistical Methods:}
#' \itemize{
#'   \item \strong{Detection proportions:} Exact Clopper-Pearson CIs via \code{binom.test()}
#'   \item \strong{Paired difference (Delta):} Exact CI via \code{mcnemarExactDP()} - the primary effect size
#'   \item \strong{Matched-pair OR:} Exact central method via \code{mcnemar.exact()} - secondary effect size
#'   \item \strong{Transition rates:} Exact binomial CIs for upward (c/N) and downward (b/N) rates
#'   \item \strong{Agreement:} Cohen's kappa for time-to-time reliability
#'   \item \strong{Multiplicity:} BH-FDR correction across cytokines
#' }
#'
#' \strong{Contingency Table Structure:}
#' \tabular{lcc}{
#'                \tab Comparison=1 \tab Comparison=0 \cr
#'   Baseline=1   \tab a            \tab b (loss)     \cr
#'   Baseline=0   \tab c (gain)     \tab d
#' }
#'
#' The test focuses on discordant pairs (b + c). Under the null hypothesis,
#' gains and losses are equally likely. The exact methods avoid asymptotic
#' approximations and handle zero cells without pseudo-count corrections.
#'
#' \strong{Replicate aggregation:}
#' If multiple rows share the same \code{(subject, cytokine, timepoint)} key
#' (e.g., technical replicates), they are collapsed to a single detection
#' status according to \code{replicate_agg}. The default
#' (\code{"majority_vote"}) preserves the original behavior:
#' \code{as.logical(round(mean(detected, na.rm = TRUE)))}, where R's IEC-60559
#' banker's rounding sends an even-count 50/50 split to \code{FALSE}
#' (non-detected). \code{"any_detected"} and \code{"all_detected"} offer the
#' two natural alternatives. For full caller control, aggregate upstream and
#' pass the already-aggregated data with one row per triple.
#'
#' \strong{Effect Sizes:}
#' \itemize{
#'   \item \strong{Delta = (c - b)/N:} Interpretable as percentage-point change
#'   \item \strong{MPOR = c/b:} Odds of detection at T2 vs T1 among discordants
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with default column names
#' results <- mcnemar_detection(
#'   data = my_data,
#'   baseline = "Enrollment",
#'   comparison = "Delivery"
#' )
#'
#' # Custom column names
#' results <- mcnemar_detection(
#'   data = my_data,
#'   subject_col = "patient_id",
#'   cytokine_col = "analyte",
#'   timepoint_col = "visit",
#'   censoring_col = "below_lod",
#'   baseline = "Visit1",
#'   comparison = "Visit2",
#'   detection_change_threshold = 10
#' )
#'
#' # If detection status already computed
#' results <- mcnemar_detection(
#'   data = my_data,
#'   detection_col = "is_detected",
#'   censoring_col = NULL,
#'   baseline = "Pre",
#'   comparison = "Post"
#' )
#' }
#'
#' @export
mcnemar_detection <- function(data,
                              subject_col = "subject_id",
                              cytokine_col = "cytokine",
                              timepoint_col = "timepoint",
                              detection_col = NULL,
                              censoring_col = "cens_lod",
                              baseline = NULL,
                              comparison = NULL,
                              conf_level = 0.95,
                              detection_change_threshold = 15,
                              plate_col = NULL,
                              replicate_agg = c("majority_vote",
                                                "any_detected",
                                                "all_detected"),
                              fdr_method = "BH",
                              method = c("exact", "midp", "noCC"),
                              quiet = FALSE) {

  method <- match.arg(method)
  replicate_agg <- match.arg(replicate_agg)

  # Validate fdr_method up front rather than letting stats::p.adjust emit its
  # terse "'arg' should be one of ..." from deep in the stack. match.arg also
  # expands partial matches (e.g. "bonf" -> "bonferroni") for convenience.
  fdr_method <- match.arg(fdr_method, choices = stats::p.adjust.methods)

  # Check required packages
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required. Please install it.", call. = FALSE)
  }
  if (method == "exact" && !requireNamespace("exact2x2", quietly = TRUE)) {
    stop("Package 'exact2x2' is required for method = 'exact'. ",
         "Install it, or choose method = 'midp' or 'noCC'.",
         call. = FALSE)
  }
  
  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame", call. = FALSE)
  }
  
  required_cols <- c(subject_col, cytokine_col, timepoint_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  
  # Determine detection column
  if (is.null(detection_col)) {
    if (is.null(censoring_col) || !censoring_col %in% names(data)) {
      stop("Must provide either 'detection_col' or 'censoring_col'", call. = FALSE)
    }
    if (!quiet) {
      cat("Computing detection status from censoring column:", censoring_col, "\n")
    }
    # Create detection column: detected = !censored
    data$detected <- !data[[censoring_col]]
    detection_col <- "detected"
  } else {
    if (!detection_col %in% names(data)) {
      stop("Detection column '", detection_col, "' not found in data", call. = FALSE)
    }
  }
  
  # Determine timepoints
  timepoints <- unique(data[[timepoint_col]])
  if (length(timepoints) < 2) {
    stop("Need at least 2 timepoints for paired analysis", call. = FALSE)
  }

  baseline_defaulted <- is.null(baseline)
  comparison_defaulted <- is.null(comparison)
  is_ordered <- is.ordered(data[[timepoint_col]])
  is_factor <- is.factor(data[[timepoint_col]])
  sort_basis <- if (is_ordered) {
    "factor level order"
  } else if (is_factor) {
    "factor level order"
  } else {
    "alphabetical order"
  }

  if (baseline_defaulted) {
    baseline <- as.character(sort(timepoints)[1])
    msg <- paste0("Using baseline timepoint: ", baseline,
                  " (default picked by ", sort_basis, ")")
    if (!is_ordered) {
      msg <- paste0(msg,
                    "; pass `baseline` explicitly or supply an ordered factor",
                    " to make this deterministic.")
    }
    message(msg)
  }

  if (comparison_defaulted) {
    comparison <- as.character(setdiff(sort(timepoints), baseline)[1])
    msg <- paste0("Using comparison timepoint: ", comparison,
                  " (default picked by ", sort_basis, ")")
    if (!is_ordered) {
      msg <- paste0(msg,
                    "; pass `comparison` explicitly or supply an ordered factor",
                    " to make this deterministic.")
    }
    message(msg)
  }

  # When defaults are used and >2 timepoints exist, the silent pick drops the rest.
  if ((baseline_defaulted || comparison_defaulted) && length(timepoints) > 2) {
    dropped <- setdiff(as.character(timepoints), c(baseline, comparison))
    warning(
      "More than 2 timepoints found; defaults compared only baseline='",
      baseline, "' vs comparison='", comparison, "'. Dropped: ",
      paste(dropped, collapse = ", "),
      ". Pass `baseline` and `comparison` explicitly to silence this warning.",
      call. = FALSE
    )
  }
  
  if (!baseline %in% timepoints) {
    stop("Baseline '", baseline, "' not found in timepoint column", call. = FALSE)
  }
  if (!comparison %in% timepoints) {
    stop("Comparison '", comparison, "' not found in timepoint column", call. = FALSE)
  }

  # Plate-vs-timepoint confounding check (advisory only). A 2x2 McNemar test
  # has no free parameter for plate drift, so if plate is perfectly confounded
  # with timepoint on the (baseline, comparison) slice, plate drift and the
  # biological change are indistinguishable.
  if (!is.null(plate_col)) {
    if (!plate_col %in% names(data)) {
      stop("plate_col '", plate_col, "' not found in data", call. = FALSE)
    }
    plate_sub <- data[data[[timepoint_col]] %in% c(baseline, comparison),
                      c(plate_col, timepoint_col), drop = FALSE]
    plate_sub <- plate_sub[!is.na(plate_sub[[plate_col]]) &
                             !is.na(plate_sub[[timepoint_col]]), , drop = FALSE]
    if (nrow(plate_sub) > 0L) {
      tp_per_plate <- tapply(
        as.character(plate_sub[[timepoint_col]]),
        as.character(plate_sub[[plate_col]]),
        function(x) length(unique(x))
      )
      # Confounded when there are >= 2 plates AND every plate sits at a single
      # timepoint. A single plate means there is nothing to confound with.
      if (length(tp_per_plate) >= 2L &&
          all(tp_per_plate == 1L, na.rm = TRUE)) {
        warning(
          "Plate is perfectly confounded with timepoint (plate_col='",
          plate_col, "', timepoint_col='", timepoint_col, "'). ",
          "McNemar cannot separate plate drift from biological change; ",
          "consider residualizing against plate upstream.",
          call. = FALSE
        )
      }
    }
  }

  if (!quiet) {
    cat("Performing McNemar's test for paired detection analysis\n")
    cat("  Baseline:", baseline, "\n")
    cat("  Comparison:", comparison, "\n")
    cat("  Number of cytokines:", length(unique(data[[cytokine_col]])), "\n")
    cat("  Number of subjects:", length(unique(data[[subject_col]])), "\n\n")
  }
  
  # Replicate-aggregation rule: collapses multiple rows sharing the same
  # (subject, cytokine, timepoint) key to a single boolean detection call.
  agg_fn <- switch(
    replicate_agg,
    majority_vote = function(x) {
      # mean(NA, na.rm = TRUE) -> NaN -> round -> NaN -> as.logical -> NA,
      # so the explicit guard is redundant here but keeps the three rules
      # visibly symmetric on all-NA input.
      x <- x[!is.na(x)]
      if (length(x) == 0) NA else as.logical(round(mean(x)))
    },
    any_detected  = function(x) {
      # Without the guard, any(logical(0), na.rm = TRUE) returns FALSE,
      # which disagrees with majority_vote and all_detected (both NA).
      # Aligning all three rules on "all-NA group -> NA detection call"
      # so the downstream pivot_wider / filter drops the subject rather
      # than silently coding them as below-LOD.
      x <- x[!is.na(x)]
      if (length(x) == 0) NA else any(x)
    },
    all_detected  = function(x) {
      x <- x[!is.na(x)]
      if (length(x) == 0) NA else all(x)
    }
  )

  # Prepare detection data.
  # Internal column names are prefixed with .mcn_ so the select-rename cannot
  # silently clobber an unrelated user column that happens to be named
  # "subject_id" / "cytokine" / "timepoint" / "detected". The collision was
  # previously undefined-behavior-by-dplyr: the user's column would vanish,
  # replaced by whatever subject_col referenced. Prefixed names are reserved
  # per R naming convention and won't appear in the output schema.
  detection_data <- data %>%
    dplyr::select(
      .mcn_subject   = !!rlang::sym(subject_col),
      .mcn_cytokine  = !!rlang::sym(cytokine_col),
      .mcn_timepoint = !!rlang::sym(timepoint_col),
      .mcn_detected  = !!rlang::sym(detection_col)
    ) %>%
    dplyr::filter(.data$.mcn_timepoint %in% c(baseline, comparison)) %>%
    dplyr::group_by(.data$.mcn_subject, .data$.mcn_cytokine, .data$.mcn_timepoint) %>%
    dplyr::summarise(
      .mcn_detected = agg_fn(.data$.mcn_detected),
      .groups = "drop"
    )

  # Create paired detection data (wide format). After pivot_wider, columns are:
  # .mcn_subject, .mcn_cytokine, <baseline value>, <comparison value>.
  paired_detection <- detection_data %>%
    tidyr::pivot_wider(names_from = ".mcn_timepoint", values_from = ".mcn_detected") %>%
    dplyr::filter(!is.na(.data[[baseline]]) & !is.na(.data[[comparison]]))
  
  # Check if we have any paired data
  if (nrow(paired_detection) == 0) {
    stop("No complete pairs found. Check timepoint names and data structure.", call. = FALSE)
  }
  
  if (!quiet) {
    cat("Found", nrow(paired_detection), "complete paired observations\n")
    cat("Running McNemar's test for each cytokine...\n")
  }
  
  # Perform exact McNemar's test for each cytokine using exact2x2 package
  # Process each cytokine separately
  cytokine_list <- unique(paired_detection$.mcn_cytokine)

  mcnemar_results_list <- lapply(cytokine_list, function(cy) {
    cy_data <- paired_detection %>% dplyr::filter(.data$.mcn_cytokine == cy)
    
    # Build 2x2 contingency table
    # Table structure:
    #           Comparison=TRUE  Comparison=FALSE
    # Baseline=TRUE      a              b (loss)
    # Baseline=FALSE     c (gain)       d
    a <- sum(cy_data[[baseline]] & cy_data[[comparison]], na.rm = TRUE)
    b <- sum(cy_data[[baseline]] & !cy_data[[comparison]], na.rm = TRUE)  # loss
    c <- sum(!cy_data[[baseline]] & cy_data[[comparison]], na.rm = TRUE)  # gain
    d <- sum(!cy_data[[baseline]] & !cy_data[[comparison]], na.rm = TRUE)
    
    n_pairs <- a + b + c + d
    n_discordant <- b + c
    n_baseline_detect <- a + b
    n_comparison_detect <- a + c
    
    # 1) Detection proportions with exact Clopper-Pearson CIs
    prop_baseline <- 100 * n_baseline_detect / n_pairs
    prop_comparison <- 100 * n_comparison_detect / n_pairs
    
    bt_baseline <- binom.test(n_baseline_detect, n_pairs, conf.level = conf_level)
    prop_baseline_ci <- 100 * bt_baseline$conf.int
    
    bt_comparison <- binom.test(n_comparison_detect, n_pairs, conf.level = conf_level)
    prop_comparison_ci <- 100 * bt_comparison$conf.int
    
    # 2) McNemar p-value and matched-pair odds ratio (method-dependent)
    contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    delta_detection <- prop_comparison - prop_baseline

    if (n_discordant == 0) {
      # No discordant pairs - perfect stability: p=1, MPOR undefined, delta CI = 0.
      p_mcnemar <- 1.0
      mpor <- NA_real_
      mpor_ci_lo <- NA_real_
      mpor_ci_hi <- NA_real_
      delta_ci_lo <- 0
      delta_ci_hi <- 0
    } else if (method == "exact") {
      # Call exact2x2::exact2x2 directly rather than mcnemar.exact, because
      # mcnemar.exact uses do.call("exact2x2", ..., envir = parent.frame())
      # which fails when the package is only namespace-loaded (not attached).
      p_mcnemar <- NA_real_
      mpor <- NA_real_
      mpor_ci_lo <- NA_real_
      mpor_ci_hi <- NA_real_
      tryCatch({
        mcnemar_result <- exact2x2::exact2x2(
          contingency_table,
          alternative = "two.sided",
          tsmethod = "central",
          paired = TRUE,
          conf.level = conf_level
        )
        p_mcnemar <- mcnemar_result$p.value
        mpor <- unname(mcnemar_result$estimate)
        mpor_ci_lo <- mcnemar_result$conf.int[1]
        mpor_ci_hi <- mcnemar_result$conf.int[2]
      }, error = function(e) {
        warning("exact2x2::exact2x2 (mcnemar) failed for ", cy, ": ", e$message,
                "; returning NA for p-value and MPOR.", call. = FALSE)
      })

      delta_ci_lo <- NA_real_
      delta_ci_hi <- NA_real_
      tryCatch({
        # mcnemarExactDP(n, m, x): n=total, m=discordant pairs, x=gains (c).
        # Returns CI for (2x - m)/n = (c - b)/n = delta on proportion scale.
        dp_result <- exact2x2::mcnemarExactDP(n = n_pairs, m = n_discordant, x = c,
                                              conf.level = conf_level)
        delta_ci_lo <- 100 * dp_result$conf.int[1]
        delta_ci_hi <- 100 * dp_result$conf.int[2]
      }, error = function(e) {
        warning("exact2x2::mcnemarExactDP failed for ", cy, ": ", e$message,
                "; returning NA for delta CI.", call. = FALSE)
      })
    } else if (method == "midp") {
      # Mid-p: exact two-sided binomial p minus point mass at observed c.
      p_exact <- stats::binom.test(c, n_discordant, p = 0.5)$p.value
      p_mcnemar <- max(p_exact - stats::dbinom(c, n_discordant, 0.5), 0)
      # Haldane-corrected MPOR (always finite); no exact CI for midp -> NA.
      mpor <- (c + 0.5) / (b + 0.5)
      mpor_ci_lo <- NA_real_
      mpor_ci_hi <- NA_real_
      wald <- .paired_wald_delta_ci(b, c, n_pairs, conf_level)
      delta_ci_lo <- wald$ci_lo
      delta_ci_hi <- wald$ci_hi
    } else { # method == "noCC"
      # Asymptotic McNemar chi-squared without continuity correction.
      chi2 <- (b - c)^2 / (b + c)
      p_mcnemar <- stats::pchisq(chi2, df = 1, lower.tail = FALSE)
      mpor <- (c + 0.5) / (b + 0.5)
      mpor_ci_lo <- NA_real_
      mpor_ci_hi <- NA_real_
      wald <- .paired_wald_delta_ci(b, c, n_pairs, conf_level)
      delta_ci_lo <- wald$ci_lo
      delta_ci_hi <- wald$ci_hi
    }
    
    # 4) Transition rates with exact CIs
    rate_upward <- 100 * c / n_pairs
    bt_upward <- binom.test(c, n_pairs, conf.level = conf_level)
    rate_upward_ci <- 100 * bt_upward$conf.int
    
    rate_downward <- 100 * b / n_pairs
    bt_downward <- binom.test(b, n_pairs, conf.level = conf_level)
    rate_downward_ci <- 100 * bt_downward$conf.int
    
    # 5) Cohen's kappa for agreement
    p_observed <- (a + d) / n_pairs
    p_expected <- ((a + b) * (a + c) + (c + d) * (b + d)) / (n_pairs^2)
    cohen_kappa <- if (p_expected < 1) {
      (p_observed - p_expected) / (1 - p_expected)
    } else {
      1.0  # Perfect agreement when expected = 1
    }
    
    # Return results as a data frame
    data.frame(
      cytokine = cy,
      n_pairs = n_pairs,
      both_detect = a,
      loss = b,
      gain = c,
      neither_detect = d,
      n_discordant = n_discordant,
      n_baseline_detect = n_baseline_detect,
      n_comparison_detect = n_comparison_detect,
      prop_baseline = round(prop_baseline, 1),
      prop_baseline_ci_lo = round(prop_baseline_ci[1], 1),
      prop_baseline_ci_hi = round(prop_baseline_ci[2], 1),
      prop_comparison = round(prop_comparison, 1),
      prop_comparison_ci_lo = round(prop_comparison_ci[1], 1),
      prop_comparison_ci_hi = round(prop_comparison_ci[2], 1),
      delta_detection = round(delta_detection, 1),
      delta_ci_lo = round(delta_ci_lo, 1),
      delta_ci_hi = round(delta_ci_hi, 1),
      rate_upward = round(rate_upward, 1),
      rate_upward_ci_lo = round(rate_upward_ci[1], 1),
      rate_upward_ci_hi = round(rate_upward_ci[2], 1),
      rate_downward = round(rate_downward, 1),
      rate_downward_ci_lo = round(rate_downward_ci[1], 1),
      rate_downward_ci_hi = round(rate_downward_ci[2], 1),
      mpor = mpor,
      mpor_ci_lo = mpor_ci_lo,
      mpor_ci_hi = mpor_ci_hi,
      p_mcnemar = p_mcnemar,
      cohen_kappa = round(cohen_kappa, 3),
      stringsAsFactors = FALSE
    )
  })
  
  # Combine results
  mcnemar_results <- dplyr::bind_rows(mcnemar_results_list) %>%
    dplyr::mutate(
      # FDR correction
      q_mcnemar = p.adjust(p_mcnemar, method = fdr_method),
      # Detection pattern categorization based on delta_detection
      detection_pattern = dplyr::case_when(
        delta_detection > detection_change_threshold ~ "Increased detection",
        delta_detection < -detection_change_threshold ~ "Decreased detection",
        TRUE ~ "Stable detection"
      ),
      # Significance annotation
      mcnemar_significance = dplyr::case_when(
        q_mcnemar < 0.001 ~ "***",
        q_mcnemar < 0.01 ~ "**",
        q_mcnemar < 0.05 ~ "*",
        TRUE ~ ""
      )
    ) %>%
    dplyr::arrange(q_mcnemar, p_mcnemar)
  
  # Summary statistics
  if (!quiet) {
    cat("\nResults summary:\n")
    cat("  Total cytokines analyzed:", nrow(mcnemar_results), "\n")
    cat("  Significant changes (FDR < 0.05):", sum(mcnemar_results$q_mcnemar < 0.05, na.rm = TRUE), "\n")
    cat("  Increased detection:", sum(mcnemar_results$detection_pattern == "Increased detection", na.rm = TRUE), "\n")
    cat("  Decreased detection:", sum(mcnemar_results$detection_pattern == "Decreased detection", na.rm = TRUE), "\n")
    cat("  Stable detection:", sum(mcnemar_results$detection_pattern == "Stable detection", na.rm = TRUE), "\n")
    cat(sprintf("  Mean Cohen's kappa: %.3f\n", mean(mcnemar_results$cohen_kappa, na.rm = TRUE)))
    
    # Show top significant results with comprehensive metrics
    sig_results <- mcnemar_results %>%
      dplyr::filter(q_mcnemar < 0.05) %>%
      dplyr::slice_head(n = 5)
    
    if (nrow(sig_results) > 0) {
      cat("\n  Top significant changes (with 95% CIs):\n")
      for (i in 1:nrow(sig_results)) {
        row <- sig_results[i, ]
        cat(sprintf("    %s:\n", row$cytokine))
        cat(sprintf("      Delta = %+.1f%% [%.1f%%, %.1f%%]; MPOR = %.2f [%.2f, %.2f]\n",
                    row$delta_detection, row$delta_ci_lo, row$delta_ci_hi,
                    row$mpor, row$mpor_ci_lo, row$mpor_ci_hi))
        cat(sprintf("      Detection: %.1f%% -> %.1f%%; Gain=%d, Loss=%d; q=%.4f\n",
                    row$prop_baseline, row$prop_comparison, row$gain, row$loss, row$q_mcnemar))
      }
    }
    cat("\n")
  }

  class(mcnemar_results) <- c("mcnemar_detection", class(mcnemar_results))
  # Tag the return with the FDR method used so downstream consumers (and the
  # print method) can distinguish BH from bonferroni q-values, which would
  # otherwise look identical at the column-name level.
  attr(mcnemar_results, "fdr_method") <- fdr_method
  return(mcnemar_results)
}


#' Print method for McNemar detection results
#'
#' @param x Result from mcnemar_detection()
#' @param n_show Number of top results to show (default: 10)
#' @param ... Additional arguments (not used)
#'
#' @export
print.mcnemar_detection <- function(x, n_show = 10, ...) {
  if (!inherits(x, "data.frame")) {
    cat("McNemar Detection Results\n")
    cat("=========================\n\n")
    print(x)
    return(invisible(x))
  }
  
  cat("McNemar Detection Analysis Results\n")
  cat("==================================\n\n")
  cat("Total cytokines analyzed:", nrow(x), "\n")
  fdr_used <- attr(x, "fdr_method")
  if (is.null(fdr_used)) fdr_used <- "unknown"
  cat(sprintf("Significant changes (q < 0.05, %s-adjusted): %d\n\n",
              fdr_used,
              sum(x$q_mcnemar < 0.05, na.rm = TRUE)))
  
  cat("Detection patterns:\n")
  patterns <- table(x$detection_pattern)
  for (pattern in names(patterns)) {
    cat(sprintf("  %s: %d\n", pattern, patterns[pattern]))
  }
  
  cat("\nTop", min(n_show, nrow(x)), "results by significance:\n")
  top_results <- head(x, n_show)
  
  # Select key columns for display
  display_cols <- c("cytokine", "prop_baseline", "prop_comparison", 
                   "delta_detection", "delta_ci_lo", "delta_ci_hi",
                   "mpor", "p_mcnemar", "q_mcnemar", "mcnemar_significance")
  
  # Only include columns that exist
  display_cols <- intersect(display_cols, names(top_results))
  
  print(as.data.frame(top_results[, display_cols]))
  
  invisible(x)
}

