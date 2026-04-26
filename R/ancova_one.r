# Core single-analyte rank-based ANCOVA engine for immunoPlex
# Fits one rank-based ANCOVA with IQR outlier removal, conditional covariate
# inclusion, effect sizes via partial omega-squared, and optional log-scale validation

#' Fit a Single Rank-Based ANCOVA
#'
#' Fits one rank-based ANCOVA for a single analyte. Applies IQR outlier removal
#' on both outcome and covariate, rank-transforms outcome and primary covariate,
#' fits an OLS model, computes partial omega-squared effect sizes, and optionally
#' fits a log-scale validation model for concordance checking.
#'
#' @importFrom stats lm anova as.formula quantile summary.lm
#' @importFrom stats qt p.adjust coef confint residuals fitted
#'
#' @param data A data frame (one row per subject) containing the outcome, group,
#'   and covariate columns.
#' @param outcome Character. Column name for the outcome variable (ranked internally).
#' @param group Character. Column name for the group factor (2 or k levels).
#' @param covariate Character. Column name for the primary covariate (ranked internally).
#' @param covariates Character vector. Additional fixed covariates to include in
#'   the model (not rank-transformed). Default \code{NULL}.
#' @param conditional_covariates Character vector. Binary covariates to include only
#'   if their prevalence (proportion == 1) meets \code{prevalence_threshold}.
#'   Default \code{NULL}.
#' @param prevalence_threshold Numeric. Minimum prevalence (proportion == 1) for
#'   conditional covariates to be included. Default from config (0.03).
#' @param outlier_removal Logical. If \code{TRUE} (default), apply IQR outlier
#'   removal on both outcome and covariate columns.
#' @param iqr_multiplier Numeric. IQR multiplier for outlier fencing. Default
#'   from config (3).
#' @param min_n Integer. Minimum number of observations required after outlier
#'   removal. Returns a warning if below this threshold. Default from config (12).
#' @param conf_level Numeric. Confidence level for CIs. Default 0.95.
#' @param log_validate Logical. If \code{TRUE}, also fit a log-scale (untransformed)
#'   model for concordance checking. Default \code{FALSE}.
#' @param rep_col Optional character. Column name identifying technical
#'   replicates within a subject. When supplied, the fitter switches from
#'   \code{lm()} to \code{glmmTMB::glmmTMB()} (Gaussian) and appends a
#'   **nested** random intercept \code{(1 | subject_id:rep_col)}. The
#'   nested form is used because the same \code{rep_col} value across
#'   different subjects does not refer to the same replicate. For k-level
#'   groups the omnibus test switches from Type-I F to a likelihood-ratio
#'   test (\code{omnibus_f} is then \code{NA}; \code{omnibus_p} carries the
#'   LRT p-value). Default \code{NULL} (standard fixed-effects \code{lm()}
#'   fit).
#'   See \code{\link{aggregate_replicates}} for a pre-aggregation alternative.
#' @param plate_col Optional character. Column name identifying plates (a
#'   grouping that crosses subjects). When supplied, the fitter switches
#'   to \code{glmmTMB::glmmTMB()} and appends a **flat** random intercept
#'   \code{(1 | plate_col)}. Same k-level caveat as \code{rep_col}.
#'   Default \code{NULL}.
#'
#' @return An S3 object of class \code{"immuno_ancova"} (inherits \code{"immuno_model"}),
#'   a list containing:
#'   \describe{
#'     \item{model}{The fitted \code{lm} object (rank-based)}
#'     \item{log_model}{The fitted \code{lm} object on log scale, or \code{NULL}}
#'     \item{analyte}{Character or \code{NULL}}
#'     \item{group_var}{Character: name of the group variable}
#'     \item{n_groups}{Integer: number of group levels}
#'     \item{group_levels}{Character vector of group level names}
#'     \item{n_obs}{Integer: observations after outlier removal}
#'     \item{n_removed}{Integer: outliers removed}
#'     \item{formula_used}{Character: model formula as string}
#'     \item{covariates_included}{Character vector of covariates kept}
#'     \item{covariates_dropped}{Character vector of covariates dropped}
#'     \item{group_effect}{Numeric: group coefficient (binary models)}
#'     \item{group_se}{Numeric: SE of group coefficient (binary models)}
#'     \item{group_p}{Numeric: p-value for group effect (binary models)}
#'     \item{group_ci}{Numeric vector of length 2: CI for group effect (binary models)}
#'     \item{contrasts}{Data frame of pairwise contrasts (k-level models)}
#'     \item{omnibus_f}{Numeric: overall F statistic (k-level models)}
#'     \item{omnibus_p}{Numeric: overall F p-value (k-level models)}
#'     \item{r_squared}{Numeric: model R-squared}
#'     \item{omega_sq_partial}{Numeric: partial omega-squared for group}
#'     \item{omega_sq_ci}{Numeric vector of length 2: CI for partial omega-squared}
#'     \item{omega_magnitude}{Character: effect magnitude label}
#'     \item{log_effect}{Numeric or \code{NULL}: log-scale group coefficient}
#'     \item{log_fold_change}{Numeric or \code{NULL}: exp(log_effect)}
#'     \item{log_p}{Numeric or \code{NULL}: log-scale p-value}
#'     \item{warnings}{Character vector of warnings encountered}
#'   }
#'
#' @details
#' The rank-based ANCOVA approach:
#' \enumerate{
#'   \item Removes outliers using IQR fencing on both outcome and covariate
#'   \item Rank-transforms outcome and primary covariate (ties averaged)
#'   \item Fits OLS: \code{outcome_rank ~ group + covariate_rank + covariates}
#'   \item Computes partial omega-squared via \code{effectsize::omega_squared()}
#'   \item For k-level groups, extracts pairwise contrasts by releveling
#' }
#'
#' @examples
#' \dontrun{
#' result <- ancova_one(
#'   data = my_data,
#'   outcome = "cord_value",
#'   group = "covid_status",
#'   covariate = "maternal_value",
#'   covariates = c("age", "gest_delivery"),
#'   conditional_covariates = c("HIV", "malaria")
#' )
#' print(result)
#' summary(result)
#' }
#'
#' @export
ancova_one <- function(data,
                       outcome,
                       group,
                       covariate,
                       covariates = NULL,
                       conditional_covariates = NULL,
                       prevalence_threshold = getOption("immunoplex.ancova_prevalence_threshold", 0.03),
                       outlier_removal = TRUE,
                       iqr_multiplier = getOption("immunoplex.ancova_iqr_multiplier", 3),
                       min_n = getOption("immunoplex.ancova_min_n", 12L),
                       conf_level = 0.95,
                       log_validate = FALSE,
                       rep_col = NULL,
                       plate_col = NULL) {

  # --- Input validation ---
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame", call. = FALSE)
  }

  required_cols <- c(outcome, group, covariate)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  # Validate rep_col / plate_col up front so the error surfaces before the
  # downstream switch to the mixed-model fitter.
  if (!is.null(rep_col)) {
    if (!is.character(rep_col) || length(rep_col) != 1L) {
      stop("`rep_col` must be a single column name (character).", call. = FALSE)
    }
    if (!rep_col %in% names(data)) {
      stop("`rep_col` not found in data: '", rep_col, "'", call. = FALSE)
    }
  }
  if (!is.null(plate_col)) {
    if (!is.character(plate_col) || length(plate_col) != 1L) {
      stop("`plate_col` must be a single column name (character).", call. = FALSE)
    }
    if (!plate_col %in% names(data)) {
      stop("`plate_col` not found in data: '", plate_col, "'", call. = FALSE)
    }
  }
  use_re <- !is.null(rep_col) || !is.null(plate_col)
  if (use_re && !requireNamespace("glmmTMB", quietly = TRUE)) {
    stop("glmmTMB package required when rep_col or plate_col is supplied.",
         call. = FALSE)
  }

  # Drop rows with NA in outcome, group, or covariate

  data <- data[complete.cases(data[, c(outcome, group, covariate)]), , drop = FALSE]

  # Ensure group is a factor
  if (!is.factor(data[[group]])) {
    data[[group]] <- factor(data[[group]])
  }

  group_levels <- levels(data[[group]])
  n_groups <- length(group_levels)
  warn_msgs <- character(0)

  n_before <- nrow(data)

  # --- Outlier removal ---
  if (outlier_removal) {
    if (iqr_multiplier < 2 && n_before < 50) {
      warning(
        sprintf("iqr_multiplier=%g at n=%d is not recommended: B3 benchmarks show negative bias up to -0.112 and coverage below 0.93 for small samples with tight fences.",
                iqr_multiplier, n_before),
        call. = FALSE, immediate. = TRUE
      )
    }
    outlier_result <- .ancova_remove_outliers(data, outcome, covariate, iqr_multiplier)
    data <- outlier_result$data
    n_removed <- outlier_result$n_removed
  } else {
    n_removed <- 0L
  }

  n_obs <- nrow(data)

  # Check minimum sample size
  if (n_obs < min_n) {
    warn_msgs <- c(warn_msgs, sprintf(
      "Only %d observations after outlier removal (min_n = %d)", n_obs, min_n
    ))
  }

  # --- Conditional covariate inclusion ---
  covariates_included <- character(0)
  covariates_dropped <- character(0)

  # Add fixed covariates (check they exist in data)
  if (!is.null(covariates)) {
    present <- covariates[covariates %in% names(data)]
    absent <- setdiff(covariates, present)
    if (length(absent) > 0) {
      warn_msgs <- c(warn_msgs, paste("Covariates not in data:", paste(absent, collapse = ", ")))
    }
    covariates_included <- present
  }

  # Check conditional covariates for prevalence
  if (!is.null(conditional_covariates)) {
    for (cc in conditional_covariates) {
      if (cc %in% names(data)) {
        if (.ancova_check_prevalence(data, cc, prevalence_threshold)) {
          covariates_included <- c(covariates_included, cc)
        } else {
          covariates_dropped <- c(covariates_dropped, cc)
        }
      } else {
        covariates_dropped <- c(covariates_dropped, cc)
      }
    }
  }

  # --- Rank transformation ---
  data$outcome_rank <- rank(data[[outcome]], ties.method = "average")
  data$covariate_rank <- rank(data[[covariate]], ties.method = "average")

  # --- Build formula ---
  rhs_terms <- c(group, "covariate_rank", covariates_included)
  rhs_re <- character(0)
  if (!is.null(rep_col))   rhs_re <- c(rhs_re, paste0("(1 | subject_id:", rep_col, ")"))
  if (!is.null(plate_col)) rhs_re <- c(rhs_re, paste0("(1 | ", plate_col, ")"))
  formula_str <- paste("outcome_rank ~", paste(c(rhs_terms, rhs_re), collapse = " + "))
  model_formula <- as.formula(formula_str)

  # --- Fit rank-based ANCOVA ---
  # When rep_col/plate_col is supplied, switch to glmmTMB Gaussian so the
  # random-effect terms in the formula are honored; otherwise keep lm() so
  # rep_col = NULL reproduces the original behavior exactly.
  fit <- if (use_re) {
    tryCatch(
      suppressMessages(suppressWarnings(
        glmmTMB::glmmTMB(model_formula, data = data, family = stats::gaussian())
      )),
      error = function(e) {
        warn_msgs <<- c(warn_msgs, paste("Model fitting failed:", e$message))
        NULL
      }
    )
  } else {
    tryCatch(lm(model_formula, data = data), error = function(e) {
      warn_msgs <<- c(warn_msgs, paste("Model fitting failed:", e$message))
      NULL
    })
  }

  if (is.null(fit)) {
    return(.ancova_empty_result(
      analyte = NULL, group_var = group, n_groups = n_groups,
      group_levels = group_levels, n_obs = n_obs, n_removed = n_removed,
      formula_used = formula_str, covariates_included = covariates_included,
      covariates_dropped = covariates_dropped, warnings = warn_msgs
    ))
  }

  # --- R-squared ---
  # For glmmTMB, there is no single r.squared; leave it NA for the mixed case.
  fit_summary <- summary(fit)
  r_squared <- if (use_re) NA_real_ else fit_summary$r.squared

  # --- Effect sizes: partial omega-squared ---
  omega_result <- .ancova_compute_omega(fit, group, conf_level)

  # --- Extract group effects ---
  if (use_re) {
    cc_table <- fit_summary$coefficients$cond
    p_col <- "Pr(>|z|)"
    # Wald CIs using the standard normal critical value
    t_crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  } else {
    cc_table <- fit_summary$coefficients
    p_col <- "Pr(>|t|)"
    df_resid <- fit$df.residual
    t_crit <- qt(1 - (1 - conf_level) / 2, df = df_resid)
  }

  # Binary vs k-level
  group_effect <- group_se <- group_p <- NA_real_
  group_ci <- c(NA_real_, NA_real_)
  contrasts_df <- NULL
  omnibus_f <- omnibus_p <- NA_real_


  if (n_groups == 2) {
    # Binary model: extract coefficient for the non-reference level
    group_term <- paste0(group, group_levels[2])
    if (group_term %in% rownames(cc_table)) {
      group_effect <- as.numeric(cc_table[group_term, "Estimate"])
      group_se <- as.numeric(cc_table[group_term, "Std. Error"])
      group_p <- as.numeric(cc_table[group_term, p_col])
      group_ci <- c(group_effect - t_crit * group_se,
                     group_effect + t_crit * group_se)
    } else {
      warn_msgs <- c(warn_msgs, paste("Group term", group_term, "not found in model"))
    }
  } else {
    # k-level: pairwise contrasts + omnibus test
    contrasts_df <- .ancova_pairwise_contrasts(fit, group, data, model_formula,
                                               conf_level, use_re = use_re)

    if (use_re) {
      # Omnibus via LRT: refit without the group term and compare via anova().
      tryCatch({
        null_terms <- c("covariate_rank", covariates_included)
        null_formula <- as.formula(paste(
          "outcome_rank ~",
          paste(c(null_terms, rhs_re), collapse = " + ")
        ))
        fit_null <- suppressMessages(suppressWarnings(
          glmmTMB::glmmTMB(null_formula, data = data, family = stats::gaussian())
        ))
        lrt <- stats::anova(fit_null, fit)
        p_col_lrt <- grep("^Pr", colnames(lrt), value = TRUE)[1]
        if (!is.na(p_col_lrt)) omnibus_p <- as.numeric(lrt[[p_col_lrt]][2])
      }, error = function(e) {
        warn_msgs <<- c(warn_msgs, paste("Omnibus LRT extraction failed:", e$message))
      })
    } else {
      # Omnibus F from anova() on the lm fit
      tryCatch({
        aov_result <- anova(fit)
        grp_row <- which(rownames(aov_result) == group)
        if (length(grp_row) > 0) {
          omnibus_f <- aov_result$`F value`[grp_row]
          omnibus_p <- aov_result$`Pr(>F)`[grp_row]
        }
      }, error = function(e) {
        warn_msgs <<- c(warn_msgs, paste("Omnibus F extraction failed:", e$message))
      })
    }
  }

  # --- Log-scale validation ---
  log_model <- NULL
  log_effect <- log_fold_change <- log_p <- NULL

  if (log_validate) {
    log_result <- .ancova_log_validate(data, outcome, group, covariate,
                                        covariates_included, group_levels, n_groups)
    log_model <- log_result$model
    log_effect <- log_result$effect
    log_fold_change <- log_result$fold_change
    log_p <- log_result$p_value
  }

  # --- Construct result ---
  result <- structure(list(
    model = fit,
    log_model = log_model,
    analyte = NULL,
    group_var = group,
    n_groups = n_groups,
    group_levels = group_levels,
    n_obs = n_obs,
    n_removed = n_removed,
    formula_used = formula_str,
    covariates_included = covariates_included,
    covariates_dropped = covariates_dropped,
    group_effect = group_effect,
    group_se = group_se,
    group_p = group_p,
    group_ci = group_ci,
    contrasts = contrasts_df,
    omnibus_f = omnibus_f,
    omnibus_p = omnibus_p,
    r_squared = r_squared,
    omega_sq_partial = omega_result$omega,
    omega_sq_ci = c(omega_result$ci_low, omega_result$ci_high),
    omega_magnitude = .ancova_effect_label(omega_result$omega),
    log_effect = log_effect,
    log_fold_change = log_fold_change,
    log_p = log_p,
    warnings = warn_msgs
  ), class = c("immuno_ancova", "immuno_model"))

  result
}


# ============================================================================
# Internal helpers (not exported)
# ============================================================================

#' IQR outlier removal on outcome and covariate
#' @noRd
.ancova_remove_outliers <- function(data, outcome_col, covariate_col, iqr_multiplier) {
  n_before <- nrow(data)

  # Outcome outliers
  vals <- data[[outcome_col]]
  q1 <- quantile(vals, 0.25, na.rm = TRUE)
  q3 <- quantile(vals, 0.75, na.rm = TRUE)
  iqr_val <- q3 - q1
  keep_outcome <- vals >= (q1 - iqr_multiplier * iqr_val) &
                  vals <= (q3 + iqr_multiplier * iqr_val)

  # Covariate outliers
  cov_vals <- data[[covariate_col]]
  q1c <- quantile(cov_vals, 0.25, na.rm = TRUE)
  q3c <- quantile(cov_vals, 0.75, na.rm = TRUE)
  iqr_cov <- q3c - q1c
  keep_cov <- cov_vals >= (q1c - iqr_multiplier * iqr_cov) &
              cov_vals <= (q3c + iqr_multiplier * iqr_cov)

  data <- data[keep_outcome & keep_cov, , drop = FALSE]
  list(data = data, n_removed = n_before - nrow(data))
}


#' Check prevalence of a binary covariate
#' @noRd
.ancova_check_prevalence <- function(data, covariate_col, threshold) {
  mean(data[[covariate_col]] == 1, na.rm = TRUE) >= threshold
}


#' Compute effect magnitude label from omega-squared
#' @noRd
.ancova_effect_label <- function(omega_sq) {
  if (is.na(omega_sq)) return("Unknown")
  # Floor negative values at 0 (bias correction can produce negatives)
  omega_sq <- max(omega_sq, 0)
  if (omega_sq >= 0.14) return("Large")
  if (omega_sq >= 0.06) return("Medium")
  if (omega_sq >= 0.01) return("Small")
  "Negligible"
}


#' Compute partial omega-squared for the group variable
#' @noRd
.ancova_compute_omega <- function(fit, group_var, conf_level) {
  tryCatch({
    if (!requireNamespace("effectsize", quietly = TRUE)) {
      stop("Package 'effectsize' is required for omega_squared()", call. = FALSE)
    }
    es <- effectsize::omega_squared(fit, partial = TRUE, ci = conf_level,
                                     alternative = "two.sided")
    grp_row <- which(es$Parameter == group_var)
    if (length(grp_row) > 0) {
      list(
        omega = es$Omega2_partial[grp_row],
        ci_low = es$CI_low[grp_row],
        ci_high = es$CI_high[grp_row]
      )
    } else {
      list(omega = NA_real_, ci_low = NA_real_, ci_high = NA_real_)
    }
  }, error = function(e) {
    list(omega = NA_real_, ci_low = NA_real_, ci_high = NA_real_)
  })
}


#' Pairwise contrasts for k-level groups by releveling and refitting
#' @noRd
.ancova_pairwise_contrasts <- function(fit, group_var, data, model_formula,
                                       conf_level, use_re = FALSE) {
  group_levels <- levels(data[[group_var]])
  pairs <- utils::combn(group_levels, 2, simplify = FALSE)

  if (use_re) {
    # Wald CIs / z-based p-values from glmmTMB summary
    t_crit <- stats::qnorm(1 - (1 - conf_level) / 2)
    p_col  <- "Pr(>|z|)"
    refitter <- function(f, d) {
      suppressMessages(suppressWarnings(
        glmmTMB::glmmTMB(f, data = d, family = stats::gaussian())
      ))
    }
    coef_table <- function(m) summary(m)$coefficients$cond
  } else {
    t_crit <- qt(1 - (1 - conf_level) / 2, df = fit$df.residual)
    p_col  <- "Pr(>|t|)"
    refitter <- function(f, d) lm(f, data = d)
    coef_table <- function(m) summary(m)$coefficients
  }

  contrast_list <- lapply(pairs, function(pair) {
    # Relevel to use first element as reference
    data[[group_var]] <- relevel(factor(data[[group_var]], levels = group_levels),
                                 ref = pair[1])
    refit <- tryCatch(refitter(model_formula, data), error = function(e) NULL)
    if (is.null(refit)) {
      return(data.frame(
        term = paste(pair[2], "vs", pair[1]),
        estimate = NA_real_, se = NA_real_, p = NA_real_,
        ci_lower = NA_real_, ci_upper = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    cc_table <- coef_table(refit)
    term_name <- paste0(group_var, pair[2])

    if (term_name %in% rownames(cc_table)) {
      est <- as.numeric(cc_table[term_name, "Estimate"])
      se <- as.numeric(cc_table[term_name, "Std. Error"])
      p_val <- as.numeric(cc_table[term_name, p_col])
      data.frame(
        term = paste(pair[2], "vs", pair[1]),
        estimate = est, se = se, p = p_val,
        ci_lower = est - t_crit * se,
        ci_upper = est + t_crit * se,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        term = paste(pair[2], "vs", pair[1]),
        estimate = NA_real_, se = NA_real_, p = NA_real_,
        ci_lower = NA_real_, ci_upper = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  })

  do.call(rbind, contrast_list)
}


#' Fit log-scale validation model (on untransformed values)
#' @noRd
.ancova_log_validate <- function(data, outcome, group, covariate,
                                  covariates_included, group_levels, n_groups) {
  # Build formula using raw (log-scale) values, not ranked
  rhs_terms <- c(group, covariate, covariates_included)
  log_formula <- as.formula(paste(outcome, "~", paste(rhs_terms, collapse = " + ")))

  log_fit <- tryCatch(lm(log_formula, data = data), error = function(e) NULL)

  if (is.null(log_fit)) {
    return(list(model = NULL, effect = NULL, fold_change = NULL, p_value = NULL))
  }

  cc_table <- summary(log_fit)$coefficients
  log_effect <- log_p <- NULL
  log_fold_change <- NULL

  if (n_groups == 2) {
    term_name <- paste0(group, group_levels[2])
    if (term_name %in% rownames(cc_table)) {
      log_effect <- as.numeric(cc_table[term_name, "Estimate"])
      log_p <- as.numeric(cc_table[term_name, "Pr(>|t|)"])
      log_fold_change <- exp(log_effect)
    }
  }

  list(model = log_fit, effect = log_effect, fold_change = log_fold_change, p_value = log_p)
}


#' Construct empty result for failed models
#' @noRd
.ancova_empty_result <- function(analyte, group_var, n_groups, group_levels,
                                  n_obs, n_removed, formula_used,
                                  covariates_included, covariates_dropped, warnings) {
  structure(list(
    model = NULL,
    log_model = NULL,
    analyte = analyte,
    group_var = group_var,
    n_groups = n_groups,
    group_levels = group_levels,
    n_obs = n_obs,
    n_removed = n_removed,
    formula_used = formula_used,
    covariates_included = covariates_included,
    covariates_dropped = covariates_dropped,
    group_effect = NA_real_,
    group_se = NA_real_,
    group_p = NA_real_,
    group_ci = c(NA_real_, NA_real_),
    contrasts = NULL,
    omnibus_f = NA_real_,
    omnibus_p = NA_real_,
    r_squared = NA_real_,
    omega_sq_partial = NA_real_,
    omega_sq_ci = c(NA_real_, NA_real_),
    omega_magnitude = "Unknown",
    log_effect = NULL,
    log_fold_change = NULL,
    log_p = NULL,
    warnings = warnings
  ), class = c("immuno_ancova", "immuno_model"))
}
