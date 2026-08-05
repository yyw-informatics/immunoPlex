# Model comparison wrapper function for immunoPlex
# Fits multiple model families and provides comprehensive comparison analysis

#' Fit a set of candidate models
#'
#' @description
#' Wrapper that iterates over multiple families and returns an
#' `immuno_model_set` for model comparison and selection.
#' @importFrom utils capture.output
#'
#' @inheritParams fit_one
#' @param families Character vector of families to evaluate.
#'   Options: "gamma", "tobit", "tobit_censreg", "aft", "gaussian"
#' @param return_best Logical. If TRUE, returns only the best model by AIC.
#'   If FALSE, returns all fitted models.
#' @param quiet Logical. If TRUE, suppresses progress messages and warnings.
#'   Useful for automated testing or batch processing.
#' @param debug Logical. If TRUE, provides additional debugging information
#'   when models fail to fit, including data summaries and censoring details.
#'
#' @return An object of class `immuno_model_set` containing fitted models,
#'   comparison statistics, and the best model. The comparison table includes
#'   an 'estimand' column indicating whether each family targets ratio_of_means
#'   (gamma, tobit, tobit_censreg) or ratio_of_medians (aft).
#'
#' @details
#' **AIC comparability across families.** The comparison table's `aic` and
#' `delta_aic` columns rank families by a single information criterion, but
#' the underlying likelihoods are not on a common response scale: Gamma is
#' fit raw-scale with a log link, Tobit is fit on `log(value)` (Gaussian
#' log-scale), and AFT uses a lognormal internally. Cross-family AIC
#' comparisons across a candidate set that spans scales are therefore not
#' strictly valid, and "best by AIC" should be interpreted with caution
#' whenever the candidate set mixes raw-scale and log-scale likelihoods.
#' When this matters operationally, restrict `families` to a single scale
#' (e.g. `c("tobit", "aft")` or `c("gamma")` alone).
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data("immunoplex_example", package = "immunoPlex")
#' 
#' # Prepare data for one cytokine
#' dat <- prepare_cytokine_data("Cytokine_01") 
#' 
#' # Fit multiple models
#' models <- fit_models(dat, families = c("gamma", "tobit", "aft"))
#' print(models)
#' 
#' # View comparison table with estimand information
#' models$comparison
#' # Shows: family, estimand, aic, bic, loglik, n_cens_lod, n_cens_ulod, delta_aic
#' 
#' summary(models)
#' plot(models)
#' }
#' @export
fit_models <- function(dat,
                       families = c("gamma", "tobit", "aft"),
                       fixed = "timepoint*disease + age",
                       random = "(1|subject_id)",
                       ulod = FALSE,
                       dispformula = ~1,
                       rep_col = NULL,
                       plate_col = NULL,
                       return_best = FALSE,
                       quiet = FALSE,
                       debug = FALSE,
                       ...) {

  # Validate inputs - allow non-dataframes to proceed to get detailed error messages
  stopifnot(length(families) > 0)
  supported_families <- c("gamma", "tobit", "tobit_censreg", "aft", "gaussian")
  families <- match.arg(families, supported_families, several.ok = TRUE)
  
  # Fit models for each family
  models <- list()
  converged_models <- character()
  failures <- character()

  for (fam in families) {
    if (!quiet) cat("Fitting", fam, "model...\n")

    model <- tryCatch(
      fit_one(dat, family = fam, fixed = fixed, random = random, ulod = ulod,
              dispformula = dispformula, rep_col = rep_col, plate_col = plate_col,
              ...),
      fit_one_missing_cols = function(e) {
        if (!quiet) {
          missing_cols <- setdiff(c("value", "cens_lod", "cens_ulod"), names(dat))
          available_cols <- paste(names(dat), collapse = ", ")
          message("Failed to fit ", fam, " model: Missing required columns")
          message("  -> Missing: ", paste(missing_cols, collapse = ", "))
          message("  -> Available: ", available_cols)
          message("  -> Fix: Ensure data has 'value', 'cens_lod', 'cens_ulod' columns")
        }
        failures[[fam]] <<- conditionMessage(e)
        NULL
      },
      fit_one_not_dataframe = function(e) {
        if (!quiet) {
          message("Failed to fit ", fam, " model: Input must be a data.frame")
          message("  -> Issue: Input is not a data.frame (got: ", class(dat)[1], ")")
          message("  -> Fix: Pass a valid data.frame to fit_models()")
        }
        failures[[fam]] <<- conditionMessage(e)
        NULL
      },
      fit_one_formula_vars = function(e) {
        if (!quiet) {
          message("Failed to fit ", fam, " model: Variable not found")
          message("  -> Issue: ", conditionMessage(e))
          message("  -> Fix: Check that all variables in formula exist in data")
        }
        failures[[fam]] <<- conditionMessage(e)
        NULL
      },
      fit_one_missing_package = function(e) {
        if (!quiet) {
          message("Failed to fit ", fam, " model: Missing package dependency")
          message("  -> Issue: ", conditionMessage(e))
          message("  -> Fix: Install required package")
        }
        failures[[fam]] <<- conditionMessage(e)
        NULL
      },
      fit_one_missing_lod = function(e) {
        if (!quiet) {
          message("Failed to fit ", fam, " model: LOD data incomplete")
          message("  -> Issue: Left-censored observations missing LOD values")
          message("  -> Fix: Ensure 'lod' column has values where cens_lod = TRUE")
        }
        failures[[fam]] <<- conditionMessage(e)
        NULL
      },
      fit_one_missing_ulod = function(e) {
        if (!quiet) {
          message("Failed to fit ", fam, " model: ULOD data incomplete")
          message("  -> Issue: Right-censored observations missing ULOD values")
          message("  -> Fix: Ensure 'ulod' column has values where cens_ulod = TRUE")
        }
        failures[[fam]] <<- conditionMessage(e)
        NULL
      },
      fit_one_censreg_row_lod = function(e) {
        if (!quiet) {
          message("Failed to fit ", fam, " model: censReg incompatible with row-specific LODs")
          message("  -> Issue: censReg package requires a single LOD value")
          message("  -> Fix: Use 'tobit' family instead for row-specific LODs")
        }
        failures[[fam]] <<- conditionMessage(e)
        NULL
      },
      fit_one_invalid_bounds = function(e) {
        if (!quiet) {
          message("Failed to fit ", fam, " model: Invalid censoring intervals")
          message("  -> Issue: Left bound > right bound in censoring intervals")
          message("  -> Fix: Check LOD/ULOD values and censoring flags")
        }
        failures[[fam]] <<- conditionMessage(e)
        NULL
      },
      error = function(e) {
        if (!quiet) {
          message("Failed to fit ", fam, " model: ", e$message)
          # Safe dimension reporting for any input type
          if (is.data.frame(dat)) {
            if (debug) {
              message("Data dimensions: ", nrow(dat), " rows")
            } else {
              message("  -> Data dimensions: ", nrow(dat), " rows x ", ncol(dat), " cols")
            }
          } else {
            if (debug) {
              message("Data dimensions: input is not a data.frame (", class(dat)[1], ")")
            } else {
              message("  -> Data dimensions: input is not a data.frame (", class(dat)[1], ")")
            }
          }
          message("  -> Formula: ", fixed)
          if (nzchar(random)) message("  -> Random: ", random)
          if (debug && is.data.frame(dat)) {
            message("  -> Censoring: ", sum(dat$cens_lod, na.rm = TRUE), " left, ",
                   sum(dat$cens_ulod, na.rm = TRUE), " right")
            message("  -> Data summary: ", paste(capture.output(summary(dat)), collapse = "; "))
          }
        }
        failures[[fam]] <<- conditionMessage(e)
        NULL
      }
    )

    if (!is.null(model) && model$converged) {
      models[[fam]] <- model
      converged_models <- c(converged_models, fam)
      if (!quiet) cat("[OK]", fam, "converged, AIC =", round(model$aic, 2), "\n")
    } else {
      if (!quiet) cat("[FAIL]", fam, "failed to converge\n")
      # Show debug information for non-converged models when debug=TRUE
      if (debug && !is.null(model) && !quiet) {
        if (is.data.frame(dat)) {
          message("Data dimensions: ", nrow(dat), " rows")
        } else {
          message("Data dimensions: input is not a data.frame (", class(dat)[1], ")")
        }
      }
    }
  }

  if (length(models) == 0) {
    cond <- structure(
      list(
        message = "No models converged successfully",
        call    = sys.call(),
        failures = failures
      ),
      class = c("fit_models_all_failed", "error", "condition")
    )
    stop(cond)
  }
  
  # Create model comparison table
  comparison <- data.frame(
    family = names(models),
    estimand = sapply(models, function(x) x$estimand),
    aic = sapply(models, function(x) x$aic),
    bic = sapply(models, function(x) x$bic),
    loglik = sapply(models, function(x) x$logLik),
    n_cens_lod = sapply(models, function(x) x$n_cens_lod),
    n_cens_ulod = sapply(models, function(x) x$n_cens_ulod),
    stringsAsFactors = FALSE
  )
  comparison$delta_aic <- comparison$aic - min(comparison$aic)
  comparison <- comparison[order(comparison$aic), ]

  if (!quiet && length(unique(comparison$estimand)) > 1) {
    warning("fit_models comparison mixes estimands (",
            paste(unique(comparison$estimand), collapse = ", "),
            "); the AIC-best model may not be estimand-stable.",
            call. = FALSE)
  }

  # Identify best model
  best_family <- comparison$family[1]
  best_model <- models[[best_family]]
  
  # Create result object
  result <- list(
    models = models,
    comparison = comparison,
    best_model = best_model,
    best_family = best_family,
    call = match.call(),
    n_models = length(models),
    option_state = list(
      min_subjects = getOption("fit_one.min_subjects"),
      min_reps     = getOption("fit_one.min_reps")
    ),
    failures = failures
  )
  
  class(result) <- "immuno_model_set"
  
  if (return_best) {
    attr(best_model, "comparison") <- comparison
    return(best_model)
  } else {
    return(result)
  }
}
