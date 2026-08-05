# LOD handling comparison function for immunoPlex
# Compares preprocessing strategies versus censoring-aware models for limit-of-detection data

#' Compare LOD-handling model families
#'
#' @description
#' High-level helper that compares different approaches to handling 
#' limit-of-detection (LOD) values: preprocessing strategies vs. 
#' censoring-aware models.
#'
#' @param dat Data frame from prepare_cytokine_data() or similar, with
#'   columns `value`, `cens_lod` and `lod`. Unlike \code{\link{fit_one}},
#'   this function expects `value` on the **raw concentration scale**,
#'   because the LOD-substitution strategies it tests are defined in raw
#'   units (e.g. LOD/2). It log-transforms internally where the fitted
#'   family requires it, so pass raw values and do not pre-log.
#' @param families Character vector of model families to compare.
#'   Options: "gamma", "tobit", "tobit_censreg", "aft"
#' @param lod_methods Character vector of LOD preprocessing strategies to test
#'   with gamma models. Options: "half", "zero", "sqrt", "lod", "halfmin"
#' @param include_log_transform Logical. Whether to include log-transformed
#'   versions of preprocessing methods.
#' @param fixed Fixed effects formula string
#' @param random Random effects formula string
#' @param ... Additional arguments passed to fit_one()
#'
#' @return An object of class `immuno_lod_comparison` containing:
#'   - Preprocessed model results (gamma with different LOD methods)
#'   - Censoring-aware model results (tobit, AFT)
#'   - Overall comparison table
#'   - Best model recommendation
#'
#' @examples
#' \dontrun{
#' # Load and prepare data
#' data("immunoplex_example", package = "immunoPlex")
#' dat <- prepare_cytokine_data("Cytokine_01")
#' 
#' # Compare LOD handling approaches
#' comparison <- compare_lod_models(
#'   dat, 
#'   families = c("gamma", "tobit", "tobit_censreg"),
#'   lod_methods = c("half", "halfmin"),
#'   include_log_transform = TRUE
#' )
#' 
#' print(comparison)
#' plot(comparison)
#' }
#' @export
compare_lod_models <- function(dat,
                               families = c("gamma", "tobit", "aft"),
                               lod_methods = c("half", "halfmin"),
                               include_log_transform = TRUE,
                               fixed = "timepoint*disease + age",
                               random = "(1|subject_id)",
                               ...) {
  
  # Validate inputs
  stopifnot(is.data.frame(dat))
  stopifnot(all(c("value", "cens_lod", "lod") %in% names(dat)))
  
  results <- list()
  all_models <- list()
  
  # 1. Test gamma models with different LOD preprocessing methods
  if ("gamma" %in% families) {
    cat("=== Testing LOD Preprocessing Strategies ===\n")
    
    for (method in lod_methods) {
      cat("Testing gamma with LOD method:", method, "\n")
      
      # Apply LOD preprocessing
      dat_processed <- .apply_lod_method(dat, method)
      
      # Fit gamma model
      model_name <- paste0("gamma_", method)
      model <- tryCatch({
        fit_one(dat_processed, family = "gamma", fixed = fixed, 
                random = random, ...)
      }, error = function(e) NULL)
      
      if (!is.null(model) && model$converged) {
        all_models[[model_name]] <- model
        cat("[OK]", model_name, "AIC =", round(model$aic, 2), "\n")
      }
      
      # Test log-transformed version if requested
      if (include_log_transform) {
        log_model_name <- paste0("gamma_log_", method)
        dat_log <- dat_processed
        dat_log$value <- log(pmax(dat_processed$value, min(dat_processed$lod)/2))
        
        model_log <- tryCatch({
          fit_one(dat_log, family = "gamma", fixed = fixed, 
                  random = random, ...)
        }, error = function(e) NULL)
        
        if (!is.null(model_log) && model_log$converged) {
          all_models[[log_model_name]] <- model_log
          cat("[OK]", log_model_name, "AIC =", round(model_log$aic, 2), "\n")
        }
      }
    }
  }
  
  # 2. Test censoring-aware models
  cat("\n=== Testing Censoring-Aware Models ===\n")
  
  # `dat$value` arrives on the raw concentration scale (that is what the
  # LOD-substitution branch above requires), but fit_one() expects `value`
  # on the log scale for every censoring-aware family: it builds Tobit
  # bounds as log(LOD), and AFT exponentiates `value` back to raw. Convert
  # here so both branches are fed the scale they document.
  dat_cens <- dat
  dat_cens$value <- log(pmax(dat$value, min(dat$lod, na.rm = TRUE) / 2))

  censoring_families <- setdiff(families, "gamma")
  for (fam in censoring_families) {
    cat("Testing", fam, "model...\n")

    model <- tryCatch({
      fit_one(dat_cens, family = fam, fixed = fixed, random = random, ...)
    }, error = function(e) NULL)
    
    if (!is.null(model) && model$converged) {
      all_models[[fam]] <- model
      cat("[OK]", fam, "AIC =", round(model$aic, 2), "\n")
    }
  }
  
  if (length(all_models) == 0) {
    stop("No models converged successfully")
  }
  
  # 3. Create comprehensive comparison
  comparison <- data.frame(
    model = names(all_models),
    family = sapply(all_models, function(x) x$family),
    lod_method = sapply(names(all_models), function(x) {
      if (grepl("gamma_", x)) {
        parts <- strsplit(x, "_")[[1]]
        if (length(parts) > 2) paste(parts[2:length(parts)], collapse = "_")
        else parts[2]
      } else "native_censoring"
    }),
    aic = sapply(all_models, function(x) x$aic),
    bic = sapply(all_models, function(x) x$bic),
    loglik = sapply(all_models, function(x) x$logLik),
    stringsAsFactors = FALSE
  )
  comparison$delta_aic <- comparison$aic - min(comparison$aic)
  comparison <- comparison[order(comparison$aic), ]
  
  # Best model
  best_name <- comparison$model[1]
  best_model <- all_models[[best_name]]
  
  # Create result object
  result <- list(
    models = all_models,
    comparison = comparison,
    best_model = best_model,
    best_name = best_name,
    lod_methods_tested = lod_methods,
    families_tested = families,
    include_log_transform = include_log_transform,
    call = match.call()
  )
  
  class(result) <- "immuno_lod_comparison"
  return(result)
}

# Helper function to apply LOD preprocessing methods
.apply_lod_method <- function(dat, method) {
  dat_out <- dat
  lod_val <- dat$lod
  
  # Apply the specified LOD method to censored values
  censored_idx <- which(dat$cens_lod)
  
  if (length(censored_idx) > 0) {
    replacement_values <- switch(method,
      "half" = lod_val[censored_idx] / 2,
      "zero" = rep(0, length(censored_idx)),
      "sqrt" = sqrt(lod_val[censored_idx]),
      "lod" = lod_val[censored_idx],
      "halfmin" = {
        min_pos <- min(dat$value[dat$value > 0 & !dat$cens_lod], na.rm = TRUE)
        rep(min_pos / 2, length(censored_idx))
      },
      stop("Unsupported LOD method: ", method)
    )
    
    dat_out$value[censored_idx] <- replacement_values
  }
  
  return(dat_out)
}
