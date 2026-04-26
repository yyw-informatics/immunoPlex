# Statistical Test Data Generation Framework
#
# Provides controlled data generation for statistical validation:
# - Synthetic data with defined statistical properties
# - Systematic data source management
# - Test environment configuration
# - Model convergence validation
#
# Core functions:
# - create_robust_test_data(): Generates data with specified statistical properties
# - setup_robust_test_data(): Implements data source selection protocol
# - test_convergence_likelihood(): Validates model estimation conditions
#
# Test environment configuration
if (requireNamespace("grDevices", quietly = TRUE)) {
  # Suppress graphics output during testing
  old_dev <- grDevices::dev.cur()
  if (old_dev > 1) {
    grDevices::dev.off()
  }
  # Use null device to prevent file generation
  grDevices::pdf(file = NULL)
}

#' Generate synthetic test data with controlled statistical properties
#' 
#' Creates synthetic data with specified distributional characteristics
#' and covariate structures for statistical model validation.
#' 
#' @param n_subjects Number of subjects (default: 20)
#' @param n_reps Number of replicates per subject (default: 3) 
#' @param censoring_rate Proportion of values to be censored (default: 0.3)
#' @param family_type Distribution type: "gamma", "tobit", or "normal" (default: "gamma")
#' @param seed Random seed for reproducibility (default: 123)
#' @param add_covariates Include structured covariate effects (default: TRUE)
#' 
#' @return Data frame with defined statistical properties for model validation
create_robust_test_data <- function(n_subjects = 20, 
                                   n_reps = 3, 
                                   censoring_rate = 0.3, 
                                   family_type = "gamma", 
                                   seed = 123,
                                   add_covariates = TRUE) {
  
  set.seed(seed)
  total_n <- n_subjects * n_reps
  
  # Generate observations from specified distribution
  if (family_type == "gamma") {
    # Gamma distribution with controlled shape parameter
    base_values <- rgamma(total_n, shape = 2, rate = 1)
  } else if (family_type == "tobit") {
    # Log-normal distribution with defined censoring characteristics
    base_values <- exp(rnorm(total_n, mean = 1, sd = 0.8))
  } else {
    # Log-normal distribution for accelerated failure time models
    base_values <- exp(rnorm(total_n, mean = 1.5, sd = 0.7))
  }
  
  # Incorporate structured covariate effects
  if (add_covariates) {
    # Define systematic covariate effects with controlled magnitudes
    timepoint_effect <- rep(c(0, 0.2), length.out = total_n)  # Temporal effect
    disease_effect <- rep(c(0, 0.15), each = total_n/2)       # Condition effect
    age_effect <- rep(as.numeric(scale(25:(25+n_subjects-1))[,1]), each = n_reps) * 0.1  # Age gradient
    
    base_values <- base_values * exp(timepoint_effect + disease_effect + age_effect)
  }
  
  # Set LOD to achieve desired censoring rate
  lod_value <- quantile(base_values, censoring_rate)
  ulod_value <- quantile(base_values, 0.95)  # 5% upper censoring
  
  # Create data frame with proper structure (avoid row name warnings)
  # Build the data frame all at once to avoid row name propagation issues
  if (add_covariates) {
    dat <- data.frame(
      subject_id = factor(rep(1:n_subjects, each = n_reps)),
      value = as.numeric(base_values),
      timepoint = factor(rep(c("Pre", "Post"), length.out = total_n)),
      disease = factor(rep(c("Healthy", "Disease"), each = total_n/2)),
      age = as.numeric(rep(25:(25+n_subjects-1), each = n_reps)),
      lod = lod_value,
      ulod = ulod_value,
      stringsAsFactors = FALSE,
      row.names = NULL  # Explicitly avoid row name warnings
    )
  } else {
    dat <- data.frame(
      subject_id = factor(rep(1:n_subjects, each = n_reps)),
      value = as.numeric(base_values),
      lod = lod_value,
      ulod = ulod_value,
      stringsAsFactors = FALSE,
      row.names = NULL  # Explicitly avoid row name warnings
    )
  }
  
  # Set censoring flags
  dat$cens_lod <- dat$value < dat$lod
  dat$cens_ulod <- dat$value > dat$ulod
  
  # Apply constrained random imputation for censored values
  dat$value[dat$cens_lod] <- runif(sum(dat$cens_lod), 
                                   min = lod_value * 0.1, 
                                   max = lod_value * 0.9)
  
  return(dat)
}

#' Data source selection protocol implementation
#' 
#' Implements a systematic protocol for test data selection,
#' prioritizing empirical data with defined fallback mechanisms.
#' 
#' @param prefer_real Prioritize empirical data selection (default: TRUE)
#' @param family_hint Distribution specification for synthetic data generation
setup_robust_test_data <- function(prefer_real = TRUE, family_hint = "gamma") {
  
  if (prefer_real) {
    # Try real data first
    tryCatch({
      data("immunoplex_example", package = "immunoPlex", envir = environment())
      
      raw_expr <- immunoplex_example$expression
      meta <- immunoplex_example$metadata  
      lod_lookup <- immunoplex_example$lod_lookup
      
      # Identify cytokines with complete LOD information
      valid_cytos <- which(!is.na(lod_lookup$lod))
      if (length(valid_cytos) == 0) stop("No cytokines with LOD information found")
      
      target_cyto <- lod_lookup$cytokine[valid_cytos[1]]
      target_lod <- lod_lookup$lod[lod_lookup$cytokine == target_cyto]
      target_ulod <- lod_lookup$ulod[lod_lookup$cytokine == target_cyto]
      
      dat <- data.frame(
        subject_id = factor(meta$subject_id[match(rownames(raw_expr), meta$sample_id)]),
        value = as.numeric(raw_expr[, target_cyto]),
        timepoint = factor(meta$timepoint[match(rownames(raw_expr), meta$sample_id)]),
        disease = factor(meta$disease[match(rownames(raw_expr), meta$sample_id)]),
        age = as.numeric(meta$age[match(rownames(raw_expr), meta$sample_id)]),
        lod = target_lod,
        ulod = if(is.na(target_ulod)) max(raw_expr[, target_cyto], na.rm = TRUE) * 2 else target_ulod,
        stringsAsFactors = FALSE,
        row.names = NULL  # Explicitly avoid row name warnings
      )
      
      dat$cens_lod <- dat$value < dat$lod
      dat$cens_ulod <- !is.na(dat$ulod) & dat$value > dat$ulod
      
      return(dat)
      
    }, error = function(e) {
      # Fall back to synthetic data if real data fails
      NULL
    })
  }
  
  # Use synthetic data (either by choice or as fallback)
  return(create_robust_test_data(family_type = family_hint))
}

#' Model estimation condition validation
#' 
#' Evaluates dataset compatibility with specified statistical models
#' 
#' @param dat Data frame for validation
#' @param families Statistical models to evaluate (default: c("gamma", "tobit"))
#' @return Named logical vector indicating estimation feasibility per model
test_convergence_likelihood <- function(dat, families = c("gamma", "tobit")) {
  
  if (!requireNamespace("immunoPlex", quietly = TRUE)) {
    return(setNames(rep(FALSE, length(families)), families))
  }
  
  results <- setNames(rep(FALSE, length(families)), families)
  
  for (family in families) {
    tryCatch({
      fit <- immunoPlex::fit_one(dat, family = family, random = "")
      results[family] <- fit$converged
    }, error = function(e) {
      results[family] <- FALSE
    })
  }
  
  return(results)
}