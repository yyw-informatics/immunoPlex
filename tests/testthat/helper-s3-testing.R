# S3 Method Testing Framework for immunoPlex
#
# Provides comprehensive testing utilities for S3 methods:
# - Test data generation optimized for S3 testing
# - Method existence and functionality verification
# - Object creation and validation
# - Dependency checking for conditional tests
#
# Core functions:
# - create_s3_test_data(): Generates test data for S3 methods
# - test_s3_method(): Verifies S3 method functionality
# - create_all_immunoplex_objects(): Creates test objects
# - validate_s3_object(): Validates object structure

#' Create test data optimized for S3 method testing
#' 
#' @param family_hint Hint for model family to optimize for
#' @param n_obs Number of observations
#' @param censoring_rate Proportion of censored observations
#' @return Data frame suitable for testing S3 methods
create_s3_test_data <- function(family_hint = "gamma", n_obs = 50, censoring_rate = 0.3) {
  
  set.seed(42)  # Reproducible for tests
  
  # Generate well-behaved data for S3 testing
  if (family_hint == "gamma") {
    base_values <- rgamma(n_obs, shape = 2, rate = 0.5)
  } else if (family_hint == "tobit") {
    base_values <- exp(rnorm(n_obs, mean = 1, sd = 0.5))
  } else {
    base_values <- exp(rnorm(n_obs, mean = 1.5, sd = 0.6))
  }
  
  # Add systematic effects for realistic modeling
  timepoint_effect <- rep(c(0, 0.2), length.out = n_obs)
  disease_effect <- rep(c(0, 0.15), each = n_obs/2)
  age_effect <- scale(rep(25:74, length.out = n_obs))[,1] * 0.1
  
  base_values <- base_values * exp(timepoint_effect + disease_effect + age_effect)
  
  # Set LOD to achieve desired censoring
  lod_value <- quantile(base_values, censoring_rate)
  ulod_value <- quantile(base_values, 0.95)
  
  # Create comprehensive data frame
  dat <- data.frame(
    subject_id = factor(rep(1:(n_obs/2), each = 2)),
    value = as.numeric(base_values),
    timepoint = factor(rep(c("Pre", "Post"), length.out = n_obs)),
    disease = factor(rep(c("Healthy", "Disease"), each = n_obs/2)),
    age = as.numeric(rep(25:74, length.out = n_obs)),
    lod = lod_value,
    ulod = ulod_value,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  # Set censoring flags
  dat$cens_lod <- dat$value < dat$lod
  dat$cens_ulod <- dat$value > dat$ulod
  
  # Apply censoring to values
  dat$value[dat$cens_lod] <- runif(sum(dat$cens_lod), 
                                   min = lod_value * 0.5, 
                                   max = lod_value * 0.9)
  
  return(dat)
}

#' Test S3 method existence and basic functionality
#' 
#' @param object Object to test
#' @param method Method name (e.g., "print", "summary", "plot")
#' @param ... Additional arguments for the method
#' @return TRUE if method works, FALSE otherwise
test_s3_method <- function(object, method, ...) {
  # If method already contains a dot, use it as-is (e.g., "summary.immuno_model")
  if (grepl("\\.", method)) {
    method_name <- method
    base_method <- sub("\\..*", "", method)
  } else {
    method_name <- paste0(method, ".", class(object)[1])
    base_method <- method
  }
  
  if (!exists(method_name)) {
    return(FALSE)
  }
  
  tryCatch({
    result <- do.call(base_method, list(object, ...))
    return(TRUE)
  }, error = function(e) {
    return(FALSE)
  })
}

#' Create objects of all immunoPlex types for comprehensive testing
#' 
#' @return List with objects of different types
create_all_immunoplex_objects <- function() {
  
  # Create base data
  dat <- create_s3_test_data()
  
  objects <- list()
  
  # Try to create immuno_fit object
  tryCatch({
    fit <- fit_one(dat, family = "gamma", random = "")
    if (fit$converged) {
      objects$immuno_fit <- fit
    }
  }, error = function(e) NULL)
  
  # Try to create immuno_model_set object  
  tryCatch({
    models <- fit_models(dat, families = c("gamma", "tobit"), random = "")
    if (length(models$models) > 0) {
      objects$immuno_model_set <- models
    }
  }, error = function(e) NULL)
  
  # Try to create immuno_preprocess object
  tryCatch({
    if (requireNamespace("immunoPlex", quietly = TRUE) && 
        exists("immunoplex_example", envir = asNamespace("immunoPlex"))) {
      
      data("immunoplex_example", package = "immunoPlex", envir = environment())
      preprocessed <- immuno_preprocess(
        immunoplex_example$expression[1:20, 1:5],  # Smaller subset for testing
        immunoplex_example$metadata[1:20, ],
        immunoplex_example$lod_lookup[1:5, ]
      )
      objects$immuno_preprocess <- preprocessed
    }
  }, error = function(e) NULL)
  
  return(objects)
}

#' Validate that an S3 object has proper structure
#' 
#' @param object Object to validate
#' @param expected_class Expected class name
#' @param required_components Required list components
#' @return TRUE if valid, error message if not
validate_s3_object <- function(object, expected_class, required_components) {
  
  # Check class
  if (!inherits(object, expected_class)) {
    stop("Object should have class '", expected_class, "' but has '", 
         paste(class(object), collapse = ", "), "'")
  }
  
  # Check required components
  missing_components <- setdiff(required_components, names(object))
  if (length(missing_components) > 0) {
    stop("Missing required components: ", paste(missing_components, collapse = ", "))
  }
  
  return(TRUE)
}

#' Skip test if plotting packages are not available
skip_if_no_plotting <- function() {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
}

#' Skip test if statistical packages are not available  
skip_if_no_stats_packages <- function() {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("survival")
}