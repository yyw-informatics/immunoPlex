# Residual computation functions for plotting
# Extracts and computes appropriate residuals for different model types and families

#' Get residuals and fitted values for plotting
#' 
#' @param x An immuno_fit object
#' @param plot_data Data frame with plotting data
#' @return List with residuals, fitted values, and residual type
#' @keywords internal
.get_residuals_and_fitted <- function(x, plot_data) {
  
  if (x$family == "gamma" && inherits(x$model, "glmmTMB")) {
    resids <- residuals(x$model, type = "pearson")
    fitted_vals <- fitted(x$model)
    resid_type <- "Pearson"
    
  } else if (x$family %in% c("tobit", "aft") && inherits(x$model, "survreg")) {
    # Try randomized quantile residuals first for censored models
    if (requireNamespace("statmod", quietly = TRUE)) {
      tryCatch({
        resids <- statmod::qresiduals(x$model)
        resid_type <- "Randomized Quantile"
      }, error = function(e) {
        resids <<- residuals(x$model, type = "deviance")
        resid_type <<- "Deviance"
      })
    } else {
      resids <- residuals(x$model, type = "deviance")
      resid_type <- "Deviance"
    }
    fitted_vals <- predict(x$model)
    
  } else if (x$family == "tobit_censreg" && inherits(x$model, "censReg")) {
    # Handle censReg models - use manual fitted values since predict() may not work
    if (requireNamespace("censReg", quietly = TRUE)) {
      # Initialize variables to avoid "object not found" errors
      fitted_vals <- NULL
      resids <- NULL
      resid_type <- "Response"
      
      tryCatch({
        # Manual calculation for censReg models
        fitted_vals <- as.numeric(x$model$fitted.values)
        if (is.null(fitted_vals) || all(is.na(fitted_vals))) {
          # Calculate manually using model matrix and coefficients
          # Handle cases where model.matrix might fail
          mm <- tryCatch({
            model.matrix(x$model)
          }, error = function(e) {
            # Fallback: create simple intercept-only matrix
            matrix(1, nrow = nrow(plot_data), ncol = 1)
          })
          
          # Get coefficients safely
          model_coefs <- tryCatch({
            coef(x$model)
          }, error = function(e) {
            # Fallback: use estimates directly
            if (!is.null(x$model$estimate)) {
              x$model$estimate[1:ncol(mm)]
            } else {
              rep(mean(plot_data$value, na.rm = TRUE), ncol(mm))
            }
          })
          
          fitted_vals <- as.numeric(mm %*% model_coefs[1:ncol(mm)])
        }
        # Use residuals if available, otherwise compute manually
        resids <- as.numeric(plot_data$value - fitted_vals)
        resid_type <- "Response"
      }, error = function(e) {
        # Fall back to simple response residuals
        fitted_vals <<- rep(mean(plot_data$value, na.rm = TRUE), nrow(plot_data))
        resids <<- as.numeric(plot_data$value - fitted_vals)
        resid_type <<- "Response (fallback)"
      })
      
      # Final check that variables are properly set
      if (is.null(fitted_vals) || is.null(resids)) {
        fitted_vals <- rep(mean(plot_data$value, na.rm = TRUE), nrow(plot_data))
        resids <- as.numeric(plot_data$value - fitted_vals)
        resid_type <- "Response (fallback)"
      }
    } else {
      message("Cannot plot: censReg package required for tobit_censreg plotting")
      return(NULL)
    }
    
  } else {
    message("Cannot plot: unsupported model type")
    return(NULL)
  }
  
  list(residuals = resids, fitted = fitted_vals, type = resid_type)
}