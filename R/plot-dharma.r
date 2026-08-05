# DHARMa computation functions for plotting
# Provides simulated residual diagnostics for various model types including censored models

#' Compute DHARMa residuals for model diagnostics
#' 
#' @importFrom stats pnorm predict model.matrix rnorm coef residuals fitted
#' @param x An immuno_fit object
#' @param nsim Number of simulations for DHARMa
#' @return DHARMa object or NULL if computation fails
#' @keywords internal
.compute_dharma <- function(x, nsim) {
  if (!requireNamespace("DHARMa", quietly = TRUE)) {
    message("DHARMa package required for simulated residuals")
    return(NULL)
  }
  
  # Try DHARMa for natively supported model classes
  supported_classes <- c("glmmTMB", "glm", "lmerMod", "glmerMod")
  
  if (!inherits(x$model, supported_classes)) {
    # For unsupported classes, try external simulation approach
    if (x$family == "tobit_censreg" && inherits(x$model, "censReg")) {
      message("Using external simulation approach for censReg model")
      return(.compute_dharma_external_censreg(x, nsim))
    } else {
      message("DHARMa diagnostics not available for ", class(x$model)[1], " models")
      message("Natively supported classes: ", paste(supported_classes, collapse = ", "))
      message("Use external simulation approach for unsupported models")
      return(NULL)
    }
  }
  
  # For natively supported models, use DHARMa directly
  if (inherits(x$model, "glmmTMB")) {
    tryCatch({
      DHARMa::simulateResiduals(x$model, n = nsim, plot = FALSE)
    }, error = function(e) {
      message("DHARMa simulation failed: ", e$message)
      NULL
    })
  } else if (inherits(x$model, "survreg")) {
    # For survival models, use randomized quantile residuals instead
    if (requireNamespace("statmod", quietly = TRUE)) {
      tryCatch({
        resids <- statmod::qresiduals(x$model)
        # Create a basic DHARMa-like object for compatibility
        list(
          scaledResiduals = pnorm(resids),
          fittedPredictedResponse = predict(x$model),
          observedResponse = x$data_used$value,
          DHARMa = TRUE
        )
      }, error = function(e) {
        message("Randomized quantile residuals failed: ", e$message)
        NULL
      })
    } else {
      message("statmod package required for randomized quantile residuals")
      NULL
    }
  } else {
    # For other model types, try basic simulation if possible
    tryCatch({
      DHARMa::simulateResiduals(x$model, n = nsim, plot = FALSE)
    }, error = function(e) {
      message("DHARMa simulation failed: ", e$message)
      NULL
    })
  }
}

#' External DHARMa simulation for censReg models
#' 
#' @param x An immuno_fit object with censReg model
#' @param nsim Number of simulations
#' @return DHARMa object or NULL if computation fails
#' @keywords internal
.compute_dharma_external_censreg <- function(x, nsim) {
  if (!inherits(x$model, "censReg")) {
    message("External simulation only available for censReg models")
    return(NULL)
  }
  
  tryCatch({
    message("Using external simulation approach for censReg model")
    
    # Get model matrix and coefficients (censReg doesn't have predict method)
    # More robust formula extraction
    original_formula <- x$model$call$formula
    if (is.null(original_formula)) {
      # Fallback: use simple intercept-only model
      data_formula <- stats::as.formula("value ~ 1")
    } else {
      # Extract variable names more safely
      formula_vars <- tryCatch({
        all.vars(original_formula)
      }, error = function(e) {
        c("1")  # Fallback to intercept-only
      })
      
      # Remove response variable (first one) and construct formula
      predictor_vars <- formula_vars[-1]
      if (length(predictor_vars) == 0) {
        predictor_vars <- "1"  # Intercept-only model
      }
      
      data_formula <- stats::as.formula(paste("value ~", 
        paste(predictor_vars, collapse = " + ")))
    }
    
    X <- model.matrix(data_formula, data = x$data_used)
    
    # Extract coefficients more robustly
    estimates <- x$model$estimate
    if (is.null(estimates)) {
      stop("No coefficient estimates available in censReg model")
    }
    
    # Get beta coefficients (first ncol(X) estimates)
    n_coef <- ncol(X)
    if (length(estimates) < n_coef) {
      stop("Insufficient coefficient estimates: need ", n_coef, " but got ", length(estimates))
    }
    
    beta <- estimates[1:n_coef]
    
    # Extract sigma parameter
    if ("logSigma" %in% names(estimates)) {
      sigma <- exp(estimates["logSigma"])
    } else if ("sigma" %in% names(estimates)) {
      sigma <- estimates["sigma"]
    } else {
      # Fallback: use last estimate as sigma
      sigma <- exp(estimates[length(estimates)])
    }
    
    # Calculate predicted values manually
    pred_values <- as.vector(X %*% beta)
    
    # Get LOD threshold from data (use original LOD for proper censoring)
    lod_threshold <- unique(x$data_used$lod)[1]
    
    # Generate simulations with proper censoring bounds
    n_obs <- length(pred_values)
    sim_matrix <- matrix(nrow = n_obs, ncol = nsim)
    
    for (i in 1:nsim) {
      # Simulate from normal distribution
      sim_raw <- rnorm(n_obs, mean = pred_values, sd = sigma)
      
      # Apply censoring based on data flags (more robust than simple threshold)
      sim_censored <- sim_raw
      if ("cens_lod" %in% names(x$data_used)) {
        # Left-censored observations get LOD value
        left_cens <- x$data_used$cens_lod & !is.na(x$data_used$cens_lod)
        sim_censored[left_cens] <- pmax(sim_raw[left_cens], x$data_used$lod[left_cens])
      }
      if ("cens_ulod" %in% names(x$data_used)) {
        # Right-censored observations get ULOD value  
        right_cens <- x$data_used$cens_ulod & !is.na(x$data_used$cens_ulod)
        sim_censored[right_cens] <- pmin(sim_raw[right_cens], x$data_used$ulod[right_cens])
      }
      
      sim_matrix[, i] <- sim_censored
    }
    
    # Create DHARMa object
    DHARMa::createDHARMa(
      simulatedResponse = sim_matrix,
      observedResponse = x$data_used$value,
      fittedPredictedResponse = pred_values,
      integerResponse = FALSE
    )
  }, error = function(e) {
    message("External DHARMa simulation failed: ", e$message)
    NULL
  })
}