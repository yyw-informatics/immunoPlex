# Statistical Model Diagnostic Framework
#
# Implements systematic model evaluation protocols:
# - Distribution family validation
# - Parameter estimation assessment
# - Residual analysis and diagnostics
# - Information criterion comparison
# - Standardized data processing
#
# Core components:
# - setup_test_data(): Standardized dataset initialization
# - fit_with_diagnostics(): Comprehensive model assessment
# - run_diagnostics(): Systematic evaluation protocol
# - inspect_results(): Structured diagnostic reporting

# Load package functions using devtools
if (!exists("fit_one")) {
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(quiet = TRUE)
  } else {
    # Fallback: source functions directly
    source(file.path("R", "fit_one.r"))
    source(file.path("R", "fit_models.r"))
  }
}

# Structured diagnostic report generation
inspect_results <- function(results) {
  cat("🔍 RESULTS SUMMARY\n", rep("=", 20), "\n")
  
  for (family_name in names(results)) {
    fit <- results[[family_name]]
    cat("\n📊", toupper(family_name), ":\n")
    cat("   Converged:", fit$converged, "\n")
    if (fit$converged) {
      cat("   AIC:", round(fit$aic, 1), "\n")
      cat("   BIC:", round(fit$bic, 1), "\n")
      cat("   LogLik:", round(fit$logLik, 1), "\n")
      cat("   Family:", fit$family, "\n")
      cat("   Censored obs:", fit$n_cens_lod, "/", fit$n_cens_lod + fit$n_cens_ulod, "\n")
    } else {
      cat("   ❌ Model failed to converge\n")
    }
  }
  cat("\n💡 Access models: results$tobit$model\n")
  cat("📋 Model summary: summary(results$tobit$model)\n")
  cat("📊 Plot residuals: plot(results$tobit)\n")
}

# Initialize standardized evaluation dataset
setup_test_data <- function() {
  data("immunoplex_example", package = "immunoPlex")
  
  raw_expr <- immunoplex_example$expression
  meta <- immunoplex_example$metadata
  lod_lookup <- immunoplex_example$lod_lookup
  
  # Select cytokine with complete LOD information
  target_cyto <- lod_lookup$cytokine[which(!is.na(lod_lookup$lod))[1]]
  target_lod <- lod_lookup$lod[lod_lookup$cytokine == target_cyto]
  
  # Initialize data frame with aligned expression and metadata
  # Maintain index alignment between expression matrix and metadata
  dat <- data.frame(
    sample_id = meta$sample_id,
    subject_id = meta$subject_id,
    value = raw_expr[, target_cyto],
    timepoint = meta$timepoint,
    disease = meta$disease,
    age = meta$age,
    lod = target_lod,
    ulod = lod_lookup$ulod[lod_lookup$cytokine == target_cyto],
    stringsAsFactors = FALSE
  )
  
  dat$cens_lod <- dat$value < dat$lod
  dat$cens_ulod <- !is.na(dat$ulod) & dat$value > dat$ulod
  
  # Report dataset characteristics
  cat("📊 DATASET:", target_cyto, "| Observations:", nrow(dat), 
      "| Censoring rate:", round(100*mean(dat$cens_lod), 1), "% (n =", sum(dat$cens_lod), ")\n")
  
  return(dat)
}

# Execute comprehensive model assessment protocol
fit_with_diagnostics <- function(dat, family_name) {
  cat("\n🔧", toupper(family_name), "MODEL:\n")
  
  fit <- fit_one(dat, family = family_name, random = "")
  
  if (!fit$converged) {
    cat("❌ Failed to converge\n")
    return(fit)
  }
  
  cat("✅ Converged | AIC:", round(fit$aic, 1), "| LogLik:", round(fit$logLik, 1), "\n")
  
  # Distribution-specific diagnostic evaluation
  if (family_name == "gamma" && inherits(fit$model, "glmmTMB")) {
    m <- fit$model
    cat("   Convergence code:", m$fit$convergence, "| Hessian OK:", m$sdr$pdHess, "\n")
    
    # Evaluate random effects parameter estimation
    theta <- m$fit$par[names(m$fit$par) == "theta"]
    if (length(theta) > 0) {
      cat("   RE variance:", round(theta^2, 6), "| Singular:", any(abs(theta) < 1e-4), "\n")
    }
    
    # Residual check
    resids <- residuals(m, type = "pearson")
    cat("   Residuals - Mean:", round(mean(resids), 3), 
        "| SD:", round(sd(resids), 3),
        "| Outliers:", sum(abs(resids) > 3), "\n")
    
  } else if (family_name %in% c("tobit", "aft") && inherits(fit$model, "survreg")) {
    m <- fit$model
    cat("   Iterations:", m$iter, "| vcov OK:", !any(is.na(vcov(m))), "\n")
    
    resids <- residuals(m, type = "deviance")
    cat("   Residuals - Mean:", round(mean(resids), 3),
        "| SD:", round(sd(resids), 3), 
        "| Outliers:", sum(abs(resids) > 3), "\n")
  }
  
  return(fit)
}

# Execute systematic model evaluation protocol
run_diagnostics <- function(plot_results = TRUE) {
  cat("🔍 immunoPlex Model Diagnostics\n", rep("=", 35), "\n")
  
  dat <- setup_test_data()
  families <- c("gamma", "tobit", "aft", "auto")
  results <- list()
  
  # Fit all models
  for (fam in families) {
    results[[fam]] <- fit_with_diagnostics(dat, fam)
  }
  
  # Compare converged models
  converged <- results[sapply(results, function(x) x$converged)]
  
  if (length(converged) > 1) {
    cat("\n📈 MODEL COMPARISON:\n")
    comparison <- data.frame(
      Family = names(converged),
      AIC = sapply(converged, function(x) round(x$aic, 1)),
      BIC = sapply(converged, function(x) round(x$bic, 1)),
      LogLik = sapply(converged, function(x) round(x$logLik, 1))
    )
    print(comparison)
    
    best <- names(converged)[which.min(sapply(converged, function(x) x$aic))]
    cat("🏆 Best model:", best, "\n")
    
    # Plot best model residuals using the new plot method
    if (plot_results && interactive()) {
      cat("\n📊 Creating residual plots for best model...\n")
      plot(converged[[best]])
    }
  }
  
  cat("\n✅ Diagnostics complete!\n")
  return(results)
}

# Interactive session initialization protocol
if (interactive()) {
  cat("💡 Protocol: results <- run_diagnostics()\n")
  cat("📊 Visualization: plot(results$gamma)  # Generates diagnostic plots\n")
  cat("🔍 Assessment: summary(results$gamma$model)\n")
  cat("👀 Evaluation: inspect_results(results)\n\n")
}
