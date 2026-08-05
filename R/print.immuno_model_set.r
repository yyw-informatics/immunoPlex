# S3 methods for immuno_model_set and immuno_lod_comparison classes
# Provides comprehensive printing and summary functionality for model comparison objects

#' @title Print method for immuno_model_set
#' @description Shows model comparison table and highlights best model
#' @param x An `immuno_model_set` object
#' @param ... Unused
#' @return `x` invisibly
#' @export
print.immuno_model_set <- function(x, ...) {
  cat("immunoPlex Model Set (immuno_model_set)\n")
  cat("Models fitted:", x$n_models, "\n")
  
  # Show best model with estimand if available
  best_estimand <- if (!is.null(x$best_model$estimand)) 
    paste0(", ", x$best_model$estimand) else ""
  cat("Best model:", x$best_family, "(AIC =", round(x$best_model$aic, 2), best_estimand, ")\n\n")
  
  cat("Model Comparison:\n")
  print(x$comparison, row.names = FALSE)
  
  invisible(x)
}

#' @title Summary method for immuno_model_set
#' @description Detailed summary of all fitted models
#' @param object An `immuno_model_set` object
#' @param ... Unused
#' @return Summary object
#' @export
summary.immuno_model_set <- function(object, ...) {
  cat("<immuno_model_set Summary>\n")
  cat("Number of models:", object$n_models, "\n")
  cat("Best model:", object$best_family, "\n\n")
  
  cat("=== Model Comparison Table ===\n")
  print(object$comparison, row.names = FALSE)
  
  cat("\n=== Best Model Details ===\n")
  print(object$best_model)
  
  invisible(object)
}

#' @title Print method for immuno_lod_comparison
#' @description Shows LOD comparison results
#' @param x An `immuno_lod_comparison` object
#' @param ... Unused
#' @return `x` invisibly
#' @export
print.immuno_lod_comparison <- function(x, ...) {
  cat("<immuno_lod_comparison>\n")
  cat("LOD methods tested:", paste(x$lod_methods_tested, collapse = ", "), "\n")
  cat("Model families tested:", paste(x$families_tested, collapse = ", "), "\n")
  cat("Log transform included:", x$include_log_transform, "\n")
  cat("Best approach:", x$best_name, "(AIC =", round(x$best_model$aic, 2), ")\n\n")
  
  cat("Comparison Results:\n")
  print(x$comparison, row.names = FALSE)
  
  invisible(x)
}

#' @title Summary method for immuno_lod_comparison
#' @description Detailed LOD comparison summary
#' @param object An `immuno_lod_comparison` object
#' @param ... Unused
#' @return Summary object
#' @export
summary.immuno_lod_comparison <- function(object, ...) {
  cat("<immuno_lod_comparison Summary>\n")
  cat("Total models fitted:", length(object$models), "\n")
  cat("Best approach:", object$best_name, "\n\n")
  
  # Separate preprocessing vs censoring approaches
  preprocessing <- object$comparison[grepl("gamma_", object$comparison$model), ]
  censoring <- object$comparison[!grepl("gamma_", object$comparison$model), ]
  
  if (nrow(preprocessing) > 0) {
    cat("=== Preprocessing Approaches (Gamma Models) ===\n")
    print(preprocessing, row.names = FALSE)
    cat("\n")
  }
  
  if (nrow(censoring) > 0) {
    cat("=== Censoring-Aware Models ===\n")
    print(censoring, row.names = FALSE)
    cat("\n")
  }
  
  cat("=== Best Model Details ===\n")
  print(object$best_model)
  
  invisible(object)
}