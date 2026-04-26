# S3 plotting methods for immuno_model and immuno_model_set objects
# Provides model comparison visualizations and diagnostic plotting for fitted models

#' Plot an immuno_model (placeholder)
#'
#' @description
#' This is a placeholder for future model comparison functionality.
#' For individual model plots, use objects from fit_one() which have class 'immuno_fit'.
#'
#' @param x   An `immuno_model` object.
#' @param ... Ignored.
#'
#' @return Nothing (placeholder).
#'
#' @examples
#' \dontrun{
#' # Use this instead:
#' fit <- fit_one(data, family = "gamma")
#' plot(fit)  # This works!
#' }
#' @export
plot.immuno_model <- function(x, ...) {
  # If this object also inherits immuno_fit, delegate to that method
  if (inherits(x, "immuno_fit")) {
    return(NextMethod("plot"))
  }
  message("immuno_model plotting is not yet implemented.")
  message("For model plots, use fit_one() which returns 'immuno_fit' objects:")
  message("   fit <- fit_one(data, family = 'gamma')")
  message("   plot(fit)  # This creates residual plots!")
  invisible(NULL)
}

#' Plot LOD comparison results
#'
#' @description
#' Creates comparison plots for different LOD handling approaches.
#'
#' @param x An `immuno_lod_comparison` object
#' @param plot_best Logical. Whether to create diagnostic plots for best approach.
#' @param ... Additional arguments passed to plot methods
#'
#' @return Invisibly returns the comparison plot
#' @export
#' @importFrom ggplot2 ggplot aes geom_col geom_text facet_wrap labs theme_minimal
#' @importFrom stats reorder
plot.immuno_lod_comparison <- function(x, plot_best = TRUE, ...) {
  
  # Prepare data for plotting
  comp_data <- x$comparison
  comp_data$approach_type <- ifelse(grepl("gamma_", comp_data$model), 
                                   "Preprocessing", "Censoring-Aware")
  
  # Create comparison plot
  p_comp <- ggplot2::ggplot(comp_data, ggplot2::aes(x = reorder(model, -aic), y = aic)) +
    ggplot2::geom_col(ggplot2::aes(fill = approach_type), alpha = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = round(aic, 0)), 
                       vjust = -0.3, size = getOption("immunoplex.text_size", 3)) +
    ggplot2::facet_wrap(~ approach_type, scales = "free_x") +
    ggplot2::labs(
      title = "LOD Handling Approach Comparison",
      subtitle = paste("Best approach:", x$best_name),
      x = "Approach",
      y = "AIC",
      fill = "Method Type"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  print(p_comp)
  
  # Plot diagnostics for best approach if requested
  if (plot_best && !is.null(x$best_model)) {
    cat("\nDiagnostic plots for best approach (", x$best_name, "):\n", sep = "")
    plot(x$best_model, plot_type = "censor_aware", save_pdf = FALSE)
  }
  
  invisible(p_comp)
}
