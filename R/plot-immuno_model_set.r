# Plotting method for immuno_model_set objects
# Comprehensive model comparison visualization with AIC analysis and residual diagnostics

#' Plot model comparison for immuno_model_set objects
#'
#' @description
#' Creates comparison plots for multiple fitted models, showing AIC comparisons
#' and optionally residual diagnostics for the best model.
#'
#' @param x An `immuno_model_set` object from `fit_models()`.
#' @param plot_type Character. Type of plot: "comparison" (default) for AIC comparison,
#'   "best_residuals" for residuals of best model, or "all_residuals" for all models.
#' @param save_pdf Logical. Whether to save plots as PDF (default TRUE).
#' @param ... Additional arguments passed to plot.immuno_fit for residual plots.
#'
#' @return Invisibly returns the plot object(s).
#'
#' @examples
#' \dontrun{
#' models <- fit_models(data, families = c("gamma", "tobit", "aft"))
#' plot(models)  # AIC comparison
#' plot(models, plot_type = "best_residuals")  # Best model residuals
#' }
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_col geom_text labs theme_minimal
#' @importFrom ggplot2 coord_flip scale_fill_viridis_d
plot.immuno_model_set <- function(x, 
                                  plot_type = c("comparison", "best_residuals", "all_residuals"),
                                  save_pdf = TRUE,
                                  ...) {
  
  if (!inherits(x, "immuno_model_set")) {
    stop("Object must be of class 'immuno_model_set'")
  }
  
  plot_type <- match.arg(plot_type)
  
  if (plot_type == "comparison") {
    # Create AIC comparison plot
    comp_data <- x$comparison
    comp_data$family <- factor(comp_data$family, levels = comp_data$family[order(comp_data$aic)])
    
    p <- ggplot2::ggplot(comp_data, ggplot2::aes(x = family, y = delta_aic, fill = family)) +
      ggplot2::geom_col(alpha = 0.8) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("AIC: %.1f", aic)), 
                        hjust = -0.1, size = 3) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = "Model Comparison by AIC",
        subtitle = paste("Optimal model:", x$best_family, "(minimum AIC)"),
        x = "Model Family",
        y = "Delta AIC (difference from best)"
      ) +
      ggplot2::scale_fill_viridis_d(guide = "none") +
      ggplot2::theme_minimal()
    
    if (save_pdf) {
      out_file <- "model_comparison.pdf"
      ggplot2::ggsave(
        filename = out_file,
        plot = p,
        width = 10,
        height = 6,
        device = "pdf"
      )
      message("Model comparison plot saved as: ", out_file)
    }
    
    print(p)
    return(invisible(p))

  } else if (plot_type == "best_residuals") {
    # Plot residuals for best model
    message("Plotting residuals for best model: ", x$best_family)
    # Call plot method directly on the fit object to avoid the immuno_model method
    if (inherits(x$best_model, "immuno_fit")) {
      plot.immuno_fit(x$best_model, save_pdf = save_pdf, ...)
    } else {
      plot(x$best_model, save_pdf = save_pdf, ...)
    }
    
  } else if (plot_type == "all_residuals") {
    # Plot residuals for all converged models
    plots <- list()
    for (family_name in names(x$models)) {
      message("Plotting residuals for ", family_name, " model")
      if (inherits(x$models[[family_name]], "immuno_fit")) {
        plots[[family_name]] <- plot.immuno_fit(x$models[[family_name]], save_pdf = save_pdf, ...)
      } else {
        plots[[family_name]] <- plot(x$models[[family_name]], save_pdf = save_pdf, ...)
      }
    }
    return(invisible(plots))
  }
}