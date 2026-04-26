# Main S3 plotting method for immuno_fit objects
# Comprehensive residual diagnostic plotting with censoring awareness and DHARMa support

#' Plot residual diagnostics for immuno_fit objects
#'
#' @description
#' Creates residual diagnostic plots for fitted models with optional censoring
#' awareness. Supports randomized quantile residuals for censored models and
#' DHARMa simulated diagnostics.
#'
#' @param x An `immuno_fit` object from `fit_one()`.
#' @param plot_type Character. Either "basic" (default) for standard plots or
#'   "censor_aware" for censoring-aware visualization.
#' @param dharma Logical. If TRUE, includes DHARMa residual diagnostics.
#' @param nsim Integer. Number of simulations for DHARMa (default from config, typically 1000).
#' @param point_size Numeric. Size of points in plots (default 1.2).
#' @param alpha_cens Numeric. Transparency for censored points (default 0.6).
#' @param save_pdf Logical. Whether to save plots as PDF (default TRUE).
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns a list containing:
#'   \itemize{
#'     \item `combined_plot`: The patchwork combined plot
#'     \item `rvf_plot`: Residuals vs fitted plot
#'     \item `qq_plot`: Q-Q plot  
#'     \item `dharma_obj`: DHARMa object (if dharma=TRUE)
#'   }
#'
#' @details
#' For censored models (tobit/aft), randomized quantile residuals are used when
#' `statmod` is available, otherwise deviance residuals are used. Censoring-aware
#' plots color points by censoring status and add interpretive elements.
#'
#' @examples
#' \dontrun{
#' fit <- fit_one(data, family = "tobit")
#' plot(fit)  # Basic plots
#' plot(fit, plot_type = "censor_aware")  # Censoring-aware
#' plot(fit, dharma = TRUE)  # Include DHARMa diagnostics
#' }
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_qq geom_qq_line
#' @importFrom ggplot2 labs theme_minimal ggsave scale_color_manual geom_smooth
#' @import patchwork
plot.immuno_fit <- function(x, 
                            plot_type = c("basic", "censor_aware"),
                            dharma = FALSE,
                            nsim = getOption("immunoplex.dharma_nsim", 1000),
                            point_size = getOption("immunoplex.point_size", 1.2),
                            alpha_cens = getOption("immunoplex.alpha_cens", 0.6),
                            save_pdf = TRUE,
                            ...) {
  
  if (!x$converged) {
    message("Cannot plot: model did not converge")
    return(invisible(NULL))
  }
  
  plot_type <- match.arg(plot_type)
  
  # Extract data used in fitting
  plot_data <- x$data_used
  if (is.null(plot_data)) {
    message("Cannot plot: no data_used in model object (likely small dataset < 10 rows)")
    return(invisible(NULL))
  }
  
  # Get residuals and fitted values
  resid_info <- .get_residuals_and_fitted(x, plot_data)
  if (is.null(resid_info)) return(invisible(NULL))
  
  plot_data$residuals <- resid_info$residuals
  plot_data$fitted <- resid_info$fitted
  resid_type <- resid_info$type
  
  # DHARMa diagnostics if requested
  dharma_obj <- NULL
  if (dharma) {
    dharma_obj <- .compute_dharma(x, nsim)
  }
  
  # Create plots based on type
  # Include estimand info in plot context
  plot_context <- list(
    family = x$family,
    estimand = x$estimand,
    resid_type = resid_type
  )
  
  if (plot_type == "censor_aware") {
    plots <- .create_censor_aware_plots(plot_data, plot_context, 
                                        point_size, alpha_cens)
  } else {
    plots <- .create_basic_plots(plot_data, plot_context, 
                                 point_size, alpha_cens)
  }
  
  # Combine plots
  combined_plot <- plots$rvf | plots$qq
  
  # Save PDF if requested
  if (save_pdf) {
    out_file <- paste0(x$family, "_model_residuals.pdf")
    ggplot2::ggsave(
      filename = out_file,
      plot = combined_plot,
      width = 12,
      height = 6,
      device = "pdf"
    )
    message("Residual plots saved as: ", out_file)
  }
  
  # Return results
  result <- list(
    combined_plot = combined_plot,
    rvf_plot = plots$rvf,
    qq_plot = plots$qq,
    dharma_obj = dharma_obj
  )
  
  invisible(result)
}