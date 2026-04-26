# Plotting helper functions
# Core plotting utilities for creating diagnostic plots with basic and censoring-aware options

#' Create basic diagnostic plots
#' 
#' @param plot_data Data frame with plotting data
#' @param plot_context List with family, estimand, resid_type info
#' @param point_size Size of points in plots
#' @param alpha Transparency of points
#' @return List with rvf and qq plots
#' @keywords internal
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_qq geom_qq_line
#' @importFrom ggplot2 labs theme_minimal
.create_basic_plots <- function(plot_data, plot_context, point_size, alpha) {
  
  # Extract context info
  family <- plot_context$family
  estimand <- plot_context$estimand
  resid_type <- plot_context$resid_type
  
  # Build title with estimand info if available
  main_title <- paste(toupper(family), "Model -", resid_type, "Residuals")
  subtitle <- if (!is.null(estimand)) paste("Estimand:", estimand) else NULL
  
  # Residuals vs Fitted plot
  p1 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = fitted, y = residuals)) +
    ggplot2::geom_point(color = "steelblue", alpha = alpha, size = point_size) +
    ggplot2::geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    ggplot2::labs(
      title = main_title,
      subtitle = subtitle,
      x = "Fitted Values",
      y = paste(resid_type, "Residuals")
    ) +
    ggplot2::theme_minimal()
  
  # Q-Q plot
  p2 <- ggplot2::ggplot(plot_data, ggplot2::aes(sample = residuals)) +
    ggplot2::geom_qq(color = "steelblue", alpha = alpha, size = point_size) +
    ggplot2::geom_qq_line(color = "red") +
    ggplot2::labs(
      title = "Q-Q Plot",
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    ggplot2::theme_minimal()
  
  list(rvf = p1, qq = p2)
}

#' Create censoring-aware diagnostic plots
#' 
#' @param plot_data Data frame with plotting data
#' @param plot_context List with family, estimand, resid_type info
#' @param point_size Size of points in plots
#' @param alpha_cens Transparency of censored points
#' @return List with rvf and qq plots
#' @keywords internal
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_qq geom_qq_line
#' @importFrom ggplot2 labs theme_minimal scale_color_manual geom_smooth
.create_censor_aware_plots <- function(plot_data, plot_context, 
                                       point_size, alpha_cens) {
  
  # Extract context info
  family <- plot_context$family
  estimand <- plot_context$estimand
  resid_type <- plot_context$resid_type
  
  # Add censoring info if available
  if ("cens_lod" %in% names(plot_data)) {
    plot_data$censor_status <- ifelse(plot_data$cens_lod, "Censored", "Observed")
  } else {
    plot_data$censor_status <- "Observed"
  }
  
  # Try to assign back to caller's environment - simplified approach
  tryCatch({
    caller_env <- parent.frame()
    assign("plot_data", plot_data, envir = caller_env)
  }, error = function(e) {
    # Silent fail - assignment not critical for function operation
  })
  
  # Colors for censoring status
  censor_colors <- c("Observed" = "steelblue", "Censored" = "red")
  
  # Build title with estimand info if available
  main_title <- paste(toupper(family), "Model -", resid_type, "Residuals")
  subtitle_parts <- c("Points colored by censoring status")
  if (!is.null(estimand)) {
    subtitle_parts <- c(paste("Estimand:", estimand), subtitle_parts)
  }
  subtitle <- paste(subtitle_parts, collapse = " | ")
  
  # Residuals vs Fitted plot with censoring awareness
  p1 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = fitted, y = residuals, 
                                                 color = censor_status)) +
    ggplot2::geom_point(alpha = alpha_cens, size = point_size) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    ggplot2::geom_smooth(method = "loess", se = FALSE, 
                         color = getOption("immunoplex.loess_color", "gray50"), 
                         linetype = "dotted", alpha = 0.7) +
    ggplot2::scale_color_manual(values = censor_colors, name = "Status") +
    ggplot2::labs(
      title = main_title,
      subtitle = subtitle,
      x = "Fitted Values",
      y = paste(resid_type, "Residuals")
    ) +
    ggplot2::theme_minimal()
  
  # Q-Q plot with censoring markers
  p2 <- ggplot2::ggplot(plot_data, ggplot2::aes(sample = residuals)) +
    ggplot2::geom_qq(ggplot2::aes(color = censor_status), 
                     alpha = alpha_cens, size = point_size) +
    ggplot2::geom_qq_line(color = "black") +
    ggplot2::scale_color_manual(values = censor_colors, name = "Status") +
    ggplot2::labs(
      title = "Q-Q Plot",
      subtitle = "Censored observations highlighted",
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    ggplot2::theme_minimal()
  
  list(rvf = p1, qq = p2)
}