# S3 plot methods for immuno_ancova and immuno_ancova_set objects

#' Plot diagnostics for a single rank-based ANCOVA
#'
#' @description
#' Two-panel diagnostic plot: residuals vs fitted and Q-Q plot for the
#' rank-based ANCOVA model.
#'
#' @param x An \code{immuno_ancova} object from \code{ancova_one()}.
#' @param save_pdf Logical. Save plots as PDF. Default \code{FALSE}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns a list with ggplot objects.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_smooth
#' @importFrom ggplot2 labs theme_minimal stat_qq stat_qq_line
#' @import patchwork
#' @export
plot.immuno_ancova <- function(x, save_pdf = FALSE, ...) {
  if (is.null(x$model)) {
    message("Cannot plot: model did not converge")
    return(invisible(NULL))
  }

  label <- if (!is.null(x$analyte)) x$analyte else "ANCOVA"

  plot_data <- data.frame(
    fitted = stats::fitted(x$model),
    residuals = stats::residuals(x$model)
  )

  rvf <- ggplot2::ggplot(plot_data, ggplot2::aes(x = fitted, y = residuals)) +
    ggplot2::geom_point(size = 1.2, alpha = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_smooth(method = "loess", se = FALSE, color = "steelblue",
                          formula = y ~ x, linewidth = 0.8) +
    ggplot2::labs(title = paste(label, "- Residuals vs Fitted"),
                  x = "Fitted values", y = "Residuals") +
    ggplot2::theme_minimal()

  qq <- ggplot2::ggplot(plot_data, ggplot2::aes(sample = residuals)) +
    ggplot2::stat_qq(size = 1.2, alpha = 0.6) +
    ggplot2::stat_qq_line(color = "steelblue", linewidth = 0.8) +
    ggplot2::labs(title = paste(label, "- Q-Q Plot"),
                  x = "Theoretical quantiles", y = "Sample quantiles") +
    ggplot2::theme_minimal()

  combined <- rvf | qq

  if (save_pdf) {
    out_file <- paste0(gsub("[^a-zA-Z0-9_-]", "_", label), "_ancova_diagnostics.pdf")
    ggplot2::ggsave(out_file, combined, width = 12, height = 6)
    message("Saved: ", out_file)
  }

  invisible(list(combined_plot = combined, rvf_plot = rvf, qq_plot = qq))
}


#' Plot effect sizes for a batch ANCOVA result set
#'
#' @description
#' Horizontal bar plot of partial omega-squared by analyte, colored by
#' effect magnitude category. Cohen threshold lines are overlaid.
#'
#' @param x An \code{immuno_ancova_set} object from \code{ancova_fit()}.
#' @param significance_col Character. Column to use for significance shading.
#'   Default \code{NULL} (auto-detected: \code{q_value} for binary,
#'   \code{q_omnibus} for k-level).
#' @param save_pdf Logical. Save plot as PDF. Default \code{FALSE}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_hline annotate coord_flip
#' @importFrom ggplot2 scale_y_continuous scale_fill_manual labs theme_bw theme
#' @importFrom ggplot2 element_text element_blank expansion unit guides guide_legend
#' @export
plot.immuno_ancova_set <- function(x, significance_col = NULL, save_pdf = FALSE, ...) {
  res <- x$results

  # Auto-detect significance column
  if (is.null(significance_col)) {
    significance_col <- if ("q_value" %in% names(res)) "q_value" else "q_omnibus"
  }
  if (!significance_col %in% names(res)) {
    message("Significance column '", significance_col, "' not found")
    return(invisible(NULL))
  }

  plot_data <- data.frame(
    analyte = res$analyte,
    omega_sq_partial = pmax(res$omega_sq_partial, 0, na.rm = TRUE),
    is_significant = res[[significance_col]] < x$fdr_threshold,
    stringsAsFactors = FALSE
  )
  plot_data$is_significant[is.na(plot_data$is_significant)] <- FALSE
  plot_data$fill_color <- ifelse(plot_data$is_significant,
                                  sprintf("Significant (FDR < %.2f)", x$fdr_threshold),
                                  "Not significant")

  # Order by effect size
  plot_data <- plot_data[order(plot_data$omega_sq_partial), ]
  plot_data$analyte <- factor(plot_data$analyte, levels = plot_data$analyte)

  sig_colors <- c(
    setNames("#2d2d2d", sprintf("Significant (FDR < %.2f)", x$fdr_threshold)),
    "Not significant" = "#b0b0b0"
  )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = analyte, y = omega_sq_partial,
                                                fill = fill_color)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_hline(yintercept = 0.01, linetype = "dashed", color = "gray60",
                         linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = 0.06, linetype = "dashed", color = "gray60",
                         linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = 0.14, linetype = "dashed", color = "gray60",
                         linewidth = 0.4) +
    ggplot2::annotate("text", x = 0.5, y = 0.01, label = "Small",
                       hjust = -0.1, vjust = -0.5, size = 3.5, color = "gray50") +
    ggplot2::annotate("text", x = 0.5, y = 0.06, label = "Medium",
                       hjust = -0.1, vjust = -0.5, size = 3.5, color = "gray50") +
    ggplot2::annotate("text", x = 0.5, y = 0.14, label = "Large",
                       hjust = -0.1, vjust = -0.5, size = 3.5, color = "gray50") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::scale_fill_manual(values = sig_colors, name = NULL) +
    ggplot2::labs(
      x = NULL,
      y = expression(paste("Partial ", omega[p]^2))
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1))

  if (save_pdf) {
    n <- nrow(plot_data)
    fig_h <- max(6, n * 0.38 + 2)
    out_file <- "ancova_effect_sizes.pdf"
    ggplot2::ggsave(out_file, p, width = 5.5, height = fig_h)
    message("Saved: ", out_file)
  }

  invisible(p)
}
