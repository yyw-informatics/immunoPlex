# PLS-DA Plotting Module for immunoPlex Package
# Visualization functions for PLS-DA model results

#' Plot PLS-DA model results
#'
#' @description
#' Generate publication-quality visualizations for PLS-DA model results.
#' Supports multiple plot types including score plots, loading plots,
#' VIP plots, and biplots.
#'
#' @param x A plsda_model object from plsda_fit()
#' @param plot_type character string specifying plot type. Options:
#'   "scores" (sample score plot), "loadings" (cytokine loadings),
#'   "vip" (Variable Importance in Projection), "biplot" (scores + loadings overlay)
#' @param components numeric vector of length 2 specifying which components to plot.
#'   Default: c(1, 2)
#' @param color_by character string naming the metadata column for coloring points.
#'   If NULL, uses the response variable. Default: NULL
#' @param shape_by character string naming the metadata column for point shapes.
#'   Default: NULL
#' @param add_ellipse logical indicating whether to add 95% confidence ellipses
#'   for score plots. Default: TRUE
#' @param vip_threshold numeric threshold for highlighting important variables
#'   in VIP plots. Default: 1.0
#' @param top_n integer specifying number of top cytokines to show in loading
#'   and VIP plots. If NULL, shows all. Default: NULL
#' @param arrow_scale numeric scaling factor for loading arrows in biplot.
#'   Default: 0.7
#' @param label_repel logical indicating whether to use ggrepel for labels
#'   in biplot. Default: TRUE
#' @param ... Additional arguments passed to specific plot functions
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' # Fit model
#' model <- plsda_fit(preprocessed, response_var = "group")
#' 
#' # Score plot
#' plot(model, plot_type = "scores", color_by = "group")
#' 
#' # Loading plot
#' plot(model, plot_type = "loadings", top_n = 20)
#' 
#' # VIP plot
#' plot(model, plot_type = "vip", vip_threshold = 1)
#' 
#' # Biplot
#' plot(model, plot_type = "biplot", color_by = "group", shape_by = "site")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline stat_ellipse
#'   theme_classic theme_bw labs scale_color_manual scale_shape_manual
#'   element_text element_blank element_rect geom_bar coord_flip geom_segment 
#'   arrow unit geom_text facet_wrap scale_fill_manual guides guide_legend 
#'   scale_x_continuous scale_y_continuous geom_jitter .data
#' @importFrom ggrepel geom_label_repel geom_text_repel
#' @importFrom dplyr group_by summarise filter arrange slice pull left_join
#' @importFrom tidyr pivot_longer
#' @export
plot.plsda_model <- function(x,
                              plot_type = c("scores", "loadings", "vip", "biplot"),
                              components = c(1, 2),
                              color_by = NULL,
                              shape_by = NULL,
                              add_ellipse = TRUE,
                              vip_threshold = 1.0,
                              top_n = NULL,
                              arrow_scale = 0.7,
                              label_repel = TRUE,
                              ...) {
  
  plot_type <- match.arg(plot_type)
  
  # Route to specific plotting function
  if (plot_type == "scores") {
    plot_pls_scores(x, components, color_by, shape_by, add_ellipse, ...)
  } else if (plot_type == "loadings") {
    plot_pls_loadings(x, components, top_n, ...)
  } else if (plot_type == "vip") {
    plot_pls_vip(x, vip_threshold, top_n, ...)
  } else if (plot_type == "biplot") {
    plot_pls_biplot(x, components, color_by, shape_by, add_ellipse,
                    arrow_scale, label_repel, ...)
  }
}


#' Create PLS-DA score plot
#'
#' @description
#' Plot sample scores on PLS-DA components with optional confidence ellipses.
#'
#' @param model A plsda_model object
#' @param components Numeric vector of length 2 for component indices
#' @param color_by Column name for coloring points
#' @param shape_by Column name for point shapes
#' @param add_ellipse Logical for adding confidence ellipses
#' @param point_size Numeric size of points. Default: 3
#' @param ... Additional arguments
#'
#' @return A ggplot object
#' @export
plot_pls_scores <- function(model,
                            components = c(1, 2),
                            color_by = NULL,
                            shape_by = NULL,
                            add_ellipse = TRUE,
                            point_size = 3,
                            ...) {
  
  # Validate input
  if (!inherits(model, "plsda_model")) {
    stop("model must be a plsda_model object")
  }
  
  # Get score data
  scores <- model$scores
  stats <- model$model_stats
  
  # Determine color variable
  if (is.null(color_by)) {
    color_by <- model$response_var
  }
  
  # Check if requested variable exists
  if (!color_by %in% names(scores)) {
    stop("color_by variable '", color_by, "' not found in model scores")
  }
  
  # Get component names
  comp_cols <- grep("^p[0-9]+$|^o[0-9]+$", names(scores), value = TRUE)
  
  if (length(components) == 1 || length(comp_cols) == 1) {
    # 1D plot with jittering
    comp_name <- comp_cols[components[1]]
    var_explained <- round(stats$R2X[components[1]] * 100, 1)
    
    p <- ggplot(scores, aes(x = .data[[comp_name]], y = 0, 
                           color = .data[[color_by]])) +
      geom_jitter(height = 0.1, width = 0, alpha = 0.7, size = point_size) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      labs(
        x = paste0("Component ", components[1], " (", var_explained, "%)"),
        y = "",
        title = paste(model$method, "Score Plot"),
        color = color_by
      ) +
      theme_classic() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
  } else {
    # 2D plot
    comp_x <- comp_cols[components[1]]
    comp_y <- comp_cols[components[2]]
    var_x <- round(stats$R2X[components[1]] * 100, 1)
    var_y <- round(stats$R2X[components[2]] * 100, 1)
    
    # Base plot
    p <- ggplot(scores, aes(x = .data[[comp_x]], y = .data[[comp_y]], 
                           color = .data[[color_by]])) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_point(size = point_size, alpha = 0.8) +
      labs(
        x = paste0("Component ", components[1], " (", var_x, "%)"),
        y = paste0("Component ", components[2], " (", var_y, "%)"),
        title = paste(model$method, "Score Plot"),
        color = color_by
      ) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    # Add ellipses if requested
    if (add_ellipse) {
      p <- p + stat_ellipse(
        aes(fill = .data[[color_by]]),
        geom = "polygon",
        alpha = 0.1,
        level = 0.95,
        show.legend = FALSE
      )
    }
    
    # Add shapes if requested
    if (!is.null(shape_by) && shape_by %in% names(scores)) {
      p <- p + aes(shape = .data[[shape_by]]) +
        labs(shape = shape_by)
    }
  }
  
  # Try to use NEJM colors if ggsci available
  if (requireNamespace("ggsci", quietly = TRUE)) {
    p <- p + ggsci::scale_color_nejm() +
      ggsci::scale_fill_nejm()
  }
  
  return(p)
}


#' Create PLS-DA loading plot
#'
#' @description
#' Bar plot showing cytokine loadings on PLS-DA components.
#'
#' @param model A plsda_model object
#' @param components Numeric vector of component indices to plot
#' @param top_n Number of top cytokines to show (by absolute loading)
#' @param ... Additional arguments
#'
#' @return A ggplot object
#' @export
plot_pls_loadings <- function(model,
                              components = c(1, 2),
                              top_n = NULL,
                              ...) {
  
  if (!inherits(model, "plsda_model")) {
    stop("model must be a plsda_model object")
  }
  
  # Get loadings
  loadings <- model$loadings
  
  # Get component columns
  comp_cols <- grep("^p[0-9]+$|^o[0-9]+$", names(loadings), value = TRUE)
  comp_names <- comp_cols[components]
  
  # Reshape to long format for plotting
  loadings_long <- tidyr::pivot_longer(
    loadings,
    cols = all_of(comp_names),
    names_to = "component",
    values_to = "loading"
  )
  
  # Filter top N if specified
  if (!is.null(top_n)) {
    # Get top cytokines by maximum absolute loading across components
    top_cyto <- loadings_long %>%
      dplyr::group_by(cytokine) %>%
      dplyr::summarise(max_abs_loading = max(abs(loading)), .groups = "drop") %>%
      dplyr::arrange(desc(max_abs_loading)) %>%
      dplyr::slice(1:min(top_n, nrow(.))) %>%
      dplyr::pull(cytokine)
    
    loadings_long <- loadings_long %>%
      dplyr::filter(cytokine %in% top_cyto)
  }
  
  # Create plot
  p <- ggplot(loadings_long, aes(x = reorder(cytokine, loading), 
                                 y = loading, 
                                 fill = loading > 0)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("TRUE" = "tomato", "FALSE" = "slateblue1"),
                      guide = "none") +
    coord_flip() +
    facet_wrap(~ component, scales = "free_x") +
    labs(
      x = "Cytokine",
      y = "Loading Value",
      title = paste(model$method, "Component Loadings")
    ) +
    theme_classic() +
    theme(
      strip.background = element_rect(fill = "lightblue"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 8)
    )
  
  return(p)
}


#' Create VIP score plot
#'
#' @description
#' Horizontal bar plot of Variable Importance in Projection scores.
#'
#' @param model A plsda_model object
#' @param vip_threshold Threshold for highlighting important variables
#' @param top_n Number of top cytokines to show
#' @param ... Additional arguments
#'
#' @return A ggplot object
#' @export
plot_pls_vip <- function(model,
                         vip_threshold = 1.0,
                         top_n = NULL,
                         ...) {
  
  if (!inherits(model, "plsda_model")) {
    stop("model must be a plsda_model object")
  }
  
  # Get VIP scores
  vip_df <- model$vip_scores
  
  # Filter top N if specified
  if (!is.null(top_n)) {
    vip_df <- vip_df %>%
      dplyr::slice(1:min(top_n, nrow(.)))
  }
  
  # Add importance category
  vip_df$importance <- ifelse(vip_df$vip_score >= vip_threshold,
                               "Important", "Less Important")
  
  # Create plot
  p <- ggplot(vip_df, aes(x = reorder(cytokine, vip_score), 
                         y = vip_score, 
                         fill = importance)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = vip_threshold, linetype = "dashed", 
               color = "red", linewidth = 0.8) +
    scale_fill_manual(
      values = c("Important" = "darkred", "Less Important" = "gray70"),
      name = paste0("VIP >= ", vip_threshold)
    ) +
    coord_flip() +
    labs(
      x = "Cytokine",
      y = "VIP Score",
      title = paste(model$method, "- Variable Importance in Projection")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 8),
      legend.position = "top"
    )
  
  return(p)
}


#' Create PLS-DA biplot
#'
#' @description
#' Biplot combining sample scores and cytokine loadings with quadrant statistics.
#'
#' @param model A plsda_model object
#' @param components Numeric vector of length 2 for component indices
#' @param color_by Column name for coloring sample points
#' @param shape_by Column name for point shapes
#' @param add_ellipse Logical for adding confidence ellipses
#' @param arrow_scale Scaling factor for loading arrows
#' @param label_repel Use ggrepel for labels
#' @param top_loadings Number of top loadings to label (NULL = all)
#' @param ... Additional arguments
#'
#' @return A ggplot object
#' @export
plot_pls_biplot <- function(model,
                            components = c(1, 2),
                            color_by = NULL,
                            shape_by = NULL,
                            add_ellipse = TRUE,
                            arrow_scale = 0.7,
                            label_repel = TRUE,
                            top_loadings = 10,
                            ...) {
  
  if (!inherits(model, "plsda_model")) {
    stop("model must be a plsda_model object")
  }
  
  # Get scores and loadings
  scores <- model$scores
  loadings <- model$loadings
  stats <- model$model_stats
  
  # Determine color variable
  if (is.null(color_by)) {
    color_by <- model$response_var
  }
  
  # Get component names
  comp_cols <- grep("^p[0-9]+$", names(scores), value = TRUE)
  comp_x <- comp_cols[components[1]]
  comp_y <- comp_cols[components[2]]
  
  # Get variance explained
  var_x <- round(stats$R2X[components[1]] * 100, 1)
  var_y <- round(stats$R2X[components[2]] * 100, 1)
  
  # Prepare loading data for arrows
  load_x <- paste0("p", components[1])
  load_y <- paste0("p", components[2])
  
  loadings_plot <- data.frame(
    cytokine = loadings$cytokine,
    Comp1 = loadings[[load_x]],
    Comp2 = loadings[[load_y]]
  )
  
  # Calculate arrow scaling
  score_range <- max(abs(c(scores[[comp_x]], scores[[comp_y]])), na.rm = TRUE)
  load_range <- max(abs(c(loadings_plot$Comp1, loadings_plot$Comp2)), na.rm = TRUE)
  scale_factor <- score_range / load_range * arrow_scale
  
  loadings_plot$Comp1_scaled <- loadings_plot$Comp1 * scale_factor
  loadings_plot$Comp2_scaled <- loadings_plot$Comp2 * scale_factor
  
  # Filter top loadings if specified
  if (!is.null(top_loadings)) {
    loadings_plot$magnitude <- sqrt(loadings_plot$Comp1^2 + loadings_plot$Comp2^2)
    loadings_plot <- loadings_plot %>%
      dplyr::arrange(desc(magnitude)) %>%
      dplyr::slice(1:min(top_loadings, nrow(.)))
  }
  
  # Calculate quadrant statistics
  quad_stats <- scores %>%
    dplyr::group_by(.data[[color_by]]) %>%
    dplyr::summarise(
      q1 = mean(.data[[comp_x]] > 0 & .data[[comp_y]] > 0) * 100,
      q2 = mean(.data[[comp_x]] < 0 & .data[[comp_y]] > 0) * 100,
      q3 = mean(.data[[comp_x]] < 0 & .data[[comp_y]] < 0) * 100,
      q4 = mean(.data[[comp_x]] > 0 & .data[[comp_y]] < 0) * 100,
      .groups = "drop"
    )
  
  # Create base plot with scores
  p <- ggplot(scores, aes(x = .data[[comp_x]], y = .data[[comp_y]], 
                         color = .data[[color_by]])) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(size = 3, alpha = 0.7, stroke = 1) +
    labs(
      x = paste0("Component ", components[1], " (", var_x, "%)"),
      y = paste0("Component ", components[2], " (", var_y, "%)"),
      title = paste(model$method, "Biplot"),
      color = color_by
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Add ellipses if requested
  if (add_ellipse) {
    p <- p + stat_ellipse(
      aes(fill = .data[[color_by]]),
      geom = "polygon",
      alpha = 0.1,
      level = 0.95
    )
  }
  
  # Add loading arrows
  p <- p + geom_segment(
    data = loadings_plot,
    aes(x = 0, y = 0, xend = Comp1_scaled, yend = Comp2_scaled),
    arrow = arrow(length = unit(0.2, "cm")),
    linewidth = 0.4,
    color = "gray30",
    inherit.aes = FALSE
  )
  
  # Add loading labels
  if (label_repel && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_label_repel(
      data = loadings_plot,
      aes(x = Comp1_scaled * 1.02, y = Comp2_scaled * 1.02, label = cytokine),
      size = 2.5,
      fontface = "bold",
      fill = "white",
      color = "black",
      box.padding = 0.2,
      point.padding = 0.1,
      segment.size = 0.2,
      segment.color = "gray70",
      inherit.aes = FALSE
    )
  } else {
    p <- p + geom_text(
      data = loadings_plot,
      aes(x = Comp1_scaled * 1.1, y = Comp2_scaled * 1.1, label = cytokine),
      size = 2.5,
      fontface = "bold",
      inherit.aes = FALSE
    )
  }
  
  # Add shapes if requested
  if (!is.null(shape_by) && shape_by %in% names(scores)) {
    p <- p + aes(shape = .data[[shape_by]]) +
      labs(shape = shape_by)
  }
  
  # Try to use NEJM colors if ggsci available
  if (requireNamespace("ggsci", quietly = TRUE)) {
    p <- p + ggsci::scale_color_nejm() +
      ggsci::scale_fill_nejm()
  }
  
  return(p)
}
