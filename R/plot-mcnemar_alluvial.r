# Suppress CMD check notes for NSE variables
utils::globalVariables(c("baseline", "comparison", "transition", "stratum", 
                         "x", "y", "label", "hjust", "vjust"))

#' Create Alluvial Plot for McNemar Detection Results (fixed labels + %)
#'
#' Generates an alluvial (Sankey) diagram showing detection status transitions
#' between baseline and comparison timepoints. Flows are colored by transition
#' type: gained detection (blue), lost detection (pink), stable detected (gray),
#' and stable undetected (white).
#'
#' @param mcnemar_result A data frame or tibble from \code{mcnemar_detection()}.
#' @param cytokine Character string. Name of the cytokine to plot. If NULL (default),
#'   plots the most significant cytokine (lowest q-value).
#' @param baseline_label Character string. Label for baseline timepoint (default: "Enrollment")
#' @param comparison_label Character string. Label for comparison timepoint (default: "Delivery")
#' @param colors Named vector of colors for transitions. Default uses:
#'   \itemize{
#'     \item "Gained" = steelblue/cornflowerblue (blue)
#'     \item "Lost" = salmon/indianred (pink/red)
#'     \item "No change (detected)" = gray60
#'     \item "No change (undetected)" = gray90/white
#'   }
#' @param show_labels Logical. Whether to show "Detected" and "Undetected" labels on the sides (default: TRUE)
#' @param width Numeric. Width of the plot in inches (default: 6)
#' @param height Numeric. Height of the plot in inches (default: 8)
#' @param title_size Numeric. Font size for title (default: 14)
#' @param axis_text_size Numeric. Font size for axis text (default: 12)
#' @param save_pdf Logical. Whether to save as PDF (default: FALSE)
#' @param filename Character string. Output filename if save_pdf = TRUE.
#'
#' @return A ggplot2 object showing the alluvial diagram.
#'
#' @details
#' The alluvial plot visualizes the 2x2 contingency table from McNemar's test:
#' \itemize{
#'   \item \strong{Left side (baseline):} Black bar segment = detected, White segment = undetected
#'   \item \strong{Right side (comparison):} Black bar segment = detected, White segment = undetected
#'   \item \strong{Blue flow:} Gained detection (undetected → detected)
#'   \item \strong{Pink/red flow:} Lost detection (detected → undetected)
#'   \item \strong{Gray flows:} No change (stable detection status)
#' }
#'
#' The plot title shows the cytokine name and the delta detection percentage.
#' Black bars are annotated with the percentage of subjects detected.
#'
#' @examples
#' \dontrun{
#' # Run McNemar test
#' results <- mcnemar_detection(data, baseline = "Enrollment", comparison = "Delivery")
#' 
#' # Plot most significant cytokine
#' plot_mcnemar_alluvial(results)
#' 
#' # Plot specific cytokine
#' plot_mcnemar_alluvial(results, cytokine = "IL12p40")
#' 
#' # Customize labels and colors
#' plot_mcnemar_alluvial(results, 
#'                       baseline_label = "Pre", 
#'                       comparison_label = "Post",
#'                       colors = c("Gained" = "#4169E1", 
#'                                 "Lost" = "#CD5C5C",
#'                                 "No change (detected)" = "#696969",
#'                                 "No change (undetected)" = "white"))
#' }
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_text scale_fill_manual scale_color_manual 
#' @importFrom ggplot2 labs theme_minimal theme element_text element_blank ggsave coord_cartesian guides guide_legend
#' @importFrom ggalluvial geom_alluvium geom_stratum
plot_mcnemar_alluvial <- function(mcnemar_result,
                                  cytokine = NULL,
                                  baseline_label = "Enrollment",
                                  comparison_label = "Delivery",
                                  colors = NULL,
                                  show_labels = TRUE,
                                  width = 6,
                                  height = 8,
                                  title_size = 14,
                                  axis_text_size = 12,
                                  save_pdf = FALSE,
                                  filename = NULL) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.", call. = FALSE)
  if (!requireNamespace("ggalluvial", quietly = TRUE)) stop("Package 'ggalluvial' is required.", call. = FALSE)
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.", call. = FALSE)
  if (!requireNamespace("scales", quietly = TRUE)) stop("Package 'scales' is required.", call. = FALSE)

  if (!is.data.frame(mcnemar_result)) stop("'mcnemar_result' must be a data frame", call. = FALSE)

  required_cols <- c("cytokine", "both_detect", "loss", "gain", "neither_detect", "delta_detection")
  missing_cols <- setdiff(required_cols, names(mcnemar_result))
  if (length(missing_cols) > 0) stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)

  # choose cytokine
  if (is.null(cytokine)) {
    if ("q_mcnemar" %in% names(mcnemar_result)) {
      cytokine <- mcnemar_result$cytokine[which.min(mcnemar_result$q_mcnemar)]
    } else if ("p_mcnemar" %in% names(mcnemar_result)) {
      cytokine <- mcnemar_result$cytokine[which.min(mcnemar_result$p_mcnemar)]
    } else cytokine <- mcnemar_result$cytokine[1]
    message("Plotting cytokine: ", cytokine)
  }
  cy_data <- mcnemar_result[mcnemar_result$cytokine == cytokine, ]
  if (nrow(cy_data) == 0) stop("Cytokine '", cytokine, "' not found", call. = FALSE)

  a <- cy_data$both_detect
  b <- cy_data$loss
  c <- cy_data$gain
  d <- cy_data$neither_detect
  delta <- cy_data$delta_detection

  # flow data (keep factor levels)
  flow_data <- data.frame(
    count = c(a, b, c, d),
    baseline = c("Detected", "Detected", "Undetected", "Undetected"),
    comparison = c("Detected", "Undetected", "Detected", "Undetected"),
    transition = c("No change (detected)", "Lost", "Gained", "No change (undetected)"),
    stringsAsFactors = FALSE
  )
  flow_data$baseline   <- factor(flow_data$baseline,   levels = c("Detected", "Undetected"))
  flow_data$comparison <- factor(flow_data$comparison, levels = c("Detected", "Undetected"))
  flow_data$transition <- factor(flow_data$transition,
                                 levels = c("Gained", "Lost", "No change (detected)", "No change (undetected)"))
  flow_data <- flow_data[flow_data$count > 0, ]
  if (nrow(flow_data) == 0) stop("No data to plot (all counts are zero)", call. = FALSE)

  # expand for alluvial
  alluvial_data <- flow_data[rep(seq_len(nrow(flow_data)), flow_data$count), ]
  alluvial_data$id <- seq_len(nrow(alluvial_data))
  alluvial_data$transition <- factor(alluvial_data$transition,
                                     levels = c("Gained", "Lost", "No change (detected)", "No change (undetected)"))

  # colors
  if (is.null(colors)) {
    colors <- c(
      "Gained" = "#6495ED",                    # cornflowerblue
      "Lost" = "#FA8072",                      # salmon
      "No change (detected)" = "#808080",      # gray
      "No change (undetected)" = "#D3D3D3"     # lightgray (darker than whitesmoke for visibility)
    )
  }
  stratum_colors <- c("Detected" = "black", "Undetected" = "white")

  # stars
  sig_stars <- ""
  if ("q_mcnemar" %in% names(cy_data)) {
    q <- cy_data$q_mcnemar
    if (!is.na(q)) sig_stars <- if (q < 0.001) "***" else if (q < 0.01) "**" else if (q < 0.05) "*" else ""
  } else if ("p_mcnemar" %in% names(cy_data)) {
    p <- cy_data$p_mcnemar
    if (!is.na(p)) sig_stars <- if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
  }
  plot_title <- sprintf("%s\n(%+.1f%%)%s", cytokine, delta, ifelse(sig_stars != "", paste0(" ", sig_stars), ""))

  # proportions for annotations
  n_total <- a + b + c + d
  prop_baseline_detected   <- (a + b) / n_total
  prop_comparison_detected <- (a + c) / n_total

  # y positions in data units (0..N), not 0..1
  eps_y      <- 0.02                 # 2% of the bar height for padding
  y_top      <- n_total * (1 - eps_y)
  y_bottom   <- n_total * eps_y
  y_top_blkB <- y_top                # top of black stratum (it is the top of the bar)

  # ---- plot ----
  bar_width <- 1/4
  left_label_x <- 1 - (bar_width/2 + 0.12)  # place further left to avoid bar overlap

  p <- ggplot2::ggplot(alluvial_data,
                       ggplot2::aes(axis1 = baseline, axis2 = comparison, y = 1)) +
    suppressWarnings(
      ggalluvial::geom_alluvium(ggplot2::aes(fill = transition),
                                width = bar_width, alpha = 0.7, color = NA)
    ) +
    ggalluvial::geom_stratum(ggplot2::aes(fill = ggplot2::after_stat(stratum)),
                             width = bar_width, color = "black", size = 1) +
    ggplot2::scale_fill_manual(
      values = c(colors, stratum_colors),
      breaks = c("Detected", "Undetected", "Gained", "Lost", 
                 "No change (detected)", "No change (undetected)"),
      labels = c("Detected", "Undetected", "Proportion Gained", "Proportion Lost",
                 "No change (Detected)", "No change (Undetected)"),
      drop   = FALSE,
      name   = NULL
    ) +
    ggplot2::scale_x_discrete(
      limits = c("baseline", "comparison"),
      labels = c(baseline_label, comparison_label),
      expand = c(0.35, 0.15)            # extra room on the left for vertical labels
    ) +
    ggplot2::labs(title = plot_title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = axis_text_size) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(hjust = 0.5, size = title_size, face = "bold"),
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = axis_text_size, face = "bold"),
      panel.grid  = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 10),
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 10)  # left margin for labels
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::guides(fill = ggplot2::guide_legend(
      override.aes = list(alpha = 1, color = c("black", "black", rep(NA, 4))),
      nrow = 3,
      byrow = TRUE
    ))

  # vertical side labels (left bar only), at fixed positions
  if (show_labels) {
    # Use fixed positions that don't change with bar heights
    # "Detected" at 75% height, "Undetected" at 25% height
    y_detected_label   <- n_total * 0.75    # fixed at 75% of total height
    y_undetected_label <- n_total * 0.25    # fixed at 25% of total height
    
    side_lab <- data.frame(
      x     = c(left_label_x, left_label_x),
      y     = c(y_detected_label, y_undetected_label),  # fixed positions
      label = c("Detected", "Undetected"),
      vjust = c(0.5, 0.5)  # Center labels vertically
    )
    if (prop_baseline_detected <= 0) side_lab <- subset(side_lab, label != "Detected")
    if (prop_baseline_detected >= 1) side_lab <- subset(side_lab, label != "Undetected")

    p <- p + ggplot2::geom_text(
      data = side_lab,
      ggplot2::aes(x = x, y = y, label = label, vjust = vjust),
      angle = 90, hjust = 0.5, size = 3.2, color = "black",
      inherit.aes = FALSE
    )
  }

  # % detected inside the black strata (baseline and comparison)
  # Top-aligned annotations inside the black bars
  # Split percentages into two lines: number on top, % symbol on bottom
  pct_vals <- c(prop_baseline_detected, prop_comparison_detected)
  pct_labels <- sprintf("%.1f\n%%", pct_vals * 100)  # e.g., "92.7\n%"
  
  pct_lab <- data.frame(
    x     = c(1, 2),
    y     = c(y_top_blkB, y_top_blkB),    # place at the very top of each black bar
    label = pct_labels
  )
  pct_lab <- pct_lab[c(prop_baseline_detected > 0, prop_comparison_detected > 0), , drop = FALSE]

  p <- p + ggplot2::geom_text(
    data = pct_lab,
    ggplot2::aes(x = x, y = y, label = label),
    vjust = 1.1,                          # nudge inside the stratum
    color = "white", fontface = "bold", size = 3.0, lineheight = 0.85,
    inherit.aes = FALSE
  )

  if (save_pdf) {
    if (is.null(filename)) filename <- paste0("mcnemar_alluvial_", gsub("[^A-Za-z0-9]", "_", cytokine), ".pdf")
    ggplot2::ggsave(filename, plot = p, width = width, height = height)
    message("Plot saved to: ", filename)
  }
  p
}


#' Create Multiple Alluvial Plots for Top Cytokines (single shared legend using facets)
#'
#' Generates a grid of alluvial plots for multiple cytokines using facet_wrap,
#' which guarantees a single shared legend.
#'
#' @param mcnemar_result A data frame or tibble from \code{mcnemar_detection()}.
#' @param n_cytokines Integer. Number of top cytokines to plot (default: 6).
#' @param sort_by Character string. How to sort cytokines: "significance" (default),
#'   "delta_abs" (absolute change), or "delta" (signed change).
#' @param ncol Integer. Number of columns in the grid (default: 2).
#' @param baseline_label Character string. Label for baseline timepoint.
#' @param comparison_label Character string. Label for comparison timepoint.
#' @param colors Named vector of colors for transitions.
#' @param show_labels Logical. Whether to show "Detected" and "Undetected" labels.
#' @param sig_threshold Numeric. FDR threshold for significance highlighting (default: 0.05).
#'   Significant cytokines will be marked with a filled star symbol.
#' @param save_pdf Logical. Whether to save as PDF (default: FALSE).
#' @param filename Character string. Output filename if save_pdf = TRUE.
#' @param width Numeric. Total width in inches (default: 12).
#' @param height Numeric. Total height in inches (default: 4 * ceiling(n_cytokines/ncol)).
#'
#' @return A ggplot2 object with facets.
#'
#' @examples
#' \dontrun{
#' results <- mcnemar_detection(data)
#' plot_mcnemar_alluvial_grid(results, n_cytokines = 6)
#' }
#'
#' @export
#' @importFrom ggplot2 ggsave facet_wrap
plot_mcnemar_alluvial_grid <- function(mcnemar_result,
                                       n_cytokines = 6,
                                       sort_by = c("significance", "delta_abs", "delta"),
                                       ncol = 2,
                                       baseline_label = "Enrollment",
                                       comparison_label = "Delivery",
                                       colors = NULL,
                                       show_labels = TRUE,
                                       sig_threshold = 0.05,
                                       save_pdf = FALSE,
                                       filename = NULL,
                                       width = 12,
                                       height = NULL) {

  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required. Please install it.", call. = FALSE)
  }

  sort_by <- match.arg(sort_by)
  if (sort_by == "significance") {
    if ("q_mcnemar" %in% names(mcnemar_result)) {
      mcnemar_result <- mcnemar_result[order(mcnemar_result$q_mcnemar), ]
    } else if ("p_mcnemar" %in% names(mcnemar_result)) {
      mcnemar_result <- mcnemar_result[order(mcnemar_result$p_mcnemar), ]
    }
  } else if (sort_by == "delta_abs") {
    mcnemar_result <- mcnemar_result[order(abs(mcnemar_result$delta_detection), decreasing = TRUE), ]
  } else if (sort_by == "delta") {
    mcnemar_result <- mcnemar_result[order(mcnemar_result$delta_detection, decreasing = TRUE), ]
  }

  n_cytokines <- min(n_cytokines, nrow(mcnemar_result))
  top_cytokines <- mcnemar_result$cytokine[1:n_cytokines]
  message("Creating faceted alluvial plot for ", n_cytokines, " cytokines...")

  # Build combined data for all cytokines
  all_data_list <- lapply(seq_along(top_cytokines), function(i) {
    cy <- top_cytokines[i]
    cy_data <- mcnemar_result[mcnemar_result$cytokine == cy, ]
    
    a <- cy_data$both_detect
    b <- cy_data$loss
    c <- cy_data$gain
    d <- cy_data$neither_detect
    delta <- cy_data$delta_detection
    
    # Significance indicator - use filled star for significant cytokines
    sig_indicator <- ""
    is_significant <- FALSE
    if ("q_mcnemar" %in% names(cy_data)) {
      q <- cy_data$q_mcnemar
      if (!is.na(q) && q < sig_threshold) {
        is_significant <- TRUE
        # Use filled star symbol for FDR-significant
        sig_indicator <- "\u2605 "  # ★ (filled star)
      }
    } else if ("p_mcnemar" %in% names(cy_data)) {
      p <- cy_data$p_mcnemar
      if (!is.na(p) && p < sig_threshold) {
        is_significant <- TRUE
        sig_indicator <- "\u2605 "  # ★ (filled star)
      }
    }
    
    # Create flow data
    flow_data <- data.frame(
      count = c(a, b, c, d),
      baseline = c("Detected", "Detected", "Undetected", "Undetected"),
      comparison = c("Detected", "Undetected", "Detected", "Undetected"),
      transition = c("No change (detected)", "Lost", "Gained", "No change (undetected)"),
      stringsAsFactors = FALSE
    )
    flow_data$baseline <- factor(flow_data$baseline, levels = c("Detected", "Undetected"))
    flow_data$comparison <- factor(flow_data$comparison, levels = c("Detected", "Undetected"))
    flow_data$transition <- factor(flow_data$transition,
                                   levels = c("Gained", "Lost", "No change (detected)", "No change (undetected)"))
    flow_data <- flow_data[flow_data$count > 0, ]
    
    # Expand
    alluvial_data <- flow_data[rep(seq_len(nrow(flow_data)), flow_data$count), ]
    alluvial_data$id <- seq_len(nrow(alluvial_data))
    alluvial_data$transition <- factor(alluvial_data$transition,
                                       levels = c("Gained", "Lost", "No change (detected)", "No change (undetected)"))
    
    # Add cytokine info with significance indicator prefix
    facet_label <- sprintf("%s%s\n(%+.1f%%)", sig_indicator, cy, delta)
    alluvial_data$cytokine <- cy
    alluvial_data$is_significant <- is_significant
    alluvial_data$facet_label <- facet_label
    alluvial_data$facet_order <- i
    
    # Add proportions for labels
    n_total <- a + b + c + d
    alluvial_data$prop_baseline_detected <- (a + b) / n_total
    alluvial_data$prop_comparison_detected <- (a + c) / n_total
    alluvial_data$n_total <- n_total
    
    alluvial_data
  })
  
  combined_data <- do.call(rbind, all_data_list)
  combined_data$facet_label <- factor(combined_data$facet_label, 
                                      levels = unique(combined_data$facet_label[order(combined_data$facet_order)]))
  
  # Colors
  if (is.null(colors)) {
    colors <- c(
      "Gained" = "#6495ED",
      "Lost" = "#FA8072",
      "No change (detected)" = "#808080",
      "No change (undetected)" = "#D3D3D3"
    )
  }
  stratum_colors <- c("Detected" = "black", "Undetected" = "white")
  
  bar_width <- 1/4
  left_label_x <- 1 - (bar_width/2 + 0.12)
  
  # Create side labels data (fixed positions per facet)
  if (show_labels) {
    side_labels_list <- lapply(unique(combined_data$facet_label), function(fl) {
      subset_data <- combined_data[combined_data$facet_label == fl, ][1, ]
      n_total <- subset_data$n_total
      data.frame(
        facet_label = fl,
        x = c(left_label_x, left_label_x),
        y = c(n_total * 0.75, n_total * 0.25),
        label = c("Detected", "Undetected")
      )
    })
    side_labels <- do.call(rbind, side_labels_list)
  }
  
  # Create percentage labels data
  pct_labels_list <- lapply(unique(combined_data$facet_label), function(fl) {
    subset_data <- combined_data[combined_data$facet_label == fl, ][1, ]
    n_total <- subset_data$n_total
    prop_baseline_detected <- subset_data$prop_baseline_detected
    prop_comparison_detected <- subset_data$prop_comparison_detected
    
    data.frame(
      facet_label = fl,
      x = c(1, 2),
      y = c(n_total * (1 - 0.02), n_total * (1 - 0.02)),
      label = sprintf("%.1f\n%%", c(prop_baseline_detected, prop_comparison_detected) * 100)
    )
  })
  pct_labels <- do.call(rbind, pct_labels_list)
  
  # Create plot
  p <- ggplot2::ggplot(combined_data,
                       ggplot2::aes(axis1 = baseline, axis2 = comparison, y = 1)) +
    suppressWarnings(
      ggalluvial::geom_alluvium(ggplot2::aes(fill = transition),
                                width = bar_width, alpha = 0.7, color = NA)
    ) +
    ggalluvial::geom_stratum(ggplot2::aes(fill = ggplot2::after_stat(stratum)),
                             width = bar_width, color = "black", size = 1) +
    ggplot2::scale_fill_manual(
      values = c(colors, stratum_colors),
      breaks = c("Detected", "Undetected", "Gained", "Lost", 
                 "No change (detected)", "No change (undetected)"),
      labels = c("Detected", "Undetected", "Proportion Gained", "Proportion Lost",
                 "No change (Detected)", "No change (Undetected)"),
      drop = FALSE,
      name = NULL
    ) +
    ggplot2::scale_x_discrete(
      limits = c("baseline", "comparison"),
      labels = c(baseline_label, comparison_label),
      expand = c(0.35, 0.15)
    ) +
    ggplot2::facet_wrap(~ facet_label, ncol = ncol, scales = "free_y") +
    ggplot2::labs(x = NULL, y = NULL,
                  caption = sprintf("\u2605 = FDR-significant (q < %.2f)", sig_threshold)) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 9, face = "bold"),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 10),
      strip.text = ggplot2::element_text(size = 11, face = "bold"),
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 10),
      plot.caption = ggplot2::element_text(hjust = 0, size = 10, face = "bold")
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::guides(fill = ggplot2::guide_legend(
      override.aes = list(alpha = 1, color = c("black", "black", rep(NA, 4)))
    ))
  
  # Add side labels
  if (show_labels) {
    p <- p + ggplot2::geom_text(
      data = side_labels,
      ggplot2::aes(x = x, y = y, label = label),
      angle = 90, hjust = 0.5, vjust = 0.5, size = 2.8, color = "black",
      inherit.aes = FALSE
    )
  }
  
  # Add percentage labels
  p <- p + ggplot2::geom_text(
    data = pct_labels,
    ggplot2::aes(x = x, y = y, label = label),
    vjust = 1.1, color = "white", fontface = "bold", size = 2.6, lineheight = 0.85,
    inherit.aes = FALSE
  )

  if (is.null(height)) height <- 4 * ceiling(n_cytokines / ncol)

  if (save_pdf) {
    if (is.null(filename)) filename <- "mcnemar_alluvial_grid.pdf"
    ggplot2::ggsave(filename, plot = p, width = width, height = height)
    message("Grid plot saved to: ", filename)
  }
  p
}

