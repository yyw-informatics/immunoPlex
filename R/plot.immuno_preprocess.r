# S3 plotting methods for immuno_preprocess objects
# Visualizes censoring profiles and data preprocessing results

#' Plot percent-censored profile for immunoPlex objects
#'
#' These S3 methods extend **graphics::plot()** so that a simple
#' `plot(obj)` call yields a bar chart of the proportion
#' of censored (<LOD) observations for each cytokine - no faceting by method.
#'
#' @param x   An object of class `immuno_preprocess` or
#'            `immuno_preprocess_compare`.
#' @param ... Ignored; included for S3 consistency.
#' @export
#' @importFrom ggplot2 ggplot aes geom_col coord_flip labs ggsave
#' @importFrom dplyr group_by summarise
#' @importFrom magrittr %>%
#' @method plot immuno_preprocess
plot.immuno_preprocess <- function(x, ...) {
  # If user passed the multi-method list, grab the first element:
  if (inherits(x, "immuno_preprocess_compare")) {
    x <- x[[1]]
  }

  # Summarise percent censored by cytokine
  summary_df <- x %>%
    dplyr::group_by(cytokine) %>%
    dplyr::summarise(
      pct_censored = mean(censored, na.rm = TRUE),
      .groups      = "drop"
    )

  # Build & save the bar chart
  p <- ggplot2::ggplot(summary_df, ggplot2::aes(cytokine, pct_censored)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Percent Censored by Cytokine",
      x     = "Cytokine",
      y     = "Percent Censored"
    )

  out_file <- "cytokine_censoring_profile.pdf"
  ggplot2::ggsave(
    filename = out_file,
    plot     = p,
    width    = 10,
    height   = 6,
    device   = "pdf"
  )
  message("Plot saved as: ", out_file)

  print(p)
  invisible(p)
}

#' @export
#' @method plot immuno_preprocess_compare
plot.immuno_preprocess_compare <- function(x, ...) {
  # Delegate to the immuno_preprocess method
  plot.immuno_preprocess(x, ...)
}
