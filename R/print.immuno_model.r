# S3 print method for immuno_model objects
# Provides minimal object representation for interactive debugging

#' @title Print method for immuno_model
#' @description Minimal representation for interactive use.
#' @param x   An `immuno_model` object.
#' @param ... Unused.
#' @return `x` invisibly.
#' @export
print.immuno_model <- function(x, ...) {
  cat("<immuno_model>\n")
  cat("Fields:", paste(names(x), collapse = ", "), "\n")
  invisible(x)
}
