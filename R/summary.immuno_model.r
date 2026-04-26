# S3 summary method for immuno_model objects
# Legacy compatibility wrapper for immuno_fit summary functionality

#' Summarize an immuno_model object
#'
#' @description
#' S3 summary method for `immuno_model` objects. Currently an alias 
#' for `summary.immuno_fit` since immuno_model is the legacy name.
#'
#' @param object An `immuno_model` object (legacy name for immuno_fit).
#' @param ... Additional arguments passed to summary.immuno_fit.
#'
#' @return Summary information for the model.
#'
#' @examples
#' \dontrun{
#' mdl <- fit_one(data, family = "gamma")
#' summary(mdl)
#' }
#' @export
summary.immuno_model <- function(object, ...) {
  # Legacy support - immuno_model is the old name for immuno_fit
  if (inherits(object, "immuno_fit")) {
    summary.immuno_fit(object, ...)
  } else {
    stop("Object must be of class 'immuno_fit' or 'immuno_model'")
  }
}
