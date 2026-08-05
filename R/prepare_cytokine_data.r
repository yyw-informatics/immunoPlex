# Data preparation utilities for immunoPlex
# Extracts and formats individual cytokine data for modeling with proper censoring flags

#' Prepare cytokine data for modeling
#'
#' @description
#' Helper function to extract and prepare data for a single cytokine 
#' from the immunoplex_example dataset for use with fit_one(), fit_models(),
#' or compare_lod_models().
#' @importFrom utils data
#'
#' @param cytokine_name Character. Name of the cytokine to extract.
#'   Use names(immunoplex_example$expression) to see available options.
#' @param data_source List. Source data with components 'expression', 
#'   'metadata', and 'lod_lookup'. Defaults to immunoplex_example.
#'
#' @return A data.frame suitable for immunoPlex modeling functions with columns:
#'   - sample_id: Sample identifiers
#'   - subject_id: Subject identifiers
#'   - value: Cytokine concentration values
#'   - timepoint: Timepoint factor
#'   - disease: Disease status factor
#'   - age: Age covariate
#'   - lod: Limit of detection value
#'   - ulod: Upper limit of detection (if available)
#'   - cens_lod: Logical flag for left-censored observations
#'   - cens_ulod: Logical flag for right-censored observations
#'   - any other columns present in `data_source$metadata` (e.g.
#'     `replicate`, `plate`, `preexist`) are passed through, so they
#'     can be supplied to `fit_one()` / `ancova_one()` via `rep_col` /
#'     `plate_col`.
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data("immunoplex_example", package = "immunoPlex")
#' 
#' # Prepare data for one cytokine
#' dat <- prepare_cytokine_data("Cytokine_01")
#' 
#' # Fit a single model
#' fit <- fit_one(dat, family = "gamma")
#' 
#' # Compare multiple models
#' models <- fit_models(dat, families = c("gamma", "tobit", "tobit_censreg"))
#' 
#' # Compare LOD handling approaches
#' comparison <- compare_lod_models(dat, lod_methods = c("half", "halfmin"))
#' }
#' 
#' @export
prepare_cytokine_data <- function(cytokine_name, data_source = NULL) {
  
  # Load default data if not provided
  if (is.null(data_source)) {
    if (!exists("immunoplex_example")) {
      data("immunoplex_example", package = "immunoPlex", envir = environment())
    }
    data_source <- immunoplex_example
  }
  
  # Validate inputs
  if (!cytokine_name %in% names(data_source$expression)) {
    available <- paste(names(data_source$expression), collapse = ", ")
    stop("Cytokine '", cytokine_name, "' not found. Available: ", available)
  }
  
  # Extract components
  raw_expr <- data_source$expression
  meta <- data_source$metadata
  lod_lookup <- data_source$lod_lookup
  
  # Get LOD information for this cytokine
  lod_info <- lod_lookup[lod_lookup$cytokine == cytokine_name, ]
  if (nrow(lod_info) == 0) {
    stop("No LOD information found for cytokine: ", cytokine_name)
  }
  
  target_lod <- lod_info$lod[1]
  target_ulod <- if ("ulod" %in% names(lod_info)) lod_info$ulod[1] else NA
  
  # Build a per-sample frame keyed on sample_id, then left-join all
  # metadata columns so callers get any extra fields they put on
  # `meta` (replicate, plate, preexist, ...) without us having to enumerate.
  per_sample <- data.frame(
    sample_id = rownames(raw_expr),
    value     = raw_expr[, cytokine_name],
    lod       = target_lod,
    ulod      = target_ulod,
    stringsAsFactors = FALSE
  )

  # Drop columns from meta that would collide with what we just built;
  # `sample_id` is the join key.
  reserved <- c("value", "lod", "ulod", "cens_lod", "cens_ulod")
  meta_cols <- setdiff(names(meta), reserved)
  meta_keep <- meta[, meta_cols, drop = FALSE]

  dat <- merge(per_sample, meta_keep, by = "sample_id",
               all.x = TRUE, sort = FALSE)

  # Set censoring flags
  dat$cens_lod  <- dat$value < dat$lod
  dat$cens_ulod <- !is.na(dat$ulod) & dat$value > dat$ulod
  
  # Add attributes for reference
  attr(dat, "cytokine") <- cytokine_name
  attr(dat, "n_censored") <- sum(dat$cens_lod)
  attr(dat, "pct_censored") <- round(100 * mean(dat$cens_lod), 1)
  
  return(dat)
}

#' List available cytokines in the example data
#'
#' @description
#' Convenience function to show available cytokines and their 
#' censoring characteristics.
#'
#' @param data_source List. Source data with immunoplex structure.
#'   Defaults to immunoplex_example.
#'
#' @return A data.frame with cytokine names, LOD values, and censoring info.
#' @export
list_cytokines <- function(data_source = NULL) {
  
  # Load default data if not provided  
  if (is.null(data_source)) {
    if (!exists("immunoplex_example")) {
      data("immunoplex_example", package = "immunoPlex", envir = environment())
    }
    data_source <- immunoplex_example
  }
  
  raw_expr <- data_source$expression
  lod_lookup <- data_source$lod_lookup
  
  # Validate that lod_lookup has required columns
  if (!("cytokine" %in% names(lod_lookup)) || !("lod" %in% names(lod_lookup))) {
    stop("lod_lookup must contain 'cytokine' and 'lod' columns")
  }
  
  # Calculate censoring statistics for each cytokine
  cytokine_info <- data.frame(
    cytokine = names(raw_expr),
    lod = lod_lookup$lod[match(names(raw_expr), lod_lookup$cytokine)],
    n_obs = sapply(raw_expr, function(x) sum(!is.na(x))),
    n_censored = sapply(names(raw_expr), function(cyto) {
      lod_val <- lod_lookup$lod[lod_lookup$cytokine == cyto]
      if (length(lod_val) == 0 || is.na(lod_val)) return(NA)
      sum(raw_expr[, cyto] < lod_val, na.rm = TRUE)
    }),
    stringsAsFactors = FALSE
  )
  
  cytokine_info$pct_censored <- round(100 * cytokine_info$n_censored / cytokine_info$n_obs, 1)
  cytokine_info <- cytokine_info[order(cytokine_info$pct_censored, decreasing = TRUE), ]
  
  return(cytokine_info)
}