#' Collapse technical replicates to one row per (subject, cytokine, timepoint)
#'
#' Utility for pre-aggregating long-format multiplex data before passing it to
#' downstream functions whose statistical model cannot consume replicate-level
#' rows (e.g., \code{\link{plsda_fit}}, \code{\link{mcnemar_detection}}). For
#' functions that \emph{can} consume replicate-level rows via random effects
#' (\code{\link{fit_models}}, \code{\link{ancova_one}}), leave replicates in
#' place and let the model handle them.
#'
#' @param data A data frame in long format with one row per replicate
#'   observation.
#' @param rep_col Character string naming the column that identifies the
#'   replicate within a \code{(subject, cytokine, timepoint)} group (e.g.,
#'   \code{"rep_id"}, \code{"well"}). Required; used for validation and
#'   dropped from the output. Not used to define the aggregation group
#'   (that is fixed by \code{subject_col}, \code{cytokine_col}, and
#'   \code{timepoint_col}).
#' @param subject_col Character string naming the subject identifier column.
#' @param cytokine_col Character string naming the cytokine/analyte column.
#' @param timepoint_col Character string naming the timepoint/condition
#'   column, or \code{NULL} for cross-sectional data. When \code{NULL}, the
#'   grouping key is \code{(subject, cytokine)} only. Default: \code{NULL}.
#' @param value_col Character string naming a continuous (numeric) column to
#'   aggregate. Mutually exclusive with \code{detection_col} and
#'   \code{censoring_col}. Default: \code{NULL}.
#' @param detection_col Character string naming a logical detection column
#'   (\code{TRUE} = detected, \code{FALSE} = below LOD). Mutually exclusive
#'   with \code{value_col} and \code{censoring_col}. Default: \code{NULL}.
#' @param censoring_col Character string naming a logical censoring column
#'   (\code{TRUE} = censored / below LOD). Aggregation is applied to
#'   \code{detected = !censored} using the same rule set as
#'   \code{detection_col}; the output column is renamed back to censoring
#'   semantics (\code{TRUE} = censored). Mutually exclusive with
#'   \code{value_col} and \code{detection_col}. Default: \code{NULL}.
#' @param rule Character string selecting the aggregation rule. Valid
#'   choices depend on the column type:
#'   \itemize{
#'     \item Continuous (\code{value_col}):
#'       \code{"mean"} (default), \code{"median"}, \code{"geomean"}.
#'       \code{"geomean"} requires strictly positive values; replicates with
#'       any non-positive value collapse to \code{NA} with a warning.
#'     \item Boolean (\code{detection_col} or \code{censoring_col}):
#'       \code{"majority_vote"} (default), \code{"any_detected"},
#'       \code{"all_detected"}. Semantics match
#'       \code{\link{mcnemar_detection}}'s \code{replicate_agg} argument:
#'       \code{"majority_vote"} uses
#'       \code{as.logical(round(mean(x, na.rm = TRUE)))} (R's banker's
#'       rounding sends a 50/50 tie to \code{FALSE});
#'       \code{"any_detected"} is \code{any(x)} restricted to non-NA values;
#'       \code{"all_detected"} is \code{all(x)} restricted to non-NA values.
#'       All three rules return \code{NA} for a group with zero non-NA
#'       replicates, rather than returning an informative default from the
#'       empty set (\code{any(logical(0))} is \code{FALSE},
#'       \code{all(logical(0))} is \code{TRUE}).
#'   }
#'   If \code{rule = NULL}, the type-appropriate default is used.
#' @param keep_sd Logical. If \code{TRUE} and \code{value_col} is set, adds a
#'   companion \code{<value_col>_sd} column carrying the within-replicate
#'   standard deviation (\code{NA} for single-replicate groups). Meaningless
#'   for boolean paths and rejected there. Default: \code{FALSE}.
#'
#' @return A data frame with one row per \code{(subject, cytokine,
#'   [timepoint])} group containing the group keys and the aggregated
#'   column (and \code{<value_col>_sd} when \code{keep_sd = TRUE}). Other
#'   input columns are not preserved; join on the group keys if you need to
#'   re-attach metadata.
#'
#' @details The degenerate case of a single observation per group collapses
#'   all rules to the identity. NA handling follows \code{na.rm = TRUE}
#'   throughout, with \code{"all_detected"} specifically guarding against
#'   the \code{all(logical(0))} -> \code{TRUE} pitfall by returning
#'   \code{NA} for fully-missing groups.
#'
#' @examples
#' \dontrun{
#' # Continuous data: mean-collapse three wells per sample
#' agg <- aggregate_replicates(
#'   data          = luminex_long,
#'   rep_col       = "well",
#'   subject_col   = "subject_id",
#'   cytokine_col  = "cytokine",
#'   timepoint_col = "visit",
#'   value_col     = "intensity",
#'   rule          = "mean",
#'   keep_sd       = TRUE
#' )
#'
#' # Boolean detection data: any-detected across replicates
#' agg_det <- aggregate_replicates(
#'   data          = detection_long,
#'   rep_col       = "rep_id",
#'   subject_col   = "subject_id",
#'   cytokine_col  = "cytokine",
#'   timepoint_col = "visit",
#'   detection_col = "detected",
#'   rule          = "any_detected"
#' )
#' }
#'
#' @seealso \code{\link{mcnemar_detection}} (replicate_agg argument),
#'   \code{\link{plsda_fit}} (rep_col argument).
#' @importFrom stats sd
#' @export
aggregate_replicates <- function(data,
                                 rep_col,
                                 subject_col,
                                 cytokine_col,
                                 timepoint_col = NULL,
                                 value_col = NULL,
                                 detection_col = NULL,
                                 censoring_col = NULL,
                                 rule = NULL,
                                 keep_sd = FALSE) {

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame", call. = FALSE)
  }

  # Exactly one target column must be supplied.
  target_flags <- c(
    value     = !is.null(value_col),
    detection = !is.null(detection_col),
    censoring = !is.null(censoring_col)
  )
  if (sum(target_flags) != 1L) {
    stop(
      "Provide exactly one of `value_col`, `detection_col`, or ",
      "`censoring_col` (got ", sum(target_flags), ").",
      call. = FALSE
    )
  }
  mode <- names(target_flags)[target_flags]
  target_col <- switch(mode,
    value     = value_col,
    detection = detection_col,
    censoring = censoring_col
  )

  # Validate presence of all named columns.
  required_cols <- c(rep_col, subject_col, cytokine_col, timepoint_col,
                     target_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns: ",
         paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  # Validate target column type per mode.
  col_vec <- data[[target_col]]
  if (mode == "value") {
    if (!is.numeric(col_vec)) {
      stop("`value_col` ('", target_col, "') must be numeric; got ",
           paste(class(col_vec), collapse = "/"), ".",
           call. = FALSE)
    }
  } else {
    if (!is.logical(col_vec)) {
      stop("`", mode, "_col` ('", target_col,
           "') must be logical; got ",
           paste(class(col_vec), collapse = "/"), ".",
           call. = FALSE)
    }
  }

  # Resolve + validate rule via match.arg so the API matches
  # mcnemar_detection()'s replicate_agg (partial matches accepted, e.g.
  # "maj" -> "majority_vote"). match.arg's built-in error lacks the
  # argument name, so wrap it to preserve the informative "`rule` must
  # be one of ..." message and flag the offending input.
  numeric_rules <- c("mean", "median", "geomean")
  boolean_rules <- c("majority_vote", "any_detected", "all_detected")
  if (mode == "value") {
    allowed <- numeric_rules
    default <- "mean"
  } else {
    allowed <- boolean_rules
    default <- "majority_vote"
  }
  if (is.null(rule)) {
    rule <- default
  }
  rule_raw <- rule
  rule <- tryCatch(match.arg(rule, choices = allowed),
                   error = function(e) NA_character_)
  if (is.na(rule)) {
    stop("`rule` must be one of: ",
         paste(allowed, collapse = ", "),
         " (got '", paste(rule_raw, collapse = ","), "').",
         call. = FALSE)
  }

  if (keep_sd && mode != "value") {
    stop("`keep_sd = TRUE` is only meaningful with `value_col`.",
         call. = FALSE)
  }

  # Build working frame with stable internal column names so we never clobber
  # a user column that happens to be named ".agg_*".
  work <- data.frame(
    .agg_subject  = data[[subject_col]],
    .agg_cytokine = data[[cytokine_col]],
    stringsAsFactors = FALSE
  )
  if (!is.null(timepoint_col)) {
    work$.agg_timepoint <- data[[timepoint_col]]
  }
  # For the censoring path, aggregate on detected = !censored, then invert.
  work$.agg_value <- if (mode == "censoring") !col_vec else col_vec

  # Split into groups.
  grp_keys <- if (is.null(timepoint_col)) {
    list(work$.agg_subject, work$.agg_cytokine)
  } else {
    list(work$.agg_subject, work$.agg_cytokine, work$.agg_timepoint)
  }
  grp_id <- do.call(paste, c(grp_keys, sep = ""))
  split_idx <- split(seq_len(nrow(work)), grp_id)

  # Pre-size output vectors. Use per-group first-row indexing on the original
  # key columns so factor class/levels survive; `vector(mode = typeof(factor))`
  # would allocate an integer vector and strip the labels, sending every
  # factor-keyed group through as its underlying integer code and returning
  # all-NA after the factor re-hydration step below.
  n_out <- length(split_idx)
  first_idx <- vapply(split_idx, `[[`, integer(1L), 1L)
  subject_out   <- work$.agg_subject[first_idx]
  cytokine_out  <- work$.agg_cytokine[first_idx]
  timepoint_out <- if (!is.null(timepoint_col)) {
    work$.agg_timepoint[first_idx]
  } else NULL
  agg_out <- if (mode == "value") {
    numeric(n_out)
  } else {
    logical(n_out)
  }
  sd_out <- if (keep_sd) numeric(n_out) else NULL

  # Accumulate count of geomean groups that hit a non-positive replicate;
  # one summary warning at the end beats N identical warnings on a large
  # panel (cytokines * subjects * timepoints can easily run into the
  # thousands).
  geomean_skipped <- 0L
  agg_fn <- switch(rule,
    mean          = function(x) mean(x, na.rm = TRUE),
    median        = function(x) stats::median(x, na.rm = TRUE),
    geomean       = function(x) {
      xx <- x[!is.na(x)]
      if (length(xx) == 0L) return(NA_real_)
      if (any(xx <= 0)) {
        geomean_skipped <<- geomean_skipped + 1L
        return(NA_real_)
      }
      exp(mean(log(xx)))
    },
    majority_vote = function(x) {
      # mean(NA, na.rm = TRUE) -> NaN -> round -> NaN -> as.logical -> NA,
      # so fully-missing groups already resolve to NA here. The explicit
      # guard mirrors the all_detected branch and makes the intent obvious.
      xx <- x[!is.na(x)]
      if (length(xx) == 0L) NA else as.logical(round(mean(xx)))
    },
    any_detected  = function(x) {
      # Without the guard, any(logical(0), na.rm = TRUE) returns FALSE,
      # which disagrees with majority_vote and all_detected (both return
      # NA for fully-missing groups). Aligning all three rules on
      # "all-NA input -> NA output".
      xx <- x[!is.na(x)]
      if (length(xx) == 0L) NA else any(xx)
    },
    all_detected  = function(x) {
      xx <- x[!is.na(x)]
      if (length(xx) == 0L) NA else all(xx)
    }
  )

  for (i in seq_along(split_idx)) {
    idx <- split_idx[[i]]
    vals <- work$.agg_value[idx]
    agg_out[i] <- agg_fn(vals)
    if (keep_sd) {
      sd_out[i] <- if (length(vals) < 2L) NA_real_ else stats::sd(vals, na.rm = TRUE)
    }
  }

  if (rule == "geomean" && geomean_skipped > 0L) {
    warning(sprintf(
      paste0("geomean requires strictly positive values; ",
             "%d of %d group(s) contained a non-positive replicate ",
             "and were collapsed to NA."),
      geomean_skipped, n_out),
      call. = FALSE)
  }

  # Invert detected -> censored for the censoring path.
  if (mode == "censoring") {
    agg_out <- !agg_out
  }

  out <- data.frame(
    x1 = subject_out,
    x2 = cytokine_out,
    stringsAsFactors = FALSE
  )
  names(out) <- c(subject_col, cytokine_col)
  if (!is.null(timepoint_col)) {
    out[[timepoint_col]] <- timepoint_out
  }
  out[[target_col]] <- agg_out
  if (keep_sd) {
    out[[paste0(target_col, "_sd")]] <- sd_out
  }

  rownames(out) <- NULL
  out
}
