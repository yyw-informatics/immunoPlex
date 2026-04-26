# Cytokine Data Preprocessing Module for immunoPlex Package
# Comprehensive preprocessing functionality for cytokine expression data with flexible LOD handling

#' Preprocess cytokine expression data with flexible LOD handling and replicate aggregation
#'
#' This function performs comprehensive preprocessing of cytokine expression data,
#' including data cleaning, reshaping from wide to long format, flexible
#' handling of values below the limit of detection (LOD), and optional aggregation
#' of technical replicates. Multiple LOD methods can be applied and compared simultaneously.
#' @importFrom dplyr group_by summarise mutate case_when across all_of first rename
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom rlang sym
#' @importFrom magrittr %>%
#' @importFrom stats runif
#' @param expr data.frame of raw cytokine measurements where rows represent samples
#'   and columns represent different cytokines. Sample identifiers can be provided either:
#'   (1) as rownames (traditional format), or (2) as a 'sample_id' column (modern format).
#'   For longitudinal data with repeated subjects, create unique IDs combining subject and timepoint
#'   (e.g., paste(subject_id, timepoint, sep="_")).
#' @param meta data.frame containing sample metadata. Must include a column that
#'   matches the sample identifiers from expr (typically "sample_id").
#' @param lod_lookup data.frame with two required columns: 'cytokine' (matching column
#'   names in expr) and 'lod' (numeric limit of detection values for each cytokine).
#' @param subject_id_col character string specifying the column name in meta that contains
#'   subject identifiers for grouping technical replicates. If NULL, no replicate 
#'   aggregation is available. Default: NULL.
#' @param aggregate_replicates logical indicating whether to aggregate technical replicates
#'   by averaging. If TRUE, requires subject_id_col to be specified. Aggregation occurs
#'   after LOD substitution. Default: FALSE.
#' @param lod_methods character vector specifying LOD substitution methods to apply.
#'   Options: "half" (LOD/2), "zero" (0), "sqrt" (square root of LOD), 
#'   "lod" (LOD value), "uniform" (random uniform between 0 and LOD),
#'   "halfmin" (half the minimum positive observed value). Default: "half".
#' @param replicate_qc logical indicating whether to run replicate quality control
#'   via \code{\link{flag_replicate_outliers}} before aggregation. Requires
#'   \code{subject_id_col} to identify biological samples. Default: FALSE.
#' @param replicate_qc_args named list of additional arguments passed to
#'   \code{\link{flag_replicate_outliers}} (e.g., \code{cv_threshold},
#'   \code{action}). Default: \code{list()}.
#' @param censor_flag logical indicating whether to add a 'censored' column flagging
#'   values that were below LOD. Default: TRUE.
#'
#' @return If single method: data.frame with class "immuno_preprocess" in long format
#'   containing columns for sample_id (contains subject IDs if aggregated), cytokine, 
#'   concentration, metadata columns, and optionally censored flag. If multiple methods: 
#'   list of such data.frames with class "immuno_preprocess_compare" for easy comparison.
#' 
#' @details
#' The function performs the following steps:
#' 1. Validates input data structures and consistency
#' 2. Merges expression data with metadata 
#' 3. Reshapes data from wide to long format
#' 4. Applies specified LOD substitution method(s)
#' 5. Optionally aggregates technical replicates by subject
#' 6. Ensures no remaining zeros when using log transforms via half-minimum method
#' 7. Optionally adds censoring flags
#' 8. Attaches appropriate classes and metadata attributes
#' 
#' Technical replicates are handled as follows:
#' - LOD substitution is applied at the individual replicate level
#' - If aggregate_replicates=TRUE, replicates are averaged after LOD substitution
#' - Censoring flags reflect the final values (post-aggregation if applicable)
#' - When aggregating, metadata is preserved for the first replicate of each subject
#' 
#' LOD substitution methods available:
#' \itemize{
#'   \item "half": Replace with LOD/2 (most common, assumes log-normal distribution)
#'   \item "zero": Replace with 0 (conservative, may bias toward null)
#'   \item "sqrt": Replace with sqrt(LOD) (intermediate approach)
#'   \item "lod": Replace with LOD value (sets to detection limit)
#'   \item "uniform": Replace with random uniform(0, LOD) (adds variability)
#'   \item "halfmin": Replace with half the minimum positive observed value across all data (useful for log transformations)
#' }
#' 
#' @examples
#' \dontrun{
#' # Load example data
#' data("immunoplex_example", package = "immunoPlex")
#' 
#' # Single method preprocessing (traditional format with rownames)
#' processed_df <- immuno_preprocess(
#'   expr = immunoplex_example$expression,
#'   meta = immunoplex_example$metadata,
#'   lod_lookup = immunoplex_example$lod_lookup,
#'   lod_methods = "half"
#' )
#' 
#' # Modern format with sample_id column
#' expr_with_id <- immunoplex_example$expression %>%
#'   tibble::rownames_to_column("sample_id")
#' processed_df <- immuno_preprocess(
#'   expr = expr_with_id,
#'   meta = immunoplex_example$metadata,
#'   lod_lookup = immunoplex_example$lod_lookup,
#'   lod_methods = "half"
#' )
#' 
#' # Longitudinal data: create unique IDs for repeated subjects
#' # Assuming 'wide_data' has subject_id and timepoint columns
#' expr_long <- wide_data %>% select(cytokine_columns_only)
#' rownames(expr_long) <- paste(wide_data$subject_id, wide_data$timepoint, sep="_")
#' meta_long <- wide_data %>% select(metadata_columns)
#' meta_long$sample_id <- paste(wide_data$subject_id, wide_data$timepoint, sep="_")
#' 
#' processed_df <- immuno_preprocess(
#'   expr = expr_long,
#'   meta = meta_long,
#'   lod_lookup = lod_lookup,
#'   lod_methods = "half"
#' )
#' 
#' # Compare multiple methods including halfmin for log-transform compatibility
#' processed_df <- immuno_preprocess(
#'   expr = immunoplex_example$expression,
#'   meta = immunoplex_example$metadata,
#'   lod_lookup = immunoplex_example$lod_lookup,
#'   lod_methods = c("half", "halfmin")
#' )
#' 
#' # With technical replicate aggregation
#' processed_df <- immuno_preprocess(
#'   expr = immunoplex_example$expression,
#'   meta = immunoplex_example$metadata,
#'   lod_lookup = immunoplex_example$lod_lookup,
#'   subject_id_col = "subject_id",
#'   aggregate_replicates = TRUE,
#'   lod_methods = "halfmin"
#' )
#' }
#' 
#' @export
immuno_preprocess <- function(expr, meta, lod_lookup,
                               subject_id_col = NULL,
                               aggregate_replicates = FALSE,
                               lod_methods = "half",
                               replicate_qc = FALSE,
                               replicate_qc_args = list(),
                               censor_flag = TRUE) {
  
  # ===== SECTION 1: INPUT VALIDATION =====
  # Validate basic data structure requirements
  stopifnot(is.data.frame(expr), is.data.frame(meta), is.data.frame(lod_lookup))
  
  # Validate required columns efficiently
  if (!all(c("cytokine", "lod") %in% names(lod_lookup))) {
    stop("lod_lookup must contain columns 'cytokine' and 'lod'")
  }
  
  # Ensure lod column is numeric (may be character in CSV files)
  lod_lookup$lod <- suppressWarnings(as.numeric(as.character(lod_lookup$lod)))
  
  if (!"sample_id" %in% names(meta)) {
    stop("meta data.frame must contain a 'sample_id' column")
  }
  
  # Validate replicate aggregation parameters
  if (aggregate_replicates) {
    if (is.null(subject_id_col)) {
      stop("aggregate_replicates=TRUE requires subject_id_col to be specified")
    }
    if (!subject_id_col %in% names(meta)) {
      stop("subject_id_col '", subject_id_col, "' not found in meta data.frame")
    }
  }
  
  # Handle expr input format (with or without sample_id column)
  expr_has_sample_id <- "sample_id" %in% names(expr)
  
  # Check for cytokines without LOD values (warn but allow processing)
  expr_cytokines <- if (expr_has_sample_id) {
    setdiff(names(expr), "sample_id")
  } else {
    names(expr)
  }
  
  missing_lods <- setdiff(expr_cytokines, lod_lookup$cytokine)
  if (length(missing_lods) > 0) {
    message("Note: ", length(missing_lods), " cytokine(s) without LOD values will be processed without substitution: ",
            paste(missing_lods, collapse = ", "))
    # Add entries for missing cytokines with NA LOD values
    missing_lod_entries <- data.frame(
      cytokine = missing_lods,
      lod = NA_real_,
      stringsAsFactors = FALSE
    )
    lod_lookup <- rbind(lod_lookup, missing_lod_entries)
  }

  # ===== SECTION 2: DATA PREPARATION =====
  
  # Handle different expr input formats robustly
  if (expr_has_sample_id) {
    # Case 1: expr already has sample_id column (modern format)
    merged <- expr %>%
      dplyr::left_join(meta, by = "sample_id")
  } else {
    # Case 2: expr has sample_ids as rownames (traditional format)
    # Check for valid rownames
    if (is.null(rownames(expr)) || all(rownames(expr) == as.character(1:nrow(expr)))) {
      stop("expr must have either:\n",
           "  1) A 'sample_id' column, OR\n",
           "  2) Sample identifiers as rownames (not just 1, 2, 3...)\n",
           "Hint: For longitudinal data, create unique IDs like paste(subject_id, timepoint, sep='_')")
    }
    
    # Check for duplicate rownames
    if (anyDuplicated(rownames(expr))) {
      dup_ids <- rownames(expr)[duplicated(rownames(expr))]
      stop("Duplicate rownames detected in expr: ", 
           paste(head(dup_ids, 5), collapse = ", "),
           if (length(dup_ids) > 5) " ... and more" else "",
           "\nFor longitudinal data, create unique IDs combining subject and timepoint.",
           "\nExample: rownames(expr) <- paste(subject_id, timepoint, sep='_')")
    }
    
    merged <- expr %>%
      tibble::rownames_to_column("sample_id") %>%
      dplyr::left_join(meta, by = "sample_id")
  }
  
  # Convert all cytokine columns to numeric to handle mixed types
  # (Some columns may be character due to "<LOD" notation, others numeric)
  # Track censoring information before converting to numeric
  cytokine_cols_in_merged <- setdiff(names(merged), c("sample_id", names(meta)))
  
  # Create a parallel data frame to track censoring flags
  censoring_flags <- data.frame(sample_id = merged$sample_id)
  
  for (col in cytokine_cols_in_merged) {
    if (is.character(merged[[col]])) {
      # Detect censoring markers before removing them
      censoring_flags[[paste0(col, "_censored_below")]] <- grepl("^<", merged[[col]], perl = TRUE)
      censoring_flags[[paste0(col, "_censored_above")]] <- grepl("^>", merged[[col]], perl = TRUE)
      # Extract numeric values from strings like "<0.5" or ">1000"
      merged[[col]] <- suppressWarnings(as.numeric(gsub("^[<>]", "", merged[[col]])))
    } else {
      # For numeric columns, no censoring markers present
      censoring_flags[[paste0(col, "_censored_below")]] <- FALSE
      censoring_flags[[paste0(col, "_censored_above")]] <- FALSE
    }
  }
  
  # Pivot both data and censoring flags to long format
  long_df <- merged %>%
    tidyr::pivot_longer(-c(sample_id, dplyr::one_of(names(meta))),
                        names_to = "cytokine", values_to = "concentration") %>%
    dplyr::left_join(lod_lookup, by = "cytokine")
  
  # Pivot censoring flags
  censoring_long <- censoring_flags %>%
    tidyr::pivot_longer(-sample_id,
                        names_to = c("cytokine", "flag_type"),
                        names_pattern = "(.+)_censored_(below|above)",
                        values_to = "flag_value") %>%
    tidyr::pivot_wider(names_from = flag_type, values_from = flag_value,
                       names_prefix = "censored_")
  
  # Join censoring information with concentration data
  long_df <- long_df %>%
    dplyr::left_join(censoring_long, by = c("sample_id", "cytokine"))

  results <- vector("list", length(lod_methods))
  names(results) <- lod_methods

  # Handle cases where censored_below might not exist (if all columns were numeric)
  if (!"censored_below" %in% names(long_df)) {
    long_df$censored_below <- FALSE
    long_df$censored_above <- FALSE
  }
  
  # Pre-compute censoring status before LOD substitution
  # A value is censored if it's explicitly flagged OR below the LOD threshold
  long_df <- long_df %>%
    dplyr::mutate(
      # Ensure all components handle NA properly
      censored_below = tidyr::replace_na(censored_below, FALSE),
      censored_above = tidyr::replace_na(censored_above, FALSE),
      was_censored = !is.na(lod) & !is.na(concentration) & 
                     (censored_below | concentration < lod)
    )
  
  for (method in lod_methods) {
    if (method == "uniform") {
      df <- long_df %>% dplyr::rowwise() %>%
        dplyr::mutate(
          concentration = ifelse(was_censored, runif(1, 0, lod), concentration)
        ) %>%
        dplyr::ungroup()

    } else if (method == "halfmin") {
      min_pos <- min(long_df$concentration[long_df$concentration > 0], na.rm = TRUE)
      df <- long_df %>%
        dplyr::mutate(
          concentration = ifelse(was_censored, min_pos/2, concentration)
        )

    } else {
      # Determine substitution value based on method (scalar, not vectorized)
      # Store original concentration before modification
      df <- long_df %>%
        dplyr::mutate(original_conc = concentration)
      
      if (method == "half") {
        df <- df %>%
          dplyr::mutate(concentration = dplyr::if_else(was_censored & !is.na(lod), lod/2, original_conc))
      } else if (method == "zero") {
        df <- df %>%
          dplyr::mutate(concentration = dplyr::if_else(was_censored, 0, original_conc))
      } else if (method == "sqrt") {
        df <- df %>%
          dplyr::mutate(concentration = dplyr::if_else(was_censored & !is.na(lod), sqrt(lod), original_conc))
      } else if (method == "lod") {
        df <- df %>%
          dplyr::mutate(concentration = dplyr::if_else(was_censored & !is.na(lod), lod, original_conc))
      }
      
      df <- df %>% dplyr::select(-original_conc)
    }

    # Replicate QC: flag/handle discordant replicates before aggregation
    if (replicate_qc && !is.null(subject_id_col)) {
      df <- do.call(flag_replicate_outliers, c(
        list(data = df, sample_id_col = subject_id_col,
             cytokine_col = "cytokine", value_col = "concentration"),
        replicate_qc_args))
    }

    if (aggregate_replicates) {
      df <- df %>%
        dplyr::group_by(across(all_of(c(subject_id_col, "cytokine")))) %>%
        dplyr::summarise(
          concentration = mean(concentration, na.rm = TRUE),
          lod           = dplyr::first(lod),
          dplyr::across(dplyr::all_of(setdiff(names(meta),
                          c("sample_id", subject_id_col))), dplyr::first),
          .groups = "drop"
        ) %>%
        dplyr::rename(sample_id = !!rlang::sym(subject_id_col))
    }

    if (censor_flag) {
      # Use the pre-computed censoring status
      df <- df %>% dplyr::mutate(censored = was_censored)
    }
    
    # Clean up temporary columns
    df <- df %>% dplyr::select(-censored_below, -censored_above, -was_censored)

    class(df) <- c("immuno_preprocess", "data.frame")
    attr(df, "lod_method")        <- method
    attr(df, "percent_censored")  <- if (censor_flag) mean(df$censored)*100 else NA_real_
    attr(df, "aggregated_replicates") <- aggregate_replicates
    if (aggregate_replicates) attr(df, "subject_id_col") <- subject_id_col

    results[[method]] <- df
  }

  if (length(results) == 1) return(results[[1]])
  class(results) <- "immuno_preprocess_compare"
  results
}