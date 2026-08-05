# PLS-DA Model Fitting Module for immunoPlex Package
# Wrapper around ropls package for convenient PLS-DA analysis

#' Recommend cross-validation method based on sample size
#'
#' @description
#' Provides guidance on selecting appropriate cross-validation methods for
#' PLS-DA based on total sample size and minimum group size. This helps users
#' make informed decisions, particularly for small sample studies.
#'
#' @param n_total integer, total number of samples
#' @param n_min_group integer, minimum number of samples in any group
#' @param verbose logical, whether to print recommendation details. Default: FALSE
#' @param loocv_n_threshold integer, threshold for total samples below which LOOCV is recommended. Default: 50
#' @param loocv_group_threshold integer, threshold for minimum group size below which LOOCV is recommended. Default: 15
#' @param moderate_n_threshold integer, threshold for total samples below which moderate k-fold (7-fold) is recommended. Default: 100
#' @param moderate_k integer, number of folds for moderate sample sizes. Default: 7
#' @param large_k integer, number of folds for large sample sizes. Default: 10
#'
#' @return A list containing:
#'   \item{method}{Recommended CV method (integer for k-fold or "loocv")}
#'   \item{description}{Human-readable description}
#'   \item{rationale}{Explanation for the recommendation}
#'   \item{considerations}{Additional considerations for the user}
#'
#' @details
#' Default guidelines used:
#' \itemize{
#'   \item LOOCV is recommended when:
#'     - Total samples < 50 (loocv_n_threshold), OR
#'     - Any group has < 15 samples (loocv_group_threshold)
#'   \item 7-fold CV (moderate_k) is recommended when:
#'     - Total samples 50-100 (loocv_n_threshold to moderate_n_threshold) AND all groups have >= 15 samples
#'   \item 10-fold CV (large_k) is recommended when:
#'     - Total samples > 100 (moderate_n_threshold) AND all groups meet minimum size
#' }
#'
#' These are conservative guidelines based on statistical best practices. LOOCV provides 
#' more stable estimates for small samples but is computationally expensive. k-fold CV is 
#' faster and preferred for larger datasets. All thresholds can be customized via parameters.
#'
#' @examples
#' \dontrun{
#' # Small sample study (uses defaults)
#' recommend_cv_method(n_total = 40, n_min_group = 12)
#' # Recommended: LOOCV
#'
#' # Moderate sample study
#' recommend_cv_method(n_total = 75, n_min_group = 20)
#' # Recommended: 7-fold
#'
#' # Large sample study
#' recommend_cv_method(n_total = 150, n_min_group = 40)
#' # Recommended: 10-fold
#'
#' # Custom thresholds: more conservative LOOCV use
#' recommend_cv_method(n_total = 60, n_min_group = 20, 
#'                     loocv_n_threshold = 75,
#'                     loocv_group_threshold = 20)
#' # Recommended: LOOCV (with custom thresholds)
#' }
#'
#' @export
recommend_cv_method <- function(n_total, n_min_group, verbose = FALSE,
                                loocv_n_threshold = 50,
                                loocv_group_threshold = 15,
                                moderate_n_threshold = 100,
                                moderate_k = 7,
                                large_k = 10) {
  
  # Decision logic based on sample size (using configurable thresholds)
  if (n_total < loocv_n_threshold || n_min_group < loocv_group_threshold) {
    # Small sample: use LOOCV
    method <- "loocv"
    description <- "Leave-One-Out Cross-Validation (LOOCV)"
    rationale <- sprintf(
      "LOOCV recommended: n = %d (total), min group = %d. LOOCV provides stable estimates for small samples.",
      n_total, n_min_group
    )
    considerations <- c(
      "- LOOCV maximizes training data in each fold",
      "- Computationally intensive but provides most stable Q2 estimates",
      "- Less sensitive to random fold assignment",
      if (n_min_group < 10) "! Warning: Very small groups may still yield unstable models"
    )
    
  } else if (n_total < moderate_n_threshold) {
    # Moderate sample: use configurable k-fold (default 7-fold)
    method <- moderate_k
    description <- sprintf("%d-fold Cross-Validation", moderate_k)
    rationale <- sprintf(
      "%d-fold CV recommended: n = %d (total), min group = %d. Good balance of computational efficiency and variance reduction.",
      moderate_k, n_total, n_min_group
    )
    considerations <- c(
      sprintf("- %d-fold is computationally efficient", moderate_k),
      "- Provides good variance-bias tradeoff",
      "- Consider repeating CV with different random seeds for robustness"
    )
    
  } else {
    # Large sample: use configurable k-fold (default 10-fold)
    method <- large_k
    description <- sprintf("%d-fold Cross-Validation", large_k)
    rationale <- sprintf(
      "%d-fold CV recommended: n = %d (total), min group = %d. Standard choice for well-powered studies.",
      large_k, n_total, n_min_group
    )
    considerations <- c(
      sprintf("- %d-fold is the standard for large samples", large_k),
      "- Fast and provides reliable estimates",
      "- Stratified by default to maintain group proportions"
    )
  }
  
  result <- list(
    method = method,
    description = description,
    rationale = rationale,
    considerations = considerations
  )
  
  if (verbose) {
    cat("\n========================================\n")
    cat("Cross-Validation Method Recommendation\n")
    cat("========================================\n\n")
    cat("Recommended:", description, "\n\n")
    cat("Rationale:", rationale, "\n\n")
    cat("Considerations:\n")
    for (consideration in considerations) {
      cat(" ", consideration, "\n")
    }
    cat("\n")
  }
  
  return(result)
}

#' Fit a PLS-DA model to preprocessed cytokine data
#'
#' @description
#' Fits a Partial Least Squares Discriminant Analysis (PLS-DA) or 
#' Orthogonal PLS-DA (OPLS-DA) model using the ropls package. This function
#' provides a convenient wrapper with structured output compatible with
#' immunoPlex plotting functions.
#'
#' @param preprocessed_data A plsda_preprocessed object from plsda_preprocess()
#' @param response_var character string naming the column in metadata to use
#'   as the response variable (e.g., "group", "disease")
#' @param n_components integer specifying number of predictive components.
#'   Default: 2
#' @param method character string specifying model type. Options: "PLS-DA" 
#'   (standard PLS-DA) or "OPLS-DA" (orthogonal PLS-DA). Default: "PLS-DA"
#' @param cross_validation integer specifying number of cross-validation folds,
#'   or character "loo"/"loocv" for Leave-One-Out Cross-Validation.
#'   Use LOOCV for small samples (n < 50) or when groups have < 15 samples.
#'   Use k-fold (5-10) for larger samples. If NULL, automatically selects based
#'   on sample size. Default: NULL (auto-select)
#' @param permutations integer specifying number of permutation tests for 
#'   assessing model significance. Set to 0 to skip. Default: 100
#' @param scale_method character string for scaling method passed to ropls.
#'   Options: "none", "center", "pareto", "standard". Since plsda_preprocess
#'   typically scales data, "none" is often appropriate. Default: "none"
#' @param verbose logical indicating whether to print ropls progress messages.
#'   Default: TRUE
#' @param loocv_n_threshold integer, threshold for total samples below which LOOCV is 
#'   auto-selected. Default: 50
#' @param loocv_group_threshold integer, threshold for minimum group size below which 
#'   LOOCV is auto-selected. Default: 15
#' @param moderate_n_threshold integer, threshold for total samples below which moderate 
#'   k-fold CV is auto-selected. Default: 100
#' @param moderate_k integer, number of folds for moderate sample sizes. Default: 7
#' @param large_k integer, number of folds for large sample sizes. Default: 10
#' @param rep_col character string naming a column in
#'   \code{preprocessed_data$metadata} that identifies technical replicates
#'   within a biological sample. When non-NULL, \code{plsda_fit} calls
#'   \code{\link{aggregate_replicates}} on the preprocessed expression
#'   matrix and collapses rows that share identical metadata-except-rep_col
#'   values to one row per biological sample, using the rule in
#'   \code{replicate_agg}. When NULL (default), the preprocessed data is
#'   fit as-is. PLS-DA cannot consume replicate-level rows directly (it
#'   has no random-effect machinery); pre-aggregate via this argument or
#'   call \code{aggregate_replicates} upstream. Default: NULL.
#'
#'   \strong{Cross-stage caveat.} Aggregation via \code{rep_col} runs on
#'   the \emph{already-preprocessed} expression matrix (after LOD
#'   substitution, log transformation, and scaling via
#'   \code{\link{plsda_preprocess}}). Because mean / geomean / median do
#'   not commute with log- and z-scoring, collapsing here is
#'   \strong{not} identical to aggregating replicates in the raw data
#'   before \code{plsda_preprocess}: the scale/center parameters are
#'   estimated on the replicate-level matrix and inherited by the
#'   aggregated rows. If biological-sample-level scaling is desired,
#'   call \code{\link{aggregate_replicates}} on the raw data and pass
#'   the aggregated frame into \code{plsda_preprocess} instead.
#' @param replicate_agg character, aggregation rule for
#'   \code{rep_col}. Passed through to \code{aggregate_replicates(rule =)}
#'   on the numeric path. One of \code{"mean"} (default), \code{"median"},
#'   \code{"geomean"}. Ignored when \code{rep_col = NULL}.
#' @param ... Additional arguments passed to ropls::opls()
#'
#' @return A list of class "plsda_model" containing:
#'   \item{model}{The original ropls S4 model object}
#'   \item{scores}{data.frame with sample scores and metadata}
#'   \item{loadings}{data.frame with cytokine loadings}
#'   \item{vip_scores}{data.frame with Variable Importance in Projection scores}
#'   \item{model_stats}{data.frame with model statistics (R2X, R2Y, Q2, etc.) 
#'     For OPLS-DA models, cumulative statistics are extracted from the summaryDF 
#'     slot. Includes permutation p-values with both technical names (perm_pval_R2Y, 
#'     perm_pval_Q2) and user-friendly aliases (permutation_pvalue, 
#'     permutation_pvalue_Q2)}
#'   \item{response_var}{Name of the response variable}
#'   \item{response_levels}{Levels of the response factor}
#'   \item{n_components}{Number of components in the model}
#'   \item{method}{Model type (PLS-DA or OPLS-DA)}
#'   \item{cv_method}{Cross-validation method used}
#'   \item{cv_folds}{Number of CV folds used}
#'   \item{preprocessing_info}{Preprocessing information from input data}
#'
#' @details
#' PLS-DA is a supervised multivariate method that maximizes separation between
#' groups while explaining variation in the predictor variables (cytokines).
#' OPLS-DA extends this by separating systematic variation into predictive and
#' orthogonal (uncorrelated with response) components.
#'
#' Model diagnostics:
#' \itemize{
#'   \item R2X: Fraction of X variance explained by the model
#'   \item R2Y: Fraction of Y variance explained by the model
#'   \item Q2: Cross-validated predictive ability (Q2 > 0.5 is good)
#'   \item Permutation p-value: Tests if model performs better than random
#' }
#'
#' The function handles OPLS-DA models robustly by extracting cumulative statistics
#' from the ropls summaryDF slot when modelDF contains NAs. Permutation p-values
#' are stored with multiple field names for convenience: `permutation_pvalue` 
#' (recommended for user scripts) and `perm_pval_R2Y` (technical name).
#'
#' The function requires the ropls package to be installed. For Bioconductor
#' installation: BiocManager::install("ropls")
#'
#' @examples
#' \dontrun{
#' # Preprocess data
#' preprocessed <- plsda_preprocess(
#'   data = data,
#'   cytokine_cols = cytokines,
#'   metadata_cols = metadata,
#'   lod_lookup = lod_values,
#'   lod_method = "half"
#' )
#' 
#' # Fit PLS-DA model with auto-selected CV (recommended)
#' model <- plsda_fit(
#'   preprocessed_data = preprocessed,
#'   response_var = "group",
#'   n_components = 2,
#'   method = "PLS-DA"
#'   # cross_validation = NULL by default (auto-selects based on sample size)
#' )
#' 
#' # For small samples: explicitly use LOOCV
#' model_loocv <- plsda_fit(
#'   preprocessed_data = preprocessed,
#'   response_var = "group",
#'   n_components = 2,
#'   cross_validation = "loocv"  # or "loo"
#' )
#' 
#' # For larger samples: use k-fold
#' model_kfold <- plsda_fit(
#'   preprocessed_data = preprocessed,
#'   response_var = "disease",
#'   n_components = 2,
#'   cross_validation = 10  # 10-fold CV
#' )
#' 
#' # Fit OPLS-DA model with custom settings
#' opls_model <- plsda_fit(
#'   preprocessed_data = preprocessed,
#'   response_var = "disease",
#'   n_components = 1,
#'   method = "OPLS-DA",
#'   cross_validation = "loocv",
#'   permutations = 1000
#' )
#' 
#' # Get CV recommendation before fitting
#' cv_rec <- recommend_cv_method(n_total = 47, n_min_group = 15, verbose = TRUE)
#' 
#' # Use custom CV thresholds (e.g., more conservative LOOCV use)
#' model_custom <- plsda_fit(
#'   preprocessed_data = preprocessed,
#'   response_var = "group",
#'   n_components = 2,
#'   loocv_n_threshold = 75,        # Use LOOCV up to 75 samples (instead of 50)
#'   loocv_group_threshold = 20,    # Require 20+ per group for k-fold (instead of 15)
#'   moderate_k = 5                 # Use 5-fold for moderate samples (instead of 7)
#' )
#' 
#' # View model statistics
#' print(model)
#' summary(model)
#' }
#'
#' @export
plsda_fit <- function(preprocessed_data,
                      response_var,
                      n_components = 2,
                      method = c("PLS-DA", "OPLS-DA"),
                      cross_validation = NULL,
                      permutations = 100,
                      scale_method = "none",
                      verbose = TRUE,
                      loocv_n_threshold = 50,
                      loocv_group_threshold = 15,
                      moderate_n_threshold = 100,
                      moderate_k = 7,
                      large_k = 10,
                      rep_col = NULL,
                      replicate_agg = c("mean", "median", "geomean"),
                      ...) {

  replicate_agg <- match.arg(replicate_agg)

  # ===== SECTION 1: INPUT VALIDATION =====
  if (verbose) cat("Validating inputs...\n")
  
  # Check if ropls is available
  if (!requireNamespace("ropls", quietly = TRUE)) {
    stop("Package 'ropls' is required but not installed.\n",
         "Install from Bioconductor with:\n",
         "  if (!requireNamespace('BiocManager', quietly = TRUE))\n",
         "    install.packages('BiocManager')\n",
         "  BiocManager::install('ropls')")
  }
  
  # Validate preprocessed_data
  if (!inherits(preprocessed_data, "plsda_preprocessed")) {
    stop("preprocessed_data must be a plsda_preprocessed object from plsda_preprocess()")
  }
  
  # Match method argument
  method <- match.arg(method)
  
  # Validate response variable exists
  if (!response_var %in% names(preprocessed_data$metadata)) {
    available <- paste(names(preprocessed_data$metadata), collapse = ", ")
    stop("response_var '", response_var, "' not found in metadata.\n",
         "Available: ", available)
  }
  
  # Validate numeric parameters
  stopifnot(
    is.numeric(n_components), n_components >= 1,
    is.numeric(permutations), permutations >= 0
  )
  
  # Validate cross_validation parameter
  if (!is.null(cross_validation)) {
    if (is.character(cross_validation)) {
      if (!tolower(cross_validation) %in% c("loo", "loocv")) {
        stop("cross_validation must be an integer >= 3, or 'loo'/'loocv' for Leave-One-Out")
      }
    } else if (is.numeric(cross_validation)) {
      if (cross_validation < 3 && cross_validation != nrow(preprocessed_data$expression)) {
        stop("cross_validation must be >= 3 or equal to n (for LOOCV)")
      }
    }
  }
  
  # ===== SECTION 2: DATA PREPARATION =====
  if (verbose) cat("Preparing data for modeling...\n")
  
  # Extract expression matrix and metadata
  expr_mat <- preprocessed_data$expression
  meta_df <- preprocessed_data$metadata
  sample_ids <- preprocessed_data$sample_ids

  # ===== SECTION 2a: PRE-AGGREGATE TECHNICAL REPLICATES =====
  # PLS-DA cannot consume replicate-level rows (no random-effect machinery).
  # When rep_col is supplied, collapse rows that share identical
  # metadata-except-rep_col to one row per biological sample, using
  # aggregate_replicates() on each cytokine column.
  if (!is.null(rep_col)) {
    if (!is.character(rep_col) || length(rep_col) != 1L) {
      stop("`rep_col` must be a single character string.", call. = FALSE)
    }
    if (!rep_col %in% names(meta_df)) {
      stop("rep_col '", rep_col, "' not found in preprocessed_data$metadata.",
           call. = FALSE)
    }
    if (identical(rep_col, response_var)) {
      stop("rep_col and response_var cannot be the same column ('",
           rep_col, "').", call. = FALSE)
    }
    grp_cols <- setdiff(names(meta_df), rep_col)
    if (length(grp_cols) == 0L) {
      stop("Cannot aggregate replicates: metadata has no columns besides ",
           "rep_col. Add at least one biological-sample identifier ",
           "(e.g., subject_id) to preprocessed_data$metadata.",
           call. = FALSE)
    }
    if (verbose) {
      cat("Pre-aggregating replicates by rep_col='", rep_col,
          "' using rule='", replicate_agg, "'\n", sep = "")
    }

    # Synthetic subject key = concatenation of all non-rep_col metadata.
    # Use a control-character separator so it cannot collide with user strings.
    sep <- ""
    grp_key <- do.call(paste,
                       c(lapply(meta_df[, grp_cols, drop = FALSE], as.character),
                         sep = sep))

    # Long-format frame: one row per (sample, cytokine). aggregate_replicates
    # collapses by (subject_key, cytokine_name), ignoring rep_col for grouping.
    n_samp <- nrow(expr_mat)
    cyto_names <- colnames(expr_mat)
    long_df <- data.frame(
      .grp_key = rep(grp_key, times = ncol(expr_mat)),
      .cyto    = rep(cyto_names, each = n_samp),
      .rep     = rep(meta_df[[rep_col]], times = ncol(expr_mat)),
      .value   = as.vector(expr_mat),
      stringsAsFactors = FALSE
    )

    agg_long <- aggregate_replicates(
      data          = long_df,
      rep_col       = ".rep",
      subject_col   = ".grp_key",
      cytokine_col  = ".cyto",
      timepoint_col = NULL,
      value_col     = ".value",
      rule          = replicate_agg
    )

    # Pivot back to wide. Use fixed cytokine order from cyto_names to keep
    # column order identical to the pre-aggregation matrix.
    unique_keys <- unique(grp_key)
    expr_new <- matrix(
      NA_real_,
      nrow = length(unique_keys),
      ncol = length(cyto_names),
      dimnames = list(NULL, cyto_names)
    )
    row_idx <- match(agg_long$.grp_key, unique_keys)
    col_idx <- match(agg_long$.cyto,    cyto_names)
    expr_new[cbind(row_idx, col_idx)] <- agg_long$.value

    # Collapse metadata to one row per unique_key (first row wins; replicates
    # are assumed to share biological metadata).
    #
    # Note on "first-row-wins" semantics: grp_key is built from ALL non-rep_col
    # metadata columns, so within a grp_key group every non-rep_col column is
    # constant by construction. A replicate that diverges from its peers on,
    # say, response_var becomes its own grp_key group rather than getting
    # silently averaged into the majority label. That's defensible but also
    # unobvious -- users wanting strict "replicates share biology" checking
    # must validate metadata upstream (e.g., by grouping on a subject_id
    # column alone and verifying response_var is constant per subject).
    first_idx <- match(unique_keys, grp_key)
    meta_new <- meta_df[first_idx, grp_cols, drop = FALSE]
    rownames(meta_new) <- NULL

    # Sample IDs: first-in-group from the original sample_ids vector, when
    # available. Falls back to grp_key if sample_ids was NULL.
    sample_ids_new <- if (!is.null(sample_ids)) {
      sample_ids[first_idx]
    } else {
      unique_keys
    }

    if (verbose) {
      cat("  Collapsed ", n_samp, " rows to ", length(unique_keys),
          " biological samples\n", sep = "")
    }

    expr_mat   <- expr_new
    meta_df    <- meta_new
    sample_ids <- sample_ids_new

    # Reflect the collapse in preprocessing_info so downstream consumers
    # (print.plsda_preprocessed, diagnostics, reports) see the post-
    # aggregation row count instead of the stale replicate-level n_samples
    # from plsda_preprocess(). n_cytokines is unchanged.
    if (!is.null(preprocessed_data$preprocessing_info)) {
      pi <- preprocessed_data$preprocessing_info
      pi$n_samples          <- nrow(expr_mat)
      pi$replicate_averaged <- TRUE
      pi$replicate_col      <- rep_col
      pi$replicate_agg      <- replicate_agg
      preprocessed_data$preprocessing_info <- pi
    }
  }

  # Extract and prepare response variable
  response <- meta_df[[response_var]]
  
  # Convert to factor if not already
  if (!is.factor(response)) {
    if (verbose) cat("  Converting response variable to factor\n")
    response <- as.factor(response)
  }
  
  # Get response levels
  response_levels <- levels(response)
  n_levels <- length(response_levels)
  
  if (verbose) {
    cat("  Response variable: ", response_var, "\n")
    cat("  Number of levels: ", n_levels, "\n")
    cat("  Levels: ", paste(response_levels, collapse = ", "), "\n")
    cat("  Samples per level:\n")
    print(table(response))
  }
  
  # Check minimum samples per group
  min_samples <- min(table(response))
  total_samples <- nrow(expr_mat)
  
  if (min_samples < 3) {
    warning("Some groups have fewer than 3 samples. Model may be unstable.")
  }
  
  # Auto-select cross-validation method if not specified
  if (is.null(cross_validation)) {
    cv_recommendation <- recommend_cv_method(
      n_total = total_samples, 
      n_min_group = min_samples, 
      verbose = verbose,
      loocv_n_threshold = loocv_n_threshold,
      loocv_group_threshold = loocv_group_threshold,
      moderate_n_threshold = moderate_n_threshold,
      moderate_k = moderate_k,
      large_k = large_k
    )
    cross_validation <- cv_recommendation$method
    if (verbose) {
      cat("  Auto-selected cross-validation: ", cv_recommendation$description, "\n")
      cat("  Rationale: ", cv_recommendation$rationale, "\n")
    }
  }
  
  # Convert LOOCV character to numeric
  cv_numeric <- cross_validation
  if (is.character(cross_validation) && tolower(cross_validation) %in% c("loo", "loocv")) {
    cv_numeric <- total_samples  # ropls uses n for LOOCV
    if (verbose) cat("  Using Leave-One-Out Cross-Validation (n = ", total_samples, ")\n")
  } else if (is.numeric(cross_validation)) {
    if (cross_validation == total_samples) {
      if (verbose) cat("  Using Leave-One-Out Cross-Validation (n = ", total_samples, ")\n")
    } else {
      if (verbose) cat("  Using ", cross_validation, "-fold Cross-Validation\n")
    }
  }
  
  # ===== SECTION 3: FIT PLS-DA MODEL =====
  if (verbose) cat("\nFitting ", method, " model...\n")
  
  # Set up ropls parameters
  if (method == "PLS-DA") {
    # Standard PLS-DA: all components are predictive
    predI <- n_components
    orthoI <- 0
  } else {
    # OPLS-DA: 1 predictive + orthogonal components
    predI <- 1
    orthoI <- max(0, n_components - 1)
  }
  
  # Fit the model
  if (verbose) {
    pls_model <- ropls::opls(
      x = expr_mat,
      y = response,
      predI = predI,
      orthoI = orthoI,
      crossvalI = cv_numeric,
      permI = permutations,
      scaleC = scale_method,
      ...
    )
  } else {
    # Suppress ropls output
    pls_model <- suppressMessages(
      ropls::opls(
        x = expr_mat,
        y = response,
        predI = predI,
        orthoI = orthoI,
        crossvalI = cv_numeric,
        permI = permutations,
        scaleC = scale_method,
        ...
      )
    )
  }
  
  # ===== SECTION 4: EXTRACT MODEL COMPONENTS =====
  if (verbose) cat("Extracting model components...\n")
  
  # Extract scores (sample positions in PLS space)
  score_mat <- as.matrix(pls_model@scoreMN)
  
  # Add orthogonal scores if OPLS-DA
  if (method == "OPLS-DA" && orthoI > 0) {
    ortho_scores <- as.matrix(pls_model@orthoScoreMN)
    score_mat <- cbind(score_mat, ortho_scores)
    
    # Rename columns to distinguish predictive and orthogonal
    colnames(score_mat) <- c(
      paste0("p", 1:predI),
      paste0("o", 1:orthoI)
    )
  } else {
    colnames(score_mat) <- paste0("p", 1:ncol(score_mat))
  }
  
  # Combine scores with metadata
  scores_df <- as.data.frame(score_mat)
  scores_df <- cbind(scores_df, meta_df)
  rownames(scores_df) <- sample_ids
  
  # Extract loadings (cytokine contributions)
  loading_mat <- as.matrix(pls_model@loadingMN)
  loadings_df <- as.data.frame(loading_mat)
  loadings_df$cytokine <- rownames(loading_mat)
  
  # Add orthogonal loadings if OPLS-DA
  if (method == "OPLS-DA" && orthoI > 0) {
    ortho_loadings <- as.matrix(pls_model@orthoLoadingMN)
    ortho_df <- as.data.frame(ortho_loadings)
    colnames(ortho_df) <- paste0("o", 1:ncol(ortho_loadings))
    loadings_df <- cbind(loadings_df, ortho_df)
  }
  
  # Reorder columns for clarity
  loadings_df <- loadings_df[, c("cytokine", setdiff(names(loadings_df), "cytokine"))]
  
  # Extract VIP scores (Variable Importance in Projection)
  vip_values <- ropls::getVipVn(pls_model)
  vip_df <- data.frame(
    cytokine = names(vip_values),
    vip_score = as.numeric(vip_values),
    importance = ifelse(vip_values >= 1, "Important", "Less Important"),
    stringsAsFactors = FALSE
  )
  vip_df <- vip_df[order(vip_df$vip_score, decreasing = TRUE), ]
  rownames(vip_df) <- NULL
  
  # Extract model statistics
  # For OPLS-DA, @modelDF often has NAs for cumulative stats, so extract from @summaryDF
  model_summary <- pls_model@modelDF
  stats_df <- as.data.frame(model_summary)
  
  # For OPLS-DA or when @modelDF has NAs, get cumulative stats from @summaryDF
  if (method == "OPLS-DA" || any(is.na(stats_df$R2X)) || any(is.na(stats_df$R2Y)) || any(is.na(stats_df$Q2))) {
    summary_df <- pls_model@summaryDF
    
    # Find the "Total" row or last row with cumulative statistics
    if ("Total" %in% rownames(summary_df)) {
      total_row <- summary_df["Total", ]
    } else {
      # Use the last row which typically contains cumulative stats
      total_row <- summary_df[nrow(summary_df), ]
    }
    
    # Replace NA values in stats_df with cumulative values from summaryDF
    # Create a row for cumulative statistics
    cum_stats <- data.frame(
      R2X = as.numeric(total_row["R2X(cum)"]),
      R2Y = as.numeric(total_row["R2Y(cum)"]),
      Q2 = as.numeric(total_row["Q2(cum)"]),
      RMSEE = as.numeric(total_row["RMSEE"]),
      stringsAsFactors = FALSE
    )
    
    # If stats_df has all NAs, replace with cumulative stats
    if (all(is.na(stats_df$R2X)) || all(is.na(stats_df$R2Y))) {
      stats_df <- cum_stats
      rownames(stats_df) <- "Total"
    }
  }
  
  # Add permutation results if available
  if (permutations > 0) {
    perm_pval_R2Y <- pls_model@summaryDF$`pR2Y`[1]
    perm_pval_Q2 <- pls_model@summaryDF$`pQ2`[1]
    
    # Add to stats_df with both technical and user-friendly names
    stats_df$perm_pval_R2Y <- perm_pval_R2Y
    stats_df$perm_pval_Q2 <- perm_pval_Q2
    
    # Add user-friendly aliases for easier access
    stats_df$permutation_pvalue <- perm_pval_R2Y  # Most commonly used
    stats_df$permutation_pvalue_Q2 <- perm_pval_Q2
  }

  # Warn when the model has no cross-validated predictive ability.
  #
  # Q2(cum) < 0 means the model predicts held-out samples worse than simply
  # using the response mean, so it does not generalise. The usual cause is more
  # components than the available signal supports. This is worth an explicit
  # warning because VIP scores can still look convincing when it happens: a user
  # ranking analytes by VIP would otherwise get no indication that the model
  # behind that ranking fails cross-validation.
  #
  # Read Q2 from the model's own summary rather than stats_df, which is either
  # the per-component modelDF or a substituted cumulative frame depending on the
  # method and on whether modelDF held NAs.
  q2_cum <- tryCatch(
    suppressWarnings(as.numeric(pls_model@summaryDF[["Q2(cum)"]][1])),
    error = function(e) NA_real_
  )
  if (isTRUE(is.finite(q2_cum)) && q2_cum < 0) {
    warning(sprintf(
      paste0("Model has no cross-validated predictive ability (Q2 = %.3f). ",
             "It predicts held-out samples worse than the response mean, so it ",
             "does not generalise.%s VIP scores from this model are not ",
             "validated by cross-validation and should be interpreted with ",
             "caution."),
      q2_cum,
      if (is.numeric(n_components) && !is.na(n_components) && n_components > 1) {
        sprintf(" This often means too many components for the signal present; %d were requested, so try fewer.",
                as.integer(n_components))
      } else {
        ""
      }
    ), call. = FALSE)
  }

  # ===== SECTION 5: PREPARE OUTPUT =====
  if (verbose) cat("Preparing output...\n")
  
  result <- list(
    model = pls_model,
    scores = scores_df,
    loadings = loadings_df,
    vip_scores = vip_df,
    model_stats = stats_df,
    response_var = response_var,
    response_levels = response_levels,
    n_components = n_components,
    method = method,
    cv_method = cross_validation,
    cv_folds = cv_numeric,
    preprocessing_info = preprocessed_data$preprocessing_info
  )
  
  # Add S3 class
  class(result) <- c("plsda_model", "list")
  
  # Add attributes
  attr(result, "method") <- method
  attr(result, "n_components") <- n_components
  attr(result, "response_var") <- response_var
  
  if (verbose) {
    cat("\nModel fitting complete!\n")
    cat("  Method: ", method, "\n")
    cat("  Components: ", n_components, "\n")
    
    # Safely display R2X, R2Y, Q2 (handling NAs)
    r2x_val <- sum(stats_df$R2X, na.rm = TRUE)
    r2y_val <- sum(stats_df$R2Y, na.rm = TRUE)
    q2_val <- sum(stats_df$Q2, na.rm = TRUE)
    
    cat("  R2X: ", if (!is.na(r2x_val)) round(r2x_val, 3) else "NA", "\n")
    cat("  R2Y: ", if (!is.na(r2y_val)) round(r2y_val, 3) else "NA", "\n")
    cat("  Q2: ", if (!is.na(q2_val)) round(q2_val, 3) else "NA", "\n")
    
    if (permutations > 0 && exists("perm_pval_R2Y")) {
      cat("  Permutation p-value (R2Y): ", 
          format.pval(perm_pval_R2Y, digits = 3), "\n")
    }
  }
  
  return(result)
}


#' Print method for plsda_model objects
#'
#' @param x A plsda_model object
#' @param ... Additional arguments (not used)
#' @export
print.plsda_model <- function(x, ...) {
  cat("PLS-DA Model\n")
  cat("============\n\n")
  cat("Method:", x$method, "\n")
  cat("Components:", x$n_components, "\n")
  cat("Response variable:", x$response_var, "\n")
  cat("Response levels:", paste(x$response_levels, collapse = ", "), "\n")
  cat("Samples:", nrow(x$scores), "\n")
  cat("Cytokines:", nrow(x$loadings), "\n")
  
  # Display CV method
  if (!is.null(x$cv_method)) {
    if (is.character(x$cv_method) || x$cv_folds == nrow(x$scores)) {
      cat("Cross-validation: Leave-One-Out (LOOCV)\n")
    } else {
      cat("Cross-validation:", x$cv_folds, "-fold\n")
    }
  }
  cat("\n")
  
  cat("Model Performance:\n")
  stats <- x$model_stats
  
  # Safely handle NA values in model statistics
  r2x_sum <- sum(stats$R2X, na.rm = TRUE)
  r2y_sum <- sum(stats$R2Y, na.rm = TRUE)
  q2_sum <- sum(stats$Q2, na.rm = TRUE)
  
  cat("  R2X (cumulative):", if (!is.na(r2x_sum)) round(r2x_sum, 3) else "NA", "\n")
  cat("  R2Y (cumulative):", if (!is.na(r2y_sum)) round(r2y_sum, 3) else "NA", "\n")
  cat("  Q2 (cumulative):", if (!is.na(q2_sum)) round(q2_sum, 3) else "NA", "\n")
  
  # Check for permutation p-value under multiple possible names
  if ("permutation_pvalue" %in% names(stats) && !is.na(stats$permutation_pvalue[1])) {
    cat("  Permutation p-value:", format.pval(stats$permutation_pvalue[1], digits = 3), "\n")
  } else if ("perm_pval_R2Y" %in% names(stats) && !is.na(stats$perm_pval_R2Y[1])) {
    cat("  Permutation p-value:", format.pval(stats$perm_pval_R2Y[1], digits = 3), "\n")
  }
  
  cat("\nImportant cytokines (VIP >= 1):", 
      sum(x$vip_scores$vip_score >= 1), "out of", nrow(x$vip_scores), "\n")
  
  invisible(x)
}


#' Summary method for plsda_model objects
#'
#' @param object A plsda_model object
#' @param ... Additional arguments (not used)
#' @export
summary.plsda_model <- function(object, ...) {
  cat("PLS-DA Model Summary\n")
  cat("====================\n\n")
  
  # Print basic info
  print(object)
  
  # Detailed statistics per component
  cat("\n\nPer-Component Statistics:\n")
  print(object$model_stats, row.names = TRUE)
  
  # Top VIP scores
  cat("\n\nTop 10 Cytokines by VIP Score:\n")
  top_vip <- head(object$vip_scores, 10)
  print(top_vip, row.names = FALSE)
  
  # Group statistics
  cat("\n\nSample Distribution:\n")
  response_col <- object$response_var
  if (response_col %in% names(object$scores)) {
    print(table(object$scores[[response_col]]))
  }
  
  invisible(object)
}


#' Extract model predictions from plsda_model
#'
#' @param object A plsda_model object
#' @param newdata Optional new data to predict (must be preprocessed same way)
#' @param ... Additional arguments (not used)
#' @return A vector of predicted class labels
#' @export
predict.plsda_model <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    # Return fitted predictions
    predictions <- predict(object$model, object$model@suppLs$xModelMN)
  } else {
    # Predict on new data
    predictions <- predict(object$model, newdata)
  }
  return(predictions)
}
