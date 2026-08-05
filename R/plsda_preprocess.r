# PLS-DA Preprocessing Module for immunoPlex Package
# Preprocessing functionality specifically for PLS-DA analysis with flexible LOD handling

#' Adjust expression data for covariates (internal helper)
#'
#' @description
#' Internal function to regress out covariate effects from cytokine expression data.
#' Performs within-subject imputation for missing covariate values when subject ID
#' is provided, then fits linear models for each cytokine and returns residuals.
#'
#' @param expr_mat Matrix of cytokine expression (samples x cytokines)
#' @param metadata Data frame containing covariate columns
#' @param covariates Character vector of covariate column names
#' @param subject_id_col Character string naming the subject ID column for 
#'   within-subject imputation of missing covariates. If NULL, no imputation. Default: NULL
#' @param verbose Logical for progress messages. Default: TRUE
#'
#' @return List containing:
#'   \item{adjusted_matrix}{Matrix of residuals (covariate-adjusted expression)}
#'   \item{samples_kept}{Indices of samples retained after removing missing covariates}
#'   \item{n_removed}{Number of samples removed due to missing covariates}
#'   \item{imputed_info}{Summary of imputation performed}
#'
#' @noRd
adjust_for_covariates_internal <- function(expr_mat, metadata, 
                                          covariates, 
                                          subject_id_col = NULL,
                                          verbose = TRUE) {
  
  # Validate covariates exist in metadata
  missing_covs <- setdiff(covariates, names(metadata))
  if (length(missing_covs) > 0) {
    stop("Covariates not found in metadata: ", paste(missing_covs, collapse = ", "))
  }
  
  # Extract covariate data
  covar_data <- metadata[, covariates, drop = FALSE]
  n_original <- nrow(covar_data)
  n_imputed <- 0
  
  # Check for single-level covariates (nothing to adjust)
  single_level_covariates <- c()
  for (cov in covariates) {
    n_levels <- length(unique(covar_data[[cov]][!is.na(covar_data[[cov]])]))
    if (n_levels <= 1) {
      single_level_covariates <- c(single_level_covariates, cov)
      if (verbose) {
        cat("  Covariate '", cov, "' has only ", n_levels, 
            " level - skipping adjustment (nothing to adjust)\n", sep = "")
      }
    }
  }
  
  # Remove single-level covariates from adjustment list
  if (length(single_level_covariates) > 0) {
    covariates <- setdiff(covariates, single_level_covariates)
    if (length(covariates) == 0) {
      if (verbose) {
        cat("  All covariates have single levels - no adjustment needed\n")
      }
      # Return original data unchanged
      return(list(
        adjusted_matrix = expr_mat,
        samples_kept = rep(TRUE, nrow(expr_mat)),
        n_removed = 0,
        imputed_info = list(
          n_imputed = 0,
          imputation_method = "none"
        ),
        single_level_covariates = single_level_covariates
      ))
    }
  }
  
  # Within-subject imputation for missing covariate values
  if (!is.null(subject_id_col)) {
    if (!subject_id_col %in% names(metadata)) {
      warning("subject_id_col '", subject_id_col, "' not found in metadata. Skipping imputation.")
    } else {
      subject_ids <- metadata[[subject_id_col]]
      
      for (cov in covariates) {
        # Find missing values
        missing_idx <- is.na(covar_data[[cov]])
        
        if (any(missing_idx)) {
          # For each subject with missing values, impute from other samples of same subject
          missing_subjects <- unique(subject_ids[missing_idx])
          
          for (subj in missing_subjects) {
            subj_idx <- subject_ids == subj
            subj_values <- covar_data[[cov]][subj_idx]
            
            # If subject has any non-missing values, use the first non-NA value
            non_na_values <- subj_values[!is.na(subj_values)]
            if (length(non_na_values) > 0) {
              # Impute with first non-NA value for this subject
              impute_value <- non_na_values[1]
              covar_data[[cov]][subj_idx & is.na(covar_data[[cov]])] <- impute_value
              n_imputed <- n_imputed + sum(subj_idx & missing_idx)
            }
          }
        }
      }
      
      if (verbose && n_imputed > 0) {
        cat("  Imputed", n_imputed, "missing covariate values within subjects\n")
      }
    }
  }
  
  # Check for remaining missing values
  complete_cases <- complete.cases(covar_data)
  n_removed <- sum(!complete_cases)
  
  if (n_removed > 0) {
    if (verbose) {
      cat("  Removing", n_removed, "samples with missing covariate values\n")
    }
    
    # Keep only complete cases
    expr_mat <- expr_mat[complete_cases, , drop = FALSE]
    covar_data <- covar_data[complete_cases, , drop = FALSE]
  }
  
  # Adjust each cytokine for covariates
  if (verbose) {
    cat("  Fitting linear models for", ncol(expr_mat), "cytokines\n")
  }
  
  adjusted_matrix <- apply(expr_mat, 2, function(cytokine_values) {
    # Build formula
    covar_names <- paste(covariates, collapse = " + ")
    model_formula <- as.formula(paste("cytokine_values ~", covar_names))
    
    # Fit linear model
    lm_fit <- lm(model_formula, data = covar_data)
    
    # Return residuals (covariate-adjusted values)
    return(residuals(lm_fit))
  })
  
  # Preserve row and column names
  rownames(adjusted_matrix) <- rownames(expr_mat)
  colnames(adjusted_matrix) <- colnames(expr_mat)
  
  # Return results
  list(
    adjusted_matrix = adjusted_matrix,
    samples_kept = complete_cases,
    n_removed = n_removed,
    imputed_info = list(
      n_imputed = n_imputed,
      imputation_method = ifelse(is.null(subject_id_col), "none", "within_subject")
    ),
    single_level_covariates = if (exists("single_level_covariates")) single_level_covariates else c(),
    covariates_adjusted = covariates
  )
}


#' Preprocess cytokine data for PLS-DA analysis
#'
#' @description
#' Prepares merged cytokine expression data for PLS-DA analysis by handling 
#' values below the limit of detection (LOD), applying log transformation,
#' and standardizing (z-score scaling) the data. This function accepts data
#' in merged format where expression and metadata are combined in a single
#' data frame.
#'
#' @param data data.frame containing both cytokine measurements and metadata
#'   in merged format (e.g., as in data.csv)
#' @param cytokine_cols character vector specifying column names containing
#'   cytokine measurements
#' @param metadata_cols character vector specifying column names containing
#'   sample metadata (e.g., "SubjectID", "group", "typeofsample")
#' @param lod_lookup data.frame with columns 'cytokine' and 'lod' containing
#'   limit of detection values for each cytokine
#' @param lod_method character string specifying LOD substitution method.
#'   Options: "half" (LOD/2, default), "zero" (0), "sqrt" (sqrt(LOD)),
#'   "lod" (LOD value), "uniform" (random uniform \[0, LOD\]), 
#'   "halfmin" (half of minimum positive value, good for log transforms)
#' @param log_transform logical indicating whether to apply log2(x + 1)
#'   transformation. Default: TRUE
#' @param transform character, optional overall preprocessing transform. When
#'   non-NULL this overrides \code{log_transform} and \code{scale_data}.
#'   Valid values: \code{"log2p1_zscore"} (log2(x+1) then z-score),
#'   \code{"cuberoot_zscore"} (sign(x)*|x|^(1/3) then z-score),
#'   \code{"pareto"} (log2(x+1) then Pareto scaling, i.e. mean-center then
#'   divide by sqrt(sd)). Default: NULL.
#' @param scale_data logical indicating whether to apply z-score scaling
#'   (standardization) across cytokines. Recommended for PLS-DA. Default: TRUE
#' @param replicate_col character string naming the column to use for grouping
#'   technical replicates (e.g., "SubjectID"). If NULL, no averaging. Default: NULL
#' @param replicate_fn function to use for aggregating replicates. Default: mean
#' @param covariates character vector specifying covariate column names to regress out
#'   from cytokine expression (e.g., c("age", "gest_collect")). Adjustment occurs after
#'   log transformation but before z-score scaling. If NULL, no adjustment. Default: NULL
#' @param subject_id_col character string naming the subject ID column for within-subject
#'   imputation of missing covariate values. If provided, missing covariate values will be
#'   imputed using non-missing values from the same subject. If NULL, no imputation. Default: NULL
#' @param response_var character string naming the response variable column for intelligent
#'   covariate adjustment. If provided along with covariates and check_balance = TRUE, the
#'   function will test group balance and only adjust for imbalanced covariates. Default: NULL
#' @param check_balance logical indicating whether to check group balance before adjusting.
#'   Only used when both covariates and response_var are provided. If TRUE, performs t-tests
#'   to identify imbalanced covariates (p < balance_threshold) and adjusts only those.
#'   If FALSE, adjusts all specified covariates. Default: TRUE
#' @param balance_threshold numeric p-value threshold for determining covariate imbalance.
#'   Covariates with t-test p < threshold are considered imbalanced and adjusted. Default: 0.05
#' @param verbose logical indicating whether to print progress messages. Default: TRUE
#'
#' @return A list of class "plsda_preprocessed" containing:
#'   \item{expression}{matrix of preprocessed cytokine expression (samples x cytokines)}
#'   \item{metadata}{data.frame of sample metadata}
#'   \item{preprocessing_info}{list of preprocessing parameters used, including covariate adjustment details}
#'   \item{sample_ids}{character vector of sample identifiers}
#'
#' @details
#' The preprocessing pipeline consists of the following steps:
#' \enumerate{
#'   \item Extract cytokine measurements and metadata from merged data
#'   \item Apply LOD substitution to handle censored values
#'   \item Optionally average technical replicates
#'   \item Apply log2(x + 1) transformation if requested
#'   \item Optionally adjust for covariates (regress out covariate effects)
#'   \item Apply z-score scaling (mean=0, sd=1) if requested
#' }
#'
#' Covariate adjustment (if covariates parameter is provided):
#' \itemize{
#'   \item Performed AFTER log transformation but BEFORE z-score scaling
#'   \item For each cytokine, fits linear model: cytokine ~ covariate1 + covariate2 + ...
#'   \item Residuals (adjusted values) replace original expression values
#'   \item Missing covariate values can be imputed within subjects if subject_id_col is provided
#'   \item Samples with remaining missing covariates are removed with warning
#'   \item Makes PLS-DA results comparable to regression models that control for covariates
#'   \item Adjustment information is stored in preprocessing_info for reproducibility
#' }
#'
#' Intelligent covariate adjustment (if response_var provided):
#' \itemize{
#'   \item When response_var is provided, the function can check group balance
#'   \item For each covariate, performs t-test comparing groups
#'   \item Only adjusts for covariates where groups significantly differ (p < balance_threshold)
#'   \item This is data-driven: adjusts only when statistically justified
#'   \item Set check_balance = FALSE to force adjustment of all specified covariates
#'   \item Balance check results are stored in preprocessing_info
#'   \item Particularly useful for case-control studies where groups may differ on demographics
#' }
#'
#' LOD substitution methods:
#' \itemize{
#'   \item "half": Replace with LOD/2 (most common, assumes log-normal)
#'   \item "zero": Replace with 0 (conservative)
#'   \item "sqrt": Replace with sqrt(LOD) (intermediate)
#'   \item "lod": Replace with LOD value (sets to detection limit)
#'   \item "uniform": Replace with random uniform(0, LOD) (adds variability)
#'   \item "halfmin": Replace with half of minimum positive value (ensures no zeros for log)
#' }
#'
#' Fallback for cytokines without LOD values:
#' \itemize{
#'   \item If a cytokine has no LOD value in lod_lookup, the minimum positive 
#'         value for that cytokine in the dataset will be used as an empirical LOD
#'   \item This allows consistent handling of missing/censored values across all cytokines
#'   \item Values below this empirical LOD (including NAs set to 0) will be substituted 
#'         using the selected lod_method
#' }
#'
#' Z-score scaling is recommended for PLS-DA because:
#' \itemize{
#'   \item PLS-DA is sensitive to variable scales
#'   \item Cytokines have vastly different concentration ranges
#'   \item Without standardization, high-variance cytokines dominate the model
#'   \item Standard practice in multivariate metabolomics/proteomics
#' }
#'
#' @examples
#' \dontrun{
#' # Load merged data
#' data <- read.csv("data.csv")
#' lod_values <- read.csv("cytokine_LOD_values.csv")
#' 
#' # Define cytokine columns
#' cytokines <- c("EGF", "Eotaxin", "GCSF", "GMCSF", "GROa")
#' metadata <- c("SubjectID", "group", "typeofsample", "maternal_age")
#' 
#' # Preprocess for PLS-DA
#' preprocessed <- plsda_preprocess(
#'   data = data,
#'   cytokine_cols = cytokines,
#'   metadata_cols = metadata,
#'   lod_lookup = lod_values,
#'   lod_method = "half",
#'   log_transform = TRUE,
#'   scale_data = TRUE
#' )
#' 
#' # With replicate averaging
#' preprocessed <- plsda_preprocess(
#'   data = data,
#'   cytokine_cols = cytokines,
#'   metadata_cols = metadata,
#'   lod_lookup = lod_values,
#'   lod_method = "halfmin",
#'   replicate_col = "SubjectID"
#' )
#' 
#' # With intelligent covariate adjustment (auto-checks balance)
#' preprocessed <- plsda_preprocess(
#'   data = data,
#'   cytokine_cols = cytokines,
#'   metadata_cols = c(metadata, "age", "gest_collect", "group"),
#'   lod_lookup = lod_values,
#'   covariates = c("age", "gest_collect"),
#'   response_var = "group",        # Enable balance checking
#'   check_balance = TRUE,           # Only adjust if imbalanced
#'   subject_id_col = "SubjectID"
#' )
#' 
#' # Force covariate adjustment (skip balance check)
#' preprocessed <- plsda_preprocess(
#'   data = data,
#'   cytokine_cols = cytokines,
#'   metadata_cols = c(metadata, "age"),
#'   lod_lookup = lod_values,
#'   covariates = c("age"),
#'   check_balance = FALSE           # Always adjust
#' )
#' }
#'
#' @importFrom stats runif sd
#' @export
plsda_preprocess <- function(data,
                              cytokine_cols,
                              metadata_cols,
                              lod_lookup,
                              lod_method = "half",
                              log_transform = TRUE,
                              scale_data = TRUE,
                              transform = NULL,
                              replicate_col = NULL,
                              replicate_fn = mean,
                              covariates = NULL,
                              subject_id_col = NULL,
                              response_var = NULL,
                              check_balance = TRUE,
                              balance_threshold = 0.05,
                              verbose = TRUE) {

  # Validate `transform` (overrides log_transform + scale_data when non-NULL)
  valid_transforms <- c("log2p1_zscore", "cuberoot_zscore", "pareto")
  if (!is.null(transform)) {
    if (!is.character(transform) || length(transform) != 1L ||
        !(transform %in% valid_transforms)) {
      stop("transform must be one of: ",
           paste(valid_transforms, collapse = ", "))
    }
    # Disable the separate log_transform / scale_data paths; SECTION 5/6 will
    # apply the chosen transform end-to-end.
    log_transform <- FALSE
    scale_data    <- FALSE
  }
  
  # ===== SECTION 1: INPUT VALIDATION =====
  if (verbose) cat("Validating inputs...\n")
  
  # Validate basic structure
  stopifnot(
    is.data.frame(data),
    is.character(cytokine_cols),
    is.character(metadata_cols),
    is.data.frame(lod_lookup),
    is.logical(log_transform),
    is.logical(scale_data)
  )
  
  # Validate LOD lookup structure
  if (!all(c("cytokine", "lod") %in% names(lod_lookup))) {
    stop("lod_lookup must contain columns 'cytokine' and 'lod'")
  }
  
  # Validate column existence
  missing_cytokines <- setdiff(cytokine_cols, names(data))
  if (length(missing_cytokines) > 0) {
    stop("Cytokine columns not found in data: ", paste(missing_cytokines, collapse = ", "))
  }
  
  missing_metadata <- setdiff(metadata_cols, names(data))
  if (length(missing_metadata) > 0) {
    stop("Metadata columns not found in data: ", paste(missing_metadata, collapse = ", "))
  }
  
  # Validate LOD method
  valid_methods <- c("half", "zero", "sqrt", "lod", "uniform", "halfmin")
  if (!lod_method %in% valid_methods) {
    stop("lod_method must be one of: ", paste(valid_methods, collapse = ", "))
  }
  
  # Check for cytokines without LOD values (warn but allow processing)
  missing_lods <- setdiff(cytokine_cols, lod_lookup$cytokine)
  if (length(missing_lods) > 0) {
    if (verbose) {
      cat("Note: ", length(missing_lods), " cytokine(s) without LOD values will be processed without substitution:\n")
      cat("  ", paste(missing_lods, collapse = ", "), "\n")
    }
  }
  
  # Validate replicate column if specified
  if (!is.null(replicate_col)) {
    if (!replicate_col %in% names(data)) {
      stop("replicate_col '", replicate_col, "' not found in data")
    }
  }
  
  # ===== SECTION 2: DATA EXTRACTION =====
  if (verbose) cat("Extracting cytokine measurements and metadata...\n")
  
  # Extract cytokine expression data
  expr_raw <- data[, cytokine_cols, drop = FALSE]
  
  # Create LOD lookup vector for quick access (convert to numeric)
  lod_vec <- setNames(suppressWarnings(as.numeric(lod_lookup$lod)), lod_lookup$cytokine)
  
  # Clean cytokine data (handle string values like "< 12.8", "> 6000")
  # The key insight: these string values represent censored data
  # We parse the numeric value but mark as censored based on LOD comparison
  for (cyto in cytokine_cols) {
    col_data <- expr_raw[[cyto]]
    cyto_lod <- lod_vec[cyto]
    
    # If character, clean it
    if (is.character(col_data)) {
      # Store original for detecting symbols
      original_values <- col_data
      
      # Remove '<', '>', and spaces, keep only numbers
      col_data <- gsub("[<>\\s]", "", col_data)
      # Convert to numeric
      col_data <- suppressWarnings(as.numeric(col_data))
      
      # For values that were originally "< X" or empty/NA, 
      # we should treat them as being below LOD
      # Set them to 0 temporarily - they'll be handled by LOD substitution
      needs_substitution <- grepl("^<|^\\s*$", original_values) | is.na(col_data)
      col_data[needs_substitution] <- 0
      
      expr_raw[[cyto]] <- col_data
    } else if (!is.numeric(col_data)) {
      # Try to convert to numeric
      expr_raw[[cyto]] <- suppressWarnings(as.numeric(as.character(col_data)))
      # Replace NAs with 0 (will be handled by LOD substitution)
      expr_raw[[cyto]][is.na(expr_raw[[cyto]])] <- 0
    } else {
      # Already numeric, but replace NAs with 0 (will be handled by LOD substitution)
      expr_raw[[cyto]][is.na(col_data)] <- 0
    }
  }
  
  # Convert to matrix
  expr_raw <- as.matrix(expr_raw)
  
  # Extract metadata (convert to data.frame to avoid tibble rownames deprecation)
  meta_raw <- as.data.frame(data[, metadata_cols, drop = FALSE])
  
  # Create sample IDs (use row names if available, otherwise create)
  if (!is.null(rownames(data)) && !all(rownames(data) == as.character(1:nrow(data)))) {
    sample_ids <- rownames(data)
  } else {
    sample_ids <- paste0("Sample_", 1:nrow(data))
  }
  rownames(expr_raw) <- sample_ids
  rownames(meta_raw) <- sample_ids
  
  # ===== SECTION 3: LOD SUBSTITUTION =====
  if (verbose) cat("Applying LOD substitution method: ", lod_method, "\n")
  
  expr_lod <- expr_raw
  
  # Apply LOD substitution based on method
  # Note: lod_vec was already created in Section 2
  for (cyto in cytokine_cols) {
    # Get LOD value, handling missing cytokines in lookup
    if (cyto %in% names(lod_vec)) {
      cyto_lod <- unname(lod_vec[cyto])  # Remove name attribute
    } else {
      cyto_lod <- NA_real_
    }
    
    # Fallback: If no LOD value, use minimum positive value as empirical LOD
    if (is.na(cyto_lod)) {
      min_pos <- min(expr_lod[expr_lod[, cyto] > 0, cyto], na.rm = TRUE)
      if (is.finite(min_pos) && min_pos > 0) {
        cyto_lod <- as.numeric(min_pos)  # Ensure numeric
        if (verbose) cat("  Cytokine", cyto, "has no LOD value - using minimum positive value as LOD:", 
                         round(cyto_lod, 2), "\n")
      } else {
        if (verbose) cat("  Cytokine", cyto, "has no LOD value and no positive values - skipping substitution\n")
        next
      }
    }
    
    # Identify censored values (below LOD)
    censored_idx <- expr_lod[, cyto] < cyto_lod
    n_censored <- sum(censored_idx, na.rm = TRUE)
    
    if (n_censored > 0) {
      if (lod_method == "half") {
        expr_lod[censored_idx, cyto] <- cyto_lod / 2
      } else if (lod_method == "zero") {
        expr_lod[censored_idx, cyto] <- 0
      } else if (lod_method == "sqrt") {
        expr_lod[censored_idx, cyto] <- sqrt(cyto_lod)
      } else if (lod_method == "lod") {
        expr_lod[censored_idx, cyto] <- cyto_lod
      } else if (lod_method == "uniform") {
        expr_lod[censored_idx, cyto] <- runif(n_censored, 0, cyto_lod)
      } else if (lod_method == "halfmin") {
        # Use half of minimum positive value across ALL data
        min_pos <- min(expr_lod[expr_lod > 0], na.rm = TRUE)
        if (is.finite(min_pos)) {
          expr_lod[censored_idx, cyto] <- min_pos / 2
        } else {
          # Fallback to LOD/2 if no positive values
          expr_lod[censored_idx, cyto] <- cyto_lod / 2
        }
      }
    }
  }
  
  # ===== SECTION 4: REPLICATE AVERAGING =====
  if (!is.null(replicate_col)) {
    if (verbose) cat("Averaging technical replicates by: ", replicate_col, "\n")
    
    # Get unique replicate groups
    replicate_groups <- meta_raw[[replicate_col]]
    
    # Calculate means for each replicate group
    unique_groups <- unique(replicate_groups)
    n_orig <- nrow(expr_lod)
    n_final <- length(unique_groups)
    
    expr_avg <- matrix(NA, nrow = n_final, ncol = ncol(expr_lod))
    colnames(expr_avg) <- colnames(expr_lod)
    rownames(expr_avg) <- unique_groups
    
    # Create as data.frame directly with rownames to avoid tibble issues
    meta_avg <- as.data.frame(matrix(NA, nrow = n_final, ncol = ncol(meta_raw)))
    colnames(meta_avg) <- colnames(meta_raw)
    rownames(meta_avg) <- unique_groups
    
    for (grp in unique_groups) {
      grp_idx <- which(replicate_groups == grp)
      
      # Average expression values
      if (length(grp_idx) == 1) {
        expr_avg[as.character(grp), ] <- expr_lod[grp_idx, ]
      } else {
        expr_avg[as.character(grp), ] <- apply(expr_lod[grp_idx, , drop = FALSE], 2, replicate_fn, na.rm = TRUE)
      }
      
      # Take first replicate's metadata (assuming consistent within group)
      meta_avg[as.character(grp), ] <- meta_raw[grp_idx[1], ]
    }
    
    expr_lod <- expr_avg
    meta_raw <- meta_avg
    sample_ids <- unique_groups
    
    if (verbose) cat("  Reduced from ", n_orig, " to ", n_final, " samples\n")
  }
  
  # ===== SECTION 5: TRANSFORMATION =====
  if (!is.null(transform)) {
    if (verbose) cat("Applying transform: ", transform, "\n", sep = "")
    if (transform == "log2p1_zscore") {
      if (any(expr_lod < 0, na.rm = TRUE)) {
        warning("Negative values detected before log transformation. Setting to 0.")
        expr_lod[expr_lod < 0] <- 0
      }
      expr_log <- log2(expr_lod + 1)
    } else if (transform == "cuberoot_zscore") {
      expr_log <- sign(expr_lod) * abs(expr_lod)^(1 / 3)
    } else if (transform == "pareto") {
      if (any(expr_lod < 0, na.rm = TRUE)) {
        warning("Negative values detected before log transformation. Setting to 0.")
        expr_lod[expr_lod < 0] <- 0
      }
      expr_log <- log2(expr_lod + 1)
    }
    if (any(!is.finite(expr_log))) {
      warning("Non-finite values after transform")
    }
  } else if (log_transform) {
    if (verbose) cat("Applying log2(x + 1) transformation...\n")

    # Check for negative values
    if (any(expr_lod < 0, na.rm = TRUE)) {
      warning("Negative values detected before log transformation. Setting to 0.")
      expr_lod[expr_lod < 0] <- 0
    }

    expr_log <- log2(expr_lod + 1)

    # Check for infinities or NaNs
    if (any(!is.finite(expr_log))) {
      warning("Non-finite values after log transformation")
    }
  } else {
    expr_log <- expr_lod
  }
  
  # ===== SECTION 5.5: COVARIATE ADJUSTMENT =====
  covariates_adjusted <- FALSE
  n_samples_removed_covariates <- 0
  covariate_imputation_info <- list(n_imputed = 0, imputation_method = "none")
  covariates_to_adjust <- NULL
  balance_check_performed <- FALSE
  balance_results <- NULL
  covariates_tested <- NULL
  
  if (!is.null(covariates)) {
    # Validate that covariates are in metadata_cols
    missing_cov_meta <- setdiff(covariates, colnames(meta_raw))
    if (length(missing_cov_meta) > 0) {
      stop("Covariates not found in metadata: ", paste(missing_cov_meta, collapse = ", "),
           "\nEnsure covariate columns are included in metadata_cols parameter")
    }
    
    # Intelligent covariate adjustment: check group balance first
    if (!is.null(response_var) && check_balance) {
      if (verbose) cat("Checking group balance for covariates...\n")
      
      # Validate response variable exists
      if (!response_var %in% colnames(meta_raw)) {
        warning("response_var '", response_var, "' not found in metadata. Adjusting all covariates.")
        covariates_to_adjust <- covariates
      } else {
        # Extract response groups
        response_groups <- meta_raw[[response_var]]
        unique_groups <- unique(response_groups)
        
        if (length(unique_groups) < 2) {
          warning("Response variable has < 2 groups. Cannot check balance. Adjusting all covariates.")
          covariates_to_adjust <- covariates
        } else {
          # Check balance for each covariate
          balance_results <- data.frame(
            covariate = character(),
            group1 = character(),
            group2 = character(),
            mean_group1 = numeric(),
            mean_group2 = numeric(),
            p_value = numeric(),
            imbalanced = logical(),
            stringsAsFactors = FALSE
          )
          
          imbalanced_covariates <- c()
          covariates_tested <- covariates
          
          for (cov in covariates) {
            # For binary response, compare two groups
            if (length(unique_groups) == 2) {
              group1_vals <- meta_raw[[cov]][response_groups == unique_groups[1]]
              group2_vals <- meta_raw[[cov]][response_groups == unique_groups[2]]
              
              # Remove NAs
              group1_vals <- group1_vals[!is.na(group1_vals)]
              group2_vals <- group2_vals[!is.na(group2_vals)]
              
              if (length(group1_vals) > 0 && length(group2_vals) > 0) {
                # Check if covariate is categorical/factor
                if (is.factor(group1_vals) || is.character(group1_vals)) {
                  # For categorical variables, use chi-square test
                  contingency_table <- table(
                    c(rep(as.character(unique_groups[1]), length(group1_vals)),
                      rep(as.character(unique_groups[2]), length(group2_vals))),
                    c(as.character(group1_vals), as.character(group2_vals))
                  )
                  chi_result <- chisq.test(contingency_table)
                  p_val <- chi_result$p.value
                  
                  # Calculate mode (most common category) instead of mean
                  mode1 <- names(sort(table(group1_vals), decreasing = TRUE))[1]
                  mode2 <- names(sort(table(group2_vals), decreasing = TRUE))[1]
                  mean1 <- NA  # Not applicable for categorical
                  mean2 <- NA  # Not applicable for categorical
                  
                } else {
                  # For numeric variables, use t-test
                  t_result <- t.test(group1_vals, group2_vals)
                  p_val <- t_result$p.value
                  mean1 <- mean(group1_vals)
                  mean2 <- mean(group2_vals)
                }
                
                # Store results
                balance_results <- rbind(balance_results, data.frame(
                  covariate = cov,
                  group1 = as.character(unique_groups[1]),
                  group2 = as.character(unique_groups[2]),
                  mean_group1 = mean1,
                  mean_group2 = mean2,
                  p_value = p_val,
                  imbalanced = p_val < balance_threshold,
                  stringsAsFactors = FALSE
                ))
                
                if (verbose) {
                  cat(sprintf("  %s:\n", cov))
                  if (is.na(mean1)) {
                    # Categorical variable - show distribution
                    cat(sprintf("    %s: mode = %s (n = %d)\n", 
                               unique_groups[1], mode1, length(group1_vals)))
                    cat(sprintf("    %s: mode = %s (n = %d)\n", 
                               unique_groups[2], mode2, length(group2_vals)))
                    cat(sprintf("    chi-square test p = %.4f\n", p_val))
                  } else {
                    # Numeric variable - show means
                    cat(sprintf("    %s: mean = %.2f (n = %d)\n", 
                               unique_groups[1], mean1, length(group1_vals)))
                    cat(sprintf("    %s: mean = %.2f (n = %d)\n", 
                               unique_groups[2], mean2, length(group2_vals)))
                    cat(sprintf("    t-test p = %.4f\n", p_val))
                  }
                  
                  if (p_val < balance_threshold) {
                    cat(sprintf("    ! IMBALANCED (p < %.2f) - will adjust\n", balance_threshold))
                    imbalanced_covariates <- c(imbalanced_covariates, cov)
                  } else {
                    cat(sprintf("    v Balanced (p >= %.2f) - no adjustment\n", balance_threshold))
                  }
                }
              }
            } else {
              # For multi-group, use ANOVA (simplified for now - just adjust all)
              if (verbose) cat(sprintf("  %s: Multi-group response - adjusting without balance check\n", cov))
              imbalanced_covariates <- c(imbalanced_covariates, cov)
            }
          }
          
          balance_check_performed <- TRUE
          covariates_to_adjust <- imbalanced_covariates
          
          if (verbose) {
            cat("\n=== Balance Check Summary ===\n")
            if (length(imbalanced_covariates) > 0) {
              cat("  Covariates to adjust: ", paste(imbalanced_covariates, collapse = ", "), "\n")
            } else {
              cat("  No imbalanced covariates - skipping adjustment\n")
            }
            cat("=============================\n\n")
          }
        }
      }
    } else {
      # No balance check - adjust all covariates
      covariates_to_adjust <- covariates
      if (verbose) {
        cat("Adjusting for covariates: ", paste(covariates_to_adjust, collapse = ", "), "\n")
        cat("  (Balance check disabled or no response variable provided)\n")
      }
    }
    
    # Perform adjustment if any covariates selected
    if (!is.null(covariates_to_adjust) && length(covariates_to_adjust) > 0) {
      if (verbose && balance_check_performed) {
        cat("Performing covariate adjustment...\n")
      }
      
      # Call internal adjustment function
      adjustment_result <- adjust_for_covariates_internal(
        expr_mat = expr_log,
        metadata = meta_raw,
        covariates = covariates_to_adjust,
        subject_id_col = subject_id_col,
        verbose = verbose
      )
      
      # Extract adjusted data
      expr_adjusted <- adjustment_result$adjusted_matrix
      samples_kept <- adjustment_result$samples_kept
      n_samples_removed_covariates <- adjustment_result$n_removed
      covariate_imputation_info <- adjustment_result$imputed_info
      
      # Check if any covariates were skipped due to single level
      single_level_skipped <- if (!is.null(adjustment_result$single_level_covariates)) {
        adjustment_result$single_level_covariates
      } else {
        c()
      }
      
      # Update covariates_to_adjust to reflect what was actually adjusted
      if (!is.null(adjustment_result$covariates_adjusted)) {
        covariates_to_adjust <- adjustment_result$covariates_adjusted
      }
      
      # Update metadata and sample IDs to match retained samples
      if (n_samples_removed_covariates > 0) {
        meta_raw <- meta_raw[samples_kept, , drop = FALSE]
        sample_ids <- sample_ids[samples_kept]
      }
      
      # Use adjusted expression for downstream processing
      expr_log <- expr_adjusted
      covariates_adjusted <- TRUE
      
      if (verbose) {
        cat("  Covariate adjustment complete\n")
        if (length(covariates_to_adjust) > 0) {
          cat("  Adjusted for:", paste(covariates_to_adjust, collapse = ", "), "\n")
        }
        if (length(single_level_skipped) > 0) {
          cat("  Skipped (single level):", paste(single_level_skipped, collapse = ", "), "\n")
        }
        if (n_samples_removed_covariates > 0) {
          cat("  Samples retained:", nrow(expr_log), "(removed", n_samples_removed_covariates, "due to missing covariates)\n")
        }
      }
    } else if (!is.null(covariates) && balance_check_performed) {
      if (verbose) cat("  No covariate adjustment performed (all covariates balanced)\n")
    }
  }
  
  # ===== SECTION 6: SCALING =====
  if (!is.null(transform)) {
    if (transform %in% c("log2p1_zscore", "cuberoot_zscore")) {
      expr_scaled <- scale(expr_log)
      expr_final <- as.matrix(expr_scaled)
      rownames(expr_final) <- rownames(expr_log)
      colnames(expr_final) <- colnames(expr_log)
    } else if (transform == "pareto") {
      # Pareto scaling: center by mean, divide by sqrt(sd)
      col_means <- colMeans(expr_log, na.rm = TRUE)
      col_sds   <- apply(expr_log, 2, sd, na.rm = TRUE)
      col_sds[col_sds == 0 | !is.finite(col_sds)] <- 1
      expr_final <- sweep(expr_log, 2, col_means, "-")
      expr_final <- sweep(expr_final, 2, sqrt(col_sds), "/")
      expr_final <- as.matrix(expr_final)
      rownames(expr_final) <- rownames(expr_log)
      colnames(expr_final) <- colnames(expr_log)
    }
  } else if (scale_data) {
    if (verbose) cat("Applying z-score scaling (standardization)...\n")

    # Scale by column (each cytokine)
    expr_scaled <- scale(expr_log)

    # Convert back to matrix and preserve row/col names
    expr_final <- as.matrix(expr_scaled)
    rownames(expr_final) <- rownames(expr_log)
    colnames(expr_final) <- colnames(expr_log)

    # Check scaling
    if (verbose) {
      means <- colMeans(expr_final, na.rm = TRUE)
      sds <- apply(expr_final, 2, sd, na.rm = TRUE)
      if (any(abs(means) > 1e-10, na.rm = TRUE)) {
        warning("Scaling may not have centered data properly")
      }
    }
  } else {
    expr_final <- expr_log
  }
  
  # ===== SECTION 7: PREPARE OUTPUT =====
  if (verbose) cat("Finalizing preprocessed data...\n")
  
  # Create preprocessing info
  preprocess_info <- list(
    lod_method = lod_method,
    log_transformed = log_transform,
    scaled = scale_data,
    replicate_averaged = !is.null(replicate_col),
    replicate_col = replicate_col,
    covariates_adjusted = covariates_adjusted,
    covariates = if (covariates_adjusted) covariates_to_adjust else NULL,
    covariates_specified = if (!is.null(covariates)) covariates else NULL,
    covariate_balance_checked = balance_check_performed,
    covariates_tested = covariates_tested,
    balance_threshold = if (balance_check_performed) balance_threshold else NULL,
    balance_results = if (balance_check_performed) balance_results else NULL,
    covariate_imputation_method = covariate_imputation_info$imputation_method,
    n_imputed_covariate_values = covariate_imputation_info$n_imputed,
    n_samples_removed_for_covariates = n_samples_removed_covariates,
    subject_id_col = if (!is.null(subject_id_col)) subject_id_col else NULL,
    n_samples = nrow(expr_final),
    n_cytokines = ncol(expr_final),
    cytokine_names = colnames(expr_final),
    metadata_cols = colnames(meta_raw)
  )
  
  # Create result object
  result <- list(
    expression = expr_final,
    metadata = meta_raw,
    preprocessing_info = preprocess_info,
    sample_ids = sample_ids
  )
  
  # Add S3 class
  class(result) <- c("plsda_preprocessed", "list")
  
  # Add attributes for easy access
  attr(result, "lod_method") <- lod_method
  attr(result, "log_transformed") <- log_transform
  attr(result, "scaled") <- scale_data
  
  if (verbose) {
    cat("\nPreprocessing complete!\n")
    cat("  Samples: ", nrow(expr_final), "\n")
    cat("  Cytokines: ", ncol(expr_final), "\n")
    cat("  LOD method: ", lod_method, "\n")
    cat("  Log transformed: ", log_transform, "\n")
    if (covariates_adjusted) {
      cat("  Covariates adjusted: ", paste(covariates, collapse = ", "), "\n")
      if (covariate_imputation_info$n_imputed > 0) {
        cat("    Imputed values: ", covariate_imputation_info$n_imputed, "\n")
      }
      if (n_samples_removed_covariates > 0) {
        cat("    Samples removed: ", n_samples_removed_covariates, "\n")
      }
    }
    cat("  Scaled: ", scale_data, "\n")
  }
  
  return(result)
}


#' Print method for plsda_preprocessed objects
#'
#' @param x A plsda_preprocessed object
#' @param ... Additional arguments (not used)
#' @export
print.plsda_preprocessed <- function(x, ...) {
  cat("PLS-DA Preprocessed Data\n")
  cat("========================\n\n")
  cat("Samples:", x$preprocessing_info$n_samples, "\n")
  cat("Cytokines:", x$preprocessing_info$n_cytokines, "\n\n")
  cat("Preprocessing steps:\n")
  cat("  LOD method:", x$preprocessing_info$lod_method, "\n")
  cat("  Log transformed:", x$preprocessing_info$log_transformed, "\n")
  if (x$preprocessing_info$covariates_adjusted) {
    cat("  Covariates adjusted:", paste(x$preprocessing_info$covariates, collapse = ", "), "\n")
    if (x$preprocessing_info$covariate_balance_checked) {
      cat("    Balance check: Performed (threshold p <", x$preprocessing_info$balance_threshold, ")\n")
      if (!is.null(x$preprocessing_info$balance_results) && nrow(x$preprocessing_info$balance_results) > 0) {
        n_imbalanced <- sum(x$preprocessing_info$balance_results$imbalanced)
        n_tested <- nrow(x$preprocessing_info$balance_results)
        cat("    Tested:", n_tested, "covariates |", n_imbalanced, "imbalanced\n")
      }
    }
    if (x$preprocessing_info$n_imputed_covariate_values > 0) {
      cat("    Imputation:", x$preprocessing_info$covariate_imputation_method, 
          "(", x$preprocessing_info$n_imputed_covariate_values, "values )\n")
    }
    if (x$preprocessing_info$n_samples_removed_for_covariates > 0) {
      cat("    Samples removed:", x$preprocessing_info$n_samples_removed_for_covariates, "\n")
    }
  } else if (!is.null(x$preprocessing_info$covariates_specified) && x$preprocessing_info$covariate_balance_checked) {
    cat("  Covariates checked but balanced - no adjustment performed\n")
    if (!is.null(x$preprocessing_info$balance_results)) {
      cat("    Tested:", paste(x$preprocessing_info$covariates_tested, collapse = ", "), "\n")
    }
  }
  cat("  Scaled:", x$preprocessing_info$scaled, "\n")
  if (x$preprocessing_info$replicate_averaged) {
    cat("  Replicates averaged by:", x$preprocessing_info$replicate_col, "\n")
  }
  cat("\nMetadata columns:", paste(x$preprocessing_info$metadata_cols, collapse = ", "), "\n")
  invisible(x)
}


#' Summary method for plsda_preprocessed objects
#'
#' @param object A plsda_preprocessed object
#' @param ... Additional arguments (not used)
#' @export
summary.plsda_preprocessed <- function(object, ...) {
  cat("PLS-DA Preprocessed Data Summary\n")
  cat("=================================\n\n")
  
  print(object)
  
  cat("\nExpression matrix summary:\n")
  cat("  Range: [", 
      round(min(object$expression, na.rm = TRUE), 3), ", ",
      round(max(object$expression, na.rm = TRUE), 3), "]\n", sep = "")
  cat("  Mean:", round(mean(object$expression, na.rm = TRUE), 3), "\n")
  cat("  SD:", round(sd(object$expression, na.rm = TRUE), 3), "\n")
  cat("  Missing values:", sum(is.na(object$expression)), "\n")
  
  invisible(object)
}
