# Replicate Quality Control for Luminex/Multiplex Immunoassay Data
# Detects discordant technical replicates using CV and population-relative metrics

# Internal helper: default tuning parameters for sample-level well detection
well_params_defaults <- function() {
  list(
    well_threshold         = 0.3,   # Fraction of flagged analytes to classify as well failure
    sign_bias_threshold    = 0.7,   # |mean(sign(z))| above which failure is directional
    lod_fraction_threshold = 0.5,   # LOD fraction above which a well is definitively dropout
    score_ratio_threshold  = 1.5    # well_score ratio below which sign_bias is tiebreaker
  )
}

#' Flag discordant technical replicates in multiplex immunoassay data
#'
#' For data with technical replicates (typically duplicate Luminex wells),
#' identifies discordant replicate pairs using coefficient of variation (CV)
#' and/or population-relative deviation metrics. With n=2 replicates, the pair
#' itself is assessed for discordance rather than individual outliers.
#'
#' @importFrom dplyr group_by summarise mutate left_join filter ungroup across
#'   all_of first n if_else case_when select arrange
#' @importFrom stats median sd mad
#' @importFrom rlang sym .data
#'
#' @param data data.frame in long format with one row per replicate measurement.
#' @param sample_id_col Character name of column identifying biological samples
#'   (replicates share this ID). Default: \code{"original_sample_id"}.
#' @param cytokine_col Character name of column identifying the analyte.
#'   Default: \code{"cytokine"}.
#' @param value_col Character name of column containing concentration values
#'   (raw scale, not log-transformed). Default: \code{"concentration"}.
#' @param replicate_id_col Character name of column uniquely identifying each
#'   physical replicate well (e.g., \code{"sample_id"}). Required when
#'   \code{well_detection = "sample"}. Default: \code{NULL}.
#' @param group_cols Optional character vector of additional grouping columns
#'   (e.g., \code{"timepoint"}) for computing population-level statistics
#'   used in MAD-based flagging. Default: \code{NULL}.
#' @param cv_threshold Numeric CV threshold (%) above which a replicate pair
#'   is flagged. Standard Luminex QC uses 15--25%. Default: 25.
#' @param mad_threshold Numeric multiplier for MAD-based population-relative
#'   flagging. Pairs whose within-pair absolute difference exceeds
#'   \code{mad_threshold * MAD} of all within-pair differences for that analyte
#'   are flagged. Default: 3.
#' @param flag_rule Character specifying how to combine CV and MAD criteria.
#'   One of \code{"cv"} (CV only), \code{"mad"} (MAD only), \code{"either"}
#'   (flag if either triggers), or \code{"both"} (flag only if both trigger).
#'   Default: \code{"cv"}.
#' @param action Character specifying what to do with flagged pairs. One of
#'   \code{"flag"} (annotate only), \code{"drop_both"} (remove both replicates),
#'   \code{"drop_farther"} (remove the replicate farther from the population
#'   median), or \code{"winsorize"} (replace the farther replicate with the
#'   population median). Default: \code{"flag"}.
#' @param well_detection Character specifying the detection scope. \code{"analyte"}
#'   uses existing per-analyte logic. \code{"sample"} enables sample-level well
#'   failure detection: when a large fraction of analytes are flagged for a single
#'   sample, the function identifies which physical well failed and consistently
#'   applies the action across all flagged analytes for that sample.
#'   Default: \code{"analyte"}.
#' @param well_params Named list of tuning parameters for sample-level detection.
#'   Elements are merged onto defaults via \code{modifyList(well_params_defaults(), well_params)}.
#'   Available elements: \code{well_threshold} (fraction of analytes flagged to
#'   classify as well failure, default 0.3), \code{sign_bias_threshold} (directional
#'   bias cutoff, default 0.7), \code{lod_fraction_threshold} (LOD fraction cutoff
#'   for dropout, default 0.5), \code{score_ratio_threshold} (score ratio below which
#'   sign_bias is used as tiebreaker, default 1.5). Default: \code{list()}.
#' @param lod_col Character name of column containing LOD values per analyte.
#'   Optional; enables LOD-fraction scoring for dropout detection when
#'   \code{well_detection = "sample"}. Default: \code{NULL}.
#' @param min_reps Integer; samples with fewer than this many replicates are
#'   not evaluated for discordance. Default: 2L.
#' @param verbose Logical; if \code{TRUE}, prints a summary of flagged pairs.
#'   Default: \code{TRUE}.
#'
#' @return The input data.frame with additional columns:
#'   \describe{
#'     \item{replicate_n}{Number of replicates for this sample-analyte pair}
#'     \item{replicate_cv}{Within-pair CV (%) for each observation}
#'     \item{replicate_discordant}{Logical; \code{TRUE} for discordant pairs}
#'     \item{replicate_action}{Character describing the action taken
#'       (\code{"keep"}, \code{"dropped"}, or \code{"winsorized"})}
#'     \item{well_failure}{Logical; \code{TRUE} for rows in samples classified
#'       as well-level failures. Only populated when \code{well_detection = "sample"}.}
#'     \item{well_failure_type}{Character; \code{"inflation"}, \code{"dropout"},
#'       or \code{"mixed"}. \code{NA} for non-well-failure samples.}
#'     \item{well_score}{Numeric; per-replicate \code{median(abs(z))} across
#'       analytes. \code{NA} for non-well-failure samples.}
#'     \item{well_sign_bias}{Numeric; per-replicate \code{abs(mean(sign(z)))}.
#'       1.0 = perfectly systematic, 0.0 = random. \code{NA} when not applicable.}
#'     \item{well_lod_fraction}{Numeric; fraction of analytes at/below LOD per
#'       replicate. Only populated when \code{lod_col} is provided.}
#'   }
#'   When \code{action != "flag"}, rows may be removed or values modified.
#'
#' @details
#' \strong{CV metric}: For each replicate pair, \code{CV = SD / mean * 100}.
#' For n=2 this simplifies to \code{|x1 - x2| / mean(x1, x2) * 100}.
#' Pairs where the mean is effectively zero (\code{< 1e-10}) are excluded
#' from CV-based flagging to avoid division-by-zero artifacts.
#'
#' \strong{MAD metric}: The absolute difference within each pair is compared
#' to the population distribution of within-pair differences for that analyte
#' (optionally stratified by \code{group_cols}). A pair is flagged if its
#' difference exceeds \code{median + mad_threshold * MAD}.
#'
#' \strong{Actions for flagged pairs}:
#' \itemize{
#'   \item \code{"flag"}: Adds annotation columns only; no data modification
#'   \item \code{"drop_both"}: Removes all replicates for flagged pairs
#'   \item \code{"drop_farther"}: Keeps the replicate closest to the
#'     population median for that analyte-group; drops the other
#'   \item \code{"winsorize"}: Replaces the farther replicate's value with
#'     the population median for that analyte-group
#' }
#'
#' \strong{Sample-level well detection} (\code{well_detection = "sample"}):
#' When a large fraction of analytes are flagged for a single sample (exceeding
#' \code{well_params$well_threshold}), this indicates a systematic well failure
#' rather than analyte-level noise. The function computes per-well composite
#' scores (robust z-scores, directional sign bias, optional LOD fraction) to
#' identify the failed well. For these samples, the action is applied consistently
#' across all flagged analytes. Failure types are classified as \code{"inflation"}
#' (systematically high), \code{"dropout"} (systematically low), or \code{"mixed"}.
#'
#' \strong{Detection sensitivity}: Benchmark B6 shows that the default
#' \code{cv_threshold = 25} detects near-all gross contamination (replicate
#' offset \eqn{\geq} 1 SD) but only 0.74--0.96 of subtle outliers
#' (0.5 SD offset). If small-amplitude contamination is a concern, inspect the
#' \code{replicate_cv} distribution for heavy tails, consider a tighter
#' threshold (e.g., 15--20), or combine with MAD-based flagging via
#' \code{flag_rule = "either"}.
#'
#' @examples
#' \dontrun{
#' # Flag discordant duplicate wells
#' flagged <- flag_replicate_outliers(
#'   data = long_data,
#'   sample_id_col = "original_sample_id",
#'   value_col = "concentration",
#'   cv_threshold = 25
#' )
#'
#' # Sample-level well detection with drop_farther
#' cleaned <- flag_replicate_outliers(
#'   data = long_data,
#'   sample_id_col = "original_sample_id",
#'   replicate_id_col = "sample_id",
#'   value_col = "concentration",
#'   group_cols = "timepoint",
#'   cv_threshold = 25,
#'   action = "drop_farther",
#'   well_detection = "sample"
#' )
#' }
#'
#' @export
flag_replicate_outliers <- function(data,
                                    sample_id_col = "original_sample_id",
                                    cytokine_col = "cytokine",
                                    value_col = "concentration",
                                    replicate_id_col = NULL,
                                    group_cols = NULL,
                                    cv_threshold = 25,
                                    mad_threshold = 3,
                                    flag_rule = c("cv", "mad", "either", "both"),
                                    action = c("flag", "drop_both",
                                               "drop_farther", "winsorize"),
                                    well_detection = c("analyte", "sample"),
                                    well_params = list(),
                                    lod_col = NULL,
                                    min_reps = 2L,
                                    verbose = TRUE) {

  # --- Input validation ---
  flag_rule <- match.arg(flag_rule)
  action <- match.arg(action)
  well_detection <- match.arg(well_detection)

  stopifnot(is.data.frame(data))
  required_cols <- c(sample_id_col, cytokine_col, value_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  if (!is.null(group_cols)) {
    missing_group <- setdiff(group_cols, names(data))
    if (length(missing_group) > 0) {
      stop("Missing group column(s): ", paste(missing_group, collapse = ", "))
    }
  }
  stopifnot(is.numeric(cv_threshold), cv_threshold > 0)
  stopifnot(is.numeric(mad_threshold), mad_threshold > 0)
  stopifnot(is.numeric(min_reps), min_reps >= 2L)

  # Well detection validation
  if (well_detection == "sample") {
    if (is.null(replicate_id_col)) {
      stop("replicate_id_col is required when well_detection = 'sample'")
    }
    if (!replicate_id_col %in% names(data)) {
      stop("replicate_id_col '", replicate_id_col, "' not found in data")
    }
    if (!is.null(lod_col) && !lod_col %in% names(data)) {
      stop("lod_col '", lod_col, "' not found in data")
    }
    # Resolve well_params: merge user overrides onto defaults
    wp <- modifyList(well_params_defaults(), well_params)
    stopifnot(is.numeric(wp$well_threshold),
              wp$well_threshold > 0, wp$well_threshold <= 1)
    stopifnot(is.numeric(wp$sign_bias_threshold),
              wp$sign_bias_threshold > 0, wp$sign_bias_threshold <= 1)
    stopifnot(is.numeric(wp$lod_fraction_threshold),
              wp$lod_fraction_threshold > 0, wp$lod_fraction_threshold <= 1)
    stopifnot(is.numeric(wp$score_ratio_threshold),
              wp$score_ratio_threshold > 1)
  }

  # --- Step 1: Compute per-pair statistics ---
  pair_group_cols <- c(sample_id_col, cytokine_col)

  pair_stats <- data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(pair_group_cols))) %>%
    dplyr::summarise(
      pair_n = dplyr::n(),
      pair_mean = mean(.data[[value_col]], na.rm = TRUE),
      pair_sd = stats::sd(.data[[value_col]], na.rm = TRUE),
      pair_abs_diff = if (all(is.na(.data[[value_col]]))) NA_real_
                      else max(.data[[value_col]], na.rm = TRUE) -
                           min(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      # CV: guard against near-zero mean
      pair_cv = dplyr::if_else(
        pair_mean > 1e-10 & !is.na(pair_sd),
        pair_sd / pair_mean * 100,
        NA_real_
      ),
      # Replace NaN from single-rep sd with 0
      pair_sd = dplyr::if_else(is.na(pair_sd), 0, pair_sd),
      pair_abs_diff = dplyr::if_else(is.na(pair_abs_diff), 0, pair_abs_diff)
    )

  # --- Step 2: Join group_cols to pair_stats, then compute population stats ---
  pop_group_cols <- c(cytokine_col, group_cols)

  # Attach group_cols from original data before computing population stats
  if (!is.null(group_cols)) {
    pair_groups <- data %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(pair_group_cols))) %>%
      dplyr::summarise(
        dplyr::across(dplyr::all_of(group_cols), dplyr::first),
        .groups = "drop"
      )
    pair_stats <- pair_stats %>%
      dplyr::left_join(pair_groups, by = pair_group_cols)
  }

  # Population median of within-pair absolute differences (for MAD flagging)
  pop_diff_stats <- pair_stats %>%
    dplyr::filter(pair_n >= min_reps) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(pop_group_cols))) %>%
    dplyr::summarise(
      pop_median_diff = stats::median(pair_abs_diff, na.rm = TRUE),
      pop_mad_diff = stats::mad(pair_abs_diff, na.rm = TRUE),
      .groups = "drop"
    )

  # Population median of values (for drop_farther / winsorize)
  # When well_detection = "sample", also compute population_mad for z-scores
  if (well_detection == "sample") {
    pop_val_stats <- data %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(pop_group_cols))) %>%
      dplyr::summarise(
        population_median = stats::median(.data[[value_col]], na.rm = TRUE),
        population_mad = stats::mad(.data[[value_col]], na.rm = TRUE),
        population_sd = stats::sd(.data[[value_col]], na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    pop_val_stats <- data %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(pop_group_cols))) %>%
      dplyr::summarise(
        population_median = stats::median(.data[[value_col]], na.rm = TRUE),
        .groups = "drop"
      )
  }

  # --- Step 3: Join population stats to pair stats ---
  pair_stats <- pair_stats %>%
    dplyr::left_join(pop_diff_stats, by = pop_group_cols) %>%
    dplyr::left_join(pop_val_stats, by = pop_group_cols)

  # --- Step 4: Apply flagging rules ---
  pair_stats <- pair_stats %>%
    dplyr::mutate(
      # CV criterion: flag if CV exceeds threshold (skip when CV is NA/undefined)
      flag_cv = !is.na(pair_cv) & pair_cv > cv_threshold,
      # MAD criterion: flag if abs_diff exceeds population threshold
      # When MAD is 0 (all pairs identical), only flag if diff > 0
      mad_cutoff = dplyr::if_else(
        pop_mad_diff > 1e-10,
        pop_median_diff + mad_threshold * pop_mad_diff,
        pop_median_diff + 1e-10
      ),
      flag_mad = pair_abs_diff > mad_cutoff,
      # Combine based on rule
      replicate_discordant = dplyr::case_when(
        pair_n < min_reps ~ FALSE,
        flag_rule == "cv"     ~ flag_cv,
        flag_rule == "mad"    ~ flag_mad,
        flag_rule == "either" ~ flag_cv | flag_mad,
        flag_rule == "both"   ~ flag_cv & flag_mad,
        TRUE ~ FALSE
      )
    )

  # --- Step 4b: Sample-level well classification ---
  # Classify each sample as "well_failure" or "analyte_noise" based on
  # what fraction of its analytes are flagged
  if (well_detection == "sample") {
    sample_flag_summary <- pair_stats %>%
      dplyr::filter(pair_n >= min_reps) %>%
      dplyr::group_by(.data[[sample_id_col]]) %>%
      dplyr::summarise(
        n_total_analytes = dplyr::n(),
        n_flagged_analytes = sum(replicate_discordant, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        frac_flagged = n_flagged_analytes / n_total_analytes,
        well_failure = frac_flagged >= wp$well_threshold & n_total_analytes > 1
      )

    pair_stats <- pair_stats %>%
      dplyr::left_join(
        sample_flag_summary %>%
          dplyr::select(dplyr::all_of(sample_id_col), well_failure,
                        n_flagged_analytes, n_total_analytes, frac_flagged),
        by = sample_id_col
      ) %>%
      dplyr::mutate(
        well_failure = dplyr::if_else(is.na(well_failure), FALSE, well_failure)
      )
  }

  # Select columns to join back
  join_cols <- c(pair_group_cols,
                 "pair_n", "pair_cv", "replicate_discordant", "population_median")
  if (well_detection == "sample") {
    join_cols <- c(join_cols, "well_failure")
  }
  pair_join <- pair_stats %>%
    dplyr::select(dplyr::all_of(join_cols))

  # --- Step 5: Join flags back to original data ---
  result <- data %>%
    dplyr::left_join(pair_join, by = pair_group_cols) %>%
    dplyr::rename(replicate_n = pair_n, replicate_cv = pair_cv) %>%
    dplyr::mutate(
      replicate_discordant = dplyr::if_else(
        is.na(replicate_discordant), FALSE, replicate_discordant
      ),
      replicate_action = "keep"
    )

  if (well_detection == "sample") {
    result <- result %>%
      dplyr::mutate(
        well_failure = dplyr::if_else(is.na(well_failure), FALSE, well_failure)
      )
  }

  # --- Step 5b: Compute well scores for well-failure samples ---
  # This runs for ALL actions when well_detection = "sample", to populate

  # diagnostic columns (well_score, well_sign_bias, well_failure_type, etc.)
  well_scores <- NULL
  good_wells <- NULL
  bad_well_info <- NULL
  tied_samples <- character(0)

  if (well_detection == "sample") {
    well_failure_samples <- unique(
      result[[sample_id_col]][result$well_failure]
    )

    if (length(well_failure_samples) > 0) {
      # Compute robust z-scores for well-failure sample rows
      # result already has population_median from Step 5; only join mad/sd
      pop_scale_stats <- pop_val_stats %>%
        dplyr::select(dplyr::all_of(pop_group_cols), population_mad, population_sd)

      wf_rows <- result %>%
        dplyr::filter(.data[[sample_id_col]] %in% well_failure_samples) %>%
        dplyr::left_join(pop_scale_stats, by = pop_group_cols) %>%
        dplyr::mutate(
          .scale = dplyr::if_else(
            population_mad > 1e-10, population_mad,
            dplyr::if_else(population_sd > 1e-10, population_sd, NA_real_)
          ),
          .z = dplyr::if_else(
            is.na(.scale), 0,
            (.data[[value_col]] - population_median) / .scale
          )
        )

      # Compute per-well composite scores
      well_scores <- wf_rows %>%
        dplyr::group_by(.data[[sample_id_col]], .data[[replicate_id_col]]) %>%
        dplyr::summarise(
          well_score = stats::median(abs(.z), na.rm = TRUE),
          sign_bias = abs(mean(sign(.z), na.rm = TRUE)),
          mean_sign_z = mean(sign(.z), na.rm = TRUE),
          .groups = "drop"
        )

      # Optionally compute LOD fraction
      if (!is.null(lod_col)) {
        lod_fractions <- wf_rows %>%
          dplyr::group_by(.data[[sample_id_col]], .data[[replicate_id_col]]) %>%
          dplyr::summarise(
            lod_fraction = mean(.data[[value_col]] <= .data[[lod_col]],
                                na.rm = TRUE),
            .groups = "drop"
          )
        well_scores <- well_scores %>%
          dplyr::left_join(lod_fractions,
                           by = c(sample_id_col, replicate_id_col))
      } else {
        well_scores <- well_scores %>%
          dplyr::mutate(lod_fraction = NA_real_)
      }

      # Determine good well per sample
      good_wells <- well_scores %>%
        dplyr::group_by(.data[[sample_id_col]]) %>%
        dplyr::summarise(
          good_well = {
            scores <- well_score
            biases <- sign_bias
            lod_fracs <- lod_fraction
            well_ids <- .data[[replicate_id_col]]
            n_w <- length(well_ids)

            if (n_w < 2) {
              well_ids[1]
            } else {
              sr <- if (min(scores) > 1e-10) max(scores) / min(scores) else Inf

              if (!is.null(lod_col) && !all(is.na(lod_fracs))) {
                high_lod <- lod_fracs > wp$lod_fraction_threshold
                low_lod <- !high_lod
                if (any(high_lod, na.rm = TRUE) && any(low_lod, na.rm = TRUE)) {
                  well_ids[which.min(lod_fracs)]
                } else if (sr >= wp$score_ratio_threshold) {
                  well_ids[which.min(scores)]
                } else if (length(unique(round(biases, 6))) > 1) {
                  well_ids[which.min(biases)]
                } else {
                  NA_character_
                }
              } else if (sr >= wp$score_ratio_threshold) {
                well_ids[which.min(scores)]
              } else if (length(unique(round(biases, 6))) > 1) {
                well_ids[which.min(biases)]
              } else {
                NA_character_
              }
            }
          },
          .groups = "drop"
        )

      # Classify failure type based on the bad well's z-scores
      bad_well_info <- well_scores %>%
        dplyr::left_join(good_wells, by = sample_id_col) %>%
        dplyr::filter(!is.na(good_well) &
                      .data[[replicate_id_col]] != good_well) %>%
        dplyr::group_by(.data[[sample_id_col]]) %>%
        dplyr::summarise(
          well_failure_type = dplyr::case_when(
            sign_bias[1] > wp$sign_bias_threshold &
              mean_sign_z[1] > 0 ~ "inflation",
            sign_bias[1] > wp$sign_bias_threshold &
              mean_sign_z[1] < 0 ~ "dropout",
            TRUE ~ "mixed"
          ),
          .groups = "drop"
        )

      tied_samples <- good_wells %>%
        dplyr::filter(is.na(good_well)) %>%
        dplyr::pull(.data[[sample_id_col]])

      # Join well diagnostic columns back to result
      well_diag <- well_scores %>%
        dplyr::select(dplyr::all_of(c(sample_id_col, replicate_id_col)),
                      well_score, well_sign_bias = sign_bias,
                      well_lod_fraction = lod_fraction)

      result <- result %>%
        dplyr::left_join(well_diag,
                         by = c(sample_id_col, replicate_id_col)) %>%
        dplyr::left_join(bad_well_info, by = sample_id_col) %>%
        dplyr::left_join(good_wells, by = sample_id_col)

      # For tied samples, set failure type to "mixed"
      if (length(tied_samples) > 0) {
        result <- result %>%
          dplyr::mutate(
            well_failure_type = dplyr::if_else(
              .data[[sample_id_col]] %in% tied_samples & is.na(well_failure_type),
              "mixed", well_failure_type
            )
          )
      }
    } else {
      # No well-failure samples: add empty diagnostic columns
      result <- result %>%
        dplyr::mutate(
          well_failure_type = NA_character_,
          well_score = NA_real_,
          well_sign_bias = NA_real_,
          well_lod_fraction = NA_real_,
          good_well = NA_character_
        )
    }
  }

  # --- Step 6: Apply action on flagged pairs ---
  if (action == "drop_both") {
    result <- result %>%
      dplyr::mutate(
        replicate_action = dplyr::if_else(replicate_discordant, "dropped", "keep")
      ) %>%
      dplyr::filter(!replicate_discordant)

  } else if (action == "drop_farther" || action == "winsorize") {

    has_well_failures <- well_detection == "sample" &&
      exists("well_failure_samples") && length(well_failure_samples) > 0

    if (has_well_failures) {
      # --- Well-failure samples: consistent well-based action ---
      # Build a boolean mask for resolvable well-failure rows
      is_resolvable_wf <- result$well_failure &
        !is.na(result$good_well) &
        !(result[[sample_id_col]] %in% tied_samples)

      if (any(is_resolvable_wf)) {
        # Store the mask as a column for use inside dplyr
        result$.is_rwf <- is_resolvable_wf

        if (action == "drop_farther") {
          result <- result %>%
            dplyr::mutate(
              replicate_action = dplyr::case_when(
                .is_rwf & !replicate_discordant ~ "keep",
                .is_rwf & replicate_discordant &
                  .data[[replicate_id_col]] == good_well ~ "keep",
                .is_rwf & replicate_discordant ~ "dropped",
                TRUE ~ replicate_action
              )
            )
        } else if (action == "winsorize") {
          is_bad_well_flagged <- is_resolvable_wf &
            result$replicate_discordant &
            result[[replicate_id_col]] != result$good_well

          result[[value_col]][is_bad_well_flagged] <-
            result$population_median[is_bad_well_flagged]

          result$replicate_action[is_bad_well_flagged] <- "winsorized"
        }
      } else {
        is_resolvable_wf <- rep(FALSE, nrow(result))
        result$.is_rwf <- FALSE
      }

      # Non-resolvable rows: per-analyte logic
      if (any(!is_resolvable_wf & result$replicate_discordant)) {
        if (action == "drop_farther") {
          result <- result %>%
            dplyr::mutate(
              dist_to_median = abs(.data[[value_col]] - population_median)
            ) %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(pair_group_cols))) %>%
            dplyr::mutate(
              replicate_action = dplyr::case_when(
                .is_rwf ~ replicate_action,
                !replicate_discordant ~ "keep",
                dist_to_median == min(dist_to_median) &
                  cumsum(dist_to_median == min(dist_to_median)) == 1 ~ "keep",
                TRUE ~ "dropped"
              )
            ) %>%
            dplyr::ungroup() %>%
            dplyr::select(-dist_to_median)
        } else if (action == "winsorize") {
          result <- result %>%
            dplyr::mutate(
              dist_to_median = abs(.data[[value_col]] - population_median)
            ) %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(pair_group_cols))) %>%
            dplyr::mutate(
              is_farther = !.is_rwf & replicate_discordant &
                dist_to_median == max(dist_to_median) &
                cumsum(dist_to_median == max(dist_to_median) &
                       !.is_rwf & replicate_discordant) == 1
            ) %>%
            dplyr::ungroup()

          result[[value_col]][result$is_farther] <-
            result$population_median[result$is_farther]

          result$replicate_action[result$is_farther] <- "winsorized"

          result <- result %>%
            dplyr::select(-dist_to_median, -is_farther)
        }
      }

      result <- result %>% dplyr::select(-`.is_rwf`)

      # Filter dropped rows for drop_farther
      if (action == "drop_farther") {
        result <- result %>%
          dplyr::filter(replicate_action != "dropped")
      }

    } else {
      # --- Per-analyte logic (well_detection = "analyte" or no well failures) ---
      if (action == "drop_farther") {
        result <- result %>%
          dplyr::mutate(
            dist_to_median = abs(.data[[value_col]] - population_median)
          ) %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(pair_group_cols))) %>%
          dplyr::mutate(
            replicate_action = dplyr::case_when(
              !replicate_discordant ~ "keep",
              dist_to_median == min(dist_to_median) &
                cumsum(dist_to_median == min(dist_to_median)) == 1 ~ "keep",
              TRUE ~ "dropped"
            )
          ) %>%
          dplyr::ungroup() %>%
          dplyr::filter(replicate_action != "dropped") %>%
          dplyr::select(-dist_to_median)

      } else if (action == "winsorize") {
        result <- result %>%
          dplyr::mutate(
            dist_to_median = abs(.data[[value_col]] - population_median)
          ) %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(pair_group_cols))) %>%
          dplyr::mutate(
            is_farther = replicate_discordant &
              dist_to_median == max(dist_to_median) &
              cumsum(dist_to_median == max(dist_to_median) &
                     replicate_discordant) == 1
          ) %>%
          dplyr::ungroup()

        result[[value_col]][result$is_farther] <-
          result$population_median[result$is_farther]

        result <- result %>%
          dplyr::mutate(
            replicate_action = dplyr::if_else(is_farther, "winsorized", "keep")
          ) %>%
          dplyr::select(-dist_to_median, -is_farther)
      }
    }
  }

  # Clean up temporary columns
  result <- result %>%
    dplyr::select(-dplyr::any_of(c("population_median", "population_mad",
                                    "population_sd", "good_well")))

  # --- Step 7: Report ---
  if (verbose) {
    n_pairs <- pair_stats %>%
      dplyr::filter(pair_n >= min_reps) %>%
      nrow()
    n_flagged <- sum(pair_stats$replicate_discordant, na.rm = TRUE)

    message(sprintf("Replicate QC: %d/%d pairs flagged (%.1f%%)",
                    n_flagged, n_pairs,
                    if (n_pairs > 0) 100 * n_flagged / n_pairs else 0))

    if (n_flagged > 0) {
      # Per-cytokine summary
      per_cyt <- pair_stats %>%
        dplyr::filter(pair_n >= min_reps) %>%
        dplyr::group_by(.data[[cytokine_col]]) %>%
        dplyr::summarise(
          n_pairs = dplyr::n(),
          n_flagged = sum(replicate_discordant),
          median_cv = round(stats::median(pair_cv, na.rm = TRUE), 1),
          .groups = "drop"
        ) %>%
        dplyr::filter(n_flagged > 0) %>%
        dplyr::arrange(dplyr::desc(n_flagged))

      message(sprintf("  Action: %s | Rule: %s | CV threshold: %.0f%% | MAD threshold: %.1f",
                      action, flag_rule, cv_threshold, mad_threshold))
      for (i in seq_len(nrow(per_cyt))) {
        message(sprintf("  %s: %d/%d flagged (median CV: %.1f%%)",
                        per_cyt[[cytokine_col]][i],
                        per_cyt$n_flagged[i],
                        per_cyt$n_pairs[i],
                        per_cyt$median_cv[i]))
      }
    }

    # Well-failure summary
    if (well_detection == "sample" && exists("sample_flag_summary") &&
        any(sample_flag_summary$well_failure, na.rm = TRUE)) {
      wf_samples <- sample_flag_summary %>%
        dplyr::filter(well_failure)

      message(sprintf("\n  Well-level failures detected: %d samples",
                      nrow(wf_samples)))

      # Get well scores and failure types for verbose output
      if (exists("well_scores") && exists("good_wells") &&
          exists("bad_well_info")) {
        for (s in wf_samples[[sample_id_col]]) {
          nf <- wf_samples$n_flagged_analytes[wf_samples[[sample_id_col]] == s]
          nt <- wf_samples$n_total_analytes[wf_samples[[sample_id_col]] == s]

          ft <- bad_well_info$well_failure_type[
            bad_well_info[[sample_id_col]] == s]
          if (length(ft) == 0) ft <- "mixed"

          gw <- good_wells$good_well[good_wells[[sample_id_col]] == s]

          # Get scores for this sample
          s_scores <- well_scores %>%
            dplyr::filter(.data[[sample_id_col]] == s)

          if (!is.na(gw) && nrow(s_scores) >= 2) {
            good_row <- s_scores[s_scores[[replicate_id_col]] == gw, ]
            bad_row <- s_scores[s_scores[[replicate_id_col]] != gw, ]

            bad_info <- sprintf("score: %.1f, bias: %.2f",
                                bad_row$well_score[1], bad_row$sign_bias[1])
            if (!is.null(lod_col) && !is.na(bad_row$lod_fraction[1])) {
              bad_info <- sprintf("%s, LOD: %.0f%%",
                                 bad_info, bad_row$lod_fraction[1] * 100)
            }

            message(sprintf(
              "    %s: %d/%d analytes flagged [%s] - good well: %s (score: %.1f, bias: %.2f), bad well: %s (%s)",
              s, nf, nt, ft, gw,
              good_row$well_score[1], good_row$sign_bias[1],
              bad_row[[replicate_id_col]][1], bad_info
            ))
          } else {
            message(sprintf("    %s: %d/%d analytes flagged [%s] - tied wells, per-analyte fallback",
                            s, nf, nt, ft))
          }
        }
      }
    }
  }

  result
}
