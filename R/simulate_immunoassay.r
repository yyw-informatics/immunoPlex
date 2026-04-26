# Simulation engine for immunoPlex: S3 class, methods, and analyte library
#
# immuno_sim objects hold $data (observed tibble), $truth (ground-truth tibble),
# and $meta (simulation parameters + realized summaries).

# ---- Built-in analyte library ------------------------------------------------

#' @noRd
analyte_library <- data.frame(
 analyte   = c("IL-1b", "IL-2", "IL-4", "IL-6", "IL-8", "IL-10", "IL-12p70",
                "IL-13", "IL-17A", "TNF-a", "IFN-g", "MCP-1", "MIP-1a",
                "MIP-1b", "RANTES", "IP-10", "Eotaxin", "G-CSF", "VEGF", "EGF"),
  log_mean = c(1.5, 1.0, 1.2, 2.5, 4.5, 2.0, 1.8, 1.5, 1.8, 2.8,
               2.0, 5.0, 3.5, 4.0, 6.0, 5.5, 4.0, 2.5, 3.5, 3.0),
  log_sd   = c(1.2, 1.5, 1.4, 1.3, 0.8, 1.2, 1.3, 1.1, 1.0, 1.0,
               1.5, 0.7, 1.0, 0.9, 0.6, 0.8, 0.7, 1.2, 1.0, 1.1),
  lod      = c(1.0, 2.0, 1.5, 0.5, 1.0, 1.0, 2.0, 1.0, 1.5, 0.5,
               2.0, 1.0, 2.0, 1.0, 5.0, 2.0, 3.0, 1.0, 5.0, 10.0),
  category = c("Pro-inflammatory", "Th1", "Th2", "Pro-inflammatory",
               "Chemokine", "Anti-inflammatory", "Th1", "Th2", "Th17",
               "Pro-inflammatory", "Th1", "Chemokine", "Chemokine",
               "Chemokine", "Chemokine", "Chemokine", "Chemokine",
               "Growth factor", "Growth factor", "Growth factor"),
  stringsAsFactors = FALSE
)


# ---- S3 constructor ----------------------------------------------------------

#' Create an immuno_sim object
#'
#' Constructor for the \code{immuno_sim} S3 class. Validates inputs and
#' assembles the three-component list (\code{$data}, \code{$truth},
#' \code{$meta}).
#'
#' @param data A tibble (or data.frame) of observed data in long format. Must
#'   contain at least columns \code{subject_id}, \code{cytokine}, and
#'   \code{value}.
#' @param truth A tibble (or data.frame) of ground-truth values with the same
#'   number of rows as \code{data}. Must contain at least column
#'   \code{true_latent}.
#' @param meta A named list of simulation metadata. Must contain at least
#'   \code{design} (character) and \code{n_subjects} (integer/numeric).
#'
#' @return An object of class \code{"immuno_sim"}.
#'
#' @details
#' The constructor enforces:
#' \itemize{
#'   \item \code{data} and \code{truth} must be data frames with identical row
#'     counts.
#'   \item \code{data} must contain columns \code{subject_id}, \code{cytokine},
#'     and \code{value}.
#'   \item \code{truth} must contain column \code{true_latent}.
#'   \item \code{meta} must be a list containing \code{design} (character) and
#'     \code{n_subjects} (positive integer).
#' }
#'
#' @examples
#' \dontrun{
#' sim <- new_immuno_sim(data = sim_data, truth = sim_truth, meta = sim_meta)
#' }
#'
#' @export
new_immuno_sim <- function(data, truth, meta) {
  # --- Validate data ----------------------------------------------------------

if (!is.data.frame(data)) {
    stop("`data` must be a data.frame or tibble.", call. = FALSE)
  }
  required_data_cols <- c("subject_id", "cytokine", "value")
  missing_data <- setdiff(required_data_cols, names(data))
  if (length(missing_data) > 0) {
    stop("`data` is missing required columns: ",
         paste(missing_data, collapse = ", "), call. = FALSE)
  }

  # --- Validate truth ---------------------------------------------------------
  if (!is.data.frame(truth)) {
    stop("`truth` must be a data.frame or tibble.", call. = FALSE)
  }
  if (!"true_latent" %in% names(truth)) {
    stop("`truth` must contain column `true_latent`.", call. = FALSE)
  }

  # --- Row alignment ----------------------------------------------------------
  if (nrow(data) != nrow(truth)) {
    stop("`data` and `truth` must have the same number of rows (got ",
         nrow(data), " vs ", nrow(truth), ").", call. = FALSE)
  }

  # --- Validate meta ----------------------------------------------------------
  if (!is.list(meta) || is.data.frame(meta)) {
    stop("`meta` must be a named list.", call. = FALSE)
  }
  if (is.null(meta$design) || !is.character(meta$design) ||
      length(meta$design) != 1) {
    stop("`meta$design` must be a single character string.", call. = FALSE)
  }
  if (is.null(meta$n_subjects) || !is.numeric(meta$n_subjects) ||
      length(meta$n_subjects) != 1 || meta$n_subjects < 1) {
    stop("`meta$n_subjects` must be a positive integer.", call. = FALSE)
  }

  # Coerce to tibble if plain data.frame
  if (!inherits(data, "tbl_df")) {
    data <- tibble::as_tibble(data)
  }
  if (!inherits(truth, "tbl_df")) {
    truth <- tibble::as_tibble(truth)
  }

  structure(
    list(data = data, truth = truth, meta = meta),
    class = "immuno_sim"
  )
}


# ---- Print method ------------------------------------------------------------

#' @title Print method for immuno_sim
#' @description One-line summary of a simulated immunoassay dataset.
#' @param x An \code{immuno_sim} object.
#' @param ... Unused.
#' @return \code{x} invisibly.
#' @export
print.immuno_sim <- function(x, ...) {
  n_sub <- length(unique(x$data$subject_id))
  n_cyt <- length(unique(x$data$cytokine))
  n_obs <- nrow(x$data)

  n_signal <- sum(x$truth$signal_analyte, na.rm = TRUE)
  signal_str <- if ("signal_analyte" %in% names(x$truth)) {
    n_unique_signal <- length(unique(
      x$data$cytokine[x$truth$signal_analyte]
    ))
    sprintf(", %d signal analytes", n_unique_signal)
  } else {
    ""
  }

  cat(sprintf("immuno_sim: %s design | %d subjects x %d analytes (%d obs%s)\n",
              x$meta$design, n_sub, n_cyt, n_obs, signal_str))
  invisible(x)
}


# ---- Summary method ----------------------------------------------------------

#' @title Summary method for immuno_sim
#' @description Detailed summary of a simulated immunoassay dataset including
#'   design, dimensions, censoring, and realized vs. target summaries.
#' @param object An \code{immuno_sim} object.
#' @param ... Unused.
#' @return Invisibly returns a list of summary statistics.
#' @export
summary.immuno_sim <- function(object, ...) {
  d <- object$data
  tr <- object$truth
  m <- object$meta

  n_sub <- length(unique(d$subject_id))
  n_cyt <- length(unique(d$cytokine))
  n_tp  <- if ("timepoint" %in% names(d)) length(unique(d$timepoint)) else 1L
  n_obs <- nrow(d)

  cat("Simulated Immunoassay Summary\n")
  cat("=============================\n\n")
  cat("Design:      ", m$design, "\n")
  cat("Subjects:    ", n_sub, "\n")
  cat("Analytes:    ", n_cyt, "\n")
  cat("Timepoints:  ", n_tp, "\n")
  cat("Observations:", n_obs, "\n")

  # Groups
  if ("group" %in% names(d)) {
    grp_tbl <- table(d$group[!duplicated(d$subject_id)])
    cat("Groups:      ", paste(names(grp_tbl), paste0("(n=", grp_tbl, ")"),
                               collapse = ", "), "\n")
  }

  # Signal analytes
  if ("signal_analyte" %in% names(tr)) {
    sig_analytes <- unique(d$cytokine[tr$signal_analyte])
    cat("Signal:      ", length(sig_analytes), "analytes",
        if (length(sig_analytes) > 0 && length(sig_analytes) <= 10)
          paste0("(", paste(sig_analytes, collapse = ", "), ")")
        else "",
        "\n")
  }

  cat("\n")

  # Censoring summary
  if ("cens_lod" %in% names(d)) {
    overall_lod <- mean(d$cens_lod, na.rm = TRUE)
    cat(sprintf("Censoring (LOD):  %.1f%% overall\n", overall_lod * 100))
    if ("cens_ulod" %in% names(d)) {
      overall_ulod <- mean(d$cens_ulod, na.rm = TRUE)
      cat(sprintf("Censoring (ULOD): %.1f%% overall\n", overall_ulod * 100))
    }

    # Per-analyte censoring
    if (!is.null(m$realized_censoring)) {
      cat("\nPer-analyte censoring (LOD):\n")
      rc <- sort(m$realized_censoring, decreasing = TRUE)
      rc_display <- rc[rc > 0]
      if (length(rc_display) > 0) {
        for (nm in names(rc_display)) {
          cat(sprintf("  %-12s %.1f%%\n", nm, rc_display[nm] * 100))
        }
        n_zero <- sum(rc == 0)
        if (n_zero > 0) cat(sprintf("  (%d analytes with 0%% censoring)\n", n_zero))
      } else {
        cat("  All analytes: 0%\n")
      }
    }
  }

  # Seed and timestamp
  cat("\n")
  if (!is.null(m$seed)) cat("Seed:        ", m$seed, "\n")
  if (!is.null(m$timestamp)) cat("Created:     ", format(m$timestamp), "\n")

  invisible(list(
    design = m$design,
    n_subjects = n_sub,
    n_analytes = n_cyt,
    n_timepoints = n_tp,
    n_obs = n_obs
  ))
}


# ---- Subset method -----------------------------------------------------------

#' @title Subset an immuno_sim object
#' @description Subsets \code{$data} and \code{$truth} in sync, preserving
#'   alignment and class.
#' @param x An \code{immuno_sim} object.
#' @param i Row indices (integer, logical, or missing).
#' @param j Column indices for \code{$data} (integer, character, or missing).
#'   \code{$truth} is always subset by rows only.
#' @param ... Unused.
#' @param drop Ignored (always FALSE for tibbles).
#' @return A new \code{immuno_sim} object with subset data and truth.
#' @export
`[.immuno_sim` <- function(x, i, j, ..., drop = FALSE) {
  d <- x$data
  tr <- x$truth
  m <- x$meta

  if (!missing(i)) {
    d  <- d[i, , drop = FALSE]
    tr <- tr[i, , drop = FALSE]
  }

  if (!missing(j)) {
    d <- d[, j, drop = FALSE]
  }

  # Update meta to reflect subset
  m$n_subjects <- length(unique(d$subject_id))
  if ("cytokine" %in% names(d)) {
    m$n_analytes <- length(unique(d$cytokine))
  }

  structure(
    list(data = d, truth = tr, meta = m),
    class = "immuno_sim"
  )
}


# ---- Stage 1: Latent biology generator ---------------------------------------

#' @noRd
simulate_latent_biology <- function(
    n_subjects = 60,
    n_timepoints = 2,
    n_analytes = 20,
    design = "pre_post",
    group_levels = c("control", "treatment"),
    group_allocation = NULL,
    group_effects = 0.5,
    time_effects = NULL,
    interaction_effects = NULL,
    signal_analytes = 5,
    effect_direction = "mixed",
    re_intercept_sd = 0.5,
    re_slope_sd = 0,
    residual_sd = 1.0,
    residual_ar1 = 0,
    heteroscedastic = FALSE,
    analyte_means = NULL,
    analyte_sds = NULL,
    analyte_correlation = NULL,
    block_rho = 0.5,
    covariates = NULL,
    covariate_group_cor = 0,
    dropout_rate = 0,
    dropout_informative = FALSE,
    visit_missingness = 0,
    seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  ## ---- Analyte setup --------------------------------------------------------
  lib <- analyte_library
  if (n_analytes <= nrow(lib)) {
    analyte_names <- lib$analyte[seq_len(n_analytes)]
    lib_idx <- seq_len(n_analytes)
  } else {
    analyte_names <- paste0("Analyte_", seq_len(n_analytes))
    lib_idx <- ((seq_len(n_analytes) - 1L) %% nrow(lib)) + 1L
  }

  mu <- if (is.null(analyte_means)) {
    stats::setNames(lib$log_mean[lib_idx], analyte_names)
  } else {
    stats::setNames(analyte_means[seq_len(n_analytes)], analyte_names)
  }

  sigma <- if (is.null(analyte_sds)) {
    stats::setNames(lib$log_sd[lib_idx], analyte_names)
  } else {
    stats::setNames(analyte_sds[seq_len(n_analytes)], analyte_names)
  }

  ## ---- Signal analyte selection ---------------------------------------------
  if (is.logical(signal_analytes)) {
    signal_mask <- signal_analytes
  } else if (is.character(signal_analytes)) {
    signal_mask <- analyte_names %in% signal_analytes
  } else {
    n_sig <- min(as.integer(signal_analytes), n_analytes)
    signal_mask <- rep(FALSE, n_analytes)
    if (n_sig > 0) {
      signal_mask[sort(sample.int(n_analytes, n_sig))] <- TRUE
    }
  }
  n_signal <- sum(signal_mask)
  signal_idx <- which(signal_mask)

  ## ---- Group assignment -----------------------------------------------------
  n_groups <- length(group_levels)
  if (is.null(group_allocation)) {
    group_allocation <- rep(1, n_groups)
  }
  group_probs <- group_allocation / sum(group_allocation)
  group_sizes <- rep(0L, n_groups)
  remaining <- n_subjects
  for (g in seq_len(n_groups - 1L)) {
    group_sizes[g] <- round(n_subjects * group_probs[g])
    remaining <- remaining - group_sizes[g]
  }
  group_sizes[n_groups] <- remaining

  subject_ids <- factor(paste0("S", sprintf("%03d", seq_len(n_subjects))))
  subject_groups <- factor(
    rep(group_levels, times = group_sizes),
    levels = group_levels
  )
  group_int <- as.integer(subject_groups)

  ## ---- Timepoint setup ------------------------------------------------------
  tp_labels <- factor(paste0("T", seq_len(n_timepoints)),
                      levels = paste0("T", seq_len(n_timepoints)))

  ## ---- Effect vectors -------------------------------------------------------
  group_eff_full <- rep(0, n_analytes)
  time_eff_full  <- rep(0, n_analytes)
  int_eff_full   <- rep(0, n_analytes)

  if (n_signal > 0) {
    g_eff <- if (length(group_effects) == 1) {
      rep(group_effects, n_signal)
    } else {
      group_effects[seq_len(n_signal)]
    }

    # Direction
    if (effect_direction == "mixed" && n_signal > 1) {
      n_down <- floor(n_signal / 2)
      dirs <- c(rep(-1, n_down), rep(1, n_signal - n_down))
      dirs <- sample(dirs)
      g_eff <- abs(g_eff) * dirs
    } else if (effect_direction == "down") {
      g_eff <- -abs(g_eff)
    }

    # Time effects
    t_eff <- if (!is.null(time_effects)) {
      if (length(time_effects) == 1) rep(time_effects, n_signal)
      else time_effects[seq_len(n_signal)]
    } else {
      rep(0, n_signal)
    }

    # Interaction effects
    i_eff <- if (!is.null(interaction_effects)) {
      if (length(interaction_effects) == 1) rep(interaction_effects, n_signal)
      else interaction_effects[seq_len(n_signal)]
    } else {
      rep(0, n_signal)
    }

    group_eff_full[signal_idx] <- g_eff
    time_eff_full[signal_idx]  <- t_eff
    int_eff_full[signal_idx]   <- i_eff
  }

  ## ---- Covariates (subject-level) -------------------------------------------
  cov_values  <- list()
  cov_effects <- list()
  if (!is.null(covariates) && length(covariates) > 0) {
    for (nm in names(covariates)) {
      spec <- covariates[[nm]]
      eff  <- if (is.null(spec$effect)) 0 else spec$effect
      tp   <- if (is.null(spec$type))   "continuous" else spec$type
      cov_effects[[nm]] <- eff

      z <- stats::rnorm(n_subjects)
      if (covariate_group_cor != 0 && n_groups >= 2) {
        g_sd <- stats::sd(group_int)
        if (g_sd > 0) {
          g_std <- (group_int - mean(group_int)) / g_sd
        } else {
          g_std <- rep(0, n_subjects)
        }
        latent <- covariate_group_cor * g_std +
          sqrt(1 - covariate_group_cor^2) * z
      } else {
        latent <- z
      }

      if (tp == "binary") {
        cov_values[[nm]] <- as.integer(latent > 0)
      } else {
        cov_values[[nm]] <- latent
      }
    }
  }

  ## ---- Analyte correlation matrix -------------------------------------------
  cor_mat <- diag(n_analytes)
  if (!is.null(analyte_correlation)) {
    if (is.character(analyte_correlation) && length(analyte_correlation) == 1) {
      if (analyte_correlation == "block") {
        cats <- lib$category[lib_idx]
        for (cat in unique(cats)) {
          idx <- which(cats == cat)
          if (length(idx) > 1) {
            for (a in idx) for (b in idx) {
              if (a != b) cor_mat[a, b] <- block_rho
            }
          }
        }
      } else if (analyte_correlation == "toeplitz") {
        for (a in seq_len(n_analytes)) for (b in seq_len(n_analytes)) {
          cor_mat[a, b] <- block_rho^abs(a - b)
        }
      }
    } else if (is.matrix(analyte_correlation)) {
      cor_mat <- analyte_correlation
    }
  }
  chol_cor <- chol(cor_mat)

  ## ---- Residual SD vector ---------------------------------------------------
  if (length(residual_sd) == 1) {
    resid_sd <- rep(residual_sd, n_analytes)
  } else {
    resid_sd <- residual_sd[seq_len(n_analytes)]
  }
  names(resid_sd) <- analyte_names

  ## ---- Random intercepts (n_subjects x n_analytes, correlated) --------------
  re_int_mat <- matrix(stats::rnorm(n_subjects * n_analytes),
                       nrow = n_subjects, ncol = n_analytes)
  re_int_mat <- re_int_mat %*% chol_cor * re_intercept_sd

  ## ---- Random slopes (optional) ---------------------------------------------
  re_slp_mat <- NULL
  if (re_slope_sd > 0) {
    re_slp_mat <- matrix(stats::rnorm(n_subjects * n_analytes),
                         nrow = n_subjects, ncol = n_analytes)
    re_slp_mat <- re_slp_mat %*% chol_cor * re_slope_sd
  }

  ## ---- Residuals (3-D array: subject x timepoint x analyte) -----------------
  eps <- array(stats::rnorm(n_subjects * n_timepoints * n_analytes),
               dim = c(n_subjects, n_timepoints, n_analytes))

  # AR(1) within subject-analyte across timepoints
  if (residual_ar1 != 0 && n_timepoints > 1) {
    innov_sd <- sqrt(1 - residual_ar1^2)
    for (t in 2:n_timepoints) {
      eps[, t, ] <- residual_ar1 * eps[, t - 1, ] + innov_sd * eps[, t, ]
    }
  }

  # Analyte correlation via Cholesky (per subject-timepoint slice)
  for (t in seq_len(n_timepoints)) {
    eps[, t, ] <- eps[, t, ] %*% chol_cor
  }

  # Scale by per-analyte residual SD
  for (k in seq_len(n_analytes)) {
    eps[, , k] <- eps[, , k] * resid_sd[k]
  }

  # Heteroscedastic: inflate variance for non-reference groups
  if (heteroscedastic) {
    non_ref <- which(group_int > 1)
    eps[non_ref, , ] <- eps[non_ref, , ] * 1.5
  }

  ## ---- Build grid and compute latent values ---------------------------------
  grid <- expand.grid(
    subject_idx   = seq_len(n_subjects),
    timepoint_idx = seq_len(n_timepoints),
    analyte_idx   = seq_len(n_analytes),
    KEEP.OUT.ATTRS = FALSE
  )
  n_obs <- nrow(grid)
  si <- grid$subject_idx
  ti <- grid$timepoint_idx
  ki <- grid$analyte_idx
  gi <- group_int[si]

  # Fixed effects (vectorised)
  fe_group <- (gi - 1) * group_eff_full[ki]
  fe_time  <- (ti - 1) * time_eff_full[ki]
  fe_int   <- (gi - 1) * (ti - 1) * int_eff_full[ki]

  # Covariate contribution
  fe_cov <- rep(0, n_obs)
  for (nm in names(cov_values)) {
    fe_cov <- fe_cov + cov_effects[[nm]] * cov_values[[nm]][si]
  }

  # Random effects (index into matrices / arrays)
  re  <- re_int_mat[cbind(si, ki)]
  rs  <- if (!is.null(re_slp_mat)) re_slp_mat[cbind(si, ki)] * (ti - 1) else 0
  epsilon <- eps[cbind(si, ti, ki)]

  # Latent log-concentration
  true_latent <- unname(mu[ki]) + fe_group + fe_time + fe_int +
    fe_cov + re + rs + epsilon

  ## ---- Missingness ----------------------------------------------------------
  observed       <- true_latent
  missing_reason <- rep(NA_character_, n_obs)

  # Monotone dropout
  if (dropout_rate > 0 && n_timepoints > 1) {
    for (s in seq_len(n_subjects)) {
      dropped <- FALSE
      for (tp in 2:n_timepoints) {
        if (dropped) {
          rows <- si == s & ti == tp
          observed[rows] <- NA_real_
          missing_reason[rows] <- "dropout"
        } else {
          p_drop <- dropout_rate
          if (dropout_informative) {
            prev_vals <- true_latent[si == s & ti == (tp - 1)]
            p_drop <- stats::plogis(stats::qlogis(dropout_rate) -
                                      0.5 * mean(prev_vals))
          }
          if (stats::runif(1) < p_drop) {
            dropped <- TRUE
            rows <- si == s & ti == tp
            observed[rows] <- NA_real_
            missing_reason[rows] <- "dropout"
          }
        }
      }
    }
  }

  # Non-monotone visit missingness
  if (visit_missingness > 0) {
    for (s in seq_len(n_subjects)) {
      for (tp in seq_len(n_timepoints)) {
        rows <- si == s & ti == tp
        if (any(!is.na(observed[rows])) && stats::runif(1) < visit_missingness) {
          observed[rows] <- NA_real_
          missing_reason[rows] <- ifelse(
            is.na(missing_reason[rows]), "missed_visit", missing_reason[rows]
          )
        }
      }
    }
  }

  ## ---- Assemble output ------------------------------------------------------
  latent_data <- tibble::tibble(
    subject_id     = subject_ids[si],
    timepoint      = tp_labels[ti],
    cytokine       = analyte_names[ki],
    group          = subject_groups[si],
    value          = observed,
    missing_reason = missing_reason
  )
  # Append covariate columns
  for (nm in names(cov_values)) {
    latent_data[[nm]] <- cov_values[[nm]][si]
  }

  truth_data <- tibble::tibble(
    true_latent       = true_latent,
    signal_analyte    = signal_mask[ki],
    true_group_effect = fe_group,
    true_time_effect  = fe_time,
    true_interaction  = fe_int
  )

  ## ---- Realized summaries ---------------------------------------------------
  realized_icc <- stats::setNames(rep(NA_real_, n_analytes), analyte_names)
  for (k in seq_len(n_analytes)) {
    var_re  <- stats::var(re_int_mat[, k])
    var_eps <- stats::var(as.vector(eps[, , k]))
    realized_icc[k] <- var_re / (var_re + var_eps)
  }

  if (re_intercept_sd > 0 && n_analytes > 1) {
    realized_cor <- stats::cor(re_int_mat)
  } else {
    realized_cor <- diag(n_analytes)
  }
  dimnames(realized_cor) <- list(analyte_names, analyte_names)

  realized_dropout <- sum(missing_reason %in% "dropout") / n_obs

  params <- list(
    design              = design,
    n_subjects          = n_subjects,
    n_timepoints        = n_timepoints,
    n_analytes          = n_analytes,
    group_levels        = group_levels,
    group_allocation    = group_allocation,
    group_effects       = group_eff_full,
    time_effects        = time_eff_full,
    interaction_effects = int_eff_full,
    signal_analytes     = signal_mask,
    signal_analyte_names = analyte_names[signal_mask],
    effect_direction    = effect_direction,
    re_intercept_sd     = re_intercept_sd,
    re_slope_sd         = re_slope_sd,
    residual_sd         = resid_sd,
    residual_ar1        = residual_ar1,
    heteroscedastic     = heteroscedastic,
    analyte_names       = analyte_names,
    analyte_means       = mu,
    analyte_sds         = sigma,
    analyte_correlation = analyte_correlation,
    cor_matrix          = cor_mat,
    covariates          = covariates,
    covariate_group_cor = covariate_group_cor,
    cov_effects         = cov_effects,
    dropout_rate        = dropout_rate,
    dropout_informative = dropout_informative,
    visit_missingness   = visit_missingness,
    seed                = seed,
    realized_icc        = realized_icc,
    realized_correlation = realized_cor,
    realized_dropout    = realized_dropout
  )

  list(
    latent_data = latent_data,
    truth_data  = truth_data,
    params      = params
  )
}


# ---- Stage 2: Assay process layer -------------------------------------------

#' @noRd
apply_assay_process <- function(
    latent_output,
    lod_values = NULL,
    ulod_values = NULL,
    lod_quantile = NULL,
    n_plates = 1,
    plate_sd = 0,
    n_lots = 1,
    lot_sd = 0,
    bg_shape = 0,
    bg_rate = 1,
    n_replicates = 1,
    replicate_sd = 0,
    well_failure_rate = 0,
    cv_threshold = 0.25,
    covariate_missingness = 0) {

  d      <- latent_output$latent_data
  tr     <- latent_output$truth_data
  params <- latent_output$params

  n_obs         <- nrow(d)
  analyte_names <- params$analyte_names
  n_analytes    <- params$n_analytes
  subjects      <- unique(d$subject_id)
  n_subj        <- length(subjects)

  ## ---- LOD / ULOD setup ------------------------------------------------------
  if (!is.null(lod_values)) {
    if (!is.null(names(lod_values))) {
      base_lod <- lod_values[analyte_names]
    } else {
      base_lod <- stats::setNames(lod_values, analyte_names)
    }
  } else if (!is.null(lod_quantile)) {
    base_lod <- stats::setNames(vapply(analyte_names, function(nm) {
      rv <- exp(d$value[d$cytokine == nm])
      rv <- rv[!is.na(rv)]
      if (length(rv) > 0) stats::quantile(rv, lod_quantile, names = FALSE)
      else 1.0
    }, numeric(1)), analyte_names)
  } else {
    lib <- analyte_library
    base_lod <- stats::setNames(vapply(analyte_names, function(nm) {
      idx <- match(nm, lib$analyte)
      if (!is.na(idx)) lib$lod[idx] else 1.0
    }, numeric(1)), analyte_names)
  }

  base_ulod <- if (!is.null(ulod_values)) {
    if (!is.null(names(ulod_values))) ulod_values[analyte_names]
    else stats::setNames(ulod_values, analyte_names)
  } else {
    stats::setNames(rep(Inf, n_analytes), analyte_names)
  }

  ## ---- Assign subjects to plates and lots ------------------------------------
  plate_labels <- paste0("plate_", ((seq_len(n_subj) - 1L) %% n_plates) + 1L)
  lot_labels   <- paste0("lot_", ((seq_len(n_subj) - 1L) %% n_lots) + 1L)
  plate_map <- stats::setNames(plate_labels, as.character(subjects))
  lot_map   <- stats::setNames(lot_labels, as.character(subjects))

  d$plate_id <- factor(plate_map[as.character(d$subject_id)])
  d$lot_id   <- factor(lot_map[as.character(d$subject_id)])

  plate_ids <- unique(plate_labels)
  lot_ids   <- unique(lot_labels)

  ## ---- Generate plate effects ------------------------------------------------
  plate_eff_vec <- stats::setNames(
    if (plate_sd > 0) stats::rnorm(length(plate_ids), 0, plate_sd)
    else rep(0, length(plate_ids)),
    plate_ids
  )

  ## ---- Generate lot effects (per lot x per analyte) --------------------------
  lot_eff_mat <- matrix(0, nrow = length(lot_ids), ncol = n_analytes,
                        dimnames = list(lot_ids, analyte_names))
  if (lot_sd > 0) {
    lot_eff_mat[] <- stats::rnorm(length(lot_ids) * n_analytes, 0, lot_sd)
  }

  ## ---- Lot-specific LODs -----------------------------------------------------
  lot_lod_mat  <- matrix(NA_real_, nrow = length(lot_ids), ncol = n_analytes,
                         dimnames = list(lot_ids, analyte_names))
  lot_ulod_mat <- matrix(NA_real_, nrow = length(lot_ids), ncol = n_analytes,
                         dimnames = list(lot_ids, analyte_names))
  for (l in lot_ids) {
    lot_lod_mat[l, ]  <- base_lod  * exp(lot_eff_mat[l, ])
    lot_ulod_mat[l, ] <- base_ulod * exp(lot_eff_mat[l, ])
  }

  ## ---- Convert to raw scale and apply pre-replicate effects ------------------
  raw_value <- exp(d$value)            # NAs from missingness preserved
  not_na    <- !is.na(raw_value)

  # Pre-expansion row indices
  plate_idx <- match(as.character(d$plate_id), plate_ids)
  lot_idx   <- match(as.character(d$lot_id), lot_ids)
  cyt_idx   <- match(d$cytokine, analyte_names)

  # 1. Plate effects (multiplicative on raw scale)
  pe_per_row <- plate_eff_vec[plate_idx]
  raw_value[not_na] <- raw_value[not_na] * exp(pe_per_row[not_na])

  # 2. Lot effects (per-analyte, multiplicative)
  le_per_row <- lot_eff_mat[cbind(lot_idx, cyt_idx)]
  raw_value[not_na] <- raw_value[not_na] * exp(le_per_row[not_na])

  # 3. Background binding (additive, per analyte per plate)
  if (bg_shape > 0) {
    bg_mat <- matrix(stats::rgamma(length(plate_ids) * n_analytes,
                                   shape = bg_shape, rate = bg_rate),
                     nrow = length(plate_ids), ncol = n_analytes,
                     dimnames = list(plate_ids, analyte_names))
    bg_per_row <- bg_mat[cbind(plate_idx, cyt_idx)]
    raw_value[not_na] <- raw_value[not_na] + bg_per_row[not_na]
  }

  # Record truth (these apply to all rows, including missing)
  obs_plate_effect <- pe_per_row
  obs_lot_effect   <- le_per_row
  true_raw         <- exp(tr$true_latent)

  ## ---- Replicate expansion ---------------------------------------------------
  if (n_replicates > 1) {
    rep_idx <- rep(seq_len(n_obs), each = n_replicates)
    d                <- d[rep_idx, ]
    tr               <- tr[rep_idx, ]
    raw_value        <- raw_value[rep_idx]
    true_raw         <- true_raw[rep_idx]
    obs_plate_effect <- obs_plate_effect[rep_idx]
    obs_lot_effect   <- obs_lot_effect[rep_idx]
    not_na           <- !is.na(raw_value)
    d$replicate_id   <- rep(seq_len(n_replicates), times = n_obs)
    n_obs_expanded   <- nrow(d)
  } else {
    d$replicate_id <- 1L
    n_obs_expanded <- n_obs
  }

  # 4. Replicate noise (applied even with 1 replicate if replicate_sd > 0)
  if (replicate_sd > 0) {
    rep_noise <- stats::rnorm(sum(not_na), 0, replicate_sd)
    raw_value[not_na] <- raw_value[not_na] * exp(rep_noise)
  }

  ## ---- 5. Censoring ----------------------------------------------------------
  # Per-row LOD/ULOD (lot-specific, on expanded data)
  lot_idx_exp <- match(as.character(d$lot_id), lot_ids)
  cyt_idx_exp <- match(d$cytokine, analyte_names)
  obs_lod  <- lot_lod_mat[cbind(lot_idx_exp, cyt_idx_exp)]
  obs_ulod <- lot_ulod_mat[cbind(lot_idx_exp, cyt_idx_exp)]

  cens_lod_flag  <- not_na & raw_value < obs_lod
  cens_ulod_flag <- not_na & is.finite(obs_ulod) & raw_value > obs_ulod

  raw_value[cens_lod_flag]  <- obs_lod[cens_lod_flag]
  raw_value[cens_ulod_flag] <- obs_ulod[cens_ulod_flag]

  ## ---- 6. QC flags -----------------------------------------------------------
  qc_flag <- rep("pass", n_obs_expanded)

  # Well failures
  if (well_failure_rate > 0) {
    n_fail <- round(n_obs_expanded * well_failure_rate)
    if (n_fail > 0) {
      fail_idx <- sample.int(n_obs_expanded, n_fail)
      raw_value[fail_idx] <- NA_real_
      qc_flag[fail_idx] <- "failed_well"
      d$missing_reason[fail_idx] <- "failed_well"
    }
  }

  # High-CV flag (replicates only)
  if (n_replicates > 1) {
    orig_id <- rep(seq_len(n_obs), each = n_replicates)
    cv_vals <- tapply(raw_value, orig_id, function(x) {
      xv <- x[!is.na(x)]
      if (length(xv) >= 2) stats::sd(xv) / mean(xv) else NA_real_
    })
    high_cv_obs <- !is.na(cv_vals) & cv_vals > cv_threshold
    row_high_cv <- high_cv_obs[orig_id]
    qc_flag[row_high_cv & qc_flag == "pass"] <- "high_cv"
  }

  ## ---- Convert back to log scale ---------------------------------------------
  log_value <- log(raw_value)

  ## ---- 7. Covariate missingness ----------------------------------------------
  if (covariate_missingness > 0 && !is.null(params$covariates)) {
    cov_cols <- intersect(names(d), names(params$covariates))
    for (col in cov_cols) {
      n_miss <- round(n_obs_expanded * covariate_missingness)
      if (n_miss > 0) {
        miss_idx <- sample.int(n_obs_expanded, n_miss)
        d[[col]][miss_idx] <- NA
      }
    }
  }

  ## ---- Assemble output -------------------------------------------------------
  d$value     <- log_value
  d$value_raw <- raw_value
  d$cens_lod  <- cens_lod_flag
  d$cens_ulod <- cens_ulod_flag
  d$lod       <- obs_lod
  d$ulod      <- obs_ulod
  d$qc_flag   <- qc_flag

  tr$true_raw      <- true_raw
  tr$plate_effect  <- obs_plate_effect
  tr$lot_effect    <- obs_lot_effect

  ## ---- Update params ---------------------------------------------------------
  realized_censoring <- stats::setNames(vapply(analyte_names, function(nm) {
    mask <- d$cytokine == nm & !is.na(log_value)
    if (sum(mask) > 0) mean(cens_lod_flag[mask]) else 0
  }, numeric(1)), analyte_names)

  params$assay <- list(
    lod_values           = base_lod,
    ulod_values          = base_ulod,
    lod_quantile         = lod_quantile,
    n_plates             = n_plates,
    plate_sd             = plate_sd,
    n_lots               = n_lots,
    lot_sd               = lot_sd,
    bg_shape             = bg_shape,
    bg_rate              = bg_rate,
    n_replicates         = n_replicates,
    replicate_sd         = replicate_sd,
    well_failure_rate    = well_failure_rate,
    cv_threshold         = cv_threshold,
    covariate_missingness = covariate_missingness,
    plate_effects        = plate_eff_vec,
    lot_effects          = lot_eff_mat,
    lot_lod              = lot_lod_mat,
    lot_ulod             = lot_ulod_mat
  )
  params$realized_censoring <- realized_censoring

  list(
    latent_data = d,
    truth_data  = tr,
    params      = params
  )
}


# ---- Stage 3: Format for immunoPlex -----------------------------------------

#' @noRd
format_for_immunoplex <- function(assay_output, format = "long") {
  d  <- assay_output$latent_data
  tr <- assay_output$truth_data
  p  <- assay_output$params

  # Ensure proper column types for downstream functions
  if (!is.factor(d$subject_id)) d$subject_id <- factor(d$subject_id)
  if (!is.factor(d$timepoint))  d$timepoint  <- factor(d$timepoint)
  if (!is.factor(d$group))      d$group      <- factor(d$group, levels = p$group_levels)
  d$cytokine <- as.character(d$cytokine)
  if (!is.factor(d$plate_id))   d$plate_id   <- factor(d$plate_id)
  if (!is.factor(d$lot_id))     d$lot_id     <- factor(d$lot_id)

  # Coerce to tibble
  if (!inherits(d, "tbl_df"))  d  <- tibble::as_tibble(d)
  if (!inherits(tr, "tbl_df")) tr <- tibble::as_tibble(tr)

  # Build meta list

  meta <- list(
    design               = p$design,
    dgp_params           = p,
    n_subjects           = as.integer(p$n_subjects),
    n_analytes           = as.integer(p$n_analytes),
    realized_censoring   = p$realized_censoring,
    realized_icc         = p$realized_icc,
    realized_correlation = p$realized_correlation,
    realized_dropout     = p$realized_dropout,
    seed                 = p$seed,
    timestamp            = Sys.time()
  )

  new_immuno_sim(data = d, truth = tr, meta = meta)
}


# ---- Exported entry point: simulate_immunoassay() ----------------------------

#' Simulate immunoassay data
#'
#' Single entry point for generating simulated multiplex immunoassay data.
#' Calls \code{simulate_latent_biology()} to generate latent biological
#' concentrations, \code{apply_assay_process()} to apply measurement effects
#' (plates, lots, censoring, replicates, QC failures), and
#' \code{format_for_immunoplex()} to format the output for downstream
#' immunoPlex analysis functions.
#'
#' @param n_subjects Integer. Total number of subjects. Default 60.
#' @param n_timepoints Integer. Number of time points. Default 2.
#' @param n_analytes Integer. Number of analytes. Default 20.
#' @param design Character. One of \code{"pre_post"}, \code{"time_course"},
#'   \code{"cross_sectional"}, \code{"paired_exposure"}. Default
#'   \code{"pre_post"}.
#' @param group_levels Character vector. Group names. Default
#'   \code{c("control", "treatment")}.
#' @param group_allocation Numeric vector or NULL. Allocation ratio.
#'   Default NULL (balanced).
#' @param group_effects Numeric. Log-scale group effect per signal analyte.
#'   Default 0.5.
#' @param time_effects Numeric vector or NULL. Log-scale time effects.
#' @param interaction_effects Numeric vector or NULL. Log-scale interaction
#'   effects.
#' @param signal_analytes Integer, logical vector, or character vector. Which
#'   analytes carry true signal. Integer = count (random selection). Default 5.
#' @param effect_direction Character. \code{"up"}, \code{"down"}, or
#'   \code{"mixed"}. Default \code{"mixed"}.
#' @param re_intercept_sd Numeric. Between-subject SD (log). Default 0.5.
#' @param re_slope_sd Numeric. Random slope SD. Default 0 (no random slopes).
#' @param residual_sd Numeric or named numeric. Within-subject SD (log).
#'   Default 1.0.
#' @param residual_ar1 Numeric. AR(1) correlation for within-subject residuals.
#'   Default 0.
#' @param heteroscedastic Logical. Allow variance to differ by group? Default
#'   FALSE.
#' @param analyte_means Named numeric or NULL. Per-analyte baseline log-means.
#' @param analyte_sds Named numeric or NULL. Per-analyte log-SDs.
#' @param analyte_correlation Matrix, character (\code{"block"}/\code{"toeplitz"}),
#'   or NULL. Default NULL (independent).
#' @param block_rho Numeric. Within-block correlation when
#'   \code{analyte_correlation = "block"}. Default 0.5.
#' @param covariates List or NULL. Named list of covariate specifications,
#'   each a list with elements \code{type} (\code{"continuous"} or
#'   \code{"binary"}) and \code{effect} (numeric).
#' @param covariate_group_cor Numeric. Correlation between covariates and
#'   group assignment (for confounding). Default 0.
#' @param dropout_rate Numeric. Per-timepoint dropout probability (monotone).
#'   Default 0.
#' @param dropout_informative Logical. Informative dropout? Default FALSE.
#' @param visit_missingness Numeric. Per-visit random missingness probability
#'   (non-monotone). Default 0.
#' @param lod_values Named numeric or NULL. Per-analyte LOD in pg/mL.
#' @param ulod_values Named numeric or NULL. Per-analyte ULOD in pg/mL.
#' @param lod_quantile Numeric or NULL. Set LOD at this quantile of the
#'   marginal distribution.
#' @param n_plates Integer. Number of plates/batches. Default 1.
#' @param plate_sd Numeric. Plate effect SD (log-scale). Default 0.
#' @param n_lots Integer. Number of reagent lots. Default 1.
#' @param lot_sd Numeric. Lot effect SD (log-scale). Default 0.
#' @param bg_shape Numeric. Background binding Gamma shape (0 = none).
#'   Default 0.
#' @param bg_rate Numeric. Background binding Gamma rate. Default 1.
#' @param n_replicates Integer. Technical replicates per observation. Default 1.
#' @param replicate_sd Numeric. Replicate-to-replicate SD (log-scale).
#'   Default 0.
#' @param well_failure_rate Numeric. Fraction of wells that fail. Default 0.
#' @param cv_threshold Numeric. CV threshold for high-CV QC flag. Default 0.25.
#' @param covariate_missingness Numeric. Fraction of covariate values set to
#'   NA. Default 0.
#' @param format Character. Output format: \code{"long"} (default).
#' @param seed Integer or NULL. RNG seed for reproducibility.
#'
#' @return An \code{\link{new_immuno_sim}} S3 object with components
#'   \code{$data} (observed tibble), \code{$truth} (ground-truth tibble),
#'   and \code{$meta} (simulation parameters and realized summaries).
#'
#' @details
#' \strong{Artifact-induced bias cannot be fixed downstream.} Benchmark B0
#' shows that when plate, lot, and background-binding artifacts are all
#' present simultaneously (\code{plate_sd > 0}, \code{lot_sd > 0},
#' \code{bg_shape > 0}), every downstream statistical family -- Tobit,
#' Gaussian glmmTMB, LOD/2 substitution, and Gamma -- exhibits negative
#' bias of roughly -0.09 to -0.17 in the estimated group effect, even
#' when Tobit is superior on clean data. Adding plate fixed effects inside
#' the model does not recover the full effect because plate indicators alone
#' do not capture lot effects or background binding. Experimental design and
#' preprocessing (lot harmonization, background subtraction, plate
#' balance) must address these artifacts upstream; no statistical method in
#' this package can fully correct for them.
#'
#' @examples
#' \dontrun{
#' # Simple cross-sectional simulation
#' sim <- simulate_immunoassay(
#'   n_subjects = 40, n_timepoints = 1, n_analytes = 10,
#'   design = "cross_sectional", signal_analytes = 3, seed = 42
#' )
#' print(sim)
#'
#' # Pre-post design with plate effects and censoring
#' sim <- simulate_immunoassay(
#'   n_subjects = 60, design = "pre_post",
#'   n_plates = 3, plate_sd = 0.1, seed = 123
#' )
#' summary(sim)
#' }
#'
#' @export
simulate_immunoassay <- function(
    n_subjects = 60,
    n_timepoints = 2,
    n_analytes = 20,
    design = "pre_post",
    # --- Biology parameters ---
    group_levels = c("control", "treatment"),
    group_allocation = NULL,
    group_effects = 0.5,
    time_effects = NULL,
    interaction_effects = NULL,
    signal_analytes = 5,
    effect_direction = "mixed",
    re_intercept_sd = 0.5,
    re_slope_sd = 0,
    residual_sd = 1.0,
    residual_ar1 = 0,
    heteroscedastic = FALSE,
    analyte_means = NULL,
    analyte_sds = NULL,
    analyte_correlation = NULL,
    block_rho = 0.5,
    covariates = NULL,
    covariate_group_cor = 0,
    dropout_rate = 0,
    dropout_informative = FALSE,
    visit_missingness = 0,
    # --- Assay parameters ---
    lod_values = NULL,
    ulod_values = NULL,
    lod_quantile = NULL,
    n_plates = 1,
    plate_sd = 0,
    n_lots = 1,
    lot_sd = 0,
    bg_shape = 0,
    bg_rate = 1,
    n_replicates = 1,
    replicate_sd = 0,
    well_failure_rate = 0,
    cv_threshold = 0.25,
    covariate_missingness = 0,
    # --- Output ---
    format = "long",
    seed = NULL) {

  # Stage 1: Latent biology
  latent <- simulate_latent_biology(
    n_subjects          = n_subjects,
    n_timepoints        = n_timepoints,
    n_analytes          = n_analytes,
    design              = design,
    group_levels        = group_levels,
    group_allocation    = group_allocation,
    group_effects       = group_effects,
    time_effects        = time_effects,
    interaction_effects = interaction_effects,
    signal_analytes     = signal_analytes,
    effect_direction    = effect_direction,
    re_intercept_sd     = re_intercept_sd,
    re_slope_sd         = re_slope_sd,
    residual_sd         = residual_sd,
    residual_ar1        = residual_ar1,
    heteroscedastic     = heteroscedastic,
    analyte_means       = analyte_means,
    analyte_sds         = analyte_sds,
    analyte_correlation = analyte_correlation,
    block_rho           = block_rho,
    covariates          = covariates,
    covariate_group_cor = covariate_group_cor,
    dropout_rate        = dropout_rate,
    dropout_informative = dropout_informative,
    visit_missingness   = visit_missingness,
    seed                = seed
  )

  # Stage 2: Assay process
  assay <- apply_assay_process(
    latent,
    lod_values            = lod_values,
    ulod_values           = ulod_values,
    lod_quantile          = lod_quantile,
    n_plates              = n_plates,
    plate_sd              = plate_sd,
    n_lots                = n_lots,
    lot_sd                = lot_sd,
    bg_shape              = bg_shape,
    bg_rate               = bg_rate,
    n_replicates          = n_replicates,
    replicate_sd          = replicate_sd,
    well_failure_rate     = well_failure_rate,
    cv_threshold          = cv_threshold,
    covariate_missingness = covariate_missingness
  )

  # Stage 3: Format for immunoPlex
  format_for_immunoplex(assay, format = format)
}
