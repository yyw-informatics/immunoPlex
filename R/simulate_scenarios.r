# Scenario wrappers and specialized generators for simulate_immunoassay()
#
# Each function calls simulate_immunoassay() with biologically motivated
# defaults and returns an immuno_sim object.


# ---- Scenario Wrappers (exported) -------------------------------------------

#' Simulate a vaccine/treatment trial
#'
#' Pre-post design with responder/non-responder groups and upregulated signal
#' analytes. Uses built-in analyte library LODs for realistic heterogeneous
#' per-analyte censoring. Suitable for testing \code{fit_one()},
#' \code{mcnemar_detection()}, and PLS-DA workflows.
#'
#' @param n_subjects Integer. Total number of subjects. Default 80.
#' @param n_analytes Integer. Number of analytes. Default 20.
#' @param signal_analytes Integer, logical, or character vector. Number or
#'   identity of signal analytes. Default 5.
#' @param group_effects Numeric. Log-scale group effect for signal analytes.
#'   Default 0.6.
#' @param time_effects Numeric. Log-scale time effect for signal analytes.
#'   Default 0.3.
#' @param interaction_effects Numeric. Log-scale interaction effect. Default 0.4.
#' @param seed Integer or NULL. RNG seed.
#' @param ... Additional arguments passed to \code{\link{simulate_immunoassay}}.
#'
#' @return An \code{\link{new_immuno_sim}} object.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_vaccine_trial(n_subjects = 60, seed = 42)
#' print(sim)
#' }
#'
#' @export
simulate_vaccine_trial <- function(
    n_subjects = 80,
    n_analytes = 20,
    signal_analytes = 5,
    group_effects = 0.6,
    time_effects = 0.3,
    interaction_effects = 0.4,
    seed = NULL,
    ...) {

  simulate_immunoassay(
    n_subjects          = n_subjects,
    n_timepoints        = 2,
    n_analytes          = n_analytes,
    design              = "pre_post",
    group_levels        = c("responder", "non-responder"),
    group_effects       = group_effects,
    time_effects        = time_effects,
    interaction_effects = interaction_effects,
    signal_analytes     = signal_analytes,
    effect_direction    = "up",
    re_intercept_sd     = 0.5,
    residual_sd         = 1.0,
    seed                = seed,
    ...
  )
}


#' Simulate an infection time course
#'
#' Time-course design with 5 timepoints and mild/severe groups. Features
#' severity-by-time interaction, heteroscedastic residuals, random slopes,
#' AR(1) correlation within subjects, and per-timepoint dropout.
#'
#' @param n_subjects Integer. Total number of subjects. Default 60.
#' @param n_timepoints Integer. Number of time points. Default 5.
#' @param n_analytes Integer. Number of analytes. Default 20.
#' @param signal_analytes Integer, logical, or character vector. Default 5.
#' @param group_effects Numeric. Log-scale severity effect. Default 0.8.
#' @param time_effects Numeric. Log-scale time effect. Default 0.4.
#' @param interaction_effects Numeric. Log-scale severity x time interaction.
#'   Default 0.5.
#' @param re_slope_sd Numeric. Random slope SD over time. Default 0.3.
#' @param residual_ar1 Numeric. AR(1) correlation. Default 0.4.
#' @param dropout_rate Numeric. Per-timepoint dropout probability (monotone
#'   missingness). Default 0.08, producing ~30% incomplete cases over 5
#'   timepoints.
#' @param seed Integer or NULL. RNG seed.
#' @param ... Additional arguments passed to \code{\link{simulate_immunoassay}}.
#'
#' @return An \code{\link{new_immuno_sim}} object.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_infection_timecourse(n_subjects = 40, seed = 42)
#' print(sim)
#' }
#'
#' @export
simulate_infection_timecourse <- function(
    n_subjects = 60,
    n_timepoints = 5,
    n_analytes = 20,
    signal_analytes = 5,
    group_effects = 0.8,
    time_effects = 0.4,
    interaction_effects = 0.5,
    re_slope_sd = 0.3,
    residual_ar1 = 0.4,
    dropout_rate = 0.08,
    seed = NULL,
    ...) {

  simulate_immunoassay(
    n_subjects          = n_subjects,
    n_timepoints        = n_timepoints,
    n_analytes          = n_analytes,
    design              = "time_course",
    group_levels        = c("mild", "severe"),
    group_effects       = group_effects,
    time_effects        = time_effects,
    interaction_effects = interaction_effects,
    signal_analytes     = signal_analytes,
    effect_direction    = "mixed",
    heteroscedastic     = TRUE,
    re_intercept_sd     = 0.5,
    re_slope_sd         = re_slope_sd,
    residual_sd         = 1.0,
    residual_ar1        = residual_ar1,
    dropout_rate        = dropout_rate,
    seed                = seed,
    ...
  )
}


#' Simulate a severity comparison study
#'
#' Cross-sectional design with 4 severity groups and graded effects.
#' Features correlated analyte blocks and optional class imbalance.
#'
#' @param n_subjects Integer. Total number of subjects. Default 120.
#' @param n_analytes Integer. Number of analytes. Default 20.
#' @param signal_analytes Integer, logical, or character vector. Default 5.
#' @param group_effects Numeric vector of length 4. Graded log-scale effects
#'   per severity group. Default \code{c(0, 0.3, 0.6, 1.0)}.
#' @param block_rho Numeric. Within-block analyte correlation. Default 0.6.
#' @param class_imbalance Logical. If TRUE, uses imbalanced allocation
#'   (4:3:2:1). Default FALSE.
#' @param seed Integer or NULL. RNG seed.
#' @param ... Additional arguments passed to \code{\link{simulate_immunoassay}}.
#'
#' @return An \code{\link{new_immuno_sim}} object.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_severity_comparison(n_subjects = 80, seed = 42)
#' print(sim)
#' }
#'
#' @export
simulate_severity_comparison <- function(
    n_subjects = 120,
    n_analytes = 20,
    signal_analytes = 5,
    group_effects = c(0, 0.3, 0.6, 1.0),
    block_rho = 0.6,
    class_imbalance = FALSE,
    seed = NULL,
    ...) {

  alloc <- if (class_imbalance) c(4, 3, 2, 1) else NULL

  # For the graded effect design, we need to translate 4-group effects to

  # the simulate_immunoassay interface. The DGP uses (group_int - 1) * effect,
  # so group 1 (none) gets 0, group 2 gets 1*effect, group 3 gets 2*effect,
  # group 4 gets 3*effect. We set the per-analyte effect so that group 4
  # achieves the maximum effect (group_effects[4]).
  # The linear coding means effect_per_analyte = max_effect / (n_groups - 1).
  max_eff <- max(group_effects)
  per_analyte_eff <- if (max_eff > 0) max_eff / 3 else 0

  simulate_immunoassay(
    n_subjects          = n_subjects,
    n_timepoints        = 1,
    n_analytes          = n_analytes,
    design              = "cross_sectional",
    group_levels        = c("none", "mild", "moderate", "severe"),
    group_allocation    = alloc,
    group_effects       = per_analyte_eff,
    signal_analytes     = signal_analytes,
    effect_direction    = "up",
    re_intercept_sd     = 0.5,
    residual_sd         = 1.0,
    analyte_correlation = "block",
    block_rho           = block_rho,
    seed                = seed,
    ...
  )
}


#' Simulate an exposure cohort study
#'
#' Paired-exposure design with 3 groups (unexposed, low, high) and household
#' clustering. Uses built-in analyte library LODs for realistic heterogeneous
#' per-analyte censoring, with LODs scaled up to produce heavier censoring
#' in unexposed groups. Features confounding covariates (age, sex).
#'
#' @param n_subjects Integer. Total number of subjects. Default 90.
#' @param n_analytes Integer. Number of analytes. Default 20.
#' @param signal_analytes Integer, logical, or character vector. Default 5.
#' @param group_effects Numeric. Log-scale exposure effect. Default 0.7.
#' @param covariate_group_cor Numeric. Correlation between covariates and group,
#'   creating genuine confounding. Default 0.3.
#' @param lod_scale Numeric. Multiplier applied to built-in analyte library LODs
#'   to increase overall censoring. Default 2.5 (produces heavy censoring in
#'   unexposed group while preserving per-analyte heterogeneity).
#' @param visit_missingness Numeric. Per-visit random missingness probability
#'   from scheduling/logistics. Default 0.05 (~5% missing visits).
#' @param seed Integer or NULL. RNG seed.
#' @param ... Additional arguments passed to \code{\link{simulate_immunoassay}}.
#'
#' @return An \code{\link{new_immuno_sim}} object.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_exposure_cohort(n_subjects = 60, seed = 42)
#' print(sim)
#' }
#'
#' @export
simulate_exposure_cohort <- function(
    n_subjects = 90,
    n_analytes = 20,
    signal_analytes = 5,
    group_effects = 0.7,
    covariate_group_cor = 0.3,
    lod_scale = 2.5,
    visit_missingness = 0.05,
    seed = NULL,
    ...) {

  # Scale up library LODs to produce heavier censoring while preserving

  # per-analyte heterogeneity (unlike lod_quantile which gives uniform rates)
  lib <- analyte_library
  n_use <- min(n_analytes, nrow(lib))
  scaled_lods <- stats::setNames(
    lib$lod[seq_len(n_use)] * lod_scale,
    lib$analyte[seq_len(n_use)]
  )
  # Pad if n_analytes > library size
  if (n_analytes > nrow(lib)) {
    extra <- stats::setNames(
      rep(1.0 * lod_scale, n_analytes - nrow(lib)),
      paste0("Analyte_", seq(nrow(lib) + 1, n_analytes))
    )
    scaled_lods <- c(scaled_lods, extra)
  }

  simulate_immunoassay(
    n_subjects          = n_subjects,
    n_timepoints        = 2,
    n_analytes          = n_analytes,
    design              = "paired_exposure",
    group_levels        = c("unexposed", "low_exposure", "high_exposure"),
    group_effects       = group_effects,
    signal_analytes     = signal_analytes,
    effect_direction    = "up",
    re_intercept_sd     = 0.6,
    residual_sd         = 1.0,
    covariates          = list(
      age = list(type = "continuous", effect = 0.2),
      sex = list(type = "binary", effect = 0.15)
    ),
    covariate_group_cor = covariate_group_cor,
    lod_values          = scaled_lods,
    visit_missingness   = visit_missingness,
    seed                = seed,
    ...
  )
}


# ---- Specialized Generators (exported) --------------------------------------

#' Simulate McNemar detection data
#'
#' Generates paired detection data parameterized directly by the
#' discordant-pair rate (the power-relevant quantity for McNemar's test).
#' Uses a latent threshold model so that analyte detections can be
#' correlated across analytes.
#'
#' @param n_subjects Integer. Number of paired subjects. Default 50.
#' @param n_analytes Integer. Number of analytes. Default 10.
#' @param baseline_detection Numeric in (0, 1). Baseline detection probability.
#'   Default 0.5.
#' @param discordant_rate Numeric in [0, 1). Target fraction of discordant
#'   pairs (detected at one timepoint but not the other). Default 0.15.
#' @param signal_analytes Integer, logical, or character vector. Which analytes
#'   have true detection changes (discordant pairs). Default 3.
#' @param analyte_rho Numeric in [0, 1). Correlation between analyte
#'   detection latents (0 = independent). Default 0.
#' @param seed Integer or NULL. RNG seed.
#'
#' @return An \code{\link{new_immuno_sim}} object with pre-post design. The
#'   \code{$data} includes \code{cens_lod} flags encoding detection status.
#'   The \code{$meta$mcnemar} element stores McNemar-specific parameters
#'   including realized discordant rates.
#'
#' @details
#' The generator works by constructing latent Gaussian variables for each
#' subject-analyte-timepoint combination, then thresholding to determine
#' detection. The latent model ensures:
#' \itemize{
#'   \item The marginal detection probability matches \code{baseline_detection}
#'   \item The fraction of discordant pairs (pre-detected / post-undetected or
#'     vice versa) matches \code{discordant_rate} for signal analytes
#'   \item Analyte detections can be correlated (shared latent factor) when
#'     \code{analyte_rho > 0}
#' }
#'
#' For non-signal analytes, detections are stable across timepoints
#' (discordant rate near 0).
#'
#' @examples
#' \dontrun{
#' sim <- simulate_mcnemar_data(
#'   n_subjects = 100, discordant_rate = 0.20, seed = 42
#' )
#' print(sim)
#' }
#'
#' @export
simulate_mcnemar_data <- function(
    n_subjects = 50,
    n_analytes = 10,
    baseline_detection = 0.5,
    discordant_rate = 0.15,
    signal_analytes = 3,
    analyte_rho = 0,
    seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  ## ---- Analyte setup ----------------------------------------------------------
  lib <- analyte_library
  if (n_analytes <= nrow(lib)) {
    analyte_names <- lib$analyte[seq_len(n_analytes)]
  } else {
    analyte_names <- paste0("Analyte_", seq_len(n_analytes))
  }

  ## ---- Signal analyte selection -----------------------------------------------
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

  ## ---- Latent threshold model -------------------------------------------------
  # Baseline threshold for desired marginal detection probability
  threshold <- stats::qnorm(1 - baseline_detection)

  # For discordant pairs: compute within-subject correlation across timepoints.
  # P(discordant) = P(X1 > t, X2 < t) + P(X1 < t, X2 > t)
  #               = 2 * P(X1 > t) * P(X2 < t | X1 > t)
  # With bivariate normal (rho_t = temporal correlation):
  #   P(disc) = 2 * [Phi(-t) - Phi2(-t, -t; rho_t)]
  # We need to solve for rho_t given discordant_rate.
  # For non-signal analytes, rho_t is high (stable detection).

  # Compute temporal correlation for signal analytes
  rho_t_signal <- .solve_temporal_rho(baseline_detection, discordant_rate)
  rho_t_null   <- 0.99  # very stable across timepoints for null analytes

  subject_ids <- factor(paste0("S", sprintf("%03d", seq_len(n_subjects))))
  n_tp <- 2L
  tp_labels <- factor(c("T1", "T2"), levels = c("T1", "T2"))

  # Build analyte correlation matrix for shared latent factor
  cor_analyte <- diag(n_analytes)
  if (analyte_rho > 0) {
    cor_analyte[lower.tri(cor_analyte)] <- analyte_rho
    cor_analyte[upper.tri(cor_analyte)] <- analyte_rho
  }
  chol_analyte <- chol(cor_analyte)

  # Generate latent values: n_subjects x 2 timepoints x n_analytes
  # For each analyte, (Z_pre, Z_post) have correlation rho_t
  latent_pre  <- matrix(stats::rnorm(n_subjects * n_analytes),
                         nrow = n_subjects, ncol = n_analytes)
  latent_pre  <- latent_pre %*% chol_analyte

  latent_post <- matrix(stats::rnorm(n_subjects * n_analytes),
                         nrow = n_subjects, ncol = n_analytes)
  latent_post <- latent_post %*% chol_analyte

  # Apply temporal correlation (mix pre and innovation for post)
  for (k in seq_len(n_analytes)) {
    rho_t <- if (signal_mask[k]) rho_t_signal else rho_t_null
    latent_post[, k] <- rho_t * latent_pre[, k] +
      sqrt(1 - rho_t^2) * latent_post[, k]
  }

  # Threshold to get detection
  detected_pre  <- latent_pre  > threshold
  detected_post <- latent_post > threshold

  # For non-signal analytes with rho_t_null = 0.99, a small residual discordance
  # rate is expected - that's fine and realistic.

  ## ---- Build data in long format -----------------------------------------------
  # Use the analyte library means for plausible concentration values
  lib_idx <- ((seq_len(n_analytes) - 1L) %% nrow(lib)) + 1L
  mu <- lib$log_mean[lib_idx]
  sigma <- lib$log_sd[lib_idx]
  base_lod <- lib$lod[lib_idx]

  grid <- expand.grid(
    subject_idx   = seq_len(n_subjects),
    timepoint_idx = seq_len(n_tp),
    analyte_idx   = seq_len(n_analytes),
    KEEP.OUT.ATTRS = FALSE
  )
  si <- grid$subject_idx
  ti <- grid$timepoint_idx
  ki <- grid$analyte_idx

  # True latent concentrations (log-scale)
  true_latent <- mu[ki] + stats::rnorm(nrow(grid), 0, sigma[ki] * 0.5)

  # Detection status based on the threshold model
  is_detected <- ifelse(ti == 1, detected_pre[cbind(si, ki)],
                         detected_post[cbind(si, ki)])

  # For undetected observations, set value at LOD (censored)
  value <- true_latent
  lod_per_row <- log(base_lod[ki])
  value[!is_detected] <- lod_per_row[!is_detected]

  data <- tibble::tibble(
    subject_id = subject_ids[si],
    timepoint  = tp_labels[ti],
    cytokine   = analyte_names[ki],
    group      = factor("all", levels = "all"),
    value      = value,
    value_raw  = exp(value),
    cens_lod   = !is_detected,
    cens_ulod  = FALSE,
    lod        = base_lod[ki],
    ulod       = Inf,
    plate_id   = factor("plate_1"),
    lot_id     = factor("lot_1"),
    replicate_id = 1L,
    qc_flag    = "pass",
    missing_reason = NA_character_
  )

  truth <- tibble::tibble(
    true_latent       = true_latent,
    true_raw          = exp(true_latent),
    signal_analyte    = signal_mask[ki],
    true_group_effect = 0,
    true_time_effect  = 0,
    true_interaction  = 0,
    plate_effect      = 0,
    lot_effect        = 0
  )

  ## ---- Realized discordant rates per analyte ----------------------------------
  realized_discordant <- stats::setNames(vapply(seq_len(n_analytes), function(k) {
    pre  <- detected_pre[, k]
    post <- detected_post[, k]
    mean(pre != post)
  }, numeric(1)), analyte_names)

  meta <- list(
    design               = "pre_post",
    dgp_params           = list(
      n_subjects         = n_subjects,
      n_analytes         = n_analytes,
      baseline_detection = baseline_detection,
      discordant_rate    = discordant_rate,
      signal_analytes    = signal_mask,
      signal_analyte_names = analyte_names[signal_mask],
      analyte_rho        = analyte_rho
    ),
    n_subjects           = as.integer(n_subjects),
    n_analytes           = as.integer(n_analytes),
    realized_censoring   = stats::setNames(vapply(seq_len(n_analytes), function(k) {
      mean(data$cens_lod[data$cytokine == analyte_names[k]])
    }, numeric(1)), analyte_names),
    realized_icc         = NULL,
    realized_correlation = NULL,
    realized_dropout     = 0,
    seed                 = seed,
    timestamp            = Sys.time(),
    mcnemar              = list(
      baseline_detection    = baseline_detection,
      target_discordant     = discordant_rate,
      realized_discordant   = realized_discordant,
      rho_t_signal          = rho_t_signal,
      rho_t_null            = rho_t_null,
      analyte_rho           = analyte_rho
    )
  )

  new_immuno_sim(data = data, truth = truth, meta = meta)
}


# Solve for temporal correlation rho_t that yields the target discordant rate
# given a marginal detection probability p.
# P(discordant) = 2 * [p - Phi2(-t, -t; rho_t)] where t = qnorm(1-p).
# Uses bisection search.
#' @noRd
.solve_temporal_rho <- function(p, disc_rate, tol = 1e-4) {
  if (disc_rate <= 0) return(0.99)

  threshold <- stats::qnorm(1 - p)

  # P(both detected) = Phi2(-t, -t; rho)
  # P(discordant) = 2 * (p - P(both detected))
  # For a given rho, compute discordant rate using mvtnorm if available,

  # otherwise use the approximation.

  disc_fn <- function(rho) {
    # P(X1 > t, X2 > t) with correlation rho
    p_both <- .bvnorm_upper(threshold, threshold, rho)
    2 * (p - p_both)
  }

  # Bisection: disc_fn is decreasing in rho (higher correlation = fewer discordant)
  lo <- -0.99
  hi <- 0.99
  for (i in 1:100) {
    mid <- (lo + hi) / 2
    d_mid <- disc_fn(mid)
    if (abs(d_mid - disc_rate) < tol) return(mid)
    if (d_mid > disc_rate) {
      lo <- mid  # need higher rho to reduce discordance
    } else {
      hi <- mid  # need lower rho to increase discordance
    }
  }
  mid
}


# Upper-tail probability of bivariate normal: P(X > a, Y > b) with
# correlation rho.
#' @noRd
.bvnorm_upper <- function(a, b, rho) {
  # P(X > a, Y > b) = P(X < -a, Y < -b) with same rho
  # Use the identity with standard bivariate normal CDF
  .bvnorm_cdf(-a, -b, rho)
}


# Bivariate normal CDF: P(X < a, Y < b) with correlation rho.
# Implements Drezner & Wesolowsky (1990) for moderate accuracy.
#' @noRd
.bvnorm_cdf <- function(a, b, rho) {
  if (abs(rho) < 1e-10) {
    return(stats::pnorm(a) * stats::pnorm(b))
  }

  # Use Gauss-Legendre quadrature via the formula:
  # Phi2(a, b; rho) = Phi(a)*Phi(b) + integral
  # For simplicity and reliability, use the tetrachoric series expansion

  # Simple numerical integration via transformation
  # P(X<a, Y<b; rho) using conditional:
  # = integral from -inf to a of Phi((b - rho*x)/sqrt(1-rho^2)) * phi(x) dx
  # Use Gauss-Hermite quadrature with 20 points

  sq <- sqrt(1 - rho^2)
  # 20-point Gauss-Hermite nodes and weights
  gh <- .gauss_hermite_20()
  nodes <- gh$nodes
  weights <- gh$weights

  # Transform: x = sqrt(2) * node, phi(x) dx -> exp(-node^2)/sqrt(pi) * sqrt(2) dnode
  # So integral = sum w_i * Phi((b - rho*sqrt(2)*node_i) / sq) where w_i already
  # include the 1/sqrt(pi) factor

  integrand <- vapply(seq_along(nodes), function(i) {
    x <- sqrt(2) * nodes[i]
    weights[i] * stats::pnorm((b - rho * x) / sq)
  }, numeric(1))

  # But we need P(X < a, ...) not the full integral to infinity
  # Actually use: P(X<a, Y<b; rho) = integral_{-inf}^{a} phi(x) Phi((b-rho*x)/sq) dx
  # With change of variable x = sqrt(2)*u, phi(x)dx = exp(-u^2)/sqrt(pi) du
  # P = integral_{-inf}^{a/sqrt(2)} exp(-u^2)/sqrt(pi) * Phi((b-rho*sqrt(2)*u)/sq) du

  # Use a different approach: truncated GH quadrature is messy.
  # Instead use the exact formula for bivariate normal via pmvnorm if possible,
  # or fall back to a simple approximation.


  # Most reliable simple approach: use the decomposition
  # P(X<a, Y<b; rho) = P(Y<b|X) integrated over X<a
  # = E[Phi((b - rho*X)/sq) | X < a] * Phi(a)
  # Use conditional sampling approximation with many points

  # Actually, let's use the classic Drezner formula which is compact:
  if (a <= 0 && b <= 0 && rho <= 0) {
    ap <- a / sqrt(2 * (1 - rho^2))
    bp <- b / sqrt(2 * (1 - rho^2))
    A <- c(0.3253030, 0.4211071, 0.1334425, 0.006374323)
    B <- c(0.1337764, 0.6243247, 1.3425378, 2.2626645)
    s <- 0
    for (i in 1:4) {
      for (j in 1:4) {
        s <- s + A[i] * A[j] *
          exp(ap * (2 * B[i] - ap) + bp * (2 * B[j] - bp) +
                2 * rho * (B[i] - ap) * (B[j] - bp))
      }
    }
    return(sqrt(1 - rho^2) / pi * s)
  }

  # Reduce other cases to the base case
  if (a * b * rho <= 0) {
    if (a <= 0 && b >= 0 && rho >= 0) {
      return(stats::pnorm(a) - .bvnorm_cdf(a, -b, -rho))
    } else if (a >= 0 && b <= 0 && rho >= 0) {
      return(stats::pnorm(b) - .bvnorm_cdf(-a, b, -rho))
    } else if (a >= 0 && b >= 0 && rho <= 0) {
      return(max(0, stats::pnorm(a) + stats::pnorm(b) - 1 +
                   .bvnorm_cdf(-a, -b, rho)))
    }
  }

  # a*b*rho > 0 case
  if (a * b * rho > 0) {
    sgn_a <- sign(a)
    sgn_b <- sign(b)
    rho_ab <- (rho * sgn_a * sgn_b - sqrt(a^2 - 2 * rho * a * b + b^2) *
                 sqrt((1 - rho^2) * 0 + 0)) # placeholder

    # Use the formula from Drezner (1978) for the general case
    delta <- sqrt(a^2 - 2 * rho * a * b + b^2)
    rho1 <- (rho * a - b) * sign(a) / delta
    rho2 <- (rho * b - a) * sign(b) / delta

    bvn <- .bvnorm_cdf(a, 0, rho1) + .bvnorm_cdf(b, 0, rho2) -
      (if (sign(a) * sign(b) < 0) max(0, stats::pnorm(a) + stats::pnorm(b) - 1) else 0)
    # Clamp
    return(max(0, min(1, bvn)))
  }

  # Fallback: independence approximation (should not reach here normally)
  stats::pnorm(a) * stats::pnorm(b)
}


#' @noRd
.gauss_hermite_20 <- function() {
  # 20-point Gauss-Hermite quadrature nodes and weights
  # (probabilist's convention: weight function exp(-x^2))
  nodes <- c(
    -5.387480890011233, -4.603682449550744, -3.944764040115625,
    -3.347854567383216, -2.788806058428131, -2.254974002089276,
    -1.738537712116586, -1.234076215395323, -0.737473728545394,
    -0.245340708300901,
     0.245340708300901,  0.737473728545394,  1.234076215395323,
     1.738537712116586,  2.254974002089276,  2.788806058428131,
     3.347854567383216,  3.944764040115625,  4.603682449550744,
     5.387480890011233
  )
  weights <- c(
    2.229393645534e-13, 4.399340992273e-10, 1.086069370769e-07,
    7.802556478532e-06, 2.283386360164e-04, 3.243773342238e-03,
    2.481052088746e-02, 1.090172060200e-01, 2.866755053628e-01,
    4.622436696006e-01,
    4.622436696006e-01, 2.866755053628e-01, 1.090172060200e-01,
    2.481052088746e-02, 3.243773342238e-03, 2.283386360164e-04,
    7.802556478532e-06, 1.086069370769e-07, 4.399340992273e-10,
    2.229393645534e-13
  )
  # Normalize weights to sum to sqrt(pi)
  weights <- weights / sum(weights) * sqrt(pi)
  list(nodes = nodes, weights = weights)
}


#' Simulate PLS-DA classification data
#'
#' Generates data with explicit discriminatory vs. noise analytes for PLS-DA
#' validation. Supports block correlation, class imbalance, and control over
#' signal strength.
#'
#' @param n_subjects Integer. Total subjects. Default 60.
#' @param n_analytes Integer. Total analytes (discriminatory + noise). Default 20.
#' @param n_discriminatory Integer. Number of truly discriminatory analytes.
#'   Default 5.
#' @param group_effects Numeric. Log-scale effect size for discriminatory
#'   analytes. Default 0.7.
#' @param analyte_correlation Character or NULL. \code{"block"} for block
#'   correlation, NULL for independent. Default \code{"block"}.
#' @param block_rho Numeric. Within-block correlation. Default 0.5.
#' @param class_imbalance Numeric or NULL. If numeric, specifies the ratio of
#'   group 1 to group 2 (e.g., 2 for 2:1). NULL for balanced. Default NULL.
#' @param lod_quantile Numeric or NULL. LOD quantile for censoring. Default NULL
#'   (use library LOD values).
#' @param seed Integer or NULL. RNG seed.
#' @param ... Additional arguments passed to \code{\link{simulate_immunoassay}}.
#'
#' @return An \code{\link{new_immuno_sim}} object.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_plsda_data(n_subjects = 50, n_discriminatory = 5, seed = 42)
#' print(sim)
#' }
#'
#' @export
simulate_plsda_data <- function(
    n_subjects = 60,
    n_analytes = 20,
    n_discriminatory = 5,
    group_effects = 0.7,
    analyte_correlation = "block",
    block_rho = 0.5,
    class_imbalance = NULL,
    lod_quantile = NULL,
    seed = NULL,
    ...) {

  alloc <- if (!is.null(class_imbalance)) c(class_imbalance, 1) else NULL

  simulate_immunoassay(
    n_subjects          = n_subjects,
    n_timepoints        = 1,
    n_analytes          = n_analytes,
    design              = "cross_sectional",
    group_levels        = c("case", "control"),
    group_allocation    = alloc,
    group_effects       = group_effects,
    signal_analytes     = n_discriminatory,
    effect_direction    = "mixed",
    re_intercept_sd     = 0.5,
    residual_sd         = 1.0,
    analyte_correlation = analyte_correlation,
    block_rho           = block_rho,
    lod_quantile        = lod_quantile,
    seed                = seed,
    ...
  )
}


#' Simulate ANCOVA data with confounding
#'
#' Generates group comparison data with covariates correlated with group
#' assignment, enabling tests of confounding adjustment. Supports nonlinear
#' covariate effects, heteroscedastic residuals, and outlier injection.
#'
#' @param n_subjects Integer. Total subjects. Default 60.
#' @param n_analytes Integer. Number of analytes. Default 10.
#' @param signal_analytes Integer, logical, or character vector. Default 3.
#' @param group_effects Numeric. Log-scale group effect. Default 0.5.
#' @param covariate_group_cor Numeric. Correlation between covariates and group,
#'   creating genuine confounding. Default 0.4.
#' @param heteroscedastic Logical. If TRUE, non-reference groups have inflated
#'   variance. Default FALSE.
#' @param outlier_rate Numeric in [0, 1). Fraction of observations to replace
#'   with outliers (shifted by 3 SD). Default 0.
#' @param nonlinear_cov Logical. If TRUE, adds a quadratic term to the
#'   covariate effect. Default FALSE.
#' @param seed Integer or NULL. RNG seed.
#' @param ... Additional arguments passed to \code{\link{simulate_immunoassay}}.
#'
#' @return An \code{\link{new_immuno_sim}} object. The \code{$meta$ancova}
#'   element stores ANCOVA-specific parameters including outlier indices
#'   and nonlinear covariate specification.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_ancova_data(
#'   n_subjects = 80, covariate_group_cor = 0.4,
#'   outlier_rate = 0.05, seed = 42
#' )
#' print(sim)
#' }
#'
#' @export
simulate_ancova_data <- function(
    n_subjects = 60,
    n_analytes = 10,
    signal_analytes = 3,
    group_effects = 0.5,
    covariate_group_cor = 0.4,
    heteroscedastic = FALSE,
    outlier_rate = 0,
    nonlinear_cov = FALSE,
    seed = NULL,
    ...) {

  # Build covariate spec with nonlinear option
  cov_spec <- list(
    age = list(type = "continuous", effect = 0.3),
    bmi = list(type = "continuous", effect = 0.2)
  )

  sim <- simulate_immunoassay(
    n_subjects          = n_subjects,
    n_timepoints        = 1,
    n_analytes          = n_analytes,
    design              = "cross_sectional",
    group_levels        = c("control", "treatment"),
    group_effects       = group_effects,
    signal_analytes     = signal_analytes,
    effect_direction    = "up",
    re_intercept_sd     = 0.5,
    residual_sd         = 1.0,
    heteroscedastic     = heteroscedastic,
    covariates          = cov_spec,
    covariate_group_cor = covariate_group_cor,
    seed                = seed,
    ...
  )

  # Post-processing: add nonlinear covariate effects and outliers
  # These must be applied after simulate_immunoassay since they aren't
  # native DGP parameters.

  if (!is.null(seed)) set.seed(seed + 1000L)

  outlier_idx <- integer(0)

  # Inject outliers
  if (outlier_rate > 0) {
    n_obs <- nrow(sim$data)
    n_outlier <- round(n_obs * outlier_rate)
    if (n_outlier > 0) {
      outlier_idx <- sort(sample.int(n_obs, n_outlier))
      # Shift outliers by +/- 3 SD (randomly)
      shift_sign <- sample(c(-1, 1), n_outlier, replace = TRUE)
      outlier_shift <- shift_sign * 3 * 1.0  # 3 * residual_sd
      sim$data$value[outlier_idx] <- sim$data$value[outlier_idx] + outlier_shift
      sim$data$value_raw[outlier_idx] <- exp(sim$data$value[outlier_idx])
    }
  }

  # Add nonlinear covariate effect (quadratic)
  if (nonlinear_cov && "age" %in% names(sim$data)) {
    quad_effect <- 0.1 * (sim$data$age^2 - 1)  # centered quadratic
    sim$data$value <- sim$data$value + quad_effect
    sim$data$value_raw <- exp(sim$data$value)
  }

  # Store ANCOVA-specific metadata
  sim$meta$ancova <- list(
    covariate_group_cor = covariate_group_cor,
    heteroscedastic     = heteroscedastic,
    outlier_rate        = outlier_rate,
    outlier_idx         = outlier_idx,
    nonlinear_cov       = nonlinear_cov
  )

  sim
}
