# Test data generator for Gaussian glmmTMB tests
# Creates repeated-measures data with known effect sizes for validation

#' Generate Gaussian test data for glmmTMB model validation
#'
#' @param n_subjects Number of subjects (default: 30)
#' @param n_timepoints Number of timepoints (default: 3)
#' @param seed Random seed (default: 42)
#' @return Data frame with subject_id, timepoint, value, age, disease,
#'   cens_lod, cens_ulod, lod, ulod columns
create_gaussian_test_data <- function(n_subjects = 30,
                                      n_timepoints = 3,
                                      seed = 42) {
  set.seed(seed)
  n_obs <- n_subjects * n_timepoints

  subject_ids <- rep(seq_len(n_subjects), each = n_timepoints)
  timepoints <- rep(paste0("T", seq_len(n_timepoints)), times = n_subjects)

  # Random intercepts per subject
  re_intercept <- rnorm(n_subjects, sd = 0.3)

  # Fixed effects: intercept + timepoint + disease + age
  disease <- rep(sample(c("Control", "Case"), n_subjects, replace = TRUE),
                 each = n_timepoints)
  age <- rep(round(rnorm(n_subjects, mean = 30, sd = 5)), each = n_timepoints)

  # Known effects on log scale
  intercept <- 5.0
  tp_effects <- seq(0, 0.5, length.out = n_timepoints)
  disease_effect <- ifelse(disease == "Case", 0.3, 0)
  age_effect <- (age - 30) * 0.02

  mu <- intercept +
    tp_effects[match(timepoints, paste0("T", seq_len(n_timepoints)))] +
    disease_effect +
    age_effect +
    re_intercept[subject_ids]

  value <- rnorm(n_obs, mean = mu, sd = 0.5)

  dat <- data.frame(
    subject_id = factor(subject_ids),
    timepoint = factor(timepoints),
    disease = factor(disease),
    age = age,
    value = value,
    cens_lod = FALSE,
    cens_ulod = FALSE,
    lod = 0.1,
    ulod = 100,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  dat
}

#' Generate small Gaussian test data (for edge cases)
#'
#' @param n_subjects Number of subjects (default: 10)
#' @param seed Random seed (default: 99)
create_small_gaussian_data <- function(n_subjects = 10, seed = 99) {
  create_gaussian_test_data(n_subjects = n_subjects, n_timepoints = 2,
                            seed = seed)
}
