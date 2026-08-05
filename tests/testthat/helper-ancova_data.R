# Test data generator for ANCOVA tests
# Creates synthetic data with controlled effects for ancova_one / ancova_fit validation

#' Generate synthetic ANCOVA test data
#'
#' @param n_per_group Number of subjects per group (default: 30)
#' @param n_groups Number of group levels: 2 or 3 (default: 2)
#' @param n_analytes Number of analytes for batch testing (default: 1)
#' @param effect_size Numeric effect on outcome for group 2 vs 1 (default: 0.5)
#' @param include_conditional Logical. Add binary conditional covariates (default: TRUE)
#' @param high_prevalence Numeric prevalence for the high-prevalence binary covariate (default: 0.15)
#' @param low_prevalence Numeric prevalence for the low-prevalence binary covariate (default: 0.01)
#' @param seed Random seed (default: 42)
#'
#' @return A data frame suitable for ancova_one() / ancova_fit()
create_ancova_test_data <- function(n_per_group = 30,
                                     n_groups = 2,
                                     n_analytes = 1,
                                     effect_size = 0.5,
                                     include_conditional = TRUE,
                                     high_prevalence = 0.15,
                                     low_prevalence = 0.01,
                                     seed = 42) {
  set.seed(seed)

  group_labels <- if (n_groups == 2) {
    c("Control", "Treatment")
  } else {
    c("Control", "GroupA", "GroupB")
  }

  n_total <- n_per_group * n_groups
  analyte_names <- if (n_analytes == 1) "IL-6" else paste0("analyte_", seq_len(n_analytes))

  rows <- vector("list", n_analytes)
  for (ai in seq_len(n_analytes)) {
    # Generate per-analyte data with slightly different effects
    analyte_effect <- effect_size * (1 + (ai - 1) * 0.3)

    group_vec <- factor(rep(group_labels, each = n_per_group), levels = group_labels)

    # Continuous outcome and covariate (log-scale values)
    maternal_value <- rnorm(n_total, mean = 5, sd = 1.5)
    cord_value <- 0.4 * maternal_value +
      rnorm(n_total, mean = 3, sd = 1) +
      analyte_effect * (as.numeric(group_vec) - 1)

    # Continuous covariates
    age <- round(rnorm(n_total, mean = 28, sd = 5))
    gest_delivery <- round(rnorm(n_total, mean = 38, sd = 2), 1)

    d <- data.frame(
      subject_id = paste0("S", seq_len(n_total)),
      cytokine = analyte_names[ai],
      cord_value = cord_value,
      maternal_value = maternal_value,
      group = group_vec,
      age = age,
      gest_delivery = gest_delivery,
      stringsAsFactors = FALSE
    )

    # Conditional binary covariates
    if (include_conditional) {
      d$HIV <- rbinom(n_total, 1, prob = high_prevalence)
      d$malaria <- rbinom(n_total, 1, prob = low_prevalence)
    }

    rows[[ai]] <- d
  }

  do.call(rbind, rows)
}
