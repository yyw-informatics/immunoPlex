# Tests for ancova_one() — single-analyte rank-based ANCOVA engine

describe("ancova_one()", {

  test_that("returns immuno_ancova class", {
    d <- create_ancova_test_data(n_per_group = 30, n_groups = 2)
    result <- ancova_one(d, outcome = "cord_value", group = "group",
                          covariate = "maternal_value")
    expect_s3_class(result, "immuno_ancova")
    expect_s3_class(result, "immuno_model")
  })

  test_that("binary model populates group_effect, group_p, omega_sq_partial", {
    d <- create_ancova_test_data(n_per_group = 30, n_groups = 2)
    result <- ancova_one(d, outcome = "cord_value", group = "group",
                          covariate = "maternal_value")

    expect_true(!is.na(result$group_effect))
    expect_true(!is.na(result$group_p))
    expect_true(!is.na(result$omega_sq_partial))
    expect_equal(result$n_groups, 2L)
    expect_equal(result$group_levels, c("Control", "Treatment"))
    expect_true(length(result$group_ci) == 2)
    expect_true(!is.na(result$r_squared))
  })

  test_that("k-level model populates contrasts, omnibus_f, omnibus_p", {
    d <- create_ancova_test_data(n_per_group = 25, n_groups = 3)
    result <- ancova_one(d, outcome = "cord_value", group = "group",
                          covariate = "maternal_value")

    expect_equal(result$n_groups, 3L)
    expect_true(!is.null(result$contrasts))
    expect_true(is.data.frame(result$contrasts))
    # 3 choose 2 = 3 pairwise contrasts
    expect_equal(nrow(result$contrasts), 3)
    expect_true(!is.na(result$omnibus_f))
    expect_true(!is.na(result$omnibus_p))
    expect_true(all(c("term", "estimate", "se", "p", "ci_lower", "ci_upper") %in%
                     names(result$contrasts)))
  })

  test_that("outlier removal reduces n_obs", {
    d <- create_ancova_test_data(n_per_group = 30, n_groups = 2)
    # Inject extreme outliers
    d$cord_value[1] <- 1000
    d$cord_value[2] <- -500

    result <- ancova_one(d, outcome = "cord_value", group = "group",
                          covariate = "maternal_value",
                          outlier_removal = TRUE, iqr_multiplier = 3)
    expect_true(result$n_removed > 0)
    expect_true(result$n_obs < nrow(d))
  })

  test_that("outlier_removal = FALSE keeps all observations", {
    d <- create_ancova_test_data(n_per_group = 30, n_groups = 2)
    d$cord_value[1] <- 1000

    result <- ancova_one(d, outcome = "cord_value", group = "group",
                          covariate = "maternal_value",
                          outlier_removal = FALSE)
    expect_equal(result$n_removed, 0L)
    expect_equal(result$n_obs, nrow(d))
  })

  test_that("conditional covariates: low-prevalence excluded, high-prevalence included", {
    d <- create_ancova_test_data(n_per_group = 40, n_groups = 2,
                                  include_conditional = TRUE,
                                  high_prevalence = 0.20,
                                  low_prevalence = 0.00)
    result <- ancova_one(d, outcome = "cord_value", group = "group",
                          covariate = "maternal_value",
                          conditional_covariates = c("HIV", "malaria"),
                          prevalence_threshold = 0.03)

    # HIV should be included (20% prevalence > 3%)
    expect_true("HIV" %in% result$covariates_included)
    # malaria should be dropped (0% prevalence < 3%)
    expect_true("malaria" %in% result$covariates_dropped)
  })

  test_that("min_n warning for small data", {
    d <- create_ancova_test_data(n_per_group = 5, n_groups = 2)
    result <- ancova_one(d, outcome = "cord_value", group = "group",
                          covariate = "maternal_value",
                          min_n = 20)
    expect_true(any(grepl("min_n", result$warnings)))
  })

  test_that("log_validate = TRUE populates log fields", {
    d <- create_ancova_test_data(n_per_group = 30, n_groups = 2)
    result <- ancova_one(d, outcome = "cord_value", group = "group",
                          covariate = "maternal_value",
                          log_validate = TRUE)

    expect_true(!is.null(result$log_effect))
    expect_true(!is.null(result$log_fold_change))
    expect_true(!is.null(result$log_p))
    expect_true(!is.null(result$log_model))
  })

  test_that("log_validate = FALSE leaves log fields NULL", {
    d <- create_ancova_test_data(n_per_group = 30, n_groups = 2)
    result <- ancova_one(d, outcome = "cord_value", group = "group",
                          covariate = "maternal_value",
                          log_validate = FALSE)

    expect_null(result$log_effect)
    expect_null(result$log_model)
  })

  test_that("effect size labels are correct", {
    expect_equal(.ancova_effect_label(0.20), "Large")
    expect_equal(.ancova_effect_label(0.14), "Large")
    expect_equal(.ancova_effect_label(0.10), "Medium")
    expect_equal(.ancova_effect_label(0.06), "Medium")
    expect_equal(.ancova_effect_label(0.03), "Small")
    expect_equal(.ancova_effect_label(0.01), "Small")
    expect_equal(.ancova_effect_label(0.005), "Negligible")
    expect_equal(.ancova_effect_label(NA), "Unknown")
    # Negative values (from bias correction) should be treated as Negligible
    expect_equal(.ancova_effect_label(-0.01), "Negligible")
  })

  test_that("errors on missing columns", {
    d <- create_ancova_test_data(n_per_group = 30, n_groups = 2)
    expect_error(
      ancova_one(d, outcome = "nonexistent", group = "group",
                  covariate = "maternal_value"),
      "Missing required columns"
    )
  })

  test_that("errors on non-data.frame input", {
    expect_error(
      ancova_one("not a data frame", outcome = "x", group = "g", covariate = "c"),
      "must be a data frame"
    )
  })

  test_that("additional covariates are included in formula", {
    d <- create_ancova_test_data(n_per_group = 30, n_groups = 2)
    result <- ancova_one(d, outcome = "cord_value", group = "group",
                          covariate = "maternal_value",
                          covariates = c("age", "gest_delivery"))
    expect_true(all(c("age", "gest_delivery") %in% result$covariates_included))
    expect_true(grepl("age", result$formula_used))
    expect_true(grepl("gest_delivery", result$formula_used))
  })

  test_that("formula includes rank-transformed terms", {
    d <- create_ancova_test_data(n_per_group = 30, n_groups = 2)
    result <- ancova_one(d, outcome = "cord_value", group = "group",
                          covariate = "maternal_value")
    expect_true(grepl("outcome_rank", result$formula_used))
    expect_true(grepl("covariate_rank", result$formula_used))
  })

  test_that("handles NA values in outcome/group/covariate", {
    d <- create_ancova_test_data(n_per_group = 30, n_groups = 2)
    d$cord_value[1:3] <- NA
    d$group[4] <- NA
    result <- ancova_one(d, outcome = "cord_value", group = "group",
                          covariate = "maternal_value")
    expect_true(result$n_obs <= nrow(d) - 4)
    expect_s3_class(result, "immuno_ancova")
  })
})


# ---------------------------------------------------------------------------
# rep_col / plate_col: switch from lm() to glmmTMB mixed model with RE terms
# ---------------------------------------------------------------------------

# Inline replicate-level DGP. Each subject contributes n_rep × n_tp rows so
# that the nested (subject_id:rep_id) random effect is actually identifiable
# (multiple observations per (subject, rep) cell).
.make_ancova_repdat <- function(n_subj = 40, n_rep = 2, n_tp = 3,
                                n_plates = 4, effect = 0.5,
                                sigma_subj = 0.5, sigma_rep = 0.3,
                                sigma_plate = 0.25, sigma_err = 0.3,
                                seed = 2468) {
  set.seed(seed)
  nr <- n_subj * n_rep * n_tp
  subj <- rep(seq_len(n_subj), each = n_rep * n_tp)
  rep_lab <- rep(rep(paste0("R", seq_len(n_rep)), each = n_tp), times = n_subj)
  group <- factor(rep(rep(c("Control", "Treatment"), length.out = n_subj),
                      each = n_rep * n_tp))
  plate <- factor(rep(paste0("P", sample(seq_len(n_plates), n_subj, TRUE)),
                      each = n_rep * n_tp))

  subj_re <- rnorm(n_subj, 0, sigma_subj)[subj]
  rep_key <- paste0(subj, ":", rep_lab)
  rep_re <- setNames(rnorm(length(unique(rep_key)), 0, sigma_rep),
                     unique(rep_key))[rep_key]
  plate_re <- setNames(rnorm(length(levels(plate)), 0, sigma_plate),
                       levels(plate))[as.character(plate)]

  maternal_value <- rnorm(nr, 5, 1.5)
  cord_value <- 0.4 * maternal_value + rnorm(nr, 3, sigma_err) +
                ifelse(group == "Treatment", effect, 0) +
                subj_re + rep_re + plate_re

  data.frame(
    subject_id = paste0("S", subj),
    cord_value = cord_value,
    maternal_value = maternal_value,
    group = group,
    rep_id = rep_lab,
    plate = plate,
    stringsAsFactors = FALSE
  )
}

describe("ancova_one() rep_col / plate_col", {

  test_that("rep_col = NULL reproduces existing lm-based behavior exactly", {
    d <- .make_ancova_repdat()
    base <- ancova_one(d, outcome = "cord_value", group = "group",
                       covariate = "maternal_value",
                       outlier_removal = FALSE)
    new  <- ancova_one(d, outcome = "cord_value", group = "group",
                       covariate = "maternal_value",
                       outlier_removal = FALSE,
                       rep_col = NULL, plate_col = NULL)
    expect_s3_class(base$model, "lm")
    expect_s3_class(new$model,  "lm")
    expect_equal(base$group_effect, new$group_effect)
    expect_equal(base$group_p,       new$group_p)
    expect_equal(base$r_squared,     new$r_squared)
  })

  test_that("rep_col switches to glmmTMB and estimates nested (subject:rep) variance > 0", {
    skip_if_not_installed("glmmTMB")
    d <- .make_ancova_repdat()
    res <- ancova_one(d, outcome = "cord_value", group = "group",
                      covariate = "maternal_value",
                      outlier_removal = FALSE, rep_col = "rep_id")
    expect_s3_class(res$model, "glmmTMB")
    expect_true(grepl("subject_id:rep_id", res$formula_used))
    expect_true(!is.na(res$group_effect))
    expect_true(!is.na(res$group_p))

    vc <- glmmTMB::VarCorr(res$model)$cond
    expect_true("subject_id:rep_id" %in% names(vc))
    sigma_rep <- attr(vc[["subject_id:rep_id"]], "stddev")
    expect_true(is.numeric(sigma_rep) && sigma_rep > 0)
  })

  test_that("plate_col adds flat plate variance and preserves the group point estimate direction", {
    skip_if_not_installed("glmmTMB")
    d <- .make_ancova_repdat()

    base <- ancova_one(d, outcome = "cord_value", group = "group",
                       covariate = "maternal_value",
                       outlier_removal = FALSE)
    res  <- ancova_one(d, outcome = "cord_value", group = "group",
                       covariate = "maternal_value",
                       outlier_removal = FALSE, plate_col = "plate")

    expect_s3_class(res$model, "glmmTMB")
    expect_true(grepl("\\(1 \\| plate\\)", res$formula_used))

    vc <- glmmTMB::VarCorr(res$model)$cond
    expect_true("plate" %in% names(vc))

    # Group point estimate should stay on the same side as the lm() version —
    # plate RE soaks up plate-level noise, it should not invert the sign.
    expect_true(sign(base$group_effect) == sign(res$group_effect))
  })

  test_that("rep_col + plate_col both appended and both variance components present", {
    skip_if_not_installed("glmmTMB")
    d <- .make_ancova_repdat()
    res <- ancova_one(d, outcome = "cord_value", group = "group",
                      covariate = "maternal_value",
                      outlier_removal = FALSE,
                      rep_col = "rep_id", plate_col = "plate")
    expect_s3_class(res$model, "glmmTMB")
    vc <- glmmTMB::VarCorr(res$model)$cond
    expect_true(all(c("subject_id:rep_id", "plate") %in% names(vc)))
  })

  test_that("rep_col / plate_col errors on nonexistent column name", {
    d <- .make_ancova_repdat()
    expect_error(
      ancova_one(d, outcome = "cord_value", group = "group",
                 covariate = "maternal_value", rep_col = "no_such_col"),
      "`rep_col` not found in data"
    )
    expect_error(
      ancova_one(d, outcome = "cord_value", group = "group",
                 covariate = "maternal_value", plate_col = "no_such_col"),
      "`plate_col` not found in data"
    )
  })

  test_that("non-character rep_col / plate_col errors clearly", {
    d <- .make_ancova_repdat()
    expect_error(
      ancova_one(d, outcome = "cord_value", group = "group",
                 covariate = "maternal_value", rep_col = 1),
      "must be a single column name"
    )
    expect_error(
      ancova_one(d, outcome = "cord_value", group = "group",
                 covariate = "maternal_value", plate_col = c("a", "b")),
      "must be a single column name"
    )
  })
})


describe("ancova_fit() rep_col / plate_col forwarding", {

  test_that("rep_col propagates into each per-analyte glmmTMB fit", {
    skip_if_not_installed("glmmTMB")
    d1 <- .make_ancova_repdat(n_subj = 30)
    d2 <- .make_ancova_repdat(n_subj = 30, seed = 1357)
    d1$cytokine <- "IL-6"
    d2$cytokine <- "IL-8"
    d <- rbind(d1, d2)

    res <- ancova_fit(d, outcome = "cord_value", group = "group",
                      covariate = "maternal_value",
                      analyte_col = "cytokine",
                      rep_col = "rep_id",
                      outlier_removal = FALSE, quiet = TRUE)
    expect_s3_class(res, "immuno_ancova_set")
    expect_equal(res$n_analytes, 2L)
    for (m in res$models) {
      expect_s3_class(m$model, "glmmTMB")
      expect_true(grepl("subject_id:rep_id", m$formula_used))
    }
  })

  test_that("ancova_fit surfaces rep_col validation errors from ancova_one", {
    d <- .make_ancova_repdat(n_subj = 30)
    d$cytokine <- "IL-6"
    expect_error(
      ancova_fit(d, outcome = "cord_value", group = "group",
                 covariate = "maternal_value",
                 analyte_col = "cytokine", rep_col = "no_such_col",
                 quiet = TRUE),
      "`rep_col` not found in data"
    )
  })
})
