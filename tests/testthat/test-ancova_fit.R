# Tests for ancova_fit() — batch ANCOVA + S3 methods

describe("ancova_fit()", {

  test_that("returns immuno_ancova_set class with results tibble", {
    d <- create_ancova_test_data(n_per_group = 25, n_groups = 2, n_analytes = 3)
    result <- ancova_fit(d, outcome = "cord_value", group = "group",
                          covariate = "maternal_value",
                          analyte_col = "cytokine", quiet = TRUE)

    expect_s3_class(result, "immuno_ancova_set")
    expect_true(is.data.frame(result$results))
    expect_equal(result$n_analytes, 3L)
    expect_equal(nrow(result$results), 3)
  })

  test_that("FDR q-value column present with BH correction", {
    d <- create_ancova_test_data(n_per_group = 25, n_groups = 2, n_analytes = 4)
    result <- ancova_fit(d, outcome = "cord_value", group = "group",
                          covariate = "maternal_value",
                          analyte_col = "cytokine", quiet = TRUE)

    expect_true("q_value" %in% names(result$results))
    # q-values should be >= raw p-values (BH adjustment)
    valid <- !is.na(result$results$group_p) & !is.na(result$results$q_value)
    expect_true(all(result$results$q_value[valid] >= result$results$group_p[valid]))
  })

  test_that("per-analyte rows in results; named models list", {
    d <- create_ancova_test_data(n_per_group = 25, n_groups = 2, n_analytes = 3)
    result <- ancova_fit(d, outcome = "cord_value", group = "group",
                          covariate = "maternal_value",
                          analyte_col = "cytokine", quiet = TRUE)

    analytes <- unique(d$cytokine)
    expect_equal(sort(result$results$analyte), sort(analytes))
    expect_equal(sort(names(result$models)), sort(analytes))
    for (m in result$models) {
      expect_s3_class(m, "immuno_ancova")
    }
  })

  test_that("k-level model produces q_omnibus and pairwise q-values", {
    d <- create_ancova_test_data(n_per_group = 20, n_groups = 3, n_analytes = 3)
    result <- ancova_fit(d, outcome = "cord_value", group = "group",
                          covariate = "maternal_value",
                          analyte_col = "cytokine", quiet = TRUE)

    expect_true("q_omnibus" %in% names(result$results))
    # Should have pairwise contrast q-value columns
    q_cols <- grep("^q_", names(result$results), value = TRUE)
    expect_true(length(q_cols) >= 1)
  })

  test_that("quiet = TRUE suppresses messages", {
    d <- create_ancova_test_data(n_per_group = 20, n_groups = 2, n_analytes = 2)
    expect_message(
      ancova_fit(d, outcome = "cord_value", group = "group",
                  covariate = "maternal_value",
                  analyte_col = "cytokine", quiet = TRUE),
      NA  # no messages expected
    )
  })

  test_that("quiet = FALSE produces messages", {
    d <- create_ancova_test_data(n_per_group = 20, n_groups = 2, n_analytes = 2)
    expect_message(
      ancova_fit(d, outcome = "cord_value", group = "group",
                  covariate = "maternal_value",
                  analyte_col = "cytokine", quiet = FALSE),
      "ancova_fit"
    )
  })

  test_that("config integration: defaults from get_immunoplex_config()", {
    cfg <- get_immunoplex_config()
    expect_equal(cfg$ancova$iqr_multiplier, 3)
    expect_equal(cfg$ancova$min_n, 12L)
    expect_equal(cfg$ancova$fdr_method, "BH")
    expect_equal(cfg$ancova$fdr_threshold, 0.05)
    expect_equal(cfg$ancova$prevalence_threshold, 0.03)
  })

  test_that("log_validate = TRUE adds log columns to results", {
    d <- create_ancova_test_data(n_per_group = 25, n_groups = 2, n_analytes = 2)
    result <- ancova_fit(d, outcome = "cord_value", group = "group",
                          covariate = "maternal_value",
                          analyte_col = "cytokine",
                          log_validate = TRUE, quiet = TRUE)
    expect_true("log_effect" %in% names(result$results))
    expect_true("log_fold_change" %in% names(result$results))
    expect_true("log_p" %in% names(result$results))
    expect_true("q_log" %in% names(result$results))
  })

  test_that("errors on missing analyte column", {
    d <- create_ancova_test_data(n_per_group = 20, n_groups = 2)
    expect_error(
      ancova_fit(d, outcome = "cord_value", group = "group",
                  covariate = "maternal_value",
                  analyte_col = "nonexistent", quiet = TRUE),
      "not found"
    )
  })
})


describe("S3 methods execute without error", {

  test_that("print.immuno_ancova works for binary model", {
    d <- create_ancova_test_data(n_per_group = 25, n_groups = 2)
    result <- ancova_one(d, outcome = "cord_value", group = "group",
                          covariate = "maternal_value")
    result$analyte <- "IL-6"
    expect_output(print(result), "ANCOVA: IL-6")
  })

  test_that("print.immuno_ancova works for k-level model", {
    d <- create_ancova_test_data(n_per_group = 20, n_groups = 3)
    result <- ancova_one(d, outcome = "cord_value", group = "group",
                          covariate = "maternal_value")
    result$analyte <- "TNF-a"
    expect_output(print(result), "ANCOVA: TNF-a")
    expect_output(print(result), "contrasts")
  })

  test_that("print.immuno_ancova_set works", {
    d <- create_ancova_test_data(n_per_group = 20, n_groups = 2, n_analytes = 3)
    result <- ancova_fit(d, outcome = "cord_value", group = "group",
                          covariate = "maternal_value",
                          analyte_col = "cytokine", quiet = TRUE)
    expect_output(print(result), "Rank-based ANCOVA")
    expect_output(print(result), "3 analytes")
  })

  test_that("summary.immuno_ancova works", {
    d <- create_ancova_test_data(n_per_group = 25, n_groups = 2)
    result <- ancova_one(d, outcome = "cord_value", group = "group",
                          covariate = "maternal_value")
    expect_output(summary(result), "Rank-Based ANCOVA Summary")
  })

  test_that("summary.immuno_ancova_set works", {
    d <- create_ancova_test_data(n_per_group = 20, n_groups = 2, n_analytes = 2)
    result <- ancova_fit(d, outcome = "cord_value", group = "group",
                          covariate = "maternal_value",
                          analyte_col = "cytokine", quiet = TRUE)
    expect_output(summary(result), "Batch Summary")
  })

  test_that("plot.immuno_ancova returns ggplot objects", {
    d <- create_ancova_test_data(n_per_group = 25, n_groups = 2)
    result <- ancova_one(d, outcome = "cord_value", group = "group",
                          covariate = "maternal_value")
    plots <- plot(result, save_pdf = FALSE)
    expect_true(!is.null(plots))
    expect_true(inherits(plots$rvf_plot, "gg"))
    expect_true(inherits(plots$qq_plot, "gg"))
  })

  test_that("plot.immuno_ancova_set returns ggplot object", {
    d <- create_ancova_test_data(n_per_group = 20, n_groups = 2, n_analytes = 3)
    result <- ancova_fit(d, outcome = "cord_value", group = "group",
                          covariate = "maternal_value",
                          analyte_col = "cytokine", quiet = TRUE)
    p <- plot(result, save_pdf = FALSE)
    expect_true(inherits(p, "gg"))
  })
})
