# Tests for PLS-DA plotting functions

# Helper function to create a test model
create_test_model <- function() {
  test_obj <- list(
    data = data.frame(
      SubjectID = paste0("S", 1:30),
      group = rep(c("Control", "Treatment"), each = 15),
      site = rep(c("A", "B", "C"), length.out = 30),
      stringsAsFactors = FALSE
    ),
    lod = data.frame(
      cytokine = paste0("cyt", 1:10),
      lod = rep(10, 10)
    )
  )
  
  # Add cytokines with some group separation
  for (i in 1:10) {
    control_vals <- rnorm(15, mean = 100, sd = 20)
    treatment_vals <- rnorm(15, mean = 120, sd = 20)
    test_obj$data[[paste0("cyt", i)]] <- c(control_vals, treatment_vals)
  }
  
  preprocessed <- plsda_preprocess(
    data = test_obj$data,
    cytokine_cols = paste0("cyt", 1:10),
    metadata_cols = c("SubjectID", "group", "site"),
    lod_lookup = test_obj$lod,
    verbose = FALSE
  )
  
  model <- plsda_fit(
    preprocessed_data = preprocessed,
    response_var = "group",
    n_components = 2,
    permutations = 0,
    verbose = FALSE
  )
  
  return(model)
}


test_that("plot.plsda_model validates plot_type argument", {
  skip_if_not_installed("ropls")
  
  model <- create_test_model()
  
  # Valid plot types should work
  expect_no_error(plot(model, plot_type = "scores"))
  expect_no_error(plot(model, plot_type = "loadings"))
  expect_no_error(plot(model, plot_type = "vip"))
  expect_no_error(plot(model, plot_type = "biplot"))
  
  # Invalid plot type should error
  expect_error(
    plot(model, plot_type = "invalid"),
    "'arg' should be one of"
  )
})


test_that("plot_pls_scores creates valid ggplot", {
  skip_if_not_installed("ropls")
  
  model <- create_test_model()
  
  p <- plot_pls_scores(model, components = c(1, 2))
  
  # Should return a ggplot object
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Should have correct labels
  expect_true(grepl("Component 1", p$labels$x))
  expect_true(grepl("Component 2", p$labels$y))
})


test_that("plot_pls_scores handles 1D case", {
  skip_if_not_installed("ropls")
  
  model <- create_test_model()
  
  # 1D plot (only component 1)
  p <- plot_pls_scores(model, components = 1)
  
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})


test_that("plot_pls_scores respects color_by argument", {
  skip_if_not_installed("ropls")
  
  model <- create_test_model()
  
  # Should work with default (response variable)
  p1 <- plot_pls_scores(model)
  expect_s3_class(p1, "ggplot")
  
  # Should work with custom color variable
  p2 <- plot_pls_scores(model, color_by = "site")
  expect_s3_class(p2, "ggplot")
  
  # Should error with non-existent variable
  expect_error(
    plot_pls_scores(model, color_by = "nonexistent"),
    "not found in model scores"
  )
})


test_that("plot_pls_loadings creates valid ggplot", {
  skip_if_not_installed("ropls")
  
  model <- create_test_model()
  
  p <- plot_pls_loadings(model, components = c(1, 2))
  
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Should have facets for components
  expect_true(!is.null(p$facet))
})


test_that("plot_pls_loadings respects top_n argument", {
  skip_if_not_installed("ropls")
  
  model <- create_test_model()
  
  p <- plot_pls_loadings(model, top_n = 5)
  
  expect_s3_class(p, "ggplot")
  
  # Build plot to check data
  built <- ggplot2::ggplot_build(p)
  # Should have fewer points than total cytokines (5 cytokines × 2 components = 10 bars)
  expect_lte(nrow(built$data[[1]]), 10)
})


test_that("plot_pls_vip creates valid ggplot", {
  skip_if_not_installed("ropls")
  
  model <- create_test_model()
  
  p <- plot_pls_vip(model, vip_threshold = 1.0)
  
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Should have horizontal bars
  expect_true(!is.null(p$coordinates))
})


test_that("plot_pls_vip respects top_n argument", {
  skip_if_not_installed("ropls")
  
  model <- create_test_model()
  
  p <- plot_pls_vip(model, top_n = 5)
  
  expect_s3_class(p, "ggplot")
  
  # Build plot to check data
  built <- ggplot2::ggplot_build(p)
  # Should have at most 5 bars
  expect_lte(nrow(built$data[[1]]), 5)
})


test_that("plot_pls_biplot creates valid ggplot", {
  skip_if_not_installed("ropls")
  skip_if_not_installed("ggrepel")
  
  model <- create_test_model()
  
  p <- plot_pls_biplot(model, components = c(1, 2), top_loadings = 5)
  
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Should have multiple layers (points, arrows, labels)
  expect_gte(length(p$layers), 3)
})


test_that("plot_pls_biplot respects color_by and shape_by", {
  skip_if_not_installed("ropls")
  skip_if_not_installed("ggrepel")
  
  model <- create_test_model()
  
  # With color_by
  p1 <- plot_pls_biplot(model, color_by = "group")
  expect_s3_class(p1, "ggplot")
  
  # With both color_by and shape_by
  p2 <- plot_pls_biplot(model, color_by = "group", shape_by = "site")
  expect_s3_class(p2, "ggplot")
})


test_that("all plot types work through generic plot() method", {
  skip_if_not_installed("ropls")
  skip_if_not_installed("ggrepel")
  
  model <- create_test_model()
  
  # Test all plot types through generic
  p1 <- plot(model, plot_type = "scores")
  p2 <- plot(model, plot_type = "loadings")
  p3 <- plot(model, plot_type = "vip")
  p4 <- plot(model, plot_type = "biplot")
  
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
})
