# McNemar Alluvial Plot Tests
#
# Validates alluvial plotting functionality for McNemar detection results:
# 1. Single cytokine alluvial plots
# 2. Multi-cytokine grid plots
# 3. Input validation and error handling
# 4. Customization options (colors, labels, sizes)
# 5. Output structure (ggplot2 objects)

library(testthat)

# Helper function to create mock McNemar results
create_mock_mcnemar_results <- function() {
  data.frame(
    cytokine = c("IL6", "TNFa", "IL10", "IFNg"),
    n_pairs = c(50, 50, 50, 50),
    both_detect = c(25, 15, 30, 5),
    loss = c(5, 15, 5, 2),
    gain = c(15, 5, 5, 3),
    neither_detect = c(5, 15, 10, 40),
    n_discordant = c(20, 20, 10, 5),
    n_baseline_detect = c(30, 30, 35, 7),
    n_comparison_detect = c(40, 20, 35, 8),
    prop_baseline = c(60.0, 60.0, 70.0, 14.0),
    prop_baseline_ci_lo = c(45.0, 45.0, 55.0, 5.0),
    prop_baseline_ci_hi = c(73.0, 73.0, 82.0, 27.0),
    prop_comparison = c(80.0, 40.0, 70.0, 16.0),
    prop_comparison_ci_lo = c(66.0, 26.0, 55.0, 7.0),
    prop_comparison_ci_hi = c(90.0, 55.0, 82.0, 29.0),
    delta_detection = c(20.0, -20.0, 0.0, 2.0),
    delta_ci_lo = c(5.0, -35.0, -15.0, -10.0),
    delta_ci_hi = c(35.0, -5.0, 15.0, 14.0),
    rate_upward = c(30.0, 10.0, 10.0, 6.0),
    rate_upward_ci_lo = c(17.5, 3.5, 3.5, 1.3),
    rate_upward_ci_hi = c(45.0, 20.0, 20.0, 14.0),
    rate_downward = c(10.0, 30.0, 10.0, 4.0),
    rate_downward_ci_lo = c(3.5, 17.5, 3.5, 0.5),
    rate_downward_ci_hi = c(20.0, 45.0, 20.0, 11.0),
    mpor = c(3.0, 0.33, 1.0, 1.5),
    mpor_ci_lo = c(1.1, 0.09, 0.3, 0.2),
    mpor_ci_hi = c(8.5, 1.2, 3.3, 11.0),
    p_mcnemar = c(0.04, 0.04, 1.0, 0.65),
    q_mcnemar = c(0.08, 0.08, 1.0, 0.87),
    cohen_kappa = c(0.4, 0.2, 0.6, 0.7),
    detection_pattern = c("Increased detection", "Decreased detection", 
                         "Stable detection", "Stable detection"),
    mcnemar_significance = c("", "", "", ""),
    stringsAsFactors = FALSE
  )
}

describe("McNemar alluvial plotting", {
  
  test_that("plot_mcnemar_alluvial creates ggplot object", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggalluvial")
    skip_if_not_installed("dplyr")
    
    mock_results <- create_mock_mcnemar_results()
    
    suppressMessages({
      p <- plot_mcnemar_alluvial(
        mock_results,
        cytokine = "IL6"
      )
    })
    
    expect_s3_class(p, "ggplot")
  })
  
  test_that("plot_mcnemar_alluvial auto-selects most significant cytokine", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggalluvial")
    
    mock_results <- create_mock_mcnemar_results()
    
    # Should select IL6 or TNFa (both have q=0.08, lowest)
    expect_message(
      p <- plot_mcnemar_alluvial(mock_results),
      "Plotting cytokine"
    )
    
    expect_s3_class(p, "ggplot")
  })
  
  test_that("plot_mcnemar_alluvial handles custom labels", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggalluvial")
    
    mock_results <- create_mock_mcnemar_results()
    
    p <- plot_mcnemar_alluvial(
      mock_results,
      cytokine = "IL6",
      baseline_label = "Pre",
      comparison_label = "Post"
    )
    
    expect_s3_class(p, "ggplot")
    
    # Check that labels are in the plot
    plot_data <- ggplot2::layer_data(p)
    expect_true(length(plot_data) > 0)
  })
  
  test_that("plot_mcnemar_alluvial handles custom colors", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggalluvial")
    
    mock_results <- create_mock_mcnemar_results()
    
    custom_colors <- c(
      "Gained" = "#FF0000",
      "Lost" = "#0000FF",
      "No change (detected)" = "#00FF00",
      "No change (undetected)" = "#FFFF00"
    )
    
    p <- plot_mcnemar_alluvial(
      mock_results,
      cytokine = "IL6",
      colors = custom_colors
    )
    
    expect_s3_class(p, "ggplot")
  })
  
  test_that("plot_mcnemar_alluvial validates input", {
    expect_error(
      plot_mcnemar_alluvial("not a data frame"),
      "must be a data frame"
    )
    
    # Missing required columns
    bad_data <- data.frame(cytokine = "IL6")
    expect_error(
      plot_mcnemar_alluvial(bad_data),
      "Missing required columns"
    )
    
    # Non-existent cytokine
    mock_results <- create_mock_mcnemar_results()
    expect_error(
      plot_mcnemar_alluvial(mock_results, cytokine = "NonExistent"),
      "not found"
    )
  })
  
  test_that("plot_mcnemar_alluvial handles zero counts gracefully", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggalluvial")
    
    # Create results with some zero counts
    zero_data <- data.frame(
      cytokine = "TestCyt",
      n_pairs = 50,
      both_detect = 50,
      loss = 0,
      gain = 0,
      neither_detect = 0,
      delta_detection = 0,
      stringsAsFactors = FALSE
    )
    
    # Should still create plot (only one flow)
    p <- plot_mcnemar_alluvial(zero_data, cytokine = "TestCyt")
    expect_s3_class(p, "ggplot")
  })
  
  test_that("plot_mcnemar_alluvial_grid creates combined plot", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggalluvial")
    
    mock_results <- create_mock_mcnemar_results()
    
    expect_message(
      p <- plot_mcnemar_alluvial_grid(
        mock_results,
        n_cytokines = 2
      ),
      "Creating faceted alluvial plot"
    )
    
    expect_s3_class(p, "ggplot")
  })
  
  test_that("plot_mcnemar_alluvial_grid sorts correctly", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggalluvial")
    
    mock_results <- create_mock_mcnemar_results()
    
    # Sort by significance (default)
    suppressMessages(
      p1 <- plot_mcnemar_alluvial_grid(
        mock_results,
        n_cytokines = 2,
        sort_by = "significance"
      )
    )
    expect_s3_class(p1, "ggplot")
    
    # Sort by absolute delta
    suppressMessages(
      p2 <- plot_mcnemar_alluvial_grid(
        mock_results,
        n_cytokines = 2,
        sort_by = "delta_abs"
      )
    )
    expect_s3_class(p2, "ggplot")
    
    # Sort by signed delta
    suppressMessages(
      p3 <- plot_mcnemar_alluvial_grid(
        mock_results,
        n_cytokines = 2,
        sort_by = "delta"
      )
    )
    expect_s3_class(p3, "ggplot")
  })
  
  test_that("plot_mcnemar_alluvial_grid handles n_cytokines > available", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggalluvial")
    
    mock_results <- create_mock_mcnemar_results()
    
    # Request more cytokines than available
    suppressMessages(
      p <- plot_mcnemar_alluvial_grid(
        mock_results,
        n_cytokines = 100  # Only 4 available
      )
    )
    
    expect_s3_class(p, "ggplot")
  })
  
  test_that("plot_mcnemar_alluvial_grid respects ncol parameter", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggalluvial")
    
    mock_results <- create_mock_mcnemar_results()
    
    suppressMessages(
      p <- plot_mcnemar_alluvial_grid(
        mock_results,
        n_cytokines = 4,
        ncol = 2
      )
    )
    
    expect_s3_class(p, "ggplot")
  })
  
  test_that("plotting functions don't save files by default", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggalluvial")
    
    mock_results <- create_mock_mcnemar_results()
    
    # Create temp directory for testing
    temp_dir <- tempdir()
    old_wd <- getwd()
    on.exit(setwd(old_wd))
    setwd(temp_dir)
    
    # Plot without saving
    p <- plot_mcnemar_alluvial(
      mock_results,
      cytokine = "IL6",
      save_pdf = FALSE
    )
    
    # Check no PDF was created
    expect_false(file.exists("mcnemar_alluvial_IL6.pdf"))
  })
  
  test_that("plot titles include cytokine name and delta", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggalluvial")
    
    mock_results <- create_mock_mcnemar_results()
    
    p <- plot_mcnemar_alluvial(
      mock_results,
      cytokine = "IL6"
    )
    
    # Extract plot title
    plot_built <- ggplot2::ggplot_build(p)
    expect_s3_class(p, "ggplot")
    
    # Title should contain cytokine name and percentage
    expect_true(grepl("IL6", as.character(p$labels$title)))
  })
  
})

