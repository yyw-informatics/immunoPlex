library(testthat)
library(tibble)

source(test_path("..", "..", "inst", "simulations", "helpers_benchmark.R"))

describe("Benchmark summary merging", {

test_that("existing summaries without cell_id append cleanly", {
  summary_file <- tempfile(fileext = ".rds")
  saveRDS(
    tibble(
      method = "baseline",
      estimate = 0.12
    ),
    summary_file
  )

  new_summary <- tibble(
    cell_id = 20001L,
    method = "phase-b",
    estimate = 0.34
  )

  expect_message(
    merge_with_existing_summary(new_summary, summary_file, verbose = TRUE),
    "usable cell_id unavailable from existing summary"
  )

  merged <- readRDS(summary_file)
  expect_equal(nrow(merged), 2)
  expect_true("cell_id" %in% names(merged))
  expect_equal(merged$method, c("baseline", "phase-b"))
})

test_that("overlapping cell_ids are replaced during merge", {
  summary_file <- tempfile(fileext = ".rds")
  saveRDS(
    tibble(
      cell_id = c(1L, 1L, 2L),
      method = c("a", "b", "a"),
      estimate = c(0.1, 0.2, 0.3)
    ),
    summary_file
  )

  new_summary <- tibble(
    cell_id = c(1L, 1L, 3L),
    method = c("a", "b", "a"),
    estimate = c(1.1, 1.2, 1.3)
  )

  merge_with_existing_summary(new_summary, summary_file, verbose = FALSE)
  merged <- readRDS(summary_file)

  expect_equal(nrow(merged), 4)
  expect_equal(sum(merged$cell_id == 1L), 2)
  expect_equal(merged$estimate[merged$cell_id == 2L], 0.3)
  expect_true(all(c(1.1, 1.2, 1.3) %in% merged$estimate))
})

})
