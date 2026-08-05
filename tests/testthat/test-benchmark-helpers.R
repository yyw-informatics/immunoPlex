library(testthat)
library(tibble)

# helpers_benchmark.R ships under inst/simulations/, which R CMD check flattens
# to <pkg>/simulations/ in the installed tree. Resolve it through system.file()
# so the test works both from the source tree and from an installed package;
# fall back to the source layout for devtools::test() before installation.
.helpers_path <- system.file("simulations", "helpers_benchmark.R",
                             package = "immunoPlex")
if (!nzchar(.helpers_path)) {
  .helpers_path <- test_path("..", "..", "inst", "simulations",
                             "helpers_benchmark.R")
}
if (!file.exists(.helpers_path)) {
  skip("helpers_benchmark.R not available in this installation")
}
source(.helpers_path)

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
    "usable key .* unavailable from existing summary"
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
