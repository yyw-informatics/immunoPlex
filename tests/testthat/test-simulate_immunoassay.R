# Tests for immuno_sim S3 class and methods

# ---- Helpers: build minimal valid components ---------------------------------

make_data <- function(n_sub = 5, n_cyt = 3) {
  grid <- expand.grid(
    subject_id = factor(paste0("S", seq_len(n_sub))),
    cytokine   = paste0("CYT", seq_len(n_cyt)),
    stringsAsFactors = FALSE
  )
  grid$value      <- rnorm(nrow(grid), 2, 1)
  grid$group      <- factor(ifelse(as.integer(grid$subject_id) <= n_sub / 2,
                                   "ctrl", "trt"))
  grid$timepoint  <- factor("T1")
  grid$cens_lod   <- grid$value < 0.5
  grid$cens_ulod  <- rep(FALSE, nrow(grid))
  tibble::as_tibble(grid)
}

make_truth <- function(data) {
  tibble::tibble(
    true_latent      = data$value + rnorm(nrow(data), 0, 0.1),
    true_raw         = exp(data$value),
    signal_analyte   = data$cytokine == "CYT1",
    true_group_effect = ifelse(data$cytokine == "CYT1", 0.5, 0),
    true_time_effect  = 0,
    true_interaction  = 0,
    plate_effect      = 0,
    lot_effect        = 0
  )
}

make_meta <- function(n_sub = 5, n_cyt = 3) {
  list(
    design     = "pre_post",
    n_subjects = as.integer(n_sub),
    n_analytes = as.integer(n_cyt),
    seed       = 42L,
    timestamp  = Sys.time()
  )
}

make_sim <- function(n_sub = 5, n_cyt = 3) {
  d <- make_data(n_sub, n_cyt)
  new_immuno_sim(data = d, truth = make_truth(d), meta = make_meta(n_sub, n_cyt))
}


# ---- Construction tests ------------------------------------------------------

test_that("new_immuno_sim creates valid object", {
  sim <- make_sim()
  expect_s3_class(sim, "immuno_sim")
  expect_named(sim, c("data", "truth", "meta"))
  expect_s3_class(sim$data, "tbl_df")
  expect_s3_class(sim$truth, "tbl_df")
  expect_true(is.list(sim$meta) && !is.data.frame(sim$meta))
})

test_that("constructor coerces plain data.frame to tibble", {
  d <- as.data.frame(make_data())
  tr <- as.data.frame(make_truth(make_data()))
  sim <- new_immuno_sim(d, tr, make_meta())
  expect_s3_class(sim$data, "tbl_df")
  expect_s3_class(sim$truth, "tbl_df")
})

test_that("constructor accepts tibble inputs directly", {
  d <- make_data()
  tr <- make_truth(d)
  sim <- new_immuno_sim(d, tr, make_meta())
  expect_s3_class(sim, "immuno_sim")
})


# ---- Validation error tests --------------------------------------------------

test_that("data must be a data.frame", {
  expect_error(
    new_immuno_sim(data = list(a = 1), truth = make_truth(make_data()),
                   meta = make_meta()),
    "data.*must be a data.frame"
  )
})

test_that("data must contain required columns", {
  d <- make_data()
  d$value <- NULL
  expect_error(
    new_immuno_sim(d, make_truth(make_data()), make_meta()),
    "missing required columns.*value"
  )
})

test_that("truth must be a data.frame", {
  expect_error(
    new_immuno_sim(make_data(), truth = "not_a_df", meta = make_meta()),
    "truth.*must be a data.frame"
  )
})

test_that("truth must contain true_latent", {
  d <- make_data()
  tr <- make_truth(d)
  tr$true_latent <- NULL
  expect_error(
    new_immuno_sim(d, tr, make_meta()),
    "true_latent"
  )
})

test_that("data and truth must have same nrow", {
  d <- make_data(n_sub = 5)
  tr <- make_truth(make_data(n_sub = 10))
  expect_error(
    new_immuno_sim(d, tr, make_meta()),
    "same number of rows"
  )
})

test_that("meta must be a list", {
  expect_error(
    new_immuno_sim(make_data(), make_truth(make_data()), meta = "not_list"),
    "meta.*must be a named list"
  )
})

test_that("meta must not be a data.frame", {
  expect_error(
    new_immuno_sim(make_data(), make_truth(make_data()),
                   meta = data.frame(design = "x", n_subjects = 1)),
    "meta.*must be a named list"
  )
})

test_that("meta$design must be a character string", {
  m <- make_meta()
  m$design <- 123
  expect_error(
    new_immuno_sim(make_data(), make_truth(make_data()), m),
    "meta\\$design.*character"
  )
})

test_that("meta$n_subjects must be a positive number", {
  m <- make_meta()
  m$n_subjects <- -1
  expect_error(
    new_immuno_sim(make_data(), make_truth(make_data()), m),
    "meta\\$n_subjects.*positive"
  )

  m$n_subjects <- NULL
  expect_error(
    new_immuno_sim(make_data(), make_truth(make_data()), m),
    "meta\\$n_subjects.*positive"
  )
})


# ---- Print method tests ------------------------------------------------------

test_that("print.immuno_sim produces expected output", {
  sim <- make_sim(n_sub = 10, n_cyt = 4)
  out <- capture.output(print(sim))
  expect_length(out, 1)
  expect_match(out, "immuno_sim:")
  expect_match(out, "pre_post")
  expect_match(out, "10 subjects")
  expect_match(out, "4 analytes")
  expect_match(out, "signal analyte")
})

test_that("print returns object invisibly", {
  sim <- make_sim()
  expect_invisible(print(sim))
  ret <- print(sim)
  expect_identical(ret, sim)
})


# ---- Summary method tests ----------------------------------------------------

test_that("summary.immuno_sim produces structured output", {
  sim <- make_sim(n_sub = 10, n_cyt = 4)
  out <- capture.output(summary(sim))
  expect_true(any(grepl("Simulated Immunoassay Summary", out)))
  expect_true(any(grepl("Design:", out)))
  expect_true(any(grepl("Subjects:", out)))
  expect_true(any(grepl("Analytes:", out)))
  expect_true(any(grepl("Signal:", out)))
})

test_that("summary returns list invisibly", {
  sim <- make_sim()
  ret <- summary(sim)
  expect_type(ret, "list")
  expect_named(ret, c("design", "n_subjects", "n_analytes",
                       "n_timepoints", "n_obs"))
  expect_equal(ret$design, "pre_post")
})

test_that("summary handles censoring info when present", {
  sim <- make_sim(n_sub = 20, n_cyt = 3)
  sim$meta$realized_censoring <- c(CYT1 = 0.3, CYT2 = 0.0, CYT3 = 0.1)
  out <- capture.output(summary(sim))
  expect_true(any(grepl("Censoring", out)))
  expect_true(any(grepl("CYT1", out)))
})


# ---- Subset method tests -----------------------------------------------------

test_that("[.immuno_sim subsets data and truth in sync", {
  sim <- make_sim(n_sub = 10, n_cyt = 3)
  idx <- 1:5
  sub <- sim[idx, ]

  expect_s3_class(sub, "immuno_sim")
  expect_equal(nrow(sub$data), 5)
  expect_equal(nrow(sub$truth), 5)
  expect_equal(sub$data$value, sim$data$value[idx])
  expect_equal(sub$truth$true_latent, sim$truth$true_latent[idx])
})

test_that("[.immuno_sim with logical index", {
  sim <- make_sim(n_sub = 10, n_cyt = 2)
  keep <- sim$data$cytokine == "CYT1"
  sub <- sim[keep, ]

  expect_equal(nrow(sub$data), sum(keep))
  expect_equal(nrow(sub$truth), sum(keep))
  expect_true(all(sub$data$cytokine == "CYT1"))
})

test_that("[.immuno_sim column subset", {
  sim <- make_sim()
  sub <- sim[, c("subject_id", "cytokine", "value")]

  expect_equal(ncol(sub$data), 3)
  expect_equal(ncol(sub$truth), ncol(sim$truth))  # truth unchanged
})

test_that("[.immuno_sim updates meta counts", {
  sim <- make_sim(n_sub = 10, n_cyt = 4)
  keep <- sim$data$cytokine %in% c("CYT1", "CYT2")
  sub <- sim[keep, ]

  expect_equal(sub$meta$n_analytes, 2)
  expect_true(sub$meta$n_subjects <= sim$meta$n_subjects)
})

test_that("[.immuno_sim preserves class after subsetting", {
  sim <- make_sim(n_sub = 6, n_cyt = 2)
  sub <- sim[1:3, ]
  expect_s3_class(sub, "immuno_sim")

  # Can still print/summary without error
  expect_output(print(sub))
  expect_output(summary(sub))
})

test_that("[.immuno_sim with no row index returns all rows", {
  sim <- make_sim(n_sub = 5, n_cyt = 2)
  sub <- sim[, c("subject_id", "value")]
  expect_equal(nrow(sub$data), nrow(sim$data))
})


# ---- Edge cases --------------------------------------------------------------

test_that("single-row immuno_sim works", {
  d <- make_data(n_sub = 1, n_cyt = 1)
  tr <- make_truth(d)
  m <- make_meta(n_sub = 1, n_cyt = 1)
  sim <- new_immuno_sim(d, tr, m)

  expect_s3_class(sim, "immuno_sim")
  expect_equal(nrow(sim$data), 1)
  expect_output(print(sim), "1 subjects")
  expect_output(print(sim), "1 analytes")
})

test_that("analyte_library has expected structure", {
  lib <- immunoPlex:::analyte_library
  expect_s3_class(lib, "data.frame")
  expect_equal(nrow(lib), 20)
  expect_named(lib, c("analyte", "log_mean", "log_sd", "lod", "category"))
  expect_true(all(lib$log_mean > 0))
  expect_true(all(lib$log_sd > 0))
  expect_true(all(lib$lod > 0))
})
