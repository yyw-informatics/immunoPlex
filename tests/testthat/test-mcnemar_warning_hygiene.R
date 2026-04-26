# Warning-hygiene test for mcnemar_detection().
#
# The fix: when exact2x2 backend fits fail, the wrapper must surface the
# failure via warning(), NOT via cat() gated on `!quiet`. Otherwise users
# running with quiet = TRUE silently receive NA results.
#
# This test isolates the binding patch in its own file so it cannot leak
# into other test files via testthat's mocking mechanism.

library(testthat)

create_one_cytokine_data <- function(n_subjects = 20) {
  set.seed(42)
  expand.grid(
    subject_id = seq_len(n_subjects),
    cytokine = "IL6",
    timepoint = c("Baseline", "Follow-up")
  ) %>%
    dplyr::mutate(
      cens_lod = dplyr::case_when(
        timepoint == "Baseline" ~ runif(dplyr::n()) > 0.4,
        timepoint == "Follow-up" ~ runif(dplyr::n()) > 0.7
      )
    )
}

test_that("backend errors surface as warning() even under quiet = TRUE", {
  skip_if_not_installed("exact2x2")

  test_data <- create_one_cytokine_data()

  # Patch exact2x2::exact2x2 (the underlying engine the wrapper calls) to
  # always error, then restore on exit.
  ns <- asNamespace("exact2x2")
  orig <- get("exact2x2", envir = ns)
  on.exit({
    unlockBinding("exact2x2", ns)
    assign("exact2x2", orig, envir = ns)
    lockBinding("exact2x2", ns)
  }, add = TRUE)

  unlockBinding("exact2x2", ns)
  assign("exact2x2",
         function(...) stop("forced backend failure"),
         envir = ns)
  lockBinding("exact2x2", ns)

  warnings_seen <- character()
  withCallingHandlers(
    immunoPlex::mcnemar_detection(
      data = test_data,
      baseline = "Baseline",
      comparison = "Follow-up",
      quiet = TRUE
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(
    any(grepl("exact2x2::exact2x2 \\(mcnemar\\) failed", warnings_seen)),
    info = paste(warnings_seen, collapse = " | ")
  )
  expect_true(any(grepl("forced backend failure", warnings_seen)))
})
