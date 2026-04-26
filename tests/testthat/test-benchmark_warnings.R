# Phase 1 benchmark-derived runtime warnings.
# See inst/simulations/reports/{B1,B2,B3,B5}*.md for the empirical basis.

library(testthat)

.make_censored_dat <- function(n = 30, cens_frac = 0.4, seed = 1) {
  set.seed(seed)
  val <- exp(rnorm(n, mean = 1, sd = 0.8))
  lod <- quantile(val, probs = cens_frac, names = FALSE)
  cens_lod <- val <= lod
  data.frame(
    value     = log(pmax(val, lod)),
    cens_lod  = cens_lod,
    cens_ulod = FALSE,
    lod       = rep(lod, n),
    ulod      = rep(Inf, n),
    timepoint = factor(rep(c("pre", "post"), length.out = n)),
    disease   = factor(rep(c("ctrl", "case"), length.out = n)),
    age       = rnorm(n, 40, 10),
    subject_id = factor(seq_len(n))
  )
}

describe("Gamma + censored data emits recommendation warning (B1/B5)", {

  test_that("warning fires when any censoring is present", {
    skip_if_not_installed("glmmTMB")
    dat <- .make_censored_dat(n = 30, cens_frac = 0.4)
    expect_warning(
      fit_one(dat, family = "gamma", random = "",
              fixed = "disease + age"),
      regexp = "family='gamma' is not recommended"
    )
  })

  test_that("warning does NOT fire when no censoring is present", {
    skip_if_not_installed("glmmTMB")
    dat <- .make_censored_dat(n = 30, cens_frac = 0)
    dat$cens_lod <- FALSE
    expect_silent(
      suppressMessages(
        fit_one(dat, family = "gamma", random = "",
                fixed = "disease + age")
      )
    )
  })

  test_that("auto-selector never falls back to gamma for censored data", {
    skip_if_not_installed("glmmTMB")
    dat <- .make_censored_dat(n = 30, cens_frac = 0.4)
    fit <- suppressWarnings(suppressMessages(
      fit_one(dat, family = "auto", random = "",
              fixed = "disease + age")
    ))
    expect_false(identical(fit$family, "gamma"))
  })
})

describe("Tobit / AFT low-n uncensored warning (B1/B5)", {

  test_that("Tobit warns when uncensored count < 10", {
    skip_if_not_installed("survival")
    dat <- .make_censored_dat(n = 15, cens_frac = 0.7)
    expect_warning(
      fit_one(dat, family = "tobit", random = "",
              fixed = "disease + age"),
      regexp = "Tobit fit has only [0-9]+ uncensored"
    )
  })

  test_that("Tobit does NOT warn when uncensored count >= 10", {
    skip_if_not_installed("survival")
    dat <- .make_censored_dat(n = 40, cens_frac = 0.3)
    # 28 uncensored — above threshold
    warns <- capture_warnings(
      suppressMessages(
        fit_one(dat, family = "tobit", random = "",
                fixed = "disease + age")
      )
    )
    expect_false(any(grepl("uncensored observation", warns)))
  })

  test_that("AFT warns when uncensored count < 10", {
    skip_if_not_installed("survival")
    dat <- .make_censored_dat(n = 15, cens_frac = 0.7)
    expect_warning(
      fit_one(dat, family = "aft", random = "",
              fixed = "disease + age"),
      regexp = "AFT fit has only [0-9]+ uncensored"
    )
  })
})

describe("ANCOVA tight-fence guard (B3)", {

  test_that("iqr_multiplier < 2 at n < 50 emits warning", {
    set.seed(42)
    n <- 30
    dat <- data.frame(
      outcome   = rnorm(n),
      covariate = rnorm(n),
      grp       = factor(rep(c("A", "B"), length.out = n))
    )
    expect_warning(
      ancova_one(
        data = dat, outcome = "outcome", group = "grp",
        covariate = "covariate", iqr_multiplier = 1.5
      ),
      regexp = "iqr_multiplier=1.5.*not recommended"
    )
  })

  test_that("default iqr_multiplier = 3 does NOT warn", {
    set.seed(42)
    n <- 30
    dat <- data.frame(
      outcome   = rnorm(n),
      covariate = rnorm(n),
      grp       = factor(rep(c("A", "B"), length.out = n))
    )
    warns <- capture_warnings(
      ancova_one(data = dat, outcome = "outcome", group = "grp",
                 covariate = "covariate")
    )
    expect_false(any(grepl("iqr_multiplier", warns)))
  })

  test_that("iqr_multiplier < 2 at n >= 50 does NOT warn (large-sample OK)", {
    set.seed(42)
    n <- 60
    dat <- data.frame(
      outcome   = rnorm(n),
      covariate = rnorm(n),
      grp       = factor(rep(c("A", "B"), length.out = n))
    )
    warns <- capture_warnings(
      ancova_one(data = dat, outcome = "outcome", group = "grp",
                 covariate = "covariate", iqr_multiplier = 1.5)
    )
    expect_false(any(grepl("iqr_multiplier.*not recommended", warns)))
  })
})
