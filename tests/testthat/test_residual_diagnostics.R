test_that("compute_residual_diagnostics returns non-NA for survreg with statmod", {
  skip_if_not_installed("survival")
  skip_if_not_installed("statmod")

  set.seed(123)
  n <- 120
  age <- rnorm(n, 30, 5)
  x <- rnorm(n)
  # Lognormal AFT generating process
  mu <- 0.5 + 0.2 * x + 0.03 * age
  sigma <- 0.5
  y <- rlnorm(n, meanlog = mu, sdlog = sigma)

  lod <- quantile(y, 0.1)
  ulod <- quantile(y, 0.99)
  cens_lod <- y < lod
  cens_ulod <- y > ulod
  y_obs <- pmin(pmax(y, lod), ulod)

  df <- data.frame(
    value = y_obs,
    cens_lod = cens_lod,
    cens_ulod = cens_ulod,
    lod = rep(lod, n),
    ulod = rep(ulod, n),
    x = x,
    age = age
  )

  fit <- immunoPlex::fit_one(df, family = "aft", fixed = "x + age", random = "", ulod = TRUE)
  expect_true(inherits(fit$model, "survreg"))

  diag <- immunoPlex::compute_residual_diagnostics(fit, nsim = 200)
  expect_true(is.list(diag))
  expect_false(all(is.na(unlist(diag))))
  expect_true(is.finite(diag$p_uniform) || is.na(diag$p_uniform))
})

test_that("compute_residual_diagnostics returns non-NA for gamma glmmTMB with DHARMa", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("DHARMa")

  set.seed(456)
  n <- 150
  subject_id <- rep(1:30, each = 5)
  age <- rnorm(n, 30, 5)
  tp <- factor(rep(c("Enrollment","Delivery","Convalescent"), length.out = n))
  eta <- 0.3 + 0.2 * (tp == "Delivery") + 0.1 * (tp == "Convalescent") + 0.02 * age + rnorm(30)[subject_id] * 0.0
  mu <- exp(eta)
  value <- rgamma(n, shape = 5, rate = 5 / mu)

  df <- data.frame(
    value = value,
    cens_lod = FALSE,
    cens_ulod = FALSE,
    lod = NA_real_,
    ulod = NA_real_,
    timepoint = tp,
    age = age,
    subject_id = factor(subject_id)
  )

  fit <- immunoPlex::fit_one(df, family = "gamma", fixed = "timepoint + age", random = "(1|subject_id)")
  expect_true(inherits(fit$model, "glmmTMB"))

  diag <- immunoPlex::compute_residual_diagnostics(fit, nsim = 200)
  expect_true(is.list(diag))
  expect_false(all(is.na(unlist(diag))))
})


