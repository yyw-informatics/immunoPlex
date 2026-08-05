# Utility diagnostics for model fitness: residual uniformity, CV metrics, and PI coverage

#' Compute residual diagnostics (uniformity p-value and dispersion proxy)
#'
#' Returns a list with fields `p_uniform` and `dispersion_proxy`.
#' Uses DHARMa for glmmTMB if available; uses randomized quantile residuals for survreg;
#' falls back to NA on failure.
#' @importFrom stats Gamma predict simulate
#' @param fit an `immuno_fit` object
#' @param nsim integer, simulations for DHARMa when applicable
#' @export
compute_residual_diagnostics <- function(fit, nsim = getOption("immunoplex.dharma_nsim", 1000L)) {
  out <- list(p_uniform = NA_real_, dispersion_proxy = NA_real_)
  if (is.null(fit) || !isTRUE(fit$converged) || is.null(fit$model)) return(out)

  # glmmTMB path
  if (inherits(fit$model, "glmmTMB")) {
    if (requireNamespace("DHARMa", quietly = TRUE)) {
      dh <- try(DHARMa::simulateResiduals(fit$model, n = as.integer(max(100L, nsim)), plot = FALSE), silent = TRUE)
      if (!inherits(dh, "try-error")) {
        tu <- try(DHARMa::testUniformity(dh), silent = TRUE)
        if (!inherits(tu, "try-error")) out$p_uniform <- suppressWarnings(tu$p.value)
        # Dispersion test; if fails, compute variance proxy of scaled residuals
        td <- try(DHARMa::testDispersion(dh), silent = TRUE)
        if (!inherits(td, "try-error")) {
          out$dispersion_proxy <- suppressWarnings(td$p.value)
        } else {
          out$dispersion_proxy <- tryCatch(stats::var(dh$scaledResiduals, na.rm = TRUE), error = function(e) NA_real_)
        }
      }
    }
    return(out)
  }

  # survreg path (tobit/aft)
  if (inherits(fit$model, "survreg")) {
    # 1) Try statmod RQR if available (guard against unsupported classes)
    if (requireNamespace("statmod", quietly = TRUE)) {
      rqr_try <- try(statmod::qresiduals(fit$model), silent = TRUE)
      if (!inherits(rqr_try, "try-error")) {
        u_try <- try(stats::pnorm(rqr_try), silent = TRUE)
        if (!inherits(u_try, "try-error")) {
          ks <- try(suppressWarnings(stats::ks.test(u_try, "punif")), silent = TRUE)
          if (!inherits(ks, "try-error")) out$p_uniform <- as.numeric(ks$p.value)
          out$dispersion_proxy <- tryCatch(stats::var(rqr_try, na.rm = TRUE), error = function(e) NA_real_)
          return(out)
        }
      }
    }

    # 2) Deterministic fallback: Dunn-Smyth-style RQR via mid-interval probability
    manual <- try({
      dat <- try(fit$data_used, silent = TRUE)
      n <- if (inherits(dat, "try-error") || is.null(dat)) 0L else nrow(dat)
      if (n <= 0L) stop("no data_used attached to fit")

      # Reconstruct censoring intervals [left, right]
      if (all(c("left", "right") %in% names(dat))) {
        left  <- dat$left;  right <- dat$right
      } else if (all(c("ltime", "rtime") %in% names(dat))) {
        left  <- ifelse(is.na(dat$ltime), -Inf, dat$ltime)
        right <- ifelse(is.na(dat$rtime),  Inf,  dat$rtime)
      } else if (all(c("value", "cens_lod", "cens_ulod") %in% names(dat))) {
        left  <- ifelse(dat$cens_lod,  -Inf, dat$value)
        right <- ifelse(dat$cens_ulod,  Inf,  dat$value)
      } else {
        stop("unable to reconstruct censoring intervals from data_used")
      }

      # Linear predictor and scale
      mu    <- as.numeric(stats::predict(fit$model, type = "lp"))
      sigma <- as.numeric(fit$model$scale)
      fam   <- tolower(if (!is.null(fit$family)) fit$family else "")

      if (identical(fam, "aft")) {
        # Log-normal AFT: operate on log scale, clip non-positive
        eps <- .Machine$double.xmin
        left_z  <- ifelse(is.finite(left),  (log(pmax(left,  eps)) - mu) / sigma,  left)
        right_z <- ifelse(is.finite(right), (log(pmax(right, eps)) - mu) / sigma, right)
      } else {
        # Tobit Gaussian on identity scale
        left_z  <- (left  - mu) / sigma
        right_z <- (right - mu) / sigma
      }

      FL <- stats::pnorm(left_z)
      FR <- stats::pnorm(right_z)
      u_mid <- pmin(pmax((FL + FR) / 2, 0), 1)
      rqr   <- stats::qnorm(pmin(pmax(u_mid, .Machine$double.eps), 1 - .Machine$double.eps))
      list(u = u_mid, rqr = rqr)
    }, silent = TRUE)

    if (!inherits(manual, "try-error")) {
      ks <- try(suppressWarnings(stats::ks.test(manual$u, "punif")), silent = TRUE)
      if (!inherits(ks, "try-error")) out$p_uniform <- as.numeric(ks$p.value)
      out$dispersion_proxy <- tryCatch(stats::var(manual$rqr, na.rm = TRUE), error = function(e) NA_real_)
    }
    return(out)
  }

  # censReg external residuals unsupported here; return NA
  out
}


#' K-fold CV metrics (RMSE/MAE) for the best model family
#'
#' Attempts to refit the same family and formula on K folds and compute predictive RMSE/MAE.
#' Returns a list with numeric fields `cv_rmse` and `cv_mae`; NA if not feasible.
#' @param fit an `immuno_fit` object (best model)
#' @param k integer folds
#' @param seed integer RNG seed
#' @export
kfold_cv_metrics <- function(fit, k = 5L, seed = 123L) {
  res <- list(cv_rmse = NA_real_, cv_mae = NA_real_)
  if (is.null(fit) || !isTRUE(fit$converged) || is.null(fit$model)) return(res)
  dat <- fit$data_used
  if (!is.data.frame(dat) || nrow(dat) < (k + 3)) return(res)

  # Extract formula for refitting
  form <- try(stats::formula(fit$model), silent = TRUE)
  if (inherits(form, "try-error") || is.null(form)) {
    # glmmTMB stores separate formulas; attempt to reconstruct from call
    if (inherits(fit$model, "glmmTMB")) {
      form <- try(fit$model$call$formula, silent = TRUE)
      if (inherits(form, "try-error")) return(res)
    } else {
      return(res)
    }
  }

  # Group-wise K-fold if subject_id available; adapt k down to avoid tiny folds
  set.seed(as.integer(seed))
  if ("subject_id" %in% names(dat)) {
    groups <- unique(dat$subject_id)
    ng <- length(groups)
    k_use <- max(2L, min(as.integer(k), as.integer(ng)))
    if (ng >= k_use) {
      folds <- sample(rep(1:k_use, length.out = ng))
      fold_id <- folds[match(dat$subject_id, groups)]
    } else {
      fold_id <- sample(rep(1:k_use, length.out = nrow(dat)))
    }
  } else {
    k_use <- max(2L, min(as.integer(k), as.integer(floor(nrow(dat) / 5))))
    if (!is.finite(k_use) || k_use < 2L) k_use <- 2L
    fold_id <- sample(rep(1:k_use, length.out = nrow(dat)))
  }

  preds <- rep(NA_real_, nrow(dat))
  obs   <- dat$value

  for (i in 1:k) {
    train <- dat[fold_id != i, , drop = FALSE]
    test  <- dat[fold_id == i, , drop = FALSE]
    if (nrow(test) == 0L || nrow(train) < 5L) next

    fam <- fit$family
    fitted_obj <- try({
      if (inherits(fit$model, "glmmTMB")) {
        ctrl <- try(glmmTMB::glmmTMBControl(optimizer = stats::nlminb,
                                            optCtrl = list(iter.max = 5000, eval.max = 5000, rel.tol = 1e-8)), silent = TRUE)
        if (inherits(ctrl, "try-error")) ctrl <- glmmTMB::glmmTMBControl()
        suppressMessages(suppressWarnings(
          glmmTMB::glmmTMB(formula = form, data = train, family = Gamma(link = "log"),
                           ziformula = ~0, dispformula = ~1, control = ctrl)
        ))
      } else if (inherits(fit$model, "survreg")) {
        # Determine distribution from original model
        dist <- try(fit$model$dist, silent = TRUE)
        if (inherits(dist, "try-error") || is.null(dist)) dist <- "gaussian"
        ctrl <- try(survival::survreg.control(maxiter = 200, rel.tolerance = 1e-8), silent = TRUE)
        suppressWarnings(survival::survreg(formula = form, data = train, dist = dist, control = ctrl))
      } else {
        stop("Unsupported model class")
      }
    }, silent = TRUE)

    if (inherits(fitted_obj, "try-error")) next
    p <- try(stats::predict(fitted_obj, newdata = test, type = "response"), silent = TRUE)
    if (!inherits(p, "try-error") && length(p) == nrow(test)) {
      preds[fold_id == i] <- as.numeric(p)
    }
  }

  keep <- is.finite(preds) & is.finite(obs)
  if (!any(keep)) return(res)
  err <- preds[keep] - obs[keep]
  res$cv_rmse <- sqrt(mean(err^2))
  res$cv_mae  <- mean(abs(err))
  res
}


#' Predictive interval coverage and mean width via simulation
#'
#' @param fit an `immuno_fit` object
#' @param nsim integer number of simulations
#' @param level numeric interval level (e.g., 0.95)
#' @export
simulate_pi_coverage <- function(fit, nsim = 200L, level = 0.95) {
  out <- list(pi_coverage = NA_real_, pi_mean_width = NA_real_)
  if (is.null(fit) || !isTRUE(fit$converged) || is.null(fit$model)) return(out)
  dat <- fit$data_used
  if (!is.data.frame(dat) || nrow(dat) < 3) return(out)

  alpha <- (1 - level) / 2
  # Align observed values to the rows actually used by the fitted model
  obs <- dat$value
  mf <- try(stats::model.frame(fit$model), silent = TRUE)
  if (!inherits(mf, "try-error") && is.data.frame(mf)) {
    idx <- suppressWarnings(as.integer(rownames(mf)))
    if (all(is.finite(idx))) {
      idx <- idx[idx >= 1 & idx <= nrow(dat)]
      if (length(idx)) obs <- dat$value[idx]
    } else if (nrow(mf) == nrow(dat)) {
      obs <- dat$value
    } else {
      # Fallback: truncate to common length
      len <- min(length(obs), nrow(mf))
      if (len >= 1) obs <- obs[seq_len(len)]
    }
  }

  # glmmTMB gamma: use simulate to generate nsim reps
  if (inherits(fit$model, "glmmTMB")) {
    sims <- try(simulate(fit$model, nsim = as.integer(max(50L, nsim))), silent = TRUE)
    if (inherits(sims, "try-error")) return(out)
    sim_mat <- as.matrix(as.data.frame(sims))
    # Ensure obs length matches rows of simulation matrix
    if (length(obs) != nrow(sim_mat)) {
      len <- min(length(obs), nrow(sim_mat))
      obs <- obs[seq_len(len)]
      sim_mat <- sim_mat[seq_len(len), , drop = FALSE]
    }
    lo <- apply(sim_mat, 1, stats::quantile, probs = alpha, na.rm = TRUE)
    hi <- apply(sim_mat, 1, stats::quantile, probs = 1 - alpha, na.rm = TRUE)
    cover <- mean(obs >= lo & obs <= hi, na.rm = TRUE)
    width <- mean(hi - lo, na.rm = TRUE)
    out$pi_coverage <- as.numeric(cover)
    out$pi_mean_width <- as.numeric(width)
    return(out)
  }

  # survreg: simulate from parametric distribution at linear predictor
  if (inherits(fit$model, "survreg")) {
    lp <- try(stats::predict(fit$model, type = "lp"), silent = TRUE)
    if (inherits(lp, "try-error")) return(out)
    sigma <- try(fit$model$scale, silent = TRUE)
    if (inherits(sigma, "try-error") || is.null(sigma)) sigma <- 1
    dist <- try(fit$model$dist, silent = TRUE)
    if (inherits(dist, "try-error") || is.null(dist)) dist <- "gaussian"

    n <- length(lp)
    # Ensure obs length matches n
    if (length(obs) != n) {
      len <- min(length(obs), n)
      obs <- obs[seq_len(len)]
      n <- len
    }
    sim_mat <- matrix(NA_real_, nrow = n, ncol = as.integer(max(50L, nsim)))
    for (j in 1:ncol(sim_mat)) {
      if (identical(dist, "lognormal")) {
        sim_mat[, j] <- stats::rlnorm(n, meanlog = lp, sdlog = sigma)
      } else {
        sim_mat[, j] <- stats::rnorm(n, mean = lp, sd = sigma)
      }
    }
    lo <- apply(sim_mat, 1, stats::quantile, probs = alpha, na.rm = TRUE)
    hi <- apply(sim_mat, 1, stats::quantile, probs = 1 - alpha, na.rm = TRUE)
    cover <- mean(obs >= lo & obs <= hi, na.rm = TRUE)
    width <- mean(hi - lo, na.rm = TRUE)
    out$pi_coverage <- as.numeric(cover)
    out$pi_mean_width <- as.numeric(width)
    return(out)
  }

  out
}


