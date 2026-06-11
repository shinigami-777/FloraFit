# Complex non-linear least squares fitting of the floral circuit.
#
# Stage 1: fit with tau_V = R_V * C_T (multi-start CNLS)
# Stage 2: fix tau_V, refine remaining parameters
# Stage 3: 1-D R_V refine (bounded) with C_T = tau_V / R_V
# Stage 4: band-limited 1-D refinements for R_S and C_M (vacuole branch fixed)

to_search <- function(theta) {
  alpha <- min(max(theta[3], 0.5001), 0.9999)
  x_alpha <- log((alpha - 0.5) / (1 - alpha))
  tau_V <- max(theta[8] * theta[7], 1e-18)
  log_pos <- log10(c(
    theta[1], theta[2], theta[4], theta[5], theta[6],
    tau_V, theta[8], theta[9]
  ))
  c(log_pos, x_alpha)
}

from_search <- function(x) {
  log_vals <- x[seq_len(length(x) - 1)]
  x_alpha  <- x[length(x)]
  alpha    <- 0.5 + 0.5 * (1 / (1 + exp(-x_alpha)))
  v        <- 10^log_vals
  RV       <- v[7]
  CT       <- v[6] / RV
  theta    <- c(v[1], v[2], alpha, v[3], v[4], v[5], CT, RV, v[8])
  setNames(theta, PARAM_NAMES)
}

to_search_fixed_tau <- function(theta, tau_V) {
  alpha <- min(max(theta[3], 0.5001), 0.9999)
  x_alpha <- log((alpha - 0.5) / (1 - alpha))
  log_pos <- log10(c(
    theta[1], theta[2], theta[4], theta[5], theta[6],
    theta[8], theta[9]
  ))
  c(log_pos, x_alpha)
}

from_search_fixed_tau <- function(x, tau_V) {
  log_vals <- x[seq_len(length(x) - 1)]
  x_alpha  <- x[length(x)]
  alpha    <- 0.5 + 0.5 * (1 / (1 + exp(-x_alpha)))
  v        <- 10^log_vals
  RV       <- v[6]
  CT       <- tau_V / RV
  theta    <- c(v[1], v[2], alpha, v[3], v[4], v[5], CT, RV, v[7])
  setNames(theta, PARAM_NAMES)
}

clip_search <- function(x, lower_log, upper_log) {
  n <- length(lower_log)
  x[seq_len(n)] <- pmin(pmax(x[seq_len(n)], lower_log), upper_log)
  x
}

arc_frequency_weights <- function(freq, f_c, sigma_log = 0.65) {
  w <- exp(-0.5 * (log10(pmax(freq, f_c / 1000) / f_c) / sigma_log)^2)
  w / mean(w)
}

residuals_fn <- function(x, omega, Z_meas, tau_V_fixed = NULL) {
  theta <- if (!is.null(tau_V_fixed)) {
    from_search_fixed_tau(x, tau_V_fixed)
  } else {
    from_search(x)
  }
  Z_model <- floral_impedance(omega, theta)
  w <- 1 / pmax(Mod(Z_meas), 1)
  diff <- (Z_meas - Z_model) * w
  c(Re(diff), Im(diff))
}

fit_rss <- function(Z_meas, Z_model) {
  sum(Mod(Z_meas - Z_model)^2)
}

default_log_bounds <- function(freq) {
  tau_rng <- tau_v_bounds(freq)
  tau_v_hi <- min(tau_rng[2], 0.2)
  list(
    lower = c(
      RS = 0, Qc = -9, RE = 1, CM = -11, RCYT = 0,
      tau_V = log10(tau_rng[1]), RV = 2.7, Rvasc = 3.5
    ),
    upper = c(
      RS = 4, Qc = -2, RE = 7, CM = -4, RCYT = 6,
      tau_V = log10(tau_v_hi), RV = 4.3, Rvasc = 6.5
    ),
    lower_fixed_tau = c(
      RS = 0, Qc = -9, RE = 1, CM = -11, RCYT = 0,
      RV = 2.7, Rvasc = 3.5
    ),
    upper_fixed_tau = c(
      RS = 4, Qc = -2, RE = 7, CM = -4, RCYT = 6,
      RV = 4.3, Rvasc = 6.5
    )
  )
}

estimate_rv_prior <- function(freq, Z, theta) {
  omega <- 2 * pi * freq
  Zc <- z_cpe(omega, theta["Qc"], theta["alpha_c"])
  Zs0 <- theta["RE"] + 1 / (1i * omega * theta["CM"]) + theta["RCYT"]
  Zp <- z_parallel(Zs0, theta["R_vasc"])
  Z_bg <- theta["RS"] + Zc + Zp
  R <- Z - Z_bg
  tau_V <- max(theta["RV"] * theta["CT"], 1e-12)
  f_c <- 1 / (2 * pi * tau_V)
  band <- which(freq >= f_c / 5 & freq <= f_c * 20)
  rv_sem <- if (length(band) >= 3) {
    2 * max(-Im(R[band]), na.rm = TRUE)
  } else {
    NA_real_
  }
  rv_scale <- 2.5 * theta["RE"]
  rv <- median(c(rv_sem, rv_scale, theta["RV"]), na.rm = TRUE)
  min(max(rv, 1.2 * theta["RE"]), 8 * theta["RE"])
}

vacuole_rv_bounds <- function(theta, rv_prior = NULL,
                              ratio = c(2.0, 2.6),
                              RE_floor = 1500) {
  RE <- max(theta["RE"], RE_floor)
  lo <- ratio[1] * RE
  hi <- ratio[2] * RE
  if (!is.null(rv_prior) && is.finite(rv_prior)) {
    lo <- max(lo, 0.75 * rv_prior)
    hi <- min(hi, 1.15 * rv_prior)
    if (lo >= hi) {
      lo <- 0.85 * rv_prior
      hi <- 1.15 * rv_prior
    }
  }
  c(lo, hi)
}

fit_once <- function(x0, omega, Z_meas, lower, upper,
                     maxiter = 500, tau_V_fixed = NULL) {
  ctrl <- minpack.lm::nls.lm.control(maxiter = maxiter, ftol = 1e-12,
                                     ptol = 1e-12)
  lm_args <- list(
    par = x0, fn = residuals_fn,
    lower = lower, upper = upper,
    omega = omega, Z_meas = Z_meas,
    control = ctrl
  )
  if (!is.null(tau_V_fixed)) {
    lm_args$tau_V_fixed <- tau_V_fixed
  }
  fit <- do.call(minpack.lm::nls.lm, lm_args)
  theta_hat <- if (!is.null(tau_V_fixed)) {
    from_search_fixed_tau(fit$par, tau_V_fixed)
  } else {
    from_search(fit$par)
  }
  Z_hat <- floral_impedance(omega, theta_hat)
  list(
    par       = fit$par,
    theta     = theta_hat,
    Z_fit     = Z_hat,
    rss       = fit_rss(Z_meas, Z_hat),
    niter     = fit$niter,
    converged = fit$info %in% 1:4,
    message   = fit$message
  )
}

multi_start_points <- function(x0, lower_log, upper_log, n_extra = 3L) {
  starts <- list(clip_search(x0, lower_log, upper_log))
  if (n_extra < 1L) {
    return(starts)
  }
  deltas <- c(-0.25, 0.25, 0.5)[seq_len(min(n_extra, 3L))]
  for (d in deltas) {
    xs <- x0
    xs[6] <- x0[6] + d
    xs[7] <- x0[7] + 0.3 * d
    starts[[length(starts) + 1L]] <- clip_search(xs, lower_log, upper_log)
  }
  starts
}

pick_best_start <- function(results, rss_ref, min_improve = 0.01) {
  best <- results[[1]]
  for (r in results[-1]) {
    if ((rss_ref - r$rss) / rss_ref >= min_improve && r$rss < best$rss) {
      best <- r
    }
  }
  best
}

refine_fixed_tau_v <- function(freq, Z, theta, bounds, maxiter = 300) {
  omega <- 2 * pi * freq
  tau_V <- theta["RV"] * theta["CT"]
  x0 <- clip_search(
    to_search_fixed_tau(theta, tau_V),
    bounds$lower_fixed_tau,
    bounds$upper_fixed_tau
  )
  lower <- c(bounds$lower_fixed_tau, -6)
  upper <- c(bounds$upper_fixed_tau,  6)
  fit_once(x0, omega, Z, lower, upper, maxiter, tau_V_fixed = tau_V)
}

refine_vacuole_rv <- function(freq, Z, theta, rv_re_ratio = c(2.0, 2.6)) {
  omega <- 2 * pi * freq
  tau_V <- theta["RV"] * theta["CT"]
  rv_prior <- estimate_rv_prior(freq, Z, theta)
  rv_rng <- vacuole_rv_bounds(theta, rv_prior, ratio = rv_re_ratio)

  rss_rv <- function(log_rv) {
    RV <- 10^log_rv
    th <- theta
    th["RV"] <- RV
    th["CT"] <- tau_V / RV
    Zm <- floral_impedance(omega, th)
    w <- (1 / pmax(Mod(Z), 1)) * arc_frequency_weights(
      freq, 1 / (2 * pi * tau_V)
    )
    sum(w * Mod(Z - Zm)^2)
  }
  opt <- optimize(rss_rv, interval = log10(rv_rng), maximum = FALSE)
  RV_hat <- 10^opt$minimum
  theta_out <- theta
  theta_out["RV"] <- RV_hat
  theta_out["CT"] <- tau_V / RV_hat
  Z_hat <- floral_impedance(omega, theta_out)
  list(
    theta     = theta_out,
    Z_fit     = Z_hat,
    rss       = fit_rss(Z, Z_hat),
    niter     = 1L,
    converged = TRUE,
    message   = "Vacuole R_V 1-D refine"
  )
}

refine_param_1d <- function(freq, Z, theta, param, log_interval, idx = NULL) {
  omega <- 2 * pi * freq
  if (log_interval[1] >= log_interval[2]) {
    return(list(theta = theta, rss = fit_rss(Z, floral_impedance(omega, theta))))
  }
  rss <- function(log_v) {
    th <- theta
    th[param] <- 10^log_v
    ii <- if (is.null(idx)) seq_along(freq) else idx
    sum(Mod(Z[ii] - floral_impedance(omega[ii], th))^2)
  }
  opt <- optimize(rss, log_interval)
  theta_out <- theta
  theta_out[param] <- 10^opt$minimum
  list(
    theta = theta_out,
    rss   = rss(opt$minimum)
  )
}

refine_post_vacuole <- function(freq, Z, theta, rss_tol = 1.012) {
  omega <- 2 * pi * freq
  base_rss <- fit_rss(Z, floral_impedance(omega, theta))

  accept <- function(res) {
    is.finite(res$rss) && res$rss <= base_rss * rss_tol
  }

  mem <- which(freq > 800 & freq < 30000)
  r_re <- refine_param_1d(
    freq, Z, theta, "RE",
    log10(c(0.85 * theta["RE"], 1.15 * theta["RE"])),
    mem
  )
  if (accept(r_re)) {
    theta <- r_re$theta
    base_rss <- r_re$rss
  }

  r_cm <- refine_param_1d(
    freq, Z, theta, "CM",
    log10(c(5e-9, 3e-7)),
    mem
  )
  if (accept(r_cm)) {
    theta <- r_cm$theta
    base_rss <- r_cm$rss
  }

  hf <- which(freq > 5000)
  r_rs <- refine_param_1d(
    freq, Z, theta, "RS",
    log10(c(max(theta["RS"] * 0.6, 80), min(theta["RS"] * 1.8, 400))),
    hf
  )
  if (accept(r_rs)) {
    theta <- r_rs$theta
    base_rss <- r_rs$rss
  }

  r_cyt <- refine_param_1d(
    freq, Z, theta, "RCYT",
    log10(c(0.2 * theta["RE"], 0.35 * theta["RE"]))
  )
  if (accept(r_cyt)) {
    theta <- r_cyt$theta
    base_rss <- r_cyt$rss
  }

  Z_hat <- floral_impedance(omega, theta)
  list(
    theta     = theta,
    Z_fit     = Z_hat,
    rss       = fit_rss(Z, Z_hat),
    niter     = 0L,
    converged = TRUE,
    message   = "Post-vacuole band 1-D refine (RE, CM, RS, RCYT)"
  )
}

fit_floral_circuit <- function(freq, Z, theta_init = NULL,
                               maxiter = 500, multi_start = TRUE,
                               refine_tau = TRUE,
                               rv_re_ratio = c(2.0, 2.6)) {
  if (is.null(theta_init)) theta_init <- initial_guess(freq, Z)
  omega <- 2 * pi * freq
  bounds <- default_log_bounds(freq)
  lower_log <- bounds$lower
  upper_log <- bounds$upper
  lower <- c(lower_log, -6)
  upper <- c(upper_log,  6)

  x0 <- clip_search(to_search(theta_init), lower_log, upper_log)

  starts <- if (isTRUE(multi_start)) {
    multi_start_points(x0, lower_log, upper_log)
  } else {
    list(x0)
  }

  results <- lapply(starts, function(x_start) {
    tryCatch(
      fit_once(x_start, omega, Z, lower, upper, maxiter),
      error = function(e) NULL
    )
  })
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0L) {
    stop("All fit attempts failed")
  }

  best <- if (length(results) == 1L) {
    results[[1]]
  } else {
    pick_best_start(results, results[[1]]$rss)
  }

  if (isTRUE(refine_tau)) {
    refined <- tryCatch(
      refine_fixed_tau_v(freq, Z, best$theta, bounds, maxiter),
      error = function(e) NULL
    )
    if (!is.null(refined) && refined$rss <= best$rss * 1.001) {
      best <- refined
    }

    vac <- tryCatch(
      refine_vacuole_rv(freq, Z, best$theta, rv_re_ratio = rv_re_ratio),
      error = function(e) NULL
    )
    if (!is.null(vac)) {
      best <- vac
    }

    post <- tryCatch(
      refine_post_vacuole(freq, Z, best$theta),
      error = function(e) NULL
    )
    if (!is.null(post)) {
      best <- post
    }
  }

  Z_hat <- best$Z_fit
  rss   <- best$rss
  tss   <- sum(Mod(Z - mean(Z))^2)
  r2    <- 1 - rss / tss
  rmse  <- sqrt(mean(Mod(Z - Z_hat)^2))
  tau_V <- best$theta["RV"] * best$theta["CT"]
  tau_M <- best$theta["RE"] * best$theta["CM"]

  list(
    theta      = best$theta,
    theta_init = theta_init,
    Z_fit      = Z_hat,
    r2         = r2,
    rmse       = rmse,
    niter      = best$niter,
    converged  = best$converged,
    message    = best$message,
    n_starts   = length(starts),
    tau_V      = tau_V,
    tau_M      = tau_M,
    refined    = isTRUE(refine_tau)
  )
}

# ---- generic model fitter used by Shiny app --------------------------------
#
# Strategy: log-space Levenberg-Marquardt with multi-start.
# All positive parameters (R, C, tau, Q) are optimized in log10 space so a
# single LM step can move them across many orders of magnitude. Parameters
# bounded on (0, 1) (e.g. Cole alpha) use a logit transform. Multi-start
# perturbs the initial guess in log-space and keeps the lowest-RSS result,
# which makes the fitter robust to poor user guesses on any of the circuits.

is_bounded_unit <- function(param_name) {
  param_name == "alpha"
}

param_to_x <- function(theta, lower, upper, param_names) {
  vapply(param_names, function(p) {
    val <- as.numeric(theta[[p]])
    if (is_bounded_unit(p)) {
      lo <- max(as.numeric(lower[[p]]), 1e-6)
      hi <- min(as.numeric(upper[[p]]), 1 - 1e-6)
      val <- min(max(val, lo + 1e-6), hi - 1e-6)
      u <- (val - lo) / (hi - lo)
      log(u / (1 - u))
    } else {
      val <- max(val, 1e-15)
      log10(val)
    }
  }, numeric(1))
}

x_to_param <- function(x, lower, upper, param_names) {
  out <- vapply(seq_along(param_names), function(i) {
    p <- param_names[i]
    if (is_bounded_unit(p)) {
      lo <- max(as.numeric(lower[[p]]), 1e-6)
      hi <- min(as.numeric(upper[[p]]), 1 - 1e-6)
      u <- 1 / (1 + exp(-x[i]))
      lo + (hi - lo) * u
    } else {
      10^x[i]
    }
  }, numeric(1))
  names(out) <- param_names
  out
}

fit_impedance_model <- function(freq, Z, model_id, theta_init = NULL,
                                voigt_n = 2L, maxiter = 500,
                                n_starts = 12L, seed = 1L) {
  spec <- get_model_spec(model_id, voigt_n = voigt_n)
  omega <- 2 * pi * freq

  if (is.null(theta_init)) {
    theta_init <- default_initial_guess_by_model(freq, Z, model_id,
                                                 voigt_n = voigt_n)
  }
  theta_init <- theta_init[spec$param_names]

  lower <- spec$lower[spec$param_names]
  upper <- spec$upper[spec$param_names]

  theta_init_clipped <- vapply(spec$param_names, function(p) {
    v <- as.numeric(theta_init[[p]])
    if (is_bounded_unit(p)) {
      min(max(v, as.numeric(lower[[p]])), as.numeric(upper[[p]]))
    } else {
      min(max(v, as.numeric(lower[[p]])), as.numeric(upper[[p]]))
    }
  }, numeric(1))
  names(theta_init_clipped) <- spec$param_names

  weights <- 1 / pmax(Mod(Z), 1)

  residual_fn <- function(x) {
    theta <- x_to_param(x, lower, upper, spec$param_names)
    Z_fit <- spec$impedance(omega, theta)
    if (any(!is.finite(Z_fit))) {
      return(rep(1e6, 2 * length(omega)))
    }
    diff <- (Z - Z_fit) * weights
    c(Re(diff), Im(diff))
  }

  run_lm <- function(x0) {
    ctrl <- minpack.lm::nls.lm.control(maxiter = maxiter, ftol = 1e-12,
                                       ptol = 1e-12, factor = 100)
    tryCatch({
      fit <- minpack.lm::nls.lm(par = x0, fn = residual_fn, control = ctrl)
      theta_hat <- x_to_param(fit$par, lower, upper, spec$param_names)
      Z_hat <- spec$impedance(omega, theta_hat)
      if (any(!is.finite(Z_hat))) return(NULL)
      list(
        theta = theta_hat,
        Z_fit = Z_hat,
        rss = fit_rss(Z, Z_hat),
        niter = fit$niter,
        converged = fit$info %in% 1:4,
        message = fit$message
      )
    }, error = function(e) NULL)
  }

  x0_base <- param_to_x(theta_init_clipped, lower, upper, spec$param_names)
  n_starts <- max(1L, as.integer(n_starts))

  set.seed(seed)
  starts <- vector("list", n_starts)
  starts[[1]] <- x0_base
  if (n_starts > 1L) {
    for (i in seq.int(2L, n_starts)) {
      perturb <- stats::rnorm(length(x0_base), mean = 0,
                              sd = ifelse(is_bounded_unit(spec$param_names),
                                          0.5, 0.7))
      starts[[i]] <- x0_base + perturb
    }
  }

  results <- Filter(Negate(is.null), lapply(starts, run_lm))
  if (length(results) == 0L) {
    stop("All fit attempts failed for model: ", spec$id)
  }
  rss_vec <- vapply(results, function(r) r$rss, numeric(1))
  best <- results[[which.min(rss_vec)]]

  theta <- best$theta
  Z_hat <- best$Z_fit
  rss <- best$rss
  tss <- sum(Mod(Z - mean(Z))^2)
  r2 <- 1 - rss / tss
  rmse <- sqrt(mean(Mod(Z - Z_hat)^2))

  extra <- character(0)
  if (identical(spec$id, "single_shell")) {
    extra <- c(extra, sprintf("tau_m [s]:  %s  (Ri * Cm)",
                              signif(theta["Ri"] * theta["Cm"], 4)))
  } else if (identical(spec$id, "double_shell")) {
    extra <- c(
      extra,
      sprintf("tau_m [s]:  %s  (Ri * Cm)",
              signif(theta["Ri"] * theta["Cm"], 4)),
      sprintf("tau_v [s]:  %s  (Rv * Cv)",
              signif(theta["Rv"] * theta["Cv"], 4))
    )
  } else if (identical(spec$id, "cole")) {
    extra <- c(extra, sprintf("alpha:      %s", signif(theta["alpha"], 4)))
  } else if (identical(spec$id, "voigt")) {
    tau_lines <- vapply(seq_len(as.integer(voigt_n)), function(k) {
      rk <- theta[paste0("R", k)]
      ck <- theta[paste0("C", k)]
      sprintf("tau_%d [s]:  %s  (R%d * C%d)", k,
              signif(rk * ck, 4), k, k)
    }, character(1))
    extra <- c(extra, tau_lines)
  }

  list(
    model_id = spec$id,
    model_label = spec$label,
    equation = spec$equation,
    param_names = spec$param_names,
    param_units = spec$param_units,
    theta = theta,
    theta_init = theta_init_clipped,
    Z_fit = Z_hat,
    r2 = r2,
    rmse = rmse,
    niter = best$niter,
    converged = best$converged,
    message = best$message,
    n_starts = length(results),
    extra_summary = extra
  )
}
