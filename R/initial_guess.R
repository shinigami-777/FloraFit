initial_guess <- function(freq, Z) {
  omega <- 2 * pi * freq
  Zre   <- Re(Z)
  Zim   <- Im(Z)

  hf_intercept <- min(Zre)
  lf_intercept <- max(Zre)

  RS_init    <- max(hf_intercept * 0.1, 50)
  RE_init    <- max(hf_intercept - RS_init, 1)
  RCYT_init  <- max(RE_init * 0.3, 1)
  Rvasc_init <- max(lf_intercept - RS_init, RE_init * 10, 1e4)

  upper_half <- which(freq > median(freq))
  if (length(upper_half) > 2) {
    j <- upper_half[which.max(-Zim[upper_half])]
    CM_init <- 1 / (omega[j] * RE_init)
  } else {
    CM_init <- 1e-8
  }
  CM_init <- min(max(CM_init, 1e-11), 1e-5)

  f_lo <- max(min(freq) * 2, 20)
  f_hi <- min(max(freq) / 50, 150)
  mid_band <- which(freq >= f_lo & freq <= f_hi)
  if (length(mid_band) > 2) {
    j <- mid_band[which.max(-Zim[mid_band])]
    tau_V_init <- 1 / (2 * pi * freq[j])
  } else {
    tau_V_init <- 1 / (2 * pi * median(freq))
  }
  tau_rng <- tau_v_bounds(freq)
  tau_V_init <- min(max(tau_V_init, tau_rng[1]), tau_rng[2])

  RV_init <- max(2.5 * RE_init, 100)
  CT_init <- min(max(tau_V_init / RV_init, 1e-11), 1e-4)

  alpha_init <- 0.85
  Qc_init    <- 1e-6

  setNames(
    c(RS_init, Qc_init, alpha_init, RE_init, CM_init,
      RCYT_init, CT_init, RV_init, Rvasc_init),
    PARAM_NAMES
  )
}

find_imag_peaks <- function(freq, Z, n_peaks = 1L) {
  ord <- order(freq)
  f <- freq[ord]
  y <- -Im(Z)[ord]
  if (length(y) < 3) {
    return(rep(median(f[f > 0]), n_peaks))
  }
  is_peak <- logical(length(y))
  for (i in seq.int(2L, length(y) - 1L)) {
    if (y[i] > y[i - 1] && y[i] > y[i + 1]) is_peak[i] <- TRUE
  }
  peak_idx <- which(is_peak)
  if (length(peak_idx) == 0) {
    peak_idx <- which.max(y)
  }
  peak_idx <- peak_idx[order(-y[peak_idx])]
  picked <- peak_idx[seq_len(min(n_peaks, length(peak_idx)))]
  if (length(picked) < n_peaks) {
    need <- n_peaks - length(picked)
    pool <- setdiff(seq_along(f), picked)
    pool <- pool[order(-y[pool])]
    picked <- c(picked, pool[seq_len(min(need, length(pool)))])
  }
  f_peaks <- f[picked]
  if (length(f_peaks) < n_peaks) {
    f_peaks <- c(f_peaks,
                 rep(median(f[f > 0]), n_peaks - length(f_peaks)))
  }
  f_peaks[seq_len(n_peaks)]
}

default_initial_guess_by_model <- function(freq, Z, model_id, voigt_n = 2L) {
  spec <- get_model_spec(model_id, voigt_n = voigt_n)
  Zre <- Re(Z)
  hf_intercept <- max(min(Zre, na.rm = TRUE), 1e-3)
  lf_intercept <- max(max(Zre, na.rm = TRUE), hf_intercept + 1)
  span <- max(lf_intercept - hf_intercept, 1)
  f_pos <- freq[freq > 0]
  f_min <- if (length(f_pos)) min(f_pos) else 1
  f_max <- if (length(f_pos)) max(f_pos) else 1e5

  if (identical(spec$id, "voigt")) {
    n <- as.integer(voigt_n)
    n <- max(1L, min(n, 6L))
    peaks <- find_imag_peaks(freq, Z, n)
    peaks <- sort(peaks, decreasing = TRUE)
    r_each <- span / n
    vals <- c(R0 = hf_intercept)
    for (k in seq_len(n)) {
      rk <- max(r_each, 1e-3)
      fk <- max(peaks[k], f_min / 2)
      ck <- 1 / (2 * pi * fk * rk)
      vals[paste0("R", k)] <- rk
      vals[paste0("C", k)] <- min(max(ck, 1e-12), 1)
    }
    return(vals[spec$param_names])
  }

  if (identical(spec$id, "single_shell")) {
    f_peak <- find_imag_peaks(freq, Z, 1L)[1]
    Ri <- max(span, 1e-3)
    Cm <- min(max(1 / (2 * pi * f_peak * Ri), 1e-12), 1)
    return(setNames(c(hf_intercept, Ri, Cm), spec$param_names))
  }

  if (identical(spec$id, "double_shell")) {
    peaks <- sort(find_imag_peaks(freq, Z, 2L), decreasing = TRUE)
    f_hi <- peaks[1]
    f_lo <- peaks[2]
    Ri <- max(0.5 * span, 1e-3)
    Rv <- max(0.5 * span, 1e-3)
    Cm <- min(max(1 / (2 * pi * f_hi * Ri), 1e-12), 1)
    Cv <- min(max(1 / (2 * pi * f_lo * Rv), 1e-12), 1)
    return(setNames(c(hf_intercept, Cm, Ri, Rv, Cv), spec$param_names))
  }

  if (identical(spec$id, "cole")) {
    f_peak <- find_imag_peaks(freq, Z, 1L)[1]
    tau <- 1 / (2 * pi * f_peak)
    tau <- min(max(tau, 1e-9), 1e3)
    return(setNames(c(hf_intercept, lf_intercept, tau, 0.8), spec$param_names))
  }

  stop("No default guess for model: ", model_id)
}
