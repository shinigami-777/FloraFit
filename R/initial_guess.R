# Graphical initial-guess routine.
#
# In the paper's circuit (Eq. 8/9) the membrane capacitance CM is in
# series, so Z_shell -> infinity at DC. The DC limit of the full network
# is therefore set by the vascular bypass: Re(Z) at low f ~ R_S + R_vasc.
# (The paper's text on p. 4 calls this "R_S + R_E + R_CYT + R_V"; that
# matches a standard double-shell with CM in parallel, not the equation
# given. We follow the equation, which is what the model uses.)

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
