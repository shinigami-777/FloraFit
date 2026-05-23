# Floral equivalent circuit impedance model.
#
# Implements the petal-specific circuit from Oduniyi (2026), Fig. 2 / Eq. 9:
#
#   Z_total(omega) = R_S + Z_CPE(omega) + [ Z_shell(omega) || R_vasc ]
#
# where
#   Z_CPE(omega)   = 1 / ( Q_c * (j*omega)^alpha_c )           # cuticle CPE (Eq. 6)
#   Z_shell(omega) = R_E + 1/(j*omega*C_M) + R_CYT             # double-shell (Eq. 8)
#                      + R_V / (1 + j*omega*R_V*C_T)
#
# The parameter vector theta packs all nine fitted values:
#   theta = c(RS, Qc, alpha_c, RE, CM, RCYT, CT, RV, R_vasc)

z_cpe <- function(omega, Qc, alpha_c) {
  1 / (Qc * (1i * omega)^alpha_c)
}

z_shell <- function(omega, RE, CM, RCYT, CT, RV) {
  s <- 1i * omega
  RE + 1 / (s * CM) + RCYT + RV / (1 + s * RV * CT)
}

z_parallel <- function(Z1, Z2) {
  (Z1 * Z2) / (Z1 + Z2)
}

# Full floral impedance as a complex vector evaluated at the supplied
# angular frequencies. `theta` is a named or positional length-9 numeric.
floral_impedance <- function(omega, theta) {
  RS      <- theta[1]
  Qc      <- theta[2]
  alpha_c <- theta[3]
  RE      <- theta[4]
  CM      <- theta[5]
  RCYT    <- theta[6]
  CT      <- theta[7]
  RV      <- theta[8]
  Rvasc   <- theta[9]

  Zc <- z_cpe(omega, Qc, alpha_c)
  Zs <- z_shell(omega, RE, CM, RCYT, CT, RV)
  Zp <- z_parallel(Zs, Rvasc)
  RS + Zc + Zp
}

# Names of parameters, kept in one place so UI/fitting/plots agree.
PARAM_NAMES <- c("RS", "Qc", "alpha_c", "RE", "CM", "RCYT", "CT", "RV", "R_vasc")
PARAM_UNITS <- c("Ohm", "S*s^a", "-", "Ohm", "F", "Ohm", "F", "Ohm", "Ohm")

# Allowed range for vacuole time constant tau_V = R_V * C_T given measured frequencies.
tau_v_bounds <- function(freq) {
  f_min <- min(freq, na.rm = TRUE)
  f_max <- max(freq, na.rm = TRUE)
  if (f_min <= 0 || f_max <= 0) {
    return(c(1e-8, 1e2))
  }
  c(
    0.02 / (2 * pi * f_max),
    50 / (2 * pi * f_min)
  )
}
