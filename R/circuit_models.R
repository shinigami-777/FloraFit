# Floral equivalent circuit impedance model.
#
# Implements the petal-specific circuit
#   Z_total(omega) = R_S + Z_CPE(omega) + [ Z_shell(omega) || R_vasc ]
#
# where
#   Z_CPE(omega)   = 1 / ( Q_c * (j*omega)^alpha_c )                                 cuticle CPE
#   Z_shell(omega) = R_E + 1/(j*omega*C_M) + R_CYT  +  R_V / (1 + j*omega*R_V*C_T)   double shell
# The parameter vector theta packs all nine fitted values:
# theta = c(RS, Qc, alpha_c, RE, CM, RCYT, CT, RV, R_vasc)

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

# Full floral impedance as a complex vector evaluated at the supplied angular frequencies. 
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

# ---- generic model support used by the Shiny app ----------------------------

MODEL_CHOICES <- c(
  "Voigt" = "voigt",
  "Single-Shell" = "single_shell",
  "Double-Shell" = "double_shell",
  "Cole" = "cole"
)

`%||%` <- function(x, y) if (is.null(x)) y else x

voigt_impedance <- function(omega, theta, n_elements) {
  Z <- rep(theta["R0"], length(omega))
  for (k in seq_len(n_elements)) {
    rk <- theta[paste0("R", k)]
    ck <- theta[paste0("C", k)]
    Z <- Z + rk / (1 + 1i * omega * rk * ck)
  }
  Z
}

single_shell_impedance <- function(omega, theta) {
  Re <- theta["Re"]
  Ri <- theta["Ri"]
  Cm <- theta["Cm"]
  Re + Ri / (1 + 1i * omega * Ri * Cm)
}

double_shell_impedance <- function(omega, theta) {
  Re <- theta["Re"]
  Cm <- theta["Cm"]
  Ri <- theta["Ri"]
  Rv <- theta["Rv"]
  Cv <- theta["Cv"]
  Z_cm <- 1 / (1i * omega * Cm)
  Z_v <- Rv / (1 + 1i * omega * Rv * Cv)
  Re + z_parallel(Z_cm, Ri + Z_v)
}

cole_impedance <- function(omega, theta) {
  R_inf <- theta["R_inf"]
  R_0   <- theta["R_0"]
  tau   <- theta["tau"]
  alpha <- theta["alpha"]
  R_inf + (R_0 - R_inf) / (1 + (1i * omega * tau)^alpha)
}

get_model_spec <- function(model_id, voigt_n = 2L) {
  model_id <- as.character(model_id %||% "voigt")
  voigt_n <- as.integer(voigt_n %||% 2L)
  voigt_n <- max(1L, min(voigt_n, 6L))

  if (identical(model_id, "voigt")) {
    param_names <- c("R0", as.vector(rbind(paste0("R", seq_len(voigt_n)),
                                           paste0("C", seq_len(voigt_n)))))
    param_units <- c("Ohm", rep(c("Ohm", "F"), voigt_n))
    lower <- setNames(c(1e-3, rep(c(1e-3, 1e-12), voigt_n)), param_names)
    upper <- setNames(c(1e8, rep(c(1e8, 1), voigt_n)), param_names)
    return(list(
      id = "voigt",
      label = "Voigt",
      param_names = param_names,
      param_units = param_units,
      lower = lower,
      upper = upper,
      equation = "Z(w)=R0+sum_k Rk/(1+j*w*Rk*Ck)",
      impedance = function(omega, theta) voigt_impedance(omega, theta, voigt_n)
    ))
  }

  if (identical(model_id, "single_shell")) {
    param_names <- c("Re", "Ri", "Cm")
    return(list(
      id = "single_shell",
      label = "Single-Shell",
      param_names = param_names,
      param_units = c("Ohm", "Ohm", "F"),
      lower = setNames(c(1e-3, 1e-3, 1e-12), param_names),
      upper = setNames(c(1e8, 1e8, 1), param_names),
      equation = "Z(w)=Re+Ri/(1+j*w*Ri*Cm)",
      impedance = function(omega, theta) single_shell_impedance(omega, theta)
    ))
  }

  if (identical(model_id, "double_shell")) {
    param_names <- c("Re", "Cm", "Ri", "Rv", "Cv")
    return(list(
      id = "double_shell",
      label = "Double-Shell",
      param_names = param_names,
      param_units = c("Ohm", "F", "Ohm", "Ohm", "F"),
      lower = setNames(c(1e-3, 1e-12, 1e-3, 1e-3, 1e-12), param_names),
      upper = setNames(c(1e8, 1, 1e8, 1e8, 1), param_names),
      equation = "Z(w)=Re+[(1/(j*w*Cm))||(Ri+Rv/(1+j*w*Rv*Cv))]",
      impedance = function(omega, theta) double_shell_impedance(omega, theta)
    ))
  }

  if (identical(model_id, "cole")) {
    param_names <- c("R_inf", "R_0", "tau", "alpha")
    return(list(
      id = "cole",
      label = "Cole",
      param_names = param_names,
      param_units = c("Ohm", "Ohm", "s", "-"),
      lower = setNames(c(1e-3, 1e-3, 1e-9, 0.1), param_names),
      upper = setNames(c(1e8, 1e8, 1e3, 0.999), param_names),
      equation = "Z(w)=Rinf+(R0-Rinf)/(1+(j*w*tau)^alpha)",
      impedance = function(omega, theta) cole_impedance(omega, theta)
    ))
  }

  stop("Unknown model id: ", model_id)
}
