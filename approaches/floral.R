install.packages(c(
  "pracma",
  "ggplot2",
  "minpack.lm",
  "dplyr",
  "tidyr"
))

library(pracma)
library(ggplot2)
library(minpack.lm)
library(dplyr)
library(tidyr)

z_r <- function(R) {
  R
}
z_c <- function(w, C) {
  1 / (1i * w * C)
}
z_cpe <- function(w, Q, alpha) {
  1 / (Q * (1i * w)^alpha)
}
parallel_z <- function(z1, z2) {
  (z1 * z2) / (z1 + z2)
}

z_shell <- function(w, PE, CM, RCYT, RV, CT) {
  membrane <- z_c(w, CM)
  vacuole <- RV / (1 + 1i * w * RV * CT)
  RE + membrane + RCYT + vacuole
}

z_total <- function(w, RS, Qc, alpha_c, RE, CM, RCYT, RV, CT, Rvasc){
  zcut <- z_cpe(w, Qc, alpha_c)
  zsh <- z_shell(w, RE, CM, RCYT, RV, CT)
  zparallel <- parallel_z(zsh, Rvasc)
  RS + zcut + zparallel
}

freq <- 10^seq(0, 5, length.out = 100)
w <- 2 * pi * freq

true_params <- list(
  RS = 200,
  Qc = 1e-5,
  alpha_c = 0.82,
  RE = 10000,
  CM = 1e-8,
  RCYT = 5000,
  RV = 3000,
  CT = 1e-7,
  Rvasc = 500
)

Z_clean <- z_total(
  w,
  true_params$RS,
  true_params$Qc,
  true_params$alpha_c,
  true_params$RE,
  true_params$CM,
  true_params$RCYT,
  true_params$RV,
  true_params$CT,
  true_params$Rvasc
)

# Adding some white noise
noise_level <- 0.02
noise_mag <- noise_level * Mod(Z_clean)

Z_noisy <- Z_clean +
  rnorm(length(Z_clean), 0, noise_mag) +
  1i * rnorm(length(Z_clean), 0, noise_mag)

# Nyquist plot
nyquist_df <- data.frame(
  Re = Re(Z_noisy),
  Im = -Im(Z_noisy)
)

ggplot(nyquist_df, aes(Re, Im)) +
  geom_point(size = 2) +
  labs(
    title = "Nyquist Plot",
    x = "Z' (Ohms)",
    y = "-Z'' (Ohms)"
  ) +
  theme_minimal()


# Bode plot
bode_df <- data.frame(
  freq = freq,
  magnitude = Mod(Z_noisy),
  phase = Arg(Z_noisy) * 180 / pi
)

ggplot(bode_df, aes(freq, magnitude)) +
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Bode Magnitude",
    x = "Frequency (Hz)",
    y = "|Z|"
  ) +
  theme_minimal()

ggplot(bode_df, aes(freq, phase)) +
  geom_line() +
  scale_x_log10() +
  labs(
    title = "Bode Phase",
    x = "Frequency (Hz)",
    y = "Phase (deg)"
  ) +
  theme_minimal()

# Parameter fitting
residual_function <- function(p, w, Zmeas) {
  Zmodel <- z_total(
    w,
    p["RS"],
    p["Qc"],
    p["alpha_c"],
    p["RE"],
    p["CM"],
    p["RCYT"],
    p["RV"],
    p["CT"],
    p["Rvasc"]
  )

  c (Re(Zmodel - Zmeas),Im(Zmodel - Zmeas))
}

initial_guess <- c(
  RS = 1000,
  Qc = 1e-6,
  alpha_c = 0.9,
  RE = 5000,
  CM = 1e-9,
  RCYT = 2000,
  RV = 1000,
  CT = 1e-8,
  Rvasc = 1000
)

# fit using levenberg-marquardt
fit <- nls.lm(
  par = initial_guess,
  fn = residual_function,
  w = w,
  Zmeas = Z_noisy,
  control = nls.lm.control(
    maxiter = 200
  )
)
fitted_params <- fit$par
print(fitted_params)

# Compare fitting
Z_fit <- z_total(
  w,
  fitted_params["RS"],
  fitted_params["Qc"],
  fitted_params["alpha_c"],
  fitted_params["RE"],
  fitted_params["CM"],
  fitted_params["RCYT"],
  fitted_params["RV"],
  fitted_params["CT"],
  fitted_params["Rvasc"]
)

# Compare nyquist
compare_df <- data.frame(
  Re = c(Re(Z_noisy), Re(Z_fit)),
  Im = c(-Im(Z_noisy), -Im(Z_fit)),
  Type = rep(
    c("Measured", "Fitted"),
    each = length(Z_fit)
  )
)

ggplot(compare_df, aes(Re, Im, color = Type)) +
  geom_point() +
  geom_line() +
  labs(
    title = "Measured vs Fitted Nyquist Plot",
    x = "Z' (Ohms)",
    y = "-Z'' (Ohms)"
  ) +
  theme_minimal()


# Compare RMSE
rmse <- sqrt(
  mean(Mod(Z_fit - Z_noisy)^2)
)

print(rmse)
