# Generate synthetic EIS spectra for each fitting model in the app.
# Each spectrum is written to data/sample_<model>.csv with the columns
# expected by the Shiny app (frequency_Hz, Z_real, Z_imag).
#
# Run from the project root:
#   Rscript scripts/make_model_samples.R

source("R/circuit_models.R")

set.seed(123)

freq <- 10^seq(0, 5, length.out = 80)
omega <- 2 * pi * freq

add_noise <- function(Z, rel_sd = 0.02) {
  sd <- rel_sd * Mod(Z)
  Z + complex(
    real      = rnorm(length(Z), 0, sd),
    imaginary = rnorm(length(Z), 0, sd)
  )
}

write_sample <- function(name, Z, params) {
  df <- data.frame(
    frequency_Hz = freq,
    Z_real       = Re(Z),
    Z_imag       = Im(Z)
  )
  if (!dir.exists("data")) dir.create("data")
  path <- file.path("data", sprintf("sample_%s.csv", name))
  write.csv(df, path, row.names = FALSE)
  cat(sprintf("Wrote %s (%d rows)\n", path, nrow(df)))
  cat("  true parameters:\n")
  print(signif(params, 4))
  cat("\n")
}

# ---- 1. Single-Shell sample ------------------------------------------------
#   Z(w) = Re + Ri / (1 + j*w*Ri*Cm)
single_shell_params <- setNames(
  c(Re = 500, Ri = 5000, Cm = 1e-7),
  c("Re", "Ri", "Cm")
)
Z_ss <- single_shell_impedance(omega, single_shell_params)
write_sample("single_shell", add_noise(Z_ss), single_shell_params)

# ---- 2. Double-Shell sample ------------------------------------------------
#   Z(w) = Re + [(1 / (j*w*Cm)) || (Ri + Rv / (1 + j*w*Rv*Cv))]
double_shell_params <- setNames(
  c(Re = 200, Cm = 1e-7, Ri = 1000, Rv = 5000, Cv = 1e-6),
  c("Re", "Cm", "Ri", "Rv", "Cv")
)
Z_ds <- double_shell_impedance(omega, double_shell_params)
write_sample("double_shell", add_noise(Z_ds), double_shell_params)

# ---- 3. Cole sample --------------------------------------------------------
#   Z(w) = R_inf + (R_0 - R_inf) / (1 + (j*w*tau)^alpha)
cole_params <- setNames(
  c(R_inf = 100, R_0 = 5000, tau = 1e-3, alpha = 0.8),
  c("R_inf", "R_0", "tau", "alpha")
)
Z_co <- cole_impedance(omega, cole_params)
write_sample("cole", add_noise(Z_co), cole_params)

# ---- 4. Voigt (N=2) sample -------------------------------------------------
#   Z(w) = R0 + sum_k Rk / (1 + j*w*Rk*Ck)
voigt2_params <- setNames(
  c(R0 = 100, R1 = 500, C1 = 1e-6, R2 = 2000, C2 = 1e-4),
  c("R0", "R1", "C1", "R2", "C2")
)
Z_v2 <- voigt_impedance(omega, voigt2_params, n_elements = 2)
write_sample("voigt_n2", add_noise(Z_v2), voigt2_params)

# ---- 5. Voigt (N=3) sample -------------------------------------------------
voigt3_params <- setNames(
  c(R0 = 100, R1 = 300, C1 = 1e-7, R2 = 800, C2 = 1e-5, R3 = 2500, C3 = 1e-3),
  c("R0", "R1", "C1", "R2", "C2", "R3", "C3")
)
Z_v3 <- voigt_impedance(omega, voigt3_params, n_elements = 3)
write_sample("voigt_n3", add_noise(Z_v3), voigt3_params)

cat("Done. Use these CSVs in the app to verify each model fits its own data.\n")
