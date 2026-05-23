# Generate a synthetic EIS spectrum from known floral-circuit parameters
# Gaussian noise 5%. added, and write data/sample_eis.csv.
#
# Run from the project root:    Rscript scripts/make_sample_data.R

source("R/circuit_models.R")

set.seed(42)

theta_true <- setNames(
  c(RS = 200, Qc = 1e-6, alpha_c = 0.85,
    RE = 2000, CM = 1e-7,
    RCYT = 500, CT = 1e-6,
    RV = 5000, R_vasc = 5e4),
  PARAM_NAMES
)

freq  <- 10^seq(0, 5, length.out = 60)    # 1 Hz to 100 kHz, log-spaced
omega <- 2 * pi * freq
Z     <- floral_impedance(omega, theta_true)

# Adding som enoise to the data
noise_sd <- 0.05 * Mod(Z)
Z_noisy <- Z + complex(real = rnorm(length(Z), 0, noise_sd),
                       imaginary = rnorm(length(Z), 0, noise_sd))

out <- data.frame(
  frequency_Hz = freq,
  Z_real       = Re(Z_noisy),
  Z_imag       = Im(Z_noisy)
)

dir.create("data")
write.csv(out, "data/sample_eis.csv", row.names = FALSE)

cat("Wrote data/sample_eis.csv with", nrow(out), "rows\n")
cat("True parameters used:\n")
print(theta_true)
