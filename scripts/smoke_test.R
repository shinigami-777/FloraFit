# Smoke test: load the synthetic CSV, run the fitter, print recovered
# parameters vs the truth used in make_sample_data.R.

source("R/circuit_models.R")
source("R/initial_guess.R")
source("R/fitting.R")

df <- read.csv("data/sample_eis.csv")
freq <- df$frequency_Hz
Z    <- complex(real = df$Z_real, imaginary = df$Z_imag)

theta_true <- setNames(
  c(200, 1e-6, 0.85, 2000, 1e-7, 500, 1e-6, 5000, 5e4),
  PARAM_NAMES
)

t0 <- Sys.time()
fit <- fit_floral_circuit(freq, Z)
dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cat(sprintf("\nFit converged: %s  in %d iter  (%.2f s)\n",
            fit$converged, fit$niter, dt))
cat(sprintf("R^2 = %.4f   RMSE = %.3g Ohm   starts = %d\n",
            fit$r2, fit$rmse, fit$n_starts))
cat(sprintf("tau_V = %.4g s  (true %.4g s)\n", fit$tau_V, theta_true["RV"] * theta_true["CT"]))
cat(sprintf("tau_M = %.4g s  (true %.4g s)\n", fit$tau_M, theta_true["RE"] * theta_true["CM"]))

cmp <- data.frame(
  param = PARAM_NAMES,
  true  = signif(theta_true, 4),
  init  = signif(fit$theta_init, 4),
  fit   = signif(fit$theta, 4),
  rel_err_pct = signif(100 * (fit$theta - theta_true) / theta_true, 3)
)
print(cmp, row.names = FALSE)
