# FloraFit

An R application for fitting equivalent-circuit models to electrochemical impedance spectroscopy (EIS) data, with a focus on flower petals and other thin plant tissues. The tool ships four general-purpose reference circuits (Voigt, single-shell,
double-shell, Cole), a data-driven initial-guess generator and a robust
multi-stage non-linear least-squares fitter.

<img width="1899" height="954" alt="image" src="https://github.com/user-attachments/assets/a2edc329-2312-4f3c-a29d-e5447e012e24" />


---

## Table of contents

1. [Features](#features)
2. [How to use](#how-to-use)
3. [Input data format](#input-data-format)
4. [Available models](#available-models)
5. [Programmatic use from R](#programmatic-use-from-r)
6. [Repository layout](#repository-layout)
7. [Reproducing the bundled samples](#reproducing-the-bundled-samples)
8. [License](#license)

---

## Features

- Interactive Shiny GUI with Nyquist, Bode-magnitude and Bode-phase tabs.
- Five equivalent-circuit families, switchable from the sidebar.
- Data-driven initial guesses (intercepts, peak frequencies, vacuole time-constant bounds) computed directly from the uploaded spectrum.
- Robust fitter: log-space Levenberg–Marquardt with multi-start, plus a four-stage refinement for the floral circuit.
- Fit summary reports $R^2$, RMSE, derived time constants and per-parameter values with units.
- Bundled synthetic datasets for every model, regenerable from scripts.

## How to use

1. Have R on your machine and install the following dependencies.
    ```r
    install.packages(c("shiny", "ggplot2", "DT", "minpack.lm"))
    ```

2. Clone the repo.
    ```bash
    git clone https://github.com/shinigami-777/FloraFit.git
    cd FloraFit
    ```

3. From the project root run

    ```bash
    Rscript -e 'shiny::runApp(".")'
    ```

The console will print a local URL. Open it
in a browser. The app starts with the bundled `data/sample_eis.csv` already
selected.

## Input data format

Upload a CSV with **exactly three columns** and a header row:

| Column         | Description                                     | Unit   |
|----------------|-------------------------------------------------|--------|
| `frequency_Hz` | Excitation frequency                            | Hz     |
| `Z_real`       | Real part of complex impedance, $\mathrm{Re}(Z)$ | Ω      |
| `Z_imag`       | Imaginary part of complex impedance, $\mathrm{Im}(Z)$ | Ω  |

Rows may be in any order — the app sorts by frequency on load. Use the
convention $Z = Z' + jZ''$ . Capacitive tissues will have negative `Z_imag`.

Example (first lines of `data/sample_eis.csv`):

```csv
frequency_Hz,Z_real,Z_imag
1,114726.73,-209638.76
1.215,86015.47,-172799.92
1.477,88316.13,-143651.90
```

## Available models

All models share the same fitting pipeline; switching between them only changes
the impedance function, parameter set and bounds.

### 1. Floral (default for petal data)

Petal-specific double-shell circuit with a constant-phase element (CPE) for the
cuticle and a parallel vascular shunt:

$$
Z(\omega) = R_S + Z_{\mathrm{CPE}}(\omega) + \bigl[\,Z_{\mathrm{shell}}(\omega) \,\|\, R_{\mathrm{vasc}}\bigr]
$$

where,

$$
Z_{\mathrm{CPE}}(\omega) = \frac{1}{Q_c\,(j\omega)^{\alpha_c}}, \qquad
Z_{\mathrm{shell}}(\omega) = R_E + \frac{1}{j\omega C_M} + R_{\mathrm{CYT}} + \frac{R_V}{1 + j\omega R_V C_T}.
$$

| Symbol           | Meaning                                          | Unit          |
|------------------|--------------------------------------------------|---------------|
| $R_S$          | Series/electrolyte resistance                    | Ω             |
| $Q_c, \alpha_c$| Cuticle CPE magnitude and exponent               | S·sᵅ, –       |
| $R_E$          | Extracellular (apoplast) resistance              | Ω             |
| $C_M$          | Membrane capacitance                             | F             |
| $R_{\mathrm{CYT}}$ | Cytoplasmic resistance                       | Ω             |
| $C_T, R_V$     | Tonoplast capacitance, vacuole resistance        | F, Ω          |
| $R_{\mathrm{vasc}}$ | Vascular shunt resistance                   | Ω             |

Derived: $\tau_M = R_E C_M$, $\tau_V = R_V C_T$.

### 2. Voigt (1–6 RC branches)

$$
Z(\omega) = R_0 + \sum_{k=1}^{N} \frac{R_k}{1 + j\omega R_k C_k}
$$

Branch count $N$ is set in the sidebar.

### 3. Single-shell

$$
Z(\omega) = R_e + \frac{R_i}{1 + j\omega R_i C_m}
$$

### 4. Double-shell

$$
Z(\omega) = R_e + \left[\frac{1}{j\omega C_m} \,\|\, \left(R_i + \frac{R_v}{1 + j\omega R_v C_v}\right)\right]
$$

### 5. Cole

$$
Z(\omega) = R_\infty + \frac{R_0 - R_\infty}{1 + (j\omega\tau)^{\alpha}}
$$

## Programmatic use from R

The fitter can be called directly, without launching Shiny:

```r
source("R/circuit_models.R")
source("R/initial_guess.R")
source("R/fitting.R")

df   <- read.csv("data/sample_eis.csv")
freq <- df$frequency_Hz
Z    <- complex(real = df$Z_real, imaginary = df$Z_imag)

fit <- fit_impedance_model(freq, Z, model_id = "voigt", voigt_n = 2)

fit$theta          # named vector of recovered parameters
fit$r2             # coefficient of determination
fit$rmse           # RMSE in Ohm
fit$Z_fit          # fitted complex impedance at the input frequencies
fit$extra_summary  # derived quantities (time constants, etc.)
```

For the petal-specific floral circuit, use the dedicated multi-stage fitter:

```r
fit <- fit_floral_circuit(freq, Z)
print(fit$theta)
cat(sprintf("tau_M = %.3g s, tau_V = %.3g s\n", fit$tau_M, fit$tau_V))
```

### Expected output on `sample_eis.csv`

```
Converged:  TRUE
Iterations: 0
R^2:        0.99357
RMSE [Ohm]: 4266.3
tau_V [s]:  0.004272   (R_V * C_T)
tau_M [s]:  0.0001787  (R_E * C_M)
LM message: Post-vacuole band 1-D refine (RE, CM, RS, RCYT)
```

## Repository layout

```
Floral_Circuit_Analysis/
├── app.R                       Shiny UI + server
├── R/
│   ├── circuit_models.R        Impedance functions for all 5 circuits
│   ├── initial_guess.R         Data-driven starting values
│   ├── fitting.R               Multi-stage CNLS estimators
│   └── plotting.R              Nyquist and Bode plots
├── scripts/
│   ├── make_sample_data.R      Synthetic floral-circuit spectrum
│   ├── make_model_samples.R    One synthetic spectrum per circuit
│   └── smoke_test.R            End-to-end recovery check
├── data/
│   ├── sample_eis.csv         
│   ├── sample_single_shell.csv
│   ├── sample_double_shell.csv
│   ├── sample_cole.csv
│   ├── sample_voigt_n2.csv
│   └── sample_voigt_n3.csv
├── README.md                  
└── LICENSE                     MIT
```

## Reproducing the bundled samples

All sample CSVs in `data/` are regenerable from known ground-truth parameters:

```bash
Rscript scripts/make_sample_data.R     # regenerate data/sample_eis.csv
Rscript scripts/make_model_samples.R   # regenerate the four model-specific CSVs
Rscript scripts/smoke_test.R           # fit sample_eis.csv and compare to truth
```

`smoke_test.R` prints a side-by-side comparison of true, initial and fitted
parameters, and is the recommended way to verify a fresh install.


## License

Released under the MIT License — see [`LICENSE`](LICENSE) for the full text.
