# Floral Circuit Analysis - Shiny app
#
# Upload an EIS spectrum (frequency, Re(Z), Im(Z)), fit the floral
# equivalent circuit from Oduniyi (2026), and visualise the result.
#
# Run from the project root:    Rscript -e 'shiny::runApp(".")'

library(shiny)
library(ggplot2)

source("R/circuit_models.R")
source("R/initial_guess.R")
source("R/fitting.R")
source("R/plotting.R")

# ---- helpers ---------------------------------------------------------------

read_eis_csv <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  needed <- c("frequency_Hz", "Z_real", "Z_imag")
  if (!all(needed %in% names(df))) {
    stop("CSV must have columns: frequency_Hz, Z_real, Z_imag")
  }
  df <- df[order(df$frequency_Hz), ]
  list(
    freq = df$frequency_Hz,
    Z    = complex(real = df$Z_real, imaginary = df$Z_imag)
  )
}

format_params <- function(theta) {
  data.frame(
    Parameter = PARAM_NAMES,
    Value     = signif(theta, 4),
    Unit      = PARAM_UNITS,
    stringsAsFactors = FALSE
  )
}

# ---- UI --------------------------------------------------------------------

ui <- fluidPage(
  titlePanel("Floral Circuit Analysis - Petal EIS Fitter"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      fileInput("file", "Upload EIS CSV",
                accept = c(".csv", "text/csv")),
      helpText("CSV columns required: frequency_Hz, Z_real, Z_imag."),
      checkboxInput("use_sample", "Use bundled sample (data/sample_eis.csv)", TRUE),
      actionButton("fit", "Fit floral circuit", class = "btn-primary"),
      hr(),
      h4("Fit summary"),
      verbatimTextOutput("summary"),
      hr(),
      h4("Fitted parameters"),
      DT::DTOutput("params")
    ),
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Nyquist",  plotOutput("nyquist", height = "520px")),
        tabPanel("Bode |Z|", plotOutput("bode_mag", height = "320px")),
        tabPanel("Bode phase", plotOutput("bode_phase", height = "320px")),
        tabPanel("Raw data", DT::DTOutput("raw"))
      )
    )
  )
)

# ---- server ----------------------------------------------------------------

server <- function(input, output, session) {

  data_in <- reactive({
    if (isTRUE(input$use_sample) && is.null(input$file)) {
      if (!file.exists("data/sample_eis.csv")) {
        showNotification(
          "Sample file missing - run scripts/make_sample_data.R first.",
          type = "error")
        return(NULL)
      }
      return(read_eis_csv("data/sample_eis.csv"))
    }
    req(input$file)
    tryCatch(read_eis_csv(input$file$datapath),
             error = function(e) {
               showNotification(conditionMessage(e), type = "error")
               NULL
             })
  })

  fit_result <- eventReactive(input$fit, {
    d <- data_in()
    req(d)
    withProgress(message = "Fitting circuit...", value = 0.3, {
      fit_floral_circuit(d$freq, d$Z)
    })
  }, ignoreNULL = FALSE)

  output$summary <- renderPrint({
    fit <- fit_result()
    if (is.null(fit)) {
      cat("Press 'Fit floral circuit' after loading data.\n")
      return(invisible())
    }
    cat("Converged: ", fit$converged, "\n", sep = "")
    cat("Iterations: ", fit$niter, "\n", sep = "")
    cat("R^2:        ", signif(fit$r2, 5), "\n", sep = "")
    cat("RMSE [Ohm]: ", signif(fit$rmse, 5), "\n", sep = "")
    cat("tau_V [s]:  ", signif(fit$tau_V, 4),
        "  (R_V * C_T)\n", sep = "")
    cat("tau_M [s]:  ", signif(fit$tau_M, 4),
        "  (R_E * C_M)\n", sep = "")
    cat("LM message: ", fit$message, "\n", sep = "")
  })

  output$params <- DT::renderDT({
    fit <- fit_result()
    if (is.null(fit)) return(NULL)
    DT::datatable(format_params(fit$theta),
                  options = list(dom = "t", paging = FALSE),
                  rownames = FALSE)
  })

  output$nyquist <- renderPlot({
    d <- data_in(); req(d)
    fit <- tryCatch(fit_result(), error = function(e) NULL)
    nyquist_plot(d$freq, d$Z, if (!is.null(fit)) fit$Z_fit else NULL)
  })

  output$bode_mag <- renderPlot({
    d <- data_in(); req(d)
    fit <- tryCatch(fit_result(), error = function(e) NULL)
    bode_plot(d$freq, d$Z, if (!is.null(fit)) fit$Z_fit else NULL)$magnitude
  })

  output$bode_phase <- renderPlot({
    d <- data_in(); req(d)
    fit <- tryCatch(fit_result(), error = function(e) NULL)
    bode_plot(d$freq, d$Z, if (!is.null(fit)) fit$Z_fit else NULL)$phase
  })

  output$raw <- DT::renderDT({
    d <- data_in(); req(d)
    DT::datatable(data.frame(
      frequency_Hz = d$freq,
      Z_real       = Re(d$Z),
      Z_imag       = Im(d$Z),
      mag          = Mod(d$Z),
      phase_deg    = Arg(d$Z) * 180 / pi
    ), options = list(pageLength = 15))
  })
}

shinyApp(ui, server)
