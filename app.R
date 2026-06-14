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
    Parameter = names(theta),
    Value     = signif(as.numeric(theta), 4),
    Unit      = "",
    stringsAsFactors = FALSE
  )
}

format_model_params <- function(theta, units) {
  data.frame(
    Parameter = names(theta),
    Value     = signif(as.numeric(theta), 4),
    Unit      = units,
    stringsAsFactors = FALSE
  )
}

# ---- UI --------------------------------------------------------------------

app_css <- "
  body { background-color: #f7f8fa; }
  .app-title {
    padding: 14px 20px;
    margin-bottom: 16px;
    background: linear-gradient(135deg, #2c3e50 0%, #4a6785 100%);
    color: #fff;
    border-radius: 6px;
    box-shadow: 0 2px 6px rgba(0,0,0,0.08);
  }
  .app-title h2 { margin: 0; font-weight: 500; }
  .app-title .subtitle { font-size: 13px; opacity: 0.85; margin-top: 4px; }
  .panel-card {
    background: #ffffff;
    padding: 16px 18px;
    border-radius: 6px;
    border: 1px solid #e5e7eb;
    box-shadow: 0 1px 2px rgba(0,0,0,0.03);
    margin-bottom: 16px;
  }
  .panel-card h4 {
    margin-top: 0;
    margin-bottom: 12px;
    padding-bottom: 8px;
    border-bottom: 1px solid #eef0f3;
    font-weight: 500;
    color: #2c3e50;
  }
  .btn-primary { background-color: #2c7be5; border-color: #2c7be5; }
  .btn-primary:hover { background-color: #1a68d1; border-color: #1a68d1; }
  #fit { width: 100%; margin-top: 6px; }
  #reset_guess { width: 100%; margin-top: 6px; }
  #summary {
    background-color: #f8f9fb;
    border: 1px solid #e5e7eb;
    border-radius: 4px;
    font-size: 13px;
    padding: 10px 12px;
    max-height: 220px;
  }
"

ui <- fluidPage(
  tags$head(tags$style(HTML(app_css))),
  div(class = "app-title",
      h2("Floral Circuit Analysis"),
      div(class = "subtitle", "Petal EIS Fitter")),
  fluidRow(
    column(
      width = 3,
      div(class = "panel-card",
        h4("Data"),
        fileInput("file", "Upload EIS CSV",
                  accept = c(".csv", "text/csv")),
        helpText("Columns required: frequency_Hz, Z_real, Z_imag."),
        checkboxInput("use_sample", "Use bundled sample", TRUE)
      ),
      div(class = "panel-card",
        h4("Model"),
        selectInput("model_id", "Fitting model",
                    choices = MODEL_CHOICES, selected = "voigt"),
        conditionalPanel(
          condition = "input.model_id == 'voigt'",
          numericInput("voigt_n", "Voigt RC branches (N)",
                       value = 2, min = 1, max = 6, step = 1)
        )
      ),
      div(class = "panel-card",
        h4("Initial guesses"),
        uiOutput("initial_guess_ui"),
        actionButton("reset_guess", "Reset initial guesses"),
        actionButton("fit", "Fit circuit", class = "btn-primary")
      )
    ),
    column(
      width = 6,
      div(class = "panel-card",
        tabsetPanel(
          tabPanel("Nyquist",    plotOutput("nyquist", height = "520px")),
          tabPanel("Bode |Z|",   plotOutput("bode_mag", height = "320px")),
          tabPanel("Bode phase", plotOutput("bode_phase", height = "320px")),
          tabPanel("Raw data",   DT::DTOutput("raw"))
        )
      )
    ),
    column(
      width = 3,
      div(class = "panel-card",
        h4("Fit summary"),
        verbatimTextOutput("summary")
      ),
      div(class = "panel-card",
        h4("Fitted parameters"),
        DT::DTOutput("params")
      )
    )
  )
)

# ---- server ----------------------------------------------------------------

server <- function(input, output, session) {

  model_spec <- reactive({
    get_model_spec(input$model_id, voigt_n = input$voigt_n %||% 2L)
  })

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

  default_guess <- reactive({
    d <- data_in()
    req(d)
    spec <- model_spec()
    default_initial_guess_by_model(
      d$freq, d$Z,
      model_id = spec$id,
      voigt_n = input$voigt_n %||% 2L
    )
  })

  guess_values <- reactiveVal(NULL)

  observeEvent(default_guess(), {
    guess_values(default_guess())
  }, ignoreInit = FALSE)

  observeEvent(input$reset_guess, {
    guess_values(default_guess())
  })

  output$initial_guess_ui <- renderUI({
    spec <- model_spec()
    dflt <- guess_values()
    if (is.null(dflt)) dflt <- default_guess()
    dflt <- dflt[spec$param_names]

    controls <- lapply(spec$param_names, function(p) {
      val <- as.numeric(dflt[[p]])
      lower <- as.numeric(spec$lower[[p]])
      upper <- as.numeric(spec$upper[[p]])
      step <- if (val >= 1) 1 else 10^floor(log10(max(val, 1e-12)))
      numericInput(
        inputId = paste0("guess_", p),
        label = sprintf("%s [%s]", p, spec$param_units[match(p, spec$param_names)]),
        value = signif(val, 6),
        min = lower,
        max = upper,
        step = step
      )
    })
    do.call(tagList, controls)
  })

  user_guess <- reactive({
    spec <- model_spec()
    dflt <- default_guess()
    vals <- vapply(spec$param_names, function(p) {
      v <- input[[paste0("guess_", p)]]
      if (is.null(v) || !is.finite(v)) dflt[[p]] else as.numeric(v)
    }, numeric(1))
    names(vals) <- spec$param_names
    vals
  })

  fit_result <- eventReactive(input$fit, {
    d <- data_in()
    req(d)
    spec <- model_spec()
    theta_init <- user_guess()
    withProgress(message = "Fitting circuit...", value = 0.3, {
      fit_impedance_model(
        d$freq, d$Z,
        model_id = spec$id,
        theta_init = theta_init,
        voigt_n = input$voigt_n %||% 2L
      )
    })
  }, ignoreNULL = FALSE)

  output$summary <- renderPrint({
    fit <- fit_result()
    if (is.null(fit)) {
      cat("Press 'Fit circuit' after loading data.\n")
      return(invisible())
    }
    cat("Model:      ", fit$model_label, "\n", sep = "")
    cat("Converged: ", fit$converged, "\n", sep = "")
    cat("Iterations: ", fit$niter, "\n", sep = "")
    cat("R^2:        ", signif(fit$r2, 5), "\n", sep = "")
    cat("RMSE [Ohm]: ", signif(fit$rmse, 5), "\n", sep = "")
    if (length(fit$extra_summary) > 0) {
      cat(paste0(fit$extra_summary, collapse = "\n"), "\n", sep = "")
    }
  })

  output$params <- DT::renderDT({
    fit <- fit_result()
    if (is.null(fit)) return(NULL)
    DT::datatable(format_model_params(fit$theta, fit$param_units),
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
