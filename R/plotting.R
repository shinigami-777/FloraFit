# Nyquist and Bode plots for measured vs fitted EIS spectra.

nyquist_plot <- function(freq, Z_meas, Z_fit = NULL) {
  df <- data.frame(
    Zre = Re(Z_meas),
    negZim = -Im(Z_meas),
    series = "Measured",
    freq = freq
  )
  if (!is.null(Z_fit)) {
    df <- rbind(df, data.frame(
      Zre = Re(Z_fit),
      negZim = -Im(Z_fit),
      series = "Fit",
      freq = freq
    ))
  }
  ggplot2::ggplot(df, ggplot2::aes(Zre, negZim, color = series)) +
    ggplot2::geom_point(data = subset(df, series == "Measured"), size = 2) +
    ggplot2::geom_path(data = subset(df, series == "Fit"), linewidth = 1) +
    ggplot2::scale_color_manual(values = c(Measured = "#1f77b4", Fit = "#d62728")) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = "Re(Z) [Ohm]", y = "-Im(Z) [Ohm]",
                  title = "Nyquist plot", color = NULL) +
    ggplot2::theme_minimal(base_size = 13)
}

bode_plot <- function(freq, Z_meas, Z_fit = NULL) {
  make_df <- function(Z, label) {
    data.frame(
      freq = freq,
      mag  = Mod(Z),
      phase = Arg(Z) * 180 / pi,
      series = label
    )
  }
  df <- make_df(Z_meas, "Measured")
  if (!is.null(Z_fit)) df <- rbind(df, make_df(Z_fit, "Fit"))

  p_mag <- ggplot2::ggplot(df, ggplot2::aes(freq, mag, color = series)) +
    ggplot2::geom_point(data = subset(df, series == "Measured"), size = 1.8) +
    ggplot2::geom_line(data = subset(df, series == "Fit"), linewidth = 1) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::scale_color_manual(values = c(Measured = "#1f77b4", Fit = "#d62728")) +
    ggplot2::labs(x = "Frequency [Hz]", y = "|Z| [Ohm]",
                  title = "Bode magnitude", color = NULL) +
    ggplot2::theme_minimal(base_size = 13)

  p_phase <- ggplot2::ggplot(df, ggplot2::aes(freq, phase, color = series)) +
    ggplot2::geom_point(data = subset(df, series == "Measured"), size = 1.8) +
    ggplot2::geom_line(data = subset(df, series == "Fit"), linewidth = 1) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_color_manual(values = c(Measured = "#1f77b4", Fit = "#d62728")) +
    ggplot2::labs(x = "Frequency [Hz]", y = "Phase [deg]",
                  title = "Bode phase", color = NULL) +
    ggplot2::theme_minimal(base_size = 13)

  list(magnitude = p_mag, phase = p_phase)
}
