# Plotting utilities --------------------------------------------------------
# ggplot2-based trajectory and structure visualization.
# DiagrammeR is used for graph/network diagrams (see catlab/R/graphviz.R).

#' Plot a Petri net as a bipartite graph via DiagrammeR
#'
#' Species are shown as blue circles, transitions as orange boxes.
#' Accepts LabelledPetriNet, TypedPetriNet, or Open Petri nets.
#'
#' @param pn A Petri net ACSet, TypedPetriNet, or Open Petri net
#' @returns A DiagrammeR `htmlwidget`
#' @examples
#' sir <- labelled_petri_net(
#'   c("S", "I", "R"),
#'   inf = c("S", "I") %=>% c("I", "I"),
#'   rec = "I" %=>% "R"
#' )
#' if (requireNamespace("DiagrammeR", quietly = TRUE)) {
#'   plot_petri(sir)
#' }
#' @export
plot_petri <- function(pn) {
  if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
    cli::cli_abort(c(
      "The {.pkg DiagrammeR} package is required for rendering.",
      i = "Install it with {.code install.packages('DiagrammeR')}."
    ))
  }
  if (S7::S7_inherits(pn, TypedPetriNet)) pn <- pn@pn
  if (S7::S7_inherits(pn, Open)) pn <- pn@apex
  DiagrammeR::grViz(catlab::petri_to_dot(pn))
}

#' Plot an undirected wiring diagram via DiagrammeR
#'
#' Boxes are yellow rectangles, junctions are grey circles,
#' outer ports are small black diamonds.
#'
#' @param w A UWD ACSet
#' @returns A DiagrammeR `htmlwidget`
#' @examples
#' w <- uwd(
#'   outer = c("s", "i", "r"),
#'   infection = c("s", "i"),
#'   recovery = c("i", "r")
#' )
#' if (requireNamespace("DiagrammeR", quietly = TRUE)) {
#'   plot_uwd(w)
#' }
#' @export
plot_uwd <- function(w) {
  if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
    cli::cli_abort(c(
      "The {.pkg DiagrammeR} package is required for rendering.",
      i = "Install it with {.code install.packages('DiagrammeR')}."
    ))
  }
  DiagrammeR::grViz(catlab::uwd_to_dot(w))
}

#' Plot ODE solution trajectories with ggplot2
#'
#' @param sol A deSolve solution matrix or data.frame with a `time` column
#' @param vars Character vector of variable names to plot (default: all non-time)
#' @param title Optional plot title
#' @returns A ggplot2 object
#' @examples
#' sol <- data.frame(time = seq(0, 10, 0.1),
#'                   S = 990 * exp(-0.3 * seq(0, 10, 0.1)),
#'                   I = 10 * exp(0.1 * seq(0, 10, 0.1)))
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_trajectory(sol, title = "SIR Trajectory")
#' }
#' @export
plot_trajectory <- function(sol, vars = NULL, title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort(c(
      "The {.pkg ggplot2} package is required for plotting.",
      i = "Install it with {.code install.packages('ggplot2')}."
    ))
  }

  # Convert deSolve matrix to data.frame if needed
  if (is.matrix(sol)) {
    sol <- as.data.frame(sol)
  }

  # Find time column
  time_col <- intersect(c("time", "Time", "t"), names(sol))[1]
  if (is.na(time_col)) {
    cli::cli_abort("No time column found. Expected 'time', 'Time', or 't'.")
  }

  if (is.null(vars)) {
    vars <- setdiff(names(sol), time_col)
  }

  # Reshape to long format
  long_df <- data.frame(
    time = rep(sol[[time_col]], length(vars)),
    variable = rep(vars, each = nrow(sol)),
    value = unlist(lapply(vars, function(v) sol[[v]]))
  )
  long_df$variable <- factor(long_df$variable, levels = vars)

  p <- ggplot2::ggplot(long_df,
    ggplot2::aes(x = .data$time, y = .data$value, color = .data$variable)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::labs(x = "Time", y = "Count", color = "Compartment") +
    ggplot2::theme_minimal()

  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }
  p
}

#' Plot a comparison of two ODE solutions (e.g., R vs Julia)
#'
#' @param sol1 First solution (data.frame with time column)
#' @param sol2 Second solution (data.frame with time column)
#' @param vars Character vector of variable names to compare
#' @param labels Two-element character vector labeling each solution
#' @param title Optional plot title
#' @returns A ggplot2 object
#' @examples
#' t <- seq(0, 10, 0.1)
#' sol1 <- data.frame(time = t, S = 990 * exp(-0.3 * t), I = 10 * exp(0.1 * t))
#' sol2 <- data.frame(time = t, S = 980 * exp(-0.3 * t), I = 20 * exp(0.1 * t))
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_comparison(sol1, sol2, labels = c("Run 1", "Run 2"))
#' }
#' @export
plot_comparison <- function(sol1, sol2, vars = NULL,
                            labels = c("Model 1", "Model 2"),
                            title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort(c(
      "The {.pkg ggplot2} package is required for plotting.",
      i = "Install it with {.code install.packages('ggplot2')}."
    ))
  }

  time_col1 <- intersect(c("time", "Time", "t"), names(sol1))[1]
  time_col2 <- intersect(c("time", "Time", "t"), names(sol2))[1]

  if (is.null(vars)) {
    vars <- intersect(
      setdiff(names(sol1), time_col1),
      setdiff(names(sol2), time_col2)
    )
  }

  # Build long-format data from both solutions
  make_long <- function(sol, time_col, label) {
    data.frame(
      time = rep(sol[[time_col]], length(vars)),
      variable = rep(vars, each = nrow(sol)),
      value = unlist(lapply(vars, function(v) sol[[v]])),
      source = label
    )
  }

  long_df <- rbind(
    make_long(sol1, time_col1, labels[1]),
    make_long(sol2, time_col2, labels[2])
  )
  long_df$variable <- factor(long_df$variable, levels = vars)

  p <- ggplot2::ggplot(long_df,
    ggplot2::aes(x = .data$time, y = .data$value,
                 color = .data$variable, linetype = .data$source)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::labs(x = "Time", y = "Count",
                  color = "Compartment", linetype = "Source") +
    ggplot2::theme_minimal()

  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }
  p
}

#' Plot absolute differences between two ODE solutions
#'
#' @param sol1 First solution (data.frame)
#' @param sol2 Second solution (data.frame)
#' @param vars Character vector of variable names
#' @param title Optional plot title
#' @returns A ggplot2 object
#' @examples
#' t <- seq(0, 10, 0.1)
#' sol1 <- data.frame(time = t, S = 990 * exp(-0.3 * t), I = 10 * exp(0.1 * t))
#' sol2 <- data.frame(time = t, S = 989 * exp(-0.3 * t), I = 11 * exp(0.1 * t))
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_diff(sol1, sol2)
#' }
#' @export
plot_diff <- function(sol1, sol2, vars = NULL, title = "Absolute Difference") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort(c(
      "The {.pkg ggplot2} package is required for plotting.",
      i = "Install it with {.code install.packages('ggplot2')}."
    ))
  }

  time_col1 <- intersect(c("time", "Time", "t"), names(sol1))[1]
  time_col2 <- intersect(c("time", "Time", "t"), names(sol2))[1]

  if (is.null(vars)) {
    vars <- intersect(
      setdiff(names(sol1), time_col1),
      setdiff(names(sol2), time_col2)
    )
  }

  # Interpolate to common time grid
  common_times <- intersect(
    round(sol1[[time_col1]], 4),
    round(sol2[[time_col2]], 4)
  )
  idx1 <- match(common_times, round(sol1[[time_col1]], 4))
  idx2 <- match(common_times, round(sol2[[time_col2]], 4))

  diff_df <- data.frame(
    time = rep(common_times, length(vars)),
    variable = rep(vars, each = length(common_times)),
    value = unlist(lapply(vars, function(v) {
      abs(sol1[[v]][idx1] - sol2[[v]][idx2])
    }))
  )
  diff_df$variable <- factor(diff_df$variable, levels = vars)

  p <- ggplot2::ggplot(diff_df,
    ggplot2::aes(x = .data$time, y = .data$value, color = .data$variable)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Time", y = "|Difference|", color = "Variable") +
    ggplot2::theme_minimal()

  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }
  p
}
