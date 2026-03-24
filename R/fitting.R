# Data fitting with dust2 particle filters and monty MCMC ----------------

#' Create an observation specification
#'
#' Specifies how a model variable relates to observed data for use in
#' data comparison (likelihood computation). Used with the `compare`
#' argument of code generation functions.
#'
#' @param dist Distribution name (e.g., "Poisson", "Normal", "NegativeBinomial")
#' @param transition Name of a Petri net transition to track incidence
#' @param flow Name of a stock-flow flow to track incidence
#' @param state Name of a state variable to observe directly (prevalence)
#' @param expr Custom expression string for the distribution parameter
#' @param noise Noise parameter name (added to mean to prevent zero lambda)
#' @param sd Standard deviation parameter name (for Normal distribution)
#' @param zero_every Reset interval for incidence accumulators (default: 1)
#' @returns An observe_spec object
#' @export
observe <- function(dist, transition = NULL, flow = NULL, state = NULL,
                    expr = NULL, noise = NULL, sd = NULL,
                    zero_every = 1) {
  sources <- c(!is.null(transition), !is.null(flow),
               !is.null(state), !is.null(expr))
  if (sum(sources) != 1) {
    cli::cli_abort(paste(
      "Exactly one of {.arg transition}, {.arg flow},",
      "{.arg state}, or {.arg expr} must be specified"))
  }
  structure(list(
    dist = dist,
    transition = transition,
    flow = flow,
    state = state,
    expr = expr,
    noise = noise,
    sd = sd,
    zero_every = zero_every
  ), class = "observe_spec")
}

# Internal: generate odin2 comparison code from compare specs
generate_compare_code <- function(compare, model_type = c("petri", "stockflow")) {
  if (is.null(compare)) return(character())
  model_type <- match.arg(model_type)

  lines <- character()
  incidence_lines <- character()
  param_lines <- character()
  noise_lines <- character()

  for (obs_var in names(compare)) {
    spec <- compare[[obs_var]]

    # Normalize old simple format: list(cases = "Poisson") → observe_spec
    if (is.character(spec)) {
      spec <- observe(spec, expr = obs_var)
    }
    if (!inherits(spec, "observe_spec")) {
      cli::cli_abort(
        "compare[[\"{obs_var}\"]] must be a string or {.fn observe} spec")
    }

    # Determine model expression for comparison
    if (!is.null(spec$transition)) {
      inc_var <- paste0("incidence_", spec$transition)
      incidence_lines <- c(incidence_lines,
        sprintf("initial(%s, zero_every = %d) <- 0",
                inc_var, spec$zero_every),
        sprintf("update(%s) <- %s + n_%s",
                inc_var, inc_var, spec$transition))
      model_expr <- inc_var

    } else if (!is.null(spec$flow)) {
      inc_var <- paste0("incidence_", spec$flow)
      incidence_lines <- c(incidence_lines,
        sprintf("initial(%s, zero_every = %d) <- 0",
                inc_var, spec$zero_every),
        sprintf("update(%s) <- %s + n_%s",
                inc_var, inc_var, spec$flow))
      model_expr <- inc_var

    } else if (!is.null(spec$state)) {
      model_expr <- spec$state

    } else {
      model_expr <- spec$expr
    }

    # Handle noise parameter
    if (!is.null(spec$noise)) {
      param_lines <- c(param_lines,
        sprintf("%s <- parameter(1e6)", spec$noise))
      noise_var <- paste0("noise_", obs_var)
      noise_lines <- c(noise_lines,
        sprintf("%s <- %s + 1 / %s", noise_var, model_expr, spec$noise))
      model_expr <- noise_var
    }

    # Build comparison statement
    if (spec$dist == "Normal" && !is.null(spec$sd)) {
      param_lines <- c(param_lines,
        sprintf("%s <- parameter()", spec$sd))
      lines <- c(lines,
        sprintf("%s <- data()", obs_var),
        sprintf("%s ~ Normal(%s, %s)", obs_var, model_expr, spec$sd))
    } else {
      lines <- c(lines,
        sprintf("%s <- data()", obs_var),
        sprintf("%s ~ %s(%s)", obs_var, spec$dist, model_expr))
    }
  }

  result <- character()
  if (length(param_lines) > 0)
    result <- c(result, "", "## Observation parameters", unique(param_lines))
  if (length(incidence_lines) > 0)
    result <- c(result, "", "## Incidence tracking", incidence_lines)
  if (length(noise_lines) > 0)
    result <- c(result, "", "## Observation noise", noise_lines)
  if (length(lines) > 0)
    result <- c(result, "", "## Comparison to data", lines)

  result
}


#' Create a dust2 particle filter from a model
#'
#' Compiles a Petri net or stock-flow model to odin2 with data comparison,
#' then wraps it in a dust2 particle filter (or deterministic unfilter).
#'
#' @param model A LabelledPetriNet, ReactionNet, or StockFlowModel
#' @param data A data.frame with a `time` column and observed data columns
#' @param compare Named list of observation specifications (see [observe()])
#' @param type Model type: "stochastic" (default), "discrete", or "ode"
#' @param n_particles Number of particles (default 200, ignored for ODE)
#' @param initial Optional named list of initial state values
#' @param dt Time step for discrete models (passed to dust2)
#' @param time_start Start time (default: 0)
#' @param seed Random seed (optional)
#' @param ... Additional arguments passed to dust2 filter creation
#' @returns A dust2 filter object (dust_likelihood)
#' @export
create_filter <- function(model, data, compare = NULL,
                          type = "stochastic",
                          n_particles = 200L, initial = NULL,
                          dt = NULL, time_start = 0, seed = NULL,
                          ...) {
  if (!requireNamespace("odin2", quietly = TRUE))
    cli::cli_abort("Package {.pkg odin2} required for filtering")
  if (!requireNamespace("dust2", quietly = TRUE))
    cli::cli_abort("Package {.pkg dust2} required for filtering")

  # Generate odin code with comparison
  if (S7::S7_inherits(model, StockFlowModel)) {
    code <- sf_to_odin(model, type = type, initial = initial,
                       compare = compare)
  } else {
    code <- pn_to_odin(model, type = type, compare = compare)
  }

  gen <- odin2::odin(code)

  # ODE models use deterministic unfilter; discrete use particle filter
  if (type == "ode") {
    dust2::dust_unfilter_create(
      gen(), time_start = time_start, data = data,
      dt = dt, ...)
  } else {
    args <- list(generator = gen(), time_start = time_start,
                 data = data, n_particles = n_particles)
    if (!is.null(dt)) args$dt <- dt
    if (!is.null(seed)) args$seed <- seed
    do.call(dust2::dust_filter_create, c(args, list(...)))
  }
}


#' Fit a model to data using MCMC
#'
#' Convenience wrapper around dust2 particle filtering and monty MCMC
#' sampling. Handles model compilation, filter creation, packer setup,
#' and MCMC execution.
#'
#' @param model A LabelledPetriNet, ReactionNet, or StockFlowModel
#' @param data A data.frame with `time` column and observed data columns
#' @param pars Character vector of parameter names to estimate
#' @param fixed Named list of fixed parameter values
#' @param prior A monty_model for the prior (e.g., from monty_dsl or
#'   [build_prior()])
#' @param compare Named list of observation specifications (see [observe()])
#' @param type Model type: "stochastic" (default) or "discrete"
#' @param n_particles Number of particles (default 200)
#' @param n_steps Number of MCMC steps (default 1000)
#' @param n_chains Number of MCMC chains (default 4)
#' @param vcv Proposal variance-covariance matrix. If NULL, uses
#'   diagonal with 0.01 variance.
#' @param initial Initial parameter values (numeric vector or named list).
#'   If NULL, uses midpoint of prior domain or 0.1 for each parameter.
#' @param initial_state Optional named list of initial state values
#' @param save_trajectories Logical; save particle trajectories?
#' @param dt Time step for discrete models
#' @param time_start Start time (default: 0)
#' @param filter Pre-built dust2 filter (if provided, model/data/compare
#'   are ignored)
#' @param ... Additional arguments passed to [monty::monty_sample()]
#' @returns A monty_samples object with posterior samples
#' @export
fit_mcmc <- function(model = NULL, data = NULL, pars, fixed = list(),
                     prior = NULL, compare = NULL,
                     type = "stochastic",
                     n_particles = 200L,
                     n_steps = 1000L, n_chains = 4L,
                     vcv = NULL, initial = NULL,
                     initial_state = NULL,
                     save_trajectories = FALSE,
                     dt = NULL, time_start = 0,
                     filter = NULL, ...) {
  if (!requireNamespace("monty", quietly = TRUE))
    cli::cli_abort("Package {.pkg monty} required for MCMC fitting")
  if (!requireNamespace("dust2", quietly = TRUE))
    cli::cli_abort("Package {.pkg dust2} required for filtering")

  # Build or use provided filter
  if (is.null(filter)) {
    if (is.null(model) || is.null(data))
      cli::cli_abort("Must provide either {.arg filter} or both {.arg model} and {.arg data}")
    filter <- create_filter(model, data, compare = compare,
                            type = type, n_particles = n_particles,
                            initial = initial_state, dt = dt,
                            time_start = time_start)
  }

  # Create packer
  packer <- monty::monty_packer(pars, fixed = fixed)

  # Create likelihood model
  likelihood <- dust2::dust_likelihood_monty(
    filter, packer,
    save_trajectories = save_trajectories)

  # Combine with prior
  posterior <- if (!is.null(prior)) likelihood + prior else likelihood

  # Proposal variance-covariance
  if (is.null(vcv)) {
    vcv <- diag(length(pars)) * 0.01
  }
  sampler <- monty::monty_sampler_random_walk(vcv)

  # Initial parameter values
  if (is.null(initial)) {
    initial <- rep(0.1, length(pars))
  } else if (is.list(initial)) {
    initial <- unlist(initial[pars])
  }

  monty::monty_sample(posterior, sampler, n_steps,
                      initial = initial, n_chains = n_chains, ...)
}


#' Build a prior model from specifications
#'
#' Constructs a monty_model for use as a prior in [fit_mcmc()].
#' Each argument is a named prior specification created by
#' [prior_exp()], [prior_uniform()], or [prior_normal()].
#'
#' @param ... Named prior specifications
#' @returns A monty_model object
#' @export
#' @examples
#' \dontrun{
#' p <- build_prior(
#'   beta = prior_exp(mean = 0.3),
#'   gamma = prior_exp(mean = 0.1)
#' )
#' }
build_prior <- function(...) {
  if (!requireNamespace("monty", quietly = TRUE))
    cli::cli_abort("Package {.pkg monty} required for prior specification")
  specs <- list(...)
  param_names <- names(specs)
  if (is.null(param_names) || any(param_names == ""))
    cli::cli_abort("All prior specifications must be named")

  # Build density function
  density <- function(x) {
    ll <- 0
    for (i in seq_along(specs)) {
      ll <- ll + specs[[i]]$log_density(x[i])
    }
    ll
  }

  # Build domain matrix
  domain_mat <- do.call(rbind, lapply(specs, function(s) s$domain))
  rownames(domain_mat) <- param_names

  monty::monty_model(list(
    parameters = param_names,
    density = density,
    domain = domain_mat
  ))
}


#' Exponential prior
#' @param mean Mean of the exponential distribution (default 1)
#' @returns A prior specification for [build_prior()]
#' @export
prior_exp <- function(mean = 1) {
  rate <- 1 / mean
  list(
    log_density = function(x) {
      if (x <= 0) return(-Inf)
      stats::dexp(x, rate = rate, log = TRUE)
    },
    domain = c(0, Inf),
    sample = function() stats::rexp(1, rate = rate)
  )
}


#' Uniform prior
#' @param min Lower bound (default 0)
#' @param max Upper bound (default 1)
#' @returns A prior specification for [build_prior()]
#' @export
prior_uniform <- function(min = 0, max = 1) {
  list(
    log_density = function(x) stats::dunif(x, min = min, max = max, log = TRUE),
    domain = c(min, max),
    sample = function() stats::runif(1, min = min, max = max)
  )
}


#' Normal prior
#' @param mean Mean (default 0)
#' @param sd Standard deviation (default 1)
#' @returns A prior specification for [build_prior()]
#' @export
prior_normal <- function(mean = 0, sd = 1) {
  list(
    log_density = function(x) stats::dnorm(x, mean = mean, sd = sd, log = TRUE),
    domain = c(-Inf, Inf),
    sample = function() stats::rnorm(1, mean = mean, sd = sd)
  )
}


#' Log-normal prior
#' @param meanlog Mean of the log (default 0)
#' @param sdlog Standard deviation of the log (default 1)
#' @returns A prior specification for [build_prior()]
#' @export
prior_lognormal <- function(meanlog = 0, sdlog = 1) {
  list(
    log_density = function(x) {
      if (x <= 0) return(-Inf)
      stats::dlnorm(x, meanlog = meanlog, sdlog = sdlog, log = TRUE)
    },
    domain = c(0, Inf),
    sample = function() stats::rlnorm(1, meanlog = meanlog, sdlog = sdlog)
  )
}


#' Plot MCMC trace plots
#'
#' @param samples A monty_samples object from [fit_mcmc()]
#' @param pars Optional character vector of parameter names to plot
#' @param burnin Number of initial samples to shade as burn-in
#' @returns A ggplot object
#' @export
plot_traces <- function(samples, pars = NULL, burnin = 0) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    cli::cli_abort("Package {.pkg ggplot2} required for plotting")

  # Extract parameter names and samples
  if (is.null(pars)) {
    pars <- rownames(samples$pars)
    if (is.null(pars)) pars <- paste0("p", seq_len(nrow(samples$pars)))
  }
  n_chains <- dim(samples$pars)[3]
  n_steps <- dim(samples$pars)[2]

  df_list <- list()
  for (i in seq_along(pars)) {
    for (ch in seq_len(n_chains)) {
      df_list[[length(df_list) + 1]] <- data.frame(
        step = seq_len(n_steps),
        value = samples$pars[i, , ch],
        parameter = pars[i],
        chain = factor(ch)
      )
    }
  }
  df <- do.call(rbind, df_list)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$step, y = .data$value,
                                         colour = .data$chain)) +
    ggplot2::geom_line(alpha = 0.6, linewidth = 0.3) +
    ggplot2::facet_wrap(~ parameter, scales = "free_y") +
    ggplot2::labs(title = "MCMC trace plots", x = "Step", y = "Value") +
    ggplot2::theme_minimal()

  if (burnin > 0) {
    p <- p + ggplot2::geom_vline(xintercept = burnin,
                                  linetype = "dashed", alpha = 0.5)
  }
  p
}


#' Plot posterior density estimates
#'
#' @param samples A monty_samples object from [fit_mcmc()]
#' @param pars Optional character vector of parameter names to plot
#' @param burnin Number of initial samples to discard
#' @param true_values Optional named list of true parameter values
#'   (shown as vertical lines)
#' @returns A ggplot object
#' @export
plot_posterior <- function(samples, pars = NULL, burnin = 0,
                           true_values = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    cli::cli_abort("Package {.pkg ggplot2} required for plotting")

  if (is.null(pars)) {
    pars <- rownames(samples$pars)
    if (is.null(pars)) pars <- paste0("p", seq_len(nrow(samples$pars)))
  }
  n_chains <- dim(samples$pars)[3]
  n_steps <- dim(samples$pars)[2]
  start <- burnin + 1

  df_list <- list()
  for (i in seq_along(pars)) {
    vals <- as.vector(samples$pars[i, start:n_steps, ])
    df_list[[i]] <- data.frame(
      value = vals,
      parameter = pars[i]
    )
  }
  df <- do.call(rbind, df_list)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_density(fill = "steelblue", alpha = 0.5) +
    ggplot2::facet_wrap(~ parameter, scales = "free") +
    ggplot2::labs(title = "Posterior distributions", x = "Value",
                  y = "Density") +
    ggplot2::theme_minimal()

  if (!is.null(true_values)) {
    tv_df <- data.frame(
      parameter = names(true_values),
      value = unlist(true_values)
    )
    p <- p + ggplot2::geom_vline(data = tv_df,
                                  ggplot2::aes(xintercept = .data$value),
                                  colour = "red", linetype = "dashed",
                                  linewidth = 0.8)
  }
  p
}


#' Plot model fit against data
#'
#' Plots observed data with model predictions from MCMC posterior.
#' Requires that `save_trajectories = TRUE` was used in [fit_mcmc()].
#'
#' @param samples A monty_samples object with saved trajectories
#' @param data The observed data (data.frame with time and data columns)
#' @param obs_var Name of the observed variable column in data
#' @param state_var Name of the state variable to plot from trajectories
#' @param burnin Number of initial samples to discard
#' @param quantiles Quantile levels for credible intervals (default:
#'   50% and 95%)
#' @returns A ggplot object
#' @export
plot_fit <- function(samples, data, obs_var, state_var = NULL,
                     burnin = 0, quantiles = c(0.025, 0.25, 0.75, 0.975)) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    cli::cli_abort("Package {.pkg ggplot2} required for plotting")

  if (is.null(samples$observations) ||
      is.null(samples$observations$trajectories)) {
    cli::cli_abort(
      "No trajectories saved. Use {.code save_trajectories = TRUE} in {.fn fit_mcmc}")
  }

  traj <- samples$observations$trajectories
  # traj dimensions: [state × time × steps × chains]
  n_states <- dim(traj)[1]
  n_times <- dim(traj)[2]
  n_steps <- dim(traj)[3]
  n_chains <- dim(traj)[4]

  # Find state index
  state_idx <- if (!is.null(state_var)) {
    state_names <- dimnames(traj)[[1]]
    if (!is.null(state_names)) {
      which(state_names == state_var)
    } else 1
  } else 1

  start <- burnin + 1
  times <- data$time

  # Extract all trajectory values for the state at each time
  vals_mat <- matrix(NA, nrow = n_times,
                     ncol = (n_steps - start + 1) * n_chains)
  col <- 1
  for (ch in seq_len(n_chains)) {
    for (s in start:n_steps) {
      vals_mat[, col] <- traj[state_idx, , s, ch]
      col <- col + 1
    }
  }

  # Compute quantiles
  q_df <- data.frame(
    time = times,
    median = apply(vals_mat, 1, stats::median, na.rm = TRUE)
  )
  for (q in quantiles) {
    q_df[[paste0("q", q)]] <- apply(vals_mat, 1, stats::quantile,
                                     probs = q, na.rm = TRUE)
  }

  p <- ggplot2::ggplot(q_df, ggplot2::aes(x = .data$time)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Model fit", x = "Time", y = "Value")

  # Add credible intervals (outer first)
  if (length(quantiles) >= 4) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data[[paste0("q", quantiles[1])]],
                   ymax = .data[[paste0("q", quantiles[4])]]),
      fill = "steelblue", alpha = 0.2)
  }
  if (length(quantiles) >= 4) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data[[paste0("q", quantiles[2])]],
                   ymax = .data[[paste0("q", quantiles[3])]]),
      fill = "steelblue", alpha = 0.4)
  }

  # Median line
  p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$median),
                               colour = "steelblue", linewidth = 0.8)

  # Observed data
  p <- p + ggplot2::geom_point(
    data = data,
    ggplot2::aes(x = .data$time, y = .data[[obs_var]]),
    colour = "black", size = 1.5)

  p
}
