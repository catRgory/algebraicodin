# Fit a model to data using MCMC

Convenience wrapper around dust2 particle filtering and monty MCMC
sampling. Handles model compilation, filter creation, packer setup, and
MCMC execution.

## Usage

``` r
fit_mcmc(
  model = NULL,
  data = NULL,
  pars,
  fixed = list(),
  prior = NULL,
  compare = NULL,
  type = "stochastic",
  n_particles = 200L,
  n_steps = 1000L,
  n_chains = 4L,
  vcv = NULL,
  initial = NULL,
  initial_state = NULL,
  save_trajectories = FALSE,
  dt = NULL,
  time_start = 0,
  filter = NULL,
  ...
)
```

## Arguments

- model:

  A LabelledPetriNet, ReactionNet, or StockFlowModel

- data:

  A data.frame with `time` column and observed data columns

- pars:

  Character vector of parameter names to estimate

- fixed:

  Named list of fixed parameter values

- prior:

  A monty_model for the prior (e.g., from monty_dsl or
  [`build_prior()`](https://catrgory.github.io/algebraicodin/reference/build_prior.md))

- compare:

  Named list of observation specifications (see
  [`observe()`](https://catrgory.github.io/algebraicodin/reference/observe.md))

- type:

  Model type: "stochastic" (default) or "discrete"

- n_particles:

  Number of particles (default 200)

- n_steps:

  Number of MCMC steps (default 1000)

- n_chains:

  Number of MCMC chains (default 4)

- vcv:

  Proposal variance-covariance matrix. If NULL, uses diagonal with 0.01
  variance.

- initial:

  Initial parameter values (numeric vector or named list). If NULL, uses
  midpoint of prior domain or 0.1 for each parameter.

- initial_state:

  Optional named list of initial state values

- save_trajectories:

  Logical; save particle trajectories?

- dt:

  Time step for discrete models

- time_start:

  Start time (default: 0)

- filter:

  Pre-built dust2 filter (if provided, model/data/compare are ignored)

- ...:

  Additional arguments passed to
  [`monty::monty_sample()`](https://mrc-ide.github.io/monty/reference/monty_sample.html)

## Value

A monty_samples object with posterior samples
