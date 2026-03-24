# Create a dust2 particle filter from a model

Compiles a Petri net or stock-flow model to odin2 with data comparison,
then wraps it in a dust2 particle filter (or deterministic unfilter).

## Usage

``` r
create_filter(
  model,
  data,
  compare = NULL,
  type = "stochastic",
  n_particles = 200L,
  initial = NULL,
  dt = NULL,
  time_start = 0,
  seed = NULL,
  ...
)
```

## Arguments

- model:

  A LabelledPetriNet, ReactionNet, or StockFlowModel

- data:

  A data.frame with a `time` column and observed data columns

- compare:

  Named list of observation specifications (see
  [`observe()`](https://catrgory.github.io/algebraicodin/reference/observe.md))

- type:

  Model type: "stochastic" (default), "discrete", or "ode"

- n_particles:

  Number of particles (default 200, ignored for ODE)

- initial:

  Optional named list of initial state values

- dt:

  Time step for discrete models (passed to dust2)

- time_start:

  Start time (default: 0)

- seed:

  Random seed (optional)

- ...:

  Additional arguments passed to dust2 filter creation

## Value

A dust2 filter object (dust_likelihood)
