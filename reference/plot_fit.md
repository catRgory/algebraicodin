# Plot model fit against data

Plots observed data with model predictions from MCMC posterior. Requires
that `save_trajectories = TRUE` was used in
[`fit_mcmc()`](https://catrgory.github.io/algebraicodin/reference/fit_mcmc.md).

## Usage

``` r
plot_fit(
  samples,
  data,
  obs_var,
  state_var = NULL,
  burnin = 0,
  quantiles = c(0.025, 0.25, 0.75, 0.975)
)
```

## Arguments

- samples:

  A monty_samples object with saved trajectories

- data:

  The observed data (data.frame with time and data columns)

- obs_var:

  Name of the observed variable column in data

- state_var:

  Name of the state variable to plot from trajectories

- burnin:

  Number of initial samples to discard

- quantiles:

  Quantile levels for credible intervals (default: 50% and 95%)

## Value

A ggplot object
