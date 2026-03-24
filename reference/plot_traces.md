# Plot MCMC trace plots

Plot MCMC trace plots

## Usage

``` r
plot_traces(samples, pars = NULL, burnin = 0)
```

## Arguments

- samples:

  A monty_samples object from
  [`fit_mcmc()`](https://catrgory.github.io/algebraicodin/reference/fit_mcmc.md)

- pars:

  Optional character vector of parameter names to plot

- burnin:

  Number of initial samples to shade as burn-in

## Value

A ggplot object
