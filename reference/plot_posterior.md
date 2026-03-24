# Plot posterior density estimates

Plot posterior density estimates

## Usage

``` r
plot_posterior(samples, pars = NULL, burnin = 0, true_values = NULL)
```

## Arguments

- samples:

  A monty_samples object from
  [`fit_mcmc()`](https://catrgory.github.io/algebraicodin/reference/fit_mcmc.md)

- pars:

  Optional character vector of parameter names to plot

- burnin:

  Number of initial samples to discard

- true_values:

  Optional named list of true parameter values (shown as vertical lines)

## Value

A ggplot object
