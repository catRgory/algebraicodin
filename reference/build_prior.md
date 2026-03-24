# Build a prior model from specifications

Constructs a monty_model for use as a prior in
[`fit_mcmc()`](https://catrgory.github.io/algebraicodin/reference/fit_mcmc.md).
Each argument is a named prior specification created by
[`prior_exp()`](https://catrgory.github.io/algebraicodin/reference/prior_exp.md),
[`prior_uniform()`](https://catrgory.github.io/algebraicodin/reference/prior_uniform.md),
or
[`prior_normal()`](https://catrgory.github.io/algebraicodin/reference/prior_normal.md).

## Usage

``` r
build_prior(...)
```

## Arguments

- ...:

  Named prior specifications

## Value

A monty_model object

## Examples

``` r
if (FALSE) { # \dontrun{
p <- build_prior(
  beta = prior_exp(mean = 0.3),
  gamma = prior_exp(mean = 0.1)
)
} # }
```
