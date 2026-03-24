# Generate array-based odin2 code from a stratified stock-flow model

Produces compact odin2 code using arrays indexed by strata, instead of
separate scalar variables per stratum. This scales much better for
models with many strata (\>5). The model must have been created with
[`sf_stratify()`](https://catrgory.github.io/algebraicodin/reference/sf_stratify.md),
which attaches the strata metadata needed.

## Usage

``` r
sf_to_odin_array(sf, type = "ode", initial = NULL, compare = NULL)
```

## Arguments

- sf:

  A stratified StockFlowModel (from
  [`sf_stratify()`](https://catrgory.github.io/algebraicodin/reference/sf_stratify.md))

- type:

  One of "ode" (default), "stochastic", "discrete"

- initial:

  Optional named list of initial values per base stock (vectors of
  length n_strata)

- compare:

  Optional named list of observation specifications

## Value

Character string of odin2 code using arrays
