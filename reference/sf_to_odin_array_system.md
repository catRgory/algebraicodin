# Generate an array-based odin2 system from a stratified model

Generate an array-based odin2 system from a stratified model

## Usage

``` r
sf_to_odin_array_system(sf, type = "ode", initial = NULL, compare = NULL)
```

## Arguments

- sf:

  A stratified StockFlowModel (from
  [`sf_stratify()`](https://catrgory.github.io/algebraicodin/reference/sf_stratify.md))

- type:

  One of "ode", "stochastic", "discrete"

- initial:

  Optional named list of initial values

- compare:

  Optional named list of observation specifications

## Value

A function that returns an odin2 generator
