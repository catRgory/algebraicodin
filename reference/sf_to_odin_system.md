# Generate an odin2 system from a stock-flow model

Returns a function that, when called, creates an odin2 system generator.

## Usage

``` r
sf_to_odin_system(sf, type = "ode", initial = NULL, compare = NULL)
```

## Arguments

- sf:

  A StockFlowModel

- type:

  One of "ode", "stochastic", "discrete"

- initial:

  Optional named list of initial values

- compare:

  Optional named list of observation specifications

## Value

A function that returns an odin2 generator
