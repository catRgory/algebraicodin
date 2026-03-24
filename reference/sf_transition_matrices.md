# Compute inflow/outflow transition matrices for a stock-flow model

Compute inflow/outflow transition matrices for a stock-flow model

## Usage

``` r
sf_transition_matrices(sf)
```

## Arguments

- sf:

  A StockFlowModel

## Value

A list with `inflow` and `outflow` matrices (rows = flows, cols =
stocks)

## Examples

``` r
sir_sf <- stock_and_flow(
  stocks = c("S", "I", "R"),
  flows = list(
    inf = flow(from = "S", to = "I", rate = quote(beta * S * I / N)),
    rec = flow(from = "I", to = "R", rate = quote(gamma * I))
  ),
  params = c("beta", "gamma"),
  sums = list(N = c("S", "I", "R"))
)
tm <- sf_transition_matrices(sir_sf)
tm$inflow  # flows x stocks inflow matrix
#>     S I R
#> inf 0 1 0
#> rec 0 0 1
tm$outflow # flows x stocks outflow matrix
#>     S I R
#> inf 1 0 0
#> rec 0 1 0
```
