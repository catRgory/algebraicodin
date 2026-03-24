# Create an open stock-flow model

Wraps a StockFlowModel with interface specification for composition. The
feet are StockFlow0 ACSets specifying which stocks and sum variables are
exposed.

## Usage

``` r
open_stock_flow(sf, ...)
```

## Arguments

- sf:

  A StockFlowModel

- ...:

  Character vectors of exposed stock names (one per foot/leg)

## Value

An OpenStockFlow object

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
open_sf <- open_stock_flow(sir_sf, c("S", "I", "R"))
```
