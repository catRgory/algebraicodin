# Build a stock-and-flow model using a declarative interface

Build a stock-and-flow model using a declarative interface

## Usage

``` r
stock_and_flow(stocks, flows, params = character(0), sums = list())
```

## Arguments

- stocks:

  Character vector of stock names

- flows:

  Named list of flow specifications. Each flow is a list with `from`
  (source stock or NULL for external), `to` (target stock or NULL for
  external), and `rate` (an R expression or function).

- params:

  Character vector of parameter names

- sums:

  Named list of sum variable definitions: `list(N = c("S", "I", "R"))`
  means N = S + I + R

## Value

A StockFlowModel

## Examples

``` r
# SIR stock-flow model
sir_sf <- stock_and_flow(
  stocks = c("S", "I", "R"),
  flows = list(
    inf = flow(from = "S", to = "I", rate = quote(beta * S * I / N)),
    rec = flow(from = "I", to = "R", rate = quote(gamma * I))
  ),
  params = c("beta", "gamma"),
  sums = list(N = c("S", "I", "R"))
)
sf_snames(sir_sf) # c("S", "I", "R")
#> [1] "S" "I" "R"
sf_fnames(sir_sf) # c("inf", "rec")
#> [1] "inf" "rec"
```
