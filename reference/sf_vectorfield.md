# Generate a vectorfield function from a stock-flow model

Generate a vectorfield function from a stock-flow model

## Usage

``` r
sf_vectorfield(sf)
```

## Arguments

- sf:

  A StockFlowModel

## Value

A function with signature (t, state, parms) suitable for deSolve

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
vf <- sf_vectorfield(sir_sf)
state <- c(S = 990, I = 10, R = 0)
parms <- list(beta = 0.3, gamma = 0.1)
dstate <- vf(0, state, parms)
dstate[[1]] # derivatives: dS/dt, dI/dt, dR/dt
#>     S     I     R 
#> -2.97  1.97  1.00 
```
