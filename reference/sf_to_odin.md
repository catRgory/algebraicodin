# Generate odin2 code from a stock-flow model

Generate odin2 code from a stock-flow model

## Usage

``` r
sf_to_odin(sf, type = "ode", initial = NULL, compare = NULL)
```

## Arguments

- sf:

  A StockFlowModel

- type:

  One of "ode", "stochastic", "discrete"

- initial:

  Named numeric vector of initial conditions (optional)

- compare:

  Optional named list of observation specifications (see
  [`observe()`](https://catrgory.github.io/algebraicodin/reference/observe.md))

## Value

Character string of odin2 code

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
code <- sf_to_odin(sir_sf, type = "ode")
cat(code)
#> beta <- parameter()
#> gamma <- parameter()
#> 
#> S0 <- parameter(0)
#> initial(S) <- S0
#> I0 <- parameter(0)
#> initial(I) <- I0
#> R0 <- parameter(0)
#> initial(R) <- R0
#> 
#> N <- S + I + R
#> 
#> v_inf <- beta * S * I/N
#> v_rec <- gamma * I
#> 
#> deriv(S) <- - v_inf
#> deriv(I) <- v_inf - v_rec
#> deriv(R) <- v_rec
```
