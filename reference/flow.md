# Create a flow specification

Create a flow specification

## Usage

``` r
flow(from = NULL, to = NULL, rate)
```

## Arguments

- from:

  Source stock name (NULL or "cloud" for external inflow)

- to:

  Target stock name (NULL or "cloud" for external outflow)

- rate:

  Rate expression (quoted) or function(u, p, t)

## Value

A flow spec list

## Examples

``` r
# Internal flow between stocks
f1 <- flow(from = "S", to = "I", rate = quote(beta * S * I / N))

# External inflow (birth)
f2 <- flow(from = "cloud", to = "S", rate = quote(mu * N))

# External outflow (death)
f3 <- flow(from = "I", to = "cloud", rate = quote(mu * I))
```
