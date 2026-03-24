# Generate a pure-R mass-action vector field function

Generate a pure-R mass-action vector field function

## Usage

``` r
vectorfield(pn, rate_names = NULL)
```

## Arguments

- pn:

  A Petri net ACSet

- rate_names:

  Optional named vector mapping transition names to rate parameter names

## Value

A function(t, state, parms) suitable for deSolve::ode

## Examples

``` r
sir <- labelled_petri_net(
  c("S", "I", "R"),
  inf = c("S", "I") %=>% c("I", "I"),
  rec = "I" %=>% "R"
)
vf <- vectorfield(sir)
state <- c(S = 990, I = 10, R = 0)
parms <- list(inf = 0.001, rec = 0.1)
# Evaluate derivatives at t = 0
dstate <- vf(0, state, parms)
dstate[[1]] # named numeric vector of derivatives
#>    S    I    R 
#> -9.9  8.9  1.0 
```
