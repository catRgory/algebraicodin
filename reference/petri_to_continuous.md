# Convert a Petri net to a ContinuousResourceSharer (ODE)

Convert a Petri net to a ContinuousResourceSharer (ODE)

## Usage

``` r
petri_to_continuous(pn)
```

## Arguments

- pn:

  A Petri net ACSet

## Value

A ResourceSharer with system_type "continuous"

## Examples

``` r
sir <- labelled_petri_net(
  c("S", "I", "R"),
  inf = c("S", "I") %=>% c("I", "I"),
  rec = "I" %=>% "R"
)
rs <- petri_to_continuous(sir)
rs@system_type   # "continuous"
#> [1] "continuous"
rs@state_names   # c("S", "I", "R")
#> [1] "S" "I" "R"
```
