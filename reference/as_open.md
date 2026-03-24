# Convert a catlab StructuredCospan back to an Open Petri net

Extracts the apex and recovers each leg as an integer vector of species
indices from the ACSetTransformation components.

## Usage

``` r
as_open(sc)
```

## Arguments

- sc:

  A catlab::StructuredCospan

## Value

An Open Petri net

## Examples

``` r
sir <- labelled_petri_net(
  c("S", "I", "R"),
  inf = c("S", "I") %=>% c("I", "I"),
  rec = "I" %=>% "R"
)
sc <- open_petri_net(sir, 1:3)
open_sir <- as_open(sc)
species_names(apex(open_sir))
#> [1] "S" "I" "R"
```
