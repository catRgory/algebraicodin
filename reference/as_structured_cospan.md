# Convert an Open Petri net to a catlab StructuredCospan

Each integer-vector leg becomes an ACSetTransformation from a discrete
Petri net foot (with parts only in "S") to the apex.

## Usage

``` r
as_structured_cospan(open_pn)
```

## Arguments

- open_pn:

  An Open Petri net

## Value

A catlab::StructuredCospan with interface_ob = "S"

## Examples

``` r
sir <- labelled_petri_net(
  c("S", "I", "R"),
  inf = c("S", "I") %=>% c("I", "I"),
  rec = "I" %=>% "R"
)
open_sir <- Open(sir)
sc <- as_structured_cospan(open_sir)
# Round-trip back to Open
open_again <- as_open(sc)
species_names(apex(open_again))
#> [1] "S" "I" "R"
```
