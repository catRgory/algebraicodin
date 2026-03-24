# Open Petri net: a Petri net with exposed "legs" (species accessible from outside for composition).

Open Petri net: a Petri net with exposed "legs" (species accessible from
outside for composition).

## Usage

``` r
Open(pn, legs = NULL)
```

## Arguments

- pn:

  A Petri net ACSet

- legs:

  List of integer vectors mapping outer junctions to species IDs

## Examples

``` r
sir <- labelled_petri_net(
  c("S", "I", "R"),
  inf = c("S", "I") %=>% c("I", "I"),
  rec = "I" %=>% "R"
)
# Expose all species as a single leg
open_sir <- Open(sir)
species_names(apex(open_sir)) # c("S", "I", "R")
#> [1] "S" "I" "R"

# Expose specific species
open_si <- Open(sir, legs = list(c(1L, 2L)))
```
