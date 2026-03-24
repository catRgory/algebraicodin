# Compose two Open Petri nets by matching species names

Creates a UWD where species with the same name are identified (shared
junctions), then composes via `oapply`.

## Usage

``` r
compose_open(a, b)
```

## Arguments

- a, b:

  Open Petri nets

## Value

An Open Petri net

## Examples

``` r
inf <- exposure_petri("S", "I", "I", "inf")
rec <- spontaneous_petri("I", "R", "rec")
sir <- compose_open(inf, rec)
species_names(apex(sir)) # c("S", "I", "R")
#> [1] "S" "I" "R"
```
