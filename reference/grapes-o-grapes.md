# Compose two Open Petri nets by matching species names

Automatically generates a UWD that connects the two components by
identifying species with the same name, then calls `oapply`.

## Usage

``` r
a %o% b
```

## Arguments

- a, b:

  Open Petri nets, or objects for `base::%o%`

## Value

An Open Petri net (if both are Open), or outer product otherwise

## Details

When applied to non-Open arguments, falls through to `base::%o%` (outer
product).

## Examples

``` r
inf <- exposure_petri("S", "I", "I", "inf")
rec <- spontaneous_petri("I", "R", "rec")
# Compose by matching species "I"
sir <- inf %o% rec
species_names(apex(sir)) # c("S", "I", "R")
#> [1] "S" "I" "R"
```
