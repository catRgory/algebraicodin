# Extract the closed Petri net from an Open Petri net

Extract the closed Petri net from an Open Petri net

## Usage

``` r
apex(open_pn)
```

## Arguments

- open_pn:

  An Open Petri net

## Value

A Petri net ACSet

## Examples

``` r
pn <- labelled_petri_net(c("S", "I"), inf = "S" %=>% "I")
open_pn <- Open(pn)
closed <- apex(open_pn)
species_names(closed) # c("S", "I")
#> [1] "S" "I"
```
