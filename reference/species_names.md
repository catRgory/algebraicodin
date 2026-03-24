# Get species names from a Petri net

Get species names from a Petri net

## Usage

``` r
species_names(pn)
```

## Arguments

- pn:

  A Petri net ACSet (or TypedPetriNet)

## Value

Character vector of species names

## Examples

``` r
sir <- labelled_petri_net(
  c("S", "I", "R"),
  inf = c("S", "I") %=>% c("I", "I"),
  rec = "I" %=>% "R"
)
species_names(sir)    # c("S", "I", "R")
#> [1] "S" "I" "R"
transition_names(sir) # c("inf", "rec")
#> [1] "inf" "rec"
```
