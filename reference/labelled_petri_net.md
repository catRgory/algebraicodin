# Build a LabelledPetriNet from species names and transition specs

Build a LabelledPetriNet from species names and transition specs

## Usage

``` r
labelled_petri_net(species, ...)
```

## Arguments

- species:

  Character vector of species names

- ...:

  Named transition specs: `name = inputs %=>% outputs` where
  inputs/outputs are character vectors of species names

## Value

A LabelledPetriNet ACSet

## Examples

``` r
# SIR model as a Petri net
sir <- labelled_petri_net(
  c("S", "I", "R"),
  inf = c("S", "I") %=>% c("I", "I"),
  rec = "I" %=>% "R"
)
species_names(sir)    # c("S", "I", "R")
#> [1] "S" "I" "R"
transition_names(sir) # c("inf", "rec")
#> [1] "inf" "rec"

# SIS model (no recovery compartment needed)
sis <- labelled_petri_net(
  c("S", "I"),
  inf = c("S", "I") %=>% c("I", "I"),
  rec = "I" %=>% "S"
)
```
