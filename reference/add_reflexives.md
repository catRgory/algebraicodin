# Add reflexive (self-loop) transitions for stratification

Reflexive transitions are identity transitions that map a species to
itself. They are required so that `typed_product` can pair transitions
from one model with the corresponding species of the other.

## Usage

``` r
add_reflexives(typed_pn, reflexives)
```

## Arguments

- typed_pn:

  A `TypedPetriNet`

- reflexives:

  Named list mapping species names to character vectors of transition
  type names, e.g. `list(S = "strata", I = "strata", R = "strata")`

## Value

A new `TypedPetriNet` with reflexive transitions added

## Examples

``` r
sir <- labelled_petri_net(
  c("S", "I", "R"),
  inf = c("S", "I") %=>% c("I", "I"),
  rec = "I" %=>% "R"
)
onto <- infectious_ontology()
tp <- typed_petri(sir, onto,
  species_types = c(S = "Pop", I = "Pop", R = "Pop"),
  transition_types = c(inf = "infect", rec = "disease")
)
# Add strata-type reflexives for each species
tp_refl <- add_reflexives(tp, list(
  S = "strata", I = "strata", R = "strata"
))
transition_names(tp_refl)
#> [1] "inf"           "rec"           "refl_S_strata" "refl_I_strata"
#> [5] "refl_R_strata"
```
