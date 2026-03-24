# Compute the typed product (stratification) of typed Petri nets

This is the core stratification operation. Given two typed Petri nets
sharing the same type system, the typed product creates a new Petri net
whose species are pairs of species with matching types, and whose
transitions are pairs of transitions with matching types.

## Usage

``` r
typed_product(tp1, tp2)
```

## Arguments

- tp1, tp2:

  `TypedPetriNet` objects with the same type system

## Value

A `TypedPetriNet` (the stratified model)

## Examples

``` r
# Stratify SIR by two age groups
onto <- infectious_ontology()

sir <- labelled_petri_net(
  c("S", "I", "R"),
  inf = c("S", "I") %=>% c("I", "I"),
  rec = "I" %=>% "R"
)
tp_sir <- typed_petri(sir, onto,
  species_types = c(S = "Pop", I = "Pop", R = "Pop"),
  transition_types = c(inf = "infect", rec = "disease")
)
tp_sir <- add_reflexives(tp_sir, list(
  S = "strata", I = "strata", R = "strata"
))

ages <- labelled_petri_net(c("young", "old"), aging = "young" %=>% "old")
tp_ages <- typed_petri(ages, onto,
  species_types = c(young = "Pop", old = "Pop"),
  transition_types = c(aging = "strata")
)
tp_ages <- add_reflexives(tp_ages, list(
  young = c("infect", "disease"), old = c("infect", "disease")
))

result <- typed_product(tp_sir, tp_ages)
species_names(result) # e.g. S_young, S_old, I_young, ...
#> [1] "S_young" "S_old"   "I_young" "I_old"   "R_young" "R_old"  
```
