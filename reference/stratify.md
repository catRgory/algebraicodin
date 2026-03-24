# User-friendly stratification wrapper

Computes the typed product of two TypedPetriNets. Equivalent to
`typed_product(base, strata)`.

## Usage

``` r
stratify(base, strata)
```

## Arguments

- base:

  A `TypedPetriNet` (the base model, e.g. SIR)

- strata:

  A `TypedPetriNet` (the stratification model, e.g. age groups)

## Value

A `TypedPetriNet` (the stratified model)

## Examples

``` r
# Same as typed_product(); see that function for a full example
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

result <- stratify(tp_sir, tp_ages)
species_names(result)
#> [1] "S_young" "S_old"   "I_young" "I_old"   "R_young" "R_old"  
```
