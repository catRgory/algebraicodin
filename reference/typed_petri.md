# Create a typed Petri net from type assignments

Create a typed Petri net from type assignments

## Usage

``` r
typed_petri(pn, type_system, species_types = NULL, transition_types = NULL)
```

## Arguments

- pn:

  A LabelledPetriNet

- type_system:

  A LabelledPetriNet (the ontology)

- species_types:

  Named character vector mapping species names to type system species
  names, e.g. `c(S = "Pop", I = "Pop")`

- transition_types:

  Named character vector mapping transition names to type system
  transition names, e.g. `c(inf = "infect", rec = "disease")`

## Value

A `TypedPetriNet`

## Examples

``` r
# Type an SIR model over the infectious disease ontology
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
species_names(tp) # c("S", "I", "R")
#> [1] "S" "I" "R"
```
