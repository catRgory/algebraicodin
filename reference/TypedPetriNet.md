# Typed Petri net: a Petri net with a morphism to a type system.

Typed Petri net: a Petri net with a morphism to a type system.

## Usage

``` r
TypedPetriNet(
  pn,
  type_system,
  species_map,
  transition_map,
  input_map,
  output_map,
  refl_species = character(0)
)
```

## Arguments

- pn:

  A LabelledPetriNet ACSet

- type_system:

  A LabelledPetriNet ACSet (the ontology)

- species_map:

  Integer vector mapping PN species -\> type system species

- transition_map:

  Integer vector mapping PN transitions -\> type system transitions

- input_map:

  Integer vector mapping PN input arcs -\> type system input arcs

- output_map:

  Integer vector mapping PN output arcs -\> type system output arcs

- refl_species:

  Character vector: for each transition, the species name it loops on
  (non-empty only for reflexive transitions added by `add_reflexives`).
