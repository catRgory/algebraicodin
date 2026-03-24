# Standard infectious disease ontology

Returns a type system (LabelledPetriNet) with:

- One species type: `Pop`

- Transition types: `infect` (Pop,Pop -\> Pop,Pop), `disease` (Pop -\>
  Pop), `strata` (Pop -\> Pop)

## Usage

``` r
infectious_ontology()
```

## Value

A LabelledPetriNet ACSet

## Examples

``` r
onto <- infectious_ontology()
species_names(onto)    # "Pop"
#> [1] "Pop"
transition_names(onto) # c("infect", "disease", "strata")
#> [1] "infect"  "disease" "strata" 
```
