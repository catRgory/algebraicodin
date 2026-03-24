# Compose multiple Open Petri nets

Convenience function for composing more than two Open Petri nets by
iterative pairwise composition.

## Usage

``` r
compose(...)
```

## Arguments

- ...:

  Open Petri nets

## Value

An Open Petri net

## Examples

``` r
d <- epi_dict()
# Build SEIR from three building blocks
seir <- compose(d$exposure, d$progression, d$recovery)
species_names(apex(seir)) # c("S", "I", "E", "R")
#> [1] "S" "I" "E" "R"
```
