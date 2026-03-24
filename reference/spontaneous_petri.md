# Spontaneous transition: A -\> B

Spontaneous transition: A -\> B

## Usage

``` r
spontaneous_petri(input, output, tname = NULL)
```

## Arguments

- input:

  Name of the input species

- output:

  Name of the output species

- tname:

  Optional transition name

## Value

An Open Petri net

## Examples

``` r
# Recovery: I -> R
rec <- spontaneous_petri("I", "R", "rec")
species_names(apex(rec)) # c("I", "R")
#> [1] "I" "R"

# Waning immunity: R -> S
wan <- spontaneous_petri("R", "S", "wan")
```
