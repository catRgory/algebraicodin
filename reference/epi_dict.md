# Pre-built epidemiology dictionary for SIR-family models

Returns a function that maps box names to Open Petri nets.

## Usage

``` r
epi_dict()
```

## Examples

``` r
d <- epi_dict()
names(d) # available building blocks
#> [1] "infection"       "recovery"        "vaccination"     "waning"         
#> [5] "death"           "exposure"        "progression"     "hospitalization"
# Each entry is an Open Petri net
species_names(apex(d$infection)) # c("S", "I")
#> [1] "S" "I"
species_names(apex(d$recovery))  # c("I", "R")
#> [1] "I" "R"
```
