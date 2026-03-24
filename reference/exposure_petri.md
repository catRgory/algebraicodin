# Exposure/contact transition: A + B -\> C + B (B catalyzes A -\> C)

Exposure/contact transition: A + B -\> C + B (B catalyzes A -\> C)

## Usage

``` r
exposure_petri(susceptible, infectious, output, tname = NULL)
```

## Arguments

- susceptible:

  Name of the susceptible species

- infectious:

  Name of the infectious species

- output:

  Name of the output species

- tname:

  Optional transition name

## Value

An Open Petri net

## Examples

``` r
# Standard infection: S + I -> I + I
inf <- exposure_petri("S", "I", "I", "inf")
species_names(apex(inf)) # c("S", "I")
#> [1] "S" "I"

# Exposure (SEIR): S + I -> E + I
exp_petri <- exposure_petri("S", "I", "E", "exp")
species_names(apex(exp_petri)) # c("S", "I", "E")
#> [1] "S" "I" "E"
```
