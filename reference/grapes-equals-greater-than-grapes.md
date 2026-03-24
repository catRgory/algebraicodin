# Define a Petri net transition: inputs to outputs

Define a Petri net transition: inputs to outputs

## Usage

``` r
inputs %=>% outputs
```

## Arguments

- inputs:

  Character vector of input species names

- outputs:

  Character vector of output species names

## Value

A list with `inputs` and `outputs`

## Examples

``` r
# Infection: S + I -> I + I
inf <- c("S", "I") %=>% c("I", "I")
inf$inputs  # c("S", "I")
#> [1] "S" "I"
inf$outputs # c("I", "I")
#> [1] "I" "I"

# Recovery: I -> R
rec <- "I" %=>% "R"
```
