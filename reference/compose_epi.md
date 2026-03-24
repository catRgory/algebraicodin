# Compose an epidemiological model from a UWD and dictionary

Compose an epidemiological model from a UWD and dictionary

## Usage

``` r
compose_epi(w, dict, box_names)
```

## Arguments

- w:

  A UWD (wiring diagram)

- dict:

  A named list of Open Petri nets, keyed by box names

- box_names:

  Character vector mapping box indices to dict keys

## Value

An Open Petri net

## Examples

``` r
# Build SIR from epi building blocks
d <- epi_dict()
w <- catlab::uwd(
  outer = c("S", "I", "R"),
  infection = c("S", "I"),
  recovery = c("I", "R")
)
sir <- compose_epi(w, d, c("infection", "recovery"))
species_names(apex(sir)) # c("S", "I", "R")
#> [1] "S" "I" "R"
```
