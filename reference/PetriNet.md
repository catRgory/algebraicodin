# PetriNet ACSet type (unlabelled)

PetriNet ACSet type (unlabelled)

## Usage

``` r
PetriNet(...)
```

## Arguments

- ...:

  Arguments passed to the ACSet constructor

## Examples

``` r
pn <- PetriNet()
s1 <- acsets::add_part(pn, "S")
t1 <- acsets::add_part(pn, "T")
acsets::add_part(pn, "I", is = s1, it = t1)
#> [1] 1
acsets::add_part(pn, "O", os = s1, ot = t1)
#> [1] 1
acsets::nparts(pn, "S") # 1
#> [1] 1
```
