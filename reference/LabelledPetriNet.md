# LabelledPetriNet ACSet type

LabelledPetriNet ACSet type

## Usage

``` r
LabelledPetriNet(...)
```

## Arguments

- ...:

  Arguments passed to the ACSet constructor

## Examples

``` r
pn <- LabelledPetriNet()
s1 <- acsets::add_part(pn, "S", sname = "S")
s2 <- acsets::add_part(pn, "S", sname = "I")
t1 <- acsets::add_part(pn, "T", tname = "inf")
acsets::add_part(pn, "I", is = s1, it = t1)
#> [1] 1
acsets::add_part(pn, "O", os = s2, ot = t1)
#> [1] 1
species_names(pn) # c("S", "I")
#> [1] "S" "I"
```
