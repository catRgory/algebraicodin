# ReactionNet ACSet type (labelled with rates and concentrations)

ReactionNet ACSet type (labelled with rates and concentrations)

## Usage

``` r
ReactionNet(...)
```

## Arguments

- ...:

  Arguments passed to the ACSet constructor

## Examples

``` r
rn <- ReactionNet()
s1 <- acsets::add_part(rn, "S", sname = "S", concentration = 990)
s2 <- acsets::add_part(rn, "S", sname = "I", concentration = 10)
t1 <- acsets::add_part(rn, "T", tname = "inf", rate = 0.001)
acsets::add_part(rn, "I", is = s1, it = t1)
#> [1] 1
acsets::add_part(rn, "O", os = s2, ot = t1)
#> [1] 1
```
