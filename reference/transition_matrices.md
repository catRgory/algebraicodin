# Extract input/output stoichiometry matrices

Extract input/output stoichiometry matrices

## Usage

``` r
transition_matrices(pn)
```

## Arguments

- pn:

  A Petri net ACSet

## Value

A list with `input` and `output` matrices (species × transitions)

## Examples

``` r
sir <- labelled_petri_net(
  c("S", "I", "R"),
  inf = c("S", "I") %=>% c("I", "I"),
  rec = "I" %=>% "R"
)
tm <- transition_matrices(sir)
tm$input         # species x transitions input matrix
#>   inf rec
#> S   1   0
#> I   1   1
#> R   0   0
tm$output        # species x transitions output matrix
#>   inf rec
#> S   0   0
#> I   2   0
#> R   0   1
tm$stoichiometry # net change: output - input
#>   inf rec
#> S  -1   0
#> I   1  -1
#> R   0   1
```
