# Open a Petri net at specified species, returning a StructuredCospan

Convenience wrapper around catlab::open_acset with interface_ob = "S".

## Usage

``` r
open_petri_net(pn, ...)
```

## Arguments

- pn:

  A Petri net ACSet

- ...:

  Integer vectors, each specifying a leg (species indices)

## Value

A catlab::StructuredCospan with interface_ob = "S"

## Examples

``` r
sir <- labelled_petri_net(
  c("S", "I", "R"),
  inf = c("S", "I") %=>% c("I", "I"),
  rec = "I" %=>% "R"
)
sc <- open_petri_net(sir, c(1L, 2L), c(2L, 3L))
```
