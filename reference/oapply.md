# Compose open Petri nets over a UWD

The UWD specifies how boxes (= open PNs) connect via junctions. Each
box's ports connect to junctions; shared junctions mean species are
identified (glued together).

## Usage

``` r
oapply(w, components)
```

## Arguments

- w:

  A UWD ACSet

- components:

  List of Open Petri nets, one per box

## Value

An Open Petri net (the composed system)

## Examples

``` r
# Build SIR by composing infection and recovery
infection <- exposure_petri("S", "I", "I", "inf")
recovery <- spontaneous_petri("I", "R", "rec")

# Create a UWD with shared junction "I"
w <- catlab::uwd(
  outer = c("S", "I", "R"),
  infection = c("S", "I"),
  recovery = c("I", "R")
)
sir <- oapply(w, list(infection, recovery))
species_names(apex(sir)) # c("S", "I", "R")
#> [1] "S" "I" "R"
```
