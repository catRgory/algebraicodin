# Convert a Petri net to a DelayResourceSharer (DDE)

Creates a DDE system where specified transitions use delayed state
values.

## Usage

``` r
petri_to_delay(pn, delays)
```

## Arguments

- pn:

  A Petri net ACSet

- delays:

  Named list mapping transition names to delay specs. Each spec is a
  list with `tau` (delay time) and optionally `species` (which input
  species are delayed; default: all inputs).
