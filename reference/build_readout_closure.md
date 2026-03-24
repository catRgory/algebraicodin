# Build a readout closure from symbolic expressions for a Machine

Build a readout closure from symbolic expressions for a Machine

## Usage

``` r
build_readout_closure(state_names, readout_expr)
```

## Arguments

- state_names:

  Character vector of state variable names

- readout_expr:

  List of quoted expressions (one per output port)

## Value

A function(u, p, t) -\> y
