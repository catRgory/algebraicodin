# Build a dynamics closure from symbolic expressions for a Machine

Build a dynamics closure from symbolic expressions for a Machine

## Usage

``` r
build_machine_closure(state_names, ninputs, dynamics_expr)
```

## Arguments

- state_names:

  Character vector of state variable names

- ninputs:

  Number of input ports

- dynamics_expr:

  Named list of quoted expressions

## Value

A function(u, x, p, t) -\> du
