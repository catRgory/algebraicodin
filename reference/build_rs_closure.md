# Build a dynamics closure from symbolic expressions

Build a dynamics closure from symbolic expressions

## Usage

``` r
build_rs_closure(state_names, dynamics_expr, param_names)
```

## Arguments

- state_names:

  Character vector of state variable names

- dynamics_expr:

  Named list of quoted expressions

- param_names:

  Character vector of parameter names

## Value

A function(u, p, t) -\> du suitable for ResourceSharer
