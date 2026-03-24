# Generate odin2 code from a Machine

Converts a (possibly composed) Machine to compilable odin2 code.
Requires symbolic expressions (`dynamics_expr` and optionally
`readout_expr`). External inputs are declared as odin2 parameters.

## Usage

``` r
machine_to_odin(m, initial, params = NULL, input_params = NULL, type = "ode")
```

## Arguments

- m:

  A Machine with dynamics_expr

- initial:

  Named numeric vector of initial conditions

- params:

  Character vector of parameter names. If NULL, uses m@params.

- input_params:

  Named list mapping external input variable names to default values. If
  NULL, external inputs are declared as parameters with no default.

- type:

  One of "ode", "discrete"

## Value

Character string of compilable odin2 code
