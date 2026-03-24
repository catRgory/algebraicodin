# Generate odin2 code from a ResourceSharer

Converts a (possibly composed) ResourceSharer to compilable odin2 code.
Requires that the ResourceSharer was created with `dynamics_expr`
(symbolic expressions). If only a closure is available, code generation
is not possible and an error is thrown.

## Usage

``` r
rs_to_odin(rs, initial, params = NULL, type = "ode")
```

## Arguments

- rs:

  A ResourceSharer with dynamics_expr

- initial:

  Named numeric vector of initial conditions

- params:

  Character vector of parameter names. If NULL, uses rs@params.

- type:

  One of "ode" (continuous ODE), "discrete" (Euler step), "stochastic"
  (Euler-Maruyama with Binomial draws)

## Value

Character string of compilable odin2 code
