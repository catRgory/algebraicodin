# Create a ContinuousMachine

Create a ContinuousMachine

## Usage

``` r
continuous_machine(
  ninputs,
  nstates = NULL,
  noutputs = NULL,
  dynamics = NULL,
  readout = NULL,
  state_names = NULL,
  dynamics_expr = list(),
  readout_expr = list(),
  params = character(0)
)
```

## Arguments

- ninputs:

  Number of input ports

- nstates:

  Number of internal states (inferred from dynamics_expr if omitted)

- noutputs:

  Number of output ports (inferred from readout_expr if omitted)

- dynamics:

  Function (u, x, p, t) -\> du. Auto-built from dynamics_expr if
  omitted.

- readout:

  Function (u, p, t) -\> y. Auto-built from readout_expr if omitted.

- state_names:

  Optional state names

- dynamics_expr:

  Optional named list of quoted expressions per state. Use `input_1`,
  `input_2`, etc. for input port references.

- readout_expr:

  Optional list of quoted expressions, one per output port

- params:

  Character vector of parameter names

## Examples

``` r
# Damped harmonic oscillator with external forcing
osc <- continuous_machine(
  ninputs = 1,
  dynamics_expr = list(
    x = quote(v),
    v = quote(-k * x - c * v + input_1)
  ),
  readout_expr = list(quote(x)),
  params = c("k", "c")
)
osc@ninputs    # 1
#> [1] 1
osc@nstates    # 2
#> [1] 2
osc@noutputs   # 1
#> [1] 1
```
