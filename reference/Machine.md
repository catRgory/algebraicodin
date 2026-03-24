# A machine: open dynamical system with directed inputs and outputs.

Corresponds to AlgebraicDynamics.jl's Machine. Has explicit input ports
(receiving signals) and output ports (emitting readouts).

## Usage

``` r
Machine(
  system_type = "continuous",
  ninputs = 0L,
  nstates = 0L,
  noutputs = 0L,
  dynamics = function(u, x, p, t) numeric(0),
  readout = function(u, p, t) u,
  state_names = character(0),
  dynamics_expr = list(),
  readout_expr = list(),
  params = character(0)
)
```

## Arguments

- system_type:

  One of "continuous", "discrete", "delay"

- ninputs:

  Number of input ports

- nstates:

  Number of internal state variables

- noutputs:

  Number of output ports

- dynamics:

  Function: `(u, x, p, t) -> du` or `(u, x, h, p, t) -> du` for delay

- readout:

  Function: `(u, p, t) -> y` or `(u, h, p, t) -> y` for delay

- state_names:

  Character vector of state variable names

- dynamics_expr:

  List of symbolic expressions for code generation

- readout_expr:

  List of symbolic readout expressions

- params:

  Character vector of parameter names
