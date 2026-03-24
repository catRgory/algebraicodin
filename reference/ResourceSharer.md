# A resource sharer: open dynamical system with symmetric ports.

Corresponds to AlgebraicDynamics.jl's ResourceSharer. The `portmap`
specifies which internal state variables are exposed at each port.

## Usage

``` r
ResourceSharer(
  system_type = "continuous",
  nstates = 0L,
  dynamics = function(u, p, t) numeric(0),
  portmap = seq_len(nstates),
  state_names = character(0),
  dynamics_expr = list(),
  params = character(0)
)
```

## Arguments

- system_type:

  One of "continuous", "discrete", "delay"

- nstates:

  Number of internal state variables

- dynamics:

  Function implementing the system evolution

- portmap:

  Integer vector mapping ports to state indices

- state_names:

  Character vector of state variable names

- dynamics_expr:

  List of symbolic expressions for code generation

- params:

  Character vector of parameter names
