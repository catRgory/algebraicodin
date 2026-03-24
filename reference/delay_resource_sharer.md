# Create a DelayResourceSharer from a DDE dynamics function

Create a DelayResourceSharer from a DDE dynamics function

## Usage

``` r
delay_resource_sharer(
  nstates,
  dynamics,
  portmap = seq_len(nstates),
  state_names = paste0("x", seq_len(nstates)),
  dynamics_expr = list(),
  params = character(0)
)
```

## Arguments

- nstates:

  Number of state variables

- dynamics:

  Function (u, h, p, t) -\> du, where h is a history function

- portmap:

  Integer vector mapping ports to states

- state_names:

  Optional state names

- dynamics_expr:

  Optional named list of expressions (for code gen)

- params:

  Parameter names

## Value

A ResourceSharer with system_type = "delay"
