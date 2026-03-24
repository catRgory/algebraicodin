# Create a DelayMachine

Create a DelayMachine

## Usage

``` r
delay_machine(
  ninputs,
  nstates,
  noutputs,
  dynamics,
  readout,
  state_names = paste0("x", seq_len(nstates)),
  dynamics_expr = list(),
  readout_expr = list(),
  params = character(0)
)
```

## Arguments

- ninputs:

  Number of input ports

- nstates:

  Number of internal states

- noutputs:

  Number of output ports

- dynamics:

  Function (u, x, h, p, t) -\> du

- readout:

  Function (u, h, p, t) -\> y

- state_names:

  Optional state names

- dynamics_expr:

  Optional named list of expressions (for code gen)

- readout_expr:

  Optional list of readout expressions

- params:

  Parameter names
