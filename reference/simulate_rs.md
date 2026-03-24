# Simulate a continuous ResourceSharer using deSolve

Simulate a continuous ResourceSharer using deSolve

## Usage

``` r
simulate_rs(rs, initial, times, params = list())
```

## Arguments

- rs:

  A continuous ResourceSharer

- initial:

  Named numeric vector of initial conditions

- times:

  Numeric vector of output times

- params:

  Parameter list (passed to dynamics as p)

## Value

A data.frame with time and state columns
