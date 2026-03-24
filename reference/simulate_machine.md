# Simulate a Machine with external inputs

Simulate a Machine with external inputs

## Usage

``` r
simulate_machine(
  m,
  initial,
  times,
  input_fn = function(t) numeric(m@ninputs),
  params = list()
)
```

## Arguments

- m:

  A continuous Machine

- initial:

  Named numeric vector of initial state

- times:

  Numeric vector of output times

- input_fn:

  Function (t) -\> numeric vector of inputs

- params:

  Parameter list

## Value

A data.frame
