# Compose Machines via a directed wiring diagram

Each box of the directed wiring diagram is filled by a Machine. Outputs
of upstream machines are wired to inputs of downstream ones.

## Usage

``` r
oapply_machine(d, machines)
```

## Arguments

- d:

  A directed wiring diagram (from
  [`dwd()`](https://catrgory.github.io/algebraicodin/reference/dwd.md))

- machines:

  A list of Machine objects, one per box

## Value

A composite Machine
