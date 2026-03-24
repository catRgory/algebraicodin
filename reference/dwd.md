# Create a directed wiring diagram for Machine composition

Create a directed wiring diagram for Machine composition

## Usage

``` r
dwd(ninputs, noutputs, boxes, wires)
```

## Arguments

- ninputs:

  Number of external input ports

- noutputs:

  Number of external output ports

- boxes:

  Named list where each element is c(ninputs, noutputs)

- wires:

  List of wire specs: c(source_box, source_port, target_box,
  target_port). Use 0 for external input box, -1 for external output
  box.

## Value

A DWD list structure
