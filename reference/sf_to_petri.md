# Convert a stock-flow model to a labelled Petri net

Stocks become species, flows become transitions. Inflows/outflows become
input/output arcs.

## Usage

``` r
sf_to_petri(sf)
```

## Arguments

- sf:

  A StockFlowModel

## Value

A LabelledPetriNet ACSet
