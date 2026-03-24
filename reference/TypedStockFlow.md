# Typed stock-flow model for stratification

Typed stock-flow model for stratification

## Usage

``` r
TypedStockFlow(
  sf,
  type_sf,
  stock_map,
  flow_map,
  param_map,
  refl_flows = character(0)
)
```

## Arguments

- sf:

  A StockFlowModel

- type_sf:

  A StockFlowModel acting as the type template

- stock_map:

  Named character: base stock name -\> type stock name

- flow_map:

  Named character: base flow name -\> type flow name

- param_map:

  Named character: base param name -\> type param name

- refl_flows:

  Character vector of reflexive flow names

## Value

A TypedStockFlow object
