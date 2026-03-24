# Create a typed stock-flow model

Maps stocks, flows, and parameters of a stock-flow model to a type
template.

## Usage

``` r
sf_typed(
  sf,
  type_sf,
  stock_types = NULL,
  flow_types = NULL,
  param_types = NULL
)
```

## Arguments

- sf:

  A StockFlowModel

- type_sf:

  A StockFlowModel acting as type template

- stock_types:

  Named character: sf stock name -\> type stock name

- flow_types:

  Named character: sf flow name -\> type flow name

- param_types:

  Named character: sf param name -\> type param name

## Value

A TypedStockFlow
