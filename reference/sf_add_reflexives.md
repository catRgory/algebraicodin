# Add reflexive (identity) flows for stock-flow stratification

Adds self-loop flows for each stock, typed to a given flow type. Needed
so the typed product can pair real flows from one model with identity
flows on corresponding stocks of the other.

## Usage

``` r
sf_add_reflexives(tsf, reflexives)
```

## Arguments

- tsf:

  A TypedStockFlow

- reflexives:

  Named list: stock_name -\> character vector of flow type names

## Value

A new TypedStockFlow with reflexive flows added
