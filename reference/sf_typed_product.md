# Compute typed product of two typed stock-flow models (stratification)

Creates all valid pairs of (base_obj, strata_obj) where both map to the
same type object, producing a stratified stock-flow model.

## Usage

``` r
sf_typed_product(tsf1, tsf2, independent_params = FALSE)
```

## Arguments

- tsf1:

  A TypedStockFlow (e.g., base disease model)

- tsf2:

  A TypedStockFlow (e.g., strata model like age groups)

- independent_params:

  If TRUE, keep parameters from both models independent (union) instead
  of pairing them via type mapping. Default FALSE.

## Value

A new StockFlowModel with stratified stocks and flows
