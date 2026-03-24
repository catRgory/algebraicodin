# Apply contact-matrix mixing to stratified flow expressions

Post-processes a stratified stock-flow model to replace per-stratum
infection expressions with contact-matrix-based force-of-infection.

## Usage

``` r
sf_apply_mixing(sf, base, strata_names, mixing)
```

## Arguments

- sf:

  A stratified StockFlowModel

- base:

  The original (unstratified) base model

- strata_names:

  Character vector of strata names

- mixing:

  Named list: base flow name -\> contact spec

## Value

Modified StockFlowModel with mixing expressions
