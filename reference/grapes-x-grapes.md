# Typed product (stratification) operator

When applied to two TypedPetriNet objects, computes their typed product
(stratification). Falls through to `base::%x%` (Kronecker product) for
other types.

## Usage

``` r
a %x% b
```

## Arguments

- a, b:

  TypedPetriNet objects, or objects for `base::%x%`

## Value

A TypedPetriNet (if both are typed), or Kronecker product otherwise
