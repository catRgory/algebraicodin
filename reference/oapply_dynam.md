# Compose ResourceSharers via an undirected wiring diagram

Implements the operad algebra for undirected composition of dynamical
systems. Each box of the UWD is filled by a resource sharer from
`sharers`. States at the same junction are identified and derivatives
(or updates) are summed.

## Usage

``` r
oapply_dynam(d, sharers)
```

## Arguments

- d:

  A UWD ACSet (from catlab)

- sharers:

  A list of ResourceSharer objects, one per box

## Value

A composite ResourceSharer
