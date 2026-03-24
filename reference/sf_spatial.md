# Create a spatial (multi-patch) model via stratification

Convenience wrapper around
[`sf_stratify`](https://catrgory.github.io/algebraicodin/reference/sf_stratify.md)
for spatial models. Each patch runs independent disease dynamics, with
optional migration between patches and optional cross-patch transmission
via a contact matrix.

## Usage

``` r
sf_spatial(
  base,
  patch_names,
  flow_types = NULL,
  migration = list(),
  contact = NULL
)
```

## Arguments

- base:

  A StockFlowModel (the within-patch disease model)

- patch_names:

  Character vector of patch names

- flow_types:

  Named character: base flow name -\> type (default all "disease")

- migration:

  Named list of migration specs. Names are migration flow names. Each
  entry is a list with:

  - `from`, `to`: patch names

  - `rate`: quoted R expression (use the source patch name as the stock
    variable, e.g., `quote(m * patch1)`)

  - `params`: character vector of parameter names used in rate

- contact:

  Optional string: contact matrix parameter prefix. When provided,
  transmission between patches uses a contact matrix
  (frequency-dependent). Contact parameters are named
  `<prefix>_<to>_<from>`.

## Value

A stratified StockFlowModel with spatial structure
