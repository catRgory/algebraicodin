# User-friendly stock-flow stratification

Stratifies a stock-flow model by strata (e.g., age groups, locations).

## Usage

``` r
sf_stratify(
  base,
  strata_names,
  flow_types = NULL,
  cross_strata_flows = list(),
  param_types = NULL,
  mixing = NULL
)
```

## Arguments

- base:

  A StockFlowModel (the base disease model)

- strata_names:

  Character vector of strata names (e.g., c("child", "adult", "senior"))

- flow_types:

  Named character: flow name -\> "disease" or "strata"

- cross_strata_flows:

  Optional list of cross-strata flow specs. Each is a list with `from`,
  `to` (strata names), `rate` (expression).

- param_types:

  Optional named character: param name -\> type param name

- mixing:

  Optional named list specifying contact-matrix mixing for flows. Names
  are base flow names (e.g., "infection"). Values are either:

  - A string: the contact matrix parameter prefix (e.g., "C"). Contact
    parameters will be named `C_<from>_<to>`.

  - A list with elements: `contact` (prefix), `type` ("frequency" or
    "density", default "frequency").

  When mixing is applied, the flow rate uses a force-of-infection that
  sums infectious pressure from all strata weighted by contact rates.

## Value

A stratified StockFlowModel
