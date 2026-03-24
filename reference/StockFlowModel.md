# A Stock-and-Flow model

Wraps a StockFlow ACSet with additional metadata: flow rate functions
(expressed as R expressions or functions) and operator information.

## Usage

``` r
StockFlowModel(acset = NULL, flow_fns = list(), var_exprs = list())
```

## Arguments

- acset:

  The underlying StockFlow ACSet

- flow_fns:

  Named list of flow rate functions: name -\> function(u, p, t)

- var_exprs:

  Named list of variable expression strings (for display)
