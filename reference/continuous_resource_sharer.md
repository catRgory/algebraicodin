# Create a ContinuousResourceSharer

Create a ContinuousResourceSharer

## Usage

``` r
continuous_resource_sharer(
  nstates = NULL,
  dynamics = NULL,
  portmap = NULL,
  state_names = NULL,
  dynamics_expr = list(),
  params = character(0)
)
```

## Arguments

- nstates:

  Number of state variables (inferred from dynamics_expr if omitted)

- dynamics:

  Function (u, p, t) -\> du. Auto-built from dynamics_expr if omitted.

- portmap:

  Integer vector mapping ports to states (default: identity)

- state_names:

  Optional character vector of state names

- dynamics_expr:

  Optional named list of quoted R expressions, one per state. Enables
  odin2 code generation via
  [`rs_to_odin()`](https://catrgory.github.io/algebraicodin/reference/rs_to_odin.md).

- params:

  Character vector of parameter names used in dynamics_expr

## Value

A ResourceSharer with system_type = "continuous"

## Examples

``` r
# Lotka-Volterra predator-prey via symbolic expressions
lv <- continuous_resource_sharer(
  dynamics_expr = list(
    prey = quote(alpha * prey - beta * prey * predator),
    predator = quote(delta * prey * predator - gamma * predator)
  ),
  params = c("alpha", "beta", "delta", "gamma")
)
lv@system_type  # "continuous"
#> [1] "continuous"
lv@state_names  # c("prey", "predator")
#> [1] "prey"     "predator"
nports(lv)      # 2
#> [1] 2
```
