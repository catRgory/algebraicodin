# Create a DiscreteResourceSharer

Create a DiscreteResourceSharer

## Usage

``` r
discrete_resource_sharer(
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

## Examples

``` r
# Discrete logistic growth
logistic <- discrete_resource_sharer(
  dynamics_expr = list(
    N = quote(N + r * N * (1 - N / K))
  ),
  params = c("r", "K")
)
logistic@system_type # "discrete"
#> [1] "discrete"
```
