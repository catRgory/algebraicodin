# Generate odin2 code and return it as a callable odin_system generator

If odin2 is available, this returns a function that when called creates
an odin system. Otherwise returns the code string.

## Usage

``` r
pn_to_odin_system(
  pn,
  type = c("ode", "dde", "stochastic", "discrete"),
  rate_style = "mass_action",
  delays = NULL,
  compare = NULL
)

to_odin_system(
  pn,
  type = c("ode", "dde", "stochastic", "discrete"),
  rate_style = "mass_action",
  delays = NULL,
  compare = NULL
)
```

## Arguments

- pn:

  A Petri net ACSet (LabelledPetriNet or ReactionNet)

- type:

  Model type:

  - `"ode"`: continuous ODE with
    [`deriv()`](https://rdrr.io/r/stats/deriv.html)

  - `"dde"`: delay differential equation with
    [`deriv()`](https://rdrr.io/r/stats/deriv.html) + `delay()`

  - `"stochastic"`: discrete-time stochastic with
    [`update()`](https://rdrr.io/r/stats/update.html) + `Binomial()`

  - `"discrete"`: discrete-time deterministic Euler with
    [`update()`](https://rdrr.io/r/stats/update.html)

- rate_style:

  "mass_action" (default) or "custom"

- delays:

  Named list for DDE: maps transition names to delay specs. Each element
  is a list with `tau` (numeric delay time) and optionally `species`
  (character vector of which input species are delayed).

- compare:

  Optional named list mapping observed variables to distribution
  families (e.g., `list(cases = "Poisson")`)

## Value

A function that creates an odin system, or a code string

## Examples

``` r
if (FALSE) { # \dontrun{
sir <- labelled_petri_net(
  c("S", "I", "R"),
  inf = c("S", "I") %=>% c("I", "I"),
  rec = "I" %=>% "R"
)
gen <- pn_to_odin_system(sir, type = "ode")
sys <- gen() # creates an odin2 system (requires odin2)
} # }
```
