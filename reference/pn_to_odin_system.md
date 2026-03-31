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

  Currently only `"mass_action"` is implemented. `"custom"` is reserved
  and will error if requested.

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
# \donttest{
sir <- labelled_petri_net(
  c("S", "I", "R"),
  inf = c("S", "I") %=>% c("I", "I"),
  rec = "I" %=>% "R"
)
gen <- pn_to_odin_system(sir, type = "ode")
sys <- gen() # creates an odin2 system (requires odin2)
#> ✔ Wrote 'DESCRIPTION'
#> ✔ Wrote 'NAMESPACE'
#> ✔ Wrote 'R/dust.R'
#> ✔ Wrote 'src/dust.cpp'
#> ✔ Wrote 'src/Makevars'
#> ℹ 13 functions decorated with [[cpp11::register]]
#> ✔ generated file cpp11.R
#> ✔ generated file cpp11.cpp
#> ℹ Re-compiling odin.system840d9d16
#> ── R CMD INSTALL ───────────────────────────────────────────────────────────────
#> * installing *source* package ‘odin.system840d9d16’ ...
#> ** this is package ‘odin.system840d9d16’ version ‘0.0.1’
#> ** using staged installation
#> ** libs
#> using C++ compiler: ‘g++ (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0’
#> g++ -std=gnu++17 -I"/opt/R/4.5.3/lib/R/include" -DNDEBUG  -I'/home/runner/work/_temp/Library/cpp11/include' -I'/home/runner/work/_temp/Library/dust2/include' -I'/home/runner/work/_temp/Library/monty/include' -I/usr/local/include   -DHAVE_INLINE -fopenmp  -fpic  -g -O2  -Wall -pedantic -fdiagnostics-color=always  -c cpp11.cpp -o cpp11.o
#> g++ -std=gnu++17 -I"/opt/R/4.5.3/lib/R/include" -DNDEBUG  -I'/home/runner/work/_temp/Library/cpp11/include' -I'/home/runner/work/_temp/Library/dust2/include' -I'/home/runner/work/_temp/Library/monty/include' -I/usr/local/include   -DHAVE_INLINE -fopenmp  -fpic  -g -O2  -Wall -pedantic -fdiagnostics-color=always  -c dust.cpp -o dust.o
#> g++ -std=gnu++17 -shared -L/opt/R/4.5.3/lib/R/lib -L/usr/local/lib -o odin.system840d9d16.so cpp11.o dust.o -fopenmp -L/opt/R/4.5.3/lib/R/lib -lR
#> installing to /tmp/RtmpYCTZRF/devtools_install_1ce444e77696/00LOCK-dust_1ce464cfb89b/00new/odin.system840d9d16/libs
#> ** checking absolute paths in shared objects and dynamic libraries
#> * DONE (odin.system840d9d16)
#> ℹ Loading odin.system840d9d16
# }
```
