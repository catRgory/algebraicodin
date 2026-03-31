# algebraicodin

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/catRgory/algebraicodin/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/catRgory/algebraicodin/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/catRgory/algebraicodin/graph/badge.svg)](https://codecov.io/gh/catRgory/algebraicodin)
<!-- badges: end -->

**Compositional epidemic modeling using category theory.**

Build Petri nets and stock-flow diagrams, compose them via open systems and
undirected wiring diagrams (UWDs), stratify models by age, geography, or risk
group, and generate [odin2](https://mrc-ide.github.io/odin2)/[dust2](https://mrc-ide.github.io/dust2)
code for high-performance simulation. An R port of
[AlgebraicPetri.jl](https://github.com/AlgebraicJulia/AlgebraicPetri.jl) and
[StockFlow.jl](https://github.com/AlgebraicJulia/StockFlow.jl) from the
[AlgebraicJulia](https://www.algebraicjulia.org/) ecosystem.

## Installation

Install from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("catRgory/algebraicodin")
```

algebraicodin depends on the **catRgory** companion packages
[acsets](https://github.com/catRgory/acsets) and
[catlab](https://github.com/catRgory/catlab), plus
[odin2](https://mrc-ide.github.io/odin2) and
[dust2](https://mrc-ide.github.io/dust2) for simulation.

## Quick example

Define an SIR model as a labelled Petri net, generate odin2 code, and simulate:

```r
library(algebraicodin)

# 1. Define an SIR Petri net
# Transition names (inf, rec) become the odin2 rate parameter names
sir <- labelled_petri_net(
  c("S", "I", "R"),
  inf = c("S", "I") %=>% c("I", "I"),
  rec = "I" %=>% "R"
)

# 2. Generate odin2 code
code <- pn_to_odin(sir, type = "ode")
cat(code)

# 3. Simulate with odin2/dust2
# Parameters: rate constants named after transitions (inf, rec)
# plus initial conditions (S0, I0, R0)
gen <- odin2::odin(code)
sys <- dust2::dust_system_create(gen,
  list(inf = 0.3, rec = 0.1, S0 = 990, I0 = 10, R0 = 0),
  n_particles = 1)
dust2::dust_system_set_state_initial(sys)
t <- seq(0, 100)
y <- dust2::dust_system_simulate(sys, t)

# 4. Plot results (y is [species x time] for n_particles = 1)
sol <- data.frame(time = t, S = y[1, ], I = y[2, ], R = y[3, ])
plot_trajectory(sol, title = "SIR model")
```

## Features

- **Petri nets** — Define compartmental models as bipartite graphs with
  `labelled_petri_net()`, `ReactionNet`, and the `%=>%` transition operator.
- **Stock-flow diagrams** — Build models with explicit stocks, flows, auxiliary
  variables, and sum variables via `stock_and_flow()` and `flow()`.
- **Open systems & UWD composition** — Expose boundary species, then compose
  models with undirected wiring diagrams: `oapply()`, `%o%`, or `compose()`.
- **Stratification** — Typed Petri nets and typed products (`%x%`) for
  age, risk, or geographic strata via `infectious_ontology()` and `stratify()`.
- **odin2 code generation** — Generate ODE, DDE, stochastic, and discrete-time
  odin2 code from Petri nets (`pn_to_odin()`) or stock-flow models (`sf_to_odin()`).
- **Compositional dynamics** — `ResourceSharer` and `Machine` open dynamical
  systems for simulation without code generation.
- **Bayesian fitting** — Observation models, dust2 particle filtering, and MCMC
  via monty with `fit_mcmc()`.
- **Epi building blocks** — Pre-built components for infection, recovery,
  vaccination, waning, and more via `epi_dict()` and `compose_epi()`.
- **Array-based codegen** — Efficient array-indexed odin2 code for large
  stratified models.
- **Spatial models** — Multi-patch dynamics with `sf_spatial()`.

## Vignettes

| Vignette | Description |
|----------|-------------|
| [Getting started](https://catrgory.github.io/algebraicodin/articles/algebraicodin.html) | Core workflow: Petri net → odin2 → simulate |
| [Compositional modeling](https://catrgory.github.io/algebraicodin/articles/composition.html) | Open systems and UWD composition |
| [Stock-and-flow models](https://catrgory.github.io/algebraicodin/articles/stockflow.html) | Stock-flow diagram interface |
| [Compositional dynamics](https://catrgory.github.io/algebraicodin/articles/dynamics.html) | ResourceSharers and Machines |
| [Stratification](https://catrgory.github.io/algebraicodin/articles/stratification.html) | Typed products for age/risk/geo strata |
| [Stock-flow stratification](https://catrgory.github.io/algebraicodin/articles/sf-stratification.html) | Stratification for stock-flow models |
| [Model examples](https://catrgory.github.io/algebraicodin/articles/examples.html) | SIR, SEIR, SIS, and more |
| [More examples](https://catrgory.github.io/algebraicodin/articles/more-examples.html) | Additional epidemiological models |
| [Advanced topics](https://catrgory.github.io/algebraicodin/articles/advanced.html) | Morphisms and complex compositions |
| [Model types](https://catrgory.github.io/algebraicodin/articles/model-types.html) | ODE, DDE, stochastic, and discrete |
| [Array codegen](https://catrgory.github.io/algebraicodin/articles/array-codegen.html) | Array-based odin2 code generation |
| [Contact matrices](https://catrgory.github.io/algebraicodin/articles/contact-matrix.html) | Age-structured contact mixing |
| [Spatial models](https://catrgory.github.io/algebraicodin/articles/spatial.html) | Multi-patch spatial dynamics |
| [Fitting](https://catrgory.github.io/algebraicodin/articles/fitting.html) | Bayesian inference with particle MCMC |

## Author

[Simon Frost](https://github.com/sdwfrost)
(ORCID: [0000-0002-5207-9879](https://orcid.org/0000-0002-5207-9879))

## License

MIT © 2024 Simon Frost. See [LICENSE](LICENSE) for details.

## Part of the catRgory ecosystem

algebraicodin is part of the [catRgory](https://github.com/catRgory) project,
bringing applied category theory to R for compositional scientific modeling.

- [acsets](https://github.com/catRgory/acsets) — Attributed C-sets
- [catlab](https://github.com/catRgory/catlab) — Categories, functors, and wiring diagrams
- **algebraicodin** — Compositional epidemiological modeling with odin2
