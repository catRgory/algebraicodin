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
sir <- labelled_petri_net(
  c("S", "I", "R"),
  inf = c("S", "I") %=>% c("I", "I"),
  rec = "I" %=>% "R"
)

# 2. Open the system at boundary species
sir_open <- open_petri_net(sir, legs = list(c("S", "I"), c("I", "R")))

# 3. Generate odin2 code
code <- pn_to_odin(sir, type = "ode",
  params = list(beta = 0.3, gamma = 0.1),
  initial = list(S = 990, I = 10, R = 0))
cat(code)

# 4. Simulate with odin2/dust2
gen <- odin2::odin(code)
sys <- dust2::dust_system_create(gen, list())
dust2::dust_system_set_state_initial(sys)
t <- seq(0, 100)
y <- dust2::dust_system_simulate(sys, t)

# 5. Plot results
plot_trajectory(y, t, species = c("S", "I", "R"))
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
| [Getting started](vignettes/algebraicodin.Rmd) | Core workflow: Petri net → odin2 → simulate |
| [Compositional modeling](vignettes/composition.Rmd) | Open systems and UWD composition |
| [Stock-and-flow models](vignettes/stockflow.Rmd) | Stock-flow diagram interface |
| [Compositional dynamics](vignettes/dynamics.Rmd) | ResourceSharers and Machines |
| [Stratification](vignettes/stratification.Rmd) | Typed products for age/risk/geo strata |
| [Stock-flow stratification](vignettes/sf-stratification.Rmd) | Stratification for stock-flow models |
| [Model examples](vignettes/examples.Rmd) | SIR, SEIR, SIS, and more |
| [More examples](vignettes/more-examples.Rmd) | Additional epidemiological models |
| [Advanced topics](vignettes/advanced.Rmd) | Morphisms and complex compositions |
| [Model types](vignettes/model-types.Rmd) | ODE, DDE, stochastic, and discrete |
| [Array codegen](vignettes/array-codegen.Rmd) | Array-based odin2 code generation |
| [Contact matrices](vignettes/contact-matrix.Rmd) | Age-structured contact mixing |
| [Spatial models](vignettes/spatial.Rmd) | Multi-patch spatial dynamics |
| [Fitting](vignettes/fitting.Rmd) | Bayesian inference with particle MCMC |

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
