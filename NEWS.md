# algebraicodin 0.1.0

## Initial CRAN release

### Core Petri Net Functionality
* Petri net types: `PetriNet`, `LabelledPetriNet`, `ReactionNet` built on the
  `acsets` ACSet framework.
* Declarative Petri net builder via `labelled_petri_net()` with the `%=>%`
  transition operator.
* Stoichiometry via `transition_matrices()` and pure-R mass-action
  `vectorfield()` compatible with **deSolve**.

### Stock-and-Flow Diagrams
* Full stock-flow schema (`SchStockFlow`) with stocks, flows, auxiliary
  variables, sum variables, parameters, and typed links.
* Declarative builder `stock_and_flow()` with `flow()` helper.
* Stock-flow vectorfield `sf_vectorfield()` and simulation with `simulate_sf()`.
* Stock-flow transition matrices via `sf_transition_matrices()`.

### Open Systems and Composition
* `Open` Petri nets with exposed legs for compositional assembly.
* UWD-based composition via `oapply()` and name-based `compose_open()`.
* Infix operators `%o%` (composition) and `%x%` (stratification).
* Open stock-flow models via `open_stock_flow()` and `sf_oapply()`.
* Bridge between `Open` and `catlab::StructuredCospan` via
 `as_structured_cospan()` / `as_open()`.

### Stratification (Typed Products)
* `TypedPetriNet` with morphisms to a type system (ontology).
* `typed_petri()` for assigning species/transition types.
* `typed_product()` / `stratify()` for categorical product in the slice
  category.
* `add_reflexives()` for identity transitions needed by stratification.
* `infectious_ontology()` providing a standard epi type system.
* Stock-flow stratification via `sf_stratify()`, `sf_typed_product()`, and
  `sf_spatial()` for spatial/contact-structured models.

### odin2 / dust2 Code Generation
* `pn_to_odin()` / `to_odin()`: generate odin2 DSL code from Petri nets
  supporting ODE, DDE, stochastic (Binomial), and discrete (Euler) model types.
* `sf_to_odin()`: generate odin2 code from stock-flow models.
* `sf_to_odin_array()`: compact array-based code generation for stratified
  stock-flow models.
* System constructors `pn_to_odin_system()`, `sf_to_odin_system()`,
  `sf_to_odin_array_system()` returning callable odin2 generators.
* Observation/comparison blocks via `observe()` for Bayesian fitting.

### Compositional Dynamics
* `ResourceSharer` for undirected open dynamical systems (UWD-composable).
* `Machine` for directed open dynamical systems (DWD-composable).
* Continuous, discrete, and delay system types.
* Composition via `oapply_dynam()` and `oapply_machine()`.
* `petri_to_continuous()`, `petri_to_discrete()`, `petri_to_delay()` converters.
* `euler_approx()` for continuous-to-discrete conversion.
* `rs_to_odin()` and `machine_to_odin()` for code generation from dynamical
  systems.

### Epidemiology Building Blocks
* `spontaneous_petri()` and `exposure_petri()` for common epi transitions.
* `epi_dict()` pre-built dictionary: infection, recovery, vaccination, waning,
  death, exposure, progression, hospitalization.
* `compose_epi()` for assembling epi models from UWDs and dictionaries.

### Bayesian Fitting
* `fit_mcmc()` for MCMC fitting via **monty** and **dust2**.
* Prior specification helpers: `prior_uniform()`, `prior_normal()`,
  `prior_lognormal()`, `prior_exp()`, and `build_prior()`.
* Diagnostic plots: `plot_traces()`, `plot_posterior()`, `plot_fit()`.

### Visualization
* `plot_petri()`, `plot_stock_flow()`, `plot_uwd()` using **DiagrammeR**.
* `plot_trajectory()`, `plot_comparison()`, `plot_diff()` for simulation output.

### Dependencies
* Imports: **acsets**, **catlab**, **cli**, **rlang**, **S7**.
* Suggests: **odin2**, **dust2**, **deSolve**, **monty**, **DiagrammeR**,
  **ggplot2**, **tidyr**.
