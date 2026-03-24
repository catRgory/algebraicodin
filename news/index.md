# Changelog

## algebraicodin 0.1.0

### Initial CRAN release

#### Core Petri Net Functionality

- Petri net types: `PetriNet`, `LabelledPetriNet`, `ReactionNet` built
  on the `acsets` ACSet framework.
- Declarative Petri net builder via
  [`labelled_petri_net()`](https://catrgory.github.io/algebraicodin/reference/labelled_petri_net.md)
  with the `%=>%` transition operator.
- Stoichiometry via
  [`transition_matrices()`](https://catrgory.github.io/algebraicodin/reference/transition_matrices.md)
  and pure-R mass-action
  [`vectorfield()`](https://catrgory.github.io/algebraicodin/reference/vectorfield.md)
  compatible with **deSolve**.

#### Stock-and-Flow Diagrams

- Full stock-flow schema (`SchStockFlow`) with stocks, flows, auxiliary
  variables, sum variables, parameters, and typed links.
- Declarative builder
  [`stock_and_flow()`](https://catrgory.github.io/algebraicodin/reference/stock_and_flow.md)
  with
  [`flow()`](https://catrgory.github.io/algebraicodin/reference/flow.md)
  helper.
- Stock-flow vectorfield
  [`sf_vectorfield()`](https://catrgory.github.io/algebraicodin/reference/sf_vectorfield.md)
  and simulation with
  [`simulate_sf()`](https://catrgory.github.io/algebraicodin/reference/simulate_sf.md).
- Stock-flow transition matrices via
  [`sf_transition_matrices()`](https://catrgory.github.io/algebraicodin/reference/sf_transition_matrices.md).

#### Open Systems and Composition

- `Open` Petri nets with exposed legs for compositional assembly.
- UWD-based composition via
  [`oapply()`](https://catrgory.github.io/algebraicodin/reference/oapply.md)
  and name-based
  [`compose_open()`](https://catrgory.github.io/algebraicodin/reference/compose_open.md).
- Infix operators `%o%` (composition) and `%x%` (stratification).
- Open stock-flow models via
  [`open_stock_flow()`](https://catrgory.github.io/algebraicodin/reference/open_stock_flow.md)
  and
  [`sf_oapply()`](https://catrgory.github.io/algebraicodin/reference/sf_oapply.md).
- Bridge between `Open` and
  [`catlab::StructuredCospan`](https://catrgory.github.io/catlab/reference/StructuredCospan.html)
  via
  [`as_structured_cospan()`](https://catrgory.github.io/algebraicodin/reference/as_structured_cospan.md)
  /
  [`as_open()`](https://catrgory.github.io/algebraicodin/reference/as_open.md).

#### Stratification (Typed Products)

- `TypedPetriNet` with morphisms to a type system (ontology).
- [`typed_petri()`](https://catrgory.github.io/algebraicodin/reference/typed_petri.md)
  for assigning species/transition types.
- [`typed_product()`](https://catrgory.github.io/algebraicodin/reference/typed_product.md)
  /
  [`stratify()`](https://catrgory.github.io/algebraicodin/reference/stratify.md)
  for categorical product in the slice category.
- [`add_reflexives()`](https://catrgory.github.io/algebraicodin/reference/add_reflexives.md)
  for identity transitions needed by stratification.
- [`infectious_ontology()`](https://catrgory.github.io/algebraicodin/reference/infectious_ontology.md)
  providing a standard epi type system.
- Stock-flow stratification via
  [`sf_stratify()`](https://catrgory.github.io/algebraicodin/reference/sf_stratify.md),
  [`sf_typed_product()`](https://catrgory.github.io/algebraicodin/reference/sf_typed_product.md),
  and
  [`sf_spatial()`](https://catrgory.github.io/algebraicodin/reference/sf_spatial.md)
  for spatial/contact-structured models.

#### odin2 / dust2 Code Generation

- [`pn_to_odin()`](https://catrgory.github.io/algebraicodin/reference/pn_to_odin.md)
  /
  [`to_odin()`](https://catrgory.github.io/algebraicodin/reference/pn_to_odin.md):
  generate odin2 DSL code from Petri nets supporting ODE, DDE,
  stochastic (Binomial), and discrete (Euler) model types.
- [`sf_to_odin()`](https://catrgory.github.io/algebraicodin/reference/sf_to_odin.md):
  generate odin2 code from stock-flow models.
- [`sf_to_odin_array()`](https://catrgory.github.io/algebraicodin/reference/sf_to_odin_array.md):
  compact array-based code generation for stratified stock-flow models.
- System constructors
  [`pn_to_odin_system()`](https://catrgory.github.io/algebraicodin/reference/pn_to_odin_system.md),
  [`sf_to_odin_system()`](https://catrgory.github.io/algebraicodin/reference/sf_to_odin_system.md),
  [`sf_to_odin_array_system()`](https://catrgory.github.io/algebraicodin/reference/sf_to_odin_array_system.md)
  returning callable odin2 generators.
- Observation/comparison blocks via
  [`observe()`](https://catrgory.github.io/algebraicodin/reference/observe.md)
  for Bayesian fitting.

#### Compositional Dynamics

- `ResourceSharer` for undirected open dynamical systems
  (UWD-composable).
- `Machine` for directed open dynamical systems (DWD-composable).
- Continuous, discrete, and delay system types.
- Composition via
  [`oapply_dynam()`](https://catrgory.github.io/algebraicodin/reference/oapply_dynam.md)
  and
  [`oapply_machine()`](https://catrgory.github.io/algebraicodin/reference/oapply_machine.md).
- [`petri_to_continuous()`](https://catrgory.github.io/algebraicodin/reference/petri_to_continuous.md),
  [`petri_to_discrete()`](https://catrgory.github.io/algebraicodin/reference/petri_to_discrete.md),
  [`petri_to_delay()`](https://catrgory.github.io/algebraicodin/reference/petri_to_delay.md)
  converters.
- [`euler_approx()`](https://catrgory.github.io/algebraicodin/reference/euler_approx.md)
  for continuous-to-discrete conversion.
- [`rs_to_odin()`](https://catrgory.github.io/algebraicodin/reference/rs_to_odin.md)
  and
  [`machine_to_odin()`](https://catrgory.github.io/algebraicodin/reference/machine_to_odin.md)
  for code generation from dynamical systems.

#### Epidemiology Building Blocks

- [`spontaneous_petri()`](https://catrgory.github.io/algebraicodin/reference/spontaneous_petri.md)
  and
  [`exposure_petri()`](https://catrgory.github.io/algebraicodin/reference/exposure_petri.md)
  for common epi transitions.
- [`epi_dict()`](https://catrgory.github.io/algebraicodin/reference/epi_dict.md)
  pre-built dictionary: infection, recovery, vaccination, waning, death,
  exposure, progression, hospitalization.
- [`compose_epi()`](https://catrgory.github.io/algebraicodin/reference/compose_epi.md)
  for assembling epi models from UWDs and dictionaries.

#### Bayesian Fitting

- [`fit_mcmc()`](https://catrgory.github.io/algebraicodin/reference/fit_mcmc.md)
  for MCMC fitting via **monty** and **dust2**.
- Prior specification helpers:
  [`prior_uniform()`](https://catrgory.github.io/algebraicodin/reference/prior_uniform.md),
  [`prior_normal()`](https://catrgory.github.io/algebraicodin/reference/prior_normal.md),
  [`prior_lognormal()`](https://catrgory.github.io/algebraicodin/reference/prior_lognormal.md),
  [`prior_exp()`](https://catrgory.github.io/algebraicodin/reference/prior_exp.md),
  and
  [`build_prior()`](https://catrgory.github.io/algebraicodin/reference/build_prior.md).
- Diagnostic plots:
  [`plot_traces()`](https://catrgory.github.io/algebraicodin/reference/plot_traces.md),
  [`plot_posterior()`](https://catrgory.github.io/algebraicodin/reference/plot_posterior.md),
  [`plot_fit()`](https://catrgory.github.io/algebraicodin/reference/plot_fit.md).

#### Visualization

- [`plot_petri()`](https://catrgory.github.io/algebraicodin/reference/plot_petri.md),
  [`plot_stock_flow()`](https://catrgory.github.io/algebraicodin/reference/plot_stock_flow.md),
  [`plot_uwd()`](https://catrgory.github.io/algebraicodin/reference/plot_uwd.md)
  using **DiagrammeR**.
- [`plot_trajectory()`](https://catrgory.github.io/algebraicodin/reference/plot_trajectory.md),
  [`plot_comparison()`](https://catrgory.github.io/algebraicodin/reference/plot_comparison.md),
  [`plot_diff()`](https://catrgory.github.io/algebraicodin/reference/plot_diff.md)
  for simulation output.

#### Dependencies

- Imports: **acsets**, **catlab**, **cli**, **rlang**, **S7**.
- Suggests: **odin2**, **dust2**, **deSolve**, **monty**,
  **DiagrammeR**, **ggplot2**, **tidyr**.
