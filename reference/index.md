# Package index

## Petri Nets

Petri net constructors and accessors

- [`PetriNet()`](https://catrgory.github.io/algebraicodin/reference/PetriNet.md)
  : PetriNet ACSet type (unlabelled)
- [`LabelledPetriNet()`](https://catrgory.github.io/algebraicodin/reference/LabelledPetriNet.md)
  : LabelledPetriNet ACSet type
- [`ReactionNet()`](https://catrgory.github.io/algebraicodin/reference/ReactionNet.md)
  : ReactionNet ACSet type (labelled with rates and concentrations)
- [`TypedPetriNet()`](https://catrgory.github.io/algebraicodin/reference/TypedPetriNet.md)
  : Typed Petri net: a Petri net with a morphism to a type system.
- [`SchPetriNet`](https://catrgory.github.io/algebraicodin/reference/SchPetriNet.md)
  : Petri net schema (unlabelled)
- [`SchLabelledPetriNet`](https://catrgory.github.io/algebraicodin/reference/SchLabelledPetriNet.md)
  : Labelled Petri net schema
- [`SchReactionNet`](https://catrgory.github.io/algebraicodin/reference/SchReactionNet.md)
  : Reaction net schema (labelled + rates + concentrations)
- [`labelled_petri_net()`](https://catrgory.github.io/algebraicodin/reference/labelled_petri_net.md)
  : Build a LabelledPetriNet from species names and transition specs
- [`open_petri_net()`](https://catrgory.github.io/algebraicodin/reference/open_petri_net.md)
  : Open a Petri net at specified species, returning a StructuredCospan
- [`species_names()`](https://catrgory.github.io/algebraicodin/reference/species_names.md)
  : Get species names from a Petri net
- [`transition_names()`](https://catrgory.github.io/algebraicodin/reference/transition_names.md)
  : Get transition names from a Petri net
- [`typed_petri()`](https://catrgory.github.io/algebraicodin/reference/typed_petri.md)
  : Create a typed Petri net from type assignments

## Stock-Flow Diagrams

Stock-and-flow model constructors and accessors

- [`StockFlow0()`](https://catrgory.github.io/algebraicodin/reference/StockFlow0.md)
  : Create StockFlow0 ACSet type (for composition feet)
- [`StockFlowModel()`](https://catrgory.github.io/algebraicodin/reference/StockFlowModel.md)
  : A Stock-and-Flow model
- [`StockFlowType()`](https://catrgory.github.io/algebraicodin/reference/StockFlowType.md)
  : Create StockFlow ACSet type
- [`TypedStockFlow()`](https://catrgory.github.io/algebraicodin/reference/TypedStockFlow.md)
  : Typed stock-flow model for stratification
- [`SchStockFlow`](https://catrgory.github.io/algebraicodin/reference/SchStockFlow.md)
  : Full Stock-and-Flow schema
- [`SchStockFlow0`](https://catrgory.github.io/algebraicodin/reference/SchStockFlow0.md)
  : Stock-and-Flow base schema (stocks + sum variables only)
- [`stock_and_flow()`](https://catrgory.github.io/algebraicodin/reference/stock_and_flow.md)
  : Build a stock-and-flow model using a declarative interface
- [`flow()`](https://catrgory.github.io/algebraicodin/reference/flow.md)
  : Create a flow specification
- [`open_stock_flow()`](https://catrgory.github.io/algebraicodin/reference/open_stock_flow.md)
  : Create an open stock-flow model
- [`add_reflexives()`](https://catrgory.github.io/algebraicodin/reference/add_reflexives.md)
  : Add reflexive (self-loop) transitions for stratification
- [`sf_add_reflexives()`](https://catrgory.github.io/algebraicodin/reference/sf_add_reflexives.md)
  : Add reflexive (identity) flows for stock-flow stratification
- [`sf_fnames()`](https://catrgory.github.io/algebraicodin/reference/sf_fnames.md)
  : Get flow names
- [`sf_pnames()`](https://catrgory.github.io/algebraicodin/reference/sf_pnames.md)
  : Get parameter names
- [`sf_snames()`](https://catrgory.github.io/algebraicodin/reference/sf_snames.md)
  : Get stock names
- [`sf_svnames()`](https://catrgory.github.io/algebraicodin/reference/sf_svnames.md)
  : Get sum variable names
- [`sf_vnames()`](https://catrgory.github.io/algebraicodin/reference/sf_vnames.md)
  : Get variable names

## Open Systems & Composition

Compose models via open systems and wiring diagrams

- [`Open()`](https://catrgory.github.io/algebraicodin/reference/Open.md)
  : Open Petri net: a Petri net with exposed "legs" (species accessible
  from outside for composition).
- [`as_open()`](https://catrgory.github.io/algebraicodin/reference/as_open.md)
  : Convert a catlab StructuredCospan back to an Open Petri net
- [`as_structured_cospan()`](https://catrgory.github.io/algebraicodin/reference/as_structured_cospan.md)
  : Convert an Open Petri net to a catlab StructuredCospan
- [`compose()`](https://catrgory.github.io/algebraicodin/reference/compose.md)
  : Compose multiple Open Petri nets
- [`compose_epi()`](https://catrgory.github.io/algebraicodin/reference/compose_epi.md)
  : Compose an epidemiological model from a UWD and dictionary
- [`compose_open()`](https://catrgory.github.io/algebraicodin/reference/compose_open.md)
  : Compose two Open Petri nets by matching species names
- [`oapply()`](https://catrgory.github.io/algebraicodin/reference/oapply.md)
  : Compose open Petri nets over a UWD
- [`oapply_dynam()`](https://catrgory.github.io/algebraicodin/reference/oapply_dynam.md)
  : Compose ResourceSharers via an undirected wiring diagram
- [`oapply_machine()`](https://catrgory.github.io/algebraicodin/reference/oapply_machine.md)
  : Compose Machines via a directed wiring diagram
- [`sf_oapply()`](https://catrgory.github.io/algebraicodin/reference/sf_oapply.md)
  : Compose open stock-flow models via a UWD
- [`nports()`](https://catrgory.github.io/algebraicodin/reference/nports.md)
  : Number of ports
- [`apex()`](https://catrgory.github.io/algebraicodin/reference/apex.md)
  : Extract the closed Petri net from an Open Petri net
- [`exposed_states()`](https://catrgory.github.io/algebraicodin/reference/exposed_states.md)
  : Get exposed state values at ports
- [`dwd()`](https://catrgory.github.io/algebraicodin/reference/dwd.md) :
  Create a directed wiring diagram for Machine composition
- [`epi_dict()`](https://catrgory.github.io/algebraicodin/reference/epi_dict.md)
  : Pre-built epidemiology dictionary for SIR-family models
- [`` `%=>%` ``](https://catrgory.github.io/algebraicodin/reference/grapes-equals-greater-than-grapes.md)
  : Define a Petri net transition: inputs to outputs
- [`` `%o%` ``](https://catrgory.github.io/algebraicodin/reference/grapes-o-grapes.md)
  : Compose two Open Petri nets by matching species names
- [`` `%x%` ``](https://catrgory.github.io/algebraicodin/reference/grapes-x-grapes.md)
  : Typed product (stratification) operator

## Stratification

Typed products for model stratification

- [`typed_product()`](https://catrgory.github.io/algebraicodin/reference/typed_product.md)
  : Compute the typed product (stratification) of typed Petri nets
- [`sf_typed()`](https://catrgory.github.io/algebraicodin/reference/sf_typed.md)
  : Create a typed stock-flow model
- [`sf_typed_product()`](https://catrgory.github.io/algebraicodin/reference/sf_typed_product.md)
  : Compute typed product of two typed stock-flow models
  (stratification)
- [`sf_stratify()`](https://catrgory.github.io/algebraicodin/reference/sf_stratify.md)
  : User-friendly stock-flow stratification
- [`stratify()`](https://catrgory.github.io/algebraicodin/reference/stratify.md)
  : User-friendly stratification wrapper
- [`infectious_ontology()`](https://catrgory.github.io/algebraicodin/reference/infectious_ontology.md)
  : Standard infectious disease ontology
- [`sf_infectious_type()`](https://catrgory.github.io/algebraicodin/reference/sf_infectious_type.md)
  : SF type template for simple infectious disease models
- [`create_filter()`](https://catrgory.github.io/algebraicodin/reference/create_filter.md)
  : Create a dust2 particle filter from a model
- [`flatten()`](https://catrgory.github.io/algebraicodin/reference/flatten.md)
  : Extract the underlying Petri net from a TypedPetriNet

## Code Generation

Generate odin2 code from models

- [`pn_to_odin()`](https://catrgory.github.io/algebraicodin/reference/pn_to_odin.md)
  [`to_odin()`](https://catrgory.github.io/algebraicodin/reference/pn_to_odin.md)
  : Generate odin2 code from a Petri net
- [`pn_to_odin_system()`](https://catrgory.github.io/algebraicodin/reference/pn_to_odin_system.md)
  [`to_odin_system()`](https://catrgory.github.io/algebraicodin/reference/pn_to_odin_system.md)
  : Generate odin2 code and return it as a callable odin_system
  generator
- [`sf_to_odin()`](https://catrgory.github.io/algebraicodin/reference/sf_to_odin.md)
  : Generate odin2 code from a stock-flow model
- [`sf_to_odin_array()`](https://catrgory.github.io/algebraicodin/reference/sf_to_odin_array.md)
  : Generate array-based odin2 code from a stratified stock-flow model
- [`sf_to_odin_array_system()`](https://catrgory.github.io/algebraicodin/reference/sf_to_odin_array_system.md)
  : Generate an array-based odin2 system from a stratified model
- [`sf_to_odin_system()`](https://catrgory.github.io/algebraicodin/reference/sf_to_odin_system.md)
  : Generate an odin2 system from a stock-flow model
- [`rs_to_odin()`](https://catrgory.github.io/algebraicodin/reference/rs_to_odin.md)
  : Generate odin2 code from a ResourceSharer
- [`machine_to_odin()`](https://catrgory.github.io/algebraicodin/reference/machine_to_odin.md)
  : Generate odin2 code from a Machine

## Dynamics

Vector fields, transition matrices, and model conversions

- [`vectorfield()`](https://catrgory.github.io/algebraicodin/reference/vectorfield.md)
  : Generate a pure-R mass-action vector field function
- [`sf_vectorfield()`](https://catrgory.github.io/algebraicodin/reference/sf_vectorfield.md)
  : Generate a vectorfield function from a stock-flow model
- [`transition_matrices()`](https://catrgory.github.io/algebraicodin/reference/transition_matrices.md)
  : Extract input/output stoichiometry matrices
- [`sf_transition_matrices()`](https://catrgory.github.io/algebraicodin/reference/sf_transition_matrices.md)
  : Compute inflow/outflow transition matrices for a stock-flow model
- [`petri_to_continuous()`](https://catrgory.github.io/algebraicodin/reference/petri_to_continuous.md)
  : Convert a Petri net to a ContinuousResourceSharer (ODE)
- [`petri_to_delay()`](https://catrgory.github.io/algebraicodin/reference/petri_to_delay.md)
  : Convert a Petri net to a DelayResourceSharer (DDE)
- [`petri_to_discrete()`](https://catrgory.github.io/algebraicodin/reference/petri_to_discrete.md)
  : Convert a Petri net to a DiscreteResourceSharer (stochastic)
- [`petri_to_sf()`](https://catrgory.github.io/algebraicodin/reference/petri_to_sf.md)
  : Convert a Petri net to a stock-flow model
- [`sf_to_petri()`](https://catrgory.github.io/algebraicodin/reference/sf_to_petri.md)
  : Convert a stock-flow model to a labelled Petri net
- [`sf_to_resource_sharer()`](https://catrgory.github.io/algebraicodin/reference/sf_to_resource_sharer.md)
  : Convert a stock-flow model to a continuous ResourceSharer
- [`sf_spatial()`](https://catrgory.github.io/algebraicodin/reference/sf_spatial.md)
  : Create a spatial (multi-patch) model via stratification
- [`euler_approx()`](https://catrgory.github.io/algebraicodin/reference/euler_approx.md)
  : Convert a continuous system to discrete via Euler's method
- [`eval_dynamics()`](https://catrgory.github.io/algebraicodin/reference/eval_dynamics.md)
  : Evaluate dynamics of a ResourceSharer
- [`eval_readout()`](https://catrgory.github.io/algebraicodin/reference/eval_readout.md)
  : Evaluate readout of a Machine

## Machines & Resource Sharers

Continuous and discrete dynamical systems

- [`Machine()`](https://catrgory.github.io/algebraicodin/reference/Machine.md)
  : A machine: open dynamical system with directed inputs and outputs.
- [`ResourceSharer()`](https://catrgory.github.io/algebraicodin/reference/ResourceSharer.md)
  : A resource sharer: open dynamical system with symmetric ports.
- [`continuous_machine()`](https://catrgory.github.io/algebraicodin/reference/continuous_machine.md)
  : Create a ContinuousMachine
- [`continuous_resource_sharer()`](https://catrgory.github.io/algebraicodin/reference/continuous_resource_sharer.md)
  : Create a ContinuousResourceSharer
- [`delay_machine()`](https://catrgory.github.io/algebraicodin/reference/delay_machine.md)
  : Create a DelayMachine
- [`delay_resource_sharer()`](https://catrgory.github.io/algebraicodin/reference/delay_resource_sharer.md)
  : Create a DelayResourceSharer from a DDE dynamics function
- [`discrete_resource_sharer()`](https://catrgory.github.io/algebraicodin/reference/discrete_resource_sharer.md)
  : Create a DiscreteResourceSharer
- [`simulate_machine()`](https://catrgory.github.io/algebraicodin/reference/simulate_machine.md)
  : Simulate a Machine with external inputs
- [`simulate_rs()`](https://catrgory.github.io/algebraicodin/reference/simulate_rs.md)
  : Simulate a continuous ResourceSharer using deSolve
- [`simulate_rs_discrete()`](https://catrgory.github.io/algebraicodin/reference/simulate_rs_discrete.md)
  : Simulate a discrete ResourceSharer
- [`simulate_sf()`](https://catrgory.github.io/algebraicodin/reference/simulate_sf.md)
  : Simulate a stock-flow model using deSolve
- [`observe()`](https://catrgory.github.io/algebraicodin/reference/observe.md)
  : Create an observation specification

## Epidemic Building Blocks

Pre-built Petri net components for epidemiological models

- [`spontaneous_petri()`](https://catrgory.github.io/algebraicodin/reference/spontaneous_petri.md)
  : Spontaneous transition: A -\> B
- [`exposure_petri()`](https://catrgory.github.io/algebraicodin/reference/exposure_petri.md)
  : Exposure/contact transition: A + B -\> C + B (B catalyzes A -\> C)

## Model Fitting

Bayesian inference with MCMC

- [`fit_mcmc()`](https://catrgory.github.io/algebraicodin/reference/fit_mcmc.md)
  : Fit a model to data using MCMC
- [`build_prior()`](https://catrgory.github.io/algebraicodin/reference/build_prior.md)
  : Build a prior model from specifications
- [`prior_exp()`](https://catrgory.github.io/algebraicodin/reference/prior_exp.md)
  : Exponential prior
- [`prior_lognormal()`](https://catrgory.github.io/algebraicodin/reference/prior_lognormal.md)
  : Log-normal prior
- [`prior_normal()`](https://catrgory.github.io/algebraicodin/reference/prior_normal.md)
  : Normal prior
- [`prior_uniform()`](https://catrgory.github.io/algebraicodin/reference/prior_uniform.md)
  : Uniform prior

## Visualization

Plotting models and results

- [`plot_petri()`](https://catrgory.github.io/algebraicodin/reference/plot_petri.md)
  : Plot a Petri net as a bipartite graph via DiagrammeR
- [`plot_stock_flow()`](https://catrgory.github.io/algebraicodin/reference/plot_stock_flow.md)
  : Plot a stock-flow diagram using DiagrammeR
- [`plot_uwd()`](https://catrgory.github.io/algebraicodin/reference/plot_uwd.md)
  : Plot an undirected wiring diagram via DiagrammeR
- [`plot_comparison()`](https://catrgory.github.io/algebraicodin/reference/plot_comparison.md)
  : Plot a comparison of two ODE solutions (e.g., R vs Julia)
- [`plot_diff()`](https://catrgory.github.io/algebraicodin/reference/plot_diff.md)
  : Plot absolute differences between two ODE solutions
- [`plot_fit()`](https://catrgory.github.io/algebraicodin/reference/plot_fit.md)
  : Plot model fit against data
- [`plot_posterior()`](https://catrgory.github.io/algebraicodin/reference/plot_posterior.md)
  : Plot posterior density estimates
- [`plot_traces()`](https://catrgory.github.io/algebraicodin/reference/plot_traces.md)
  : Plot MCMC trace plots
- [`plot_trajectory()`](https://catrgory.github.io/algebraicodin/reference/plot_trajectory.md)
  : Plot ODE solution trajectories with ggplot2
