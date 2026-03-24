# Create an observation specification

Specifies how a model variable relates to observed data for use in data
comparison (likelihood computation). Used with the `compare` argument of
code generation functions.

## Usage

``` r
observe(
  dist,
  transition = NULL,
  flow = NULL,
  state = NULL,
  expr = NULL,
  noise = NULL,
  sd = NULL,
  zero_every = 1
)
```

## Arguments

- dist:

  Distribution name (e.g., "Poisson", "Normal", "NegativeBinomial")

- transition:

  Name of a Petri net transition to track incidence

- flow:

  Name of a stock-flow flow to track incidence

- state:

  Name of a state variable to observe directly (prevalence)

- expr:

  Custom expression string for the distribution parameter

- noise:

  Noise parameter name (added to mean to prevent zero lambda)

- sd:

  Standard deviation parameter name (for Normal distribution)

- zero_every:

  Reset interval for incidence accumulators (default: 1)

## Value

An observe_spec object
