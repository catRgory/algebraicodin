# Plot ODE solution trajectories with ggplot2

Plot ODE solution trajectories with ggplot2

## Usage

``` r
plot_trajectory(sol, vars = NULL, title = NULL)
```

## Arguments

- sol:

  A deSolve solution matrix or data.frame with a `time` column

- vars:

  Character vector of variable names to plot (default: all non-time)

- title:

  Optional plot title

## Value

A ggplot2 object

## Examples

``` r
sol <- data.frame(time = seq(0, 10, 0.1),
                  S = 990 * exp(-0.3 * seq(0, 10, 0.1)),
                  I = 10 * exp(0.1 * seq(0, 10, 0.1)))
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_trajectory(sol, title = "SIR Trajectory")
}
```
