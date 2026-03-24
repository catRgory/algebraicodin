# Plot a comparison of two ODE solutions (e.g., R vs Julia)

Plot a comparison of two ODE solutions (e.g., R vs Julia)

## Usage

``` r
plot_comparison(
  sol1,
  sol2,
  vars = NULL,
  labels = c("Model 1", "Model 2"),
  title = NULL
)
```

## Arguments

- sol1:

  First solution (data.frame with time column)

- sol2:

  Second solution (data.frame with time column)

- vars:

  Character vector of variable names to compare

- labels:

  Two-element character vector labeling each solution

- title:

  Optional plot title

## Value

A ggplot2 object

## Examples

``` r
t <- seq(0, 10, 0.1)
sol1 <- data.frame(time = t, S = 990 * exp(-0.3 * t), I = 10 * exp(0.1 * t))
sol2 <- data.frame(time = t, S = 980 * exp(-0.3 * t), I = 20 * exp(0.1 * t))
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_comparison(sol1, sol2, labels = c("Run 1", "Run 2"))
}
```
