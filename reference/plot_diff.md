# Plot absolute differences between two ODE solutions

Plot absolute differences between two ODE solutions

## Usage

``` r
plot_diff(sol1, sol2, vars = NULL, title = "Absolute Difference")
```

## Arguments

- sol1:

  First solution (data.frame)

- sol2:

  Second solution (data.frame)

- vars:

  Character vector of variable names

- title:

  Optional plot title

## Value

A ggplot2 object

## Examples

``` r
t <- seq(0, 10, 0.1)
sol1 <- data.frame(time = t, S = 990 * exp(-0.3 * t), I = 10 * exp(0.1 * t))
sol2 <- data.frame(time = t, S = 989 * exp(-0.3 * t), I = 11 * exp(0.1 * t))
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_diff(sol1, sol2)
}
```
