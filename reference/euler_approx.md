# Convert a continuous system to discrete via Euler's method

Convert a continuous system to discrete via Euler's method

## Usage

``` r
euler_approx(rs, h = NULL)
```

## Arguments

- rs:

  A ResourceSharer or Machine with system_type = "continuous"

- h:

  Step size (numeric). If NULL, step size is appended to parameters.

## Value

A discrete system of the same class
