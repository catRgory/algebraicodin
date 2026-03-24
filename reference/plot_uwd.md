# Plot an undirected wiring diagram via DiagrammeR

Boxes are yellow rectangles, junctions are grey circles, outer ports are
small black diamonds.

## Usage

``` r
plot_uwd(w)
```

## Arguments

- w:

  A UWD ACSet

## Value

A DiagrammeR `htmlwidget`

## Examples

``` r
w <- uwd(
  outer = c("s", "i", "r"),
  infection = c("s", "i"),
  recovery = c("i", "r")
)
#> Error in uwd(outer = c("s", "i", "r"), infection = c("s", "i"), recovery = c("i",     "r")): could not find function "uwd"
if (requireNamespace("DiagrammeR", quietly = TRUE)) {
  plot_uwd(w)
}
#> Error: object 'w' not found
```
