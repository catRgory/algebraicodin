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
w <- catlab::uwd(
  outer = c("s", "i", "r"),
  infection = c("s", "i"),
  recovery = c("i", "r")
)
if (requireNamespace("DiagrammeR", quietly = TRUE)) {
  plot_uwd(w)
}

{"x":{"diagram":"graph UWD {\n  rankdir=LR;\n  J_1 [label=\"s\" shape=circle width=0.3 style=filled fillcolor=gray90];\n  J_2 [label=\"i\" shape=circle width=0.3 style=filled fillcolor=gray90];\n  J_3 [label=\"r\" shape=circle width=0.3 style=filled fillcolor=gray90];\n  B_1 [label=\"Box 1\" shape=box style=filled fillcolor=lightyellow];\n  B_2 [label=\"Box 2\" shape=box style=filled fillcolor=lightyellow];\n  OP_1 [label=\"\" shape=diamond width=0.2 style=filled fillcolor=black];\n  OP_1 -- J_1;\n  OP_2 [label=\"\" shape=diamond width=0.2 style=filled fillcolor=black];\n  OP_2 -- J_2;\n  OP_3 [label=\"\" shape=diamond width=0.2 style=filled fillcolor=black];\n  OP_3 -- J_3;\n  B_1 -- J_1;\n  B_1 -- J_2;\n  B_2 -- J_2;\n  B_2 -- J_3;\n}","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}
```
