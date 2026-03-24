# Plot a Petri net as a bipartite graph via DiagrammeR

Species are shown as blue circles, transitions as orange boxes. Accepts
LabelledPetriNet, TypedPetriNet, or Open Petri nets.

## Usage

``` r
plot_petri(pn)
```

## Arguments

- pn:

  A Petri net ACSet, TypedPetriNet, or Open Petri net

## Value

A DiagrammeR `htmlwidget`

## Examples

``` r
sir <- labelled_petri_net(
  c("S", "I", "R"),
  inf = c("S", "I") %=>% c("I", "I"),
  rec = "I" %=>% "R"
)
if (requireNamespace("DiagrammeR", quietly = TRUE)) {
  plot_petri(sir)
}

{"x":{"diagram":"digraph PetriNet {\n  rankdir=LR;\n  S_1 [label=\"S\" shape=circle style=filled fillcolor=lightskyblue];\n  S_2 [label=\"I\" shape=circle style=filled fillcolor=lightskyblue];\n  S_3 [label=\"R\" shape=circle style=filled fillcolor=lightskyblue];\n  T_1 [label=\"inf\" shape=box style=filled fillcolor=lightsalmon];\n  T_2 [label=\"rec\" shape=box style=filled fillcolor=lightsalmon];\n  S_1 -> T_1;\n  S_2 -> T_1;\n  S_2 -> T_2;\n  T_1 -> S_2;\n  T_1 -> S_2;\n  T_2 -> S_3;\n}","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}
```
