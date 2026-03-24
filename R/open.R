# Open Petri nets and composition -------------------------------------------
# Port of AlgebraicPetri.jl's Open and oapply

#' Open Petri net: a Petri net with exposed "legs" (species accessible
#' from outside for composition).
#' @param pn A Petri net ACSet
#' @param legs List of integer vectors mapping outer junctions to species IDs
#' @examples
#' sir <- labelled_petri_net(
#'   c("S", "I", "R"),
#'   inf = c("S", "I") %=>% c("I", "I"),
#'   rec = "I" %=>% "R"
#' )
#' # Expose all species as a single leg
#' open_sir <- Open(sir)
#' species_names(apex(open_sir)) # c("S", "I", "R")
#'
#' # Expose specific species
#' open_si <- Open(sir, legs = list(c(1L, 2L)))
#' @export
Open <- S7::new_class("Open",
  properties = list(
    apex = acsets::ACSet,   # The underlying Petri net
    legs = S7::class_list   # List of integer vectors: each leg maps
                             # outer junction → species ID in apex
  ),
  constructor = function(pn, legs = NULL) {
    if (is.null(legs)) {
      # Default: single leg exposing all species
      all_s <- acsets::parts(pn, "S")
      legs <- list(all_s)
    }
    S7::new_object(S7::S7_object(), apex = pn, legs = legs)
  }
)

#' Compose open Petri nets over a UWD
#'
#' The UWD specifies how boxes (= open PNs) connect via junctions.
#' Each box's ports connect to junctions; shared junctions mean species
#' are identified (glued together).
#'
#' @param w A UWD ACSet
#' @param components List of Open Petri nets, one per box
#' @returns An Open Petri net (the composed system)
#' @examples
#' # Build SIR by composing infection and recovery
#' infection <- exposure_petri("S", "I", "I", "inf")
#' recovery <- spontaneous_petri("I", "R", "rec")
#'
#' # Create a UWD with shared junction "I"
#' w <- catlab::uwd(
#'   outer = c("S", "I", "R"),
#'   infection = c("S", "I"),
#'   recovery = c("I", "R")
#' )
#' sir <- oapply(w, list(infection, recovery))
#' species_names(apex(sir)) # c("S", "I", "R")
#' @export
oapply <- function(w, components) {
  n_boxes <- acsets::nparts(w, "Box")
  n_junctions <- acsets::nparts(w, "Junction")
  if (length(components) != n_boxes) {
    cli::cli_abort("Number of components ({length(components)}) must match boxes ({n_boxes})")
  }

  # Start with empty labelled Petri net
  result <- LabelledPetriNet()

  # Add one species per junction
  junc_species <- integer(n_junctions)
  for (j in seq_len(n_junctions)) {
    jname <- if (acsets::has_subpart(w, "name")) {
      acsets::subpart(w, j, "name")
    } else paste0("J", j)
    junc_species[j] <- acsets::add_part(result, "S", sname = jname)
  }

  # For each box, add its transitions and internal species, mapping
  # leg species to junction species
  for (b in seq_len(n_boxes)) {
    comp <- components[[b]]
    pn <- comp@apex

    # Find which ports this box has, and their junctions
    box_ports <- acsets::incident(w, b, "box")
    port_junctions <- acsets::subpart(w, box_ports, "junction")

    # Build species map: species in component → species in result
    # Leg species map to the corresponding junction species
    n_comp_species <- acsets::nparts(pn, "S")
    species_map <- integer(n_comp_species)

    # The leg (first leg = default) maps species indices to the ports
    leg <- comp@legs[[1]]
    if (length(leg) != length(port_junctions)) {
      cli::cli_abort(c(
        "Leg length ({length(leg)}) doesn't match port count ({length(port_junctions)})",
        i = "Box {b}"
      ))
    }
    for (k in seq_along(leg)) {
      species_map[leg[k]] <- junc_species[port_junctions[k]]
    }

    # Any unmapped species (internal) get new species in result
    for (s in seq_len(n_comp_species)) {
      if (species_map[s] == 0L) {
        sname <- if (acsets::has_subpart(pn, "sname")) {
          acsets::subpart(pn, s, "sname")
        } else paste0("S", s, "_box", b)
        species_map[s] <- acsets::add_part(result, "S", sname = sname)
      }
    }

    # Add transitions
    n_trans <- acsets::nparts(pn, "T")
    for (t in seq_len(n_trans)) {
      tname <- if (acsets::has_subpart(pn, "tname")) {
        acsets::subpart(pn, t, "tname")
      } else paste0("T", t, "_box", b)
      new_t <- acsets::add_part(result, "T", tname = tname)

      # Input arcs
      input_arcs <- acsets::incident(pn, t, "it")
      for (ia in input_arcs) {
        old_s <- acsets::subpart(pn, ia, "is")
        acsets::add_part(result, "I", is = species_map[old_s], it = new_t)
      }

      # Output arcs
      output_arcs <- acsets::incident(pn, t, "ot")
      for (oa in output_arcs) {
        old_s <- acsets::subpart(pn, oa, "os")
        acsets::add_part(result, "O", os = species_map[old_s], ot = new_t)
      }
    }
  }

  # Build outer legs from outer ports
  n_op <- acsets::nparts(w, "OuterPort")
  outer_leg <- integer(n_op)
  for (op in seq_len(n_op)) {
    oj <- acsets::subpart(w, op, "outer_junction")
    outer_leg[op] <- junc_species[oj]
  }

  Open(result, legs = list(outer_leg))
}

#' Extract the closed Petri net from an Open Petri net
#' @param open_pn An Open Petri net
#' @returns A Petri net ACSet
#' @examples
#' pn <- labelled_petri_net(c("S", "I"), inf = "S" %=>% "I")
#' open_pn <- Open(pn)
#' closed <- apex(open_pn)
#' species_names(closed) # c("S", "I")
#' @export
apex <- function(open_pn) open_pn@apex

# Conversion between Open and catlab::StructuredCospan -------------------------

#' Convert an Open Petri net to a catlab StructuredCospan
#'
#' Each integer-vector leg becomes an ACSetTransformation from a discrete
#' Petri net foot (with parts only in "S") to the apex.
#'
#' @param open_pn An Open Petri net
#' @returns A catlab::StructuredCospan with interface_ob = "S"
#' @examples
#' sir <- labelled_petri_net(
#'   c("S", "I", "R"),
#'   inf = c("S", "I") %=>% c("I", "I"),
#'   rec = "I" %=>% "R"
#' )
#' open_sir <- Open(sir)
#' sc <- as_structured_cospan(open_sir)
#' # Round-trip back to Open
#' open_again <- as_open(sc)
#' species_names(apex(open_again))
#' @export
as_structured_cospan <- function(open_pn) {
  do.call(catlab::open_acset, c(list(open_pn@apex, "S"), open_pn@legs))
}

#' Convert a catlab StructuredCospan back to an Open Petri net
#'
#' Extracts the apex and recovers each leg as an integer vector of species
#' indices from the ACSetTransformation components.
#'
#' @param sc A catlab::StructuredCospan
#' @returns An Open Petri net
#' @examples
#' sir <- labelled_petri_net(
#'   c("S", "I", "R"),
#'   inf = c("S", "I") %=>% c("I", "I"),
#'   rec = "I" %=>% "R"
#' )
#' sc <- open_petri_net(sir, 1:3)
#' open_sir <- as_open(sc)
#' species_names(apex(open_sir))
#' @export
as_open <- function(sc) {
  legs <- lapply(sc@legs, function(leg) {
    leg@components[[sc@interface_ob]]
  })
  Open(sc@apex, legs = legs)
}

#' Open a Petri net at specified species, returning a StructuredCospan
#'
#' Convenience wrapper around catlab::open_acset with interface_ob = "S".
#'
#' @param pn A Petri net ACSet
#' @param ... Integer vectors, each specifying a leg (species indices)
#' @returns A catlab::StructuredCospan with interface_ob = "S"
#' @examples
#' sir <- labelled_petri_net(
#'   c("S", "I", "R"),
#'   inf = c("S", "I") %=>% c("I", "I"),
#'   rec = "I" %=>% "R"
#' )
#' sc <- open_petri_net(sir, c(1L, 2L), c(2L, 3L))
#' @export
open_petri_net <- function(pn, ...) {
  catlab::open_acset(pn, "S", ...)
}
