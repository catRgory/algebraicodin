# Petri net schemas and constructors ----------------------------------------
# Port of AlgebraicPetri.jl

#' Petri net schema (unlabelled)
#' @export
SchPetriNet <- acsets::BasicSchema(
  obs = c("S", "T", "I", "O"),
  homs = list(
    acsets::hom("is", "I", "S"), acsets::hom("it", "I", "T"),
    acsets::hom("os", "O", "S"), acsets::hom("ot", "O", "T")
  )
)

#' Labelled Petri net schema
#' @export
SchLabelledPetriNet <- acsets::BasicSchema(
  obs = c("S", "T", "I", "O"),
  homs = list(
    acsets::hom("is", "I", "S"), acsets::hom("it", "I", "T"),
    acsets::hom("os", "O", "S"), acsets::hom("ot", "O", "T")
  ),
  attrtypes = "Name",
  attrs = list(
    acsets::attr_spec("sname", "S", "Name"),
    acsets::attr_spec("tname", "T", "Name")
  )
)

#' Reaction net schema (labelled + rates + concentrations)
#' @export
SchReactionNet <- acsets::BasicSchema(
  obs = c("S", "T", "I", "O"),
  homs = list(
    acsets::hom("is", "I", "S"), acsets::hom("it", "I", "T"),
    acsets::hom("os", "O", "S"), acsets::hom("ot", "O", "T")
  ),
  attrtypes = c("Name", "Rate", "Concentration"),
  attrs = list(
    acsets::attr_spec("sname", "S", "Name"),
    acsets::attr_spec("tname", "T", "Name"),
    acsets::attr_spec("rate", "T", "Rate"),
    acsets::attr_spec("concentration", "S", "Concentration")
  )
)

# ACSet type constructors

#' PetriNet ACSet type (unlabelled)
#' @param ... Arguments passed to the ACSet constructor
#' @examples
#' pn <- PetriNet()
#' s1 <- acsets::add_part(pn, "S")
#' t1 <- acsets::add_part(pn, "T")
#' acsets::add_part(pn, "I", is = s1, it = t1)
#' acsets::add_part(pn, "O", os = s1, ot = t1)
#' acsets::nparts(pn, "S") # 1
#' @export
PetriNet <- acsets::acset_type(SchPetriNet, name = "PetriNet",
                                index = c("is", "it", "os", "ot"))

#' LabelledPetriNet ACSet type
#' @param ... Arguments passed to the ACSet constructor
#' @examples
#' pn <- LabelledPetriNet()
#' s1 <- acsets::add_part(pn, "S", sname = "S")
#' s2 <- acsets::add_part(pn, "S", sname = "I")
#' t1 <- acsets::add_part(pn, "T", tname = "inf")
#' acsets::add_part(pn, "I", is = s1, it = t1)
#' acsets::add_part(pn, "O", os = s2, ot = t1)
#' species_names(pn) # c("S", "I")
#' @export
LabelledPetriNet <- acsets::acset_type(SchLabelledPetriNet,
                                        name = "LabelledPetriNet",
                                        index = c("is", "it", "os", "ot"))

#' ReactionNet ACSet type (labelled with rates and concentrations)
#' @param ... Arguments passed to the ACSet constructor
#' @examples
#' rn <- ReactionNet()
#' s1 <- acsets::add_part(rn, "S", sname = "S", concentration = 990)
#' s2 <- acsets::add_part(rn, "S", sname = "I", concentration = 10)
#' t1 <- acsets::add_part(rn, "T", tname = "inf", rate = 0.001)
#' acsets::add_part(rn, "I", is = s1, it = t1)
#' acsets::add_part(rn, "O", os = s2, ot = t1)
#' @export
ReactionNet <- acsets::acset_type(SchReactionNet, name = "ReactionNet",
                                   index = c("is", "it", "os", "ot"))

#' Define a Petri net transition: inputs to outputs
#'
#' @param inputs Character vector of input species names
#' @param outputs Character vector of output species names
#' @returns A list with `inputs` and `outputs`
#' @examples
#' # Infection: S + I -> I + I
#' inf <- c("S", "I") %=>% c("I", "I")
#' inf$inputs  # c("S", "I")
#' inf$outputs # c("I", "I")
#'
#' # Recovery: I -> R
#' rec <- "I" %=>% "R"
#' @export
`%=>%` <- function(inputs, outputs) {
  list(inputs = inputs, outputs = outputs)
}

#' Build a LabelledPetriNet from species names and transition specs
#'
#' @param species Character vector of species names
#' @param ... Named transition specs: `name = inputs %=>% outputs`
#'   where inputs/outputs are character vectors of species names
#' @returns A LabelledPetriNet ACSet
#' @examples
#' # SIR model as a Petri net
#' sir <- labelled_petri_net(
#'   c("S", "I", "R"),
#'   inf = c("S", "I") %=>% c("I", "I"),
#'   rec = "I" %=>% "R"
#' )
#' species_names(sir)    # c("S", "I", "R")
#' transition_names(sir) # c("inf", "rec")
#'
#' # SIS model (no recovery compartment needed)
#' sis <- labelled_petri_net(
#'   c("S", "I"),
#'   inf = c("S", "I") %=>% c("I", "I"),
#'   rec = "I" %=>% "S"
#' )
#' @export
labelled_petri_net <- function(species, ...) {
  transitions <- list(...)
  pn <- LabelledPetriNet()

  # Add species
  sids <- list()
  for (s in species) {
    sids[[s]] <- acsets::add_part(pn, "S", sname = s)
  }

  # Add transitions
  for (tname in names(transitions)) {
    spec <- transitions[[tname]]
    tid <- acsets::add_part(pn, "T", tname = tname)

    # Input arcs
    for (s in spec$inputs) {
      acsets::add_part(pn, "I", is = sids[[s]], it = tid)
    }
    # Output arcs
    for (s in spec$outputs) {
      acsets::add_part(pn, "O", os = sids[[s]], ot = tid)
    }
  }
  pn
}

#' Get species names from a Petri net
#' @param pn A Petri net ACSet (or TypedPetriNet)
#' @returns Character vector of species names
#' @examples
#' sir <- labelled_petri_net(
#'   c("S", "I", "R"),
#'   inf = c("S", "I") %=>% c("I", "I"),
#'   rec = "I" %=>% "R"
#' )
#' species_names(sir)    # c("S", "I", "R")
#' transition_names(sir) # c("inf", "rec")
#' @export
species_names <- function(pn) {
  if (S7::S7_inherits(pn, TypedPetriNet)) pn <- pn@pn
  ids <- acsets::parts(pn, "S")
  if (acsets::has_subpart(pn, "sname")) {
    acsets::subpart(pn, ids, "sname")
  } else {
    paste0("S", ids)
  }
}

#' Get transition names from a Petri net
#' @param pn A Petri net ACSet (or TypedPetriNet)
#' @returns Character vector of transition names
#' @export
transition_names <- function(pn) {
  if (S7::S7_inherits(pn, TypedPetriNet)) pn <- pn@pn
  ids <- acsets::parts(pn, "T")
  if (acsets::has_subpart(pn, "tname")) {
    acsets::subpart(pn, ids, "tname")
  } else {
    paste0("T", ids)
  }
}
