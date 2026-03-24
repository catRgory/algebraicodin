# Typed Petri nets and stratification ----------------------------------------
# Port of AlgebraicPetri.jl's TypedPetri module.
#
# A typed Petri net is a morphism phi: P -> T where T is a type system
# (ontology). Stratification is computed as the categorical product in
# the slice category Petri/T (implemented via pullback).

#' Typed Petri net: a Petri net with a morphism to a type system.
#'
#' @param pn A LabelledPetriNet ACSet
#' @param type_system A LabelledPetriNet ACSet (the ontology)
#' @param species_map Integer vector mapping PN species -> type system species
#' @param transition_map Integer vector mapping PN transitions -> type system transitions
#' @param input_map Integer vector mapping PN input arcs -> type system input arcs
#' @param output_map Integer vector mapping PN output arcs -> type system output arcs
#' @param refl_species Character vector: for each transition, the species name
#'   it loops on (non-empty only for reflexive transitions added by
#'   \code{add_reflexives}).
#' @export
TypedPetriNet <- S7::new_class("TypedPetriNet",
  properties = list(
    pn = acsets::ACSet,
    type_system = acsets::ACSet,
    species_map = S7::class_integer,
    transition_map = S7::class_integer,
    input_map = S7::class_integer,
    output_map = S7::class_integer,
    refl_species = S7::class_character
  ),
  constructor = function(pn, type_system, species_map, transition_map,
                         input_map, output_map, refl_species = character(0)) {
    nt <- length(transition_map)
    if (length(refl_species) == 0L) {
      refl_species <- rep("", nt)
    }
    S7::new_object(S7::S7_object(),
      pn = pn,
      type_system = type_system,
      species_map = as.integer(species_map),
      transition_map = as.integer(transition_map),
      input_map = as.integer(input_map),
      output_map = as.integer(output_map),
      refl_species = refl_species
    )
  }
)

#' Standard infectious disease ontology
#'
#' Returns a type system (LabelledPetriNet) with:
#' \itemize{
#'   \item One species type: \code{Pop}
#'   \item Transition types: \code{infect} (Pop,Pop -> Pop,Pop),
#'     \code{disease} (Pop -> Pop), \code{strata} (Pop -> Pop)
#' }
#'
#' @returns A LabelledPetriNet ACSet
#' @export
infectious_ontology <- function() {
  pn <- LabelledPetriNet()
  pop <- acsets::add_part(pn, "S", sname = "Pop")

  # infect: Pop, Pop -> Pop, Pop (frequency-dependent transmission)
  inf <- acsets::add_part(pn, "T", tname = "infect")
  acsets::add_part(pn, "I", is = pop, it = inf)
  acsets::add_part(pn, "I", is = pop, it = inf)
  acsets::add_part(pn, "O", os = pop, ot = inf)
  acsets::add_part(pn, "O", os = pop, ot = inf)

  # disease: Pop -> Pop (spontaneous transition)
  dis <- acsets::add_part(pn, "T", tname = "disease")
  acsets::add_part(pn, "I", is = pop, it = dis)
  acsets::add_part(pn, "O", os = pop, ot = dis)

  # strata: Pop -> Pop (movement between strata)
  str <- acsets::add_part(pn, "T", tname = "strata")
  acsets::add_part(pn, "I", is = pop, it = str)
  acsets::add_part(pn, "O", os = pop, ot = str)

  pn
}

#' Create a typed Petri net from type assignments
#'
#' @param pn A LabelledPetriNet
#' @param type_system A LabelledPetriNet (the ontology)
#' @param species_types Named character vector mapping species names to
#'   type system species names, e.g. \code{c(S = "Pop", I = "Pop")}
#' @param transition_types Named character vector mapping transition names
#'   to type system transition names, e.g. \code{c(inf = "infect", rec = "disease")}
#' @returns A \code{TypedPetriNet}
#' @export
typed_petri <- function(pn, type_system,
                        species_types = NULL,
                        transition_types = NULL) {
  snames <- species_names(pn)
  tnames <- transition_names(pn)
  ts_snames <- species_names(type_system)
  ts_tnames <- transition_names(type_system)

  # Build species map
  if (is.null(species_types)) {
    species_map <- rep(1L, length(snames))
  } else {
    species_map <- integer(length(snames))
    for (i in seq_along(snames)) {
      type_name <- species_types[snames[i]]
      if (is.na(type_name)) cli::cli_abort("No type specified for species '{snames[i]}'")
      species_map[i] <- match(type_name, ts_snames)
      if (is.na(species_map[i])) cli::cli_abort("Unknown species type '{type_name}'")
    }
  }

  # Build transition map
  if (is.null(transition_types)) {
    cli::cli_abort("{.arg transition_types} must be specified")
  }
  transition_map <- integer(length(tnames))
  for (j in seq_along(tnames)) {
    type_name <- transition_types[tnames[j]]
    if (is.na(type_name)) cli::cli_abort("No type specified for transition '{tnames[j]}'")
    transition_map[j] <- match(type_name, ts_tnames)
    if (is.na(transition_map[j])) cli::cli_abort("Unknown transition type '{type_name}'")
  }

  arc_maps <- compute_arc_maps(pn, type_system, species_map, transition_map)

  TypedPetriNet(
    pn = pn,
    type_system = type_system,
    species_map = species_map,
    transition_map = transition_map,
    input_map = arc_maps$input_map,
    output_map = arc_maps$output_map
  )
}

# Compute which type system arc each PN arc maps to.
# Uses greedy matching by species type, preserving order.
compute_arc_maps <- function(pn, type_system, species_map, transition_map) {
  ni <- acsets::nparts(pn, "I")
  no <- acsets::nparts(pn, "O")
  input_map <- integer(ni)
  output_map <- integer(no)

  nt <- acsets::nparts(pn, "T")
  for (t in seq_len(nt)) {
    tt <- transition_map[t]

    # --- Input arcs ---
    ts_arcs <- acsets::incident(type_system, tt, "it")
    ts_arc_sp <- if (length(ts_arcs) > 0L) {
      vapply(ts_arcs, function(a) acsets::subpart(type_system, a, "is"), integer(1))
    } else integer(0)

    pn_arcs <- acsets::incident(pn, t, "it")
    used <- logical(length(ts_arcs))
    for (k in seq_along(pn_arcs)) {
      s <- acsets::subpart(pn, pn_arcs[k], "is")
      st <- species_map[s]
      for (m in seq_along(ts_arcs)) {
        if (!used[m] && ts_arc_sp[m] == st) {
          input_map[pn_arcs[k]] <- ts_arcs[m]
          used[m] <- TRUE
          break
        }
      }
    }

    # --- Output arcs ---
    ts_arcs <- acsets::incident(type_system, tt, "ot")
    ts_arc_sp <- if (length(ts_arcs) > 0L) {
      vapply(ts_arcs, function(a) acsets::subpart(type_system, a, "os"), integer(1))
    } else integer(0)

    pn_arcs <- acsets::incident(pn, t, "ot")
    used <- logical(length(ts_arcs))
    for (k in seq_along(pn_arcs)) {
      s <- acsets::subpart(pn, pn_arcs[k], "os")
      st <- species_map[s]
      for (m in seq_along(ts_arcs)) {
        if (!used[m] && ts_arc_sp[m] == st) {
          output_map[pn_arcs[k]] <- ts_arcs[m]
          used[m] <- TRUE
          break
        }
      }
    }
  }

  list(input_map = input_map, output_map = output_map)
}

#' Add reflexive (self-loop) transitions for stratification
#'
#' Reflexive transitions are identity transitions that map a species to
#' itself. They are required so that \code{typed_product} can pair
#' transitions from one model with the corresponding species of the other.
#'
#' @param typed_pn A \code{TypedPetriNet}
#' @param reflexives Named list mapping species names to character vectors
#'   of transition type names, e.g.
#'   \code{list(S = "strata", I = "strata", R = "strata")}
#' @returns A new \code{TypedPetriNet} with reflexive transitions added
#' @export
add_reflexives <- function(typed_pn, reflexives) {
  pn <- acsets::copy_acset(typed_pn@pn)
  ts <- typed_pn@type_system
  snames <- species_names(typed_pn@pn)
  ts_tnames <- transition_names(ts)

  new_tmap <- typed_pn@transition_map
  new_imap <- typed_pn@input_map
  new_omap <- typed_pn@output_map
  new_refl <- typed_pn@refl_species

  for (sname in names(reflexives)) {
    s <- match(sname, snames)
    if (is.na(s)) cli::cli_abort("Species '{sname}' not found in Petri net")

    ttype_names <- reflexives[[sname]]
    for (ttype_name in ttype_names) {
      tt <- match(ttype_name, ts_tnames)
      if (is.na(tt)) cli::cli_abort("Transition type '{ttype_name}' not found in type system")

      refl_name <- paste0("refl_", sname, "_", ttype_name)
      new_t <- acsets::add_part(pn, "T", tname = refl_name)
      new_tmap <- c(new_tmap, tt)
      new_refl <- c(new_refl, sname)

      # Input arcs matching type system
      for (ta in acsets::incident(ts, tt, "it")) {
        acsets::add_part(pn, "I", is = s, it = new_t)
        new_imap <- c(new_imap, ta)
      }

      # Output arcs matching type system
      for (ta in acsets::incident(ts, tt, "ot")) {
        acsets::add_part(pn, "O", os = s, ot = new_t)
        new_omap <- c(new_omap, ta)
      }
    }
  }

  TypedPetriNet(
    pn = pn,
    type_system = ts,
    species_map = typed_pn@species_map,
    transition_map = new_tmap,
    input_map = new_imap,
    output_map = new_omap,
    refl_species = new_refl
  )
}

#' Compute the typed product (stratification) of typed Petri nets
#'
#' This is the core stratification operation. Given two typed Petri nets
#' sharing the same type system, the typed product creates a new Petri net
#' whose species are pairs of species with matching types, and whose
#' transitions are pairs of transitions with matching types.
#'
#' @param tp1,tp2 \code{TypedPetriNet} objects with the same type system
#' @returns A \code{TypedPetriNet} (the stratified model)
#' @export
typed_product <- function(tp1, tp2) {
  pn1 <- tp1@pn
  pn2 <- tp2@pn
  ts <- tp1@type_system

  snames1 <- species_names(pn1)
  snames2 <- species_names(pn2)
  tnames1 <- transition_names(pn1)
  tnames2 <- transition_names(pn2)
  refl1 <- tp1@refl_species
  refl2 <- tp2@refl_species

  ns1 <- length(snames1)
  ns2 <- length(snames2)

  result <- LabelledPetriNet()

  # --- Product species ---
  sp_pair <- matrix(0L, ns1, ns2)
  result_smap <- integer(0)

  for (s1 in seq_len(ns1)) {
    for (s2 in seq_len(ns2)) {
      if (tp1@species_map[s1] == tp2@species_map[s2]) {
        name <- paste0(snames1[s1], "_", snames2[s2])
        sid <- acsets::add_part(result, "S", sname = name)
        sp_pair[s1, s2] <- sid
        result_smap <- c(result_smap, tp1@species_map[s1])
      }
    }
  }

  # --- Product transitions ---
  result_tmap <- integer(0)
  result_imap <- integer(0)
  result_omap <- integer(0)
  result_refl <- character(0)

  nt1 <- length(tnames1)
  nt2 <- length(tnames2)
  ts_tnames <- transition_names(ts)

  for (t1 in seq_len(nt1)) {
    for (t2 in seq_len(nt2)) {
      if (tp1@transition_map[t1] != tp2@transition_map[t2]) next

      # Smart naming based on reflexive status
      r1 <- refl1[t1]
      r2 <- refl2[t2]
      if (nchar(r1) > 0 && nchar(r2) == 0) {
        tname <- paste0(tnames2[t2], "_", r1)
      } else if (nchar(r1) == 0 && nchar(r2) > 0) {
        tname <- paste0(tnames1[t1], "_", r2)
      } else if (nchar(r1) > 0 && nchar(r2) > 0) {
        tname <- paste0(r1, "_", r2, "_", ts_tnames[tp1@transition_map[t1]])
      } else {
        tname <- paste0(tnames1[t1], "_", tnames2[t2])
      }

      tid <- acsets::add_part(result, "T", tname = tname)
      result_tmap <- c(result_tmap, tp1@transition_map[t1])

      if (nchar(r1) > 0 && nchar(r2) > 0) {
        result_refl <- c(result_refl, paste0(r1, "_", r2))
      } else {
        result_refl <- c(result_refl, "")
      }

      # Input arcs: pair by type system arc
      arcs1 <- acsets::incident(pn1, t1, "it")
      arcs2 <- acsets::incident(pn2, t2, "it")
      for (a1 in arcs1) {
        for (a2 in arcs2) {
          if (tp1@input_map[a1] == tp2@input_map[a2]) {
            s1 <- acsets::subpart(pn1, a1, "is")
            s2 <- acsets::subpart(pn2, a2, "is")
            ps <- sp_pair[s1, s2]
            if (ps > 0L) {
              acsets::add_part(result, "I", is = ps, it = tid)
              result_imap <- c(result_imap, tp1@input_map[a1])
            }
          }
        }
      }

      # Output arcs: pair by type system arc
      arcs1 <- acsets::incident(pn1, t1, "ot")
      arcs2 <- acsets::incident(pn2, t2, "ot")
      for (a1 in arcs1) {
        for (a2 in arcs2) {
          if (tp1@output_map[a1] == tp2@output_map[a2]) {
            s1 <- acsets::subpart(pn1, a1, "os")
            s2 <- acsets::subpart(pn2, a2, "os")
            ps <- sp_pair[s1, s2]
            if (ps > 0L) {
              acsets::add_part(result, "O", os = ps, ot = tid)
              result_omap <- c(result_omap, tp1@output_map[a1])
            }
          }
        }
      }
    }
  }

  TypedPetriNet(
    pn = result,
    type_system = ts,
    species_map = result_smap,
    transition_map = result_tmap,
    input_map = result_imap,
    output_map = result_omap,
    refl_species = result_refl
  )
}

#' Extract the underlying Petri net from a TypedPetriNet
#' @param typed_pn A TypedPetriNet
#' @returns A Petri net ACSet
#' @export
flatten <- function(typed_pn) typed_pn@pn

#' User-friendly stratification wrapper
#'
#' Computes the typed product of two TypedPetriNets.
#' Equivalent to \code{typed_product(base, strata)}.
#'
#' @param base A \code{TypedPetriNet} (the base model, e.g. SIR)
#' @param strata A \code{TypedPetriNet} (the stratification model, e.g. age groups)
#' @returns A \code{TypedPetriNet} (the stratified model)
#' @export
stratify <- function(base, strata) typed_product(base, strata)
