# Composition and stratification operators ------------------------------------

#' Compose two Open Petri nets by matching species names
#'
#' Automatically generates a UWD that connects the two components by
#' identifying species with the same name, then calls \code{oapply}.
#'
#' When applied to non-Open arguments, falls through to \code{base::\%o\%}
#' (outer product).
#'
#' @param a,b Open Petri nets, or objects for \code{base::\%o\%}
#' @returns An Open Petri net (if both are Open), or outer product otherwise
#' @export
`%o%` <- function(a, b) {
  if (S7::S7_inherits(a, Open) && S7::S7_inherits(b, Open)) {
    compose_open(a, b)
  } else {
    base::`%o%`(a, b)
  }
}

#' Typed product (stratification) operator
#'
#' When applied to two TypedPetriNet objects, computes their typed product
#' (stratification). Falls through to \code{base::\%x\%} (Kronecker product)
#' for other types.
#'
#' @param a,b TypedPetriNet objects, or objects for \code{base::\%x\%}
#' @returns A TypedPetriNet (if both are typed), or Kronecker product otherwise
#' @export
`%x%` <- function(a, b) {
  if (S7::S7_inherits(a, TypedPetriNet) && S7::S7_inherits(b, TypedPetriNet)) {
    typed_product(a, b)
  } else {
    base::`%x%`(a, b)
  }
}

#' Compose two Open Petri nets by matching species names
#'
#' Creates a UWD where species with the same name are identified (shared
#' junctions), then composes via \code{oapply}.
#'
#' @param a,b Open Petri nets
#' @returns An Open Petri net
#' @export
compose_open <- function(a, b) {
  snames_a <- species_names(apex(a))
  snames_b <- species_names(apex(b))
  leg_names_a <- snames_a[a@legs[[1]]]
  leg_names_b <- snames_b[b@legs[[1]]]

  all_names <- unique(c(leg_names_a, leg_names_b))

  box_specs <- list(leg_names_a, leg_names_b)
  names(box_specs) <- c("comp1", "comp2")

  w <- do.call(catlab::uwd, c(list(outer = all_names), box_specs))
  oapply(w, list(a, b))
}

#' Compose multiple Open Petri nets
#'
#' Convenience function for composing more than two Open Petri nets
#' by iterative pairwise composition.
#'
#' @param ... Open Petri nets
#' @returns An Open Petri net
#' @export
compose <- function(...) {
  components <- list(...)
  if (length(components) < 2L) {
    cli::cli_abort("compose() requires at least 2 components")
  }
  result <- components[[1]]
  for (i in 2:length(components)) {
    result <- compose_open(result, components[[i]])
  }
  result
}
