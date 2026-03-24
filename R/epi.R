# Epidemiology building blocks ---------------------------------------------
# Pre-built open Petri nets for common epi transitions.

#' Spontaneous transition: A -> B
#' @param input Name of the input species
#' @param output Name of the output species
#' @param tname Optional transition name
#' @returns An Open Petri net
#' @examples
#' # Recovery: I -> R
#' rec <- spontaneous_petri("I", "R", "rec")
#' species_names(apex(rec)) # c("I", "R")
#'
#' # Waning immunity: R -> S
#' wan <- spontaneous_petri("R", "S", "wan")
#' @export
spontaneous_petri <- function(input, output, tname = NULL) {
  if (is.null(tname)) tname <- paste0(input, "_to_", output)
  transitions <- list()
  transitions[[tname]] <- input %=>% output
  pn <- do.call(labelled_petri_net, c(list(c(input, output)), transitions))
  Open(pn, legs = list(seq_len(acsets::nparts(pn, "S"))))
}

#' Exposure/contact transition: A + B -> C + B (B catalyzes A -> C)
#' @param susceptible Name of the susceptible species
#' @param infectious Name of the infectious species
#' @param output Name of the output species
#' @param tname Optional transition name
#' @returns An Open Petri net
#' @examples
#' # Standard infection: S + I -> I + I
#' inf <- exposure_petri("S", "I", "I", "inf")
#' species_names(apex(inf)) # c("S", "I")
#'
#' # Exposure (SEIR): S + I -> E + I
#' exp_petri <- exposure_petri("S", "I", "E", "exp")
#' species_names(apex(exp_petri)) # c("S", "I", "E")
#' @export
exposure_petri <- function(susceptible, infectious, output, tname = NULL) {
  if (is.null(tname)) tname <- paste0(susceptible, "_", infectious, "_to_", output)
  species <- unique(c(susceptible, infectious, output))
  pn <- LabelledPetriNet()
  sids <- list()
  for (s in species) sids[[s]] <- acsets::add_part(pn, "S", sname = s)

  tid <- acsets::add_part(pn, "T", tname = tname)
  acsets::add_part(pn, "I", is = sids[[susceptible]], it = tid)
  acsets::add_part(pn, "I", is = sids[[infectious]], it = tid)
  acsets::add_part(pn, "O", os = sids[[output]], ot = tid)
  acsets::add_part(pn, "O", os = sids[[infectious]], ot = tid)

  Open(pn, legs = list(seq_len(length(species))))
}

#' Pre-built epidemiology dictionary for SIR-family models
#'
#' Returns a function that maps box names to Open Petri nets.
#' @examples
#' d <- epi_dict()
#' names(d) # available building blocks
#' # Each entry is an Open Petri net
#' species_names(apex(d$infection)) # c("S", "I")
#' species_names(apex(d$recovery))  # c("I", "R")
#' @export
epi_dict <- function() {
  list(
    infection = exposure_petri("S", "I", "I", "inf"),
    recovery = spontaneous_petri("I", "R", "rec"),
    vaccination = spontaneous_petri("S", "V", "vac"),
    waning = spontaneous_petri("R", "S", "wan"),
    death = spontaneous_petri("I", "D", "death"),
    exposure = exposure_petri("S", "I", "E", "exp"),
    progression = spontaneous_petri("E", "I", "prog"),
    hospitalization = spontaneous_petri("I", "H", "hosp")
  )
}

#' Compose an epidemiological model from a UWD and dictionary
#'
#' @param w A UWD (wiring diagram)
#' @param dict A named list of Open Petri nets, keyed by box names
#' @param box_names Character vector mapping box indices to dict keys
#' @returns An Open Petri net
#' @examples
#' # Build SIR from epi building blocks
#' d <- epi_dict()
#' w <- catlab::uwd(
#'   outer = c("S", "I", "R"),
#'   infection = c("S", "I"),
#'   recovery = c("I", "R")
#' )
#' sir <- compose_epi(w, d, c("infection", "recovery"))
#' species_names(apex(sir)) # c("S", "I", "R")
#' @export
compose_epi <- function(w, dict, box_names) {
  components <- lapply(box_names, function(nm) {
    if (!(nm %in% names(dict))) {
      cli::cli_abort("Unknown component: '{nm}'")
    }
    dict[[nm]]
  })
  oapply(w, components)
}
