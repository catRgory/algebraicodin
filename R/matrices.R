# Transition matrices and stoichiometry ------------------------------------

#' Extract input/output stoichiometry matrices
#'
#' @param pn A Petri net ACSet
#' @returns A list with `input` and `output` matrices (species × transitions)
#' @export
transition_matrices <- function(pn) {
  if (S7::S7_inherits(pn, TypedPetriNet)) pn <- pn@pn
  ns <- acsets::nparts(pn, "S")
  nt <- acsets::nparts(pn, "T")

  input_mat <- matrix(0L, nrow = ns, ncol = nt)
  output_mat <- matrix(0L, nrow = ns, ncol = nt)

  # Fill input matrix
  ni <- acsets::nparts(pn, "I")
  for (i in seq_len(ni)) {
    s <- acsets::subpart(pn, i, "is")
    t <- acsets::subpart(pn, i, "it")
    if (!is.na(s) && !is.na(t)) input_mat[s, t] <- input_mat[s, t] + 1L
  }

  # Fill output matrix
  no <- acsets::nparts(pn, "O")
  for (i in seq_len(no)) {
    s <- acsets::subpart(pn, i, "os")
    t <- acsets::subpart(pn, i, "ot")
    if (!is.na(s) && !is.na(t)) output_mat[s, t] <- output_mat[s, t] + 1L
  }

  snames <- species_names(pn)
  tnames <- transition_names(pn)
  rownames(input_mat) <- snames
  rownames(output_mat) <- snames
  colnames(input_mat) <- tnames
  colnames(output_mat) <- tnames

  list(input = input_mat, output = output_mat, stoichiometry = output_mat - input_mat)
}

#' Generate a pure-R mass-action vector field function
#'
#' @param pn A Petri net ACSet
#' @param rate_names Optional named vector mapping transition names to rate parameter names
#' @returns A function(t, state, parms) suitable for deSolve::ode
#' @export
vectorfield <- function(pn, rate_names = NULL) {
  if (S7::S7_inherits(pn, TypedPetriNet)) pn <- pn@pn
  tm <- transition_matrices(pn)
  snames <- species_names(pn)
  tnames <- transition_names(pn)
  ns <- length(snames)
  nt <- length(tnames)

  if (is.null(rate_names)) {
    rate_names <- tnames
    names(rate_names) <- tnames
  }

  input_mat <- tm$input
  stoich <- tm$stoichiometry

  function(t, state, parms) {
    # Mass-action kinetics: rate_j * prod(state_i^input_ij)
    rates <- numeric(nt)
    for (j in seq_len(nt)) {
      rate <- parms[[rate_names[j]]]
      for (i in seq_len(ns)) {
        if (input_mat[i, j] > 0L) {
          rate <- rate * state[i]^input_mat[i, j]
        }
      }
      rates[j] <- rate
    }
    dstate <- as.numeric(stoich %*% rates)
    names(dstate) <- snames
    list(dstate)
  }
}
