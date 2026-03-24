# Dynamical systems ---------------------------------------------------------
# Port of AlgebraicDynamics.jl's ContinuousMachine, DiscreteMachine,
# and DelayMachine types.
#
# These represent open dynamical systems that can be composed via
# wiring diagrams. Each system has:
#   - states (internal variables)
#   - dynamics function (how state evolves)
#   - interface (how it connects to the outside world)
#
# Three model types (matching odin2 capabilities):
#   ContinuousSystem: du/dt = f(u, x, p, t)          -> odin2 deriv()
#   DiscreteSystem:   u_{n+1} = f(u, x, p, t)        -> odin2 update()
#   DelaySystem:      du/dt = f(u, x, h, p, t)       -> odin2 deriv() + delay()

# --- ResourceSharer (undirected, UWD-composable) --------------------------

#' A resource sharer: open dynamical system with symmetric ports.
#'
#' Corresponds to AlgebraicDynamics.jl's ResourceSharer. The `portmap`
#' specifies which internal state variables are exposed at each port.
#'
#' @param system_type One of "continuous", "discrete", "delay"
#' @param nstates Number of internal state variables
#' @param dynamics Function implementing the system evolution
#' @param portmap Integer vector mapping ports → state indices
#' @export
ResourceSharer <- S7::new_class("ResourceSharer",
  properties = list(
    system_type = S7::class_character,
    nstates = S7::class_integer,
    dynamics = S7::class_function,
    portmap = S7::class_integer,
    state_names = S7::class_character,
    dynamics_expr = S7::class_list,
    params = S7::class_character
  ),
  constructor = function(system_type = "continuous", nstates = 0L,
                         dynamics = function(u, p, t) numeric(0),
                         portmap = seq_len(nstates),
                         state_names = character(0),
                         dynamics_expr = list(),
                         params = character(0)) {
    S7::new_object(S7::S7_object(),
      system_type = system_type,
      nstates = as.integer(nstates),
      dynamics = dynamics,
      portmap = as.integer(portmap),
      state_names = state_names,
      dynamics_expr = dynamics_expr,
      params = params
    )
  },
  validator = function(self) {
    if (!(self@system_type %in% c("continuous", "discrete", "delay"))) {
      return("system_type must be 'continuous', 'discrete', or 'delay'")
    }
    if (any(self@portmap < 1L | self@portmap > self@nstates)) {
      return("portmap indices must be in 1:nstates")
    }
    NULL
  }
)

#' Number of ports
#' @export
nports <- function(rs) length(rs@portmap)

#' Evaluate dynamics of a ResourceSharer
#' @export
eval_dynamics <- S7::new_generic("eval_dynamics", "rs")

S7::method(eval_dynamics, ResourceSharer) <- function(rs, u, p = NULL, t = 0, h = NULL) {
  if (rs@system_type == "delay" && !is.null(h)) {
    rs@dynamics(u, h, p, t)
  } else {
    rs@dynamics(u, p, t)
  }
}

#' Get exposed state values at ports
#' @export
exposed_states <- function(rs, u) u[rs@portmap]

# --- Machine (directed, DWD-composable) -----------------------------------

#' A machine: open dynamical system with directed inputs and outputs.
#'
#' Corresponds to AlgebraicDynamics.jl's Machine. Has explicit input ports
#' (receiving signals) and output ports (emitting readouts).
#'
#' @param system_type One of "continuous", "discrete", "delay"
#' @param ninputs Number of input ports
#' @param nstates Number of internal state variables
#' @param noutputs Number of output ports
#' @param dynamics Function: (u, x, p, t) -> du  [or (u, x, h, p, t) for delay]
#' @param readout Function: (u, p, t) -> y  [or (u, h, p, t) for delay]
#' @export
Machine <- S7::new_class("Machine",
  properties = list(
    system_type = S7::class_character,
    ninputs = S7::class_integer,
    nstates = S7::class_integer,
    noutputs = S7::class_integer,
    dynamics = S7::class_function,
    readout = S7::class_function,
    state_names = S7::class_character,
    dynamics_expr = S7::class_list,
    readout_expr = S7::class_list,
    params = S7::class_character
  ),
  constructor = function(system_type = "continuous",
                         ninputs = 0L, nstates = 0L, noutputs = 0L,
                         dynamics = function(u, x, p, t) numeric(0),
                         readout = function(u, p, t) u,
                         state_names = character(0),
                         dynamics_expr = list(),
                         readout_expr = list(),
                         params = character(0)) {
    S7::new_object(S7::S7_object(),
      system_type = system_type,
      ninputs = as.integer(ninputs),
      nstates = as.integer(nstates),
      noutputs = as.integer(noutputs),
      dynamics = dynamics,
      readout = readout,
      state_names = state_names,
      dynamics_expr = dynamics_expr,
      readout_expr = readout_expr,
      params = params
    )
  }
)

S7::method(eval_dynamics, Machine) <- function(rs, u, x = numeric(0),
                                                p = NULL, t = 0, h = NULL) {
  if (rs@system_type == "delay" && !is.null(h)) {
    rs@dynamics(u, x, h, p, t)
  } else {
    rs@dynamics(u, x, p, t)
  }
}

#' Evaluate readout of a Machine
#' @export
eval_readout <- function(m, u, p = NULL, t = 0, h = NULL) {
  if (m@system_type == "delay" && !is.null(h)) {
    m@readout(u, h, p, t)
  } else {
    m@readout(u, p, t)
  }
}

# --- Euler approximation: continuous → discrete ---------------------------

#' Convert a continuous system to discrete via Euler's method
#'
#' @param rs A ResourceSharer or Machine with system_type = "continuous"
#' @param h Step size (numeric). If NULL, step size is appended to parameters.
#' @returns A discrete system of the same class
#' @export
euler_approx <- function(rs, h = NULL) {
  if (!S7::S7_inherits(rs, ResourceSharer) && !S7::S7_inherits(rs, Machine)) {
    cli::cli_abort("euler_approx requires a ResourceSharer or Machine")
  }
  if (rs@system_type != "continuous") {
    cli::cli_abort("euler_approx only applies to continuous systems")
  }

  if (S7::S7_inherits(rs, ResourceSharer)) {
    orig_dyn <- rs@dynamics
    if (!is.null(h)) {
      h <- as.numeric(h)
      new_dyn <- function(u, p, t) u + h * orig_dyn(u, p, t)
    } else {
      new_dyn <- function(u, p, t) {
        dt <- p[["dt"]]
        u + dt * orig_dyn(u, p, t)
      }
    }
    ResourceSharer(
      system_type = "discrete",
      nstates = rs@nstates,
      dynamics = new_dyn,
      portmap = rs@portmap,
      state_names = rs@state_names
    )
  } else {
    orig_dyn <- rs@dynamics
    if (!is.null(h)) {
      h <- as.numeric(h)
      new_dyn <- function(u, x, p, t) u + h * orig_dyn(u, x, p, t)
    } else {
      new_dyn <- function(u, x, p, t) {
        dt <- p[["dt"]]
        u + dt * orig_dyn(u, x, p, t)
      }
    }
    Machine(
      system_type = "discrete",
      ninputs = rs@ninputs,
      nstates = rs@nstates,
      noutputs = rs@noutputs,
      dynamics = new_dyn,
      readout = rs@readout,
      state_names = rs@state_names
    )
  }
}

# --- Petri net → dynamical system converters ------------------------------

#' Convert a Petri net to a ContinuousResourceSharer (ODE)
#' @export
petri_to_continuous <- function(pn) {
  tm <- transition_matrices(pn)
  sn <- species_names(pn)
  tn <- transition_names(pn)
  ns <- length(sn)
  nt <- length(tn)

  dynamics <- function(u, p, t) {
    rates <- numeric(nt)
    for (j in seq_len(nt)) {
      r <- p[[tn[j]]]
      for (i in seq_len(ns)) {
        if (tm$input[i, j] > 0L) r <- r * u[i]^tm$input[i, j]
      }
      rates[j] <- r
    }
    as.numeric(tm$stoichiometry %*% rates)
  }

  ResourceSharer(
    system_type = "continuous",
    nstates = ns,
    dynamics = dynamics,
    portmap = seq_len(ns),
    state_names = sn
  )
}

#' Convert a Petri net to a DiscreteResourceSharer (stochastic)
#' @export
petri_to_discrete <- function(pn) {
  tm <- transition_matrices(pn)
  sn <- species_names(pn)
  tn <- transition_names(pn)
  ns <- length(sn)
  nt <- length(tn)

  dynamics <- function(u, p, t) {
    dt <- p[["dt"]]
    draws <- numeric(nt)
    for (j in seq_len(nt)) {
      input_sp <- which(tm$input[, j] > 0L)
      if (length(input_sp) == 0L) next
      substrate <- input_sp[1]
      # Per-capita hazard (exclude substrate from rate)
      hazard <- p[[tn[j]]]
      for (idx in input_sp[-1]) {
        hazard <- hazard * u[idx]^tm$input[idx, j]
      }
      prob <- 1 - exp(-hazard * dt)
      draws[j] <- stats::rbinom(1, size = as.integer(u[substrate]), prob = min(prob, 1))
    }
    u + as.numeric(tm$stoichiometry %*% draws)
  }

  ResourceSharer(
    system_type = "discrete",
    nstates = ns,
    dynamics = dynamics,
    portmap = seq_len(ns),
    state_names = sn
  )
}

#' Convert a Petri net to a DelayResourceSharer (DDE)
#'
#' Creates a DDE system where specified transitions use delayed state values.
#'
#' @param pn A Petri net ACSet
#' @param delays Named list mapping transition names to delay specs.
#'   Each spec is a list with `tau` (delay time) and optionally `species`
#'   (which input species are delayed; default: all inputs).
#' @export
petri_to_delay <- function(pn, delays) {
  tm <- transition_matrices(pn)
  sn <- species_names(pn)
  tn <- transition_names(pn)
  ns <- length(sn)
  nt <- length(tn)

  dynamics <- function(u, h, p, t) {
    rates <- numeric(nt)
    for (j in seq_len(nt)) {
      r <- p[[tn[j]]]
      dspec <- delays[[tn[j]]]
      for (i in seq_len(ns)) {
        if (tm$input[i, j] > 0L) {
          if (!is.null(dspec) && (is.null(dspec$species) || sn[i] %in% dspec$species)) {
            # Use delayed value
            u_delayed <- h(p, t - dspec$tau)
            r <- r * u_delayed[i]^tm$input[i, j]
          } else {
            r <- r * u[i]^tm$input[i, j]
          }
        }
      }
      rates[j] <- r
    }
    as.numeric(tm$stoichiometry %*% rates)
  }

  ResourceSharer(
    system_type = "delay",
    nstates = ns,
    dynamics = dynamics,
    portmap = seq_len(ns),
    state_names = sn
  )
}

# --- Symbolic expression helpers ------------------------------------------
# Used by rs_to_odin / machine_to_odin to generate compilable odin2 code.
# Expressions are stored as named lists of R language objects (e.g. quote(...)).

#' Rename variables in an expression
#' @param expr A language object (quoted R expression)
#' @param mapping Named character vector: old_name -> new_name
#' @returns The expression with names substituted
#' @keywords internal
rename_expr_vars <- function(expr, mapping) {
  if (is.symbol(expr)) {
    nm <- as.character(expr)
    if (nm %in% names(mapping)) return(as.symbol(mapping[[nm]]))
    return(expr)
  }
  if (is.numeric(expr) || is.character(expr) || is.logical(expr)) return(expr)
  if (is.call(expr)) {
    expr[-1] <- lapply(expr[-1], rename_expr_vars, mapping = mapping)
    return(expr)
  }
  expr
}

#' Sum a list of expressions into a single expression
#' @param exprs List of language objects to sum
#' @returns A single expression representing the sum
#' @keywords internal
sum_exprs <- function(exprs) {
  exprs <- Filter(Negate(is.null), exprs)
  if (length(exprs) == 0L) return(quote(0))
  if (length(exprs) == 1L) return(exprs[[1]])
  Reduce(function(a, b) call("+", a, b), exprs)
}

#' Build a dynamics closure from symbolic expressions
#' @param state_names Character vector of state variable names
#' @param dynamics_expr Named list of quoted expressions
#' @param param_names Character vector of parameter names
#' @returns A function(u, p, t) -> du suitable for ResourceSharer
#' @keywords internal
build_rs_closure <- function(state_names, dynamics_expr, param_names) {
  force(state_names)
  force(dynamics_expr)
  force(param_names)
  function(u, p, t) {
    env <- new.env(parent = baseenv())
    for (i in seq_along(state_names)) env[[state_names[i]]] <- u[i]
    if (is.list(p)) {
      for (nm in names(p)) env[[nm]] <- p[[nm]]
    }
    env[["t"]] <- t
    du <- numeric(length(state_names))
    for (i in seq_along(state_names)) {
      du[i] <- eval(dynamics_expr[[state_names[i]]], envir = env)
    }
    du
  }
}

#' Build a dynamics closure from symbolic expressions for a Machine
#' @param state_names Character vector of state variable names
#' @param ninputs Number of input ports
#' @param dynamics_expr Named list of quoted expressions
#' @returns A function(u, x, p, t) -> du
#' @keywords internal
build_machine_closure <- function(state_names, ninputs, dynamics_expr) {
  force(state_names)
  force(ninputs)
  force(dynamics_expr)
  function(u, x, p, t) {
    env <- new.env(parent = baseenv())
    for (i in seq_along(state_names)) env[[state_names[i]]] <- u[i]
    for (i in seq_len(ninputs)) env[[paste0("input_", i)]] <- x[i]
    if (is.list(p)) {
      for (nm in names(p)) env[[nm]] <- p[[nm]]
    }
    env[["t"]] <- t
    du <- numeric(length(state_names))
    for (i in seq_along(state_names)) {
      du[i] <- eval(dynamics_expr[[state_names[i]]], envir = env)
    }
    du
  }
}

#' Build a readout closure from symbolic expressions for a Machine
#' @param state_names Character vector of state variable names
#' @param readout_expr List of quoted expressions (one per output port)
#' @returns A function(u, p, t) -> y
#' @keywords internal
build_readout_closure <- function(state_names, readout_expr) {
  force(state_names)
  force(readout_expr)
  function(u, p, t) {
    env <- new.env(parent = baseenv())
    for (i in seq_along(state_names)) env[[state_names[i]]] <- u[i]
    if (is.list(p)) {
      for (nm in names(p)) env[[nm]] <- p[[nm]]
    }
    env[["t"]] <- t
    vapply(readout_expr, eval, numeric(1), envir = env)
  }
}

# --- oapply for ResourceSharers (UWD composition) -------------------------
# Mirrors AlgebraicDynamics.jl's oapply for resource sharers.
#
# When resource sharers are composed via a UWD:
#  1. Each box's internal states are collected into one big state vector
#     (coproduct of state spaces).
#  2. States connected to the same junction are identified (pushout).
#  3. The composite dynamics sums derivatives at shared junctions.
#
# The result is a new ResourceSharer whose portmap corresponds to
# outer ports of the UWD.

#' Compose ResourceSharers via an undirected wiring diagram
#'
#' Implements the operad algebra for undirected composition of dynamical
#' systems. Each box of the UWD is filled by a resource sharer
#' from \code{sharers}. States at the same junction are identified and
#' derivatives (or updates) are summed.
#'
#' @param d A UWD ACSet (from catlab)
#' @param sharers A list of ResourceSharer objects, one per box
#' @returns A composite ResourceSharer
#' @export
oapply_dynam <- function(d, sharers) {
  nboxes <- acsets::nparts(d, "Box")
  njunctions <- acsets::nparts(d, "Junction")

  if (length(sharers) != nboxes)
    cli::cli_abort("Need {nboxes} sharers but got {length(sharers)}")

  # Validate port counts match
  for (b in seq_len(nboxes)) {
    box_ports <- which(acsets::subpart(d, NULL, "box") == b)
    expected_nports <- length(box_ports)
    if (nports(sharers[[b]]) != expected_nports) {
      cli::cli_abort(
        "Box {b} has {expected_nports} ports in UWD but sharer has {nports(sharers[[b]])} ports"
      )
    }
  }

  # Determine system type (all must match)
  sys_type <- sharers[[1]]@system_type
  for (s in sharers) {
    if (s@system_type != sys_type)
      cli::cli_abort("All sharers must have the same system_type")
  }

  # Build coproduct: state offsets for each box
  state_counts <- vapply(sharers, function(s) s@nstates, integer(1))
  offsets <- c(0L, cumsum(state_counts))
  total_states <- offsets[nboxes + 1L]
  box_states <- function(b) (offsets[b] + 1L):offsets[b + 1L]

  # Build junction -> states mapping (pushout)
  # For each box port, the connected junction gets the corresponding state
  port_boxes <- acsets::subpart(d, NULL, "box")
  port_junctions <- acsets::subpart(d, NULL, "junction")

  # junction_states[[j]] = list of (global_state_index) values
  junction_states <- vector("list", njunctions)
  for (j in seq_len(njunctions)) junction_states[[j]] <- integer(0)

  for (p in seq_len(acsets::nparts(d, "Port"))) {
    b <- port_boxes[p]
    j <- port_junctions[p]
    # Which port number is this within box b?
    box_port_indices <- which(port_boxes == b)
    local_port_idx <- match(p, box_port_indices)
    # The portmap of sharer b maps local port -> local state
    local_state <- sharers[[b]]@portmap[local_port_idx]
    global_state <- offsets[b] + local_state
    junction_states[[j]] <- c(junction_states[[j]], global_state)
  }

  # Build the pushout: identify states at the same junction.
  # Use union-find to merge states at shared junctions.
  parent <- seq_len(total_states)
  find_root <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  union_states <- function(a, b) {
    ra <- find_root(a)
    rb <- find_root(b)
    if (ra != rb) parent[ra] <<- rb
  }

  for (j in seq_len(njunctions)) {
    sts <- unique(junction_states[[j]])
    if (length(sts) > 1L) {
      for (k in 2:length(sts)) {
        union_states(sts[1], sts[k])
      }
    }
  }

  # Compress to 1..n_composite
  roots <- vapply(seq_len(total_states), find_root, integer(1))
  unique_roots <- sort(unique(roots))
  n_composite <- length(unique_roots)
  root_to_composite <- integer(max(unique_roots))
  root_to_composite[unique_roots] <- seq_len(n_composite)
  composite_map <- root_to_composite[roots]  # old global -> new composite

  # Inverse map: for each composite state, which old global states map to it
  inverse_map <- vector("list", n_composite)
  for (i in seq_len(n_composite)) inverse_map[[i]] <- integer(0)
  for (g in seq_len(total_states)) {
    c_idx <- composite_map[g]
    inverse_map[[c_idx]] <- c(inverse_map[[c_idx]], g)
  }

  # Build composite portmap (outer ports -> composite states)
  outer_junctions <- acsets::subpart(d, NULL, "outer_junction")
  n_outer <- length(outer_junctions)
  composite_portmap <- integer(n_outer)
  for (op in seq_len(n_outer)) {
    j <- outer_junctions[op]
    sts <- junction_states[[j]]
    if (length(sts) == 0L) {
      cli::cli_abort("Outer port {op} connects to junction {j} with no states")
    }
    composite_portmap[op] <- composite_map[sts[1]]
  }

  # Build composite state names
  all_snames <- character(total_states)
  for (b in seq_len(nboxes)) {
    sn <- sharers[[b]]@state_names
    idx <- box_states(b)
    if (length(sn) == length(idx)) {
      all_snames[idx] <- sn
    } else {
      all_snames[idx] <- paste0("x", idx)
    }
  }
  # For each composite state, take the name of the first contributing state
  composite_names <- all_snames[vapply(inverse_map, `[`, integer(1), 1L)]

  # Build composite dynamics
  if (sys_type == "continuous") {
    composite_dynamics <- function(u, p, t) {
      # Expand composite state to full state vector
      u_full <- u[composite_map]
      du_full <- numeric(total_states)
      for (b in seq_len(nboxes)) {
        idx <- box_states(b)
        du_box <- sharers[[b]]@dynamics(u_full[idx], p, t)
        du_full[idx] <- du_full[idx] + du_box
      }
      # Contract: sum derivatives at identified states
      du <- numeric(n_composite)
      for (c_idx in seq_len(n_composite)) {
        du[c_idx] <- sum(du_full[inverse_map[[c_idx]]])
      }
      du
    }
  } else if (sys_type == "discrete") {
    composite_dynamics <- function(u, p, t) {
      u_full <- u[composite_map]
      delta_full <- numeric(total_states)
      for (b in seq_len(nboxes)) {
        idx <- box_states(b)
        u1_box <- sharers[[b]]@dynamics(u_full[idx], p, t)
        delta_full[idx] <- u1_box - u_full[idx]
      }
      # Sum deltas at identified states
      delta <- numeric(n_composite)
      for (c_idx in seq_len(n_composite)) {
        delta[c_idx] <- sum(delta_full[inverse_map[[c_idx]]])
      }
      u + delta
    }
  } else {
    # delay
    composite_dynamics <- function(u, h, p, t) {
      u_full <- u[composite_map]
      du_full <- numeric(total_states)
      for (b in seq_len(nboxes)) {
        idx <- box_states(b)
        h_box <- function(p2, t2) {
          h_val <- h(p2, t2)
          h_val[composite_map[idx]]
        }
        du_box <- sharers[[b]]@dynamics(u_full[idx], h_box, p, t)
        du_full[idx] <- du_full[idx] + du_box
      }
      du <- numeric(n_composite)
      for (c_idx in seq_len(n_composite)) {
        du[c_idx] <- sum(du_full[inverse_map[[c_idx]]])
      }
      du
    }
  }

  ResourceSharer(
    system_type = sys_type,
    nstates = n_composite,
    dynamics = composite_dynamics,
    portmap = composite_portmap,
    state_names = composite_names,
    dynamics_expr = merge_rs_expressions(sharers, all_snames, composite_names,
                                          composite_map, inverse_map, n_composite,
                                          nboxes, box_states),
    params = unique(unlist(lapply(sharers, function(s) s@params)))
  )
}

#' Merge symbolic expressions during oapply_dynam composition
#' @keywords internal
merge_rs_expressions <- function(sharers, all_snames, composite_names,
                                  composite_map, inverse_map, n_composite,
                                  nboxes, box_states) {
  # Check if ALL sharers have expressions
  has_exprs <- vapply(sharers, function(s) length(s@dynamics_expr) > 0L, logical(1))
  if (!all(has_exprs)) return(list())

  # Build per-global-state expressions with renaming
  total_states <- length(all_snames)
  global_exprs <- vector("list", total_states)

  for (b in seq_len(nboxes)) {
    idx <- box_states(b)
    sn <- sharers[[b]]@state_names
    dexpr <- sharers[[b]]@dynamics_expr

    # Build renaming: local state names -> composite state names
    mapping <- character(0)
    for (k in seq_along(idx)) {
      local_nm <- sn[k]
      composite_idx <- composite_map[idx[k]]
      comp_nm <- composite_names[composite_idx]
      if (local_nm != comp_nm) mapping[local_nm] <- comp_nm
    }

    for (k in seq_along(idx)) {
      local_nm <- sn[k]
      g <- idx[k]
      if (local_nm %in% names(dexpr)) {
        expr <- dexpr[[local_nm]]
        if (length(mapping) > 0L) expr <- rename_expr_vars(expr, mapping)
        global_exprs[[g]] <- expr
      }
    }
  }

  # Sum expressions for identified states
  result <- list()
  for (c_idx in seq_len(n_composite)) {
    contrib <- lapply(inverse_map[[c_idx]], function(g) global_exprs[[g]])
    result[[composite_names[c_idx]]] <- sum_exprs(contrib)
  }
  result
}

# --- oapply for Machines (DWD composition) --------------------------------

#' Compose Machines via a directed wiring diagram
#'
#' Each box of the directed wiring diagram is filled by a Machine.
#' Outputs of upstream machines are wired to inputs of downstream ones.
#'
#' @param d A directed wiring diagram (from \code{dwd()})
#' @param machines A list of Machine objects, one per box
#' @returns A composite Machine
#' @export
oapply_machine <- function(d, machines) {
  nboxes <- length(d$boxes)
  if (length(machines) != nboxes)
    cli::cli_abort("Need {nboxes} machines but got {length(machines)}")

  sys_type <- machines[[1]]@system_type

  # State offsets
  state_counts <- vapply(machines, function(m) m@nstates, integer(1))
  offsets <- c(0L, cumsum(state_counts))
  total_states <- offsets[nboxes + 1L]
  box_states_fn <- function(b) {
    if (state_counts[b] == 0L) return(integer(0))
    (offsets[b] + 1L):offsets[b + 1L]
  }

  n_total_inputs <- d$ninputs
  n_total_outputs <- d$noutputs

  # Helper: get readout for all boxes
  get_readouts <- function(states_list, p, t) {
    lapply(seq_len(nboxes), function(b) {
      eval_readout(machines[[b]], states_list[[b]], p, t)
    })
  }

  # Helper: get inputs for a box from external inputs and readouts
  get_inputs <- function(b, xs, readouts) {
    nin <- machines[[b]]@ninputs
    inputs <- numeric(nin)
    for (w in d$wires) {
      if (w$target_box == b) {
        val <- if (w$source_box == 0L) {
          xs[w$source_port]
        } else {
          readouts[[w$source_box]][w$source_port]
        }
        inputs[w$target_port] <- inputs[w$target_port] + val
      }
    }
    inputs
  }

  # Composite dynamics
  if (sys_type %in% c("continuous", "discrete")) {
    comp_dynamics <- function(u, x, p, t) {
      states_list <- lapply(seq_len(nboxes), function(b) u[box_states_fn(b)])
      readouts <- get_readouts(states_list, p, t)
      du <- numeric(total_states)
      for (b in seq_len(nboxes)) {
        inputs <- get_inputs(b, x, readouts)
        idx <- box_states_fn(b)
        if (length(idx) > 0L)
          du[idx] <- machines[[b]]@dynamics(states_list[[b]], inputs, p, t)
      }
      du
    }
  } else {
    comp_dynamics <- function(u, x, h, p, t) {
      states_list <- lapply(seq_len(nboxes), function(b) u[box_states_fn(b)])
      readouts <- get_readouts(states_list, p, t)
      du <- numeric(total_states)
      for (b in seq_len(nboxes)) {
        inputs <- get_inputs(b, x, readouts)
        idx <- box_states_fn(b)
        if (length(idx) > 0L) {
          h_box <- function(p2, t2) h(p2, t2)[idx]
          du[idx] <- machines[[b]]@dynamics(
            states_list[[b]], inputs, h_box, p, t
          )
        }
      }
      du
    }
  }

  # Composite readout: wires to output box (-1)
  comp_readout <- function(u, p, t) {
    states_list <- lapply(seq_len(nboxes), function(b) u[box_states_fn(b)])
    readouts <- get_readouts(states_list, p, t)
    out <- numeric(n_total_outputs)
    for (w in d$wires) {
      if (w$target_box == -1L) {
        val <- if (w$source_box == 0L) {
          0
        } else {
          readouts[[w$source_box]][w$source_port]
        }
        out[w$target_port] <- out[w$target_port] + val
      }
    }
    out
  }

  # State names
  all_snames <- character(0)
  for (b in seq_len(nboxes)) {
    sn <- machines[[b]]@state_names
    idx <- box_states_fn(b)
    if (length(sn) == length(idx)) {
      all_snames <- c(all_snames, sn)
    } else {
      all_snames <- c(all_snames, paste0("x", idx))
    }
  }

  Machine(
    system_type = sys_type,
    ninputs = n_total_inputs,
    nstates = total_states,
    noutputs = n_total_outputs,
    dynamics = comp_dynamics,
    readout = comp_readout,
    state_names = all_snames,
    dynamics_expr = merge_machine_expressions(d, machines, all_snames,
                                               offsets, state_counts),
    readout_expr = merge_machine_readout_expr(d, machines, all_snames,
                                               offsets),
    params = unique(unlist(lapply(machines, function(m) m@params)))
  )
}

#' Merge symbolic expressions during oapply_machine composition
#' @keywords internal
merge_machine_expressions <- function(d, machines, all_snames, offsets, state_counts) {
  nboxes <- length(machines)
  has_exprs <- vapply(machines, function(m) length(m@dynamics_expr) > 0L, logical(1))
  if (!all(has_exprs)) return(list())

  result <- list()
  for (b in seq_len(nboxes)) {
    sn <- machines[[b]]@state_names
    dexpr <- machines[[b]]@dynamics_expr
    rexpr <- machines[[b]]@readout_expr

    # Build input substitution: input_k -> wired expression
    input_subs <- list()
    for (w in d$wires) {
      if (w$target_box == b) {
        input_var <- paste0("input_", w$target_port)
        if (w$source_box == 0L) {
          # External input -> ext_input_k
          src_expr <- as.symbol(paste0("ext_input_", w$source_port))
        } else {
          # Readout of another machine
          src_m <- machines[[w$source_box]]
          if (length(src_m@readout_expr) >= w$source_port) {
            src_expr <- src_m@readout_expr[[w$source_port]]
            # Rename source machine's state names if needed
            src_sn <- src_m@state_names
            src_off <- offsets[w$source_box]
            src_map <- character(0)
            for (k in seq_along(src_sn)) {
              gidx <- src_off + k
              if (src_sn[k] != all_snames[gidx])
                src_map[src_sn[k]] <- all_snames[gidx]
            }
            if (length(src_map) > 0L) src_expr <- rename_expr_vars(src_expr, src_map)
          } else {
            src_expr <- as.symbol(paste0("readout_", w$source_box, "_", w$source_port))
          }
        }
        if (input_var %in% names(input_subs)) {
          input_subs[[input_var]] <- call("+", input_subs[[input_var]], src_expr)
        } else {
          input_subs[[input_var]] <- src_expr
        }
      }
    }

    off <- offsets[b]
    for (k in seq_along(sn)) {
      local_nm <- sn[k]
      if (local_nm %in% names(dexpr)) {
        expr <- dexpr[[local_nm]]
        # Substitute input variables
        if (length(input_subs) > 0L)
          expr <- rename_expr_vars(expr, setNames(
            vapply(input_subs, deparse, character(1), width.cutoff = 500L),
            names(input_subs)
          ))
        # Actually, substitution needs to insert expressions, not strings.
        # Use proper substitution:
        expr <- subst_exprs(expr, input_subs)
        global_nm <- all_snames[off + k]
        result[[global_nm]] <- expr
      }
    }
  }
  result
}

#' Substitute symbolic sub-expressions into an expression
#' @keywords internal
subst_exprs <- function(expr, mapping) {
  if (is.symbol(expr)) {
    nm <- as.character(expr)
    if (nm %in% names(mapping)) return(mapping[[nm]])
    return(expr)
  }
  if (is.numeric(expr) || is.character(expr) || is.logical(expr)) return(expr)
  if (is.call(expr)) {
    expr[-1] <- lapply(expr[-1], subst_exprs, mapping = mapping)
    return(expr)
  }
  expr
}

#' Merge readout expressions during oapply_machine composition
#' @keywords internal
merge_machine_readout_expr <- function(d, machines, all_snames, offsets) {
  nboxes <- length(machines)
  has_rexprs <- vapply(machines, function(m) length(m@readout_expr) > 0L, logical(1))

  result <- list()
  for (w in d$wires) {
    if (w$target_box == -1L && w$source_box > 0L) {
      src_m <- machines[[w$source_box]]
      if (length(src_m@readout_expr) >= w$source_port) {
        expr <- src_m@readout_expr[[w$source_port]]
        # Rename
        src_sn <- src_m@state_names
        src_off <- offsets[w$source_box]
        src_map <- character(0)
        for (k in seq_along(src_sn)) {
          gidx <- src_off + k
          if (src_sn[k] != all_snames[gidx])
            src_map[src_sn[k]] <- all_snames[gidx]
        }
        if (length(src_map) > 0L) expr <- rename_expr_vars(expr, src_map)
        idx <- w$target_port
        if (length(result) < idx) {
          for (fill in (length(result) + 1):idx) result[[fill]] <- NULL
        }
        result[[idx]] <- expr
      }
    }
  }
  result
}

# --- Directed wiring diagram (simple representation) ----------------------

#' Create a directed wiring diagram for Machine composition
#'
#' @param ninputs Number of external input ports
#' @param noutputs Number of external output ports
#' @param boxes Named list where each element is c(ninputs, noutputs)
#' @param wires List of wire specs: c(source_box, source_port, target_box, target_port).
#'   Use 0 for external input box, -1 for external output box.
#' @returns A DWD list structure
#' @export
dwd <- function(ninputs, noutputs, boxes, wires) {
  wire_list <- lapply(wires, function(w) {
    list(source_box = as.integer(w[1]), source_port = as.integer(w[2]),
         target_box = as.integer(w[3]), target_port = as.integer(w[4]))
  })
  list(
    ninputs = as.integer(ninputs),
    noutputs = as.integer(noutputs),
    boxes = boxes,
    wires = wire_list
  )
}

# --- Convenience constructors ---------------------------------------------

#' Create a ContinuousResourceSharer
#'
#' @param nstates Number of state variables (inferred from dynamics_expr if omitted)
#' @param dynamics Function (u, p, t) -> du. Auto-built from dynamics_expr if omitted.
#' @param portmap Integer vector mapping ports to states (default: identity)
#' @param state_names Optional character vector of state names
#' @param dynamics_expr Optional named list of quoted R expressions, one per state.
#'   Enables odin2 code generation via \code{rs_to_odin()}.
#' @param params Character vector of parameter names used in dynamics_expr
#' @returns A ResourceSharer with system_type = "continuous"
#' @export
continuous_resource_sharer <- function(nstates = NULL, dynamics = NULL,
                                        portmap = NULL,
                                        state_names = NULL,
                                        dynamics_expr = list(),
                                        params = character(0)) {
  if (length(dynamics_expr) > 0L && is.null(dynamics)) {
    if (is.null(state_names)) state_names <- names(dynamics_expr)
    if (is.null(nstates)) nstates <- length(state_names)
    dynamics <- build_rs_closure(state_names, dynamics_expr, params)
  }
  if (is.null(nstates)) cli::cli_abort("nstates is required")
  if (is.null(dynamics)) cli::cli_abort("dynamics or dynamics_expr is required")
  if (is.null(state_names)) state_names <- paste0("x", seq_len(nstates))
  if (is.null(portmap)) portmap <- seq_len(nstates)
  ResourceSharer(
    system_type = "continuous",
    nstates = as.integer(nstates),
    dynamics = dynamics,
    portmap = as.integer(portmap),
    state_names = state_names,
    dynamics_expr = dynamics_expr,
    params = params
  )
}

#' Create a DiscreteResourceSharer
#' @inheritParams continuous_resource_sharer
#' @export
discrete_resource_sharer <- function(nstates = NULL, dynamics = NULL,
                                      portmap = NULL,
                                      state_names = NULL,
                                      dynamics_expr = list(),
                                      params = character(0)) {
  if (length(dynamics_expr) > 0L && is.null(dynamics)) {
    if (is.null(state_names)) state_names <- names(dynamics_expr)
    if (is.null(nstates)) nstates <- length(state_names)
    dynamics <- build_rs_closure(state_names, dynamics_expr, params)
  }
  if (is.null(nstates)) cli::cli_abort("nstates is required")
  if (is.null(dynamics)) cli::cli_abort("dynamics or dynamics_expr is required")
  if (is.null(state_names)) state_names <- paste0("x", seq_len(nstates))
  if (is.null(portmap)) portmap <- seq_len(nstates)
  ResourceSharer(
    system_type = "discrete",
    nstates = as.integer(nstates),
    dynamics = dynamics,
    portmap = as.integer(portmap),
    state_names = state_names,
    dynamics_expr = dynamics_expr,
    params = params
  )
}

#' Create a DelayResourceSharer from a DDE dynamics function
#'
#' @param nstates Number of state variables
#' @param dynamics Function (u, h, p, t) -> du, where h is a history function
#' @param portmap Integer vector mapping ports to states
#' @param state_names Optional state names
#' @param dynamics_expr Optional named list of expressions (for code gen)
#' @param params Parameter names
#' @returns A ResourceSharer with system_type = "delay"
#' @export
delay_resource_sharer <- function(nstates, dynamics,
                                   portmap = seq_len(nstates),
                                   state_names = paste0("x", seq_len(nstates)),
                                   dynamics_expr = list(),
                                   params = character(0)) {
  ResourceSharer(
    system_type = "delay",
    nstates = as.integer(nstates),
    dynamics = dynamics,
    portmap = as.integer(portmap),
    state_names = state_names,
    dynamics_expr = dynamics_expr,
    params = params
  )
}

#' Create a ContinuousMachine
#'
#' @param ninputs Number of input ports
#' @param nstates Number of internal states (inferred from dynamics_expr if omitted)
#' @param noutputs Number of output ports (inferred from readout_expr if omitted)
#' @param dynamics Function (u, x, p, t) -> du. Auto-built from dynamics_expr if omitted.
#' @param readout Function (u, p, t) -> y. Auto-built from readout_expr if omitted.
#' @param state_names Optional state names
#' @param dynamics_expr Optional named list of quoted expressions per state.
#'   Use \code{input_1}, \code{input_2}, etc. for input port references.
#' @param readout_expr Optional list of quoted expressions, one per output port
#' @param params Character vector of parameter names
#' @export
continuous_machine <- function(ninputs, nstates = NULL, noutputs = NULL,
                                dynamics = NULL, readout = NULL,
                                state_names = NULL,
                                dynamics_expr = list(),
                                readout_expr = list(),
                                params = character(0)) {
  if (length(dynamics_expr) > 0L) {
    if (is.null(state_names)) state_names <- names(dynamics_expr)
    if (is.null(nstates)) nstates <- length(state_names)
    if (is.null(dynamics))
      dynamics <- build_machine_closure(state_names, ninputs, dynamics_expr)
  }
  if (length(readout_expr) > 0L) {
    if (is.null(noutputs)) noutputs <- length(readout_expr)
    if (is.null(readout) && !is.null(state_names))
      readout <- build_readout_closure(state_names, readout_expr)
  }
  if (is.null(nstates)) cli::cli_abort("nstates is required")
  if (is.null(noutputs)) noutputs <- nstates
  if (is.null(state_names)) state_names <- paste0("x", seq_len(nstates))
  if (is.null(dynamics)) cli::cli_abort("dynamics or dynamics_expr is required")
  if (is.null(readout)) readout <- function(u, p, t) u
  Machine(
    system_type = "continuous",
    ninputs = as.integer(ninputs),
    nstates = as.integer(nstates),
    noutputs = as.integer(noutputs),
    dynamics = dynamics,
    readout = readout,
    state_names = state_names,
    dynamics_expr = dynamics_expr,
    readout_expr = readout_expr,
    params = params
  )
}

#' Create a DelayMachine
#'
#' @param ninputs Number of input ports
#' @param nstates Number of internal states
#' @param noutputs Number of output ports
#' @param dynamics Function (u, x, h, p, t) -> du
#' @param readout Function (u, h, p, t) -> y
#' @param state_names Optional state names
#' @param dynamics_expr Optional named list of expressions (for code gen)
#' @param readout_expr Optional list of readout expressions
#' @param params Parameter names
#' @export
delay_machine <- function(ninputs, nstates, noutputs, dynamics, readout,
                           state_names = paste0("x", seq_len(nstates)),
                           dynamics_expr = list(),
                           readout_expr = list(),
                           params = character(0)) {
  Machine(
    system_type = "delay",
    ninputs = as.integer(ninputs),
    nstates = as.integer(nstates),
    noutputs = as.integer(noutputs),
    dynamics = dynamics,
    readout = readout,
    state_names = state_names,
    dynamics_expr = dynamics_expr,
    readout_expr = readout_expr,
    params = params
  )
}

# --- to_odin for ResourceSharer (generate odin2 code) --------------------

#' Generate odin2 code from a ResourceSharer
#'
#' Converts a (possibly composed) ResourceSharer to compilable odin2 code.
#' Requires that the ResourceSharer was created with \code{dynamics_expr}
#' (symbolic expressions). If only a closure is available, code generation
#' is not possible and an error is thrown.
#'
#' @param rs A ResourceSharer with dynamics_expr
#' @param initial Named numeric vector of initial conditions
#' @param params Character vector of parameter names. If NULL, uses rs@params.
#' @param type One of "ode" (continuous ODE), "discrete" (Euler step),
#'   "stochastic" (Euler-Maruyama with Binomial draws)
#' @returns Character string of compilable odin2 code
#' @export
rs_to_odin <- function(rs, initial, params = NULL, type = "ode") {
  state_names <- rs@state_names
  ns <- rs@nstates
  dexpr <- rs@dynamics_expr

  if (length(dexpr) == 0L) {
    cli::cli_abort(c(
      "Cannot generate odin2 code: no symbolic expressions available.",
      "i" = "Use {.arg dynamics_expr} when creating the ResourceSharer.",
      "i" = "Example: continuous_resource_sharer(dynamics_expr = list(S = quote(-beta * S * I)))"
    ))
  }

  if (length(initial) != ns)
    cli::cli_abort("initial must have {ns} elements, got {length(initial)}")
  if (is.null(names(initial)))
    names(initial) <- state_names

  if (is.null(params)) params <- rs@params

  lines <- c("## Auto-generated by algebraicodin (ResourceSharer)", "")

  # Parameter declarations
  if (length(params) > 0L) {
    lines <- c(lines, "## Parameters")
    for (pname in params) {
      lines <- c(lines, sprintf("%s <- parameter()", pname))
    }
    lines <- c(lines, "")
  }

  # Initial conditions
  lines <- c(lines, "## Initial conditions")
  for (i in seq_len(ns)) {
    sn <- state_names[i]
    init_name <- paste0(sn, "0")
    lines <- c(lines, sprintf("%s <- parameter(%s)", init_name, initial[i]))
  }
  lines <- c(lines, "")
  for (i in seq_len(ns)) {
    sn <- state_names[i]
    lines <- c(lines, sprintf("initial(%s) <- %s0", sn, sn))
  }
  lines <- c(lines, "")

  if (type == "ode") {
    lines <- c(lines, "## Derivatives")
    for (i in seq_len(ns)) {
      sn <- state_names[i]
      expr_str <- deparse_odin_expr(dexpr[[sn]])
      lines <- c(lines, sprintf("deriv(%s) <- %s", sn, expr_str))
    }
  } else if (type == "discrete") {
    lines <- c(lines, "## Discrete update (Euler step)")
    lines <- c(lines, "dt <- parameter(1)")
    for (i in seq_len(ns)) {
      sn <- state_names[i]
      expr_str <- deparse_odin_expr(dexpr[[sn]])
      lines <- c(lines, sprintf("update(%s) <- %s + dt * (%s)", sn, sn, expr_str))
    }
  } else if (type == "stochastic") {
    lines <- c(lines, "## Stochastic update")
    lines <- c(lines, "dt <- parameter(0.01)")
    # Attempt to decompose expressions into additive rate terms
    for (i in seq_len(ns)) {
      sn <- state_names[i]
      expr_str <- deparse_odin_expr(dexpr[[sn]])
      lines <- c(lines, sprintf("d_%s <- (%s) * dt", sn, expr_str))
    }
    for (i in seq_len(ns)) {
      sn <- state_names[i]
      lines <- c(lines, sprintf("update(%s) <- %s + d_%s", sn, sn, sn))
    }
  } else {
    cli::cli_abort("Unknown type: {type}. Use 'ode', 'discrete', or 'stochastic'.")
  }

  paste(lines, collapse = "\n")
}

#' Deparse an expression for odin2 code
#' @keywords internal
deparse_odin_expr <- function(expr) {
  if (is.null(expr)) return("0")
  s <- deparse(expr, width.cutoff = 500L)
  s <- paste(s, collapse = " ")
  # Clean up + - → -
  s <- gsub("\\+ -", "- ", s)
  s
}

# --- to_odin for Machine (generate odin2 code) ---------------------------

#' Generate odin2 code from a Machine
#'
#' Converts a (possibly composed) Machine to compilable odin2 code.
#' Requires symbolic expressions (\code{dynamics_expr} and optionally
#' \code{readout_expr}). External inputs are declared as odin2 parameters.
#'
#' @param m A Machine with dynamics_expr
#' @param initial Named numeric vector of initial conditions
#' @param params Character vector of parameter names. If NULL, uses m@params.
#' @param input_params Named list mapping external input variable names to
#'   default values. If NULL, external inputs are declared as parameters
#'   with no default.
#' @param type One of "ode", "discrete"
#' @returns Character string of compilable odin2 code
#' @export
machine_to_odin <- function(m, initial, params = NULL, input_params = NULL,
                             type = "ode") {
  state_names <- m@state_names
  ns <- m@nstates
  dexpr <- m@dynamics_expr

  if (length(dexpr) == 0L) {
    cli::cli_abort(c(
      "Cannot generate odin2 code: no symbolic expressions available.",
      "i" = "Use {.arg dynamics_expr} when creating the Machine."
    ))
  }

  if (length(initial) != ns)
    cli::cli_abort("initial must have {ns} elements, got {length(initial)}")
  if (is.null(names(initial)))
    names(initial) <- state_names

  if (is.null(params)) params <- m@params

  lines <- c("## Auto-generated by algebraicodin (Machine)", "")

  # External input declarations
  if (m@ninputs > 0L) {
    lines <- c(lines, "## External inputs")
    for (i in seq_len(m@ninputs)) {
      inp_name <- paste0("ext_input_", i)
      if (!is.null(input_params) && inp_name %in% names(input_params)) {
        lines <- c(lines, sprintf("%s <- parameter(%s)", inp_name, input_params[[inp_name]]))
      } else {
        lines <- c(lines, sprintf("%s <- parameter()", inp_name))
      }
    }
    lines <- c(lines, "")
  }

  # Parameters
  if (length(params) > 0L) {
    lines <- c(lines, "## Parameters")
    for (pname in params) {
      lines <- c(lines, sprintf("%s <- parameter()", pname))
    }
    lines <- c(lines, "")
  }

  # Initial conditions
  lines <- c(lines, "## Initial conditions")
  for (i in seq_len(ns)) {
    sn <- state_names[i]
    init_name <- paste0(sn, "0")
    lines <- c(lines, sprintf("%s <- parameter(%s)", init_name, initial[i]))
  }
  lines <- c(lines, "")
  for (i in seq_len(ns)) {
    sn <- state_names[i]
    lines <- c(lines, sprintf("initial(%s) <- %s0", sn, sn))
  }
  lines <- c(lines, "")

  if (type == "ode") {
    lines <- c(lines, "## Derivatives")
    for (i in seq_len(ns)) {
      sn <- state_names[i]
      expr_str <- deparse_odin_expr(dexpr[[sn]])
      lines <- c(lines, sprintf("deriv(%s) <- %s", sn, expr_str))
    }
  } else if (type == "discrete") {
    lines <- c(lines, "## Discrete update")
    lines <- c(lines, "dt <- parameter(1)")
    for (i in seq_len(ns)) {
      sn <- state_names[i]
      expr_str <- deparse_odin_expr(dexpr[[sn]])
      lines <- c(lines, sprintf("update(%s) <- %s + dt * (%s)", sn, sn, expr_str))
    }
  }

  # Readout declarations
  rexpr <- m@readout_expr
  if (length(rexpr) > 0L) {
    lines <- c(lines, "", "## Readout (output) variables")
    for (i in seq_along(rexpr)) {
      if (!is.null(rexpr[[i]])) {
        expr_str <- deparse_odin_expr(rexpr[[i]])
        lines <- c(lines, sprintf("output_%d <- %s", i, expr_str))
      }
    }
  }

  paste(lines, collapse = "\n")
}

# --- Simulate a ResourceSharer directly (deSolve) ------------------------

#' Simulate a continuous ResourceSharer using deSolve
#'
#' @param rs A continuous ResourceSharer
#' @param initial Named numeric vector of initial conditions
#' @param times Numeric vector of output times
#' @param params Parameter list (passed to dynamics as p)
#' @returns A data.frame with time and state columns
#' @export
simulate_rs <- function(rs, initial, times, params = list()) {
  if (rs@system_type != "continuous")
    cli::cli_abort("simulate_rs only works with continuous systems")
  if (!requireNamespace("deSolve", quietly = TRUE))
    cli::cli_abort("deSolve is required for simulate_rs")

  # deSolve expects function(t, y, parms)
  odefun <- function(t, y, parms) {
    dy <- rs@dynamics(y, parms, t)
    list(dy)
  }

  sol <- deSolve::ode(y = initial, times = times, func = odefun,
                       parms = params, method = "lsoda")
  as.data.frame(sol)
}

#' Simulate a discrete ResourceSharer
#'
#' @param rs A discrete ResourceSharer
#' @param initial Named numeric vector of initial conditions
#' @param nsteps Number of time steps
#' @param params Parameter list
#' @returns A data.frame with time and state columns
#' @export
simulate_rs_discrete <- function(rs, initial, nsteps, params = list()) {
  if (rs@system_type != "discrete")
    cli::cli_abort("simulate_rs_discrete only works with discrete systems")

  ns <- rs@nstates
  result <- matrix(NA_real_, nrow = nsteps + 1L, ncol = ns + 1L)
  colnames(result) <- c("time", rs@state_names)
  result[1, ] <- c(0, initial)

  u <- initial
  for (step in seq_len(nsteps)) {
    u <- rs@dynamics(u, params, step)
    result[step + 1L, ] <- c(step, u)
  }
  as.data.frame(result)
}

#' Simulate a Machine with external inputs
#'
#' @param m A continuous Machine
#' @param initial Named numeric vector of initial state
#' @param times Numeric vector of output times
#' @param input_fn Function (t) -> numeric vector of inputs
#' @param params Parameter list
#' @returns A data.frame
#' @export
simulate_machine <- function(m, initial, times, input_fn = function(t) numeric(m@ninputs),
                              params = list()) {
  if (m@system_type != "continuous")
    cli::cli_abort("simulate_machine only works with continuous systems")
  if (!requireNamespace("deSolve", quietly = TRUE))
    cli::cli_abort("deSolve is required for simulate_machine")

  odefun <- function(t, y, parms) {
    x <- input_fn(t)
    dy <- m@dynamics(y, x, parms, t)
    list(dy)
  }

  sol <- deSolve::ode(y = initial, times = times, func = odefun,
                       parms = params, method = "lsoda")
  as.data.frame(sol)
}
