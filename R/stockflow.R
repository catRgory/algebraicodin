# Stock-and-Flow Diagrams -----------------------------------------------------
# Port of StockFlow.jl from the AlgebraicJulia ecosystem.
#
# Stock-flow diagrams represent dynamical systems with:
#   - Stocks (S): state variables (accumulators)
#   - Flows (F): rate processes moving material between stocks
#   - Auxiliary/Dynamic Variables (V): intermediate computations
#   - Sum Variables (SV): aggregated sums of stocks (e.g., total population)
#   - Parameters (P): external constants
#   - Links: connections between components (LV, LS, LSV, LVV, LPV)
#
# The key difference from Petri nets: flows have explicit rate expressions
# built from auxiliary variables, which can depend on stocks, parameters,
# and sum variables. This makes the model more explicit and self-documenting.

# --- Schema Definitions ------------------------------------------------------

#' Stock-and-Flow base schema (stocks + sum variables only)
#'
#' Used as the "foot" type for open stock-flow composition.
#' @export
SchStockFlow0 <- acsets::BasicSchema(
  obs = c("S", "SV", "LS"),
  homs = list(
    acsets::hom("lss", "LS", "S"),
    acsets::hom("lssv", "LS", "SV")
  ),
  attrtypes = "Name",
  attrs = list(
    acsets::attr_spec("sname", "S", "Name"),
    acsets::attr_spec("svname", "SV", "Name")
  )
)

#' Full Stock-and-Flow schema
#'
#' Includes stocks, flows, inflows, outflows, auxiliary variables,
#' sum variables, parameters, and all link types.
#' @export
SchStockFlow <- acsets::BasicSchema(
  obs = c("S", "F", "I", "O", "V", "SV", "LS", "LV", "LSV", "P", "LPV"),
  homs = list(
    # Inflow/outflow connections
    acsets::hom("ifn", "I", "F"),      # inflow → flow
    acsets::hom("is", "I", "S"),       # inflow → stock (destination)
    acsets::hom("ofn", "O", "F"),      # outflow → flow
    acsets::hom("os", "O", "S"),       # outflow → stock (source)
    # Flow → auxiliary variable
    acsets::hom("fv", "F", "V"),       # flow rate defined by variable
    # Stock → variable links
    acsets::hom("lvs", "LV", "S"),     # link source stock
    acsets::hom("lvv", "LV", "V"),     # link target variable
    # Stock → sum variable links
    acsets::hom("lss", "LS", "S"),
    acsets::hom("lssv", "LS", "SV"),
    # Sum variable → variable links
    acsets::hom("lsvsv", "LSV", "SV"),
    acsets::hom("lsvv", "LSV", "V"),
    # Parameter → variable links
    acsets::hom("lpvp", "LPV", "P"),
    acsets::hom("lpvv", "LPV", "V")
  ),
  attrtypes = "Name",
  attrs = list(
    acsets::attr_spec("sname", "S", "Name"),
    acsets::attr_spec("fname", "F", "Name"),
    acsets::attr_spec("vname", "V", "Name"),
    acsets::attr_spec("svname", "SV", "Name"),
    acsets::attr_spec("pname", "P", "Name")
  )
)

#' Create StockFlow0 ACSet type (for composition feet)
#' @param ... Arguments passed to the ACSet constructor
#' @export
StockFlow0 <- acsets::acset_type(SchStockFlow0, name = "StockFlow0",
                                  index = c("lss", "lssv"))

#' Create StockFlow ACSet type
#' @param ... Arguments passed to the ACSet constructor
#' @export
StockFlowType <- acsets::acset_type(SchStockFlow, name = "StockFlow",
                                     index = c("is", "os", "ifn", "ofn", "fv",
                                               "lvs", "lvv", "lss", "lssv",
                                               "lsvsv", "lsvv", "lpvp", "lpvv"))

# --- StockFlow S7 class (wraps ACSet + flow rate functions) ------------------

#' A Stock-and-Flow model
#'
#' Wraps a StockFlow ACSet with additional metadata: flow rate functions
#' (expressed as R expressions or functions) and operator information.
#'
#' @param acset The underlying StockFlow ACSet
#' @param flow_fns Named list of flow rate functions: name -> function(u, p, t)
#' @param var_exprs Named list of variable expression strings (for display)
#' @export
StockFlowModel <- S7::new_class("StockFlowModel",
  properties = list(
    acset = S7::class_any,
    flow_fns = S7::class_list,
    var_exprs = S7::class_list
  ),
  constructor = function(acset = NULL, flow_fns = list(), var_exprs = list()) {
    if (is.null(acset)) acset <- StockFlowType()
    S7::new_object(S7::S7_object(),
      acset = acset,
      flow_fns = flow_fns,
      var_exprs = var_exprs
    )
  }
)

# --- Accessor functions ------------------------------------------------------

#' Get stock names
#' @param sf A StockFlowModel or StockFlow ACSet
#' @returns Character vector of stock names
#' @export
sf_snames <- function(sf) {
  ac <- if (S7::S7_inherits(sf, StockFlowModel)) sf@acset else sf
  ns <- acsets::nparts(ac, "S")
  if (ns == 0L) return(character(0))
  vapply(seq_len(ns), function(i) acsets::subpart(ac, i, "sname"), character(1))
}

#' Get flow names
#' @param sf A StockFlowModel or StockFlow ACSet
#' @returns Character vector of flow names
#' @export
sf_fnames <- function(sf) {
  ac <- if (S7::S7_inherits(sf, StockFlowModel)) sf@acset else sf
  nf <- acsets::nparts(ac, "F")
  if (nf == 0L) return(character(0))
  vapply(seq_len(nf), function(i) acsets::subpart(ac, i, "fname"), character(1))
}

#' Get variable names
#' @param sf A StockFlowModel or StockFlow ACSet
#' @returns Character vector of variable names
#' @export
sf_vnames <- function(sf) {
  ac <- if (S7::S7_inherits(sf, StockFlowModel)) sf@acset else sf
  nv <- acsets::nparts(ac, "V")
  if (nv == 0L) return(character(0))
  vapply(seq_len(nv), function(i) acsets::subpart(ac, i, "vname"), character(1))
}

#' Get sum variable names
#' @param sf A StockFlowModel or StockFlow ACSet
#' @returns Character vector of sum variable names
#' @export
sf_svnames <- function(sf) {
  ac <- if (S7::S7_inherits(sf, StockFlowModel)) sf@acset else sf
  nsv <- acsets::nparts(ac, "SV")
  if (nsv == 0L) return(character(0))
  vapply(seq_len(nsv), function(i) acsets::subpart(ac, i, "svname"), character(1))
}

#' Get parameter names
#' @param sf A StockFlowModel or StockFlow ACSet
#' @returns Character vector of parameter names
#' @export
sf_pnames <- function(sf) {
  ac <- if (S7::S7_inherits(sf, StockFlowModel)) sf@acset else sf
  np <- acsets::nparts(ac, "P")
  if (np == 0L) return(character(0))
  vapply(seq_len(np), function(i) acsets::subpart(ac, i, "pname"), character(1))
}

# --- DSL Builder: stock_and_flow() -------------------------------------------

#' Build a stock-and-flow model using a declarative interface
#'
#' @param stocks Character vector of stock names
#' @param flows Named list of flow specifications. Each flow is a list with
#'   \code{from} (source stock or NULL for external), \code{to} (target stock
#'   or NULL for external), and \code{rate} (an R expression or function).
#' @param params Character vector of parameter names
#' @param sums Named list of sum variable definitions:
#'   \code{list(N = c("S", "I", "R"))} means N = S + I + R
#' @returns A StockFlowModel
#' @export
stock_and_flow <- function(stocks, flows, params = character(0),
                            sums = list()) {
  ac <- StockFlowType()

  # Index lookups
  s_idx <- list()
  v_idx <- list()
  sv_idx <- list()
  p_idx <- list()
  f_idx <- list()


  # Add stocks
  for (sn in stocks) {
    s_idx[[sn]] <- acsets::add_part(ac, "S", sname = sn)
  }

  # Add parameters
  for (pn in params) {
    p_idx[[pn]] <- acsets::add_part(ac, "P", pname = pn)
  }

  # Add sum variables
  for (svn in names(sums)) {
    sv_idx[[svn]] <- acsets::add_part(ac, "SV", svname = svn)
  }

  # Link stocks to sum variables
  for (svn in names(sums)) {
    for (sn in sums[[svn]]) {
      if (sn %in% names(s_idx)) {
        acsets::add_part(ac, "LS", lss = s_idx[[sn]], lssv = sv_idx[[svn]])
      }
    }
  }

  # Process flows: extract auxiliary variables and build links
  flow_fns <- list()
  var_exprs <- list()

  for (fn in names(flows)) {
    fspec <- flows[[fn]]
    # Each flow defines an auxiliary variable with the same name prefix
    vn <- paste0("v_", fn)
    v_idx[[vn]] <- acsets::add_part(ac, "V", vname = vn)

    # Add flow
    f_idx[[fn]] <- acsets::add_part(ac, "F", fname = fn, fv = v_idx[[vn]])

    # Add inflow (flow -> target stock)
    if (!is.null(fspec$to) && fspec$to %in% names(s_idx)) {
      acsets::add_part(ac, "I", ifn = f_idx[[fn]], is = s_idx[[fspec$to]])
    }

    # Add outflow (source stock -> flow)
    if (!is.null(fspec$from) && fspec$from %in% names(s_idx)) {
      acsets::add_part(ac, "O", ofn = f_idx[[fn]], os = s_idx[[fspec$from]])
    }

    # Process rate expression to extract dependencies
    rate <- fspec$rate
    if (is.function(rate)) {
      flow_fns[[fn]] <- rate
      var_exprs[[vn]] <- "<function>"
    } else if (is.language(rate) || is.character(rate)) {
      rate_expr <- if (is.character(rate)) parse(text = rate)[[1]] else rate
      var_exprs[[vn]] <- deparse(rate_expr)

      # Extract variable names from expression
      expr_vars <- all.vars(rate_expr)

      # Link stocks referenced in expression
      for (ev in expr_vars) {
        if (ev %in% names(s_idx)) {
          acsets::add_part(ac, "LV", lvs = s_idx[[ev]], lvv = v_idx[[vn]])
        } else if (ev %in% names(sv_idx)) {
          acsets::add_part(ac, "LSV", lsvsv = sv_idx[[ev]], lsvv = v_idx[[vn]])
        } else if (ev %in% names(p_idx)) {
          acsets::add_part(ac, "LPV", lpvp = p_idx[[ev]], lpvv = v_idx[[vn]])
        }
      }

      # Build flow rate function from expression
      flow_fns[[fn]] <- make_flow_fn(rate_expr, stocks, names(sums), params)
    }
  }

  StockFlowModel(
    acset = ac,
    flow_fns = flow_fns,
    var_exprs = var_exprs
  )
}

#' Create a flow specification
#'
#' @param from Source stock name (NULL or "cloud" for external inflow)
#' @param to Target stock name (NULL or "cloud" for external outflow)
#' @param rate Rate expression (quoted) or function(u, p, t)
#' @returns A flow spec list
#' @export
flow <- function(from = NULL, to = NULL, rate) {
  if (identical(from, "cloud")) from <- NULL
  if (identical(to, "cloud")) to <- NULL
  # Capture rate as an unevaluated expression.
  # If the user passes quote(...), the substitute captures the quote() call;
  # we need to evaluate it to get the inner expression.
  rate_expr <- substitute(rate)
  if (is.call(rate_expr) && identical(rate_expr[[1]], quote(quote))) {
    rate_expr <- eval(rate_expr)
  }
  list(from = from, to = to, rate = rate_expr)
}

# Build a flow rate function from an R expression
make_flow_fn <- function(expr, stocks, sv_names, params) {
  # Force all arguments to avoid closure-in-loop capture bugs
  force(expr)
  force(stocks)
  force(sv_names)
  force(params)
  # Create a function that evaluates the expression given (u, p, t)
  # where u is a named vector of stock values, p is a named list of params
  function(u, sv, p, t) {
    env <- new.env(parent = baseenv())
    # Bind stock values
    for (sn in stocks) {
      assign(sn, u[[sn]], envir = env)
    }
    # Bind sum variable values
    for (svn in sv_names) {
      assign(svn, sv[[svn]], envir = env)
    }
    # Bind parameter values
    for (pn in params) {
      assign(pn, p[[pn]], envir = env)
    }
    assign("t", t, envir = env)
    eval(expr, envir = env)
  }
}

# --- Transition Matrices -----------------------------------------------------

#' Compute inflow/outflow transition matrices for a stock-flow model
#'
#' @param sf A StockFlowModel
#' @returns A list with \code{inflow} and \code{outflow} matrices
#'   (rows = flows, cols = stocks)
#' @export
sf_transition_matrices <- function(sf) {
  ac <- if (S7::S7_inherits(sf, StockFlowModel)) sf@acset else sf
  ns <- acsets::nparts(ac, "S")
  nf <- acsets::nparts(ac, "F")
  ni <- acsets::nparts(ac, "I")
  no <- acsets::nparts(ac, "O")

  inflow <- matrix(0L, nrow = nf, ncol = ns)
  outflow <- matrix(0L, nrow = nf, ncol = ns)

  rownames(inflow) <- sf_fnames(sf)
  colnames(inflow) <- sf_snames(sf)
  rownames(outflow) <- sf_fnames(sf)
  colnames(outflow) <- sf_snames(sf)

  for (i in seq_len(ni)) {
    f <- acsets::subpart(ac, i, "ifn")
    s <- acsets::subpart(ac, i, "is")
    inflow[f, s] <- inflow[f, s] + 1L
  }
  for (o in seq_len(no)) {
    f <- acsets::subpart(ac, o, "ofn")
    s <- acsets::subpart(ac, o, "os")
    outflow[f, s] <- outflow[f, s] + 1L
  }

  list(inflow = inflow, outflow = outflow,
       stoichiometry = t(inflow) - t(outflow))
}

# --- Vectorfield Generation --------------------------------------------------

#' Generate a vectorfield function from a stock-flow model
#'
#' @param sf A StockFlowModel
#' @returns A function with signature (t, state, parms) suitable for deSolve
#' @export
sf_vectorfield <- function(sf) {
  ac <- sf@acset
  sn <- sf_snames(sf)
  fn <- sf_fnames(sf)
  ns <- length(sn)
  nf <- length(fn)
  tm <- sf_transition_matrices(sf)
  flow_fns <- sf@flow_fns
  sv_stocks <- sf_sum_variable_stocks(sf)

  function(t, state, parms) {
    # Compute sum variables
    sv <- list()
    for (svn in names(sv_stocks)) {
      sv[[svn]] <- sum(state[sv_stocks[[svn]]])
    }

    # Compute flow rates
    rates <- numeric(nf)
    for (j in seq_len(nf)) {
      rates[j] <- flow_fns[[fn[j]]](state, sv, parms, t)
    }

    # Compute derivatives: dS/dt = sum(inflows) - sum(outflows)
    du <- numeric(ns)
    names(du) <- sn
    for (j in seq_len(nf)) {
      for (i in seq_len(ns)) {
        if (tm$inflow[j, i] == 1L) du[i] <- du[i] + rates[j]
        if (tm$outflow[j, i] == 1L) du[i] <- du[i] - rates[j]
      }
    }

    list(du)
  }
}

# Helper: get which stocks contribute to each sum variable
sf_sum_variable_stocks <- function(sf) {
  ac <- if (S7::S7_inherits(sf, StockFlowModel)) sf@acset else sf
  svn <- sf_svnames(sf)
  sn <- sf_snames(sf)
  nls <- acsets::nparts(ac, "LS")
  result <- list()
  for (sv in svn) result[[sv]] <- character(0)

  for (l in seq_len(nls)) {
    s <- acsets::subpart(ac, l, "lss")
    sv <- acsets::subpart(ac, l, "lssv")
    result[[svn[sv]]] <- c(result[[svn[sv]]], sn[s])
  }
  result
}

# --- odin2 Code Generation ---------------------------------------------------

#' Generate odin2 code from a stock-flow model
#'
#' @param sf A StockFlowModel
#' @param type One of "ode", "stochastic", "discrete"
#' @param initial Named numeric vector of initial conditions (optional)
#' @param compare Optional named list of observation specifications
#'   (see [observe()])
#' @returns Character string of odin2 code
#' @export
sf_to_odin <- function(sf, type = "ode", initial = NULL, compare = NULL) {
  ac <- sf@acset
  sn <- sf_snames(sf)
  fn <- sf_fnames(sf)
  pn <- sf_pnames(sf)
  svn <- sf_svnames(sf)
  ns <- length(sn)
  nf <- length(fn)
  tm <- sf_transition_matrices(sf)
  sv_stocks <- sf_sum_variable_stocks(sf)

  lines <- character(0)

  # Parameters
  for (p in pn) {
    lines <- c(lines, sprintf("%s <- parameter()", p))
  }
  lines <- c(lines, "")

  # Initial conditions
  for (i in seq_len(ns)) {
    init_val <- if (!is.null(initial) && sn[i] %in% names(initial)) {
      initial[sn[i]]
    } else {
      0
    }
    init_name <- paste0(sn[i], "0")
    lines <- c(lines, sprintf("%s <- parameter(%s)", init_name, init_val))
    lines <- c(lines, sprintf("initial(%s) <- %s", sn[i], init_name))
  }
  lines <- c(lines, "")

  # Sum variables
  for (svn_i in svn) {
    sum_expr <- paste(sv_stocks[[svn_i]], collapse = " + ")
    lines <- c(lines, sprintf("%s <- %s", svn_i, sum_expr))
  }
  if (length(svn) > 0) lines <- c(lines, "")

  # Auxiliary variable / flow rate expressions
  # First emit non-flow auxiliary variables (lambda, etc.)
  aux_names <- setdiff(names(sf@var_exprs),
                        paste0("v_", fn))
  for (an in aux_names) {
    expr_str <- sf@var_exprs[[an]]
    if (!is.null(expr_str) && expr_str != "<function>") {
      lines <- c(lines, sprintf("%s <- %s", an, expr_str))
    }
  }
  if (length(aux_names) > 0) lines <- c(lines, "")

  # Then emit flow rate expressions
  for (j in seq_len(nf)) {
    vn <- paste0("v_", fn[j])
    expr_str <- sf@var_exprs[[vn]]
    if (!is.null(expr_str) && expr_str != "<function>") {
      lines <- c(lines, sprintf("%s <- %s", vn, expr_str))
    }
  }
  if (nf > 0) lines <- c(lines, "")

  if (type == "ode") {
    # ODE: deriv() blocks
    for (i in seq_len(ns)) {
      terms <- character(0)
      for (j in seq_len(nf)) {
        vn <- paste0("v_", fn[j])
        if (tm$inflow[j, i] == 1L) terms <- c(terms, paste0("+ ", vn))
        if (tm$outflow[j, i] == 1L) terms <- c(terms, paste0("- ", vn))
      }
      rhs <- if (length(terms) == 0L) "0" else {
        result <- paste(terms, collapse = " ")
        # Clean up leading "+"
        sub("^\\+ ", "", result)
      }
      lines <- c(lines, sprintf("deriv(%s) <- %s", sn[i], rhs))
    }
  } else if (type == "stochastic") {
    # Stochastic: Poisson/Binomial approximation
    # dt is provided by odin2 at system creation, no need to declare it
    for (j in seq_len(nf)) {
      vn <- paste0("v_", fn[j])
      # Find source stock for outflows
      out_stocks <- which(tm$outflow[j, ] == 1L)
      if (length(out_stocks) > 0) {
        src <- sn[out_stocks[1]]
        lines <- c(lines, sprintf("n_%s <- Binomial(%s, min(1, %s / %s * dt))",
                                   fn[j], src, vn, src))
      } else {
        lines <- c(lines, sprintf("n_%s <- Poisson(%s * dt)", fn[j], vn))
      }
    }
    lines <- c(lines, "")
    for (i in seq_len(ns)) {
      terms <- character(0)
      for (j in seq_len(nf)) {
        if (tm$inflow[j, i] == 1L) terms <- c(terms, paste0("+ n_", fn[j]))
        if (tm$outflow[j, i] == 1L) terms <- c(terms, paste0("- n_", fn[j]))
      }
      rhs <- if (length(terms) == 0L) "0" else {
        sub("^\\+ ", "", paste(terms, collapse = " "))
      }
      lines <- c(lines, sprintf("update(%s) <- %s %s", sn[i], sn[i],
                                 if (nchar(rhs) > 0 && !startsWith(rhs, "-"))
                                   paste("+", rhs) else rhs))
    }
  } else if (type == "discrete") {
    # Euler discrete — dt provided by odin2 at system creation
    for (i in seq_len(ns)) {
      terms <- character(0)
      for (j in seq_len(nf)) {
        vn <- paste0("v_", fn[j])
        if (tm$inflow[j, i] == 1L) terms <- c(terms, paste0("+ ", vn))
        if (tm$outflow[j, i] == 1L) terms <- c(terms, paste0("- ", vn))
      }
      rhs <- if (length(terms) == 0L) "0" else {
        sub("^\\+ ", "", paste(terms, collapse = " "))
      }
      lines <- c(lines, sprintf("update(%s) <- %s + (%s) * dt", sn[i], sn[i], rhs))
    }
  }

  # Add comparison/likelihood code
  lines <- c(lines, generate_compare_code(compare, "stockflow"))

  paste(lines, collapse = "\n")
}

#' Generate an odin2 system from a stock-flow model
#'
#' Returns a function that, when called, creates an odin2 system generator.
#'
#' @param sf A StockFlowModel
#' @param type One of "ode", "stochastic", "discrete"
#' @param initial Optional named list of initial values
#' @param compare Optional named list of observation specifications
#' @returns A function that returns an odin2 generator
#' @export
sf_to_odin_system <- function(sf, type = "ode", initial = NULL,
                               compare = NULL) {
  code <- sf_to_odin(sf, type = type, initial = initial, compare = compare)
  if (requireNamespace("odin2", quietly = TRUE)) {
    function() odin2::odin(code)
  } else {
    cli::cli_warn("odin2 not available; returning code string")
    code
  }
}


#' Generate array-based odin2 code from a stratified stock-flow model
#'
#' Produces compact odin2 code using arrays indexed by strata, instead of
#' separate scalar variables per stratum. This scales much better for
#' models with many strata (>5). The model must have been created with
#' [sf_stratify()], which attaches the strata metadata needed.
#'
#' @param sf A stratified StockFlowModel (from [sf_stratify()])
#' @param type One of "ode" (default), "stochastic", "discrete"
#' @param initial Optional named list of initial values per base stock
#'   (vectors of length n_strata)
#' @param compare Optional named list of observation specifications
#' @returns Character string of odin2 code using arrays
#' @export
sf_to_odin_array <- function(sf, type = "ode", initial = NULL,
                              compare = NULL) {
  info <- attr(sf, "strata_info")
  if (is.null(info))
    cli::cli_abort(paste(
      "Model has no strata metadata.",
      "Use {.fn sf_stratify} to create stratified models for array codegen."))

  strata <- info$strata_names
  n_strata <- length(strata)
  base_stocks <- info$base_stocks
  base_flows <- info$base_flows
  base_params <- info$base_params
  cross_flows <- info$cross_strata_flows

  sn <- sf_snames(sf)
  fn <- sf_fnames(sf)
  pn <- sf_pnames(sf)
  tm <- sf_transition_matrices(sf)

  lines <- character()
  lines <- c(lines, "## Auto-generated by algebraicodin (array mode)")
  lines <- c(lines, "")

  # Dimensions
  lines <- c(lines, sprintf("n_strata <- %d", n_strata))
  lines <- c(lines, "")

  # Scalar parameters (shared across strata)
  for (p in base_params) {
    lines <- c(lines, sprintf("%s <- parameter()", p))
  }

  # Cross-strata flow parameters
  cross_params <- unique(unlist(lapply(cross_flows, function(x) x$params)))
  cross_params <- setdiff(cross_params, base_params)
  for (p in cross_params) {
    lines <- c(lines, sprintf("%s <- parameter()", p))
  }

  # Contact matrix parameters (if mixing was applied)
  contact_pars <- pn[grepl("^C_", pn)]
  if (length(contact_pars) > 0) {
    lines <- c(lines, "")
    lines <- c(lines, "## Contact matrix")
    lines <- c(lines, "C <- parameter()")
    lines <- c(lines, "dim(C) <- c(n_strata, n_strata)")
  }
  lines <- c(lines, "")

  # Array dimensions for state variables
  lines <- c(lines, "## State variable dimensions")
  dim_vars <- paste(base_stocks, collapse = ", ")
  lines <- c(lines, sprintf("dim(%s) <- n_strata", dim_vars))
  lines <- c(lines, "")

  # Initial conditions (array parameters)
  lines <- c(lines, "## Initial conditions")
  for (s in base_stocks) {
    init_name <- paste0(s, "0")
    lines <- c(lines, sprintf("%s <- parameter()", init_name))
    lines <- c(lines, sprintf("dim(%s) <- n_strata", init_name))
    lines <- c(lines, sprintf("initial(%s[]) <- %s[i]", s, init_name))
  }
  lines <- c(lines, "")

  # Sum variables
  svn <- sf_svnames(sf)
  base_sum_names <- unique(sub("_[^_]+$", "", svn))
  if (length(base_sum_names) > 0) {
    lines <- c(lines, "## Sum variables")
    for (sv in base_sum_names) {
      # Sum over stocks for each stratum: N[i] = S[i] + I[i] + R[i]
      svs <- sf_sum_variable_stocks(sf)
      # Find which base stocks are in this sum var (from first stratum)
      first_sv <- paste0(sv, "_", strata[1])
      if (first_sv %in% names(svs)) {
        stocks_in_sum <- svs[[first_sv]]
        # Extract base stock names
        base_in_sum <- unique(sub(paste0("_", strata[1], "$"), "", stocks_in_sum))
        sum_expr <- paste(sprintf("%s[i]", base_in_sum), collapse = " + ")
        lines <- c(lines, sprintf("dim(%s) <- n_strata", sv))
        lines <- c(lines, sprintf("%s[] <- %s", sv, sum_expr))
      }
    }
    lines <- c(lines, "")
  }

  # Lambda (force-of-infection) for contact matrix mixing
  lambda_vars <- names(sf@var_exprs)[grepl("^lambda_", names(sf@var_exprs))]
  if (length(lambda_vars) > 0) {
    lines <- c(lines, "## Force of infection (contact matrix)")
    # Group by base flow
    lambda_base <- unique(sub("_[^_]+$", "", lambda_vars))
    for (lb in lambda_base) {
      lines <- c(lines, sprintf("dim(%s) <- n_strata", lb))
      # Build lambda expression: lambda[i] = sum(C[i,j] * I[j] / N[j])
      # Parse from first stratum's expression to get the pattern
      first_expr <- sf@var_exprs[[paste0(lb, "_", strata[1])]]
      # The expression sums C_a_X * I_X / N_X terms → C[i,j] * I[j] / N[j]
      # Detect which stock appears (the infectious compartment)
      inf_stock <- NULL
      norm_var <- NULL
      for (bs in base_stocks) {
        if (grepl(paste0(bs, "_", strata[1]), first_expr)) {
          if (is.null(inf_stock)) inf_stock <- bs
        }
      }
      for (sv in base_sum_names) {
        if (grepl(paste0(sv, "_", strata[1]), first_expr)) {
          norm_var <- sv
        }
      }
      if (!is.null(inf_stock) && !is.null(norm_var)) {
        lines <- c(lines, sprintf(
          "%s[] <- sum(C[i, j] * %s[j] / %s[j])", lb, inf_stock, norm_var))
      } else if (!is.null(inf_stock)) {
        lines <- c(lines, sprintf(
          "%s[] <- sum(C[i, j] * %s[j])", lb, inf_stock))
      }
    }
    lines <- c(lines, "")
  }

  # Flow rate expressions
  lines <- c(lines, "## Flow rates")
  for (bf in base_flows) {
    vn <- paste0("v_", bf)
    # Disease flows: v_<flow>_id_<stratum>
    id_flow <- paste0(bf, "_id_", strata[1])
    flow_expr_key <- paste0("v_", id_flow)

    if (flow_expr_key %in% names(sf@var_exprs)) {
      expr <- sf@var_exprs[[flow_expr_key]]
      # Convert scalar names to array: S_a → S[i], I_a → I[i], N_a → N[i]
      arr_expr <- expr
      for (bs in base_stocks) {
        arr_expr <- gsub(paste0("\\b", bs, "_", strata[1], "\\b"),
                          paste0(bs, "[i]"), arr_expr)
      }
      for (sv in base_sum_names) {
        arr_expr <- gsub(paste0("\\b", sv, "_", strata[1], "\\b"),
                          paste0(sv, "[i]"), arr_expr)
      }
      # Replace lambda references
      for (lb in unique(sub("_[^_]+$", "", lambda_vars))) {
        arr_expr <- gsub(paste0("\\b", lb, "_", strata[1], "\\b"),
                          paste0(lb, "[i]"), arr_expr)
      }

      lines <- c(lines, sprintf("dim(%s) <- n_strata", vn))
      lines <- c(lines, sprintf("%s[] <- %s", vn, arr_expr))
    }
  }
  lines <- c(lines, "")

  # Cross-strata flow rates (aging/migration)
  if (length(cross_flows) > 0) {
    lines <- c(lines, "## Cross-strata flows")
    for (cf in cross_flows) {
      from_idx <- which(strata == cf$from)
      to_idx <- which(strata == cf$to)
      for (bs in base_stocks) {
        cf_name <- paste0(cf$name, "_", bs)
        vn <- paste0("v_", cf_name)
        if (vn %in% names(sf@var_exprs)) {
          expr <- sf@var_exprs[[vn]]
          # Replace strata-specific stock name with array ref
          arr_expr <- gsub(paste0("\\b", bs, "_", cf$from, "\\b"),
                            paste0(bs, "[", from_idx, "]"), expr)
          lines <- c(lines, sprintf("%s <- %s", vn, arr_expr))
        }
      }
    }
    lines <- c(lines, "")
  }

  # Build derivative/update equations
  if (type == "ode") {
    lines <- c(lines, "## Derivatives")
    for (bs in base_stocks) {
      # Disease flow terms (array-indexed)
      terms <- character()
      for (bf in base_flows) {
        id_flow <- paste0(bf, "_id_", strata[1])
        flow_idx <- which(fn == id_flow)
        if (length(flow_idx) > 0) {
          vn <- paste0("v_", bf)
          stock_idx <- which(sn == paste0(bs, "_", strata[1]))
          if (length(stock_idx) > 0) {
            if (tm$inflow[flow_idx, stock_idx] == 1L) terms <- c(terms, paste0("+ ", vn, "[i]"))
            if (tm$outflow[flow_idx, stock_idx] == 1L) terms <- c(terms, paste0("- ", vn, "[i]"))
          }
        }
      }
      rhs <- if (length(terms) == 0L) "0" else {
        sub("^\\+ ", "", paste(terms, collapse = " "))
      }

      # Cross-strata flow terms (scalar, add per-element)
      cross_terms <- character()
      for (cf in cross_flows) {
        cf_name <- paste0(cf$name, "_", bs)
        vn <- paste0("v_", cf_name)
        if (vn %in% names(sf@var_exprs)) {
          from_idx <- which(strata == cf$from)
          to_idx <- which(strata == cf$to)
          cross_terms <- c(cross_terms,
            sprintf("- (if (i == %d) %s else 0)", from_idx, vn),
            sprintf("+ (if (i == %d) %s else 0)", to_idx, vn))
        }
      }

      if (length(cross_terms) > 0) {
        cross_str <- paste(cross_terms, collapse = " ")
        cross_str <- sub("^\\+ ", "", cross_str)
        if (rhs == "0") rhs <- cross_str
        else rhs <- paste(rhs, cross_str, sep = " ")
      }

      lines <- c(lines, sprintf("deriv(%s[]) <- %s", bs, rhs))
    }
  } else if (type == "stochastic") {
    lines <- c(lines, "## Stochastic draws")
    for (bf in base_flows) {
      vn <- paste0("v_", bf)
      # Find outflow stock for Binomial
      id_flow <- paste0(bf, "_id_", strata[1])
      flow_idx <- which(fn == id_flow)
      if (length(flow_idx) > 0) {
        out_stocks <- which(tm$outflow[flow_idx, ] == 1L)
        if (length(out_stocks) > 0) {
          src_name <- sn[out_stocks[1]]
          src_base <- sub(paste0("_", strata[1], "$"), "", src_name)
          lines <- c(lines, sprintf("dim(n_%s) <- n_strata", bf))
          lines <- c(lines, sprintf(
            "n_%s[] <- Binomial(%s[i], min(1, %s[i] / %s[i] * dt))",
            bf, src_base, vn, src_base))
        } else {
          lines <- c(lines, sprintf("dim(n_%s) <- n_strata", bf))
          lines <- c(lines, sprintf("n_%s[] <- Poisson(%s[i] * dt)", bf, vn))
        }
      }
    }

    # Cross-strata draws (scalar)
    for (cf in cross_flows) {
      for (bs in base_stocks) {
        cf_name <- paste0(cf$name, "_", bs)
        vn <- paste0("v_", cf_name)
        if (vn %in% names(sf@var_exprs)) {
          from_idx <- which(strata == cf$from)
          lines <- c(lines, sprintf(
            "n_%s <- Binomial(%s[%d], min(1, %s / %s[%d] * dt))",
            cf_name, bs, from_idx, vn, bs, from_idx))
        }
      }
    }
    lines <- c(lines, "")

    lines <- c(lines, "## Update equations")
    for (bs in base_stocks) {
      terms <- character()
      for (bf in base_flows) {
        id_flow <- paste0(bf, "_id_", strata[1])
        flow_idx <- which(fn == id_flow)
        if (length(flow_idx) > 0) {
          stock_idx <- which(sn == paste0(bs, "_", strata[1]))
          if (length(stock_idx) > 0) {
            if (tm$inflow[flow_idx, stock_idx] == 1L)
              terms <- c(terms, paste0("+ n_", bf, "[i]"))
            if (tm$outflow[flow_idx, stock_idx] == 1L)
              terms <- c(terms, paste0("- n_", bf, "[i]"))
          }
        }
      }

      cross_terms <- character()
      for (cf in cross_flows) {
        cf_name <- paste0(cf$name, "_", bs)
        vn <- paste0("v_", cf_name)
        if (vn %in% names(sf@var_exprs)) {
          from_idx <- which(strata == cf$from)
          to_idx <- which(strata == cf$to)
          cross_terms <- c(cross_terms,
            sprintf("- (if (i == %d) n_%s else 0)", from_idx, cf_name),
            sprintf("+ (if (i == %d) n_%s else 0)", to_idx, cf_name))
        }
      }

      all_terms <- c(terms, cross_terms)
      rhs <- if (length(all_terms) == 0L) "" else {
        t <- paste(all_terms, collapse = " ")
        t <- sub("^\\+ ", "", t)
        if (!startsWith(t, "-")) paste("+", t) else t
      }
      lines <- c(lines, sprintf("update(%s[]) <- %s[i] %s", bs, bs, rhs))
    }
  } else if (type == "discrete") {
    lines <- c(lines, "## Update equations (Euler)")
    for (bs in base_stocks) {
      terms <- character()
      for (bf in base_flows) {
        vn <- paste0("v_", bf)
        id_flow <- paste0(bf, "_id_", strata[1])
        flow_idx <- which(fn == id_flow)
        if (length(flow_idx) > 0) {
          stock_idx <- which(sn == paste0(bs, "_", strata[1]))
          if (length(stock_idx) > 0) {
            if (tm$inflow[flow_idx, stock_idx] == 1L)
              terms <- c(terms, paste0("+ ", vn, "[i]"))
            if (tm$outflow[flow_idx, stock_idx] == 1L)
              terms <- c(terms, paste0("- ", vn, "[i]"))
          }
        }
      }

      cross_terms <- character()
      for (cf in cross_flows) {
        cf_name <- paste0(cf$name, "_", bs)
        vn <- paste0("v_", cf_name)
        if (vn %in% names(sf@var_exprs)) {
          from_idx <- which(strata == cf$from)
          to_idx <- which(strata == cf$to)
          cross_terms <- c(cross_terms,
            sprintf("- (if (i == %d) %s else 0)", from_idx, vn),
            sprintf("+ (if (i == %d) %s else 0)", to_idx, vn))
        }
      }

      all_terms <- c(terms, cross_terms)
      rhs <- if (length(all_terms) == 0L) "0" else {
        sub("^\\+ ", "", paste(all_terms, collapse = " "))
      }
      lines <- c(lines, sprintf("update(%s[]) <- %s[i] + (%s) * dt", bs, bs, rhs))
    }
  }

  lines <- c(lines, "")

  # Output totals (only for ODE — odin2 forbids output() in discrete time)
  if (type == "ode") {
    lines <- c(lines, "## Output totals")
    for (bs in base_stocks) {
      lines <- c(lines, sprintf("output(total_%s) <- sum(%s)", bs, bs))
    }
  }

  # Comparison code
  lines <- c(lines, generate_compare_code(compare, "stockflow"))

  paste(lines, collapse = "\n")
}


#' Generate an array-based odin2 system from a stratified model
#'
#' @param sf A stratified StockFlowModel (from [sf_stratify()])
#' @param type One of "ode", "stochastic", "discrete"
#' @param initial Optional named list of initial values
#' @param compare Optional named list of observation specifications
#' @returns A function that returns an odin2 generator
#' @export
sf_to_odin_array_system <- function(sf, type = "ode", initial = NULL,
                                     compare = NULL) {
  code <- sf_to_odin_array(sf, type = type, initial = initial,
                            compare = compare)
  if (requireNamespace("odin2", quietly = TRUE)) {
    function() odin2::odin(code)
  } else {
    cli::cli_warn("odin2 not available; returning code string")
    code
  }
}

#' Create an open stock-flow model
#'
#' Wraps a StockFlowModel with interface specification for composition.
#' The feet are StockFlow0 ACSets specifying which stocks and sum variables
#' are exposed.
#'
#' @param sf A StockFlowModel
#' @param ... Character vectors of exposed stock names (one per foot/leg)
#' @returns An OpenStockFlow object
#' @export
open_stock_flow <- function(sf, ...) {
  legs <- list(...)
  sn <- sf_snames(sf)
  svn <- sf_svnames(sf)

  leg_maps <- lapply(legs, function(exposed_stocks) {
    # Map exposed stock names to indices in the model
    s_map <- match(exposed_stocks, sn)
    if (any(is.na(s_map))) {
      bad <- exposed_stocks[is.na(s_map)]
      cli::cli_abort("Unknown stocks in leg: {paste(bad, collapse=', ')}")
    }
    s_map
  })

  list(
    apex = sf,
    legs = leg_maps,
    leg_names = legs
  )
}

#' Compose open stock-flow models via a UWD
#'
#' @param d A UWD ACSet (from catlab)
#' @param open_sfs A list of open stock-flow models
#' @returns A StockFlowModel (the composed model)
#' @export
sf_oapply <- function(d, open_sfs) {
  nboxes <- acsets::nparts(d, "Box")
  if (length(open_sfs) != nboxes)
    cli::cli_abort("Need {nboxes} open stock-flows but got {length(open_sfs)}")

  # Collect all stocks, flows, params, sum vars from all components
  all_stocks <- character(0)
  all_flows <- list()
  all_params <- character(0)
  all_sums <- list()
  all_flow_fns <- list()
  all_var_exprs <- list()

  # Track which stocks belong to which box and their local indices
  stock_offsets <- integer(nboxes + 1L)
  stock_offsets[1] <- 0L

  for (b in seq_len(nboxes)) {
    sf <- open_sfs[[b]]$apex
    sn <- sf_snames(sf)
    fn <- sf_fnames(sf)
    pn <- sf_pnames(sf)
    svn <- sf_svnames(sf)
    sv_stocks <- sf_sum_variable_stocks(sf)

    stock_offsets[b + 1L] <- stock_offsets[b] + length(sn)

    # Prefix stock names to avoid collision (will be merged at junctions)
    all_stocks <- c(all_stocks, sn)

    # Copy flows with prefixed names where needed
    for (fname in fn) {
      fspec_orig <- NULL
      # Reconstruct flow spec from the ACSet
      ac <- sf@acset
      f_id <- which(sf_fnames(sf) == fname)

      from_stock <- NULL
      to_stock <- NULL

      # Find outflow source
      no <- acsets::nparts(ac, "O")
      for (o in seq_len(no)) {
        if (acsets::subpart(ac, o, "ofn") == f_id) {
          from_stock <- sn[acsets::subpart(ac, o, "os")]
          break
        }
      }
      # Find inflow target
      ni <- acsets::nparts(ac, "I")
      for (i in seq_len(ni)) {
        if (acsets::subpart(ac, i, "ifn") == f_id) {
          to_stock <- sn[acsets::subpart(ac, i, "is")]
          break
        }
      }

      all_flows[[fname]] <- list(from = from_stock, to = to_stock)
    }

    # Copy flow functions and var expressions
    for (fname in names(sf@flow_fns)) {
      all_flow_fns[[fname]] <- sf@flow_fns[[fname]]
    }
    for (vn in names(sf@var_exprs)) {
      all_var_exprs[[vn]] <- sf@var_exprs[[vn]]
    }

    all_params <- union(all_params, pn)

    # Merge sum variables
    for (svn_i in svn) {
      if (svn_i %in% names(all_sums)) {
        all_sums[[svn_i]] <- union(all_sums[[svn_i]], sv_stocks[[svn_i]])
      } else {
        all_sums[[svn_i]] <- sv_stocks[[svn_i]]
      }
    }
  }

  # Now identify stocks at shared junctions (pushout)
  port_boxes <- acsets::subpart(d, NULL, "box")
  port_junctions <- acsets::subpart(d, NULL, "junction")
  njunctions <- acsets::nparts(d, "Junction")

  # Build junction → stock name mapping
  junction_stocks <- vector("list", njunctions)
  for (j in seq_len(njunctions)) junction_stocks[[j]] <- character(0)

  for (p in seq_len(acsets::nparts(d, "Port"))) {
    b <- port_boxes[p]
    j <- port_junctions[p]
    # Which leg port is this?
    box_port_indices <- which(port_boxes == b)
    local_port_idx <- match(p, box_port_indices)
    # Get the stock name at this leg position
    leg_idx <- open_sfs[[b]]$legs[[1]]  # Assumes single leg per box
    if (local_port_idx <= length(leg_idx)) {
      stock_idx <- leg_idx[local_port_idx]
      stock_name <- sf_snames(open_sfs[[b]]$apex)[stock_idx]
      junction_stocks[[j]] <- c(junction_stocks[[j]], stock_name)
    }
  }

  # Merge: rename duplicate stocks to canonical name
  stock_renames <- list()  # old_name -> canonical_name
  for (j in seq_len(njunctions)) {
    sts <- unique(junction_stocks[[j]])
    if (length(sts) > 1L) {
      canonical <- sts[1]
      for (s in sts[-1]) {
        stock_renames[[s]] <- canonical
      }
    }
  }

  # Apply renames: remove duplicates from stock list
  final_stocks <- character(0)
  for (sn in all_stocks) {
    canonical <- if (sn %in% names(stock_renames)) stock_renames[[sn]] else sn
    if (!(canonical %in% final_stocks)) {
      final_stocks <- c(final_stocks, canonical)
    }
  }

  # Update flow from/to references
  for (fn in names(all_flows)) {
    f <- all_flows[[fn]]
    if (!is.null(f$from) && f$from %in% names(stock_renames)) {
      all_flows[[fn]]$from <- stock_renames[[f$from]]
    }
    if (!is.null(f$to) && f$to %in% names(stock_renames)) {
      all_flows[[fn]]$to <- stock_renames[[f$to]]
    }
  }

  # Update sum variable stock references
  for (svn in names(all_sums)) {
    all_sums[[svn]] <- unique(vapply(all_sums[[svn]], function(s) {
      if (s %in% names(stock_renames)) stock_renames[[s]] else s
    }, character(1)))
  }

  # Reconstruct the composed model using the builder
  # But we need to pass rate expressions, not quoted expressions
  # So we build it manually
  ac <- StockFlowType()
  s_idx <- list()
  v_idx <- list()
  sv_idx <- list()
  p_idx <- list()
  f_idx <- list()

  for (sn in final_stocks) {
    s_idx[[sn]] <- acsets::add_part(ac, "S", sname = sn)
  }
  for (pn in all_params) {
    p_idx[[pn]] <- acsets::add_part(ac, "P", pname = pn)
  }
  for (svn in names(all_sums)) {
    sv_idx[[svn]] <- acsets::add_part(ac, "SV", svname = svn)
    for (sn in all_sums[[svn]]) {
      if (sn %in% names(s_idx)) {
        acsets::add_part(ac, "LS", lss = s_idx[[sn]], lssv = sv_idx[[svn]])
      }
    }
  }

  for (fn in names(all_flows)) {
    fspec <- all_flows[[fn]]
    vn <- paste0("v_", fn)
    v_idx[[vn]] <- acsets::add_part(ac, "V", vname = vn)
    f_idx[[fn]] <- acsets::add_part(ac, "F", fname = fn, fv = v_idx[[vn]])

    if (!is.null(fspec$to) && fspec$to %in% names(s_idx)) {
      acsets::add_part(ac, "I", ifn = f_idx[[fn]], is = s_idx[[fspec$to]])
    }
    if (!is.null(fspec$from) && fspec$from %in% names(s_idx)) {
      acsets::add_part(ac, "O", ofn = f_idx[[fn]], os = s_idx[[fspec$from]])
    }
  }

  StockFlowModel(
    acset = ac,
    flow_fns = all_flow_fns,
    var_exprs = all_var_exprs
  )
}

# --- Conversion: Stock-Flow <-> Petri Net ------------------------------------

#' Convert a stock-flow model to a labelled Petri net
#'
#' Stocks become species, flows become transitions.
#' Inflows/outflows become input/output arcs.
#'
#' @param sf A StockFlowModel
#' @returns A LabelledPetriNet ACSet
#' @export
sf_to_petri <- function(sf) {
  ac <- sf@acset
  sn <- sf_snames(sf)
  fn <- sf_fnames(sf)
  ns <- length(sn)
  nf <- length(fn)

  pn <- LabelledPetriNet()
  s_ids <- integer(ns)
  for (i in seq_len(ns)) {
    s_ids[i] <- acsets::add_part(pn, "S", sname = sn[i])
  }
  t_ids <- integer(nf)
  for (j in seq_len(nf)) {
    t_ids[j] <- acsets::add_part(pn, "T", tname = fn[j])
  }

  # Add input arcs (outflow: stock -> flow = species -> transition)
  no <- acsets::nparts(ac, "O")
  for (o in seq_len(no)) {
    f <- acsets::subpart(ac, o, "ofn")
    s <- acsets::subpart(ac, o, "os")
    acsets::add_part(pn, "I", is = s, it = f)
  }

  # Add output arcs (inflow: flow -> stock = transition -> species)
  ni <- acsets::nparts(ac, "I")
  for (i in seq_len(ni)) {
    f <- acsets::subpart(ac, i, "ifn")
    s <- acsets::subpart(ac, i, "is")
    acsets::add_part(pn, "O", os = s, ot = f)
  }

  pn
}

#' Convert a Petri net to a stock-flow model
#'
#' Creates a basic stock-flow model from a Petri net. Flow rate functions
#' default to mass-action kinetics.
#'
#' @param pn A LabelledPetriNet
#' @returns A StockFlowModel
#' @export
petri_to_sf <- function(pn) {
  sn <- species_names(pn)
  tn <- transition_names(pn)
  ns <- length(sn)
  nt <- length(tn)
  tm <- transition_matrices(pn)

  flows <- list()
  for (j in seq_len(nt)) {
    # Find source and target stocks
    in_species <- which(tm$input[, j] > 0)
    out_species <- which(tm$output[, j] > 0)

    from <- if (length(in_species) > 0) sn[in_species[1]] else NULL
    to_stocks <- setdiff(out_species, in_species)
    to <- if (length(to_stocks) > 0) sn[to_stocks[1]] else {
      if (length(out_species) > 0) sn[out_species[1]] else NULL
    }

    # Mass-action rate expression
    rate_parts <- c(tn[j])  # rate constant name
    for (s in in_species) {
      for (k in seq_len(tm$input[s, j])) {
        rate_parts <- c(rate_parts, sn[s])
      }
    }
    rate_str <- paste(rate_parts, collapse = " * ")

    flows[[tn[j]]] <- list(
      from = from,
      to = to,
      rate = parse(text = rate_str)[[1]]
    )
  }

  stock_and_flow(
    stocks = sn,
    flows = flows,
    params = tn
  )
}

# --- Stock-Flow to ResourceSharer -------------------------------------------

#' Convert a stock-flow model to a continuous ResourceSharer
#'
#' @param sf A StockFlowModel
#' @returns A ResourceSharer with system_type = "continuous"
#' @export
sf_to_resource_sharer <- function(sf) {
  vf <- sf_vectorfield(sf)
  sn <- sf_snames(sf)
  ns <- length(sn)

  dynamics <- function(u, p, t) {
    names(u) <- sn
    result <- vf(t, u, p)
    unname(result[[1]])
  }

  continuous_resource_sharer(
    nstates = ns,
    dynamics = dynamics,
    state_names = sn
  )
}

# --- Simulation convenience --------------------------------------------------

#' Simulate a stock-flow model using deSolve
#'
#' @param sf A StockFlowModel
#' @param initial Named numeric vector of initial conditions
#' @param times Numeric vector of output times
#' @param params Named list of parameter values
#' @returns A data.frame
#' @export
simulate_sf <- function(sf, initial, times, params) {
  if (!requireNamespace("deSolve", quietly = TRUE))
    cli::cli_abort("deSolve is required for simulate_sf")

  vf <- sf_vectorfield(sf)
  sol <- deSolve::ode(y = initial, times = times, func = vf,
                       parms = params, method = "lsoda")
  as.data.frame(sol)
}

# --- Visualization -----------------------------------------------------------

#' Plot a stock-flow diagram using DiagrammeR
#'
#' @param sf A StockFlowModel
#' @returns A DiagrammeR graph
#' @export
plot_stock_flow <- function(sf) {
  if (!requireNamespace("DiagrammeR", quietly = TRUE))
    cli::cli_abort("DiagrammeR is required for plot_stock_flow")

  sn <- sf_snames(sf)
  fn <- sf_fnames(sf)
  svn <- sf_svnames(sf)
  pn <- sf_pnames(sf)
  ac <- sf@acset

  # Build DOT string
  dot_lines <- c(
    "digraph StockFlow {",
    "  rankdir=LR;",
    "  node [fontname=Helvetica];",
    ""
  )

  # Stocks: rectangles
  for (s in sn) {
    dot_lines <- c(dot_lines,
      sprintf('  "%s" [shape=box, style=filled, fillcolor="#AED6F1"];', s))
  }

  # Flows: small circles
  for (f in fn) {
    dot_lines <- c(dot_lines,
      sprintf('  "%s" [shape=circle, width=0.3, style=filled, fillcolor="#F9E79F", label="%s"];', f, f))
  }

  # Sum variables: diamonds
  for (sv in svn) {
    dot_lines <- c(dot_lines,
      sprintf('  "%s" [shape=diamond, style=filled, fillcolor="#ABEBC6"];', sv))
  }

  # Parameters: ellipses
  for (p in pn) {
    dot_lines <- c(dot_lines,
      sprintf('  "%s" [shape=ellipse, style=filled, fillcolor="#D7BDE2"];', p))
  }

  dot_lines <- c(dot_lines, "")

  # Outflow edges: stock -> flow
  no <- acsets::nparts(ac, "O")
  for (o in seq_len(no)) {
    f <- acsets::subpart(ac, o, "ofn")
    s <- acsets::subpart(ac, o, "os")
    dot_lines <- c(dot_lines,
      sprintf('  "%s" -> "%s" [color="#E74C3C"];', sn[s], fn[f]))
  }

  # Inflow edges: flow -> stock
  ni <- acsets::nparts(ac, "I")
  for (i in seq_len(ni)) {
    f <- acsets::subpart(ac, i, "ifn")
    s <- acsets::subpart(ac, i, "is")
    dot_lines <- c(dot_lines,
      sprintf('  "%s" -> "%s" [color="#2ECC71"];', fn[f], sn[s]))
  }

  # LS links: stock -> sum variable (dashed)
  nls <- acsets::nparts(ac, "LS")
  for (l in seq_len(nls)) {
    s <- acsets::subpart(ac, l, "lss")
    sv <- acsets::subpart(ac, l, "lssv")
    dot_lines <- c(dot_lines,
      sprintf('  "%s" -> "%s" [style=dashed, color="#95A5A6"];', sn[s], svn[sv]))
  }

  dot_lines <- c(dot_lines, "}")
  dot_str <- paste(dot_lines, collapse = "\n")

  DiagrammeR::grViz(dot_str)
}

# --- Stock-flow stratification -----------------------------------------------
# Port of StockFlow.jl's stratification.
#
# Like Petri net stratification, this computes a pullback: given a base SF
# and a strata SF, both mapped to a common type template, the result
# contains pairs (base_obj, strata_obj) where both map to the same type.
#
# The type template defines the "shape" — e.g., one stock type (pop),
# one flow type (transition), one parameter type (rate).

#' Typed stock-flow model for stratification
#'
#' @param sf A StockFlowModel
#' @param type_sf A StockFlowModel acting as the type template
#' @param stock_map Named character: base stock name -> type stock name
#' @param flow_map Named character: base flow name -> type flow name
#' @param param_map Named character: base param name -> type param name
#' @param refl_flows Character vector of reflexive flow names
#' @returns A TypedStockFlow object
#' @export
TypedStockFlow <- S7::new_class("TypedStockFlow",
  properties = list(
    sf = StockFlowModel,
    type_sf = StockFlowModel,
    stock_map = S7::class_integer,
    flow_map = S7::class_integer,
    param_map = S7::class_integer,
    refl_flows = S7::class_character
  ),
  constructor = function(sf, type_sf, stock_map, flow_map, param_map,
                         refl_flows = character(0)) {
    nf <- length(flow_map)
    if (length(refl_flows) == 0L) refl_flows <- rep("", nf)
    S7::new_object(S7::S7_object(),
      sf = sf,
      type_sf = type_sf,
      stock_map = as.integer(stock_map),
      flow_map = as.integer(flow_map),
      param_map = as.integer(param_map),
      refl_flows = refl_flows
    )
  }
)

#' Create a typed stock-flow model
#'
#' Maps stocks, flows, and parameters of a stock-flow model to a type template.
#'
#' @param sf A StockFlowModel
#' @param type_sf A StockFlowModel acting as type template
#' @param stock_types Named character: sf stock name -> type stock name
#' @param flow_types Named character: sf flow name -> type flow name
#' @param param_types Named character: sf param name -> type param name
#' @returns A TypedStockFlow
#' @export
sf_typed <- function(sf, type_sf, stock_types = NULL, flow_types = NULL,
                     param_types = NULL) {
  sn <- sf_snames(sf)
  fn <- sf_fnames(sf)
  pn <- sf_pnames(sf)
  type_sn <- sf_snames(type_sf)
  type_fn <- sf_fnames(type_sf)
  type_pn <- sf_pnames(type_sf)

  # Build stock map
  if (is.null(stock_types)) {
    if (length(type_sn) == 1L) {
      stock_map <- rep(1L, length(sn))
    } else {
      cli::cli_abort("stock_types required when type has multiple stocks")
    }
  } else {
    stock_map <- integer(length(sn))
    for (i in seq_along(sn)) {
      tn <- stock_types[sn[i]]
      if (is.na(tn)) cli::cli_abort("No type for stock '{sn[i]}'")
      stock_map[i] <- match(tn, type_sn)
      if (is.na(stock_map[i])) cli::cli_abort("Unknown type stock '{tn}'")
    }
  }

  # Build flow map
  if (is.null(flow_types)) {
    if (length(type_fn) == 1L) {
      flow_map <- rep(1L, length(fn))
    } else {
      cli::cli_abort("flow_types required when type has multiple flows")
    }
  } else {
    flow_map <- integer(length(fn))
    for (i in seq_along(fn)) {
      tn <- flow_types[fn[i]]
      if (is.na(tn)) cli::cli_abort("No type for flow '{fn[i]}'")
      flow_map[i] <- match(tn, type_fn)
      if (is.na(flow_map[i])) cli::cli_abort("Unknown type flow '{tn}'")
    }
  }

  # Build param map
  if (is.null(param_types)) {
    if (length(type_pn) == 1L) {
      param_map <- rep(1L, length(pn))
    } else if (length(type_pn) == 0L) {
      param_map <- integer(0)
    } else {
      cli::cli_abort("param_types required when type has multiple params")
    }
  } else {
    param_map <- integer(length(pn))
    for (i in seq_along(pn)) {
      tn <- param_types[pn[i]]
      if (is.na(tn)) cli::cli_abort("No type for param '{pn[i]}'")
      param_map[i] <- match(tn, type_pn)
      if (is.na(param_map[i])) cli::cli_abort("Unknown type param '{tn}'")
    }
  }

  TypedStockFlow(sf, type_sf, stock_map, flow_map, param_map)
}

#' Add reflexive (identity) flows for stock-flow stratification
#'
#' Adds self-loop flows for each stock, typed to a given flow type.
#' Needed so the typed product can pair real flows from one model
#' with identity flows on corresponding stocks of the other.
#'
#' @param tsf A TypedStockFlow
#' @param reflexives Named list: stock_name -> character vector of flow type names
#' @returns A new TypedStockFlow with reflexive flows added
#' @export
sf_add_reflexives <- function(tsf, reflexives) {
  sf <- tsf@sf
  type_sf <- tsf@type_sf
  # Deep-copy the ACSet to avoid mutating the original model
  ac <- acsets::copy_acset(sf@acset)
  sn <- sf_snames(sf)
  type_fn <- sf_fnames(type_sf)

  new_flow_map <- tsf@flow_map
  new_refl <- tsf@refl_flows
  new_flow_fns <- sf@flow_fns
  new_var_exprs <- sf@var_exprs

  for (stock_name in names(reflexives)) {
    si <- match(stock_name, sn)
    if (is.na(si)) cli::cli_abort("Stock '{stock_name}' not found")

    for (ftype_name in reflexives[[stock_name]]) {
      fti <- match(ftype_name, type_fn)
      if (is.na(fti)) cli::cli_abort("Flow type '{ftype_name}' not found")

      refl_name <- paste0("refl_", stock_name, "_", ftype_name)
      vn <- paste0("v_", refl_name)

      # Add variable for the reflexive flow
      vi <- acsets::add_part(ac, "V", vname = vn)
      # Add flow with from=stock, to=stock (self-loop)
      fi <- acsets::add_part(ac, "F", fname = refl_name, fv = vi)
      acsets::add_part(ac, "I", ifn = fi, is = si)
      acsets::add_part(ac, "O", ofn = fi, os = si)

      new_flow_map <- c(new_flow_map, fti)
      new_refl <- c(new_refl, stock_name)
      new_var_exprs[[vn]] <- "0"
    }
  }

  new_sf <- StockFlowModel(acset = ac, flow_fns = new_flow_fns,
                             var_exprs = new_var_exprs)
  TypedStockFlow(new_sf, type_sf, tsf@stock_map, new_flow_map,
                  tsf@param_map, new_refl)
}

#' Compute typed product of two typed stock-flow models (stratification)
#'
#' Creates all valid pairs of (base_obj, strata_obj) where both map to
#' the same type object, producing a stratified stock-flow model.
#'
#' @param tsf1 A TypedStockFlow (e.g., base disease model)
#' @param tsf2 A TypedStockFlow (e.g., strata model like age groups)
#' @param independent_params If TRUE, keep parameters from both models
#'   independent (union) instead of pairing them via type mapping. Default FALSE.
#' @returns A new StockFlowModel with stratified stocks and flows
#' @export
sf_typed_product <- function(tsf1, tsf2, independent_params = FALSE) {
  sf1 <- tsf1@sf
  sf2 <- tsf2@sf
  ac1 <- sf1@acset
  ac2 <- sf2@acset

  sn1 <- sf_snames(sf1); sn2 <- sf_snames(sf2)
  fn1 <- sf_fnames(sf1); fn2 <- sf_fnames(sf2)
  pn1 <- sf_pnames(sf1); pn2 <- sf_pnames(sf2)
  refl1 <- tsf1@refl_flows; refl2 <- tsf2@refl_flows

  ns1 <- length(sn1); ns2 <- length(sn2)
  nf1 <- length(fn1); nf2 <- length(fn2)
  np1 <- length(pn1); np2 <- length(pn2)

  # --- Product stocks ---
  result_stocks <- character(0)
  sp_pair <- matrix(0L, ns1, ns2)
  for (s1 in seq_len(ns1)) {
    for (s2 in seq_len(ns2)) {
      if (tsf1@stock_map[s1] == tsf2@stock_map[s2]) {
        name <- paste0(sn1[s1], "_", sn2[s2])
        result_stocks <- c(result_stocks, name)
        sp_pair[s1, s2] <- length(result_stocks)
      }
    }
  }

  # --- Product flows ---
  result_flows <- list()
  for (f1 in seq_len(nf1)) {
    for (f2 in seq_len(nf2)) {
      if (tsf1@flow_map[f1] != tsf2@flow_map[f2]) next

      r1 <- refl1[f1]; r2 <- refl2[f2]
      if (nchar(r1) > 0 && nchar(r2) > 0) next  # skip refl-refl pairs

      # Smart naming
      if (nchar(r1) > 0 && nchar(r2) == 0) {
        fname <- paste0(fn2[f2], "_", r1)
      } else if (nchar(r1) == 0 && nchar(r2) > 0) {
        fname <- paste0(fn1[f1], "_", r2)
      } else {
        fname <- paste0(fn1[f1], "_", fn2[f2])
      }

      # Determine from/to by pairing inflow/outflow stocks
      from1 <- get_flow_source(ac1, f1)
      to1 <- get_flow_target(ac1, f1)
      from2 <- get_flow_source(ac2, f2)
      to2 <- get_flow_target(ac2, f2)

      # For reflexive flows, from/to are the reflexive stock
      if (nchar(r1) > 0) { from1 <- match(r1, sn1); to1 <- from1 }
      if (nchar(r2) > 0) { from2 <- match(r2, sn2); to2 <- from2 }

      from_pair <- if (!is.na(from1) && !is.na(from2)) sp_pair[from1, from2] else NA_integer_
      to_pair <- if (!is.na(to1) && !is.na(to2)) sp_pair[to1, to2] else NA_integer_

      # Skip invalid pairs (unmapped stock combinations)
      if ((!is.na(from_pair) && from_pair == 0L) ||
          (!is.na(to_pair) && to_pair == 0L)) next

      # Build rate expression by substitution
      expr <- build_stratified_flow_expr(
        sf1, sf2, fn1[f1], fn2[f2], r1, r2,
        sn1, sn2, pn1, pn2, sp_pair,
        independent_params = independent_params)

      result_flows[[fname]] <- list(
        from = if (is.na(from_pair)) NULL else result_stocks[from_pair],
        to = if (is.na(to_pair)) NULL else result_stocks[to_pair],
        rate_expr = expr
      )
    }
  }

  # --- Product parameters ---
  result_params <- character(0)
  if (independent_params) {
    # Keep params independent: union of both sets (no renaming)
    result_params <- unique(c(pn1, pn2))
  } else if (np1 > 0 && np2 > 0) {
    pp_pair <- matrix("", np1, np2)
    for (p1 in seq_len(np1)) {
      for (p2 in seq_len(np2)) {
        if (tsf1@param_map[p1] == tsf2@param_map[p2]) {
          name <- paste0(pn1[p1], "_", pn2[p2])
          result_params <- c(result_params, name)
          pp_pair[p1, p2] <- name
        }
      }
    }
  } else if (np1 > 0 && np2 == 0) {
    # Strata has no params: pass base params through directly
    result_params <- pn1
  } else if (np1 == 0 && np2 > 0) {
    result_params <- pn2
  }

  # --- Product sum variables ---
  svn1 <- sf_svnames(sf1); svn2 <- sf_svnames(sf2)
  result_sums <- list()
  # Each sum variable in the base gets replicated per stratum
  for (sv1 in svn1) {
    ls1 <- sf_sum_variable_stocks(sf1)[[sv1]]
    for (sv2 in svn2) {
      ls2 <- sf_sum_variable_stocks(sf2)[[sv2]]
      sv_name <- paste0(sv1, "_", sv2)
      sv_stocks <- character(0)
      for (s1_name in ls1) {
        for (s2_name in ls2) {
          s1i <- match(s1_name, sn1); s2i <- match(s2_name, sn2)
          if (!is.na(s1i) && !is.na(s2i) && sp_pair[s1i, s2i] > 0)
            sv_stocks <- c(sv_stocks, result_stocks[sp_pair[s1i, s2i]])
        }
      }
      if (length(sv_stocks) > 0)
        result_sums[[sv_name]] <- sv_stocks
    }
  }

  # If either model has sums and the other doesn't, create per-stratum sums
  if (length(svn1) > 0 && length(svn2) == 0) {
    for (sv1 in svn1) {
      ls1 <- sf_sum_variable_stocks(sf1)[[sv1]]
      for (s2 in seq_len(ns2)) {
        sv_name <- paste0(sv1, "_", sn2[s2])
        sv_stocks <- character(0)
        for (s1_name in ls1) {
          s1i <- match(s1_name, sn1)
          if (!is.na(s1i) && sp_pair[s1i, s2] > 0)
            sv_stocks <- c(sv_stocks, result_stocks[sp_pair[s1i, s2]])
        }
        if (length(sv_stocks) > 0)
          result_sums[[sv_name]] <- sv_stocks
      }
    }
  }
  if (length(svn2) > 0 && length(svn1) == 0) {
    for (sv2 in svn2) {
      ls2 <- sf_sum_variable_stocks(sf2)[[sv2]]
      for (s1 in seq_len(ns1)) {
        sv_name <- paste0(sn1[s1], "_", sv2)
        sv_stocks <- character(0)
        for (s2_name in ls2) {
          s2i <- match(s2_name, sn2)
          if (!is.na(s2i) && sp_pair[s1, s2i] > 0)
            sv_stocks <- c(sv_stocks, result_stocks[sp_pair[s1, s2i]])
        }
        if (length(sv_stocks) > 0)
          result_sums[[sv_name]] <- sv_stocks
      }
    }
  }

  # Build result model
  flow_specs <- lapply(names(result_flows), function(fname) {
    fl <- result_flows[[fname]]
    spec <- list(from = fl$from, to = fl$to, rate = fl$rate_expr)
    spec
  })
  names(flow_specs) <- names(result_flows)

  # Use stock_and_flow_raw to build from already-processed data
  sf_build_stratified(result_stocks, flow_specs, result_params, result_sums)
}

#' Get flow source stock index
#' @keywords internal
get_flow_source <- function(ac, f_idx) {
  outflows <- which(acsets::subpart(ac, NULL, "ofn") == f_idx)
  if (length(outflows) == 0L) return(NA_integer_)
  acsets::subpart(ac, outflows[1], "os")
}

#' Get flow target stock index
#' @keywords internal
get_flow_target <- function(ac, f_idx) {
  inflows <- which(acsets::subpart(ac, NULL, "ifn") == f_idx)
  if (length(inflows) == 0L) return(NA_integer_)
  acsets::subpart(ac, inflows[1], "is")
}

#' Build stratified flow rate expression
#' @keywords internal
build_stratified_flow_expr <- function(sf1, sf2, fn1, fn2, r1, r2,
                                        sn1, sn2, pn1, pn2, sp_pair,
                                        independent_params = FALSE) {
  # Determine which model provides the "real" expression
  if (nchar(r1) > 0) {
    # Model 1 is reflexive; use model 2's expression
    real_sf <- sf2; real_fn <- fn2
    refl_stock <- r1
  } else if (nchar(r2) > 0) {
    # Model 2 is reflexive; use model 1's expression
    real_sf <- sf1; real_fn <- fn1
    refl_stock <- r2
  } else {
    # Both are real flows — use model 1's expression
    vn1 <- paste0("v_", fn1)
    expr_str <- sf1@var_exprs[[vn1]]
    if (is.null(expr_str) || expr_str == "<function>") {
      return(parse(text = paste0("v_", fn1, "_", fn2))[[1]])
    }
    expr <- parse(text = expr_str)[[1]]
    # Determine which stratum this flow pair corresponds to
    from2 <- get_flow_source(sf2@acset, match(fn2, sf_fnames(sf2)))
    to2 <- get_flow_target(sf2@acset, match(fn2, sf_fnames(sf2)))
    s2i <- if (!is.na(from2)) from2 else if (!is.na(to2)) to2 else 1L
    s2_name <- sn2[s2i]

    all_vars <- all.vars(expr)
    mapping <- list()
    for (v in all_vars) {
      if (v %in% sn1) {
        mapping[[v]] <- as.symbol(paste0(v, "_", s2_name))
      } else if (!independent_params && v %in% pn1 && length(pn2) > 0) {
        # Only rename params if not independent and the other model has params
        mapping[[v]] <- as.symbol(paste0(v, "_", pn2[1]))
      } else if (v %in% sf_svnames(sf1)) {
        mapping[[v]] <- as.symbol(paste0(v, "_", s2_name))
      }
      # If independent_params or v is a param but strata has no params, leave as-is
    }
    if (length(mapping) > 0) expr <- subst_exprs(expr, mapping)
    return(expr)
  }

  # Handle reflexive case (one model is reflexive, use the other's expression)
  vn <- paste0("v_", real_fn)
  expr_str <- real_sf@var_exprs[[vn]]
  if (is.null(expr_str) || expr_str == "<function>") {
    return(parse(text = paste0("v_", real_fn, "_", refl_stock))[[1]])
  }
  expr <- parse(text = expr_str)[[1]]

  all_vars <- all.vars(expr)
  real_sn <- if (nchar(r1) > 0) sn2 else sn1
  real_pn <- if (nchar(r1) > 0) pn2 else pn1
  refl_pn <- if (nchar(r1) > 0) pn1 else pn2

  mapping <- list()
  # Determine naming order: product names are always <sn1>_<sn2>
  # When model 1 is reflexive (r1 non-empty): real_sf is model 2, so v is sn2
  #   → product name = refl_stock (from sn1) "_" v (from sn2)
  # When model 2 is reflexive (r2 non-empty): real_sf is model 1, so v is sn1
  #   → product name = v (from sn1) "_" refl_stock (from sn2)
  m1_is_refl <- nchar(r1) > 0
  for (v in all_vars) {
    if (v %in% real_sn) {
      if (m1_is_refl) {
        mapping[[v]] <- as.symbol(paste0(refl_stock, "_", v))
      } else {
        mapping[[v]] <- as.symbol(paste0(v, "_", refl_stock))
      }
    } else if (!independent_params && v %in% real_pn && length(refl_pn) > 0) {
      mapping[[v]] <- as.symbol(paste0(v, "_", refl_pn[1]))
    } else if (v %in% sf_svnames(real_sf)) {
      if (m1_is_refl) {
        mapping[[v]] <- as.symbol(paste0(refl_stock, "_", v))
      } else {
        mapping[[v]] <- as.symbol(paste0(v, "_", refl_stock))
      }
    }
    # If independent_params or v is a param but other model has no params, leave as-is
  }
  if (length(mapping) > 0) expr <- subst_exprs(expr, mapping)
  expr
}

#' Build a stratified StockFlowModel from raw components
#' @keywords internal
sf_build_stratified <- function(stocks, flow_specs, params, sums) {
  ac <- StockFlowType()
  s_idx <- list()
  f_idx <- list()
  p_idx <- list()
  sv_idx <- list()
  v_idx <- list()
  flow_fns <- list()
  var_exprs <- list()

  # Add stocks
  for (sn in stocks) {
    s_idx[[sn]] <- acsets::add_part(ac, "S", sname = sn)
  }

  # Add parameters
  for (pn in params) {
    p_idx[[pn]] <- acsets::add_part(ac, "P", pname = pn)
  }

  # Add sum variables
  for (svn in names(sums)) {
    sv_idx[[svn]] <- acsets::add_part(ac, "SV", svname = svn)
    for (sn in sums[[svn]]) {
      if (sn %in% names(s_idx))
        acsets::add_part(ac, "LS", lss = s_idx[[sn]], lssv = sv_idx[[svn]])
    }
  }

  # Add flows
  for (fname in names(flow_specs)) {
    fspec <- flow_specs[[fname]]
    vn <- paste0("v_", fname)
    v_idx[[vn]] <- acsets::add_part(ac, "V", vname = vn)
    f_idx[[fname]] <- acsets::add_part(ac, "F", fname = fname, fv = v_idx[[vn]])

    if (!is.null(fspec$to) && fspec$to %in% names(s_idx))
      acsets::add_part(ac, "I", ifn = f_idx[[fname]], is = s_idx[[fspec$to]])
    if (!is.null(fspec$from) && fspec$from %in% names(s_idx))
      acsets::add_part(ac, "O", ofn = f_idx[[fname]], os = s_idx[[fspec$from]])

    rate_expr <- fspec$rate
    if (is.language(rate_expr)) {
      var_exprs[[vn]] <- deparse(rate_expr, width.cutoff = 500L)
      # Link vars referenced in expression
      expr_vars <- all.vars(rate_expr)
      for (ev in expr_vars) {
        if (ev %in% names(s_idx))
          acsets::add_part(ac, "LV", lvs = s_idx[[ev]], lvv = v_idx[[vn]])
        else if (ev %in% names(sv_idx))
          acsets::add_part(ac, "LSV", lsvsv = sv_idx[[ev]], lsvv = v_idx[[vn]])
        else if (ev %in% names(p_idx))
          acsets::add_part(ac, "LPV", lpvp = p_idx[[ev]], lpvv = v_idx[[vn]])
      }
      flow_fns[[fname]] <- make_flow_fn(rate_expr, stocks, names(sums), params)
    }
  }

  StockFlowModel(acset = ac, flow_fns = flow_fns, var_exprs = var_exprs)
}

#' SF type template for simple infectious disease models
#'
#' One stock type (pop), with flow types for disease transitions and
#' strata transitions. Analogous to \code{infectious_ontology()}.
#'
#' @returns A StockFlowModel usable as a type template
#' @export
sf_infectious_type <- function() {
  stock_and_flow(
    stocks = c("pop"),
    flows = list(
      disease = flow(from = "pop", to = "pop", rate = quote(0)),
      strata = flow(from = "pop", to = "pop", rate = quote(0))
    ),
    params = c("rate")
  )
}

#' User-friendly stock-flow stratification
#'
#' Stratifies a stock-flow model by strata (e.g., age groups, locations).
#'
#' @param base A StockFlowModel (the base disease model)
#' @param strata_names Character vector of strata names (e.g., c("child", "adult", "senior"))
#' @param flow_types Named character: flow name -> "disease" or "strata"
#' @param cross_strata_flows Optional list of cross-strata flow specs. Each is a list
#'   with \code{from}, \code{to} (strata names), \code{rate} (expression).
#' @param param_types Optional named character: param name -> type param name
#' @param mixing Optional named list specifying contact-matrix mixing for flows.
#'   Names are base flow names (e.g., "infection"). Values are either:
#'   \itemize{
#'     \item A string: the contact matrix parameter prefix (e.g., "C").
#'       Contact parameters will be named \code{C_<from>_<to>}.
#'     \item A list with elements: \code{contact} (prefix), \code{type}
#'       ("frequency" or "density", default "frequency").
#'   }
#'   When mixing is applied, the flow rate uses a force-of-infection
#'   that sums infectious pressure from all strata weighted by contact rates.
#' @returns A stratified StockFlowModel
#' @export
sf_stratify <- function(base, strata_names, flow_types = NULL,
                         cross_strata_flows = list(),
                         param_types = NULL,
                         mixing = NULL) {
  type_sf <- sf_infectious_type()
  type_sn <- sf_snames(type_sf)
  type_fn <- sf_fnames(type_sf)
  type_pn <- sf_pnames(type_sf)

  base_sn <- sf_snames(base)
  base_fn <- sf_fnames(base)
  base_pn <- sf_pnames(base)

  # Default: all stocks map to "pop"
  base_stock_types <- setNames(rep("pop", length(base_sn)), base_sn)

  # Default: all flows map to "disease"
  if (is.null(flow_types)) {
    flow_types <- setNames(rep("disease", length(base_fn)), base_fn)
  }

  # Default: all params map to "rate"
  if (is.null(param_types)) {
    param_types <- setNames(rep("rate", length(base_pn)), base_pn)
  }

  # Type the base model
  tsf_base <- sf_typed(base, type_sf,
                        stock_types = base_stock_types,
                        flow_types = flow_types,
                        param_types = param_types)

  # Add reflexives: each stock needs a reflexive for each non-disease flow type
  # Include "strata" if cross-strata flows exist
  refl_types <- unique(flow_types[flow_types != "disease"])
  if (length(cross_strata_flows) > 0) refl_types <- unique(c(refl_types, "strata"))
  if (length(refl_types) > 0) {
    refl_spec <- setNames(
      lapply(base_sn, function(s) refl_types),
      base_sn
    )
    tsf_base <- sf_add_reflexives(tsf_base, refl_spec)
  }

  # Build strata model
  strata_flows <- list()
  strata_params <- character(0)
  strata_flow_types <- character(0)
  strata_param_types <- character(0)

  # Add identity flow for disease (so product can pair with base disease flows)
  for (sn in strata_names) {
    refl_name <- paste0("id_", sn)
    strata_flows[[refl_name]] <- list(from = sn, to = sn,
                                       rate = parse(text = "0")[[1]])
  }

  # Cross-strata flows
  for (csf in cross_strata_flows) {
    strata_flows[[csf$name]] <- list(
      from = csf$from, to = csf$to,
      rate = if (is.language(csf$rate)) csf$rate else parse(text = "0")[[1]]
    )
    if (!is.null(csf$params)) strata_params <- c(strata_params, csf$params)
    strata_flow_types <- c(strata_flow_types,
                            setNames("strata", csf$name))
  }

  strata_sf <- stock_and_flow(
    stocks = strata_names,
    flows = strata_flows,
    params = unique(strata_params)
  )

  # Type strata
  strata_stock_types <- setNames(rep("pop", length(strata_names)), strata_names)
  all_strata_fn <- sf_fnames(strata_sf)
  all_strata_flow_types <- character(length(all_strata_fn))
  for (i in seq_along(all_strata_fn)) {
    fn <- all_strata_fn[i]
    if (fn %in% names(strata_flow_types)) {
      all_strata_flow_types[i] <- strata_flow_types[fn]
    } else {
      # Identity flows are "disease" type
      all_strata_flow_types[i] <- "disease"
    }
  }
  names(all_strata_flow_types) <- all_strata_fn

  strata_pn <- sf_pnames(strata_sf)
  if (length(strata_pn) > 0) {
    strata_param_types_full <- setNames(rep("rate", length(strata_pn)), strata_pn)
  } else {
    strata_param_types_full <- NULL
  }

  tsf_strata <- sf_typed(strata_sf, type_sf,
                          stock_types = strata_stock_types,
                          flow_types = all_strata_flow_types,
                          param_types = strata_param_types_full)

  # Add reflexives to strata: each stratum needs strata-type reflexive
  if (length(refl_types) > 0) {
    refl_spec_strata <- setNames(
      lapply(strata_names, function(s) refl_types),
      strata_names
    )
    tsf_strata <- sf_add_reflexives(tsf_strata, refl_spec_strata)
  }

  result <- sf_typed_product(tsf_base, tsf_strata, independent_params = TRUE)

  # Apply contact-matrix mixing if specified
  if (!is.null(mixing)) {
    result <- sf_apply_mixing(result, base, strata_names, mixing)
  }

  # Attach stratification metadata for array-based codegen
  attr(result, "strata_info") <- list(
    strata_names = strata_names,
    base_stocks = sf_snames(base),
    base_flows = sf_fnames(base),
    base_params = sf_pnames(base),
    base_sums = names(base@var_exprs)[!grepl("^v_", names(base@var_exprs))],
    flow_types = flow_types,
    cross_strata_flows = cross_strata_flows
  )

  result
}

#' Apply contact-matrix mixing to stratified flow expressions
#'
#' Post-processes a stratified stock-flow model to replace per-stratum
#' infection expressions with contact-matrix-based force-of-infection.
#'
#' @param sf A stratified StockFlowModel
#' @param base The original (unstratified) base model
#' @param strata_names Character vector of strata names
#' @param mixing Named list: base flow name -> contact spec
#' @returns Modified StockFlowModel with mixing expressions
#' @keywords internal
sf_apply_mixing <- function(sf, base, strata_names, mixing) {
  ac <- sf@acset
  sn <- sf_snames(sf)
  fn <- sf_fnames(sf)
  pn <- sf_pnames(sf)
  svn <- sf_svnames(sf)
  sv_stocks <- sf_sum_variable_stocks(sf)

  base_sn <- sf_snames(base)
  base_fn <- sf_fnames(base)
  base_svn <- sf_svnames(base)

  new_var_exprs <- sf@var_exprs
  new_params <- pn
  n_strata <- length(strata_names)

  for (flow_name in names(mixing)) {
    spec <- mixing[[flow_name]]

    # Parse spec
    if (is.character(spec)) {
      contact_prefix <- spec
      mixing_type <- "frequency"
    } else if (is.list(spec)) {
      contact_prefix <- spec$contact
      mixing_type <- if (!is.null(spec$type)) spec$type else "frequency"
    } else {
      cli::cli_abort("Invalid mixing spec for flow '{flow_name}'")
    }

    # Find the base flow's source and target stocks
    fi_base <- match(flow_name, base_fn)
    if (is.na(fi_base)) cli::cli_abort("Flow '{flow_name}' not found in base model")

    from_stock <- NULL
    to_stock <- NULL
    base_ac <- base@acset
    out_rows <- which(acsets::subpart(base_ac, NULL, "ofn") == fi_base)
    if (length(out_rows) > 0)
      from_stock <- base_sn[acsets::subpart(base_ac, out_rows[1], "os")]
    in_rows <- which(acsets::subpart(base_ac, NULL, "ifn") == fi_base)
    if (length(in_rows) > 0)
      to_stock <- base_sn[acsets::subpart(base_ac, in_rows[1], "is")]

    if (is.null(from_stock) || is.null(to_stock))
      cli::cli_abort(
        "Flow '{flow_name}' must have both source and target stocks for mixing")

    # Find sum variable used in base expression (e.g., N)
    base_vn <- paste0("v_", flow_name)
    base_expr_str <- base@var_exprs[[base_vn]]
    norm_var <- NULL
    if (!is.null(base_expr_str) && base_expr_str != "<function>") {
      base_expr <- parse(text = base_expr_str)[[1]]
      expr_vars <- all.vars(base_expr)
      for (sv in base_svn) {
        if (sv %in% expr_vars) { norm_var <- sv; break }
      }
    }

    # Add contact matrix parameters: C_i_j for each pair
    contact_params <- character(0)
    for (si in strata_names) {
      for (sj in strata_names) {
        cp <- paste0(contact_prefix, "_", si, "_", sj)
        contact_params <- c(contact_params, cp)
      }
    }
    new_params <- unique(c(new_params, contact_params))

    # Remove original rate parameter (e.g., beta) — absorbed into contact matrix
    # Keep it as a parameter but it won't appear in mixing flow expressions

    # Build force-of-infection lambda variables for each stratum
    for (si in strata_names) {
      lambda_name <- paste0("lambda_", flow_name, "_", si)
      terms <- character(0)
      for (sj in strata_names) {
        cp <- paste0(contact_prefix, "_", si, "_", sj)
        I_j <- paste0(to_stock, "_", sj)
        if (!is.null(norm_var) && mixing_type == "frequency") {
          N_j <- paste0(norm_var, "_", sj)
          terms <- c(terms, sprintf("%s * %s / %s", cp, I_j, N_j))
        } else {
          terms <- c(terms, sprintf("%s * %s", cp, I_j))
        }
      }
      lambda_expr <- paste(terms, collapse = " + ")
      new_var_exprs[[lambda_name]] <- lambda_expr
    }

    # Rewrite flow rate expressions for each stratum
    for (si in strata_names) {
      strat_flow <- paste0(flow_name, "_id_", si)
      vn <- paste0("v_", strat_flow)

      if (vn %in% names(new_var_exprs)) {
        S_i <- paste0(from_stock, "_", si)
        lambda_name <- paste0("lambda_", flow_name, "_", si)
        new_var_exprs[[vn]] <- sprintf("%s * %s", lambda_name, S_i)
      }
    }
  }

  # Add new parameters to the ACSet
  for (cp in setdiff(new_params, pn)) {
    acsets::add_part(ac, "P", pname = cp)
  }

  # Remove parameters that are no longer referenced in any expression
  all_expr_vars <- character(0)
  for (expr_str in new_var_exprs) {
    if (!is.null(expr_str) && expr_str != "<function>") {
      tryCatch({
        all_expr_vars <- c(all_expr_vars, all.vars(parse(text = expr_str)[[1]]))
      }, error = function(e) NULL)
    }
  }
  all_expr_vars <- unique(all_expr_vars)

  # Find unreferenced params and remove from ACSet
  all_params <- sf_pnames(StockFlowModel(acset = ac))
  for (p in all_params) {
    if (!(p %in% all_expr_vars)) {
      # Remove this parameter from the ACSet
      pi <- match(p, sf_pnames(StockFlowModel(acset = ac)))
      if (!is.na(pi)) acsets::rem_part(ac, "P", pi)
    }
  }

  StockFlowModel(acset = ac, flow_fns = sf@flow_fns, var_exprs = new_var_exprs)
}

#' Create a spatial (multi-patch) model via stratification
#'
#' Convenience wrapper around \code{\link{sf_stratify}} for spatial models.
#' Each patch runs independent disease dynamics, with optional migration
#' between patches and optional cross-patch transmission via a contact matrix.
#'
#' @param base A StockFlowModel (the within-patch disease model)
#' @param patch_names Character vector of patch names
#' @param flow_types Named character: base flow name -> type (default all "disease")
#' @param migration Named list of migration specs. Names are migration flow names.
#'   Each entry is a list with:
#'   \itemize{
#'     \item \code{from}, \code{to}: patch names
#'     \item \code{rate}: quoted R expression (use the source patch name as the
#'       stock variable, e.g., \code{quote(m * patch1)})
#'     \item \code{params}: character vector of parameter names used in rate
#'   }
#' @param contact Optional string: contact matrix parameter prefix. When provided,
#'   transmission between patches uses a contact matrix (frequency-dependent).
#'   Contact parameters are named \code{<prefix>_<to>_<from>}.
#' @returns A stratified StockFlowModel with spatial structure
#' @export
sf_spatial <- function(base, patch_names, flow_types = NULL,
                       migration = list(), contact = NULL) {
  base_fn <- sf_fnames(base)
  if (is.null(flow_types))
    flow_types <- setNames(rep("disease", length(base_fn)), base_fn)

  # Convert migration list to cross_strata_flows format
  cross_strata_flows <- lapply(names(migration), function(nm) {
    spec <- migration[[nm]]
    list(name = nm, from = spec$from, to = spec$to,
         rate = spec$rate,
         params = if (!is.null(spec$params)) spec$params else character(0))
  })

  mixing <- if (!is.null(contact)) {
    # Apply contact matrix to all disease-type flows that have both src and tgt
    disease_flows <- names(flow_types[flow_types == "disease"])
    # Only apply to flows that have a source and target (not births/deaths)
    base_ac <- base@acset
    mix_flows <- character(0)
    for (fn in disease_flows) {
      fi <- match(fn, base_fn)
      has_src <- length(which(acsets::subpart(base_ac, NULL, "ofn") == fi)) > 0
      has_tgt <- length(which(acsets::subpart(base_ac, NULL, "ifn") == fi)) > 0
      if (has_src && has_tgt) mix_flows <- c(mix_flows, fn)
    }
    if (length(mix_flows) > 0) {
      setNames(lapply(mix_flows, function(f) contact), mix_flows)
    } else NULL
  } else NULL

  sf_stratify(base, patch_names,
              flow_types = flow_types,
              cross_strata_flows = cross_strata_flows,
              mixing = mixing)
}