# Tests for stock-and-flow models

if (!exists("stock_and_flow", mode = "function")) {
  library(acsets)
  library(catlab)
  library(algebraicodin)
}

# --- Schema and basic construction ---

test_that("StockFlow schemas are valid", {
  expect_true(!is.null(SchStockFlow0))
  expect_true(!is.null(SchStockFlow))
})

test_that("StockFlowType creates empty ACSet", {
  sf <- StockFlowType()
  expect_equal(acsets::nparts(sf, "S"), 0L)
  expect_equal(acsets::nparts(sf, "F"), 0L)
})

# --- SIR stock-and-flow model ---

test_that("stock_and_flow builds SIR model correctly", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      inf = list(from = "S", to = "I", rate = quote(beta * S * I / N)),
      rec = list(from = "I", to = "R", rate = quote(gamma * I))
    ),
    params = c("beta", "gamma"),
    sums = list(N = c("S", "I", "R"))
  )

  expect_true(S7::S7_inherits(sir, StockFlowModel))
  expect_equal(sf_snames(sir), c("S", "I", "R"))
  expect_equal(sf_fnames(sir), c("inf", "rec"))
  expect_equal(sf_pnames(sir), c("beta", "gamma"))
  expect_equal(sf_svnames(sir), "N")
  expect_equal(sf_vnames(sir), c("v_inf", "v_rec"))
})

test_that("sf_transition_matrices correct for SIR", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      inf = list(from = "S", to = "I", rate = quote(beta * S * I / N)),
      rec = list(from = "I", to = "R", rate = quote(gamma * I))
    ),
    params = c("beta", "gamma"),
    sums = list(N = c("S", "I", "R"))
  )

  tm <- sf_transition_matrices(sir)

  # Inflow: inf->I, rec->R
  expect_equal(tm$inflow["inf", "S"], 0L)
  expect_equal(tm$inflow["inf", "I"], 1L)
  expect_equal(tm$inflow["rec", "R"], 1L)

  # Outflow: S->inf, I->rec
  expect_equal(tm$outflow["inf", "S"], 1L)
  expect_equal(tm$outflow["rec", "I"], 1L)

  # Stoichiometry: cols=stocks, rows=stocks
  expect_equal(tm$stoichiometry["S", "inf"], -1L)
  expect_equal(tm$stoichiometry["I", "inf"], 1L)
  expect_equal(tm$stoichiometry["I", "rec"], -1L)
  expect_equal(tm$stoichiometry["R", "rec"], 1L)
})

# --- Vectorfield ---

test_that("sf_vectorfield produces correct SIR dynamics", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      inf = list(from = "S", to = "I", rate = quote(beta * S * I / N)),
      rec = list(from = "I", to = "R", rate = quote(gamma * I))
    ),
    params = c("beta", "gamma"),
    sums = list(N = c("S", "I", "R"))
  )

  vf <- sf_vectorfield(sir)
  state <- c(S = 990, I = 10, R = 0)
  parms <- list(beta = 0.4, gamma = 0.2)
  result <- vf(0, state, parms)
  du <- result[[1]]

  # N = 1000
  # inf_rate = 0.4 * 990 * 10 / 1000 = 3.96
  # rec_rate = 0.2 * 10 = 2.0
  # dS = -inf_rate = -3.96
  # dI = inf_rate - rec_rate = 1.96
  # dR = rec_rate = 2.0
  expect_equal(du[["S"]], -3.96, tolerance = 1e-10)
  expect_equal(du[["I"]], 1.96, tolerance = 1e-10)
  expect_equal(du[["R"]], 2.0, tolerance = 1e-10)
})

# --- Simulation ---

test_that("simulate_sf runs SIR model", {
  skip_if_not_installed("deSolve")

  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      inf = list(from = "S", to = "I", rate = quote(beta * S * I / N)),
      rec = list(from = "I", to = "R", rate = quote(gamma * I))
    ),
    params = c("beta", "gamma"),
    sums = list(N = c("S", "I", "R"))
  )

  result <- simulate_sf(sir,
    initial = c(S = 990, I = 10, R = 0),
    times = seq(0, 100, by = 1),
    params = list(beta = 0.4, gamma = 0.2)
  )

  expect_equal(nrow(result), 101L)
  expect_equal(ncol(result), 4L)  # time, S, I, R
  # Conservation
  totals <- result$S + result$I + result$R
  expect_equal(totals, rep(1000, 101), tolerance = 1e-5)
  # Peak I
  expect_true(max(result$I) > 10)
  # R monotonically increases
  expect_true(all(diff(result$R) >= -1e-6))
})

# --- odin code generation ---

test_that("sf_to_odin generates valid ODE code", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      inf = list(from = "S", to = "I", rate = quote(beta * S * I / N)),
      rec = list(from = "I", to = "R", rate = quote(gamma * I))
    ),
    params = c("beta", "gamma"),
    sums = list(N = c("S", "I", "R"))
  )

  code <- sf_to_odin(sir, type = "ode",
                      initial = c(S = 990, I = 10, R = 0))

  expect_true(grepl("beta <- parameter()", code, fixed = TRUE))
  expect_true(grepl("gamma <- parameter()", code, fixed = TRUE))
  expect_true(grepl("N <- S + I + R", code, fixed = TRUE))
  expect_true(grepl("deriv(S)", code, fixed = TRUE))
  expect_true(grepl("deriv(I)", code, fixed = TRUE))
  expect_true(grepl("deriv(R)", code, fixed = TRUE))
  expect_true(grepl("initial(S)", code, fixed = TRUE))
})

test_that("sf_to_odin generates stochastic code", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      inf = list(from = "S", to = "I", rate = quote(beta * S * I / N)),
      rec = list(from = "I", to = "R", rate = quote(gamma * I))
    ),
    params = c("beta", "gamma"),
    sums = list(N = c("S", "I", "R"))
  )

  code <- sf_to_odin(sir, type = "stochastic")
  expect_true(grepl("Binomial", code, fixed = TRUE))
  expect_true(grepl("update(S)", code, fixed = TRUE))
  # dt is provided by odin2 automatically, not declared in generated code
  expect_true(grepl("* dt)", code, fixed = TRUE))
})

test_that("sf_to_odin generates discrete code", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      inf = list(from = "S", to = "I", rate = quote(beta * S * I / N)),
      rec = list(from = "I", to = "R", rate = quote(gamma * I))
    ),
    params = c("beta", "gamma"),
    sums = list(N = c("S", "I", "R"))
  )

  code <- sf_to_odin(sir, type = "discrete")
  expect_true(grepl("update(S)", code, fixed = TRUE))
  expect_true(grepl("* dt", code, fixed = TRUE))
})

# --- Conversion to/from Petri net ---

test_that("sf_to_petri converts SIR correctly", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      inf = list(from = "S", to = "I", rate = quote(beta * S * I / N)),
      rec = list(from = "I", to = "R", rate = quote(gamma * I))
    ),
    params = c("beta", "gamma"),
    sums = list(N = c("S", "I", "R"))
  )

  pn <- sf_to_petri(sir)
  expect_equal(acsets::nparts(pn, "S"), 3L)
  expect_equal(acsets::nparts(pn, "T"), 2L)

  sn <- vapply(seq_len(acsets::nparts(pn, "S")),
    function(i) acsets::subpart(pn, i, "sname"), character(1))
  expect_equal(sn, c("S", "I", "R"))
  tn <- vapply(seq_len(acsets::nparts(pn, "T")),
    function(i) acsets::subpart(pn, i, "tname"), character(1))
  expect_equal(tn, c("inf", "rec"))
})

# --- Conversion to ResourceSharer ---

test_that("sf_to_resource_sharer produces correct dynamics", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      inf = list(from = "S", to = "I", rate = quote(beta * S * I / N)),
      rec = list(from = "I", to = "R", rate = quote(gamma * I))
    ),
    params = c("beta", "gamma"),
    sums = list(N = c("S", "I", "R"))
  )

  rs <- sf_to_resource_sharer(sir)
  expect_equal(rs@nstates, 3L)
  expect_equal(rs@system_type, "continuous")

  u0 <- c(S = 990, I = 10, R = 0)
  p <- list(beta = 0.4, gamma = 0.2)
  du <- rs@dynamics(u0, p, 0)

  expect_equal(du[1], -3.96, tolerance = 1e-10)
  expect_equal(du[2], 1.96, tolerance = 1e-10)
  expect_equal(du[3], 2.0, tolerance = 1e-10)
})

# --- SEIR model ---

test_that("stock_and_flow builds SEIR model", {
  seir <- stock_and_flow(
    stocks = c("S", "E", "I", "R"),
    flows = list(
      inf = list(from = "S", to = "E", rate = quote(beta * S * I / N)),
      prog = list(from = "E", to = "I", rate = quote(sigma * E)),
      rec = list(from = "I", to = "R", rate = quote(gamma * I))
    ),
    params = c("beta", "sigma", "gamma"),
    sums = list(N = c("S", "E", "I", "R"))
  )

  expect_equal(sf_snames(seir), c("S", "E", "I", "R"))
  expect_equal(length(sf_fnames(seir)), 3L)

  vf <- sf_vectorfield(seir)
  state <- c(S = 990, E = 0, I = 10, R = 0)
  parms <- list(beta = 0.4, sigma = 0.5, gamma = 0.2)
  du <- vf(0, state, parms)[[1]]

  # inf_rate = 0.4 * 990 * 10 / 1000 = 3.96
  # prog_rate = 0.5 * 0 = 0
  # rec_rate = 0.2 * 10 = 2.0
  expect_equal(du[["S"]], -3.96, tolerance = 1e-10)
  expect_equal(du[["E"]], 3.96, tolerance = 1e-10)
  expect_equal(du[["I"]], -2.0, tolerance = 1e-10)
  expect_equal(du[["R"]], 2.0, tolerance = 1e-10)
})

# --- External flows (births/deaths) ---

test_that("external inflows/outflows work", {
  sir_demo <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      birth = list(from = NULL, to = "S", rate = quote(mu * N)),
      inf = list(from = "S", to = "I", rate = quote(beta * S * I / N)),
      rec = list(from = "I", to = "R", rate = quote(gamma * I)),
      death_S = list(from = "S", to = NULL, rate = quote(mu * S)),
      death_I = list(from = "I", to = NULL, rate = quote(mu * I)),
      death_R = list(from = "R", to = NULL, rate = quote(mu * R))
    ),
    params = c("beta", "gamma", "mu"),
    sums = list(N = c("S", "I", "R"))
  )

  # Check transitions
  tm <- sf_transition_matrices(sir_demo)
  # birth: no outflow source, inflow to S
  expect_equal(tm$inflow["birth", "S"], 1L)
  expect_equal(sum(tm$outflow["birth", ]), 0L)

  # death_S: outflow from S, no inflow
  expect_equal(tm$outflow["death_S", "S"], 1L)
  expect_equal(sum(tm$inflow["death_S", ]), 0L)

  # Verify vectorfield at equilibrium-like state
  vf <- sf_vectorfield(sir_demo)
  state <- c(S = 900, I = 50, R = 50)
  parms <- list(beta = 0.4, gamma = 0.2, mu = 0.01)
  du <- vf(0, state, parms)[[1]]

  # birth = mu * N = 0.01 * 1000 = 10
  # death_S = mu * S = 0.01 * 900 = 9
  # inf = beta * S * I / N = 0.4 * 900 * 50 / 1000 = 18
  # rec = gamma * I = 0.2 * 50 = 10
  # death_I = mu * I = 0.5
  # death_R = mu * R = 0.5
  # dS = birth - inf - death_S = 10 - 18 - 9 = -17
  expect_equal(du[["S"]], -17.0, tolerance = 1e-10)
})

test_that("stock_and_flow validates malformed flow specifications", {
  expect_error(
    stock_and_flow(
      stocks = c("S", "I"),
      flows = list(bad = list(from = "S", to = "I"))
    ),
    "must define a `rate`"
  )

  expect_error(
    stock_and_flow(
      stocks = c("S", "I"),
      flows = list(bad = list(from = "X", to = "I", rate = quote(beta * S)))
    ),
    "unknown source stock"
  )

  expect_error(
    stock_and_flow(
      stocks = c("S", "I"),
      flows = list(bad = list(from = NULL, to = NULL, rate = quote(beta)))
    ),
    "must specify at least one endpoint"
  )
})

# --- Vectorfield comparison: stock-flow vs Petri net ---

test_that("sf vectorfield matches Petri net vectorfield for simple SIR", {
  skip_if_not_installed("deSolve")

  # Stock-flow SIR with mass-action
  sir_sf <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      inf = list(from = "S", to = "I", rate = quote(beta * S * I)),
      rec = list(from = "I", to = "R", rate = quote(gamma * I))
    ),
    params = c("beta", "gamma")
  )

  # Petri net SIR (mass-action by construction)
  sir_pn <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )

  # Compare dynamics at same state
  state <- c(S = 990, I = 10, R = 0)
  parms_sf <- list(beta = 0.0004, gamma = 0.2)
  parms_pn <- list(inf = 0.0004, rec = 0.2)

  vf_sf <- sf_vectorfield(sir_sf)
  du_sf <- vf_sf(0, state, parms_sf)[[1]]

  vf_pn <- vectorfield(sir_pn)
  du_pn <- vf_pn(0, state, parms_pn)[[1]]

  expect_equal(unname(du_sf), unname(du_pn), tolerance = 1e-10)
})

# --- Accessor edge cases ---

test_that("accessors work on empty model", {
  sf <- StockFlowModel()
  expect_equal(sf_snames(sf), character(0))
  expect_equal(sf_fnames(sf), character(0))
  expect_equal(sf_pnames(sf), character(0))
  expect_equal(sf_svnames(sf), character(0))
})

test_that("multiple sum variables work", {
  model <- stock_and_flow(
    stocks = c("S1", "S2", "I1", "I2"),
    flows = list(
      f1 = list(from = "S1", to = "I1", rate = quote(beta * S1 * I1 / N1)),
      f2 = list(from = "S2", to = "I2", rate = quote(beta * S2 * I2 / N2))
    ),
    params = c("beta"),
    sums = list(
      N1 = c("S1", "I1"),
      N2 = c("S2", "I2")
    )
  )

  expect_equal(sf_svnames(model), c("N1", "N2"))

  vf <- sf_vectorfield(model)
  state <- c(S1 = 90, S2 = 80, I1 = 10, I2 = 20)
  du <- vf(0, state, list(beta = 0.5))[[1]]

  # f1 = 0.5 * 90 * 10 / 100 = 4.5
  # f2 = 0.5 * 80 * 20 / 100 = 8.0
  expect_equal(du[["S1"]], -4.5, tolerance = 1e-10)
  expect_equal(du[["I1"]], 4.5, tolerance = 1e-10)
  expect_equal(du[["S2"]], -8.0, tolerance = 1e-10)
  expect_equal(du[["I2"]], 8.0, tolerance = 1e-10)
})
