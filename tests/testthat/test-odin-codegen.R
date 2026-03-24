# Tests for odin2 code generation from ResourceSharers and Machines
# These test the expression-based construction and code generation features.

if (!exists("continuous_resource_sharer", mode = "function")) {
  library(algebraicodin)
  library(catlab)
}

# === ResourceSharer expression-based construction ===

test_that("expression-based RS constructor works", {
  rs <- continuous_resource_sharer(
    dynamics_expr = list(
      S = quote(-beta * S * I),
      I = quote(beta * S * I - gamma * I)
    ),
    params = c("beta", "gamma")
  )
  expect_equal(rs@nstates, 2L)
  expect_equal(rs@state_names, c("S", "I"))
  expect_equal(rs@params, c("beta", "gamma"))
  expect_equal(length(rs@dynamics_expr), 2L)
})

test_that("expression-based RS auto-builds working closure", {
  rs <- continuous_resource_sharer(
    dynamics_expr = list(
      S = quote(-beta * S * I),
      I = quote(beta * S * I - gamma * I)
    ),
    params = c("beta", "gamma")
  )
  du <- rs@dynamics(c(990, 10), list(beta = 0.001, gamma = 0.1), 0)
  expect_equal(du, c(-9.9, 8.9), tolerance = 1e-10)
})

test_that("expression-based RS infers state_names from expr names", {
  rs <- continuous_resource_sharer(
    dynamics_expr = list(x = quote(-a * x), y = quote(a * x)),
    params = c("a")
  )
  expect_equal(rs@state_names, c("x", "y"))
  expect_equal(rs@nstates, 2L)
})

test_that("RS still works with closure-only (backward compat)", {
  rs <- continuous_resource_sharer(
    nstates = 2,
    dynamics = function(u, p, t) c(-p$beta * u[1] * u[2], p$beta * u[1] * u[2]),
    state_names = c("S", "I")
  )
  du <- rs@dynamics(c(990, 10), list(beta = 0.001), 0)
  expect_equal(du, c(-9.9, 9.9), tolerance = 1e-10)
  expect_equal(length(rs@dynamics_expr), 0L)
})

# === rs_to_odin code generation ===

test_that("rs_to_odin generates valid ODE code", {
  rs <- continuous_resource_sharer(
    dynamics_expr = list(
      S = quote(-beta * S * I),
      I = quote(beta * S * I - gamma * I),
      R = quote(gamma * I)
    ),
    params = c("beta", "gamma")
  )
  code <- rs_to_odin(rs, initial = c(S = 990, I = 10, R = 0), type = "ode")
  expect_true(grepl("deriv\\(S\\) <- -beta \\* S \\* I", code))
  expect_true(grepl("deriv\\(I\\)", code))
  expect_true(grepl("deriv\\(R\\) <- gamma \\* I", code))
  expect_true(grepl("beta <- parameter\\(\\)", code))
  expect_true(grepl("S0 <- parameter\\(990\\)", code))
  expect_true(grepl("initial\\(S\\) <- S0", code))
})

test_that("rs_to_odin generates discrete code", {
  rs <- continuous_resource_sharer(
    dynamics_expr = list(x = quote(-a * x)),
    params = c("a")
  )
  code <- rs_to_odin(rs, initial = c(x = 1), type = "discrete")
  expect_true(grepl("update\\(x\\)", code))
  expect_true(grepl("dt <- parameter\\(1\\)", code))
})

test_that("rs_to_odin generates stochastic code", {
  rs <- continuous_resource_sharer(
    dynamics_expr = list(x = quote(-a * x)),
    params = c("a")
  )
  code <- rs_to_odin(rs, initial = c(x = 1), type = "stochastic")
  expect_true(grepl("update\\(x\\)", code))
  expect_true(grepl("dt <- parameter\\(0.01\\)", code))
})

test_that("rs_to_odin errors without expressions", {
  rs <- continuous_resource_sharer(
    nstates = 1,
    dynamics = function(u, p, t) -u[1],
    state_names = c("x")
  )
  expect_error(rs_to_odin(rs, initial = c(x = 1)), "no symbolic expressions")
})

test_that("rs_to_odin auto-detects params from RS", {
  rs <- continuous_resource_sharer(
    dynamics_expr = list(x = quote(-a * x)),
    params = c("a")
  )
  code <- rs_to_odin(rs, initial = c(x = 1))
  expect_true(grepl("a <- parameter\\(\\)", code))
})

# === Composed ResourceSharers with expressions ===

test_that("oapply_dynam preserves and merges expressions", {
  infection <- continuous_resource_sharer(
    dynamics_expr = list(S = quote(-beta * S * I), I = quote(beta * S * I)),
    params = c("beta")
  )
  recovery <- continuous_resource_sharer(
    dynamics_expr = list(I = quote(-gamma * I), R = quote(gamma * I)),
    params = c("gamma")
  )
  d <- uwd(c("S", "I", "R"), c("S", "I"), c("I", "R"))
  sir <- oapply_dynam(d, list(infection, recovery))

  expect_equal(sir@nstates, 3L)
  expect_equal(sir@state_names, c("S", "I", "R"))
  expect_equal(sort(sir@params), c("beta", "gamma"))
  expect_equal(length(sir@dynamics_expr), 3L)

  # S should be: -beta * S * I (from infection only)
  expect_equal(deparse(sir@dynamics_expr$S), "-beta * S * I")
  # R should be: gamma * I (from recovery only)
  expect_equal(deparse(sir@dynamics_expr$R), "gamma * I")
  # I should be a sum: beta * S * I + -gamma * I (from both)
  I_expr <- sir@dynamics_expr$I
  expect_true(is.call(I_expr))
})

test_that("composed RS generates compilable odin code", {
  infection <- continuous_resource_sharer(
    dynamics_expr = list(S = quote(-beta * S * I), I = quote(beta * S * I)),
    params = c("beta")
  )
  recovery <- continuous_resource_sharer(
    dynamics_expr = list(I = quote(-gamma * I), R = quote(gamma * I)),
    params = c("gamma")
  )
  d <- uwd(c("S", "I", "R"), c("S", "I"), c("I", "R"))
  sir <- oapply_dynam(d, list(infection, recovery))

  code <- rs_to_odin(sir, initial = c(S = 990, I = 10, R = 0), type = "ode")
  expect_true(grepl("deriv\\(S\\)", code))
  expect_true(grepl("deriv\\(I\\)", code))
  expect_true(grepl("deriv\\(R\\)", code))
  expect_true(grepl("beta <- parameter\\(\\)", code))
  expect_true(grepl("gamma <- parameter\\(\\)", code))
})

test_that("composed RS expression dynamics matches closure dynamics", {
  infection <- continuous_resource_sharer(
    dynamics_expr = list(S = quote(-beta * S * I), I = quote(beta * S * I)),
    params = c("beta")
  )
  recovery <- continuous_resource_sharer(
    dynamics_expr = list(I = quote(-gamma * I), R = quote(gamma * I)),
    params = c("gamma")
  )
  d <- uwd(c("S", "I", "R"), c("S", "I"), c("I", "R"))
  sir <- oapply_dynam(d, list(infection, recovery))

  du <- sir@dynamics(c(990, 10, 0), list(beta = 0.001, gamma = 0.1), 0)
  expect_equal(du, c(-9.9, 8.9, 1.0), tolerance = 1e-10)
})

test_that("oapply_dynam returns empty expressions if any sharer lacks them", {
  rs_with <- continuous_resource_sharer(
    dynamics_expr = list(S = quote(-beta * S * I), I = quote(beta * S * I)),
    params = c("beta")
  )
  rs_without <- continuous_resource_sharer(
    nstates = 2,
    dynamics = function(u, p, t) c(-p$gamma * u[1], p$gamma * u[1]),
    state_names = c("I", "R")
  )
  d <- uwd(c("S", "I", "R"), c("S", "I"), c("I", "R"))
  sir <- oapply_dynam(d, list(rs_with, rs_without))

  # Should have no expressions since one sharer lacks them
  expect_equal(length(sir@dynamics_expr), 0L)
})

# === Machine expression-based construction ===

test_that("expression-based Machine constructor works", {
  m <- continuous_machine(
    ninputs = 1,
    dynamics_expr = list(x = quote(input_1)),
    readout_expr = list(quote(x)),
    state_names = c("x")
  )
  expect_equal(m@nstates, 1L)
  expect_equal(m@ninputs, 1L)
  expect_equal(m@noutputs, 1L)
  expect_equal(length(m@dynamics_expr), 1L)
  expect_equal(length(m@readout_expr), 1L)
})

test_that("expression-based Machine auto-builds working closures", {
  m <- continuous_machine(
    ninputs = 1,
    dynamics_expr = list(y = quote(-alpha * y + input_1)),
    readout_expr = list(quote(y)),
    state_names = c("y"),
    params = c("alpha")
  )
  # dynamics: dy/dt = -alpha*y + input
  du <- m@dynamics(c(5), c(10), list(alpha = 2), 0)
  expect_equal(du, c(-2 * 5 + 10), tolerance = 1e-10)

  # readout
  y <- m@readout(c(3.5), list(), 0)
  expect_equal(y, 3.5)
})

# === machine_to_odin code generation ===

test_that("machine_to_odin generates valid ODE code", {
  m <- continuous_machine(
    ninputs = 1,
    dynamics_expr = list(y = quote(-alpha * y + input_1)),
    readout_expr = list(quote(y)),
    state_names = c("y"),
    params = c("alpha")
  )
  code <- machine_to_odin(m, initial = c(y = 0), type = "ode")
  expect_true(grepl("ext_input_1 <- parameter\\(\\)", code))
  expect_true(grepl("alpha <- parameter\\(\\)", code))
  expect_true(grepl("deriv\\(y\\)", code))
  expect_true(grepl("initial\\(y\\) <- y0", code))
})

test_that("machine_to_odin errors without expressions", {
  m <- continuous_machine(
    ninputs = 1, nstates = 1, noutputs = 1,
    dynamics = function(u, x, p, t) c(x[1]),
    readout = function(u, p, t) u,
    state_names = c("x")
  )
  expect_error(machine_to_odin(m, initial = c(x = 0)), "no symbolic expressions")
})

test_that("machine_to_odin handles input_params", {
  m <- continuous_machine(
    ninputs = 1,
    dynamics_expr = list(x = quote(input_1)),
    readout_expr = list(quote(x)),
    state_names = c("x")
  )
  code <- machine_to_odin(m, initial = c(x = 0),
                           input_params = list(ext_input_1 = 1.0))
  expect_true(grepl("ext_input_1 <- parameter\\(1\\)", code))
})

# === Composed Machines with expressions ===

test_that("oapply_machine propagates expressions through wiring", {
  integrator <- continuous_machine(
    ninputs = 1,
    dynamics_expr = list(x = quote(input_1)),
    readout_expr = list(quote(x)),
    state_names = c("x")
  )
  filter <- continuous_machine(
    ninputs = 1,
    dynamics_expr = list(y = quote(-alpha * y + input_1)),
    readout_expr = list(quote(y)),
    state_names = c("y"),
    params = c("alpha")
  )
  dwg <- dwd(
    ninputs = 1, noutputs = 1,
    boxes = list(int = c(1, 1), filt = c(1, 1)),
    wires = list(
      c(0, 1, 1, 1),
      c(1, 1, 2, 1),
      c(2, 1, -1, 1)
    )
  )
  comp <- oapply_machine(dwg, list(integrator, filter))

  expect_equal(comp@nstates, 2L)
  expect_equal(comp@state_names, c("x", "y"))
  expect_equal(length(comp@dynamics_expr), 2L)

  # x's dynamics: input_1 -> ext_input_1 (external input wired to box 1)
  expect_equal(deparse(comp@dynamics_expr$x), "ext_input_1")

  # y's dynamics: -alpha * y + input_1 -> -alpha * y + x
  # (integrator's readout "x" substituted for filter's input_1)
  y_expr <- deparse(comp@dynamics_expr$y)
  expect_true(grepl("alpha", y_expr))
  expect_true(grepl("x", y_expr))
})

test_that("composed Machine generates odin code", {
  integrator <- continuous_machine(
    ninputs = 1,
    dynamics_expr = list(x = quote(input_1)),
    readout_expr = list(quote(x)),
    state_names = c("x")
  )
  filter <- continuous_machine(
    ninputs = 1,
    dynamics_expr = list(y = quote(-alpha * y + input_1)),
    readout_expr = list(quote(y)),
    state_names = c("y"),
    params = c("alpha")
  )
  dwg <- dwd(
    ninputs = 1, noutputs = 1,
    boxes = list(int = c(1, 1), filt = c(1, 1)),
    wires = list(
      c(0, 1, 1, 1),
      c(1, 1, 2, 1),
      c(2, 1, -1, 1)
    )
  )
  comp <- oapply_machine(dwg, list(integrator, filter))
  code <- machine_to_odin(comp, initial = c(x = 0, y = 0), type = "ode")

  expect_true(grepl("deriv\\(x\\) <- ext_input_1", code))
  expect_true(grepl("deriv\\(y\\)", code))
  expect_true(grepl("alpha <- parameter\\(\\)", code))
  expect_true(grepl("ext_input_1 <- parameter\\(\\)", code))
})

# === Expression helpers ===

test_that("rename_expr_vars substitutes variable names", {
  expr <- quote(a * x + b * y)
  renamed <- algebraicodin:::rename_expr_vars(expr, c(x = "S", y = "I"))
  expect_equal(deparse(renamed), "a * S + b * I")
})

test_that("sum_exprs sums list of expressions", {
  e1 <- quote(a * x)
  e2 <- quote(-b * y)
  result <- algebraicodin:::sum_exprs(list(e1, e2))
  expect_equal(deparse(result), "a * x + -b * y")
})

test_that("sum_exprs handles single expression", {
  result <- algebraicodin:::sum_exprs(list(quote(a * x)))
  expect_equal(deparse(result), "a * x")
})

test_that("sum_exprs handles empty list", {
  result <- algebraicodin:::sum_exprs(list())
  expect_equal(deparse(result), "0")
})

# === Discrete RS with expressions ===

test_that("discrete RS with expressions works", {
  rs <- discrete_resource_sharer(
    dynamics_expr = list(
      x = quote(x + a * x * (1 - x))
    ),
    params = c("a")
  )
  expect_equal(rs@system_type, "discrete")
  u_next <- rs@dynamics(c(0.1), list(a = 2.5), 0)
  expect_equal(u_next, 0.1 + 2.5 * 0.1 * 0.9, tolerance = 1e-10)
})

# === Lotka-Volterra composition example ===

test_that("Lotka-Volterra composition generates correct code", {
  rabbit_growth <- continuous_resource_sharer(
    dynamics_expr = list(x = quote(alpha * x)),
    params = c("alpha")
  )
  fox_death <- continuous_resource_sharer(
    dynamics_expr = list(y = quote(-gamma * y)),
    params = c("gamma")
  )
  predation <- continuous_resource_sharer(
    dynamics_expr = list(
      x = quote(-beta * x * y),
      y = quote(delta * x * y)
    ),
    params = c("beta", "delta")
  )

  d <- uwd(c("x", "y"), c("x"), c("y"), c("x", "y"))
  lv <- oapply_dynam(d, list(rabbit_growth, fox_death, predation))

  expect_equal(lv@nstates, 2L)
  expect_equal(sort(lv@params), c("alpha", "beta", "delta", "gamma"))

  code <- rs_to_odin(lv, initial = c(x = 10, y = 5))
  expect_true(grepl("deriv\\(x\\)", code))
  expect_true(grepl("deriv\\(y\\)", code))
  expect_true(grepl("alpha", code))
  expect_true(grepl("beta", code))

  # Verify dynamics match expected Lotka-Volterra
  du <- lv@dynamics(c(10, 5), list(alpha = 1.1, beta = 0.4, gamma = 0.4, delta = 0.1), 0)
  # dx/dt = alpha*x - beta*x*y = 1.1*10 - 0.4*10*5 = 11 - 20 = -9
  # dy/dt = -gamma*y + delta*x*y = -0.4*5 + 0.1*10*5 = -2 + 5 = 3
  expect_equal(du, c(-9, 3), tolerance = 1e-10)
})

# === odin2 compilation test ===

test_that("composed RS code compiles and runs in odin2", {
  skip_if_not_installed("odin2")
  skip_if_not_installed("dust2")

  infection <- continuous_resource_sharer(
    dynamics_expr = list(S = quote(-beta * S * I), I = quote(beta * S * I)),
    params = c("beta")
  )
  recovery <- continuous_resource_sharer(
    dynamics_expr = list(I = quote(-gamma * I), R = quote(gamma * I)),
    params = c("gamma")
  )
  d <- uwd(c("S", "I", "R"), c("S", "I"), c("I", "R"))
  sir <- oapply_dynam(d, list(infection, recovery))

  code <- rs_to_odin(sir, initial = c(S = 990, I = 10, R = 0))
  gen <- suppressMessages(odin2::odin(code))

  sys <- dust2::dust_system_create(gen, list(beta = 0.001, gamma = 0.1), n_particles = 1)
  dust2::dust_system_set_state_initial(sys)
  t <- seq(0, 100, by = 1)
  sol <- dust2::dust_system_simulate(sys, t)

  # Compare to deSolve
  sim_ds <- simulate_rs(sir, initial = c(S = 990, I = 10, R = 0),
                         times = t, params = list(beta = 0.001, gamma = 0.1))
  odin_final <- sol[, ncol(sol)]
  desolve_final <- c(tail(sim_ds$S, 1), tail(sim_ds$I, 1), tail(sim_ds$R, 1))
  expect_equal(odin_final, desolve_final, tolerance = 1e-3)
})
