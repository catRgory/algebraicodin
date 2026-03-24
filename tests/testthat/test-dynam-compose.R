# Tests for AlgebraicDynamics-style composition
# Tests oapply_dynam (ResourceSharer via UWD) and oapply_machine (Machine via DWD)
#
# When run with test_dir(), the main test file already sources all packages.
# When run standalone with test_file(), we need library() calls.

if (!exists("oapply_dynam", mode = "function")) {
  library(acsets)
  library(catlab)
  library(algebraicodin)
}

# --- Convenience constructors ---

test_that("continuous_resource_sharer creates valid RS", {
  rs <- continuous_resource_sharer(
    nstates = 2L,
    dynamics = function(u, p, t) c(-u[1], u[1]),
    state_names = c("S", "I")
  )
  expect_equal(rs@system_type, "continuous")
  expect_equal(rs@nstates, 2L)
  expect_equal(rs@portmap, 1:2)
  expect_equal(rs@state_names, c("S", "I"))
})

test_that("discrete_resource_sharer creates valid RS", {
  rs <- discrete_resource_sharer(
    nstates = 2L,
    dynamics = function(u, p, t) c(u[1] - 1, u[2] + 1),
    state_names = c("A", "B")
  )
  expect_equal(rs@system_type, "discrete")
  expect_equal(nports(rs), 2L)
})

test_that("delay_resource_sharer creates valid RS", {
  rs <- delay_resource_sharer(
    nstates = 1L,
    dynamics = function(u, h, p, t) {
      u_past <- h(NULL, t - 1)
      -u[1] + u_past[1]
    },
    state_names = c("x")
  )
  expect_equal(rs@system_type, "delay")
})

test_that("continuous_machine creates valid Machine", {
  m <- continuous_machine(
    ninputs = 1L, nstates = 1L, noutputs = 1L,
    dynamics = function(u, x, p, t) x[1] - u[1],
    readout = function(u, p, t) u,
    state_names = c("x")
  )
  expect_equal(m@system_type, "continuous")
  expect_equal(m@ninputs, 1L)
  expect_equal(m@noutputs, 1L)
})

# --- ResourceSharer oapply via UWD ---

test_that("oapply_dynam composes SIR as two resource sharers", {
  # Infection: S,I -> S,I  (both shared)
  infection <- continuous_resource_sharer(
    nstates = 2L,
    dynamics = function(u, p, t) {
      beta <- p$beta
      foi <- beta * u[1] * u[2]
      c(-foi, foi)
    },
    state_names = c("S", "I")
  )

  # Recovery: I -> I,R
  recovery <- continuous_resource_sharer(
    nstates = 2L,
    dynamics = function(u, p, t) {
      gamma <- p$gamma
      c(-gamma * u[1], gamma * u[1])
    },
    state_names = c("I", "R")
  )

  # Build UWD: 2 boxes, 3 junctions (S, I, R)
  # Box 1 (infection): ports -> junctions S, I
  # Box 2 (recovery): ports -> junctions I, R
  d <- uwd(c("S", "I", "R"),
    infection = c("S", "I"),
    recovery = c("I", "R")
  )

  sir <- oapply_dynam(d, list(infection, recovery))
  expect_equal(sir@nstates, 3L)
  expect_equal(sir@system_type, "continuous")
  expect_equal(nports(sir), 3L)

  # Evaluate dynamics at a known state
  u0 <- c(990, 10, 0)
  p <- list(beta = 0.4 / 1000, gamma = 0.2)
  du <- sir@dynamics(u0, p, 0)

  # dS/dt = -beta*S*I = -0.4/1000 * 990 * 10 = -3.96
  # dI/dt = beta*S*I - gamma*I = 3.96 - 2.0 = 1.96
  # dR/dt = gamma*I = 2.0
  expect_equal(du[1], -3.96, tolerance = 1e-10)
  expect_equal(du[2], 1.96, tolerance = 1e-10)
  expect_equal(du[3], 2.0, tolerance = 1e-10)
})

test_that("oapply_dynam composes Lotka-Volterra", {
  # Classic example from AlgebraicDynamics.jl
  # Rabbit growth: dx/dt = alpha*x (shared at junction)
  rabbit_growth <- continuous_resource_sharer(
    nstates = 1L,
    dynamics = function(u, p, t) c(p$alpha * u[1]),
    state_names = c("rabbits")
  )

  # Fox death: dy/dt = -gamma*y
  fox_death <- continuous_resource_sharer(
    nstates = 1L,
    dynamics = function(u, p, t) c(-p$gamma * u[1]),
    state_names = c("foxes")
  )

  # Predation: dx/dt = -beta*x*y, dy/dt = delta*x*y
  predation <- continuous_resource_sharer(
    nstates = 2L,
    dynamics = function(u, p, t) {
      c(-p$beta * u[1] * u[2], p$delta * u[1] * u[2])
    },
    state_names = c("rabbits", "foxes")
  )

  # UWD: 3 boxes, 2 junctions (rabbits, foxes)
  d <- uwd(c("rabbits", "foxes"),
    rabbit_growth = c("rabbits"),
    fox_death = c("foxes"),
    predation = c("rabbits", "foxes")
  )

  lv <- oapply_dynam(d, list(rabbit_growth, fox_death, predation))

  expect_equal(lv@nstates, 2L)
  expect_equal(nports(lv), 2L)

  # Evaluate at (10, 5)
  u0 <- c(10, 5)
  p <- list(alpha = 0.3, beta = 0.02, gamma = 0.1, delta = 0.01)
  du <- lv@dynamics(u0, p, 0)

  # dx/dt = alpha*x - beta*x*y = 0.3*10 - 0.02*10*5 = 3 - 1 = 2
  # dy/dt = -gamma*y + delta*x*y = -0.1*5 + 0.01*10*5 = -0.5 + 0.5 = 0
  expect_equal(du[1], 2.0, tolerance = 1e-10)
  expect_equal(du[2], 0.0, tolerance = 1e-10)
})

test_that("oapply_dynam discrete composition works", {
  # Two discrete sharers that share a state
  growth <- discrete_resource_sharer(
    nstates = 1L,
    dynamics = function(u, p, t) c(u[1] + p$a),
    state_names = c("x")
  )
  decay <- discrete_resource_sharer(
    nstates = 1L,
    dynamics = function(u, p, t) c(u[1] - p$b),
    state_names = c("x")
  )

  # UWD: 2 boxes share 1 junction
  d <- uwd(c("x"),
    growth = c("x"),
    decay = c("x")
  )

  composed <- oapply_dynam(d, list(growth, decay))
  expect_equal(composed@nstates, 1L)
  expect_equal(composed@system_type, "discrete")

  # x_next = x + (delta from growth) + (delta from decay)
  # growth: x -> x + a, delta = a
  # decay: x -> x - b, delta = -b
  # composed: x -> x + a - b
  u1 <- composed@dynamics(c(10), list(a = 3, b = 1), 0)
  expect_equal(u1[1], 12)
})

test_that("oapply_dynam validates inputs", {
  rs <- continuous_resource_sharer(1L, function(u, p, t) c(0))
  d <- uwd(c("x"),
    box1 = c("x"),
    box2 = c("x")
  )
  expect_error(oapply_dynam(d, list(rs)), "Need 2 sharers")
})

# --- Machine oapply via DWD ---

test_that("oapply_machine composes two machines in series", {
  # Machine 1: integrator  dx/dt = input, readout = x
  m1 <- continuous_machine(
    ninputs = 1L, nstates = 1L, noutputs = 1L,
    dynamics = function(u, x, p, t) x,  # dx/dt = input
    readout = function(u, p, t) u,       # output = state
    state_names = c("x1")
  )

  # Machine 2: gain  dx/dt = -x + input, readout = x
  m2 <- continuous_machine(
    ninputs = 1L, nstates = 1L, noutputs = 1L,
    dynamics = function(u, x, p, t) c(-u[1] + x[1]),
    readout = function(u, p, t) u,
    state_names = c("x2")
  )

  # Wire: ext_input -> m1 -> m2 -> ext_output
  d <- dwd(
    ninputs = 1L, noutputs = 1L,
    boxes = list(c(1, 1), c(1, 1)),
    wires = list(
      c(0, 1, 1, 1),   # external input -> m1 input 1
      c(1, 1, 2, 1),   # m1 output 1 -> m2 input 1
      c(2, 1, -1, 1)   # m2 output 1 -> external output
    )
  )

  composed <- oapply_machine(d, list(m1, m2))
  expect_equal(composed@nstates, 2L)
  expect_equal(composed@ninputs, 1L)
  expect_equal(composed@noutputs, 1L)

  # At state (3, 5) with input 1:
  # m1: dx1/dt = 1 (input)
  # m2: dx2/dt = -5 + 3 = -2 (readout of m1 is 3)
  du <- composed@dynamics(c(3, 5), c(1), list(), 0)
  expect_equal(du[1], 1.0)
  expect_equal(du[2], -2.0)

  # Readout: m2 output = x2 = 5
  y <- composed@readout(c(3, 5), list(), 0)
  expect_equal(y[1], 5.0)
})

test_that("oapply_machine handles parallel machines", {
  m1 <- continuous_machine(
    ninputs = 1L, nstates = 1L, noutputs = 1L,
    dynamics = function(u, x, p, t) c(x[1] - u[1]),
    readout = function(u, p, t) c(2 * u[1]),
    state_names = c("a")
  )
  m2 <- continuous_machine(
    ninputs = 1L, nstates = 1L, noutputs = 1L,
    dynamics = function(u, x, p, t) c(x[1] + u[1]),
    readout = function(u, p, t) c(3 * u[1]),
    state_names = c("b")
  )

  # Both get external inputs, both produce external outputs
  d <- dwd(
    ninputs = 2L, noutputs = 2L,
    boxes = list(c(1, 1), c(1, 1)),
    wires = list(
      c(0, 1, 1, 1),   # ext input 1 -> m1
      c(0, 2, 2, 1),   # ext input 2 -> m2
      c(1, 1, -1, 1),  # m1 -> ext output 1
      c(2, 1, -1, 2)   # m2 -> ext output 2
    )
  )

  composed <- oapply_machine(d, list(m1, m2))
  expect_equal(composed@nstates, 2L)
  expect_equal(composed@ninputs, 2L)
  expect_equal(composed@noutputs, 2L)

  du <- composed@dynamics(c(1, 2), c(10, 20), list(), 0)
  # m1: x[1]-u[1] = 10-1 = 9
  # m2: x[1]+u[1] = 20+2 = 22
  expect_equal(du[1], 9.0)
  expect_equal(du[2], 22.0)

  y <- composed@readout(c(1, 2), list(), 0)
  # m1: 2*1 = 2, m2: 3*2 = 6
  expect_equal(y, c(2, 6))
})

# --- simulate_rs ---

test_that("simulate_rs solves SIR composed via oapply_dynam", {
  skip_if_not_installed("deSolve")

  infection <- continuous_resource_sharer(
    nstates = 2L,
    dynamics = function(u, p, t) {
      foi <- p$beta * u[1] * u[2]
      c(-foi, foi)
    },
    state_names = c("S", "I")
  )
  recovery <- continuous_resource_sharer(
    nstates = 2L,
    dynamics = function(u, p, t) {
      c(-p$gamma * u[1], p$gamma * u[1])
    },
    state_names = c("I", "R")
  )

  d <- uwd(c("S", "I", "R"),
    infection = c("S", "I"),
    recovery = c("I", "R")
  )

  sir <- oapply_dynam(d, list(infection, recovery))

  result <- simulate_rs(sir,
    initial = c(S = 990, I = 10, R = 0),
    times = seq(0, 100, by = 1),
    params = list(beta = 0.4 / 1000, gamma = 0.2)
  )

  expect_equal(ncol(result), 4L)  # time, S, I, R
  expect_equal(nrow(result), 101L)
  # Conservation: S + I + R = 1000
  totals <- result$S + result$I + result$R
  expect_equal(totals, rep(1000, 101), tolerance = 1e-6)
  # Peak should occur
  expect_true(max(result$I) > 10)
  # R should increase monotonically
  expect_true(all(diff(result$R) >= -1e-6))
})

test_that("simulate_rs_discrete works for composed system", {
  growth <- discrete_resource_sharer(
    nstates = 1L,
    dynamics = function(u, p, t) c(u[1] * 1.1),
    state_names = c("pop")
  )
  harvest <- discrete_resource_sharer(
    nstates = 1L,
    dynamics = function(u, p, t) c(u[1] - p$h),
    state_names = c("pop")
  )

  d <- uwd(c("pop"),
    growth = c("pop"),
    harvest = c("pop")
  )

  composed <- oapply_dynam(d, list(growth, harvest))
  result <- simulate_rs_discrete(composed,
    initial = c(pop = 100),
    nsteps = 5,
    params = list(h = 5)
  )

  expect_equal(nrow(result), 6L)  # t=0..5
  expect_equal(result$time, 0:5)
  # Step 1: delta_growth = 100*1.1 - 100 = 10, delta_harvest = 100-5 - 100 = -5
  # next = 100 + 10 + (-5) = 105
  expect_equal(result$pop[2], 105, tolerance = 1e-10)
})

test_that("simulate_machine works with external input", {
  skip_if_not_installed("deSolve")

  # Simple integrator: dx/dt = input(t)
  m <- continuous_machine(
    ninputs = 1L, nstates = 1L, noutputs = 1L,
    dynamics = function(u, x, p, t) x,
    readout = function(u, p, t) u,
    state_names = c("x")
  )

  result <- simulate_machine(m,
    initial = c(x = 0),
    times = seq(0, 10, by = 0.1),
    input_fn = function(t) c(1),  # constant input
    params = list()
  )

  # x(t) = t for constant input of 1
  expect_equal(result$x[nrow(result)], 10.0, tolerance = 0.01)
})

# --- DWD constructor ---

test_that("dwd creates valid wiring diagram", {
  d <- dwd(2L, 1L,
    boxes = list(c(1, 1), c(2, 1)),
    wires = list(
      c(0, 1, 1, 1),
      c(0, 2, 2, 1),
      c(1, 1, 2, 2),
      c(2, 1, -1, 1)
    )
  )
  expect_equal(d$ninputs, 2L)
  expect_equal(d$noutputs, 1L)
  expect_equal(length(d$boxes), 2L)
  expect_equal(length(d$wires), 4L)
  expect_equal(d$wires[[1]]$source_box, 0L)
})

# --- Petri net -> ResourceSharer round-trip ---

test_that("petri_to_continuous matches direct ResourceSharer composition", {
  skip_if_not_installed("deSolve")

  # Build SIR Petri net directly
  sir_pn <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )

  # Convert to ResourceSharer
  rs_petri <- petri_to_continuous(sir_pn)

  # Direct composition (using same param names as petri)
  infection <- continuous_resource_sharer(
    nstates = 2L,
    dynamics = function(u, p, t) {
      foi <- p$inf * u[1] * u[2]
      c(-foi, foi)
    },
    state_names = c("S", "I")
  )
  recovery <- continuous_resource_sharer(
    nstates = 2L,
    dynamics = function(u, p, t) {
      c(-p$rec * u[1], p$rec * u[1])
    },
    state_names = c("I", "R")
  )
  d <- uwd(c("S", "I", "R"),
    infection = c("S", "I"),
    recovery = c("I", "R")
  )
  rs_composed <- oapply_dynam(d, list(infection, recovery))

  # Both should give same dynamics at same state
  u0 <- c(990, 10, 0)
  p <- list(inf = 0.0004, rec = 0.2)

  du_petri <- rs_petri@dynamics(u0, p, 0)
  du_composed <- rs_composed@dynamics(u0, p, 0)

  expect_equal(du_petri, du_composed, tolerance = 1e-10)
})

# --- Delay resource sharer ---

test_that("delay_resource_sharer composes via oapply_dynam", {
  # Simple delay system: dx/dt = -x(t) + x(t-tau)
  delay_rs <- delay_resource_sharer(
    nstates = 1L,
    dynamics = function(u, h, p, t) {
      u_past <- h(p, t - p$tau)
      c(-u[1] + u_past[1])
    },
    state_names = c("x")
  )

  # Compose with itself (no actual sharing, separate junctions)
  d <- uwd(c("x1", "x2"),
    box1 = c("x1"),
    box2 = c("x2")
  )

  composed <- oapply_dynam(d, list(delay_rs, delay_rs))
  expect_equal(composed@nstates, 2L)
  expect_equal(composed@system_type, "delay")
})
