# Tests for array-based codegen and model examples

if (!exists("sf_to_odin_array", mode = "function")) {
  library(algebraicodin)
  library(acsets)
  library(catlab)
}

# ============================================================================
# sf_to_odin_array() — basic array codegen
# ============================================================================

test_that("sf_to_odin_array generates array code for basic SIR", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )
  sir2 <- sf_stratify(sir, c("a", "b"),
    flow_types = c(infection = "disease", recovery = "disease"))
  code <- sf_to_odin_array(sir2, type = "ode")
  expect_true(grepl("n_strata <- 2", code, fixed = TRUE))
  expect_true(grepl("dim(S, I, R) <- n_strata", code, fixed = TRUE))
  expect_true(grepl("deriv(S[]) <-", code, fixed = TRUE))
  expect_true(grepl("v_infection[] <-", code, fixed = TRUE))
})

test_that("sf_to_odin_array errors without strata metadata", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )
  expect_error(sf_to_odin_array(sir), "strata metadata")
})

test_that("sf_to_odin_array is much smaller than scalar", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )
  sir10 <- sf_stratify(sir, paste0("g", 1:10),
    flow_types = c(infection = "disease", recovery = "disease"))
  scalar_code <- sf_to_odin(sir10, type = "ode")
  array_code <- sf_to_odin_array(sir10, type = "ode")
  expect_true(nchar(array_code) < nchar(scalar_code) / 3)
})

# ============================================================================
# Array vs scalar cross-validation
# ============================================================================

test_that("array ODE matches scalar ODE (no cross-strata)", {
  skip_if_not_installed("odin2")
  skip_if_not_installed("dust2")

  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )
  sir3 <- sf_stratify(sir, c("a", "b", "c"),
    flow_types = c(infection = "disease", recovery = "disease"))

  gen_arr <- sf_to_odin_array_system(sir3, type = "ode")
  sys_arr <- dust2::dust_system_create(gen_arr(), list(
    beta = 0.3, gamma = 0.1,
    S0 = c(990, 500, 200), I0 = c(10, 5, 2), R0 = c(0, 0, 0)
  ))
  dust2::dust_system_set_state_initial(sys_arr)
  res_arr <- dust2::dust_system_simulate(sys_arr, seq(0, 50))
  st_arr <- dust2::dust_unpack_state(sys_arr, res_arr)

  gen_sc <- sf_to_odin_system(sir3, type = "ode")
  sys_sc <- dust2::dust_system_create(gen_sc(), list(
    beta = 0.3, gamma = 0.1,
    S_a0 = 990, S_b0 = 500, S_c0 = 200,
    I_a0 = 10, I_b0 = 5, I_c0 = 2,
    R_a0 = 0, R_b0 = 0, R_c0 = 0
  ))
  dust2::dust_system_set_state_initial(sys_sc)
  res_sc <- dust2::dust_system_simulate(sys_sc, seq(0, 50))
  st_sc <- dust2::dust_unpack_state(sys_sc, res_sc)

  expect_equal(st_arr$S[1,], st_sc$S_a, tolerance = 1e-10)
  expect_equal(st_arr$S[2,], st_sc$S_b, tolerance = 1e-10)
  expect_equal(st_arr$I[3,], st_sc$I_c, tolerance = 1e-10)
})

test_that("array ODE matches scalar ODE (with cross-strata aging)", {
  skip_if_not_installed("odin2")
  skip_if_not_installed("dust2")

  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )
  sir3 <- sf_stratify(sir, c("y", "m", "o"),
    flow_types = c(infection = "disease", recovery = "disease"),
    cross_strata_flows = list(
      list(name = "age1", from = "y", to = "m",
           rate = quote(mu1 * y), params = "mu1"),
      list(name = "age2", from = "m", to = "o",
           rate = quote(mu2 * m), params = "mu2")
    ))

  gen_arr <- sf_to_odin_array_system(sir3, type = "ode")
  sys_arr <- dust2::dust_system_create(gen_arr(), list(
    beta = 0.3, gamma = 0.1, mu1 = 0.02, mu2 = 0.01,
    S0 = c(800, 400, 200), I0 = c(10, 5, 2), R0 = c(0, 0, 0)
  ))
  dust2::dust_system_set_state_initial(sys_arr)
  res_arr <- dust2::dust_system_simulate(sys_arr, seq(0, 50))
  st_arr <- dust2::dust_unpack_state(sys_arr, res_arr)

  gen_sc <- sf_to_odin_system(sir3, type = "ode")
  sys_sc <- dust2::dust_system_create(gen_sc(), list(
    beta = 0.3, gamma = 0.1, mu1 = 0.02, mu2 = 0.01,
    S_y0 = 800, S_m0 = 400, S_o0 = 200,
    I_y0 = 10, I_m0 = 5, I_o0 = 2,
    R_y0 = 0, R_m0 = 0, R_o0 = 0
  ))
  dust2::dust_system_set_state_initial(sys_sc)
  res_sc <- dust2::dust_system_simulate(sys_sc, seq(0, 50))
  st_sc <- dust2::dust_unpack_state(sys_sc, res_sc)

  expect_equal(st_arr$S[1,], st_sc$S_y, tolerance = 1e-8)
  expect_equal(st_arr$S[2,], st_sc$S_m, tolerance = 1e-8)
  expect_equal(st_arr$S[3,], st_sc$S_o, tolerance = 1e-8)
  expect_equal(st_arr$I[1,], st_sc$I_y, tolerance = 1e-8)
  expect_equal(st_arr$I[2,], st_sc$I_m, tolerance = 1e-8)
})

test_that("array stochastic compiles and runs", {
  skip_if_not_installed("odin2")
  skip_if_not_installed("dust2")

  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )
  sir2 <- sf_stratify(sir, c("a", "b"),
    flow_types = c(infection = "disease", recovery = "disease"))

  gen <- sf_to_odin_array_system(sir2, type = "stochastic")
  sys <- dust2::dust_system_create(gen(), list(
    beta = 0.3, gamma = 0.1,
    S0 = c(990, 500), I0 = c(10, 5), R0 = c(0, 0)
  ), n_particles = 1, dt = 0.25)
  dust2::dust_system_set_state_initial(sys)
  res <- dust2::dust_system_simulate(sys, seq(0, 50))
  st <- dust2::dust_unpack_state(sys, res)

  expect_true(max(st$I) > 0)
})

# ============================================================================
# Model examples: Ross-Macdonald vector-borne
# ============================================================================

test_that("Ross-Macdonald vector-borne model", {
  skip_if_not_installed("odin2")
  skip_if_not_installed("dust2")

  # Human: Sh → Ih → Rh (infected by mosquitoes)
  # Mosquito: Sv → Iv (infected by humans, no recovery)
  ross_mac <- stock_and_flow(
    stocks = c("Sh", "Ih", "Rh", "Sv", "Iv"),
    flows = list(
      human_inf = flow(from = "Sh", to = "Ih",
        rate = a * b * Iv / Nh * Sh),
      human_rec = flow(from = "Ih", to = "Rh", rate = r * Ih),
      mosq_inf = flow(from = "Sv", to = "Iv",
        rate = a * c * Ih / Nh * Sv),
      mosq_birth = flow(from = NULL, to = "Sv", rate = mu * Nm),
      mosq_death_s = flow(from = "Sv", to = NULL, rate = mu * Sv),
      mosq_death_i = flow(from = "Iv", to = NULL, rate = mu * Iv)
    ),
    sums = list(Nh = c("Sh", "Ih", "Rh"), Nm = c("Sv", "Iv")),
    params = c("a", "b", "c", "r", "mu")
  )

  code <- sf_to_odin(ross_mac, type = "ode")
  gen <- odin2::odin(code)
  sys <- dust2::dust_system_create(gen(), list(
    a = 0.3, b = 0.5, c = 0.5, r = 0.01, mu = 0.1,
    Sh0 = 990, Ih0 = 10, Rh0 = 0, Sv0 = 5000, Iv0 = 100
  ))
  dust2::dust_system_set_state_initial(sys)
  res <- dust2::dust_system_simulate(sys, seq(0, 365))
  st <- dust2::dust_unpack_state(sys, res)

  expect_true(max(st$Ih) > 10)  # infection spreads
  # Total human pop conserved
  total_h <- st$Sh[366] + st$Ih[366] + st$Rh[366]
  expect_equal(total_h, 1000, tolerance = 1)
})

# ============================================================================
# Model examples: Multi-strain SIR
# ============================================================================

test_that("Multi-strain SIR via stratification", {
  skip_if_not_installed("odin2")
  skip_if_not_installed("dust2")

  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )

  # Two strains: shared S pool, independent I and R
  sir_2strain <- sf_stratify(sir, c("strain1", "strain2"),
    flow_types = c(infection = "disease", recovery = "disease"))

  gen <- sf_to_odin_system(sir_2strain, type = "ode")
  sys <- dust2::dust_system_create(gen(), list(
    beta = 0.3, gamma = 0.1,
    S_strain10 = 990, S_strain20 = 990,
    I_strain10 = 10, I_strain20 = 5,
    R_strain10 = 0, R_strain20 = 0
  ))
  dust2::dust_system_set_state_initial(sys)
  res <- dust2::dust_system_simulate(sys, seq(0, 200))
  st <- dust2::dust_unpack_state(sys, res)

  # Both strains should cause epidemics
  expect_true(max(st$I_strain1) > 10)
  expect_true(max(st$I_strain2) > 5)
})

# ============================================================================
# Model examples: SEIR with demographics
# ============================================================================

test_that("SEIR with demographics model", {
  skip_if_not_installed("odin2")
  skip_if_not_installed("dust2")

  seir_demo <- stock_and_flow(
    stocks = c("S", "E", "I", "R"),
    flows = list(
      exposure = flow(from = "S", to = "E", rate = beta * S * I / N),
      progression = flow(from = "E", to = "I", rate = sigma * E),
      recovery = flow(from = "I", to = "R", rate = gamma * I),
      birth = flow(from = NULL, to = "S", rate = mu * N),
      death_S = flow(from = "S", to = NULL, rate = mu * S),
      death_E = flow(from = "E", to = NULL, rate = mu * E),
      death_I = flow(from = "I", to = NULL, rate = mu * I),
      death_R = flow(from = "R", to = NULL, rate = mu * R)
    ),
    sums = list(N = c("S", "E", "I", "R")),
    params = c("beta", "sigma", "gamma", "mu")
  )

  code <- sf_to_odin(seir_demo, type = "ode")
  gen <- odin2::odin(code)
  sys <- dust2::dust_system_create(gen(), list(
    beta = 0.5, sigma = 0.2, gamma = 0.1, mu = 0.01,
    S0 = 990, E0 = 0, I0 = 10, R0 = 0
  ))
  dust2::dust_system_set_state_initial(sys)
  res <- dust2::dust_system_simulate(sys, seq(0, 365))
  st <- dust2::dust_unpack_state(sys, res)

  # Epidemic should peak then settle to endemic
  expect_true(max(st$I) > 10)
  # Population stays near 1000 (births ≈ deaths)
  total <- st$S[366] + st$E[366] + st$I[366] + st$R[366]
  expect_equal(total, 1000, tolerance = 50)
})

# ============================================================================
# Model examples: SIS (no immunity)
# ============================================================================

test_that("SIS model (no immunity)", {
  skip_if_not_installed("odin2")
  skip_if_not_installed("dust2")

  sis <- stock_and_flow(
    stocks = c("S", "I"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "S", rate = gamma * I)
    ),
    sums = list(N = c("S", "I")),
    params = c("beta", "gamma")
  )

  gen <- sf_to_odin_system(sis, type = "ode")
  sys <- dust2::dust_system_create(gen(), list(
    beta = 0.5, gamma = 0.1,
    S0 = 990, I0 = 10
  ))
  dust2::dust_system_set_state_initial(sys)
  res <- dust2::dust_system_simulate(sys, seq(0, 200))
  st <- dust2::dust_unpack_state(sys, res)

  # SIS reaches endemic equilibrium: I* = N(1 - gamma/beta)
  I_star <- 1000 * (1 - 0.1 / 0.5)
  expect_equal(st$I[201], I_star, tolerance = 1)
})

# ============================================================================
# 10-strata array compilation
# ============================================================================

test_that("10 strata array model compiles", {
  skip_if_not_installed("odin2")
  skip_if_not_installed("dust2")

  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )

  sir10 <- sf_stratify(sir, paste0("g", 1:10),
    flow_types = c(infection = "disease", recovery = "disease"))

  gen <- sf_to_odin_array_system(sir10, type = "ode")
  sys <- dust2::dust_system_create(gen(), list(
    beta = 0.3, gamma = 0.1,
    S0 = rep(100, 10), I0 = c(5, rep(0, 9)), R0 = rep(0, 10)
  ))
  dust2::dust_system_set_state_initial(sys)
  res <- dust2::dust_system_simulate(sys, seq(0, 100))
  st <- dust2::dust_unpack_state(sys, res)

  # First group (seeded) should have epidemic
  expect_true(max(st$I[1,]) > 5)
  # Total population conserved
  total <- sum(st$S[, 101]) + sum(st$I[, 101]) + sum(st$R[, 101])
  expect_equal(total, 1005, tolerance = 1)
})
