# Stock-flow stratification tests

if (!exists("sf_stratify", mode = "function")) {
  library(algebraicodin)
  library(acsets)
  library(catlab)
}

# ============================================================================
# Basic sf_stratify: SIR × 2 age groups (no cross-strata flows)
# ============================================================================

test_that("sf_stratify basic SIR × 2 strata", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )

  sir2 <- sf_stratify(sir, c("young", "old"),
    flow_types = c(infection = "disease", recovery = "disease"))

  expect_equal(sort(sf_snames(sir2)),
    sort(c("S_young", "S_old", "I_young", "I_old", "R_young", "R_old")))
  expect_equal(length(sf_fnames(sir2)), 4)  # 2 flows × 2 strata
  expect_equal(sort(sf_pnames(sir2)), sort(c("beta", "gamma")))
})

test_that("sf_stratify expressions correct without cross-strata", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )

  sir2 <- sf_stratify(sir, c("young", "old"),
    flow_types = c(infection = "disease", recovery = "disease"))

  # Check infection expressions reference correct stratum stocks
  expr_inf_young <- sir2@var_exprs[["v_infection_id_young"]]
  expect_true(grepl("S_young", expr_inf_young))
  expect_true(grepl("I_young", expr_inf_young))
  expect_true(grepl("N_young", expr_inf_young))

  expr_rec_old <- sir2@var_exprs[["v_recovery_id_old"]]
  expect_true(grepl("I_old", expr_rec_old))
})

# ============================================================================
# sf_stratify with cross-strata flows (aging)
# ============================================================================

test_that("sf_stratify with cross-strata aging flows", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )

  sir_age <- sf_stratify(sir, c("young", "old"),
    flow_types = c(infection = "disease", recovery = "disease"),
    cross_strata_flows = list(
      list(name = "aging", from = "young", to = "old",
           rate = quote(aging_rate * young), params = c("aging_rate"))
    )
  )

  expect_equal(sort(sf_snames(sir_age)),
    sort(c("S_young", "S_old", "I_young", "I_old", "R_young", "R_old")))
  expect_equal(length(sf_fnames(sir_age)), 7)  # 4 disease + 3 aging
  expect_true(all(c("aging_S", "aging_I", "aging_R") %in% sf_fnames(sir_age)))
  expect_equal(sort(sf_pnames(sir_age)), sort(c("beta", "gamma", "aging_rate")))
})

test_that("sf_stratify cross-strata expressions have correct stock refs", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )

  sir_age <- sf_stratify(sir, c("young", "old"),
    flow_types = c(infection = "disease", recovery = "disease"),
    cross_strata_flows = list(
      list(name = "aging", from = "young", to = "old",
           rate = quote(aging_rate * young), params = c("aging_rate"))
    )
  )

  # aging_S should reference S_young (not young_S)
  expr_aging_S <- sir_age@var_exprs[["v_aging_S"]]
  expect_true(grepl("S_young", expr_aging_S))

  # aging_I should reference I_young
  expr_aging_I <- sir_age@var_exprs[["v_aging_I"]]
  expect_true(grepl("I_young", expr_aging_I))
})

# ============================================================================
# sf_to_odin on stratified models
# ============================================================================

test_that("sf_to_odin on stratified SIR", {
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

  code <- sf_to_odin(sir2, type = "ode")
  expect_true(grepl("deriv\\(S_a\\)", code))
  expect_true(grepl("deriv\\(S_b\\)", code))
  expect_true(grepl("N_a", code))
  expect_true(grepl("N_b", code))
})

test_that("sf_to_odin_system on stratified model compiles", {
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

  sir_age <- sf_stratify(sir, c("young", "old"),
    flow_types = c(infection = "disease", recovery = "disease"),
    cross_strata_flows = list(
      list(name = "aging", from = "young", to = "old",
           rate = quote(aging_rate * young), params = c("aging_rate"))
    )
  )

  gen <- sf_to_odin_system(sir_age, type = "ode")
  sys <- dust2::dust_system_create(gen(), list(
    beta = 0.3, gamma = 0.1, aging_rate = 0.02,
    S_young0 = 500, S_old0 = 400,
    I_young0 = 5, I_old0 = 5,
    R_young0 = 0, R_old0 = 0
  ))
  dust2::dust_system_set_state_initial(sys)
  res <- dust2::dust_system_simulate(sys, seq(0, 100, by = 1))
  state <- dust2::dust_unpack_state(sys, res)

  # Population should be conserved (closed system)
  total <- state$S_young[101] + state$I_young[101] + state$R_young[101] +
    state$S_old[101] + state$I_old[101] + state$R_old[101]
  expect_equal(total, 910, tolerance = 1e-4)
})

# ============================================================================
# 3 strata
# ============================================================================

test_that("sf_stratify with 3 strata", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )

  sir3 <- sf_stratify(sir, c("child", "adult", "senior"),
    flow_types = c(infection = "disease", recovery = "disease"),
    cross_strata_flows = list(
      list(name = "aging_ca", from = "child", to = "adult",
           rate = quote(age_rate_1 * child), params = "age_rate_1"),
      list(name = "aging_as", from = "adult", to = "senior",
           rate = quote(age_rate_2 * adult), params = "age_rate_2")
    )
  )

  expect_equal(length(sf_snames(sir3)), 9)  # 3 stocks × 3 strata
  expect_equal(length(sf_fnames(sir3)), 12) # 6 disease + 3 aging_ca + 3 aging_as
  expect_true(all(c("aging_ca_S", "aging_ca_I", "aging_ca_R") %in% sf_fnames(sir3)))
  expect_true(all(c("aging_as_S", "aging_as_I", "aging_as_R") %in% sf_fnames(sir3)))
  expect_equal(sort(sf_pnames(sir3)), sort(c("beta", "gamma", "age_rate_1", "age_rate_2")))
})

# ============================================================================
# TypedStockFlow low-level API
# ============================================================================

test_that("sf_typed creates correct TypedStockFlow", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )

  type_sf <- sf_infectious_type()
  tsf <- sf_typed(sir, type_sf,
    stock_types = c(S = "pop", I = "pop", R = "pop"),
    flow_types = c(infection = "disease", recovery = "disease"),
    param_types = c(beta = "rate", gamma = "rate"))

  expect_s3_class(tsf, "algebraicodin::TypedStockFlow")
  # Maps store indices (1-based) into type model objects
  expect_equal(length(tsf@stock_map), 3)
  expect_equal(length(tsf@flow_map), 2)
})

test_that("sf_add_reflexives adds identity flows", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )

  type_sf <- sf_infectious_type()
  tsf <- sf_typed(sir, type_sf,
    stock_types = c(S = "pop", I = "pop", R = "pop"),
    flow_types = c(infection = "disease", recovery = "disease"),
    param_types = c(beta = "rate", gamma = "rate"))

  tsf_r <- sf_add_reflexives(tsf, list(S = "strata", I = "strata", R = "strata"))
  fn <- sf_fnames(tsf_r@sf)
  expect_true(all(c("refl_S_strata", "refl_I_strata", "refl_R_strata") %in% fn))
  expect_equal(length(fn), 5)  # 2 original + 3 reflexives
})

test_that("sf_typed_product directly produces correct result", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )

  type_sf <- sf_infectious_type()

  tsf_base <- sf_typed(sir, type_sf,
    stock_types = c(S = "pop", I = "pop", R = "pop"),
    flow_types = c(infection = "disease", recovery = "disease"),
    param_types = c(beta = "rate", gamma = "rate"))
  tsf_base <- sf_add_reflexives(tsf_base,
    list(S = "strata", I = "strata", R = "strata"))

  strata_sf <- stock_and_flow(
    stocks = c("young", "old"),
    flows = list(
      id_young = list(from = "young", to = "young", rate = parse(text = "0")[[1]]),
      id_old = list(from = "old", to = "old", rate = parse(text = "0")[[1]]),
      aging = list(from = "young", to = "old", rate = quote(aging_rate))
    ),
    params = c("aging_rate")
  )

  tsf_strata <- sf_typed(strata_sf, type_sf,
    stock_types = c(young = "pop", old = "pop"),
    flow_types = c(id_young = "disease", id_old = "disease", aging = "strata"),
    param_types = c(aging_rate = "rate"))
  tsf_strata <- sf_add_reflexives(tsf_strata,
    list(young = "strata", old = "strata"))

  result <- sf_typed_product(tsf_base, tsf_strata)
  expect_equal(length(sf_snames(result)), 6)
  expect_equal(length(sf_fnames(result)), 7)
  expect_true(all(c("aging_S", "aging_I", "aging_R") %in% sf_fnames(result)))
})

test_that("sf_infectious_type has correct structure", {
  tf <- sf_infectious_type()
  expect_equal(sf_snames(tf), "pop")
  expect_true(all(c("disease", "strata") %in% sf_fnames(tf)))
  expect_equal(sf_pnames(tf), "rate")
})

# ============================================================================
# SEIR stratification
# ============================================================================

test_that("sf_stratify works on SEIR", {
  seir <- stock_and_flow(
    stocks = c("S", "E", "I", "R"),
    flows = list(
      exposure = flow(from = "S", to = "E", rate = beta * S * I / N),
      infection = flow(from = "E", to = "I", rate = sigma * E),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "E", "I", "R")),
    params = c("beta", "sigma", "gamma")
  )

  seir2 <- sf_stratify(seir, c("a", "b"),
    flow_types = c(exposure = "disease", infection = "disease",
                   recovery = "disease"))

  expect_equal(length(sf_snames(seir2)), 8)  # 4 × 2
  expect_equal(length(sf_fnames(seir2)), 6)  # 3 × 2
  expect_equal(sort(sf_pnames(seir2)), sort(c("beta", "sigma", "gamma")))
})

# ============================================================================
# Stochastic and discrete stratified models
# ============================================================================

test_that("sf_to_odin stochastic on stratified model", {
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

  code <- sf_to_odin(sir2, type = "stochastic")
  expect_true(grepl("update\\(S_a\\)", code))
  expect_true(grepl("Binomial", code))
})

# ============================================================================
# Sum variables in stratified models
# ============================================================================

test_that("sf_stratify preserves sum variables per stratum", {
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

  svn <- sf_svnames(sir2)
  expect_true(any(grepl("N_a", svn)))
  expect_true(any(grepl("N_b", svn)))
})

# ============================================================================
# Contact-matrix mixing
# ============================================================================

test_that("sf_stratify with mixing generates contact parameters", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )

  sir_mix <- sf_stratify(sir, c("a", "b"),
    flow_types = c(infection = "disease", recovery = "disease"),
    mixing = list(infection = "C"))

  pn <- sf_pnames(sir_mix)
  # Should have C_a_a, C_a_b, C_b_a, C_b_b + gamma
  expect_true(all(c("C_a_a", "C_a_b", "C_b_a", "C_b_b") %in% pn))
  expect_true("gamma" %in% pn)
  # beta should be removed (absorbed into contact matrix)
  expect_false("beta" %in% pn)
})

test_that("sf_stratify mixing generates lambda expressions", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )

  sir_mix <- sf_stratify(sir, c("a", "b"),
    flow_types = c(infection = "disease", recovery = "disease"),
    mixing = list(infection = "C"))

  # Should have lambda_infection_a and lambda_infection_b
  expect_true("lambda_infection_a" %in% names(sir_mix@var_exprs))
  expect_true("lambda_infection_b" %in% names(sir_mix@var_exprs))

  # Lambda should reference both strata
  lambda_a <- sir_mix@var_exprs[["lambda_infection_a"]]
  expect_true(grepl("C_a_a", lambda_a))
  expect_true(grepl("C_a_b", lambda_a))
  expect_true(grepl("I_a", lambda_a))
  expect_true(grepl("I_b", lambda_a))
  expect_true(grepl("N_a", lambda_a))
  expect_true(grepl("N_b", lambda_a))
})

test_that("sf_stratify mixing rewrites flow expressions", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )

  sir_mix <- sf_stratify(sir, c("a", "b"),
    flow_types = c(infection = "disease", recovery = "disease"),
    mixing = list(infection = "C"))

  # Infection flow should use lambda
  v_inf_a <- sir_mix@var_exprs[["v_infection_id_a"]]
  expect_true(grepl("lambda_infection_a", v_inf_a))
  expect_true(grepl("S_a", v_inf_a))

  # Recovery flow should be unchanged
  v_rec_a <- sir_mix@var_exprs[["v_recovery_id_a"]]
  expect_true(grepl("gamma", v_rec_a))
  expect_true(grepl("I_a", v_rec_a))
})

test_that("sf_to_odin with mixing generates valid odin code", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )

  sir_mix <- sf_stratify(sir, c("a", "b"),
    flow_types = c(infection = "disease", recovery = "disease"),
    mixing = list(infection = "C"))

  code <- sf_to_odin(sir_mix, type = "ode")
  # Lambda should appear before flow expressions
  expect_true(grepl("lambda_infection_a <-", code))
  expect_true(grepl("lambda_infection_b <-", code))
  # Contact params declared
  expect_true(grepl("C_a_a <- parameter", code))
  expect_true(grepl("C_a_b <- parameter", code))
})

test_that("sf_stratify mixing with 3 strata", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )

  sir_mix <- sf_stratify(sir, c("a", "b", "c"),
    flow_types = c(infection = "disease", recovery = "disease"),
    mixing = list(infection = "C"))

  pn <- sf_pnames(sir_mix)
  # Should have 9 contact params (3×3)
  c_params <- pn[grepl("^C_", pn)]
  expect_equal(length(c_params), 9)

  # Lambda for each stratum
  expect_true("lambda_infection_a" %in% names(sir_mix@var_exprs))
  expect_true("lambda_infection_b" %in% names(sir_mix@var_exprs))
  expect_true("lambda_infection_c" %in% names(sir_mix@var_exprs))

  # Lambda_a should sum over all 3 strata
  lambda_a <- sir_mix@var_exprs[["lambda_infection_a"]]
  expect_true(grepl("C_a_a", lambda_a))
  expect_true(grepl("C_a_b", lambda_a))
  expect_true(grepl("C_a_c", lambda_a))
})

test_that("sf_stratify mixing compiles in odin2", {
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

  sir_mix <- sf_stratify(sir, c("young", "old"),
    flow_types = c(infection = "disease", recovery = "disease"),
    mixing = list(infection = "C"))

  gen <- sf_to_odin_system(sir_mix, type = "ode")
  sys <- dust2::dust_system_create(gen(), list(
    gamma = 0.1,
    C_young_young = 0.4, C_young_old = 0.1,
    C_old_young = 0.1, C_old_old = 0.3,
    S_young0 = 500, S_old0 = 400,
    I_young0 = 5, I_old0 = 5,
    R_young0 = 0, R_old0 = 0
  ))
  dust2::dust_system_set_state_initial(sys)
  res <- dust2::dust_system_simulate(sys, seq(0, 100, by = 1))
  state <- dust2::dust_unpack_state(sys, res)

  # Population conserved
  total <- state$S_young[101] + state$I_young[101] + state$R_young[101] +
    state$S_old[101] + state$I_old[101] + state$R_old[101]
  expect_equal(total, 910, tolerance = 1e-4)

  # Both groups should get infected (cross-group mixing)
  expect_true(max(state$I_young) > 50)
  expect_true(max(state$I_old) > 50)
})

test_that("sf_stratify mixing with cross-strata flows", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )

  # Mixing + aging
  sir_mix_age <- sf_stratify(sir, c("young", "old"),
    flow_types = c(infection = "disease", recovery = "disease"),
    cross_strata_flows = list(
      list(name = "aging", from = "young", to = "old",
           rate = quote(aging_rate * young), params = c("aging_rate"))
    ),
    mixing = list(infection = "C")
  )

  # Should have both contact params and aging params
  pn <- sf_pnames(sir_mix_age)
  expect_true("aging_rate" %in% pn)
  expect_true(all(c("C_young_young", "C_young_old") %in% pn))

  # Should have both aging flows and infection flows
  fn <- sf_fnames(sir_mix_age)
  expect_true(all(c("aging_S", "aging_I", "aging_R") %in% fn))
  expect_equal(length(fn), 7)  # 4 disease + 3 aging
})

test_that("sf_stratify density mixing omits normalization", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )

  sir_dens <- sf_stratify(sir, c("a", "b"),
    flow_types = c(infection = "disease", recovery = "disease"),
    mixing = list(infection = list(contact = "C", type = "density")))

  lambda_a <- sir_dens@var_exprs[["lambda_infection_a"]]
  # Density-dependent: no N division
  expect_false(grepl("N_a", lambda_a))
  expect_false(grepl("N_b", lambda_a))
  expect_true(grepl("C_a_a \\* I_a", lambda_a))
})


# ============================================================================
# sf_spatial() convenience wrapper
# ============================================================================

test_that("sf_spatial creates 2-patch model", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )

  sir2p <- sf_spatial(sir, c("urban", "rural"),
    migration = list(
      mig_ur = list(from = "urban", to = "rural",
                    rate = quote(m * urban), params = "m"),
      mig_ru = list(from = "rural", to = "urban",
                    rate = quote(m * rural), params = "m")
    )
  )

  # 6 stocks: S, I, R × 2 patches

  expect_equal(length(sf_snames(sir2p)), 6)
  # 4 disease flows + 6 migration flows (3 compartments × 2 directions)
  expect_equal(length(sf_fnames(sir2p)), 10)
  # Params include migration rate
  expect_true("m" %in% sf_pnames(sir2p))
})

test_that("sf_spatial with contact matrix", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )

  sir2p <- sf_spatial(sir, c("a", "b"),
    migration = list(
      mig_ab = list(from = "a", to = "b",
                    rate = quote(m * a), params = "m"),
      mig_ba = list(from = "b", to = "a",
                    rate = quote(m * b), params = "m")
    ),
    contact = "C"
  )

  # Contact matrix params should exist
  pnames <- sf_pnames(sir2p)
  expect_true("C_a_a" %in% pnames)
  expect_true("C_a_b" %in% pnames)
  expect_true("C_b_a" %in% pnames)
  expect_true("C_b_b" %in% pnames)

  # Lambda expressions should exist
  expect_true("lambda_infection_a" %in% names(sir2p@var_exprs))
  expect_true("lambda_infection_b" %in% names(sir2p@var_exprs))
})

test_that("sf_spatial 3-patch chain", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )

  sir3p <- sf_spatial(sir, c("A", "B", "C"),
    migration = list(
      mig_AB = list(from = "A", to = "B",
                    rate = quote(m * A), params = "m"),
      mig_BA = list(from = "B", to = "A",
                    rate = quote(m * B), params = "m"),
      mig_BC = list(from = "B", to = "C",
                    rate = quote(m * B), params = "m"),
      mig_CB = list(from = "C", to = "B",
                    rate = quote(m * C), params = "m")
    ),
    contact = "C"
  )

  # 9 stocks, 9 contact params, 3 lambda exprs
  expect_equal(length(sf_snames(sir3p)), 9)
  pnames <- sf_pnames(sir3p)
  contact_params <- pnames[grepl("^C_", pnames)]
  expect_equal(length(contact_params), 9)  # 3×3
})

test_that("sf_spatial compiles with odin2", {
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

  sir2p <- sf_spatial(sir, c("a", "b"),
    migration = list(
      mig_ab = list(from = "a", to = "b",
                    rate = quote(m * a), params = "m"),
      mig_ba = list(from = "b", to = "a",
                    rate = quote(m * b), params = "m")
    ),
    contact = "C"
  )

  gen <- sf_to_odin_system(sir2p, type = "ode")
  sys <- dust2::dust_system_create(gen(), list(
    gamma = 0.1, m = 0.02,
    C_a_a = 0.3, C_a_b = 0.05,
    C_b_a = 0.05, C_b_b = 0.25,
    S_a0 = 990, I_a0 = 10, R_a0 = 0,
    S_b0 = 500, I_b0 = 0, R_b0 = 0
  ))
  dust2::dust_system_set_state_initial(sys)
  res <- dust2::dust_system_simulate(sys, seq(0, 100))
  state <- dust2::dust_unpack_state(sys, res)

  # Both patches should see infection
  expect_true(max(state$I_a) > 0)
  expect_true(max(state$I_b) > 0)

  # Population should be approximately conserved
  total_init <- 990 + 10 + 0 + 500 + 0 + 0
  total_end <- state$S_a[101] + state$I_a[101] + state$R_a[101] +
               state$S_b[101] + state$I_b[101] + state$R_b[101]
  expect_equal(total_end, total_init, tolerance = 1)
})

test_that("sf_spatial without contact uses within-patch transmission only", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )

  sir2p <- sf_spatial(sir, c("a", "b"),
    migration = list(
      mig_ab = list(from = "a", to = "b",
                    rate = quote(m * a), params = "m"),
      mig_ba = list(from = "b", to = "a",
                    rate = quote(m * b), params = "m")
    )
  )

  # No lambda expressions — plain within-patch transmission
  lambda_vars <- names(sir2p@var_exprs)[grepl("^lambda_", names(sir2p@var_exprs))]
  expect_equal(length(lambda_vars), 0)

  # beta should still be a parameter (not absorbed into contact matrix)
  expect_true("beta" %in% sf_pnames(sir2p))
})
