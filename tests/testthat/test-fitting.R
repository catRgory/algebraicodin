# Tests for data fitting infrastructure

if (!exists("observe", mode = "function")) {
  library(algebraicodin)
  library(acsets)
  library(catlab)
}

# ============================================================================
# observe() specification
# ============================================================================

test_that("observe() creates spec with transition", {
  spec <- observe("Poisson", transition = "inf")
  expect_s3_class(spec, "observe_spec")
  expect_equal(spec$dist, "Poisson")
  expect_equal(spec$transition, "inf")
  expect_null(spec$state)
  expect_equal(spec$zero_every, 1)
})

test_that("observe() creates spec with state", {
  spec <- observe("Normal", state = "I", sd = "sigma")
  expect_equal(spec$dist, "Normal")
  expect_equal(spec$state, "I")
  expect_equal(spec$sd, "sigma")
})

test_that("observe() creates spec with flow", {
  spec <- observe("Poisson", flow = "infection")
  expect_equal(spec$flow, "infection")
})

test_that("observe() creates spec with expr", {
  spec <- observe("Poisson", expr = "incidence + noise")
  expect_equal(spec$expr, "incidence + noise")
})

test_that("observe() errors when no source given", {
  expect_error(observe("Poisson"),
               "Exactly one of")
})

test_that("observe() errors when multiple sources given", {
  expect_error(observe("Poisson", transition = "inf", state = "I"),
               "Exactly one of")
})

test_that("observe() supports noise parameter", {
  spec <- observe("Poisson", transition = "inf", noise = "exp_noise")
  expect_equal(spec$noise, "exp_noise")
})

test_that("observe() supports zero_every", {
  spec <- observe("Poisson", transition = "inf", zero_every = 7)
  expect_equal(spec$zero_every, 7)
})

# ============================================================================
# generate_compare_code() — internal helper
# ============================================================================

test_that("generate_compare_code handles simple string spec", {
  code <- algebraicodin:::generate_compare_code(
    list(cases = "Poisson"), "petri")
  code_str <- paste(code, collapse = "\n")
  expect_true(grepl("cases <- data()", code_str, fixed = TRUE))
  expect_true(grepl("cases ~ Poisson(cases)", code_str, fixed = TRUE))
})

test_that("generate_compare_code handles transition incidence", {
  code <- algebraicodin:::generate_compare_code(
    list(cases = observe("Poisson", transition = "inf")), "petri")
  code_str <- paste(code, collapse = "\n")
  expect_true(grepl("initial(incidence_inf, zero_every = 1) <- 0",
                     code_str, fixed = TRUE))
  expect_true(grepl("update(incidence_inf) <- incidence_inf + n_inf",
                     code_str, fixed = TRUE))
  expect_true(grepl("cases ~ Poisson(incidence_inf)", code_str, fixed = TRUE))
})

test_that("generate_compare_code handles flow incidence", {
  code <- algebraicodin:::generate_compare_code(
    list(cases = observe("Poisson", flow = "infection")), "stockflow")
  code_str <- paste(code, collapse = "\n")
  expect_true(grepl("incidence_infection", code_str))
  expect_true(grepl("cases ~ Poisson(incidence_infection)",
                     code_str, fixed = TRUE))
})

test_that("generate_compare_code handles state observation", {
  code <- algebraicodin:::generate_compare_code(
    list(prevalence = observe("Normal", state = "I", sd = "sigma")), "petri")
  code_str <- paste(code, collapse = "\n")
  expect_true(grepl("sigma <- parameter()", code_str, fixed = TRUE))
  expect_true(grepl("prevalence ~ Normal(I, sigma)", code_str, fixed = TRUE))
})

test_that("generate_compare_code handles noise parameter", {
  code <- algebraicodin:::generate_compare_code(
    list(cases = observe("Poisson", transition = "inf",
                         noise = "exp_noise")), "petri")
  code_str <- paste(code, collapse = "\n")
  expect_true(grepl("exp_noise <- parameter(1e6)", code_str, fixed = TRUE))
  expect_true(grepl("noise_cases <- incidence_inf + 1 / exp_noise",
                     code_str, fixed = TRUE))
  expect_true(grepl("cases ~ Poisson(noise_cases)", code_str, fixed = TRUE))
})

test_that("generate_compare_code handles custom expr", {
  code <- algebraicodin:::generate_compare_code(
    list(deaths = observe("Poisson", expr = "mu * I")), "petri")
  code_str <- paste(code, collapse = "\n")
  expect_true(grepl("deaths ~ Poisson(mu * I)", code_str, fixed = TRUE))
})

test_that("generate_compare_code returns empty for NULL", {
  expect_equal(length(algebraicodin:::generate_compare_code(NULL)), 0)
})

# ============================================================================
# pn_to_odin() with enhanced compare
# ============================================================================

test_that("pn_to_odin with observe() transition spec", {
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  code <- pn_to_odin(sir, type = "stochastic",
    compare = list(cases = observe("Poisson", transition = "inf")))
  expect_true(grepl("initial(incidence_inf, zero_every = 1) <- 0",
                     code, fixed = TRUE))
  expect_true(grepl("update(incidence_inf) <- incidence_inf + n_inf",
                     code, fixed = TRUE))
  expect_true(grepl("cases <- data()", code, fixed = TRUE))
  expect_true(grepl("cases ~ Poisson(incidence_inf)", code, fixed = TRUE))
})

test_that("pn_to_odin with observe() and noise", {
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  code <- pn_to_odin(sir, type = "stochastic",
    compare = list(
      cases = observe("Poisson", transition = "inf", noise = "exp_noise")
    ))
  expect_true(grepl("exp_noise <- parameter(1e6)", code, fixed = TRUE))
  expect_true(grepl("noise_cases <- incidence_inf + 1 / exp_noise",
                     code, fixed = TRUE))
  expect_true(grepl("cases ~ Poisson(noise_cases)", code, fixed = TRUE))
})

test_that("pn_to_odin backward compat: simple string compare", {
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  code <- pn_to_odin(sir, type = "stochastic",
    compare = list(cases = "Poisson"))
  expect_true(grepl("cases <- data()", code, fixed = TRUE))
  expect_true(grepl("cases ~ Poisson(cases)", code, fixed = TRUE))
})

# ============================================================================
# sf_to_odin() with compare
# ============================================================================

test_that("sf_to_odin with compare parameter", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )

  code <- sf_to_odin(sir, type = "stochastic",
    compare = list(cases = observe("Poisson", flow = "infection")))
  expect_true(grepl("cases <- data()", code, fixed = TRUE))
  expect_true(grepl("incidence_infection", code))
  expect_true(grepl("cases ~ Poisson(incidence_infection)", code, fixed = TRUE))
})

test_that("sf_to_odin with prevalence observation", {
  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )

  code <- sf_to_odin(sir, type = "ode",
    compare = list(prevalence = observe("Normal", state = "I", sd = "sigma")))
  expect_true(grepl("sigma <- parameter()", code, fixed = TRUE))
  expect_true(grepl("prevalence <- data()", code, fixed = TRUE))
  expect_true(grepl("prevalence ~ Normal(I, sigma)", code, fixed = TRUE))
})

# ============================================================================
# prior helpers
# ============================================================================

test_that("prior_exp returns valid spec", {
  p <- prior_exp(mean = 0.3)
  expect_true(is.function(p$log_density))
  expect_equal(p$domain, c(0, Inf))
  # Check density at x = 0.3
  expected <- dexp(0.3, rate = 1/0.3, log = TRUE)
  expect_equal(p$log_density(0.3), expected)
  # Negative values should be -Inf
  expect_equal(p$log_density(-1), -Inf)
})

test_that("prior_uniform returns valid spec", {
  p <- prior_uniform(min = 0, max = 1)
  expect_equal(p$domain, c(0, 1))
  expect_equal(p$log_density(0.5), dunif(0.5, 0, 1, log = TRUE))
  expect_equal(p$log_density(1.5), -Inf)
})

test_that("prior_normal returns valid spec", {
  p <- prior_normal(mean = 0, sd = 1)
  expect_equal(p$domain, c(-Inf, Inf))
  expect_equal(p$log_density(0), dnorm(0, 0, 1, log = TRUE))
})

test_that("prior_lognormal returns valid spec", {
  p <- prior_lognormal(meanlog = 0, sdlog = 1)
  expect_equal(p$domain, c(0, Inf))
  expect_equal(p$log_density(1), dlnorm(1, 0, 1, log = TRUE))
  expect_equal(p$log_density(-1), -Inf)
})

# ============================================================================
# build_prior()
# ============================================================================

test_that("build_prior creates monty_model", {
  skip_if_not_installed("monty")
  p <- build_prior(
    beta = prior_exp(mean = 0.3),
    gamma = prior_exp(mean = 0.1)
  )
  expect_s3_class(p, "monty_model")
})

test_that("build_prior errors without names", {
  skip_if_not_installed("monty")
  expect_error(build_prior(prior_exp()), "must be named")
})

test_that("build_prior density matches manual", {
  skip_if_not_installed("monty")
  p <- build_prior(
    beta = prior_exp(mean = 0.3),
    gamma = prior_uniform(min = 0.05, max = 0.5)
  )
  expected <- dexp(0.2, rate = 1/0.3, log = TRUE) +
    dunif(0.1, 0.05, 0.5, log = TRUE)
  actual <- monty::monty_model_density(p, c(0.2, 0.1))
  expect_equal(actual, expected)
})

# ============================================================================
# create_filter() — compile + filter creation
# ============================================================================

test_that("create_filter works with stochastic Petri net", {
  skip_if_not_installed("odin2")
  skip_if_not_installed("dust2")

  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )

  # Synthetic data
  data <- data.frame(time = c(5, 10, 15), cases = c(10, 30, 20))

  filter <- create_filter(sir, data,
    compare = list(cases = observe("Poisson", transition = "inf")),
    type = "stochastic", n_particles = 10)
  expect_true(inherits(filter, "dust_likelihood"))
})

test_that("create_filter computes likelihood", {
  skip_if_not_installed("odin2")
  skip_if_not_installed("dust2")

  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )

  data <- data.frame(time = c(5, 10, 15, 20),
                     cases = c(5, 15, 25, 10))

  filter <- create_filter(sir, data,
    compare = list(cases = observe("Poisson", transition = "inf")),
    type = "stochastic", n_particles = 50)

  pars <- list(inf = 0.001, rec = 0.1, S0 = 990, I0 = 10, R0 = 0)
  ll <- dust2::dust_likelihood_run(filter, pars)
  expect_true(is.finite(ll))
})

test_that("create_filter works with stock-flow model", {
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

  data <- data.frame(time = c(5, 10, 15), cases = c(10, 30, 20))

  filter <- create_filter(sir, data,
    compare = list(cases = observe("Poisson", flow = "infection")),
    type = "stochastic", n_particles = 10)
  expect_true(inherits(filter, "dust_likelihood"))
})

# ============================================================================
# fit_mcmc() — end-to-end MCMC
# ============================================================================

test_that("fit_mcmc runs end-to-end with Petri net", {
  skip_if_not_installed("odin2")
  skip_if_not_installed("dust2")
  skip_if_not_installed("monty")

  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )

  # Use noise-added Poisson compare for robustness
  data <- data.frame(time = seq(5, 50, by = 5),
                     cases = c(3, 8, 15, 25, 30, 20, 12, 6, 3, 1))

  prior <- build_prior(
    inf = prior_uniform(min = 1e-7, max = 1e-4),
    rec = prior_uniform(min = 0.01, max = 0.5)
  )

  samples <- fit_mcmc(
    model = sir, data = data,
    pars = c("inf", "rec"),
    fixed = list(S0 = 990, I0 = 10, R0 = 0),
    prior = prior,
    compare = list(
      cases = observe("Poisson", transition = "inf", noise = "exp_noise")
    ),
    type = "stochastic",
    n_particles = 50,
    n_steps = 20,
    n_chains = 1,
    initial = c(5e-6, 0.1),
    vcv = diag(c(1e-13, 1e-4))
  )

  expect_true(inherits(samples, "monty_samples"))
  expect_equal(nrow(samples$pars), 2)   # 2 parameters
  expect_equal(ncol(samples$pars), 20)  # 20 steps
})

# ============================================================================
# Compilation test: compare code compiles in odin2
# ============================================================================

test_that("stochastic PN with compare compiles in odin2", {
  skip_if_not_installed("odin2")

  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )

  code <- pn_to_odin(sir, type = "stochastic",
    compare = list(
      cases = observe("Poisson", transition = "inf", noise = "exp_noise")
    ))

  # This should compile without error
  gen <- odin2::odin(code)
  expect_true(!is.null(gen))
})

test_that("stochastic SF with compare compiles in odin2", {
  skip_if_not_installed("odin2")

  sir <- stock_and_flow(
    stocks = c("S", "I", "R"),
    flows = list(
      infection = flow(from = "S", to = "I", rate = beta * S * I / N),
      recovery = flow(from = "I", to = "R", rate = gamma * I)
    ),
    sums = list(N = c("S", "I", "R")),
    params = c("beta", "gamma")
  )

  code <- sf_to_odin(sir, type = "stochastic",
    compare = list(
      cases = observe("Poisson", flow = "infection")
    ))

  gen <- odin2::odin(code)
  expect_true(!is.null(gen))
})
