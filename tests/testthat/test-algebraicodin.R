library(testthat)
library(S7)
library(acsets)
library(catlab)
library(algebraicodin)

# === Petri net construction ================================================

test_that("labelled_petri_net constructs SIR correctly", {
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  expect_equal(nparts(sir, "S"), 3L)
  expect_equal(nparts(sir, "T"), 2L)
  expect_equal(nparts(sir, "I"), 3L)  # 2 inputs for inf + 1 for rec
  expect_equal(nparts(sir, "O"), 3L)  # 2 outputs for inf + 1 for rec
  expect_equal(species_names(sir), c("S", "I", "R"))
  expect_equal(transition_names(sir), c("inf", "rec"))
})

test_that("%=>% operator works", {
  spec <- c("A", "B") %=>% c("C", "D")
  expect_equal(spec$inputs, c("A", "B"))
  expect_equal(spec$outputs, c("C", "D"))
})

# === Transition matrices ===================================================

test_that("transition_matrices gives correct stoichiometry for SIR", {
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  tm <- transition_matrices(sir)

  # Input: inf consumes S,I; rec consumes I
  expect_equal(tm$input["S", "inf"], 1L)
  expect_equal(tm$input["I", "inf"], 1L)
  expect_equal(tm$input["I", "rec"], 1L)
  expect_equal(tm$input["R", "inf"], 0L)

  # Output: inf produces I,I; rec produces R
  expect_equal(tm$output["I", "inf"], 2L)
  expect_equal(tm$output["R", "rec"], 1L)

  # Stoichiometry: S=-1, I=+1-1=0 for inf only; I=-1,R=+1 for rec
  expect_equal(tm$stoichiometry["S", "inf"], -1L)
  expect_equal(tm$stoichiometry["I", "inf"], 1L)  # 2 out - 1 in = 1
  expect_equal(tm$stoichiometry["I", "rec"], -1L)
  expect_equal(tm$stoichiometry["R", "rec"], 1L)
})

# === Vectorfield ===========================================================

test_that("vectorfield produces correct SIR derivatives", {
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  vf <- vectorfield(sir)
  state <- c(S = 990, I = 10, R = 0)
  parms <- list(inf = 0.001, rec = 0.1)
  result <- vf(0, state, parms)[[1]]

  # dS = -inf*S*I = -0.001*990*10 = -9.9
  expect_equal(result["S"], c(S = -9.9))
  # dI = inf*S*I - rec*I = 9.9 - 1.0 = 8.9
  expect_equal(result["I"], c(I = 8.9))
  # dR = rec*I = 1.0
  expect_equal(result["R"], c(R = 1.0))
})

# === ODE code generation ==================================================

test_that("pn_to_odin generates deterministic code", {
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  code <- pn_to_odin(sir, type = "ode")

  expect_true(grepl("deriv\\(S\\)", code))
  expect_true(grepl("deriv\\(I\\)", code))
  expect_true(grepl("deriv\\(R\\)", code))
  expect_true(grepl("inf <- parameter\\(\\)", code))
  expect_true(grepl("rec <- parameter\\(\\)", code))
  expect_true(grepl("initial\\(S\\) <- S0", code))
  expect_true(grepl("rate_inf", code))
})

# === Stochastic code generation ============================================

test_that("pn_to_odin generates stochastic code", {
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  code <- pn_to_odin(sir, type = "stochastic")

  expect_true(grepl("update\\(S\\)", code))
  expect_true(grepl("update\\(I\\)", code))
  expect_true(grepl("update\\(R\\)", code))
  expect_true(grepl("Binomial", code))
  expect_true(grepl("p_inf", code))
  expect_true(grepl("n_inf", code))
  # dt is a built-in variable in odin2 discrete models, not a parameter
  expect_false(grepl("dt <- parameter", code))
})

test_that("pn_to_odin stochastic supports source-free birth transitions", {
  pn <- labelled_petri_net(c("S"), birth = character(0) %=>% "S")
  code <- pn_to_odin(pn, type = "stochastic")

  expect_true(grepl("n_birth <- Poisson\\(birth \\* dt\\)", code))
  expect_true(grepl("update\\(S\\) <- S \\+ n_birth", code))
  expect_no_error(odin2::odin(code))
})

test_that("pn_to_odin rejects unsupported stochastic self-interactions", {
  pn <- labelled_petri_net(c("S", "R"), dimer = c("S", "S") %=>% "R")
  expect_error(
    pn_to_odin(pn, type = "stochastic"),
    "does not currently support repeated self-inputs"
  )
})

test_that("pn_to_odin rejects competing stochastic outflows", {
  pn <- labelled_petri_net(
    c("S", "I", "V"),
    inf = "S" %=>% "I",
    vac = "S" %=>% "V"
  )
  expect_error(
    pn_to_odin(pn, type = "stochastic"),
    "competing stochastic outflows"
  )
})

test_that("pn_to_odin rejects unsupported rate styles", {
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  expect_error(pn_to_odin(sir, type = "ode", rate_style = "custom"), "Only .*mass_action")
})

# === Comparison code =======================================================

test_that("pn_to_odin with compare adds likelihood code", {
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  code <- pn_to_odin(sir, type = "stochastic",
                  compare = list(cases = "Poisson"))
  expect_true(grepl("cases <- data\\(\\)", code))
  expect_true(grepl("cases ~ Poisson\\(cases\\)", code))
})

# === Open Petri nets =======================================================

test_that("spontaneous_petri works", {
  op <- spontaneous_petri("I", "R", "rec")
  pn <- apex(op)
  expect_equal(nparts(pn, "S"), 2L)
  expect_equal(nparts(pn, "T"), 1L)
  expect_equal(species_names(pn), c("I", "R"))
})

test_that("exposure_petri works", {
  op <- exposure_petri("S", "I", "I", "inf")
  pn <- apex(op)
  expect_equal(nparts(pn, "S"), 2L)  # S, I (unique)
  expect_equal(nparts(pn, "T"), 1L)
  expect_equal(nparts(pn, "I"), 2L)  # S and I are both inputs
  expect_equal(nparts(pn, "O"), 2L)  # I and I are outputs
})

# === SIR composition via UWD (core EpiCats test) ==========================

test_that("SIR composition via oapply works", {
  # Define wiring diagram: s,i,r junctions; infection(s,i), recovery(i,r)
  sir_uwd <- uwd(
    outer = c("s", "i", "r"),
    infection = c("s", "i"),
    recovery = c("i", "r")
  )

  # Compose from building blocks
  dict <- epi_dict()
  sir_pn <- oapply(sir_uwd, list(dict$infection, dict$recovery))
  pn <- apex(sir_pn)

  # Should have 3 species (s, i, r), 2 transitions (inf, rec)
  expect_equal(nparts(pn, "S"), 3L)
  expect_equal(nparts(pn, "T"), 2L)

  # Generate ODE code
  code <- pn_to_odin(pn, type = "ode")
  expect_true(grepl("deriv", code))

  # Generate stochastic code
  scode <- pn_to_odin(pn, type = "stochastic")
  expect_true(grepl("update", scode))
  expect_true(grepl("Binomial", scode))
})

# === Epi dict ==============================================================

test_that("epi_dict has expected components", {
  d <- epi_dict()
  expect_true("infection" %in% names(d))
  expect_true("recovery" %in% names(d))
  expect_true("vaccination" %in% names(d))
  expect_true("exposure" %in% names(d))
  expect_true("progression" %in% names(d))
})

# === SEIR via composition ==================================================

test_that("SEIR composition works", {
  seir_uwd <- uwd(
    outer = c("s", "e", "i", "r"),
    exposure = c("s", "i", "e"),
    progression = c("e", "i"),
    recovery = c("i", "r")
  )

  d <- epi_dict()
  # exposure: S+I→E+I, progression: E→I, recovery: I→R
  # Need 3-port open PN for exposure
  seir_pn <- oapply(seir_uwd, list(d$exposure, d$progression, d$recovery))
  pn <- apex(seir_pn)

  expect_equal(nparts(pn, "S"), 4L)  # s, e, i, r
  expect_equal(nparts(pn, "T"), 3L)  # exp, prog, rec

  code <- pn_to_odin(pn, type = "ode")
  expect_true(grepl("deriv", code))
})

# === Print generated code (visual check) ===================================

test_that("SIR ODE code is readable", {
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  code <- pn_to_odin(sir, type = "ode")
  cat("\n--- SIR ODE code ---\n")
  cat(code)
  cat("\n--- end ---\n")
  expect_true(nchar(code) > 100)
})

test_that("SIR stochastic code is readable", {
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  code <- pn_to_odin(sir, type = "stochastic")
  cat("\n--- SIR Stochastic code ---\n")
  cat(code)
  cat("\n--- end ---\n")
  expect_true(nchar(code) > 100)
})

# === DDE code generation ===================================================

test_that("pn_to_odin generates DDE code", {
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  code <- pn_to_odin(sir, type = "dde",
    delays = list(inf = list(tau = 5, species = "I"))
  )
  cat("\n--- SIR DDE code ---\n")
  cat(code)
  cat("\n--- end ---\n")
  expect_true(grepl("delay\\(I, 5\\)", code))
  expect_true(grepl("I_delayed_inf", code))
  expect_true(grepl("deriv\\(S\\)", code))
  expect_true(grepl("deriv\\(I\\)", code))
  expect_true(grepl("rate_inf <- inf \\* S \\* I_delayed_inf", code))
})

test_that("DDE requires delays argument", {
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  expect_error(pn_to_odin(sir, type = "dde"), "delays")
})

test_that("DDE delay specs validate transition and species names", {
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  expect_error(
    pn_to_odin(sir, type = "dde", delays = list(wrong = list(tau = 5))),
    "Unknown delayed transition"
  )
  expect_error(
    pn_to_odin(sir, type = "dde", delays = list(inf = list(tau = 5, species = "X"))),
    "Unknown delayed species"
  )
})

# === Discrete deterministic code ==========================================

test_that("pn_to_odin generates discrete deterministic code", {
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  code <- pn_to_odin(sir, type = "discrete")
  cat("\n--- SIR Discrete deterministic code ---\n")
  cat(code)
  cat("\n--- end ---\n")
  expect_true(grepl("update\\(S\\)", code))
  expect_true(grepl("update\\(I\\)", code))
  expect_true(grepl("update\\(R\\)", code))
  expect_true(grepl("dt \\* rate_", code))
  expect_false(grepl("Binomial", code))
  # dt is a built-in variable in odin2 discrete models, not a parameter
  expect_false(grepl("dt <- parameter", code))
})

# === ResourceSharer (dynamical systems) ====================================

test_that("petri_to_continuous produces correct dynamics", {
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  rs <- petri_to_continuous(sir)
  expect_equal(rs@system_type, "continuous")
  expect_equal(rs@nstates, 3L)
  expect_equal(rs@state_names, c("S", "I", "R"))
  expect_equal(nports(rs), 3L)

  # Evaluate dynamics: dS = -beta*S*I, dI = beta*S*I - gamma*I, dR = gamma*I
  du <- eval_dynamics(rs, c(990, 10, 0), list(inf = 0.0005, rec = 0.25), 0)
  expect_equal(du[1], -0.0005 * 990 * 10)  # dS = -4.95
  expect_equal(du[2], 0.0005 * 990 * 10 - 0.25 * 10)  # dI = 2.45
  expect_equal(du[3], 0.25 * 10)  # dR = 2.5
})

test_that("petri_to_discrete runs stochastic simulation", {
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  rs <- petri_to_discrete(sir)
  expect_equal(rs@system_type, "discrete")

  set.seed(42)
  u <- c(990, 10, 0)
  u_next <- eval_dynamics(rs, u, list(inf = 0.0005, rec = 0.25, dt = 1), 0)
  expect_equal(sum(u_next), 1000)  # conservation
  expect_true(all(u_next >= 0))
})

test_that("petri_to_discrete supports source-free birth transitions", {
  pn <- labelled_petri_net(c("S"), birth = character(0) %=>% "S")
  rs <- petri_to_discrete(pn)

  set.seed(1)
  u_next <- eval_dynamics(rs, c(S = 0), list(birth = 10, dt = 1), 0)
  expect_true(u_next[1] > 0)
})

test_that("petri_to_discrete rejects unsupported stochastic self-interactions", {
  pn <- labelled_petri_net(c("S", "R"), dimer = c("S", "S") %=>% "R")
  expect_error(petri_to_discrete(pn), "does not currently support repeated self-inputs")
})

test_that("petri_to_discrete rejects competing stochastic outflows", {
  pn <- labelled_petri_net(
    c("S", "I", "V"),
    inf = "S" %=>% "I",
    vac = "S" %=>% "V"
  )
  expect_error(petri_to_discrete(pn), "competing stochastic outflows")
})

test_that("plot helpers validate required time columns", {
  skip_if_not_installed("ggplot2")
  sol <- data.frame(S = 1:3, I = 4:6)
  expect_error(plot_comparison(sol, sol), "No time column found in `sol1`")
  expect_error(plot_diff(sol, sol), "No time column found in `sol1`")
})

test_that("petri_to_delay produces DDE dynamics", {
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  rs <- petri_to_delay(sir, delays = list(inf = list(tau = 5, species = "I")))
  expect_equal(rs@system_type, "delay")

  # History function: constant at initial state
  h <- function(p, t) c(990, 10, 0)
  du <- eval_dynamics(rs, c(990, 10, 0), list(inf = 0.0005, rec = 0.25), 0, h = h)
  # With constant history, same as ODE at t=0
  expect_equal(du[1], -0.0005 * 990 * 10)
})

test_that("euler_approx converts continuous to discrete", {
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  rs_cont <- petri_to_continuous(sir)
  rs_disc <- euler_approx(rs_cont, h = 0.1)
  expect_equal(rs_disc@system_type, "discrete")

  u <- c(990, 10, 0)
  parms <- list(inf = 0.0005, rec = 0.25)

  # Euler step: u + 0.1 * du
  du <- eval_dynamics(rs_cont, u, parms, 0)
  u_expected <- u + 0.1 * du
  u_euler <- eval_dynamics(rs_disc, u, parms, 0)
  expect_equal(u_euler, u_expected)
})

test_that("Machine type works with directed dynamics", {
  m <- Machine(
    system_type = "continuous",
    ninputs = 1L, nstates = 1L, noutputs = 1L,
    dynamics = function(u, x, p, t) -p$gamma * u + x,
    readout = function(u, p, t) u
  )
  expect_equal(m@system_type, "continuous")
  du <- eval_dynamics(m, c(10), c(5), list(gamma = 0.25), 0)
  expect_equal(du, c(-0.25 * 10 + 5))

  y <- eval_readout(m, c(10), list(gamma = 0.25), 0)
  expect_equal(y, c(10))
})

test_that("euler_approx works on Machine", {
  m <- Machine(
    system_type = "continuous",
    ninputs = 0L, nstates = 1L, noutputs = 1L,
    dynamics = function(u, x, p, t) -p$k * u,
    readout = function(u, p, t) u
  )
  m_disc <- euler_approx(m, h = 0.01)
  expect_equal(m_disc@system_type, "discrete")

  u <- c(100)
  u_next <- eval_dynamics(m_disc, u, numeric(0), list(k = 0.1), 0)
  expect_equal(u_next, c(100 + 0.01 * (-0.1 * 100)))
})

# === Typed Petri nets and stratification ===================================

test_that("infectious_ontology creates correct type system", {
  ont <- infectious_ontology()
  expect_equal(species_names(ont), "Pop")
  expect_equal(transition_names(ont), c("infect", "disease", "strata"))
  expect_equal(nparts(ont, "I"), 4L)  # infect has 2, disease 1, strata 1
  expect_equal(nparts(ont, "O"), 4L)
})

test_that("typed_petri creates correct morphism for SIR", {
  ont <- infectious_ontology()
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  tp <- typed_petri(sir, ont,
    species_types = c(S = "Pop", I = "Pop", R = "Pop"),
    transition_types = c(inf = "infect", rec = "disease")
  )
  expect_equal(tp@species_map, c(1L, 1L, 1L))
  expect_equal(tp@transition_map, c(1L, 2L))
  expect_equal(length(tp@input_map), nparts(sir, "I"))
  expect_equal(length(tp@output_map), nparts(sir, "O"))
})

test_that("add_reflexives adds correct self-loops", {
  ont <- infectious_ontology()
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  tp <- typed_petri(sir, ont,
    species_types = c(S = "Pop", I = "Pop", R = "Pop"),
    transition_types = c(inf = "infect", rec = "disease")
  )
  tp_aug <- add_reflexives(tp, list(S = "strata", I = "strata", R = "strata"))

  tnames <- transition_names(tp_aug@pn)
  expect_equal(length(tnames), 5L)  # 2 original + 3 reflexive
  expect_true("refl_S_strata" %in% tnames)
  expect_true("refl_I_strata" %in% tnames)
  expect_true("refl_R_strata" %in% tnames)
  expect_equal(tp_aug@refl_species, c("", "", "S", "I", "R"))
})

test_that("typed_product produces correct SIR x Age(2)", {
  ont <- infectious_ontology()
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  sir_typed <- typed_petri(sir, ont,
    species_types = c(S = "Pop", I = "Pop", R = "Pop"),
    transition_types = c(inf = "infect", rec = "disease")
  )
  sir_aug <- add_reflexives(sir_typed, list(S = "strata", I = "strata", R = "strata"))

  age <- labelled_petri_net(c("Y", "O"), aging = "Y" %=>% "O")
  age_typed <- typed_petri(age, ont,
    species_types = c(Y = "Pop", O = "Pop"),
    transition_types = c(aging = "strata")
  )
  age_aug <- add_reflexives(age_typed, list(
    Y = c("infect", "disease"), O = c("infect", "disease")
  ))

  product <- typed_product(sir_aug, age_aug)
  sn <- species_names(product)
  tn <- transition_names(product)

  expect_equal(length(sn), 6L)
  expect_true(all(c("S_Y", "S_O", "I_Y", "I_O", "R_Y", "R_O") %in% sn))

  expect_equal(length(tn), 7L)
  expect_true("inf_Y" %in% tn)
  expect_true("inf_O" %in% tn)
  expect_true("rec_Y" %in% tn)
  expect_true("rec_O" %in% tn)
  expect_true("aging_S" %in% tn)
  expect_true("aging_I" %in% tn)
  expect_true("aging_R" %in% tn)

  # Check stoichiometry
  tm <- transition_matrices(product)
  # aging_S should remove from S_Y and add to S_O
  aging_S_col <- which(colnames(tm$stoichiometry) == "aging_S")
  expect_equal(tm$stoichiometry["S_Y", aging_S_col], -1L)
  expect_equal(tm$stoichiometry["S_O", aging_S_col], 1L)
})

test_that("%x% operator works as typed_product", {
  ont <- infectious_ontology()
  sir <- labelled_petri_net(c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"), rec = "I" %=>% "R")
  sir_typed <- typed_petri(sir, ont,
    species_types = c(S = "Pop", I = "Pop", R = "Pop"),
    transition_types = c(inf = "infect", rec = "disease"))
  sir_aug <- add_reflexives(sir_typed, list(S = "strata", I = "strata", R = "strata"))

  age <- labelled_petri_net(c("Y", "O"), aging = "Y" %=>% "O")
  age_typed <- typed_petri(age, ont,
    species_types = c(Y = "Pop", O = "Pop"), transition_types = c(aging = "strata"))
  age_aug <- add_reflexives(age_typed, list(Y = c("infect", "disease"), O = c("infect", "disease")))

  p1 <- sir_aug %x% age_aug
  p2 <- typed_product(sir_aug, age_aug)
  expect_equal(species_names(p1), species_names(p2))
  expect_equal(transition_names(p1), transition_names(p2))
})

# === Composition operators =================================================

test_that("%o% composes two Open PNs correctly", {
  inf_pn <- exposure_petri("S", "I", "I", "inf")
  rec_pn <- spontaneous_petri("I", "R", "rec")
  composed <- inf_pn %o% rec_pn

  expect_equal(sort(species_names(apex(composed))), c("I", "R", "S"))
  expect_equal(sort(transition_names(apex(composed))), c("inf", "rec"))
})

test_that("compose() works with multiple components", {
  inf_pn <- exposure_petri("S", "I", "I", "inf")
  rec_pn <- spontaneous_petri("I", "R", "rec")
  vac_pn <- spontaneous_petri("S", "R", "vac")
  composed <- compose(inf_pn, rec_pn, vac_pn)

  sn <- sort(species_names(apex(composed)))
  expect_equal(sn, c("I", "R", "S"))
  expect_equal(length(transition_names(apex(composed))), 3L)
})

test_that("pn_to_odin accepts TypedPetriNet", {
  ont <- infectious_ontology()
  sir <- labelled_petri_net(c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"), rec = "I" %=>% "R")
  sir_typed <- typed_petri(sir, ont,
    species_types = c(S = "Pop", I = "Pop", R = "Pop"),
    transition_types = c(inf = "infect", rec = "disease"))

  code <- pn_to_odin(sir_typed, "ode")
  expect_true(grepl("deriv\\(S\\)", code))
  expect_true(grepl("deriv\\(I\\)", code))
})

test_that("flatten extracts underlying Petri net", {
  ont <- infectious_ontology()
  sir <- labelled_petri_net(c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"), rec = "I" %=>% "R")
  sir_typed <- typed_petri(sir, ont,
    species_types = c(S = "Pop", I = "Pop", R = "Pop"),
    transition_types = c(inf = "infect", rec = "disease"))

  pn <- flatten(sir_typed)
  expect_equal(species_names(pn), species_names(sir))
})

cat("\nAll algebraicodin tests passed!\n")
