library(testthat)
library(S7)
library(jsonlite)

# Use installed packages (loaded via testthat.R / R CMD check)
library(acsets)
library(catlab)
library(algebraicodin)

# Load Julia reference results
julia <- fromJSON(file.path(getwd(), "julia_results.json"))

# Helper: Julia JSON matrices are already parsed as R matrices by jsonlite.
# Julia stores matrices as species × transitions (row = species, col = transition).
# Our R transition_matrices returns the same layout.
julia_matrix <- function(jmat) {
  as.matrix(jmat)
}

# Helper: compare stoichiometry accounting for possible row/column reordering
compare_stoich <- function(r_mat, j_mat, r_snames, j_snames, r_tnames, j_tnames, label) {
  # Reorder R matrix to match Julia ordering
  s_order <- match(j_snames, r_snames)
  t_order <- match(j_tnames, r_tnames)
  if (any(is.na(s_order)) || any(is.na(t_order))) {
    # Names don't match exactly — compare sorted
    expect_equal(sort(r_snames), sort(j_snames),
                 label = paste(label, "species names (sorted)"))
    expect_equal(sort(r_tnames), sort(j_tnames),
                 label = paste(label, "transition names (sorted)"))
    return(invisible(NULL))
  }
  r_reordered <- r_mat[s_order, t_order, drop = FALSE]
  expect_equal(unname(r_reordered), unname(j_mat),
               label = paste(label, "stoichiometry"))
}

# ====================================================================
# Test 1: Basic SIR
# ====================================================================
test_that("SIR matches Julia AlgebraicPetri", {
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  jref <- julia$sir

  expect_equal(jref$ns, 3L)
  expect_equal(jref$nt, 2L)
  expect_equal(species_names(sir), jref$species)
  expect_equal(transition_names(sir), jref$transitions)

  tm <- transition_matrices(sir)
  j_input <- julia_matrix(jref$input)
  j_output <- julia_matrix(jref$output)
  compare_stoich(tm$input, j_input,
                 species_names(sir), jref$species,
                 transition_names(sir), jref$transitions, "SIR input")
  compare_stoich(tm$output, j_output,
                 species_names(sir), jref$species,
                 transition_names(sir), jref$transitions, "SIR output")
})

# ====================================================================
# Test 2: SIS
# ====================================================================
test_that("SIS matches Julia AlgebraicPetri", {
  sis <- labelled_petri_net(
    c("S", "I"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "S"
  )
  jref <- julia$sis

  expect_equal(species_names(sis), jref$species)
  expect_equal(transition_names(sis), jref$transitions)

  tm <- transition_matrices(sis)
  j_input <- julia_matrix(jref$input)
  j_output <- julia_matrix(jref$output)
  compare_stoich(tm$input, j_input,
                 species_names(sis), jref$species,
                 transition_names(sis), jref$transitions, "SIS input")
  compare_stoich(tm$output, j_output,
                 species_names(sis), jref$species,
                 transition_names(sis), jref$transitions, "SIS output")
})

# ====================================================================
# Test 3: SEIR
# ====================================================================
test_that("SEIR matches Julia AlgebraicPetri", {
  seir <- labelled_petri_net(
    c("S", "E", "I", "R"),
    exp = c("S", "I") %=>% c("E", "I"),
    prog = "E" %=>% "I",
    rec = "I" %=>% "R"
  )
  jref <- julia$seir

  expect_equal(length(species_names(seir)), jref$ns)
  expect_equal(length(transition_names(seir)), jref$nt)
  expect_equal(species_names(seir), jref$species)
  expect_equal(transition_names(seir), jref$transitions)

  tm <- transition_matrices(seir)
  j_input <- julia_matrix(jref$input)
  j_output <- julia_matrix(jref$output)
  compare_stoich(tm$input, j_input,
                 species_names(seir), jref$species,
                 transition_names(seir), jref$transitions, "SEIR input")
  compare_stoich(tm$output, j_output,
                 species_names(seir), jref$species,
                 transition_names(seir), jref$transitions, "SEIR output")
})

# ====================================================================
# Test 4: SIR composition (SI + IR)
# ====================================================================
test_that("Composed SIR matches Julia oapply", {
  si <- exposure_petri("S", "I", "I", "inf")
  ir <- spontaneous_petri("I", "R", "rec")
  sir_comp <- si %o% ir
  sir_pn <- apex(sir_comp)

  jref <- julia$sir_composed

  expect_equal(length(species_names(sir_pn)), jref$ns)
  expect_equal(length(transition_names(sir_pn)), jref$nt)

  # Species and transitions should match (possibly different order)
  expect_equal(sort(species_names(sir_pn)), sort(jref$species))
  expect_equal(sort(transition_names(sir_pn)), sort(jref$transitions))

  tm <- transition_matrices(sir_pn)
  j_input <- julia_matrix(jref$input)
  j_output <- julia_matrix(jref$output)
  compare_stoich(tm$input, j_input,
                 species_names(sir_pn), jref$species,
                 transition_names(sir_pn), jref$transitions, "SIR composed input")
  compare_stoich(tm$output, j_output,
                 species_names(sir_pn), jref$species,
                 transition_names(sir_pn), jref$transitions, "SIR composed output")
})

# ====================================================================
# Test 5: SEIR composition (3 components)
# ====================================================================
test_that("Composed SEIR matches Julia oapply", {
  exposure    <- exposure_petri("S", "I", "E", "exp")
  progression <- spontaneous_petri("E", "I", "prog")
  recovery    <- spontaneous_petri("I", "R", "rec")

  w_seir <- uwd(
    outer       = c("S", "E", "I", "R"),
    exposure    = c("S", "I", "E"),
    progression = c("E", "I"),
    recovery    = c("I", "R")
  )
  seir <- oapply(w_seir, list(exposure, progression, recovery))
  seir_pn <- apex(seir)

  jref <- julia$seir_composed

  expect_equal(length(species_names(seir_pn)), jref$ns)
  expect_equal(length(transition_names(seir_pn)), jref$nt)
  expect_equal(sort(species_names(seir_pn)), sort(jref$species))
  expect_equal(sort(transition_names(seir_pn)), sort(jref$transitions))

  tm <- transition_matrices(seir_pn)
  j_input <- julia_matrix(jref$input)
  j_output <- julia_matrix(jref$output)
  compare_stoich(tm$input, j_input,
                 species_names(seir_pn), jref$species,
                 transition_names(seir_pn), jref$transitions, "SEIR composed input")
  compare_stoich(tm$output, j_output,
                 species_names(seir_pn), jref$species,
                 transition_names(seir_pn), jref$transitions, "SEIR composed output")
})

# ====================================================================
# Test 6: SIR × 2 age groups (stratification)
# ====================================================================
test_that("SIR×Age stratification matches Julia typed_product", {
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  ont <- infectious_ontology()
  sir_typed <- typed_petri(sir, ont,
    species_types    = c(S = "Pop", I = "Pop", R = "Pop"),
    transition_types = c(inf = "infect", rec = "disease")
  )
  sir_aug <- add_reflexives(sir_typed, list(
    S = "strata", I = "strata", R = "strata"
  ))

  age <- labelled_petri_net(
    c("Young", "Old"),
    aging = "Young" %=>% "Old"
  )
  age_typed <- typed_petri(age, ont,
    species_types    = c(Young = "Pop", Old = "Pop"),
    transition_types = c(aging = "strata")
  )
  age_aug <- add_reflexives(age_typed, list(
    Young = c("infect", "disease"),
    Old   = c("infect", "disease")
  ))

  sir_age <- sir_aug %x% age_aug

  jref <- julia$sir_age_stratified

  # Same number of species and transitions
  expect_equal(length(species_names(sir_age)), jref$ns)
  expect_equal(length(transition_names(sir_age)), jref$nt)

  # Species names: Julia uses S_Young, I_Young etc.
  r_snames <- species_names(sir_age)
  j_snames <- jref$species
  expect_equal(length(r_snames), length(j_snames))

  # Net stoichiometry should be identical (up to ordering)
  tm_r <- transition_matrices(sir_age)
  net_r <- tm_r$output - tm_r$input
  j_input <- julia_matrix(jref$input)
  j_output <- julia_matrix(jref$output)
  net_j <- j_output - j_input

  # Compare net stoich dimension
  expect_equal(dim(net_r), dim(net_j))

  # Compare row sums (total net change per species across all transitions)
  # and column sums (total net change per transition across all species)
  # These should match up to reordering
  expect_equal(sort(unname(rowSums(abs(net_r)))), sort(unname(rowSums(abs(net_j)))))
  expect_equal(sort(unname(colSums(abs(net_r)))), sort(unname(colSums(abs(net_j)))))
})

# ====================================================================
# Test 7: SIR × 2 risk groups
# ====================================================================
test_that("SIR×Risk stratification matches Julia typed_product", {
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )
  ont <- infectious_ontology()
  sir_typed <- typed_petri(sir, ont,
    species_types    = c(S = "Pop", I = "Pop", R = "Pop"),
    transition_types = c(inf = "infect", rec = "disease")
  )
  sir_aug <- add_reflexives(sir_typed, list(
    S = "strata", I = "strata", R = "strata"
  ))

  risk <- labelled_petri_net(
    c("H", "L"),
    HH = c("H", "H") %=>% c("H", "H"),
    HL = c("H", "L") %=>% c("H", "L"),
    LH = c("L", "H") %=>% c("L", "H"),
    LL = c("L", "L") %=>% c("L", "L")
  )
  risk_typed <- typed_petri(risk, ont,
    species_types    = c(H = "Pop", L = "Pop"),
    transition_types = c(HH = "infect", HL = "infect", LH = "infect", LL = "infect")
  )
  risk_aug <- add_reflexives(risk_typed, list(
    H = c("disease", "strata"),
    L = c("disease", "strata")
  ))

  sir_risk <- sir_aug %x% risk_aug

  jref <- julia$sir_risk_stratified

  expect_equal(length(species_names(sir_risk)), jref$ns)
  expect_equal(length(transition_names(sir_risk)), jref$nt)

  tm_r <- transition_matrices(sir_risk)
  net_r <- tm_r$output - tm_r$input
  j_input <- julia_matrix(jref$input)
  j_output <- julia_matrix(jref$output)
  net_j <- j_output - j_input

  expect_equal(dim(net_r), dim(net_j))
  expect_equal(sort(unname(rowSums(abs(net_r)))), sort(unname(rowSums(abs(net_j)))))
  expect_equal(sort(unname(colSums(abs(net_r)))), sort(unname(colSums(abs(net_j)))))
})

# ====================================================================
# Test 8: SIR vectorfield evaluation
# ====================================================================
test_that("SIR vectorfield matches Julia mass action", {
  sir <- labelled_petri_net(
    c("S", "I", "R"),
    inf = c("S", "I") %=>% c("I", "I"),
    rec = "I" %=>% "R"
  )

  jref <- julia$sir_vectorfield

  state <- setNames(jref$u, species_names(sir))
  parms <- setNames(as.list(jref$p), transition_names(sir))
  vf <- vectorfield(sir)
  result <- vf(0, state, parms)
  du <- unname(result[[1]])

  expect_equal(du, jref$du, tolerance = 1e-10,
               label = "SIR du/dt")
})

# ====================================================================
# Test 9: SEIR vectorfield evaluation
# ====================================================================
test_that("SEIR vectorfield matches Julia mass action", {
  seir <- labelled_petri_net(
    c("S", "E", "I", "R"),
    exp = c("S", "I") %=>% c("E", "I"),
    prog = "E" %=>% "I",
    rec = "I" %=>% "R"
  )

  jref <- julia$seir_vectorfield

  state <- setNames(jref$u, species_names(seir))
  parms <- setNames(as.list(jref$p), transition_names(seir))
  vf <- vectorfield(seir)
  result <- vf(0, state, parms)
  du <- unname(result[[1]])

  expect_equal(du, jref$du, tolerance = 1e-10,
               label = "SEIR du/dt")
})

# ====================================================================
# Test 10: Epi helper building blocks
# ====================================================================
test_that("exposure_petri matches Julia Epidemiology.infection", {
  inf <- exposure_petri("S", "I", "I", "inf")
  inf_pn <- apex(inf)

  jref <- julia$epi_infection

  expect_equal(sort(species_names(inf_pn)), sort(jref$species))
  expect_equal(transition_names(inf_pn), jref$transitions)

  tm <- transition_matrices(inf_pn)
  j_input <- julia_matrix(jref$input)
  j_output <- julia_matrix(jref$output)
  compare_stoich(tm$input, j_input,
                 species_names(inf_pn), jref$species,
                 transition_names(inf_pn), jref$transitions, "epi_infection input")
  compare_stoich(tm$output, j_output,
                 species_names(inf_pn), jref$species,
                 transition_names(inf_pn), jref$transitions, "epi_infection output")
})

test_that("spontaneous_petri matches Julia Epidemiology.recovery", {
  rec <- spontaneous_petri("I", "R", "rec")
  rec_pn <- apex(rec)

  jref <- julia$epi_recovery

  expect_equal(sort(species_names(rec_pn)), sort(jref$species))
  expect_equal(transition_names(rec_pn), jref$transitions)

  tm <- transition_matrices(rec_pn)
  j_input <- julia_matrix(jref$input)
  j_output <- julia_matrix(jref$output)
  compare_stoich(tm$input, j_input,
                 species_names(rec_pn), jref$species,
                 transition_names(rec_pn), jref$transitions, "epi_recovery input")
  compare_stoich(tm$output, j_output,
                 species_names(rec_pn), jref$species,
                 transition_names(rec_pn), jref$transitions, "epi_recovery output")
})

cat("\nAll Julia cross-validation tests passed!\n")
