library(testthat)
library(algebraicodin)

test_that("as_structured_cospan converts Open to StructuredCospan", {
  pn <- algebraicodin::labelled_petri_net(c("S", "I"), infection = c("S", "I") %=>% c("I"))
  open_pn <- algebraicodin::Open(pn, legs = list(c(1L), c(2L)))
  sc <- algebraicodin::as_structured_cospan(open_pn)
  expect_s3_class(sc, "catlab::StructuredCospan")
  expect_equal(acsets::nparts(catlab::cospan_apex(sc), "S"), 2)
  expect_equal(catlab::nlegs(sc), 2)
})

test_that("as_open converts StructuredCospan back to Open", {
  pn <- algebraicodin::labelled_petri_net(c("S", "I"), infection = c("S", "I") %=>% c("I"))
  open_pn <- algebraicodin::Open(pn, legs = list(c(1L), c(2L)))
  sc <- algebraicodin::as_structured_cospan(open_pn)
  open_back <- algebraicodin::as_open(sc)
  expect_s3_class(open_back, "algebraicodin::Open")
  expect_equal(open_back@legs, list(c(1L), c(2L)))
})

test_that("open_petri_net creates StructuredCospan directly", {
  pn <- algebraicodin::labelled_petri_net(c("S", "I", "R"),
    infection = c("S", "I") %=>% c("I"),
    recovery = c("I") %=>% c("R"))
  sc <- algebraicodin::open_petri_net(pn, 1L, 2L, 3L)
  expect_equal(catlab::nlegs(sc), 3)
  expect_equal(catlab::foot_sizes(sc), c(1L, 1L, 1L))
})

test_that("catlab oapply_cospans gives same result as oapply for SIR", {
  infection <- algebraicodin::labelled_petri_net(c("S", "I"), infection = c("S", "I") %=>% c("I"))
  recovery <- algebraicodin::labelled_petri_net(c("I", "R"), recovery = c("I") %=>% c("R"))

  # catlab path
  sc_inf <- algebraicodin::open_petri_net(infection, 1L, 2L)
  sc_rec <- algebraicodin::open_petri_net(recovery, 1L, 2L)
  w <- catlab::uwd(c("s", "i", "r"), infection = c("s", "i"), recovery = c("i", "r"))
  result_sc <- catlab::oapply_cospans(w, list(sc_inf, sc_rec))
  result_pn <- catlab::cospan_apex(result_sc)

  # algebraicodin path (oapply expects a single leg per box mapping all ports)
  open_inf <- algebraicodin::Open(infection, legs = list(c(1L, 2L)))
  open_rec <- algebraicodin::Open(recovery, legs = list(c(1L, 2L)))
  result_open <- algebraicodin::oapply(w, list(open_inf, open_rec))

  # Both should produce 3 species
  expect_equal(acsets::nparts(result_pn, "S"), acsets::nparts(result_open@apex, "S"))
})
