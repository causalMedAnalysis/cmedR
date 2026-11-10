test_that("linmed returns the documented components", {
  dat <- sim_mediation_data(n = 500)
  fit <- linmed(data = dat, D = "D", M = "M", Y = "Y", C = "C")

  expect_type(fit, "list")
  expect_true(all(c("ATE", "NDE", "NIE", "CDE", "model_m", "model_y") %in%
                    names(fit)))
  expect_true(all(vapply(fit[c("ATE", "NDE", "NIE", "CDE")],
                         function(x) is.numeric(x) && length(x) == 1L,
                         logical(1))))
})

test_that("linmed recovers known effects and the ATE decomposes", {
  dat <- sim_mediation_data(n = 3000)
  fit <- linmed(data = dat, D = "D", M = "M", Y = "Y", C = "C")

  # Total effect = natural direct + natural indirect
  expect_equal(fit$ATE, fit$NDE + fit$NIE, tolerance = 1e-8)

  # Without interactions the CDE equals the NDE
  expect_equal(fit$CDE, fit$NDE, tolerance = 1e-8)

  # True values: NDE = 0.2, NIE = 0.2, ATE = 0.4 (loose Monte Carlo tolerance)
  expect_lt(abs(fit$NDE - 0.2), 0.1)
  expect_lt(abs(fit$NIE - 0.2), 0.1)
  expect_lt(abs(fit$ATE - 0.4), 0.1)
})

test_that("linmed is invariant to the sign of the exposure contrast", {
  dat <- sim_mediation_data(n = 500)
  up   <- linmed(data = dat, D = "D", M = "M", Y = "Y", C = "C",
                 d = 1, dstar = 0)
  down <- linmed(data = dat, D = "D", M = "M", Y = "Y", C = "C",
                 d = 0, dstar = 1)

  expect_equal(down$ATE, -up$ATE, tolerance = 1e-8)
  expect_equal(down$NDE, -up$NDE, tolerance = 1e-8)
  expect_equal(down$NIE, -up$NIE, tolerance = 1e-8)
})

test_that("linmed is deterministic without bootstrapping", {
  dat <- sim_mediation_data(n = 500)
  a <- linmed(data = dat, D = "D", M = "M", Y = "Y", C = "C")
  b <- linmed(data = dat, D = "D", M = "M", Y = "Y", C = "C")
  expect_equal(a[c("ATE", "NDE", "NIE", "CDE")],
               b[c("ATE", "NDE", "NIE", "CDE")])
})
