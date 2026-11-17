# Expected utilities were computed separately in Python (scipy): Lin et al.
# (2020) Table S4 for the independent case, and the bivariate normal
# distribution function by Owen's T function for the copula.

lin_sc1 <- list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
                p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38))
lin_sc16 <- list(p_true = c(0.25, 0.55, 0.67, 0.87, 0.95),
                 p_resp = c(0.01, 0.55, 0.55, 0.55, 0.55))

test_that("utilities and true doses follow Table S4 of Lin et al. (2020)", {
  truth <- oc_scenarios(list("Scenario 1" = lin_sc1, "Scenario 16" = lin_sc16),
                        target = 0.35, target_resp = 0.25)
  expect_s3_class(truth, "oc_scenarios")
  expect_s3_class(truth, "data.frame")
  expect_identical(levels(truth$scenario), c("Scenario 1", "Scenario 16"))
  expect_identical(names(truth), c("scenario", "dose", "p_true", "p_resp",
                                   "utility", "true_mtd", "true_obd"))
  s1 <- truth[truth$scenario == "Scenario 1", ]
  expect_equal(s1$utility, c(50.0, 56.2, 50.8, 48.2, 42.8), tolerance = 1e-12)
  expect_identical(which(s1$true_mtd), 4L)
  expect_identical(which(s1$true_obd), 2L)
  s16 <- truth[truth$scenario == "Scenario 16", ]
  expect_equal(s16$utility, c(30.6, 51.0, 46.2, 38.2, 35.0), tolerance = 1e-12)
  expect_false(any(s16$true_obd))
  expect_identical(attr(truth, "type"), "single")
  expect_identical(attr(truth, "dims"), 5L)
  expect_identical(attr(truth, "tox_limit"), 0.35)

  # Without the response limit, dose 1 of Scenario 16 is acceptable.
  no_limit <- oc_scenarios(list(lin_sc16), target = 0.35)
  expect_identical(which(no_limit$true_obd), 1L)
  expect_identical(levels(no_limit$scenario), "Scenario 1")
})

test_that("the correlation enters the utility as in the simulations", {
  scen <- list(list(p_true = c(0.10, 0.30, 0.50),
                    p_resp = c(0.60, 0.40, 0.50)))
  cop <- oc_scenarios(scen, target = 0.30, psi00 = 40, psi11 = 30,
                      resp_cor = 0.5)
  expect_equal(cop$utility, c(69.30473280549995, 46.24327939525246, 40),
               tolerance = 1e-8)
  # Pr(both) = 1 / 4 + asin(rho) / (2 pi) at the medians gives 40 exactly.
  expect_equal(cop$utility[3], 40, tolerance = 1e-9)
  bin <- oc_scenarios(scen, target = 0.30, psi00 = 40, psi11 = 30,
                      resp_cor = 0.2, resp_cor_type = "binary")
  expect_equal(bin$utility, c(69.31818369259805, 47.05300334076138, 41),
               tolerance = 1e-9)
  # A correlation of 0.3 needs Pr(response | DLT) > 1 at dose 1.
  expect_error(oc_scenarios(scen, target = 0.30, resp_cor = 0.3,
                            resp_cor_type = "binary"), "not attainable")
  # With psi00 + psi11 = 100 the utility does not depend on the correlation.
  a <- oc_scenarios(scen, target = 0.30, resp_cor = 0.6)
  b <- oc_scenarios(scen, target = 0.30)
  expect_equal(a$utility, b$utility, tolerance = 1e-9)
})

test_that("ties and given true doses are handled", {
  tied <- oc_scenarios(list(c(0.10, 0.25, 0.35, 0.50)), target = 0.30)
  expect_identical(which(tied$true_mtd), c(2L, 3L))
  expect_false("true_obd" %in% names(tied))

  given <- oc_scenarios(list(A = lin_sc1, B = lin_sc1), target = 0.35,
                        true_mtd = list(B = 3), true_obd = list(c(2, 3), 2))
  expect_identical(which(given$true_mtd[given$scenario == "A"]), 4L)
  expect_identical(which(given$true_mtd[given$scenario == "B"]), 3L)
  expect_identical(which(given$true_obd[given$scenario == "A"]), c(2L, 3L))
  none <- oc_scenarios(list(lin_sc1), target = 0.35,
                       true_obd = list(integer(0)))
  expect_false(any(none$true_obd))
})

test_that("combinations are tabulated column by column", {
  p <- matrix(c(0.05, 0.10, 0.20,
                0.10, 0.20, 0.30,
                0.20, 0.30, 0.45), nrow = 3)
  truth <- oc_scenarios(list(p), target = 0.30)
  expect_identical(attr(truth, "type"), "combination")
  expect_identical(attr(truth, "dims"), c(3L, 3L))
  expect_identical(truth$dose_a, rep(1:3, times = 3))
  expect_identical(truth$dose_b, rep(1:3, each = 3))
  expect_identical(truth$p_true, as.vector(p))
  expect_identical(cbind(truth$dose_a, truth$dose_b)[truth$true_mtd, ],
                   cbind(c(3L, 2L), c(2L, 3L)))
  given <- oc_scenarios(list(p), target = 0.30,
                        true_mtd = list(matrix(c(1, 3), ncol = 2)))
  expect_identical(which(given$true_mtd), 7L)

  pr <- matrix(0.4, 3, 3)
  with_resp <- oc_scenarios(list(list(p_true = p, p_resp = pr)),
                            target = 0.30)
  expect_true("utility" %in% names(with_resp))
  expect_identical(which(with_resp$true_obd), 1L)
})

test_that("invalid arguments are rejected", {
  expect_error(oc_scenarios(list(), target = 0.3), "non-empty list")
  expect_error(oc_scenarios(list(lin_sc1, lin_sc1$p_true), target = 0.3),
               "every scenario or none")
  expect_error(oc_scenarios(list(lin_sc1), target = 0.3, resp_cor = 1),
               "strictly between -1 and 1")
  expect_error(oc_scenarios(list(lin_sc1), target = 0.3, psi11 = 120),
               "between 0 and 100")
  expect_error(oc_scenarios(list(lin_sc1$p_true), target = 0.3,
                            true_obd = list(2)), "needs scenarios")
  expect_error(oc_scenarios(list(lin_sc1), target = 0.3,
                            true_mtd = list(6)), "between 1 and 5")
  expect_error(oc_scenarios(list(lin_sc1), target = 0.3,
                            true_mtd = list(X = 2)), "do not exist")
  expect_error(oc_scenarios(list(lin_sc1, lin_sc1), target = 0.3,
                            true_mtd = list(2)), "one element per scenario")
  expect_error(oc_scenarios(list(lin_sc1), target = 1.2), "strictly between")
})
