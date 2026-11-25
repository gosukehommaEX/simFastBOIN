test_that("uboin_utility_table reproduces Table S4 of the article", {
  # Table S4 of the supporting information: 3 patients, utilities 0, 15, 25,
  # 100. Its admissible rows are those of an efficacy prior Beta(0.05, 0.05).
  tab <- uboin_utility_table(3, tox_upper = 0.30, resp_lower = 0.20,
                             utility = c(0, 15, 25, 100),
                             prior_resp = c(0.05, 0.05))
  published <- zhou2019_table_s4()

  expect_equal(nrow(tab), choose(3 + 3, 3))
  adm <- tab[tab$admissible, ]
  key <- function(d) paste(d$n_resp, d$n_tox, d$n_resp_no_tox)
  expect_setequal(key(adm), key(published))
  at <- match(key(published), key(adm))
  # The table prints one decimal; 41.25 is printed as 41.2.
  expect_true(all(abs(adm$utility[at] - published$utility) <= 0.05 + 1e-9))
})

test_that("the default priors make doses without a response admissible", {
  # With the Beta(1, 1) prior, Pr(q < 0.2 | 0 of 3) = 1 - 0.8^4 = 0.59.
  tab <- uboin_utility_table(3, tox_upper = 0.30, resp_lower = 0.20)
  none <- tab[tab$n_resp == 0, ]
  expect_equal(nrow(none), 4L)
  expect_equal(none$admissible, c(TRUE, TRUE, TRUE, FALSE))
  expect_equal(none$toxic, c(FALSE, FALSE, FALSE, TRUE))
  expect_false(any(none$futile))
})

test_that("the counts are consistent and the utilities are posterior means", {
  tab <- uboin_utility_table(6, tox_upper = 0.30, resp_lower = 0.20)
  expect_named(tab, c("n_resp", "n_tox", "n_resp_no_tox", "utility",
                      "toxic", "futile", "admissible"))
  expect_equal(nrow(tab), choose(6 + 3, 3))
  expect_true(all(tab$n_resp_no_tox <= tab$n_resp))
  expect_true(all(tab$n_resp - tab$n_resp_no_tox <= tab$n_tox))
  expect_true(all(6 - tab$n_resp - tab$n_tox + tab$n_resp -
                    tab$n_resp_no_tox >= 0))
  expect_identical(tab$admissible, !tab$toxic & !tab$futile)
  # One response without a DLT and five patients without either outcome.
  row <- tab[tab$n_resp == 1 & tab$n_tox == 0, ]
  expect_equal(row$utility, (30 * 5.25 + 50 * 0.25 + 100 * 1.25) / 7)
  # The toxicity rule: 4 DLTs of 6 are toxic (0.971), 3 are not (0.874).
  expect_true(all(tab$toxic[tab$n_tox == 4]))
  expect_false(any(tab$toxic[tab$n_tox == 3]))
})

test_that("uboin_utility_table validates its arguments", {
  expect_error(uboin_utility_table(0, tox_upper = 0.3, resp_lower = 0.2),
               "'n' must be at least 1")
  expect_error(uboin_utility_table(3, tox_upper = 0.3, resp_lower = 0.2,
                                   utility = c(0, 50, 100)),
               "'utility' must hold four finite numbers")
  expect_error(uboin_utility_table(3, tox_upper = 0.3, resp_lower = 0.2,
                                   prior_resp = c(0, 1)),
               "'prior_resp' must hold the two positive parameters")
  expect_error(uboin_utility_table(3, tox_upper = 1.3, resp_lower = 0.2),
               "'tox_upper' must lie strictly between 0 and 1")
})
