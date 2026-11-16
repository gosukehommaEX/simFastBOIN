test_that("the scores reproduce Table 3 of Lin et al. (2020)", {
  pub <- utils::read.csv(test_path("fixtures", "lin2020-table3.csv"),
                         colClasses = c("integer", "integer", "integer",
                                        "character"))
  tab <- boin12_decision_table(target = 0.35, target_resp = 0.25)
  key_tab <- paste(tab$n, tab$n_tox, tab$n_resp)
  key_pub <- paste(pub$n, pub$n_tox, pub$n_resp)
  expect_equal(nrow(pub), 166L)
  expect_setequal(key_tab, key_pub)
  score <- ifelse(tab$admissible, as.character(tab$score), "E")
  expect_identical(score[match(key_pub, key_tab)], pub$score)
})

test_that("the boundaries and the benchmark match the articles", {
  tab <- boin12_decision_table(target = 0.35, target_resp = 0.25, n = 0)
  des <- attr(tab, "design")
  expect_lt(abs(des$lambda_e - 0.276), 5e-4)
  expect_lt(abs(des$lambda_d - 0.419), 5e-4)
  expect_lt(abs(des$u_low - 41), 1e-10)
  expect_lt(abs(des$u_benchmark - 70.5), 1e-10)
  # Zhou et al. (2022): the desirability of an untried dose is 0.295
  expect_lt(abs(tab$desirability - 0.295), 5e-4)
})

test_that("the desirabilities match Table S.3 of Zhou et al. (2022)", {
  tab <- boin12_decision_table(target = 0.35, target_resp = 0.25,
                               n = c(3, 9))
  pick <- function(n, y, x) {
    tab$desirability[tab$n == n & tab$n_tox == y & tab$n_resp == x]
  }
  expect_lt(abs(pick(3, 0, 0) - 0.113), 5e-4)
  expect_lt(abs(pick(3, 1, 1) - 0.156), 5e-4)
  expect_lt(abs(pick(3, 2, 0) - 0.022), 5e-4)
  expect_lt(abs(pick(9, 3, 5) - 0.208), 5e-4)
})

test_that("scores that do not add up to 100 need the joint count", {
  tab <- boin12_decision_table(target = 0.35, target_resp = 0.25, n = 3,
                               psi00 = 30, psi11 = 50)
  expect_true("n_tox_resp" %in% names(tab))
  # The four outcomes of three patients can be combined in choose(6, 3) ways
  expect_equal(nrow(tab), 20L)
  u_low <- 100 * 0.25 * 0.65 + 30 * 0.65 * 0.75 + 50 * 0.35 * 0.25
  ub <- (0.5 * u_low + 50) / 100
  row <- tab[tab$n_tox == 1 & tab$n_resp == 1 & tab$n_tox_resp == 1, ]
  q <- (30 * 2 + 50 * 1) / 100
  expect_lt(abs(row$desirability -
                  stats::pbeta(ub, 1 + q, 1 + 3 - q, lower.tail = FALSE)),
            1e-12)
  plain <- boin12_decision_table(target = 0.35, target_resp = 0.25, n = 3)
  expect_false("n_tox_resp" %in% names(plain))
  expect_equal(nrow(plain), 16L)
})

test_that("invalid arguments are rejected", {
  expect_error(boin12_decision_table(0.35, 0.25, n = -1), "'n'")
  expect_error(boin12_decision_table(0.35, 0.25, psi00 = 120), "psi00")
  expect_error(boin12_decision_table(0.35, 0.25, benchmark_weight = 0),
               "benchmark_weight")
  expect_error(boin12_decision_table(0.35, 0.25, utility_prior = 1),
               "utility_prior")
  expect_error(boin12_decision_table(0.35, 1.2), "target_resp")
})
