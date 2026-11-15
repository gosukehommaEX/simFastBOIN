test_that("Table 2 of Takeda et al. (2018) is reproduced", {
  tab <- boin_et_decision_table(target = 0.30, target_resp = 0.60, n = c(3, 6))
  expect_s3_class(tab, "boin_et_decision_table")

  decisions <- function(k) {
    part <- tab[tab$n == k, ]
    out <- matrix("", k + 1, k + 1)
    out[cbind(part$n_tox + 1, part$n_resp + 1)] <- part$decision
    out
  }

  # Rows are DLTs and columns responses, as transcribed from the article.
  table_3 <- rbind(
    c("E", "E", "S", "S"),
    c("E/S/D", "E/S/D", "S", "S"),
    c("D", "D", "D", "D"),
    c("D", "D", "D", "D")
  )
  esd <- c("E/S/D", "E/S/D", "E/S/D", "S", "S", "S", "S")
  table_6 <- rbind(
    c("E", "E", "E", "S", "S", "S", "S"),
    esd,
    esd,
    rep("D", 7), rep("D", 7), rep("D", 7), rep("D", 7)
  )
  dimnames(table_6) <- NULL

  expect_identical(decisions(3), table_3)
  expect_identical(decisions(6), table_6)
})

test_that("the table covers every count and agrees with the boundaries", {
  tab <- boin_et_decision_table(target = 0.25, target_resp = 0.40, n = c(1, 4, 9))
  bd <- attr(tab, "boundary")

  expect_equal(nrow(tab), sum((c(1, 4, 9) + 1)^2))
  pt <- tab$n_tox / tab$n
  pe <- tab$n_resp / tab$n
  expected <- ifelse(pt <= bd$lambda1 & pe <= bd$eta1, "E",
              ifelse(pt < bd$lambda2 & pe > bd$eta1, "S",
              ifelse(pt >= bd$lambda2, "D", "E/S/D")))
  expect_identical(tab$decision, expected)
})

test_that("the admissibility flags follow the posterior probabilities", {
  tab <- boin_et_decision_table(target = 0.30, target_resp = 0.60,
                                n = c(3, 6, 12))
  tterm <- stats::pbeta(0.30, 1 + tab$n_tox, 1 + tab$n - tab$n_tox)
  eterm <- 1 - stats::pbeta(0.36, 1 + tab$n_resp, 1 + tab$n - tab$n_resp)

  expect_identical(tab$tox_excluded, tterm < 1 - 0.95)
  expect_identical(tab$resp_excluded, eterm < 1 - 0.99)
  # 12 patients without a response exclude the dose, 9 do not.
  expect_true(tab$resp_excluded[tab$n == 12 & tab$n_resp == 0 & tab$n_tox == 0])
  tab9 <- boin_et_decision_table(target = 0.30, target_resp = 0.60, n = 9)
  expect_false(tab9$resp_excluded[tab9$n_resp == 0 & tab9$n_tox == 0])
})

test_that("invalid sample sizes are rejected", {
  expect_error(boin_et_decision_table(0.30, 0.60, n = 0), "'n'")
  expect_error(boin_et_decision_table(0.30, 0.60, n = 2.5), "'n'")
})

test_that("the elimination rules can use other probabilities", {
  # The rules of the simulations of Takeda et al. (2018).
  tab <- boin_et_decision_table(target = 0.30, target_resp = 0.60,
                                n = c(3, 6, 9), cutoff_eli = 0.90,
                                cutoff_eli_resp = 0.90, p_eli = 0.40,
                                p_eli_resp = 0.20)
  tterm <- stats::pbeta(0.40, 1 + tab$n_tox, 1 + tab$n - tab$n_tox)
  eterm <- 1 - stats::pbeta(0.20, 1 + tab$n_resp, 1 + tab$n - tab$n_resp)

  expect_identical(tab$tox_excluded, tterm < 1 - 0.90)
  expect_identical(tab$resp_excluded, eterm < 1 - 0.90)
  expect_identical(attr(tab, "p_eli"), 0.40)
  expect_error(boin_et_decision_table(0.30, 0.60, n = 3, p_eli = 1.5), "p_eli")
})
