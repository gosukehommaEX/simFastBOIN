test_that("boin_decision_table has the expected shape and vocabulary", {
  tab <- boin_decision_table(target = 0.30, max_n = 12)

  expect_equal(dim(tab), c(13L, 12L))
  expect_equal(rownames(tab), as.character(0:12))
  expect_equal(colnames(tab), as.character(1:12))
  expect_true(all(tab[!is.na(tab)] %in% c("E", "S", "D", "DE")))
})

test_that("impossible combinations are missing", {
  tab <- boin_decision_table(target = 0.30, max_n = 9)

  for (n in 1:9) {
    expect_true(all(!is.na(tab[seq_len(n + 1L), n])))
    if (n < 9) expect_true(all(is.na(tab[(n + 2L):10L, n])))
  }
})

test_that("the decision table agrees with the boundaries it is built from", {
  target <- 0.30
  max_n <- 18
  bd <- boin_boundary(target = target, max_n = max_n)
  tab <- boin_decision_table(target = target, max_n = max_n)

  for (n in seq_len(max_n)) {
    for (y in 0:n) {
      decision <- tab[y + 1L, n]
      if (!is.na(bd$b_elim[n]) && y >= bd$b_elim[n]) {
        expect_identical(decision, "DE")
      } else if (y <= bd$b_esc[n]) {
        expect_identical(decision, "E")
      } else if (y >= bd$b_deesc[n]) {
        expect_identical(decision, "D")
      } else {
        expect_identical(decision, "S")
      }
    }
  }
})

test_that("zero DLTs always escalates and all DLTs never does", {
  tab <- boin_decision_table(target = 0.30, max_n = 15)

  expect_true(all(tab[1L, ] == "E"))
  for (n in 3:15) expect_true(tab[n + 1L, n] %in% c("D", "DE"))
})
