# Example of the help page: the expected values were computed with the
# Python implementation of the rules (etc_python.py of the development
# folder). Row A1 has its isotonic estimates tied at B2 and B3, so its MTDC is
# B3; the ODC is A1B3 with the utility 62.5.
ex_pts <- matrix(c(3, 6, 0, 9, 6, 3), nrow = 2)
ex_tox <- matrix(c(0, 1, 0, 2, 2, 2), nrow = 2)
ex_resp <- matrix(c(0, 2, 0, 5, 4, 2), nrow = 2)

test_that("boin_etc_select_obd selects the example as computed in Python", {
  sel <- boin_etc_select_obd(ex_pts, ex_tox, ex_resp, target = 0.30)
  expect_equal(sel$obd_a, 1L)
  expect_equal(sel$obd_b, 3L)
  expect_equal(unname(sel$mtd_col[1, ]), c(3L, 2L))
  expect_equal(sel$reason, "selected")
  expect_equal(sel$n_tied, 1L)
  expect_equal(unname(sel$utility[1, , ]),
               matrix(c(44, 52.5, NA, 680 / 11, 62.5, 52), 2),
               tolerance = 1e-12)
  pool <- (0.05 / 0.1 + 2.05 / 9.1 + 2.05 / 6.1) / 3
  expect_equal(unname(sel$tox_est[1, , ]),
               matrix(c(0.05 / 3.1, 1.05 / 6.1, pool, pool, pool, 2.05 / 3.1),
                      2), tolerance = 1e-12)
  # An excluded ODC is not selected: the next best is A2B2.
  excl <- matrix(FALSE, 2, 3)
  excl[1, 3] <- TRUE
  sel2 <- boin_etc_select_obd(ex_pts, ex_tox, ex_resp, target = 0.30,
                              excluded = excl)
  expect_equal(c(sel2$obd_a, sel2$obd_b), c(2L, 2L))
})

test_that("tied candidates give the lowest one and are all marked with 'all'", {
  n_pts <- matrix(3, 2, 2)
  n_tox <- matrix(0, 2, 2)
  n_resp <- matrix(c(3, 1, 3, 1), 2)
  low <- boin_etc_select_obd(n_pts, n_tox, n_resp, target = 0.30)
  all <- boin_etc_select_obd(n_pts, n_tox, n_resp, target = 0.30,
                             odc_tie = "all")
  expect_equal(c(low$obd_a, low$obd_b), c(1L, 1L))
  expect_equal(c(all$obd_a, all$obd_b), c(1L, 1L))
  expect_equal(low$n_tied, 2L)
  expect_equal(unname(low$tied[1, , ]), matrix(c(TRUE, FALSE, FALSE, FALSE), 2))
  expect_equal(unname(all$tied[1, , ]), matrix(c(TRUE, FALSE, TRUE, FALSE), 2))
})

test_that("psi00 + psi11 other than 100 needs the patients with both", {
  n_pts <- matrix(c(6, 0, 0, 0), 2)
  n_tox <- matrix(c(2, 0, 0, 0), 2)
  n_resp <- matrix(c(3, 0, 0, 0), 2)
  expect_error(boin_etc_select_obd(n_pts, n_tox, n_resp, target = 0.30,
                                   psi00 = 30, psi11 = 60),
               "n_tox_resp")
  sel <- boin_etc_select_obd(n_pts, n_tox, n_resp,
                             n_tox_resp = matrix(c(2, 0, 0, 0), 2),
                             target = 0.30, psi00 = 30, psi11 = 60)
  # q = (100 * 3 + 30 * 1 - 10 * 2) / 100 = 3.1 and 100 * 4.1 / 8
  expect_equal(sel$utility[1, 1, 1], 51.25, tolerance = 1e-12)
})

test_that("the selection of the fixture trials agrees with the Python rules", {
  fx <- utils::read.csv(test_path("fixtures", "boin_etc_rules.csv"),
                        stringsAsFactors = FALSE, colClasses = "character")
  dims <- list(F1 = c(4, 4), F2 = c(4, 4), F3 = c(3, 5), F4 = c(3, 5),
               F5 = c(2, 3), F6 = c(4, 3), F7 = c(3, 4))
  split <- function(x) as.integer(strsplit(x, ";", fixed = TRUE)[[1L]])
  for (s in names(dims)) {
    rows <- fx[fx$setting == s, , drop = FALSE]
    d <- c(nrow(rows), dims[[s]])
    arr <- function(col) {
      array(t(vapply(rows[[col]], split, integer(prod(dims[[s]])))), d)
    }
    sel <- boin_etc_select_obd(
      arr("n_pts"), arr("n_tox"), arr("n_resp"), arr("n_tox_resp"),
      target = 0.30, odc_tie = "all", excluded = arr("excluded") == 1L,
      early_stop = rows$stop_reason == "lowest_excluded"
    )
    expect_identical(sel$obd_a, suppressWarnings(as.integer(rows$obd_a)))
    expect_identical(sel$obd_b, suppressWarnings(as.integer(rows$obd_b)))
    expect_identical(sel$n_tied, as.integer(rows$n_tied))
    expect_identical(array(as.integer(sel$tied), d), arr("tied"))
    for (t in which(nzchar(rows$mtd_col))) {
      expect_identical(unname(sel$mtd_col[t, ]), split(rows$mtd_col[t]))
    }
  }
})

test_that("boin_etc_select_obd checks its arguments", {
  expect_error(boin_etc_select_obd(ex_pts, ex_tox + 9, ex_resp,
                                   target = 0.30), "must not exceed")
  expect_error(boin_etc_select_obd(ex_pts, ex_tox, ex_resp, target = 0.30,
                                   psi00 = 120), "psi00")
  expect_error(boin_etc_select_obd(ex_pts, ex_tox, ex_resp, target = 0.30,
                                   utility_prior = 1), "utility_prior")
  expect_error(boin_etc_select_obd(ex_pts, ex_tox, ex_resp, target = 0.30,
                                   early_stop = NA), "early_stop")
  expect_error(boin_etc_select_obd(ex_pts, ex_tox, ex_resp, target = 0.30,
                                   excluded = matrix(1, 2, 3)), "excluded")
  stopped <- boin_etc_select_obd(ex_pts, ex_tox, ex_resp, target = 0.30,
                                 early_stop = TRUE)
  expect_true(is.na(stopped$obd_a))
  expect_equal(stopped$reason, "early_stop")
})
