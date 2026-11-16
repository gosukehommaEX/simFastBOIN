test_that("the OBD agrees with an independent implementation", {
  # Final data of 300 trials and their OBDs under the three estimates,
  # computed in Python with scipy from the rules of the article.
  fx <- utils::read.csv(test_path("fixtures", "boin12_select_obd.csv"))
  get <- function(p) as.matrix(fx[, paste0(p, 1:5)])
  methods <- c("posterior_mean", "observed", "tail_probability")
  for (m in methods) {
    sel <- boin12_select_obd(n_pts = get("n"), n_tox = get("y"),
                             n_resp = get("x"), n_tox_resp = get("z"),
                             target = 0.35, target_resp = 0.25,
                             est_method = m)
    expected <- fx[[paste0("obd_", m)]]
    expected[expected == 0L] <- NA_integer_
    expect_identical(sel$obd, expected)
    expect_identical(is.na(sel$obd), sel$reason != "selected")
  }
})

test_that("the example of Lin et al. (2020) selects dose 4", {
  args <- list(n_pts = c(3, 3, 6, 12, 6), n_tox = c(0, 1, 1, 3, 3),
               n_resp = c(0, 1, 2, 6, 3), target = 0.35, target_resp = 0.25)
  obs <- do.call(boin12_select_obd, c(args, list(est_method = "observed")))
  # Section 6 of the Supplement: 40.0, 46.7, 53.3, 60.0 and 50.0
  expect_equal(round(unname(obs$utility[1, ]), 1),
               c(40.0, 46.7, 53.3, 60.0, 50.0))
  expect_equal(obs$obd, 4L)
  post <- do.call(boin12_select_obd, args)
  expect_equal(round(unname(post$utility[1, ]), 1),
               c(44.0, 48.0, 52.5, 58.6, 50.0))
  expect_equal(post$obd, 4L)
})

test_that("early stops and inadmissible doses select no OBD", {
  sel <- boin12_select_obd(n_pts = rbind(c(3, 3, 0), c(3, 0, 0)),
                           n_tox = rbind(c(0, 1, 0), c(3, 0, 0)),
                           n_resp = rbind(c(1, 1, 0), c(1, 0, 0)),
                           target = 0.35, target_resp = 0.25,
                           early_stop = c(TRUE, FALSE))
  expect_true(all(is.na(sel$obd)))
  expect_equal(sel$reason, c("early_stop", "no_admissible_dose"))
  expect_true(is.na(sel$utility[1, 3]))
})

test_that("the joint count is needed for scores that do not add to 100", {
  args <- list(n_pts = c(3, 3), n_tox = c(0, 1), n_resp = c(1, 2),
               target = 0.35, target_resp = 0.25, psi00 = 30, psi11 = 50)
  expect_error(do.call(boin12_select_obd, args), "n_tox_resp")
  sel <- do.call(boin12_select_obd, c(args, list(n_tox_resp = c(0, 1))))
  # Dose 2: one patient with both (50), one response only (100) and one
  # with neither (30)
  expect_lt(abs(sel$utility[1, 2] - 100 * (1.8 + 1) / 5), 1e-10)
  expect_error(do.call(boin12_select_obd,
                       c(args, list(n_tox_resp = c(1, 1)))), "consistent")
})
