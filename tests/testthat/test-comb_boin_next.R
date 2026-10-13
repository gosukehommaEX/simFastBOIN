# The expected decisions of the examples below were obtained independently
# with a standalone build of the C++ core that is checked against a port of
# BOIN::get.oc.comb() to Python; the posterior probabilities come from scipy.

zero_matrix <- function() matrix(0, nrow = 3, ncol = 4)

test_that("a tie between the two candidates returns both", {
  n_pts <- zero_matrix()
  n_pts[1, 1] <- 3
  nxt <- comb_boin_next(target = 0.25, n_pts = n_pts, n_tox = zero_matrix(),
                        current = c(1, 1))

  expect_s3_class(nxt, "comb_boin_next")
  expect_equal(nxt$decision, "escalate")
  expect_true(nxt$tie)
  expect_equal(unname(nxt$next_dose), rbind(c(2L, 1L), c(1L, 2L)))
  expect_equal(nxt$candidates$score[1], nxt$candidates$score[2])
})

test_that("an escalation and a de-escalation go to the better candidate", {
  n_pts <- rbind(c(6, 3, 0, 0), c(3, 0, 0, 0), c(0, 0, 0, 0))
  n_tox <- rbind(c(1, 0, 0, 0), c(1, 0, 0, 0), c(0, 0, 0, 0))
  up <- comb_boin_next(target = 0.25, n_pts = n_pts, n_tox = n_tox,
                       current = c(1, 1))
  expect_equal(up$decision, "escalate")
  expect_false(up$tie)
  expect_equal(unname(up$next_dose), matrix(c(2L, 1L), nrow = 1))
  expect_equal(up$candidates$score, c(0.168740, 0.105060), tolerance = 1e-5)

  n_pts <- rbind(c(6, 6, 0, 0), c(3, 3, 0, 0), c(0, 0, 0, 0))
  n_tox <- rbind(c(0, 1, 0, 0), c(1, 2, 0, 0), c(0, 0, 0, 0))
  down <- comb_boin_next(target = 0.25, n_pts = n_pts, n_tox = n_tox,
                         current = c(2, 2))
  expect_equal(down$decision, "de-escalate")
  expect_equal(unname(down$next_dose), matrix(c(1L, 2L), nrow = 1))
  expect_equal(down$candidates$prob_in_interval,
               c(0.2180521404543201, 0.16724024679868177), tolerance = 1e-8)
  expect_equal(down$candidates$score, c(0.221052, 0.168740), tolerance = 1e-5)
  expect_equal(down$boundary$deescalate_at_least, 1L)
})

test_that("the trial stops for safety and for n_earlystop", {
  n_pts <- zero_matrix()
  n_pts[1, 1] <- 3
  n_tox <- zero_matrix()
  n_tox[1, 1] <- 3
  unsafe <- comb_boin_next(target = 0.25, n_pts = n_pts, n_tox = n_tox,
                           current = c(1, 1))
  expect_equal(unsafe$decision, "stop")
  expect_equal(unsafe$stop_reason, "lowest_dose_eliminated")
  expect_equal(nrow(unsafe$next_dose), 0L)
  expect_true(all(unsafe$eliminated))

  n_pts <- rbind(c(3, 3, 0, 0), c(3, 12, 0, 0), c(0, 0, 0, 0))
  n_tox <- rbind(c(0, 0, 0, 0), c(0, 3, 0, 0), c(0, 0, 0, 0))
  full <- comb_boin_next(target = 0.25, n_pts = n_pts, n_tox = n_tox,
                         current = c(2, 2), n_earlystop = 12)
  expect_equal(full$decision, "stop")
  expect_equal(full$stop_reason, "n_earlystop")
})

test_that("the dose stays between the boundaries and when no move is available", {
  n_pts <- rbind(c(3, 3, 0, 0), c(3, 9, 0, 0), c(0, 0, 0, 0))
  n_tox <- rbind(c(0, 0, 0, 0), c(0, 2, 0, 0), c(0, 0, 0, 0))
  between <- comb_boin_next(target = 0.25, n_pts = n_pts, n_tox = n_tox,
                            current = c(2, 2))
  expect_equal(between$decision, "stay")
  expect_equal(unname(between$next_dose), matrix(c(2L, 2L), nrow = 1))
  expect_null(between$candidates)
  expect_true(is.na(between$note))

  top <- comb_boin_next(target = 0.25,
                        n_pts = rbind(c(3, 3, 3, 3), c(3, 3, 3, 3),
                                      c(3, 3, 3, 6)),
                        n_tox = zero_matrix(), current = c(3, 4))
  expect_equal(top$decision, "stay")
  expect_false(any(top$candidates$available))
  expect_match(top$note, "escalation")

  n_pts <- zero_matrix()
  n_pts[1, 1] <- 3
  n_tox <- zero_matrix()
  n_tox[1, 1] <- 2
  lowest <- comb_boin_next(target = 0.25, n_pts = n_pts, n_tox = n_tox,
                           current = c(1, 1))
  expect_equal(lowest$decision, "stay")
  expect_match(lowest$note, "de-escalation")
  expect_false(any(lowest$eliminated))
})

test_that("comb_boin_next takes the decisions of the simulated trials", {
  # A trial simulated with k cohorts consumes the same random variates as the
  # first k cohorts of a longer one, so running k = 1, 2, ... reveals the
  # state after every cohort and where the next cohort went. At every state
  # comb_boin_next() must reach the decision the engine took.
  settings <- list(
    list(target = 0.25, n_earlystop = 12, extrasafe = FALSE,
         p_true = rbind(c(0.02, 0.04, 0.08, 0.14), c(0.08, 0.25, 0.42, 0.48),
                        c(0.25, 0.45, 0.50, 0.60))),
    list(target = 0.30, n_earlystop = 9, extrasafe = TRUE,
         p_true = rbind(c(0.30, 0.45, 0.60, 0.70), c(0.45, 0.55, 0.65, 0.75),
                        c(0.50, 0.60, 0.70, 0.80))),
    list(target = 0.30, n_earlystop = 12, extrasafe = FALSE,
         p_true = rbind(c(0.05, 0.10, 0.15, 0.30), c(0.10, 0.15, 0.30, 0.45),
                        c(0.15, 0.30, 0.45, 0.50)))
  )
  n_cohort <- 12
  counts <- c(move = 0, tie = 0, stay = 0, stop = 0)

  for (s in settings) {
    for (seed in 1:8) {
      states <- lapply(seq_len(n_cohort), function(k) {
        comb_boin_simulate(target = s$target, p_true = s$p_true, n_cohort = k,
                           cohort_size = 3, n_trials = 1,
                           n_earlystop = s$n_earlystop,
                           extrasafe = s$extrasafe, seed = seed)
      })
      before <- matrix(0L, 3, 4)
      for (k in seq_len(n_cohort - 1L)) {
        now <- states[[k]]
        n_now <- now$n_pts[1, , ]
        current <- which(n_now != before, arr.ind = TRUE)
        expect_equal(nrow(current), 1L)
        nxt <- comb_boin_next(target = s$target, n_pts = n_now,
                              n_tox = now$n_tox[1, , ],
                              current = as.vector(current),
                              n_earlystop = s$n_earlystop,
                              extrasafe = s$extrasafe)
        info <- paste("seed", seed, "cohort", k)
        if (now$stop_reason != "max_cohorts") {
          expect_equal(nxt$decision, "stop", info = info)
          expect_equal(nxt$stop_reason, now$stop_reason, info = info)
          counts["stop"] <- counts["stop"] + 1
          break
        }
        went <- which(states[[k + 1L]]$n_pts[1, , ] != n_now, arr.ind = TRUE)
        expect_true(nxt$decision != "stop", info = info)
        hit <- nxt$next_dose[, "dose_a"] == went[1, 1] &
          nxt$next_dose[, "dose_b"] == went[1, 2]
        expect_true(any(hit), info = info)
        if (nxt$decision == "stay") {
          counts["stay"] <- counts["stay"] + 1
        } else if (nxt$tie) {
          counts["tie"] <- counts["tie"] + 1
        } else {
          counts["move"] <- counts["move"] + 1
        }
        before <- n_now
      }
    }
  }
  # Every kind of decision was met, so the comparison is not vacuous.
  expect_true(all(counts > 0))
})

test_that("comb_boin_next agrees with BOIN::next.comb where the rules coincide", {
  skip_if_not_installed("BOIN")
  expect_true(have_package("BOIN"))

  target <- 0.30
  lambda_d <- boin_lambda(target)$lambda_d
  set.seed(2024)
  compared <- 0L

  for (draw in 1:400) {
    n <- matrix(sample(c(0, 3, 6, 9), 12, replace = TRUE,
                       prob = c(0.4, 0.3, 0.2, 0.1)), 3, 4)
    y <- matrix(vapply(n, function(k) {
      if (k == 0) 0 else sample(0:k, 1, prob = stats::dbinom(0:k, k, 0.3))
    }, numeric(1)), 3, 4)
    treated <- which(n > 0, arr.ind = TRUE)
    if (nrow(treated) == 0L) next
    cur <- as.vector(treated[sample.int(nrow(treated), 1L), ])
    ours <- comb_boin_next(target, n_pts = n, n_tox = y, current = cur,
                           n_earlystop = 100)

    # Skip the states where BOIN::next.comb() applies a different rule: a tie,
    # which it breaks at random; the lowest combination eliminated while
    # another one is current, which it does not stop for; and an escalation
    # towards a combination that its extra rule on observed DLT rates blocks.
    if (ours$tie) next
    if (ours$eliminated[1, 1] && !all(cur == c(1, 1))) next
    cand <- ours$candidates
    escalating <- !is.null(cand) &&
      any(cand$dose_a > cur[1] | cand$dose_b > cur[2], na.rm = TRUE)
    if (escalating) {
      phat <- y / pmax(n, 1)
      block_a <- cur[1] < 3 && any(phat[cur[1] + 1, seq_len(cur[2])] >= lambda_d)
      block_b <- cur[2] < 4 && any(phat[seq_len(cur[1]), cur[2] + 1] >= lambda_d)
      if ((block_a && cand$available[1]) || (block_b && cand$available[2])) next
    }

    invisible(capture.output(
      reference <- BOIN::next.comb(target = target, npts = n, ntox = y,
                                   dose.curr = cur)$next_dc
    ))
    info <- paste("draw", draw)
    if (ours$decision == "stop") {
      expect_true(all(is.na(reference)), info = info)
    } else {
      expect_equal(as.vector(reference), unname(ours$next_dose[1, ]),
                   info = info)
    }
    compared <- compared + 1L
  }
  expect_gt(compared, 150L)
})

test_that("invalid input is rejected", {
  n_pts <- zero_matrix()
  n_pts[1, 1] <- 3
  expect_error(comb_boin_next(0.25, n_pts, zero_matrix(), current = c(2, 2)),
               "no patient")
  expect_error(comb_boin_next(0.25, n_pts, zero_matrix(), current = c(4, 1)),
               "current")
  expect_error(comb_boin_next(0.25, n_pts, matrix(0, 3, 3), current = c(1, 1)),
               "same dimensions")
  expect_error(comb_boin_next(0.25, c(3, 0), c(0, 0), current = c(1, 1)),
               "matrices")
  n_tox <- zero_matrix()
  n_tox[1, 1] <- 4
  expect_error(comb_boin_next(0.25, n_pts, n_tox, current = c(1, 1)),
               "must not exceed")
})
