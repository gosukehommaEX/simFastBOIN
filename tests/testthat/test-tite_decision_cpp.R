# The simulation engine takes its decisions with tite_decide() in
# src/tite_core.h. These tests check that it agrees with
# tite_boin_decision_table() at every state and at follow-up values on both
# sides of every boundary.

# Decision of a table row at a given value of the follow-up statistic.
resolve_decision <- function(decision, esc, deesc, value) {
  out <- decision
  mixed <- grepl("/", decision, fixed = TRUE)
  first <- sub("/.*$", "", decision)
  up <- mixed & !is.na(esc) & value >= esc
  down <- mixed & !up & !is.na(deesc) & value <= deesc
  out[mixed] <- "S"
  out[up] <- first[up]
  out[down] <- "D"
  out
}

test_that("the engine takes the decisions of tite_boin_decision_table", {
  settings <- list(
    list(method = "imputation", max_pending_ratio = NULL, min_completed = NULL),
    list(method = "ess", max_pending_ratio = NULL, min_completed = NULL),
    list(method = "imputation", max_pending_ratio = 0.75, min_completed = 2),
    list(method = "ess", max_pending_ratio = 0.5, min_completed = 0)
  )
  codes <- c("E", "S", "D", "SUS", "DE")

  for (target in c(0.25, 0.30)) {
    bound <- boin_boundary(target, 15)
    b_elim <- bound$b_elim
    b_elim[is.na(b_elim)] <- 0L

    for (set in settings) {
      tab <- tite_boin_decision_table(
        target = target, max_n = 15, method = set$method,
        max_pending_ratio = set$max_pending_ratio,
        min_completed = set$min_completed
      )
      is_ess <- set$method == "ess"
      lower <- if (is_ess) tab$n - tab$n_pending else rep(0, nrow(tab))
      upper <- if (is_ess) tab$n else tab$n_pending

      # Five points spread over the range of the statistic, and points just
      # below and above each boundary that lie inside the range.
      rows <- integer(0)
      value <- numeric(0)
      for (f in c(0.001, 0.25, 0.5, 0.75, 0.999)) {
        rows <- c(rows, seq_len(nrow(tab)))
        value <- c(value, lower + (upper - lower) * f)
      }
      for (bnd in list(tab$esc_bound, tab$deesc_bound)) {
        for (shift in c(-1e-6, 1e-6)) {
          x <- bnd + shift
          keep <- !is.na(x) & x >= lower & x < upper
          rows <- c(rows, which(keep))
          value <- c(value, x[keep])
        }
      }

      expected <- resolve_decision(tab$decision[rows], tab$esc_bound[rows],
                                   tab$deesc_bound[rows], value)
      expected[tab$decision[rows] == "DE"] <- "DE"

      stft <- if (is_ess) value - lower[rows] else value
      stft[tab$n_pending[rows] == 0L] <- 0
      got <- codes[tite_decision_cpp(
        n = tab$n[rows], n_tox = tab$n_tox[rows],
        n_pending = tab$n_pending[rows], stft = stft,
        method = if (is_ess) 1L else 0L, target = target,
        lambda_e = bound$lambda_e, lambda_d = bound$lambda_d,
        max_pending_ratio = attr(tab, "max_pending_ratio"),
        min_completed = attr(tab, "min_completed"),
        b_esc = bound$b_esc, b_deesc = bound$b_deesc, b_elim = as.integer(b_elim)
      ) + 1L]

      state <- paste0("target ", target, ", ", set$method, ": n = ",
                      tab$n[rows], ", DLTs = ", tab$n_tox[rows], ", pending = ",
                      tab$n_pending[rows], ", statistic = ", signif(value, 8))
      # Guard against a vacuous comparison: every boundary is probed.
      expect_gt(length(rows), 5L * nrow(tab))
      expect_identical(state[got != expected], character(0))
    }
  }
})
