#' Doses Eliminated by the Final Data of Each Trial
#'
#' @description
#'   Internal helper. Re-derive the dose elimination from the final data, as
#'   \code{boin_select_mtd()} does: a dose given to at least three patients
#'   whose posterior probability of exceeding the target, under a uniform
#'   prior, is above the cutoff is eliminated together with every higher
#'   dose, and the stricter safety rule at the lowest dose eliminates every
#'   dose.
#'
#' @param n_pts Matrix of the patients treated, trials by doses.
#' @param n_tox Matrix of the same shape, of the DLTs or of the sums of the
#'   toxicity scores.
#' @param target Target of the design.
#' @param cutoff_eli Posterior probability cutoff for dose elimination.
#' @param extrasafe Whether to apply the safety rule at the lowest dose.
#' @param offset Amount by which \code{cutoff_eli} is relaxed for that rule.
#'
#' @return A logical matrix of the shape of \code{n_pts}.
#'
#' @importFrom stats pbeta
#'
#' @noRd
mtd_eliminated <- function(n_pts, n_tox, target, cutoff_eli = 0.95,
                           extrasafe = FALSE, offset = 0.05) {
  prob_above <- 1 - pbeta(target, n_tox + 1, n_pts - n_tox + 1)
  dim(prob_above) <- dim(n_pts)
  elim <- n_pts >= 3 & prob_above > cutoff_eli
  dim(elim) <- dim(n_pts)
  for (j in seq_len(ncol(elim))[-1L]) {
    elim[, j] <- elim[, j] | elim[, j - 1L]
  }
  if (extrasafe) {
    unsafe <- n_pts[, 1L] >= 3 & prob_above[, 1L] > cutoff_eli - offset
    elim[unsafe, ] <- TRUE
  }
  elim
}
