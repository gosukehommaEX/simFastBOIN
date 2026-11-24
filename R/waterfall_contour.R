#' True MTD Contour of a Dose Matrix
#'
#' @description
#'   Find the true MTD at every dose level of drug A, and classify the
#'   combinations as below, on or above the contour, as
#'   \code{BOIN::get.oc.comb()} with \code{mtd.contour = TRUE} does. Internal.
#'
#' @param p_true
#'   Numeric matrix of true DLT probabilities, one row per level of drug A.
#'
#' @param target
#'   Numeric scalar. Target DLT probability.
#'
#' @return
#'   A list with \code{true_mtd}, a two-column matrix of the dose levels of
#'   drug A and drug B of the true MTDs, \code{true_col}, the level of drug B
#'   of the true MTD of every row (\code{NA} for none), and the logical
#'   matrices \code{below}, \code{above} and \code{at}.
#'
#' @details
#'   The true MTD of a row is the combination whose DLT probability is closest
#'   to \code{target}, the first of them among ties, provided that probability
#'   does not exceed \code{target + 0.05}. A combination is below the contour
#'   when both of its dose levels are at most those of some true MTD, and above
#'   it when both are at least those of some true MTD; the true MTDs
#'   themselves and the remaining combinations are on the contour.
#'
#' @keywords internal
#' @noRd
waterfall_contour <- function(p_true, target) {
  n_rows <- nrow(p_true)
  n_cols <- ncol(p_true)
  true_col <- rep(NA_integer_, n_rows)
  for (a in seq_len(n_rows)) {
    b <- which.min(abs(p_true[a, ] - target))
    if (p_true[a, b] <= target + 0.05) true_col[a] <- b
  }
  rows <- which(!is.na(true_col))
  true_mtd <- cbind(dose_a = rows, dose_b = true_col[rows])

  below <- matrix(FALSE, n_rows, n_cols)
  above <- matrix(FALSE, n_rows, n_cols)
  for (a in rows) {
    b <- true_col[a]
    below[seq_len(a), seq_len(b)] <- TRUE
    above[a:n_rows, b:n_cols] <- TRUE
  }
  below[true_mtd] <- FALSE
  above[true_mtd] <- FALSE

  list(true_mtd = true_mtd, true_col = true_col, below = below,
       above = above, at = !(below | above))
}
