#' Isotonic Estimates of the DLT Probability
#'
#' @description
#'   Estimate the DLT probability at each dose under the constraint that toxicity
#'   does not decrease with dose, using the pool adjacent violators algorithm.
#'
#' @param n_pts
#'   Integer matrix with one row per trial and one column per dose, or a vector
#'   for a single trial, giving the number of patients treated.
#'
#' @param n_tox
#'   Integer matrix or vector of the same shape as \code{n_pts}, giving the number
#'   of DLTs observed.
#'
#' @param admissible
#'   Logical matrix of the same shape as \code{n_pts}, or \code{NULL}. Doses that
#'   are \code{FALSE} are left out of the estimation and returned as \code{NA}.
#'   Defaults to all doses that have treated at least one patient.
#'
#' @return
#'   A numeric matrix with one row per trial and one column per dose. Doses that
#'   did not enter the estimation are \code{NA}.
#'
#' @details
#'   Pseudo-counts of 0.05 DLTs and 0.1 patients are added before the fit, and
#'   doses are pooled with inverse variance weights. Doses left out through
#'   \code{admissible} do not influence the estimates of the remaining doses,
#'   which matters when eliminated doses must be excluded before selecting the
#'   MTD.
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @examples
#' # A single trial
#' boin_isotonic(n_pts = c(3, 6, 9, 12), n_tox = c(0, 1, 3, 4))
#'
#' # Several trials at once
#' n_pts <- matrix(c(3, 6, 9, 12,
#'                   3, 6, 9, 12,
#'                   3, 6, 9, 12), nrow = 3, byrow = TRUE)
#' n_tox <- matrix(c(0, 1, 3, 4,
#'                   0, 0, 2, 3,
#'                   1, 2, 4, 6), nrow = 3, byrow = TRUE)
#' boin_isotonic(n_pts, n_tox)
#'
#' @seealso \code{\link{boin_select_mtd}}
#'
#' @export
boin_isotonic <- function(n_pts, n_tox, admissible = NULL) {

  n_pts <- as_count_matrix(n_pts, "n_pts")
  n_tox <- as_count_matrix(n_tox, "n_tox")

  if (!identical(dim(n_pts), dim(n_tox))) {
    stop("'n_pts' and 'n_tox' must have the same dimensions", call. = FALSE)
  }
  if (any(n_tox > n_pts)) {
    stop("'n_tox' must not exceed 'n_pts' at any dose", call. = FALSE)
  }

  if (is.null(admissible)) {
    admissible <- n_pts > 0L
  } else {
    if (is.null(dim(admissible))) admissible <- matrix(admissible, nrow = 1L)
    if (!is.logical(admissible) || !identical(dim(admissible), dim(n_pts))) {
      stop("'admissible' must be a logical matrix with the same dimensions as 'n_pts'",
           call. = FALSE)
    }
  }

  out <- boin_isotonic_cpp(n_pts, n_tox, admissible)
  dimnames(out) <- dimnames(n_pts)
  out
}
