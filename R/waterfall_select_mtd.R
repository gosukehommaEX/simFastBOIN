#' Select the MTD Contour from Completed Waterfall Trials
#'
#' @description
#'   Choose the maximum tolerated dose at every dose level of drug A for each
#'   of a set of completed trials of the waterfall design, from the final
#'   patient and DLT counts and the combinations excluded during the trial.
#'
#' @param n_pts
#'   Integer matrix giving the number of patients treated at each combination of
#'   one trial, with one row per dose level of drug A and one column per dose
#'   level of drug B, or an integer array with dimensions trial by level of drug
#'   A by level of drug B for several trials, as returned by
#'   \code{\link{waterfall_simulate}}.
#'
#' @param n_tox
#'   Integer matrix or array of the same shape as \code{n_pts}, giving the number
#'   of DLTs observed.
#'
#' @param eliminated
#'   Logical matrix or array of the same shape as \code{n_pts}: whether the
#'   combination is excluded from the selection, as returned by
#'   \code{\link{waterfall_simulate}}.
#'
#' @param target
#'   Numeric scalar. Target DLT probability.
#'
#' @return
#'   A data frame with one row per trial and dose level of drug A and columns
#'   \item{trial}{Trial index.}
#'   \item{dose_a}{Dose level of drug A.}
#'   \item{dose_b}{Dose level of drug B of the MTD at that level of drug A, or \code{NA} when none is selected.}
#'
#' @details
#'   The selection follows \code{BOIN::get.oc.comb()} with
#'   \code{mtd.contour = TRUE}. The estimates \code{(y + 0.05) / (n + 0.1)},
#'   set to 1.1 at the excluded combinations, are smoothed by bivariate
#'   isotonic regression with weights \code{n + 0.1}, so that they increase
#'   along every row and column, and \code{1e-5} times the sum of the two dose
#'   levels is added. From the highest level of drug A to the lowest, the
#'   treated combination of the row that is not excluded and whose estimate is
#'   closest to \code{target} is selected. A row whose first combination is
#'   excluded, or without a treated combination that is not excluded, has no
#'   MTD. When the MTD of a row is not at a higher level of drug B than the MTD
#'   of the row above it, it is moved to the level of drug B of that MTD, so
#'   that the contour does not turn back; the combination it is moved to need
#'   not have been treated.
#'
#'   The isotonic regression is \code{Iso::biviso()}, the routine the reference
#'   implementation calls.
#'
#'   \code{BOIN::select.mtd.comb()} with \code{mtd.contour = TRUE} ends with a
#'   different rule (a row takes the last level of drug B when the row above
#'   has its MTD there, and has no MTD when both are there); the operating
#'   characteristics of \code{BOIN::get.oc.comb()} come from the rule described
#'   here.
#'
#' @references
#'   Zhang, L. and Yuan, Y. (2016). A Practical Bayesian Design to Identify the
#'   Maximum Tolerated Dose Contour for Drug Combination Trials. Statistics in
#'   Medicine, 35(27), 4924-4936.
#'
#' @examples
#' n_pts <- matrix(c(3, 3, 9,
#'                   6, 9, 3), nrow = 2, byrow = TRUE)
#' n_tox <- matrix(c(0, 0, 3,
#'                   1, 3, 2), nrow = 2, byrow = TRUE)
#' eliminated <- matrix(FALSE, nrow = 2, ncol = 3)
#' eliminated[2, 3] <- TRUE
#'
#' waterfall_select_mtd(n_pts, n_tox, eliminated, target = 0.30)
#'
#' @seealso \code{\link{waterfall_simulate}}, \code{\link{sim_waterfall}}
#'
#' @export
waterfall_select_mtd <- function(n_pts, n_tox, eliminated, target) {

  n_pts <- as_comb_count_array(n_pts, "n_pts")
  n_tox <- as_comb_count_array(n_tox, "n_tox")

  if (!identical(dim(n_pts), dim(n_tox))) {
    stop("'n_pts' and 'n_tox' must have the same dimensions", call. = FALSE)
  }
  if (any(n_tox > n_pts)) {
    stop("'n_tox' must not exceed 'n_pts' at any combination", call. = FALSE)
  }
  if (!is.logical(eliminated) || any(is.na(eliminated))) {
    stop("'eliminated' must be logical without missing values", call. = FALSE)
  }
  if (length(dim(eliminated)) == 2L) {
    eliminated <- array(eliminated, c(1L, dim(eliminated)))
  }
  if (!identical(dim(eliminated), dim(n_pts))) {
    stop("'eliminated' must have the same dimensions as 'n_pts'",
         call. = FALSE)
  }
  check_scalar_prob(target, "target")

  dims <- dim(n_pts)
  n_trials <- dims[1L]
  n_rows <- dims[2L]
  n_cols <- dims[3L]

  # A trial without patients selects nothing, and needs no regression.
  iso <- array(NA_real_, dims)
  n_by_trial <- rowSums(matrix(n_pts, nrow = n_trials))
  for (t in which(n_by_trial > 0)) {
    n <- matrix(n_pts[t, , ], n_rows, n_cols)
    y <- matrix(n_tox[t, , ], n_rows, n_cols)
    phat <- (y + 0.05) / (n + 0.1)
    phat[matrix(eliminated[t, , ], n_rows, n_cols)] <- 1.1
    iso[t, , ] <- Iso::biviso(phat, n + 0.1, warn = TRUE)
  }

  mtd <- waterfall_choose_cpp(
    iso = iso, n_pts = n_pts, eliminated = eliminated,
    target = as.numeric(target)
  )

  data.frame(
    trial = rep(seq_len(n_trials), each = n_rows),
    dose_a = rep(seq_len(n_rows), times = n_trials),
    dose_b = as.vector(t(mtd)),
    stringsAsFactors = FALSE
  )
}
