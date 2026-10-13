#' Select the MTD Combination from Completed Drug-Combination Trials
#'
#' @description
#'   Choose the maximum tolerated dose combination for each of a set of
#'   completed trials of the BOIN drug-combination design from the final
#'   patient and DLT counts.
#'
#' @param n_pts
#'   Integer matrix giving the number of patients treated at each combination of
#'   one trial, with one row per dose level of drug A and one column per dose
#'   level of drug B, or an integer array with dimensions trial by level of drug
#'   A by level of drug B for several trials, as returned by
#'   \code{\link{comb_boin_simulate}}.
#'
#' @param n_tox
#'   Integer matrix or array of the same shape as \code{n_pts}, giving the number
#'   of DLTs observed.
#'
#' @param target
#'   Numeric scalar. Target DLT probability.
#'
#' @param cutoff_eli
#'   Numeric scalar. Posterior probability cutoff for eliminating a combination.
#'   Defaults to 0.95.
#'
#' @param extrasafe
#'   Logical scalar. Apply the stricter safety rule at the lowest combination, in
#'   which case no combination is selected when it is judged too toxic.
#'   Defaults to \code{FALSE}.
#'
#' @param offset
#'   Numeric scalar between 0 and 0.5. Amount by which \code{cutoff_eli} is
#'   relaxed for the safety rule. Defaults to 0.05.
#'
#' @param bound_mtd
#'   Logical scalar. Require the isotonic estimate at the selected combination
#'   to be below the de-escalation boundary. Defaults to \code{FALSE}.
#'
#' @param p_tox
#'   Numeric scalar. Lowest DLT probability deemed overly toxic, used to derive
#'   the de-escalation boundary when \code{bound_mtd} is \code{TRUE}.
#'   Defaults to \code{1.4 * target}.
#'
#' @return
#'   A data frame with one row per trial and columns
#'   \item{trial}{Trial index.}
#'   \item{mtd_a}{Dose level of drug A in the selected combination, or \code{NA} when none is selected.}
#'   \item{mtd_b}{Dose level of drug B in the selected combination, or \code{NA}.}
#'   \item{reason}{Why a combination was or was not selected.}
#'
#' @details
#'   The selection follows the reference implementation used by
#'   \code{BOIN::get.oc.comb()}. Elimination is re-derived from the final data:
#'   within each row, the first combination with at least three patients whose
#'   posterior probability of a DLT rate above \code{target} exceeds
#'   \code{cutoff_eli} is eliminated together with the rest of its row and the
#'   rest of its column. If the lowest combination is eliminated, no
#'   combination is selected. Otherwise the estimates
#'   \code{(y + 0.05) / (n + 0.1)} are smoothed by bivariate isotonic regression
#'   with weights
#'   \code{n + 0.1}, so that they increase along every row and column, and the
#'   treated, non-eliminated combination whose estimate is closest to
#'   \code{target} is selected. A small amount, \code{1e-5} times the sum of
#'   the two dose levels, is added to every estimate first, so that ties go to
#'   the lower combination.
#'
#'   The isotonic regression is \code{Iso::biviso()}, the routine the reference
#'   implementation calls. Combinations that the regression pools together have
#'   equal estimates in exact arithmetic but can differ in the last digits, and
#'   with the small amount above, these digits decide between combinations on
#'   the same anti-diagonal. Another algorithm for the same regression would
#'   therefore select differently in a few trials out of a thousand.
#'
#'   The article of Lin and Yin (2017) estimates the DLT rates as \code{y / n}
#'   and resolves ties differently; the package follows the reference
#'   implementation so that the two can be compared trial by trial.
#'
#'   Values of \code{reason} are \code{"selected"},
#'   \code{"lowest_dose_eliminated"} and \code{"no_dose_below_bound"}.
#'
#' @references
#'   Lin, R. and Yin, G. (2017). Bayesian Optimal Interval Design for Dose
#'   Finding in Drug-Combination Trials. Statistical Methods in Medical
#'   Research, 26(5), 2155-2167.
#'
#' @examples
#' n_pts <- matrix(c(3, 6, 3,
#'                   3, 3, 0), nrow = 2, byrow = TRUE)
#' n_tox <- matrix(c(0, 1, 1,
#'                   0, 2, 0), nrow = 2, byrow = TRUE)
#'
#' comb_boin_select_mtd(n_pts, n_tox, target = 0.30)
#'
#' @seealso \code{\link{comb_boin_simulate}}, \code{\link{sim_comb_boin}}
#'
#' @export
comb_boin_select_mtd <- function(n_pts, n_tox, target, cutoff_eli = 0.95,
                                 extrasafe = FALSE, offset = 0.05,
                                 bound_mtd = FALSE, p_tox = NULL) {

  n_pts <- as_comb_count_array(n_pts, "n_pts")
  n_tox <- as_comb_count_array(n_tox, "n_tox")

  if (!identical(dim(n_pts), dim(n_tox))) {
    stop("'n_pts' and 'n_tox' must have the same dimensions", call. = FALSE)
  }
  if (any(n_tox > n_pts)) {
    stop("'n_tox' must not exceed 'n_pts' at any combination", call. = FALSE)
  }
  check_scalar_prob(cutoff_eli, "cutoff_eli")
  check_flag(extrasafe, "extrasafe")
  check_flag(bound_mtd, "bound_mtd")

  lambda_d <- boin_lambda(target, p_tox = p_tox)$lambda_d

  dims <- dim(n_pts)
  n_trials <- dims[1L]
  n_rows <- dims[2L]
  n_cols <- dims[3L]

  eliminated <- comb_boin_eliminate_cpp(
    n_pts = n_pts, n_tox = n_tox, target = as.numeric(target),
    cutoff_eli = as.numeric(cutoff_eli), extrasafe = extrasafe,
    offset = as.numeric(offset)
  )

  # The isotonic regression is needed only when the lowest combination
  # survives, which is also when the reference implementation computes it.
  iso <- array(NA_real_, dims)
  for (t in which(!eliminated[, 1L, 1L])) {
    n <- matrix(n_pts[t, , ], n_rows, n_cols)
    y <- matrix(n_tox[t, , ], n_rows, n_cols)
    iso[t, , ] <- Iso::biviso((y + 0.05) / (n + 0.1), n + 0.1, warn = TRUE)
  }

  res <- comb_boin_choose_cpp(
    iso = iso, n_pts = n_pts, eliminated = eliminated,
    target = as.numeric(target), bound_mtd = bound_mtd,
    lambda_d = as.numeric(lambda_d)
  )

  reasons <- c("selected", "lowest_dose_eliminated", "no_dose_below_bound")

  data.frame(
    trial = seq_len(n_trials),
    mtd_a = res$mtd_a,
    mtd_b = res$mtd_b,
    reason = reasons[res$reason + 1L],
    stringsAsFactors = FALSE
  )
}
