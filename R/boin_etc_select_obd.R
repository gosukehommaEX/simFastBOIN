#' Select the Optimal Dose Combination of BOIN-ETC Trials
#'
#' @description
#'   Select the optimal dose combination (ODC) of one or more trials of the
#'   BOIN-ETC design of Kakizume et al. (2024) from their final data. The
#'   DLT rates are made monotone along both drugs by bivariate isotonic
#'   regression, the maximum tolerated dose combination (MTDC) of every row
#'   bounds the candidates of that row, and the ODC is the candidate with the
#'   highest estimated utility.
#'
#' @param n_pts,n_tox,n_resp
#'   Numeric matrices for a single trial, with one row per level of drug A
#'   and one column per level of drug B, or arrays with dimensions trial by
#'   level of drug A by level of drug B, of the patients treated, the DLTs
#'   and the responses.
#'
#' @param n_tox_resp
#'   Array or matrix of the same size, or \code{NULL}. Patients with both a
#'   DLT and a response. Needed only when \code{psi00 + psi11} is not 100.
#'   Defaults to \code{NULL}.
#'
#' @param target
#'   Numeric scalar. Target DLT probability, phi.
#'
#' @param psi00,psi11
#'   Numeric scalars from 0 to 100. Utilities of a patient with neither a
#'   DLT nor a response and with both; a response without DLT scores 100 and
#'   a DLT without response 0. Default to 40 and 60, the utilities
#'   (100, 40, 60, 0) of Kakizume et al. (2024).
#'
#' @param utility_prior
#'   Numeric vector of two positive numbers. Beta prior of the utility of a
#'   combination divided by 100. Defaults to \code{c(1, 1)}.
#'
#' @param odc_tie
#'   Character scalar. How tied candidates are reported: \code{"lowest"}
#'   (the default) selects the one with the lowest sum of the two levels;
#'   \code{"all"} also marks every tied candidate in \code{tied}, which
#'   \code{\link{sim_boin_etc}} then counts as selected, as Kakizume et al.
#'   (2024) do. The ODC returned in \code{obd_a} and \code{obd_b} is the same
#'   in both cases.
#'
#' @param cutoff_eli,p_eli
#'   As in \code{\link{boin_etc_simulate}}. Used only when \code{excluded} is
#'   \code{NULL}.
#'
#' @param excluded
#'   Logical array or matrix of the same size, or \code{NULL}. Combinations
#'   excluded by the DLT rule during the trial (component \code{excluded} of
#'   \code{\link{boin_etc_simulate}}), which are never selected. Defaults to
#'   \code{NULL}: a combination is excluded when its own final DLT data fail
#'   the rule of \code{cutoff_eli}.
#'
#' @param early_stop
#'   Logical vector or \code{NULL}. Trials stopped early, which select no
#'   ODC. Defaults to \code{NULL}, no trial.
#'
#' @return
#'   A list with components
#'   \item{obd_a, obd_b}{Integer vectors, the levels of drug A and drug B of
#'     the ODC of every trial (\code{NA} when none is selected).}
#'   \item{n_tied}{Integer vector, the number of tied candidates.}
#'   \item{tied}{Logical array, the tied candidates of every trial (with
#'     \code{odc_tie = "lowest"}, only the selected ODC).}
#'   \item{mtd_col}{Integer matrix, trial by level of drug A, the level of
#'     drug B of the MTDC of every row.}
#'   \item{reason}{Character vector: \code{"selected"}, \code{"early_stop"}
#'     or \code{"no_admissible_dose"}.}
#'   \item{tox_est}{Array of the isotonic estimates of the DLT rates.}
#'   \item{utility}{Array of the estimated utilities of the tried
#'     combinations, on the scale 0 to 100.}
#'
#' @details
#'   The DLT rate of every combination, tried or not, is estimated by
#'   \code{(y + 0.05) / (n + 0.1)} (the beta(0.05, 0.05) posterior mean of
#'   Kakizume et al., 2024) and made non-decreasing in both levels by the
#'   unweighted least-squares fit of bivariate isotonic regression, as in the
#'   code of the authors (\code{Iso::biviso()} without weights), computed by
#'   the algorithm of Dykstra and Robertson (1982) until no estimate changes
#'   by \code{1e-13} in a cycle. In every row (level of drug A) the MTDC is
#'   the combination whose estimate is closest to \code{target}; among
#'   combinations within \code{1e-8} of the smallest distance, which are tied
#'   in exact arithmetic, the highest level of drug B. The candidates are the
#'   tried combinations that are not excluded and not above the MTDC of
#'   their row.
#'
#'   The utility of a candidate is estimated by
#'   \code{100 * (q + a) / (n + a + b)}, with \code{(a, b)} the
#'   \code{utility_prior} and the quasi-binomial number of events
#'   \code{q = (100 x + psi00 (n - y - x) + (psi00 + psi11 - 100) z) / 100}
#'   of \code{n} patients with \code{y} DLTs, \code{x} responses and \code{z}
#'   both, the estimate of the code of the authors (Kakizume et al., 2024, do
#'   not state it). The ODC is the candidate with the highest estimate; the
#'   candidates within \code{1e-6} of it are tied, and the one with the
#'   lowest sum of the levels, then the first in column-major order, is
#'   selected.
#'
#'   The code of the authors also lets untried combinations (utility 50) and
#'   excluded ones be selected, and counts every tied combination; on the
#'   published trials, leaving out the untried ones changes the selection of
#'   0.0 to 0.3% of the trials and leaving out the excluded ones 0.04 to
#'   0.14%. Its isotonic estimates (\code{Iso::biviso()}) leave rounding
#'   errors of about \code{1e-6}, so that which of tied MTDCs it takes
#'   depends on them; the published percentages of trials selecting a
#'   correct ODC are on average about 1 point lower than with exact ties
#'   (see \code{inst/reproduce/reproduce_Kakizume2024.R}).
#'
#' @references
#'   Kakizume, T., Takeda, K., Taguri, M. and Morita, S. (2024). BOIN-ETC: A
#'   Bayesian optimal interval design considering efficacy and toxicity to
#'   identify the optimal dose combinations. Statistical Methods in Medical
#'   Research, 33(4), 716-727. \doi{10.1177/09622802241236936}
#'
#'   Dykstra, R. L. and Robertson, T. (1982). An algorithm for isotonic
#'   regression for two or more independent variables. The Annals of
#'   Statistics, 10(3), 708-716. \doi{10.1214/aos/1176345866}
#'
#' @examples
#' n_pts <- matrix(c(3, 6, 0, 9, 6, 3), nrow = 2)
#' n_tox <- matrix(c(0, 1, 0, 2, 2, 2), nrow = 2)
#' n_resp <- matrix(c(0, 2, 0, 5, 4, 2), nrow = 2)
#' boin_etc_select_obd(n_pts, n_tox, n_resp, target = 0.30)
#'
#' @seealso \code{\link{sim_boin_etc}}, \code{\link{boin_etc_simulate}}
#'
#' @export
boin_etc_select_obd <- function(n_pts, n_tox, n_resp, n_tox_resp = NULL,
                                target, psi00 = 40, psi11 = 60,
                                utility_prior = c(1, 1),
                                odc_tie = c("lowest", "all"),
                                cutoff_eli = 0.95, p_eli = NULL,
                                excluded = NULL, early_stop = NULL) {

  odc_tie <- match.arg(odc_tie)
  check_scalar_prob(target, "target")
  boin_etc_check_utility(psi00, psi11, utility_prior)
  check_scalar_prob(cutoff_eli, "cutoff_eli")
  if (is.null(p_eli)) p_eli <- target
  check_scalar_prob(p_eli, "p_eli")

  n_pts <- as_comb_count_array(n_pts, "n_pts")
  n_tox <- as_comb_count_array(n_tox, "n_tox")
  n_resp <- as_comb_count_array(n_resp, "n_resp")
  both_weight <- psi00 + psi11 - 100
  if (abs(both_weight) < 1e-9) both_weight <- 0
  if (is.null(n_tox_resp)) {
    if (both_weight != 0) {
      stop("'n_tox_resp' is needed when psi00 + psi11 is not 100",
           call. = FALSE)
    }
    n_tox_resp <- array(0L, dim(n_pts))
  }
  n_tox_resp <- as_comb_count_array(n_tox_resp, "n_tox_resp")
  if (!identical(dim(n_pts), dim(n_tox)) ||
      !identical(dim(n_pts), dim(n_resp)) ||
      !identical(dim(n_pts), dim(n_tox_resp))) {
    stop("'n_pts', 'n_tox', 'n_resp' and 'n_tox_resp' must have the same ",
         "dimensions", call. = FALSE)
  }
  if (any(n_tox > n_pts) || any(n_resp > n_pts) ||
      any(n_tox_resp > pmin(n_tox, n_resp))) {
    stop("the counts of DLTs and responses must not exceed the patients at ",
         "any combination", call. = FALSE)
  }
  dims <- dim(n_pts)
  n_trials <- dims[1L]
  n_rows <- dims[2L]
  n_cols <- dims[3L]
  if (any(apply(n_pts, 1L, sum) == 0L)) {
    stop("every trial must have treated at least one patient", call. = FALSE)
  }
  if (is.null(early_stop)) early_stop <- rep(FALSE, n_trials)
  if (!is.logical(early_stop) || length(early_stop) != n_trials ||
      anyNA(early_stop)) {
    stop("'early_stop' must be TRUE or FALSE for every trial", call. = FALSE)
  }
  if (is.null(excluded)) {
    tab <- boin_etc_tox_table(p_eli, cutoff_eli, max(n_pts))
    excluded <- array(tab[as.vector(n_pts) * (max(n_pts) + 1L) +
                            as.vector(n_tox) + 1L] == 0L, dims)
  } else {
    if (length(dim(excluded)) == 2L) excluded <- array(excluded, dims)
    if (!is.logical(excluded) || !identical(dim(excluded), dims) ||
        anyNA(excluded)) {
      stop("'excluded' must be a logical array of the same size as 'n_pts'",
           call. = FALSE)
    }
  }

  res <- boin_etc_select_cpp(
    n_pts = n_pts, n_tox = n_tox, n_resp = n_resp, n_tox_resp = n_tox_resp,
    excluded = array(as.integer(excluded), dims), early_stop = early_stop,
    target = as.numeric(target), psi00 = as.numeric(psi00),
    both_weight = as.numeric(both_weight),
    alpha = as.numeric(utility_prior[1L]), beta = as.numeric(utility_prior[2L])
  )

  if (res$not_converged > 0L) {
    warning("the isotonic estimates of ", res$not_converged, " trial(s) ",
            "did not converge", call. = FALSE)
  }
  labels <- dimnames(n_pts)
  if (is.null(labels)) {
    labels <- list(NULL, paste0("A", seq_len(n_rows)),
                   paste0("B", seq_len(n_cols)))
  }
  tied <- array(res$tied != 0L, dims, dimnames = labels)
  if (odc_tie == "lowest") {
    tied[] <- FALSE
    ok <- !is.na(res$obd_a)
    tied[cbind(which(ok), res$obd_a[ok], res$obd_b[ok])] <- TRUE
  }
  dimnames(res$utility) <- labels
  dimnames(res$tox_est) <- labels
  dimnames(res$mtd_col) <- list(NULL, labels[[2L]])
  list(
    obd_a = res$obd_a,
    obd_b = res$obd_b,
    n_tied = res$n_tied,
    tied = tied,
    mtd_col = res$mtd_col,
    reason = c("selected", "early_stop", "no_admissible_dose")[
      res$reason + 1L],
    tox_est = res$tox_est,
    utility = res$utility
  )
}
