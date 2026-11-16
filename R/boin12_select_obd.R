#' Select the Optimal Biological Dose of BOIN12 Trials
#'
#' @description
#'   Select the optimal biological dose (OBD) from the final data of one or
#'   more trials of the BOIN12 design with the two-step procedure of Lin et
#'   al. (2020): the MTD from the isotonic DLT rates, then the admissible
#'   dose up to the MTD with the highest estimated utility.
#'
#' @inheritParams boin12_simulate
#'
#' @param n_pts
#'   Integer vector of patients treated at each dose, or a matrix with one row
#'   per trial.
#'
#' @param n_tox,n_resp
#'   DLT and response counts of the same shape as \code{n_pts}.
#'
#' @param n_tox_resp
#'   Counts of the patients with both a DLT and a response, of the same shape
#'   as \code{n_pts}, or \code{NULL}. Needed unless \code{psi00 + psi11} is
#'   100, in which case the utility depends on the margins alone. Defaults to
#'   \code{NULL}.
#'
#' @param est_method
#'   Character scalar. Estimate of the utility of a dose:
#'   \code{"posterior_mean"} (the default), the posterior mean of the
#'   quasi-binomial model given in Section 4 of the Supplement of Lin et al.
#'   (2020); \code{"observed"}, the observed mean utility, which the example
#'   of Section 6 of that Supplement reports; or \code{"tail_probability"},
#'   the desirability used during the trial, which the 'escalation' package
#'   uses.
#'
#' @param early_stop
#'   Logical vector or \code{NULL}. \code{TRUE} for the trials stopped early,
#'   which select no OBD. Defaults to \code{NULL}, no early stop.
#'
#' @return
#'   A list with components
#'   \item{obd}{Integer vector of the selected dose levels, \code{NA} when no
#'     dose is selected.}
#'   \item{reason}{Character vector: \code{"selected"}, \code{"early_stop"} or
#'     \code{"no_admissible_dose"}.}
#'   \item{utility}{Numeric matrix of the utility estimates (on the scale 0
#'     to 100; for \code{"tail_probability"}, 100 times the probability),
#'     \code{NA} for the doses not tried.}
#'
#' @details
#'   Only the doses with at least one patient are considered. The DLT rates
#'   are made monotone by isotonic regression with equal weights, and the MTD
#'   is the dose whose estimate is closest to \code{target}, ties to the
#'   higher dose as in the 'escalation' package. The OBD is the admissible
#'   dose, in the sense of \code{\link{boin12_simulate}}, at or below the MTD
#'   with the highest estimate, ties to the lower dose. When no dose up to the
#'   MTD is admissible, no OBD is selected.
#'
#' @references
#'   Lin, R., Zhou, Y., Yan, F., Li, D. and Yuan, Y. (2020). BOIN12: Bayesian
#'   optimal interval phase I/II trial design for utility-based dose finding
#'   in immunotherapy and targeted therapies. JCO Precision Oncology, 4,
#'   1393-1402. \doi{10.1200/PO.20.00257}
#'
#' @examples
#' # The example of Section 6 of the Supplement of Lin et al. (2020)
#' boin12_select_obd(
#'   n_pts = c(3, 3, 6, 12, 6), n_tox = c(0, 1, 1, 3, 3),
#'   n_resp = c(0, 1, 2, 6, 3), target = 0.35, target_resp = 0.25,
#'   est_method = "observed"
#' )
#'
#' @seealso \code{\link{sim_boin12}}, \code{\link{boin12_simulate}}
#'
#' @export
boin12_select_obd <- function(n_pts, n_tox, n_resp, n_tox_resp = NULL,
                              target, target_resp, cutoff_eli = 0.95,
                              cutoff_eli_resp = 0.90, psi00 = 40, psi11 = 60,
                              benchmark_weight = 0.5,
                              utility_prior = c(1, 1),
                              est_method = c("posterior_mean", "observed",
                                             "tail_probability"),
                              early_stop = NULL) {

  est_method <- match.arg(est_method)
  des <- boin12_check_design(target, target_resp, psi00 = psi00,
                             psi11 = psi11,
                             benchmark_weight = benchmark_weight,
                             utility_prior = utility_prior)
  check_scalar_prob(cutoff_eli, "cutoff_eli")
  check_scalar_prob(cutoff_eli_resp, "cutoff_eli_resp")

  n_pts <- as_count_matrix(n_pts, "n_pts")
  n_tox <- as_count_matrix(n_tox, "n_tox")
  n_resp <- as_count_matrix(n_resp, "n_resp")
  if (is.null(n_tox_resp)) {
    if (abs(psi00 + psi11 - 100) > 1e-9) {
      stop("'n_tox_resp' is needed unless psi00 + psi11 is 100",
           call. = FALSE)
    }
    n_tox_resp <- n_pts
    n_tox_resp[] <- pmax(0L, n_tox + n_resp - n_pts)
  }
  n_tox_resp <- as_count_matrix(n_tox_resp, "n_tox_resp")
  if (!identical(dim(n_pts), dim(n_tox)) ||
      !identical(dim(n_pts), dim(n_resp)) ||
      !identical(dim(n_pts), dim(n_tox_resp))) {
    stop("'n_pts', 'n_tox', 'n_resp' and 'n_tox_resp' must have the same ",
         "dimensions", call. = FALSE)
  }
  if (any(n_tox > n_pts) || any(n_resp > n_pts) ||
      any(n_tox_resp > pmin(n_tox, n_resp)) ||
      any(n_resp - n_tox_resp > n_pts - n_tox)) {
    stop("the counts are not consistent: 'n_tox' and 'n_resp' must not ",
         "exceed 'n_pts', and 'n_tox_resp' must fit both", call. = FALSE)
  }
  if (any(rowSums(n_pts) == 0L)) {
    stop("every trial must have treated at least one patient", call. = FALSE)
  }
  n_trials <- nrow(n_pts)
  if (is.null(early_stop)) early_stop <- rep(FALSE, n_trials)
  if (!is.logical(early_stop) || length(early_stop) != n_trials ||
      anyNA(early_stop)) {
    stop("'early_stop' must be TRUE or FALSE for every trial", call. = FALSE)
  }

  max_n <- max(n_pts)
  adm <- boin12_admissible_table(target, target_resp, cutoff_eli,
                                 cutoff_eli_resp, max_n)
  storage.mode(n_pts) <- "integer"
  storage.mode(n_tox) <- "integer"
  storage.mode(n_resp) <- "integer"
  storage.mode(n_tox_resp) <- "integer"

  res <- boin12_select_obd_cpp(
    n_pts = n_pts, n_tox = n_tox, n_resp = n_resp, n_tox_resp = n_tox_resp,
    early_stop = early_stop, target = as.numeric(target),
    psi00 = as.numeric(psi00), psi11 = as.numeric(psi11),
    alpha = des$utility_prior[1L], beta = des$utility_prior[2L],
    u_benchmark = des$u_benchmark / 100,
    est = match(est_method, c("posterior_mean", "observed",
                              "tail_probability")) - 1L,
    max_n = as.integer(max_n), tox_ok = adm$tox_ok, resp_ok = adm$resp_ok
  )

  utility <- res$utility
  colnames(utility) <- paste0("DL", seq_len(ncol(n_pts)))
  list(
    obd = res$obd,
    reason = c("selected", "early_stop", "no_admissible_dose")[
      res$reason + 1L],
    utility = utility
  )
}
