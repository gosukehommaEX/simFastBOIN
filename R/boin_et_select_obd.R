#' Select the Optimal Biological Dose of BOIN-ET Trials
#'
#' @description
#'   Select the optimal biological dose (OBD) from the final data of one or
#'   more trials of the BOIN-ET design, with the estimators and selection
#'   rules of the 'boinet' package.
#'
#' @inheritParams boin_et_decision_table
#'
#' @param n_pts
#'   Integer vector of patients treated at each dose, or a matrix with one row
#'   per trial.
#'
#' @param n_tox
#'   DLT counts of the same shape as \code{n_pts}.
#'
#' @param n_resp
#'   Response counts of the same shape as \code{n_pts}.
#'
#' @param est_method
#'   Character scalar. Estimator of the response rates: \code{"observed"}
#'   (the default) for the observed rates, \code{"multi_iso"} for the
#'   model-averaged unimodal isotonic regression of Lin and Yin (2017), or
#'   \code{"fp_logistic"} for the fractional polynomial logistic regression
#'   of Takeda et al. (2018), which needs the 'mfp' package.
#'
#' @param obd_method
#'   Character scalar. Selection rule: \code{"max_resp"} (the default) for the
#'   admissible dose with the highest estimated response rate up to the
#'   estimated MTD, or the dose that maximizes the utility
#'   \code{"utility_weighted"}, \code{"utility_truncated_linear"} or
#'   \code{"utility_scoring"} among the admissible doses.
#'
#' @param w1,w2
#'   Numeric scalars. Weights of the weighted utility
#'   \code{q - w1 * p - w2 * p * (p > p_tox)}, where \code{p} and \code{q} are
#'   the estimated DLT and response rates. Default to 0.33 and 1.09.
#'
#' @param t_low,t_upp
#'   Numeric scalars or \code{NULL}. DLT rates between which the truncated
#'   linear utility falls from one to zero. Default to \code{p_saf} and
#'   \code{p_tox}.
#'
#' @param e_low,e_upp
#'   Numeric scalars or \code{NULL}. Response rates between which the
#'   truncated linear utility rises from zero to one. Default to
#'   \code{p_saf_resp / 2} and \code{target_resp}.
#'
#' @param psi00,psi11
#'   Numeric scalars. Scores of the scoring utility for a patient with neither
#'   DLT nor response and for a patient with both, the scores of a response
#'   without DLT and of a DLT without response being 100 and 0. Default to 40
#'   and 60.
#'
#' @param early_stop
#'   Logical vector or \code{NULL}. \code{TRUE} for the trials stopped early,
#'   which select no OBD. Defaults to \code{NULL}, no early stop.
#'
#' @param boinet_compatible
#'   Logical scalar. Apply the rule of the 'boinet' package (version 1.6.0),
#'   in which a dose excluded by its DLT data does not exclude the doses
#'   above it. Defaults to \code{FALSE}; see Details.
#'
#' @return
#'   A list with components
#'   \item{obd}{Integer vector of the selected dose levels, \code{NA} when no
#'     dose is selected.}
#'   \item{reason}{Character vector: \code{"selected"}, \code{"early_stop"} or
#'     \code{"no_admissible_dose"}.}
#'   \item{probe}{Numeric matrix of the response estimates, \code{NA} for the
#'     doses not tried.}
#'
#' @details
#'   Only the doses with at least one patient are considered. A dose is
#'   admissible when the posterior probability that its DLT rate is below
#'   \code{p_eli} is at least \code{1 - cutoff_eli} and the posterior
#'   probability that its response rate is above \code{p_eli_resp} is at
#'   least \code{1 - cutoff_eli_resp}, under uniform priors, and when no
#'   lower dose is excluded by its DLT data. When no dose is admissible, no OBD
#'   is selected. Otherwise the DLT rates of the doses tried
#'   are estimated by isotonic regression of the observed rates with equal
#'   weights. Under \code{"max_resp"} the estimated MTD is the highest dose
#'   whose estimate is as close to \code{target} as that of the closest
#'   admissible dose, and the OBD is the admissible dose up to the MTD with the
#'   highest response estimate. Ties go to the lowest dose. These are the
#'   rules of \code{boinet::obd.select()}, except that in 'boinet' a dose
#'   excluded by its DLT data does not exclude the doses above it, which can
#'   then be selected; \code{boinet_compatible = TRUE} applies this rule, and
#'   \code{\link{boin_et_simulate}} explains the change.
#'
#'   The response estimates (under \code{"max_resp"}) and the utilities are
#'   compared after rounding to ten decimals, and ties go to the lowest dose.
#'   Doses whose utilities are equal as fractions then tie whatever the
#'   rounding error of the platform. This happens regularly: with
#'   \code{"utility_scoring"} and its default scores, a dose with 1 DLT and 7
#'   responses in 15 patients and a dose with 4 DLTs and 9 responses in 15
#'   patients both have the utility 196 / 3. 'boinet' compares the utilities
#'   as computed, so that rounding error decides between such doses, and
#'   versions of this package up to 2.2.0 did the same; in these trials, a
#'   few in a thousand, the selected dose can therefore differ from that of
#'   \code{boinet::obd.select()} and from that of earlier versions, also
#'   with \code{boinet_compatible = TRUE}.
#'
#'   The model-based estimates are computed in R once for every distinct set
#'   of final data; the selection itself runs in C++.
#'
#' @references
#'   Takeda, K., Taguri, M. and Morita, S. (2018). BOIN-ET: Bayesian optimal
#'   interval design for dose finding based on both efficacy and toxicity
#'   outcomes. Pharmaceutical Statistics, 17(4), 383-395.
#'   \doi{10.1002/pst.1864}
#'
#'   Yamaguchi, Y., Takeda, K., Yoshida, S. and Maruo, K. (2024). Optimal
#'   biological dose selection in dose-finding trials with model-assisted
#'   designs based on efficacy and toxicity: a simulation study. Journal of
#'   Biopharmaceutical Statistics, 34(3), 379-393.
#'   \doi{10.1080/10543406.2023.2202259}
#'
#'   Lin, R. and Yin, G. (2017). STEIN: A simple toxicity and efficacy
#'   interval design for seamless phase I/II clinical trials. Statistics in
#'   Medicine, 36(26), 4106-4120. \doi{10.1002/sim.7428}
#'
#' @examples
#' # One trial: patients, DLTs and responses at four doses
#' boin_et_select_obd(
#'   n_pts = c(3, 6, 12, 3), n_tox = c(0, 1, 3, 2), n_resp = c(0, 2, 7, 2),
#'   target = 0.30, target_resp = 0.60
#' )
#'
#' # The OBD under a utility
#' boin_et_select_obd(
#'   n_pts = c(3, 6, 12, 3), n_tox = c(0, 1, 3, 2), n_resp = c(0, 2, 7, 2),
#'   target = 0.30, target_resp = 0.60, obd_method = "utility_scoring"
#' )$obd
#'
#' @seealso \code{\link{sim_boin_et}}, \code{\link{boin_et_simulate}}
#'
#' @export
boin_et_select_obd <- function(n_pts, n_tox, n_resp, target, target_resp,
                               p_saf = NULL, p_tox = NULL, p_saf_resp = NULL,
                               cutoff_eli = 0.95, cutoff_eli_resp = 0.99,
                               p_eli = NULL, p_eli_resp = NULL,
                               est_method = c("observed", "multi_iso",
                                              "fp_logistic"),
                               obd_method = c("max_resp", "utility_weighted",
                                              "utility_truncated_linear",
                                              "utility_scoring"),
                               w1 = 0.33, w2 = 1.09, t_low = NULL,
                               t_upp = NULL, e_low = NULL, e_upp = NULL,
                               psi00 = 40, psi11 = 60, early_stop = NULL,
                               boinet_compatible = FALSE) {

  est_method <- match.arg(est_method)
  obd_method <- match.arg(obd_method)
  check_flag(boinet_compatible, "boinet_compatible")

  n_pts <- as_count_matrix(n_pts, "n_pts")
  n_tox <- as_count_matrix(n_tox, "n_tox")
  n_resp <- as_count_matrix(n_resp, "n_resp")
  if (!identical(dim(n_pts), dim(n_tox)) || !identical(dim(n_pts), dim(n_resp))) {
    stop("'n_pts', 'n_tox' and 'n_resp' must have the same dimensions",
         call. = FALSE)
  }
  if (any(n_tox > n_pts) || any(n_resp > n_pts)) {
    stop("'n_tox' and 'n_resp' must not exceed 'n_pts'", call. = FALSE)
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

  par <- boin_et_check_design(target, target_resp, p_saf, p_tox, p_saf_resp)
  check_scalar_prob(cutoff_eli, "cutoff_eli")
  check_scalar_prob(cutoff_eli_resp, "cutoff_eli_resp")
  if (is.null(t_low)) t_low <- par$p_saf
  if (is.null(t_upp)) t_upp <- par$p_tox
  if (is.null(e_low)) e_low <- par$p_saf_resp / 2
  if (is.null(e_upp)) e_upp <- par$target_resp
  for (nm in c("w1", "w2", "t_low", "t_upp", "e_low", "e_upp", "psi00",
               "psi11")) {
    value <- get(nm)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
      stop("'", nm, "' must be a single finite number", call. = FALSE)
    }
  }
  if (est_method == "fp_logistic" && !requireNamespace("mfp", quietly = TRUE)) {
    stop("Package 'mfp' is needed for est_method = \"fp_logistic\"",
         call. = FALSE)
  }

  max_n <- max(n_pts)
  eli <- boin_et_check_eli(p_eli, p_eli_resp, par)
  adm <- boin_et_admissible_table(eli$p_eli, eli$p_eli_resp, cutoff_eli,
                                  cutoff_eli_resp, max_n)

  # The model-based estimates are only needed where the rule is applied:
  # trials that did not stop early and have an admissible dose.
  tok <- matrix(adm$tox_ok[cbind(as.vector(n_tox) + 1L, as.vector(n_pts) + 1L)],
                nrow = n_trials)
  rok <- matrix(adm$resp_ok[cbind(as.vector(n_resp) + 1L, as.vector(n_pts) + 1L)],
                nrow = n_trials)
  any_adm <- rowSums((tok == 1L) & (rok == 1L) & (n_pts > 0L)) > 0L
  need <- !early_stop & any_adm

  probe <- boin_et_estimate_resp(n_pts, n_resp, est_method, need)
  probe_cpp <- probe
  probe_cpp[is.na(probe_cpp)] <- 0

  methods <- c("max_resp", "utility_weighted", "utility_truncated_linear",
               "utility_scoring")

  res <- boin_et_select_obd_cpp(
    n_pts = n_pts, n_tox = n_tox, n_resp = n_resp, probe = probe_cpp,
    early_stop = early_stop, tox_ok = as.integer(adm$tox_ok),
    resp_ok = as.integer(adm$resp_ok), max_n = as.integer(max_n),
    method = match(obd_method, methods) - 1L, target = par$target,
    tox_upper = par$p_tox, w1 = w1, w2 = w2, t_low = t_low, t_upp = t_upp,
    e_low = e_low, e_upp = e_upp, psi00 = psi00, psi11 = psi11,
    boinet_compatible = as.integer(boinet_compatible)
  )

  obd <- res$obd
  obd[obd == 0L] <- NA_integer_
  reasons <- c("selected", "early_stop", "no_admissible_dose")
  colnames(probe) <- paste0("DL", seq_len(ncol(n_pts)))

  list(obd = obd, reason = reasons[res$reason + 1L], probe = probe)
}
