#' Select the Optimal Biological Dose of gBOIN-ET Trials
#'
#' @description
#'   Select the optimal biological dose (OBD) from the final data of one or
#'   more trials of the gBOIN-ET design: the patients and the total
#'   normalized toxicity and efficacy scores at every dose.
#'
#' @inheritParams boin_et_select_obd
#' @inheritParams gboin_et_simulate
#'
#' @param n_tox
#'   Numeric vector or matrix of the same shape as \code{n_pts}. Total of the
#'   normalized toxicity scores of the patients at each dose.
#'
#' @param n_resp
#'   Numeric vector or matrix of the same shape as \code{n_pts}. Total of the
#'   normalized efficacy scores of the patients at each dose.
#'
#' @param est_method
#'   Character scalar. Estimator of the mean efficacy scores:
#'   \code{"observed"} (the default) for the observed means,
#'   \code{"fp_linear"} for the linear fractional polynomial of Takeda et al.
#'   (2022), or \code{"fp_logistic"} for the fractional polynomial logistic
#'   regression of the 'boinet' package, which needs the 'mfp' package.
#'
#' @param min_obd_sample
#'   Integer scalar. Only doses given to at least this many patients can be
#'   selected. Defaults to 1. The simulations of Takeda et al. (2022) use 6.
#'
#' @return
#'   A list with components
#'   \item{obd}{Integer vector of the selected dose levels, \code{NA} when no
#'     dose is selected.}
#'   \item{reason}{Character vector: \code{"selected"}, \code{"early_stop"} or
#'     \code{"no_admissible_dose"}.}
#'   \item{probe}{Numeric matrix of the estimates of the mean efficacy
#'     scores, \code{NA} for the doses not tried.}
#'
#' @details
#'   The rules are those of \code{\link{boin_et_select_obd}} with the mean
#'   scores in the place of the observed rates. A dose is admissible when the
#'   posterior probability that its mean toxicity score is below \code{p_eli}
#'   is at least \code{1 - cutoff_eli}, the posterior probability that its
#'   mean efficacy score is above \code{p_eli_resp} is at least
#'   \code{1 - cutoff_eli_resp}, and no lower dose is excluded by its
#'   toxicity data. The mean toxicity scores of the doses tried are estimated
#'   by isotonic regression with equal weights. Under \code{"max_resp"} the
#'   estimated MTD is the highest dose whose estimate is as close to
#'   \code{target} as that of the closest admissible dose, and the OBD is the
#'   admissible dose up to the MTD with the highest estimated efficacy score;
#'   ties go to the lowest dose. A dose with fewer than \code{min_obd_sample}
#'   patients takes part in finding the MTD but cannot be selected; when no
#'   dose can be selected the reason is \code{"no_admissible_dose"}.
#'
#'   The estimated efficacy scores (under \code{"max_resp"}) and the
#'   utilities are compared after rounding to ten decimals, and ties go to
#'   the lowest dose. Doses whose scores or utilities are equal as fractions
#'   then tie whatever the rounding error of the platform; this happens
#'   regularly with \code{"utility_scoring"}, whose default scores give
#'   different doses the same utility. \code{\link{boin_et_select_obd}}
#'   follows 'boinet' and compares the unrounded values, so that for whole
#'   numbers of DLTs and responses the two functions select the same dose
#'   under \code{"max_resp"} and can differ in such ties under a utility.
#'   The estimates in \code{probe} are rounded to ten decimals.
#'
#'   \code{est_method = "fp_linear"} fits, by least squares to the scores of
#'   the patients, the 36 models
#'   \code{b0 + b1 * dose^p1 + b2 * dose^p2} with the powers -2, -1, -0.5, 0
#'   (the logarithm), 0.5, 1, 2 and 3, where a repeated power \code{p} stands
#'   for \code{dose^p} and \code{dose^p * log(dose)} and the doses are the
#'   dose levels 1, 2 and so on, and takes the fitted values of the model
#'   with the smallest residual sum of squares. This is the estimator of
#'   Section 2.3 of Takeda et al. (2022), as their program computes it (a
#'   normal model with identity link). With three doses or fewer every model
#'   passes through the observed means, which are then used. The fit runs in
#'   C++.
#'
#'   The program of Takeda et al. (2022) differs from this function in two
#'   further points that \code{min_obd_sample = 6} and
#'   \code{boinet_compatible = TRUE} bring close: it selects only among
#'   doses with at least six patients, and it excludes no dose at the end of
#'   the trial.
#'
#' @references
#'   Takeda, K., Morita, S. and Taguri, M. (2022). gBOIN-ET: The generalized
#'   Bayesian optimal interval design for optimal dose-finding accounting for
#'   ordinal graded efficacy and toxicity in early clinical trials.
#'   Biometrical Journal, 64(7), 1178-1191. \doi{10.1002/bimj.202100263}
#'
#' @examples
#' # One trial: patients and total normalized scores at four doses
#' gboin_et_select_obd(
#'   n_pts = c(3, 6, 12, 9), n_tox = c(0, 1, 3.5, 4),
#'   n_resp = c(0.5, 2, 8, 6.5),
#'   target = 0.47 / 1.5, target_resp = 0.95 / 1.5
#' )
#'
#' # The efficacy scores estimated by the linear fractional polynomial
#' gboin_et_select_obd(
#'   n_pts = c(3, 6, 12, 9), n_tox = c(0, 1, 3.5, 4),
#'   n_resp = c(0.5, 2, 8, 6.5),
#'   target = 0.47 / 1.5, target_resp = 0.95 / 1.5, est_method = "fp_linear"
#' )$probe
#'
#' @seealso \code{\link{sim_gboin_et}}, \code{\link{gboin_et_simulate}},
#'   \code{\link{boin_et_select_obd}}
#'
#' @export
gboin_et_select_obd <- function(n_pts, n_tox, n_resp, target, target_resp,
                                p_saf = NULL, p_tox = NULL, p_saf_resp = NULL,
                                cutoff_eli = 0.95, cutoff_eli_resp = 0.99,
                                p_eli = NULL, p_eli_resp = NULL,
                                est_method = c("observed", "fp_linear",
                                               "fp_logistic"),
                                obd_method = c("max_resp", "utility_weighted",
                                               "utility_truncated_linear",
                                               "utility_scoring"),
                                min_obd_sample = 1, w1 = 0.33, w2 = 1.09,
                                t_low = NULL, t_upp = NULL, e_low = NULL,
                                e_upp = NULL, psi00 = 40, psi11 = 60,
                                early_stop = NULL,
                                boinet_compatible = FALSE) {

  est_method <- match.arg(est_method)
  obd_method <- match.arg(obd_method)
  check_flag(boinet_compatible, "boinet_compatible")
  check_count(min_obd_sample, "min_obd_sample", 1L)

  n_pts <- as_count_matrix(n_pts, "n_pts")
  n_tox <- as_score_matrix(n_tox, "n_tox")
  n_resp <- as_score_matrix(n_resp, "n_resp")
  if (!identical(dim(n_pts), dim(n_tox)) ||
      !identical(dim(n_pts), dim(n_resp))) {
    stop("'n_pts', 'n_tox' and 'n_resp' must have the same dimensions",
         call. = FALSE)
  }
  if (any(n_tox > n_pts + 1e-8) || any(n_resp > n_pts + 1e-8)) {
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
  eli <- boin_et_check_eli(p_eli, p_eli_resp, par)

  # The logistic fractional polynomial is fitted in R, once for every
  # distinct set of final data of the trials that did not stop early.
  probe <- boin_et_estimate_resp(
    n_pts, n_resp,
    if (est_method == "fp_logistic") "fp_logistic" else "observed",
    need = !early_stop
  )
  probe_cpp <- probe
  probe_cpp[is.na(probe_cpp)] <- 0

  methods <- c("max_resp", "utility_weighted", "utility_truncated_linear",
               "utility_scoring")

  res <- gboin_et_select_obd_cpp(
    n_pts = n_pts, tox_sum = n_tox, resp_sum = n_resp, probe = probe_cpp,
    early_stop = early_stop, p_eli = as.numeric(eli$p_eli),
    p_eli_resp = as.numeric(eli$p_eli_resp),
    cutoff_eli = as.numeric(cutoff_eli),
    cutoff_eli_resp = as.numeric(cutoff_eli_resp),
    min_obd_sample = as.integer(min_obd_sample),
    estimate = if (est_method == "fp_linear") 1L else 0L,
    method = match(obd_method, methods) - 1L, target = par$target,
    tox_upper = par$p_tox, w1 = w1, w2 = w2, t_low = t_low, t_upp = t_upp,
    e_low = e_low, e_upp = e_upp, psi00 = psi00, psi11 = psi11,
    boinet_compatible = as.integer(boinet_compatible)
  )

  obd <- res$obd
  obd[obd == 0L] <- NA_integer_
  reasons <- c("selected", "early_stop", "no_admissible_dose")
  probe <- res$probe
  probe[n_pts == 0L] <- NA_real_
  colnames(probe) <- paste0("DL", seq_len(ncol(n_pts)))

  list(obd = obd, reason = reasons[res$reason + 1L], probe = probe)
}
