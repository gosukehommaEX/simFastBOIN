#' Select the OBD Combination of Comb-BOIN12 Trials
#'
#' @description
#'   Select the optimal biological dose combination (OBDC) of one or more
#'   trials of the Comb-BOIN12 design of Lu et al. (2025) from their final
#'   data. The DLT rates are made monotone along both drugs by bivariate
#'   isotonic regression, the maximum tolerated dose combination (MTDC)
#'   bounds the candidates, and the OBDC is the admissible candidate with the
#'   highest estimated utility.
#'
#' @param n_pts,n_tox,n_resp
#'   Numeric matrices for a single trial, with one row per level of drug A
#'   and one column per level of drug B, or arrays with dimensions trial by
#'   level of drug A by level of drug B, of the patients treated, the DLTs
#'   and the responses.
#'
#' @param target,target_resp,cutoff_eli,cutoff_eli_resp,extrasafe,offset
#'   As in \code{\link{comb_boin12_simulate}}. The admissibility of the
#'   combinations is judged again from the final data with the same rules.
#'
#' @param psi00,psi11
#'   Numeric scalars from 0 to 100. Utilities of a patient with neither a
#'   DLT nor a response and with both, as in \code{\link{sim_boin12}}.
#'   Default to 40 and 60.
#'
#' @param est_method
#'   Character scalar. Estimate of the response probability:
#'   \code{"logistic"} (the default) the posterior mean under the logistic
#'   model of Lu et al. (2025), see Details; \code{"observed"} the observed
#'   response rate.
#'
#' @param obd_bound
#'   Character scalar. Combinations from which the OBDC is chosen:
#'   \code{"rectangle"} (the default) those with both levels at most those
#'   of the MTDC, as Lu et al. (2025) write; \code{"contour"} those, in every
#'   column (level of drug B), at most the combination of that column whose
#'   isotonic DLT rate is closest to \code{target}, the rule that reproduces
#'   the results published by Lu et al. (2025). See Details.
#'
#' @param dose_a,dose_b
#'   Numeric vectors or \code{NULL}. Dose values of the levels of drug A and
#'   of drug B in the logistic model. Default to \code{NULL}, the level
#'   numbers. Each is standardized to mean 0 and standard deviation 0.5.
#'
#' @param logistic_scale
#'   Numeric vector of two positive numbers. Scales of the Cauchy priors of
#'   the intercept and of the other coefficients of the logistic model.
#'   Defaults to \code{c(2.5, 2.5)}, as in Lu et al. (2025).
#'
#' @param n_draws
#'   Integer scalar. Number of points of the importance sampling of the
#'   logistic model. Defaults to 4000.
#'
#' @param early_stop
#'   Logical vector or \code{NULL}. Trials stopped early, which select no
#'   OBDC. Defaults to \code{NULL}, no trial.
#'
#' @return
#'   A list with components
#'   \item{obd_a, obd_b}{Integer vectors, the levels of drug A and drug B of
#'     the OBDC of every trial (\code{NA} when none is selected).}
#'   \item{mtd_a, mtd_b}{Integer vectors, the levels of the MTDC.}
#'   \item{reason}{Character vector: \code{"selected"}, \code{"early_stop"}
#'     or \code{"no_admissible_dose"}.}
#'   \item{tox_est}{Array of the isotonic estimates of the DLT rates.}
#'   \item{resp_est}{Array of the estimates of the response probabilities of
#'     the tried combinations (\code{NA} elsewhere and for the trials with no
#'     selection).}
#'   \item{utility}{Array of the estimated utilities of the tried
#'     combinations, on the scale 0 to 100.}
#'
#' @details
#'   The DLT rates are estimated by \code{(y + 0.05) / (n + 0.1)} with
#'   weights \code{n + 0.1}, as in \code{\link{comb_boin_select_mtd}}, and
#'   made non-decreasing in both levels by the weighted least-squares fit of
#'   bivariate isotonic regression, computed by the algorithm of Dykstra and
#'   Robertson (1982) until no estimate changes by \code{1e-13} in a cycle.
#'   The MTDC is the tried combination whose estimate, plus \code{1e-5}
#'   times the sum of its two levels, is closest to \code{target}: below the
#'   target the small amount favors the higher combination, above it the
#'   lower one. Among combinations whose distances differ by at most
#'   \code{1e-8}, such as two combinations with the same sum of levels in one
#'   block of the isotonic regression, the first in column-major order is
#'   the MTDC. Untried combinations are never selected.
#'
#'   The utility of a combination is estimated by
#'   \code{100 * pE * (1 - pT) + psi00 * (1 - pT) * (1 - pE) +
#'   psi11 * pT * pE}, with \code{pT} the isotonic estimate and \code{pE} the
#'   estimate of \code{est_method}. When \code{psi00 + psi11} is 100 it is a
#'   linear function of \code{pE - wT * pT}, the estimated utility of Lu et
#'   al. (2025), with \code{wT = psi00 / psi11}. The OBDC is the admissible
#'   tried candidate with the highest estimate; among the candidates within
#'   \code{1e-6} of it, the one with the lowest sum of the levels and then
#'   the first in column-major order.
#'
#'   Lu et al. (2025) state the candidates as the combinations (j, k) with
#'   j and k at most the levels of the MTDC. Applied to the trials that
#'   reproduce their allocation percentages, that rule gives selection
#'   percentages far from those published; \code{obd_bound = "contour"}, which
#'   bounds every level of drug B by its own MTD along drug A, agrees with
#'   their Figure 2 and Supplementary Figure S2 (see
#'   \code{inst/reproduce/reproduce_Lu2025.R}). It treats the two drugs
#'   differently.
#'
#'   The logistic model is \code{logit(pE) = g0 + g1 a + g2 b + g3 a^2 +
#'   g4 b^2}, with \code{a} and \code{b} the standardized doses (mean 0 and
#'   standard deviation 0.5 over the levels of each drug, the scaling of
#'   Gelman et al., 2008, for these priors; Lu et al., 2025, do not state
#'   one) and independent Cauchy priors centered at 0. Lu et al. draw the
#'   posterior with Gibbs sampling; here its mean is computed by importance
#'   sampling with \code{n_draws} Halton points: the proposal is a
#'   multivariate t distribution with 4 degrees of freedom centered at the
#'   posterior mode, with the inverse of the negative Hessian as scale
#'   matrix, and it is fitted twice more to the weighted mean and covariance
#'   of the weighted points. No random number is used, so that the selection
#'   does not change the random number stream. In checks against long Markov
#'   chain Monte Carlo runs the estimates were within 0.01.
#'
#' @references
#'   Lu, M., Zhang, J., Yuan, Y. and Lin, R. (2025). Comb-BOIN12: A
#'   utility-based Bayesian optimal interval design for dose optimization in
#'   cancer drug-combination trials. Statistics in Biopharmaceutical
#'   Research, 17(2), 266-276. \doi{10.1080/19466315.2024.2370403}
#'
#'   Gelman, A., Jakulin, A., Pittau, M. G. and Su, Y.-S. (2008). A weakly
#'   informative default prior distribution for logistic and other
#'   regression models. The Annals of Applied Statistics, 2(4), 1360-1383.
#'   \doi{10.1214/08-AOAS191}
#'
#'   Dykstra, R. L. and Robertson, T. (1982). An algorithm for isotonic
#'   regression for two or more independent variables. The Annals of
#'   Statistics, 10(3), 708-716. \doi{10.1214/aos/1176345866}
#'
#' @examples
#' n_pts <- matrix(c(3, 6, 0, 9, 12, 3), nrow = 2)
#' n_tox <- matrix(c(0, 1, 0, 1, 2, 1), nrow = 2)
#' n_resp <- matrix(c(0, 1, 0, 3, 6, 2), nrow = 2)
#' comb_boin12_select_obd(n_pts, n_tox, n_resp, target = 0.33,
#'                        target_resp = 0.2)
#'
#' @seealso \code{\link{sim_comb_boin12}}, \code{\link{comb_boin12_simulate}}
#'
#' @export
comb_boin12_select_obd <- function(n_pts, n_tox, n_resp, target, target_resp,
                                   cutoff_eli = 0.95, cutoff_eli_resp = 0.90,
                                   extrasafe = FALSE, offset = 0.05,
                                   psi00 = 40, psi11 = 60,
                                   est_method = c("logistic", "observed"),
                                   obd_bound = c("rectangle", "contour"),
                                   dose_a = NULL, dose_b = NULL,
                                   logistic_scale = c(2.5, 2.5),
                                   n_draws = 4000, early_stop = NULL) {

  est_method <- match.arg(est_method)
  obd_bound <- match.arg(obd_bound)
  boin12_check_design(target, target_resp, psi00 = psi00, psi11 = psi11)
  check_scalar_prob(cutoff_eli, "cutoff_eli")
  check_scalar_prob(cutoff_eli_resp, "cutoff_eli_resp")
  check_flag(extrasafe, "extrasafe")
  check_comb_offset(offset, cutoff_eli)
  if (!is.numeric(logistic_scale) || length(logistic_scale) != 2L ||
      any(!is.finite(logistic_scale)) || any(logistic_scale <= 0)) {
    stop("'logistic_scale' must hold two positive numbers", call. = FALSE)
  }
  check_count(n_draws, "n_draws", 100L)

  n_pts <- as_comb_count_array(n_pts, "n_pts")
  n_tox <- as_comb_count_array(n_tox, "n_tox")
  n_resp <- as_comb_count_array(n_resp, "n_resp")
  if (!identical(dim(n_pts), dim(n_tox)) ||
      !identical(dim(n_pts), dim(n_resp))) {
    stop("'n_pts', 'n_tox' and 'n_resp' must have the same dimensions",
         call. = FALSE)
  }
  if (any(n_tox > n_pts) || any(n_resp > n_pts)) {
    stop("'n_tox' and 'n_resp' must not exceed 'n_pts' at any combination",
         call. = FALSE)
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
  code_a <- comb_boin12_dose_codes(dose_a, n_rows, "dose_a")
  code_b <- comb_boin12_dose_codes(dose_b, n_cols, "dose_b")

  max_n <- max(n_pts)
  adm <- comb_boin12_tables(target, target_resp, cutoff_eli, cutoff_eli_resp,
                            extrasafe, offset, max_n)
  points <- comb_boin12_points(n_draws, df = 4)

  res <- comb_boin12_select_cpp(
    n_pts = n_pts, n_tox = n_tox, n_resp = n_resp,
    early_stop = early_stop, target = as.numeric(target),
    psi00 = as.numeric(psi00), psi11 = as.numeric(psi11),
    max_n = as.integer(max_n), tox_ok = adm$tox_ok,
    tox_ok_lowest = adm$tox_ok_lowest, resp_ok = adm$resp_ok,
    est = if (est_method == "logistic") 0L else 1L,
    bound = if (obd_bound == "contour") 1L else 0L,
    fa = rep(code_a, times = n_cols), fb = rep(code_b, each = n_rows),
    scale0 = as.numeric(logistic_scale[1L]),
    scale1 = as.numeric(logistic_scale[2L]),
    z = points$z, w = points$w, df = 4, rounds = 2L
  )

  if (res$not_converged > 0L) {
    warning("the isotonic estimates of ", res$not_converged, " trial(s) ",
            "did not converge", call. = FALSE)
  }
  dimnames(res$utility) <- dimnames(n_pts)
  dimnames(res$resp_est) <- dimnames(n_pts)
  dimnames(res$tox_est) <- dimnames(n_pts)
  list(
    obd_a = res$obd_a,
    obd_b = res$obd_b,
    mtd_a = res$mtd_a,
    mtd_b = res$mtd_b,
    reason = c("selected", "early_stop", "no_admissible_dose")[
      res$reason + 1L],
    tox_est = res$tox_est,
    resp_est = res$resp_est,
    utility = res$utility
  )
}
