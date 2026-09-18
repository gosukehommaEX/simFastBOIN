#' Escalation and De-escalation Interval Boundaries
#'
#' @description
#'   Compute the two boundaries that define the BOIN decision interval. A dose is
#'   escalated when the observed DLT rate falls at or below \code{lambda_e} and
#'   de-escalated when it reaches \code{lambda_d}.
#'
#' @param target
#'   Numeric scalar. Target DLT probability, for example 0.30.
#'
#' @param p_saf
#'   Numeric scalar. Highest DLT probability deemed subtherapeutic, so that dose
#'   escalation should be undertaken. Defaults to \code{0.6 * target}.
#'
#' @param p_tox
#'   Numeric scalar. Lowest DLT probability deemed overly toxic, so that dose
#'   de-escalation is required. Defaults to \code{1.4 * target}.
#'
#' @return
#'   A list with components \code{lambda_e} and \code{lambda_d}.
#'
#' @details
#'   The boundaries minimize the probability of incorrect dose assignment under a
#'   three-point hypothesis on the DLT probability at the current dose. Values of
#'   \code{p_saf} and \code{p_tox} close to \code{target} should be avoided,
#'   because the sample sizes of phase I trials cannot distinguish the target rate
#'   from rates close to it.
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @examples
#' boin_lambda(target = 0.30)
#'
#' boin_lambda(target = 0.25, p_saf = 0.12, p_tox = 0.40)
#'
#' @export
boin_lambda <- function(target, p_saf = NULL, p_tox = NULL) {

  if (is.null(p_saf)) p_saf <- 0.6 * target
  if (is.null(p_tox)) p_tox <- 1.4 * target

  check_thresholds(target, p_saf, p_tox)

  lambda_e <- log((1 - p_saf) / (1 - target)) /
    log(target * (1 - p_saf) / (p_saf * (1 - target)))

  lambda_d <- log((1 - target) / (1 - p_tox)) /
    log(p_tox * (1 - target) / (target * (1 - p_tox)))

  list(lambda_e = lambda_e, lambda_d = lambda_d)
}
