#' Calculate BOIN Escalation and De-escalation Boundaries
#'
#' @description
#'   Compute the escalation and de-escalation boundaries for the Bayesian Optimal
#'   Interval (BOIN) design based on the target toxicity probability. These boundaries
#'   are used to determine dose escalation and de-escalation decisions during the trial.
#'
#' @param target
#'   Numeric. The target toxicity probability (e.g., 0.30 for 30%).
#'   Should be between 0 and 1.
#'
#' @param p_saf
#'   Numeric. The toxicity probability threshold for safe (subtherapeutic) dose.
#'   Default is 0.6 * target. Doses with observed toxicity rate <= p_saf trigger escalation.
#'
#' @param p_tox
#'   Numeric. The toxicity probability threshold for overly toxic dose.
#'   Default is 1.4 * target. Doses with observed toxicity rate >= p_tox trigger de-escalation.
#'
#' @return
#'   A list containing:
#'   \item{lambda_e}{Numeric. Escalation boundary: toxicity rate threshold for escalation}
#'   \item{lambda_d}{Numeric. De-escalation boundary: toxicity rate threshold for de-escalation}
#'
#' @details
#'   Based on Bayesian theory with three-point hypothesis, these boundaries are
#'   computed to optimally balance exploration and exploitation in dose finding.
#'   The default values (p_saf = 0.6 * target, p_tox = 1.4 * target) are recommended
#'   and generally yield excellent operating characteristics.
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @examples
#' # Calculate boundaries for 30% target toxicity rate
#' boin_bound <- get_boin_boundary(target = 0.30)
#' print(boin_bound)
#' # $lambda_e
#' # [1] 0.236
#' # $lambda_d
#' # [1] 0.359
#'
#' # Calculate boundaries for 25% target with custom thresholds
#' boin_bound_custom <- get_boin_boundary(
#'   target = 0.25,
#'   p_saf = 0.12,
#'   p_tox = 0.40
#' )
#' print(boin_bound_custom)
#'
#' # Compare boundaries for different target rates
#' targets <- c(0.20, 0.25, 0.30, 0.35)
#' boundaries <- lapply(targets, get_boin_boundary)
#' names(boundaries) <- paste0("Target_", targets * 100, "%")
#' print(boundaries)
#'
#' @export
get_boin_boundary <- function(target, p_saf = NULL, p_tox = NULL) {

  # Set default thresholds if not provided
  # p_saf: safe threshold as 60% of target toxicity rate
  # p_tox: toxic threshold as 140% of target toxicity rate
  if (is.null(p_saf)) p_saf <- 0.6 * target
  if (is.null(p_tox)) p_tox <- 1.4 * target

  # Calculate escalation boundary (lambda_e)
  # This is the threshold for observed toxicity rate below which we escalate the dose
  # Derived from Bayesian likelihood ratio test between target and safe dose hypotheses
  lambda_e <- '/'(
    log((1 - p_saf) / (1 - target)),
    log(target * (1 - p_saf) / (p_saf * (1 - target)))
  )

  # Calculate de-escalation boundary (lambda_d)
  # This is the threshold for observed toxicity rate above which we de-escalate the dose
  # Derived from Bayesian likelihood ratio test between target and toxic dose hypotheses
  lambda_d <- '/'(
    log((1 - target) / (1 - p_tox)),
    log(p_tox * (1 - target) / (target * (1 - p_tox)))
  )

  # Return boundaries as a list
  return(list(
    lambda_e = lambda_e,
    lambda_d = lambda_d
  ))
}
