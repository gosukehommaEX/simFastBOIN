#' Check the Design Parameters of BOIN-ET
#'
#' @description
#'   Validate the target probabilities of the BOIN-ET design and fill in the
#'   defaults of Takeda et al. (2018) for the subtherapeutic and toxic levels.
#'   Internal.
#'
#' @param target,target_resp
#'   Numeric scalars. Target DLT and response probabilities.
#'
#' @param p_saf,p_tox,p_saf_resp
#'   Numeric scalars or \code{NULL} for the defaults \code{0.1 * target},
#'   \code{1.4 * target} and \code{0.6 * target_resp}.
#'
#' @return
#'   A list with the five probabilities.
#'
#' @keywords internal
#' @noRd
boin_et_check_design <- function(target, target_resp, p_saf = NULL,
                                 p_tox = NULL, p_saf_resp = NULL) {
  check_scalar_prob(target, "target")
  check_scalar_prob(target_resp, "target_resp")
  if (is.null(p_saf)) p_saf <- target * 0.1
  if (is.null(p_tox)) p_tox <- target * 1.4
  if (is.null(p_saf_resp)) p_saf_resp <- target_resp * 0.6
  check_scalar_prob(p_saf, "p_saf")
  check_scalar_prob(p_tox, "p_tox")
  check_scalar_prob(p_saf_resp, "p_saf_resp")
  if (!(p_saf < target && target < p_tox)) {
    stop("the design parameters must satisfy p_saf < target < p_tox",
         call. = FALSE)
  }
  if (!(p_saf_resp < target_resp)) {
    stop("the design parameters must satisfy p_saf_resp < target_resp",
         call. = FALSE)
  }
  list(target = target, target_resp = target_resp, p_saf = p_saf,
       p_tox = p_tox, p_saf_resp = p_saf_resp)
}
