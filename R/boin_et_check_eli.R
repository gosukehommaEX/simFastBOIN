#' Check the Probabilities of the Elimination Rules of BOIN-ET
#'
#' @description
#'   Fill in and validate the DLT and response probabilities used by the rules
#'   that make a dose inadmissible. Internal.
#'
#' @param p_eli,p_eli_resp
#'   Numeric scalars or \code{NULL} for \code{target} and \code{p_saf_resp}.
#'
#' @param par
#'   List with the components \code{target} and \code{p_saf_resp}, such as a
#'   \code{boin_et_boundary} object.
#'
#' @return
#'   A list with the components \code{p_eli} and \code{p_eli_resp}.
#'
#' @keywords internal
#' @noRd
boin_et_check_eli <- function(p_eli, p_eli_resp, par) {
  if (is.null(p_eli)) p_eli <- par$target
  if (is.null(p_eli_resp)) p_eli_resp <- par$p_saf_resp
  check_scalar_prob(p_eli, "p_eli")
  check_scalar_prob(p_eli_resp, "p_eli_resp")
  list(p_eli = p_eli, p_eli_resp = p_eli_resp)
}
