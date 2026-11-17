#' Probabilities of the Four Outcomes of a Patient
#'
#' @description
#'   Internal helper of \code{\link{oc_scenarios}}. The probabilities of a
#'   response without a DLT, of neither, of both and of a DLT without a
#'   response at each dose, for the two models of the correlation between the
#'   DLT and the response used by the simulation functions.
#'
#' @param p_true Numeric vector of DLT probabilities.
#' @param p_resp Numeric vector of response probabilities of the same length.
#' @param resp_cor Numeric scalar, the correlation.
#' @param resp_cor_type Character scalar, \code{"copula"} or \code{"binary"}.
#'
#' @return A list with the numeric vectors \code{resp_no_dlt},
#'   \code{neither}, \code{both} and \code{dlt_no_resp}.
#'
#' @noRd
oc_outcome_probs <- function(p_true, p_resp, resp_cor, resp_cor_type) {
  if (resp_cor_type == "binary") {
    # Pr(response | DLT) = q1 and Pr(response | no DLT) = q0
    link <- boin_et_resp_link(p_true, p_resp, resp_cor)
    both <- p_true * link$q1
    resp_no_dlt <- (1 - p_true) * link$q0
  } else {
    both <- vapply(seq_along(p_true), function(j) {
      oc_bivariate_normal(p_true[j], p_resp[j], resp_cor)
    }, numeric(1))
    resp_no_dlt <- p_resp - both
  }
  list(resp_no_dlt = resp_no_dlt, neither = 1 - p_true - resp_no_dlt,
       both = both, dlt_no_resp = p_true - both)
}
