#' Conditional Response Probabilities of BOIN-ET
#'
#' @description
#'   Response probabilities given the DLT for the binary link of the BOIN-ET
#'   simulation. With a correlation \code{resp_cor} between the binary DLT and
#'   response of a patient, as in Takeda et al. (2018) after Sato et al.
#'   (2016), Pr(response | no DLT) is (rho * s - p_resp * (1 - p_true)) /
#'   (p_true - 1) and Pr(response | DLT) is (rho * s + p_resp * p_true) /
#'   p_true, where s = sqrt(p_resp * (1 - p_resp) * p_true * (1 - p_true)). A
#'   DLT rate of 0 or 1 leaves one of them unused, and it is then set to
#'   \code{p_resp}. Internal.
#'
#' @param p_true,p_resp
#'   Numeric vectors. DLT and response probabilities per dose.
#'
#' @param resp_cor
#'   Numeric scalar. Correlation of the binary outcomes.
#'
#' @return
#'   A list with the vectors \code{q0} and \code{q1}.
#'
#' @keywords internal
#' @noRd
boin_et_resp_link <- function(p_true, p_resp, resp_cor) {
  s <- sqrt(p_resp * (1 - p_resp) * p_true * (1 - p_true))
  q0 <- p_resp
  q1 <- p_resp
  use0 <- p_true < 1
  use1 <- p_true > 0
  q0[use0] <- (resp_cor * s[use0] - p_resp[use0] * (1 - p_true[use0])) /
    (p_true[use0] - 1)
  q1[use1] <- (resp_cor * s[use1] + p_resp[use1] * p_true[use1]) / p_true[use1]
  bad <- q0 < -1e-12 | q0 > 1 + 1e-12 | q1 < -1e-12 | q1 > 1 + 1e-12
  if (any(bad)) {
    stop("'resp_cor' is not attainable for the DLT and response ",
         "probabilities of dose level(s) ",
         paste(which(bad), collapse = ", "), call. = FALSE)
  }
  list(q0 = pmin(pmax(q0, 0), 1), q1 = pmin(pmax(q1, 0), 1))
}
