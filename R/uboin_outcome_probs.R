#' Outcome Probabilities of U-BOIN from the Marginal Rates
#'
#' @description
#'   Joint probabilities of the four outcomes of a patient from the marginal
#'   DLT and response probabilities of every dose, with the Gumbel model of
#'   Zhou, Lee and Yuan (2019). Internal.
#'
#' @param p_true
#'   Numeric vector. DLT probability of every dose.
#'
#' @param p_resp
#'   Numeric vector. Response probability of every dose.
#'
#' @param resp_assoc
#'   Numeric scalar. Association parameter of the Gumbel model; 0 makes the
#'   two outcomes independent.
#'
#' @return
#'   A numeric matrix with one row per dose and one column per outcome, in the
#'   order (no response, DLT), (no response, no DLT), (response, DLT),
#'   (response, no DLT).
#'
#' @details
#'   The probability of (Y_E, Y_T) = (e, t) at a dose is
#'   \code{pe^e (1 - pe)^(1 - e) pt^t (1 - pt)^(1 - t)} \code{+ pe (1 - pe)
#'   pt (1 - pt) (-1)^(e + t) (exp(c) - 1) / (exp(c) + 1)}.
#'
#' @keywords internal
#' @noRd
uboin_outcome_probs <- function(p_true, p_resp, resp_assoc) {
  if (!is.numeric(resp_assoc) || length(resp_assoc) != 1L ||
      !is.finite(resp_assoc)) {
    stop("'resp_assoc' must be a single finite number", call. = FALSE)
  }
  shrink <- (exp(resp_assoc) - 1) / (exp(resp_assoc) + 1)
  cells <- list(c(0, 1), c(0, 0), c(1, 1), c(1, 0))
  out <- vapply(cells, function(et) {
    e <- et[1L]
    t <- et[2L]
    p_resp^e * (1 - p_resp)^(1 - e) * p_true^t * (1 - p_true)^(1 - t) +
      p_resp * (1 - p_resp) * p_true * (1 - p_true) * (-1)^(e + t) * shrink
  }, numeric(length(p_true)))
  out <- matrix(out, nrow = length(p_true))
  if (any(out < -1e-12)) {
    stop("'resp_assoc' gives negative outcome probabilities for these ",
         "DLT and response probabilities", call. = FALSE)
  }
  out[out < 0] <- 0
  colnames(out) <- uboin_outcome_names()
  out
}
