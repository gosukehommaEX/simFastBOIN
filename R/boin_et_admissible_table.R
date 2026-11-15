#' Admissibility Flags of BOIN-ET
#'
#' @description
#'   Whether the DLT data and the response data of a dose leave it admissible
#'   in the BOIN-ET design, for every count among 0 to \code{max_n} patients.
#'   With a uniform prior, a dose is excluded when Pr(DLT rate < p_eli) falls
#'   below \code{1 - cutoff_eli} or when Pr(response rate > p_eli_resp) falls
#'   below \code{1 - cutoff_eli_resp}. With \code{p_eli} the target and
#'   \code{p_eli_resp} the highest subtherapeutic response probability these
#'   are the rules of \code{boinet::boinet()}. Internal.
#'
#' @param p_eli,p_eli_resp
#'   Numeric scalars. DLT and response probabilities of the rules.
#'
#' @param cutoff_eli,cutoff_eli_resp
#'   Numeric scalars. Posterior probability cutoffs.
#'
#' @param max_n
#'   Integer scalar. Largest number of patients at one dose.
#'
#' @return
#'   A list of two integer matrices, \code{tox_ok} and \code{resp_ok},
#'   indexed by \code{[count + 1, patients + 1]}, which in column-major order
#'   is the layout read by the C++ engine. Counts above the number of patients
#'   are set to 0.
#'
#' @keywords internal
#' @noRd
boin_et_admissible_table <- function(p_eli, p_eli_resp, cutoff_eli,
                                     cutoff_eli_resp, max_n) {
  w <- max_n + 1L
  k <- rep(0:max_n, times = w)
  n <- rep(0:max_n, each = w)
  valid <- k <= n

  tterm <- stats::pbeta(p_eli, 1 + k[valid], 1 + n[valid] - k[valid])
  eterm <- 1 - stats::pbeta(p_eli_resp, 1 + k[valid], 1 + n[valid] - k[valid])

  tox_ok <- integer(w * w)
  resp_ok <- integer(w * w)
  tox_ok[valid] <- as.integer(!(tterm < (1 - cutoff_eli)))
  resp_ok[valid] <- as.integer(!(eterm < (1 - cutoff_eli_resp)))

  list(tox_ok = matrix(tox_ok, nrow = w), resp_ok = matrix(resp_ok, nrow = w))
}
