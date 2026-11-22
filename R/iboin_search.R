#' Integer Boundaries of the iBOIN Design by Search
#'
#' @description
#'   Internal helper. Among the integer boundaries \code{e < d} (escalate
#'   when the number of DLTs is at most \code{e}, de-escalate when it is at
#'   least \code{d}), find the pair with the smallest probability of an
#'   incorrect decision under the three point hypotheses. Used where the
#'   closed form of the boundaries gives an escalation boundary above the
#'   de-escalation boundary.
#'
#' @param n Number of patients treated at the dose.
#' @param prior Prior probabilities of the hypotheses that the DLT
#'   probability is \code{target}, \code{p_saf} and \code{p_tox}.
#' @param target Target DLT probability.
#' @param p_saf Highest DLT probability deemed subtherapeutic.
#' @param p_tox Lowest DLT probability deemed overly toxic.
#'
#' @return Integer vector of length two, \code{e} and \code{d}. Among pairs
#'   with the same probability the smallest \code{e} and then the smallest
#'   \code{d} is returned.
#'
#' @importFrom stats pbinom
#'
#' @noRd
iboin_search <- function(n, prior, target, p_saf, p_tox) {
  n <- as.integer(n)
  e <- rep(-1L:n, each = n + 2L)
  d <- rep(0L:(n + 1L), times = n + 2L)
  keep <- e < d
  e <- e[keep]
  d <- d[keep]
  err <- prior[1L] * (pbinom(e, n, target) +
                        pbinom(d - 1L, n, target, lower.tail = FALSE)) +
    prior[2L] * pbinom(e, n, p_saf, lower.tail = FALSE) +
    prior[3L] * pbinom(d - 1L, n, p_tox)
  best <- which.min(err)
  c(e[best], d[best])
}
