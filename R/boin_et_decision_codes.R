#' Decision Codes of BOIN-ET
#'
#' @description
#'   Decision of the BOIN-ET design for every count of DLTs and responses
#'   among up to \code{max_n} patients at the current dose: 0 escalate,
#'   1 stay, 2 de-escalate, 3 choose among the neighboring doses (rule 4).
#'   The observed rates are compared with the boundaries exactly as in
#'   \code{boinet::boinet()}. Internal.
#'
#' @param lambda1,lambda2,eta1
#'   Numeric scalars. The boundaries.
#'
#' @param max_n
#'   Integer scalar. Largest number of patients at one dose.
#'
#' @return
#'   An integer array indexed by \code{[responses + 1, DLTs + 1, patients]},
#'   which in column-major order is the layout read by the C++ engine. Entries
#'   with more DLTs or responses than patients are meaningless.
#'
#' @keywords internal
#' @noRd
boin_et_decision_codes <- function(lambda1, lambda2, eta1, max_n) {
  w <- max_n + 1L
  x <- rep(0:max_n, times = w * max_n)
  y <- rep(rep(0:max_n, each = w), times = max_n)
  n <- rep(seq_len(max_n), each = w * w)
  pt <- y / n
  pe <- x / n

  code <- rep(3L, length(n))
  deescalate <- pt >= lambda2
  stay <- (pt < lambda2) & (pe > eta1)
  escalate <- (pt <= lambda1) & (pe <= eta1)
  code[deescalate] <- 2L
  code[stay] <- 1L
  code[escalate] <- 0L

  array(code, dim = c(w, w, max_n))
}
