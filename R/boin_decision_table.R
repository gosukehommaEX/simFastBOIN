#' Decision Table for the BOIN Design
#'
#' @description
#'   Build the lookup table that maps a pair of DLT and patient counts at the
#'   current dose to the dose assignment decision.
#'
#' @param target
#'   Numeric scalar. Target DLT probability.
#'
#' @param max_n
#'   Integer scalar. Largest number of patients treated at a single dose.
#'
#' @param p_saf
#'   Numeric scalar. Highest DLT probability deemed subtherapeutic.
#'   Defaults to \code{0.6 * target}.
#'
#' @param p_tox
#'   Numeric scalar. Lowest DLT probability deemed overly toxic.
#'   Defaults to \code{1.4 * target}.
#'
#' @param cutoff_eli
#'   Numeric scalar. Posterior probability cutoff for dose elimination.
#'   Defaults to 0.95.
#'
#' @return
#'   A character matrix with \code{max_n + 1} rows, labelled by the number of
#'   DLTs from 0 to \code{max_n}, and \code{max_n} columns, labelled by the number
#'   of patients from 1 to \code{max_n}. Entries are \code{"E"} to escalate,
#'   \code{"S"} to stay, \code{"D"} to de-escalate and \code{"DE"} to de-escalate
#'   and eliminate the dose together with all higher doses. Combinations with more
#'   DLTs than patients are \code{NA}.
#'
#' @details
#'   The table is derived from \code{\link{boin_boundary}} and therefore agrees
#'   with the decisions taken by the simulation engine. Elimination takes
#'   precedence over the other rules, then escalation, then de-escalation.
#'
#' @examples
#' decisions <- boin_decision_table(target = 0.30, max_n = 12)
#' decisions[1:5, 1:12]
#'
#' @seealso \code{\link{boin_boundary}}
#'
#' @export
boin_decision_table <- function(target, max_n, p_saf = NULL, p_tox = NULL,
                                cutoff_eli = 0.95) {

  bound <- boin_boundary(target, max_n, p_saf = p_saf, p_tox = p_tox,
                         cutoff_eli = cutoff_eli)
  max_n <- as.integer(max_n)

  n_tox <- 0:max_n
  out <- matrix(NA_character_, nrow = max_n + 1L, ncol = max_n,
                dimnames = list(n_tox, seq_len(max_n)))

  for (k in seq_len(max_n)) {
    y <- 0:k
    decision <- rep("S", length(y))
    escalate <- y <= bound$b_esc[k]
    deescalate <- !escalate & (y >= bound$b_deesc[k])
    decision[escalate] <- "E"
    decision[deescalate] <- "D"
    if (!is.na(bound$b_elim[k])) decision[y >= bound$b_elim[k]] <- "DE"
    out[seq_along(y), k] <- decision
  }

  out
}
