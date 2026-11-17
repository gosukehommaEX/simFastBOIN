#' True OBD Combinations of a Comb-BOIN12 Scenario
#'
#' @param p_true,p_resp Numeric matrices of the true DLT and response
#'   probabilities.
#' @param target Numeric scalar. Highest acceptable DLT probability.
#' @param psi00,psi11 Numeric scalars. Utilities of neither outcome and of
#'   both.
#' @param resp_cor,resp_cor_type Correlation of the two outcomes, as in
#'   \code{comb_boin12_simulate()}.
#'
#' @return A logical matrix, \code{TRUE} for the combinations with the
#'   highest true utility among those with a DLT probability at most
#'   \code{target} (ties within \code{1e-10}), the definition of Lu et al.
#'   (2025) and of \code{oc_scenarios()} with \code{tox_limit = target}.
#'
#' @noRd
comb_boin12_true_obd <- function(p_true, p_resp, target, psi00, psi11,
                                 resp_cor, resp_cor_type) {
  pr <- oc_outcome_probs(as.vector(p_true), as.vector(p_resp), resp_cor,
                         resp_cor_type)
  utility <- 100 * pr$resp_no_dlt + psi00 * pr$neither + psi11 * pr$both
  ok <- as.vector(p_true) <= target + 1e-10
  obd <- logical(length(utility))
  if (any(ok)) obd <- ok & utility >= max(utility[ok]) - 1e-10
  matrix(obd, nrow(p_true), ncol(p_true))
}
