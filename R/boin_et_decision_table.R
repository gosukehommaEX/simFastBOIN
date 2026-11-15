#' Decision Table for the BOIN-ET Design
#'
#' @description
#'   Tabulate the decision of the BOIN-ET design for every number of DLTs and
#'   responses at the current dose, for the given numbers of patients, in the
#'   form of Table 2 of Takeda, Taguri and Morita (2018). The table also marks
#'   the counts for which the dose is no longer admissible.
#'
#' @inheritParams boin_et_boundary
#'
#' @param n
#'   Integer vector. Numbers of patients at the current dose to tabulate.
#'   Defaults to \code{c(3, 6, 9, 12)}.
#'
#' @param cutoff_eli
#'   Numeric scalar. A dose is no longer admissible when the posterior
#'   probability that its DLT rate is below \code{p_eli} falls below
#'   \code{1 - cutoff_eli}. Defaults to 0.95.
#'
#' @param cutoff_eli_resp
#'   Numeric scalar. A dose is no longer admissible when the posterior
#'   probability that its response rate is above \code{p_eli_resp} falls
#'   below \code{1 - cutoff_eli_resp}. Defaults to 0.99.
#'
#' @param p_eli
#'   Numeric scalar or \code{NULL}. DLT probability of the rule above.
#'   Defaults to \code{NULL}, which uses \code{target} as the 'boinet'
#'   package does. Takeda et al. (2018) used 0.40 with a cutoff of 0.90.
#'
#' @param p_eli_resp
#'   Numeric scalar or \code{NULL}. Response probability of the rule above.
#'   Defaults to \code{NULL}, which uses \code{p_saf_resp} as the 'boinet'
#'   package does. Takeda et al. (2018) used 0.20 with a cutoff of 0.90.
#'
#' @return
#'   An object of class \code{boin_et_decision_table}, which is a data frame
#'   with one row per combination of \code{n} patients, \code{n_tox} DLTs and
#'   \code{n_resp} responses, and the columns
#'   \item{n, n_tox, n_resp}{The counts.}
#'   \item{decision}{\code{"E"} to escalate, \code{"S"} to stay, \code{"D"} to
#'     de-escalate and \code{"E/S/D"} when the next dose is chosen among the
#'     current dose and its neighbors by their observed response rates.}
#'   \item{tox_excluded, resp_excluded}{Whether the DLT data or the response
#'     data make the dose inadmissible.}
#'   The boundaries are attached as the attribute \code{"boundary"}, an object
#'   of class \code{\link{boin_et_boundary}}.
#'
#' @details
#'   When the observed DLT rate is at most \code{lambda1} and the observed
#'   response rate at most \code{eta1}, the dose is escalated; when the DLT
#'   rate is below \code{lambda2} and the response rate above \code{eta1}, the
#'   next cohort stays; when the DLT rate is at least \code{lambda2}, the dose
#'   is de-escalated. Otherwise the dose is escalated when the next dose has not
#'   been tried, and the next cohort goes to the dose with the highest observed
#'   response rate among the current dose and its neighbors, with ties broken
#'   at random. An inadmissible dose is skipped by
#'   \code{\link{boin_et_simulate}}.
#'
#' @references
#'   Takeda, K., Taguri, M. and Morita, S. (2018). BOIN-ET: Bayesian optimal
#'   interval design for dose finding based on both efficacy and toxicity
#'   outcomes. Pharmaceutical Statistics, 17(4), 383-395.
#'   \doi{10.1002/pst.1864}
#'
#' @examples
#' # Table 2 of Takeda et al. (2018)
#' tab <- boin_et_decision_table(target = 0.30, target_resp = 0.60, n = c(3, 6))
#' tab
#'
#' head(as.data.frame(tab))
#'
#' @seealso \code{\link{boin_et_boundary}}, \code{\link{boin_et_simulate}}
#'
#' @export
boin_et_decision_table <- function(target, target_resp, n = c(3, 6, 9, 12),
                                   p_saf = NULL, p_tox = NULL,
                                   p_saf_resp = NULL, cutoff_eli = 0.95,
                                   cutoff_eli_resp = 0.99, p_eli = NULL,
                                   p_eli_resp = NULL,
                                   prior = rep(1 / 6, 6)) {

  if (!is.numeric(n) || length(n) < 1L || any(!is.finite(n)) ||
      any(n != round(n)) || any(n < 1)) {
    stop("'n' must contain whole numbers of at least 1", call. = FALSE)
  }
  n <- sort(unique(as.integer(n)))
  check_scalar_prob(cutoff_eli, "cutoff_eli")
  check_scalar_prob(cutoff_eli_resp, "cutoff_eli_resp")

  bd <- boin_et_boundary(target, target_resp, p_saf = p_saf, p_tox = p_tox,
                         p_saf_resp = p_saf_resp, prior = prior)
  max_n <- max(n)
  codes <- boin_et_decision_codes(bd$lambda1, bd$lambda2, bd$eta1, max_n)
  eli <- boin_et_check_eli(p_eli, p_eli_resp, bd)
  adm <- boin_et_admissible_table(eli$p_eli, eli$p_eli_resp, cutoff_eli,
                                  cutoff_eli_resp, max_n)

  labels <- c("E", "S", "D", "E/S/D")
  blocks <- lapply(n, function(k) {
    grid <- expand.grid(n_resp = 0:k, n_tox = 0:k)
    data.frame(
      n = k,
      n_tox = grid$n_tox,
      n_resp = grid$n_resp,
      decision = labels[codes[cbind(grid$n_resp + 1L, grid$n_tox + 1L, k)] + 1L],
      tox_excluded = adm$tox_ok[cbind(grid$n_tox + 1L, k + 1L)] == 0L,
      resp_excluded = adm$resp_ok[cbind(grid$n_resp + 1L, k + 1L)] == 0L,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, blocks)
  out <- out[order(out$n, out$n_tox, out$n_resp), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "boundary") <- bd
  attr(out, "cutoff_eli") <- cutoff_eli
  attr(out, "cutoff_eli_resp") <- cutoff_eli_resp
  attr(out, "p_eli") <- eli$p_eli
  attr(out, "p_eli_resp") <- eli$p_eli_resp
  class(out) <- c("boin_et_decision_table", "data.frame")
  out
}
