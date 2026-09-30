#' Decision Table for the TITE-BOIN Design
#'
#' @description
#'   Tabulate the dose assignment rule of the time-to-event Bayesian optimal
#'   interval (TITE-BOIN) design, which allows new patients to be treated while
#'   the DLT assessment of earlier patients is still pending. For every
#'   attainable combination of patients treated, DLTs observed and patients
#'   pending at the current dose, the table gives the decision or, when the
#'   decision depends on how long the pending patients have been followed, the
#'   values of the follow-up statistic at which it changes.
#'
#' @param target
#'   Numeric scalar. Target DLT probability, for example 0.30.
#'
#' @param max_n
#'   Integer scalar. Largest number of patients treated at a single dose.
#'
#' @param method
#'   Character scalar. How the pending patients enter the estimate of the DLT
#'   probability: \code{"imputation"} (the default) for the single mean
#'   imputation of Yuan et al. (2018), or \code{"ess"} for the effective sample
#'   size of Lin and Yuan (2020). See the details.
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
#' @param max_pending_ratio
#'   Numeric scalar greater than 0 and at most 1. Accrual is suspended when the
#'   proportion of the patients at the current dose whose assessment is pending
#'   exceeds this value, unless the dose is de-escalated whatever their outcomes.
#'   Defaults to 0.5 for \code{"imputation"}, following Yuan et al. (2018), and to
#'   1, which never suspends on this ground, for \code{"ess"}.
#'
#' @param min_completed
#'   Integer scalar. While some patients are pending, escalation requires at
#'   least this many patients at the current dose to have completed the
#'   assessment, and accrual is suspended instead when fewer have. Defaults to 0
#'   for \code{"imputation"} and to 2 for \code{"ess"}, following Lin and Yuan
#'   (2020).
#'
#' @return
#'   An object of class \code{tite_boin_decision_table}, which is a data frame
#'   with one row per attainable state of the current dose, ordered by
#'   \code{n}, \code{n_tox} and \code{n_pending}, and with the columns
#'   \item{n}{Number of patients treated at the current dose, from 1 to
#'     \code{max_n}.}
#'   \item{n_tox}{Number of DLTs observed, which can only come from patients who
#'     completed the assessment.}
#'   \item{n_pending}{Number of patients whose assessment is pending.}
#'   \item{decision}{Decision code. \code{"E"} escalates, \code{"S"} stays,
#'     \code{"D"} de-escalates, \code{"DE"} de-escalates and eliminates the dose
#'     together with all higher doses, and \code{"SUS"} suspends accrual until
#'     more data are available. These do not depend on the follow-up of the
#'     pending patients. The codes \code{"E/S"}, \code{"S/D"} and
#'     \code{"E/S/D"} mark states whose decision depends on the follow-up
#'     statistic through \code{esc_bound} and \code{deesc_bound}; in them
#'     \code{"SUS"} takes the place of \code{"E"} when escalation is blocked by
#'     \code{min_completed}.}
#'   \item{esc_bound}{Escalate, or suspend when escalation is blocked, when the
#'     follow-up statistic is at least this value. \code{NA} when no such
#'     boundary applies.}
#'   \item{deesc_bound}{De-escalate when the follow-up statistic is at most this
#'     value. \code{NA} when no such boundary applies.}
#'   The design parameters are stored as attributes: \code{method},
#'   \code{statistic} (\code{"STFT"} or \code{"ESS"}), \code{target},
#'   \code{p_saf}, \code{p_tox}, \code{lambda_e}, \code{lambda_d},
#'   \code{cutoff_eli}, \code{max_pending_ratio} and \code{min_completed}. The
#'   object has \code{print} and \code{plot} methods.
#'
#' @details
#'   Suppose that n patients have been treated at the current dose, that r of them
#'   have completed the DLT assessment with \code{n_tox} DLTs, and that
#'   c = n - r are pending. The standardized total follow-up time (STFT) is the
#'   sum of the follow-up times of the pending patients divided by the length of
#'   the assessment window, so that STFT lies between 0 and c. With no pending
#'   patient the decision is that of the BOIN design, taken from
#'   \code{\link{boin_boundary}}.
#'
#'   \code{method = "imputation"} replaces each pending outcome by its expected
#'   value under a uniform time to DLT. The estimate of the DLT probability is
#'   \code{(n_tox + q * (c - STFT)) / n}, where \code{q = p / (1 - p)} and
#'   \code{p = (n_tox + a) / (r + 1)} is the posterior mean under a
#'   Beta(a, 1 - a) prior with \code{a = target / 2} (Yuan et al., 2018,
#'   Supplementary Appendix A). The estimate decreases as STFT grows, so the rule
#'   becomes a pair of boundaries on STFT. As in the BOIN design, escalation is
#'   possible only while \code{n_tox / n} is below the target and de-escalation
#'   only once it reaches the target.
#'
#'   \code{method = "ess"} uses the approximated likelihood of Lin and Yuan
#'   (2020), under which the estimate is \code{n_tox / ESS} with the effective
#'   sample size ESS = r + STFT, which lies between r and n. The dose is escalated
#'   when ESS is at least \code{n_tox / lambda_e} and de-escalated when it is at
#'   most \code{n_tox / lambda_d}. The boundaries in the table refer to ESS.
#'
#'   The rules are applied in this order: elimination; de-escalation that holds
#'   whatever the outcomes of the pending patients, that is when
#'   \code{n_tox / n} is at least \code{lambda_d}; suspension for too many
#'   pending patients; and finally the boundaries on the follow-up statistic, in
#'   which an escalation blocked by \code{min_completed} becomes a suspension.
#'   Elimination requires \code{Pr(p > target | data) > cutoff_eli} computed
#'   with all n treated patients, the pending ones counted as without DLT, which
#'   is how both articles define it. It is not evaluated before three patients
#'   have been treated.
#'
#'   Neither the length of the assessment window nor the weights given to the
#'   follow-up times change the table, because both enter only through STFT. A
#'   weighted STFT, such as that of a piecewise uniform time to DLT, is compared
#'   with the same boundaries.
#'
#' @section Agreement with the published tables:
#'   With the defaults, \code{method = "imputation"} reproduces every entry of
#'   Table 1 (target 0.2) of Yuan et al. (2018) and of Table S1 (target 0.3) of
#'   its supplementary appendix. These tables allow de-escalation when
#'   \code{n_tox / n} equals the target (Table 1, 15 patients with 3 DLTs), which
#'   equation (5) of the appendix, read literally, would not. The implementation
#'   follows the tables, which agree with the verbal description of the rule in
#'   the appendix.
#'
#'   For \code{method = "ess"} no table for the BOIN design has been published;
#'   the tables of Lin and Yuan (2020) are for the keyboard and mTPI designs.
#'
#' @references
#'   Yuan, Y., Lin, R., Li, D., Nie, L. and Warren, K. E. (2018). Time-to-Event
#'   Bayesian Optimal Interval Design to Accelerate Phase I Trials. Clinical
#'   Cancer Research, 24(20), 4921-4930.
#'
#'   Lin, R. and Yuan, Y. (2020). Time-to-Event Model-Assisted Designs for
#'   Dose-Finding Trials with Delayed Toxicity. Biostatistics, 21(4), 807-824.
#'
#' @examples
#' # Table 1 of Yuan et al. (2018)
#' decisions <- tite_boin_decision_table(target = 0.2, max_n = 15)
#' print(decisions, cohort_size = 3)
#'
#' # Nine patients, one DLT and four pending: escalate once STFT reaches 2.15
#' decisions[decisions$n == 9 & decisions$n_tox == 1 & decisions$n_pending == 4, ]
#'
#' # The effective sample size of Lin and Yuan (2020)
#' decisions_ess <- tite_boin_decision_table(target = 0.3, max_n = 12, method = "ess")
#' print(decisions_ess, cohort_size = 3)
#'
#' @seealso \code{\link{print.tite_boin_decision_table}},
#'   \code{\link{plot.tite_boin_decision_table}},
#'   \code{\link{boin_decision_table}}
#'
#' @export
tite_boin_decision_table <- function(target, max_n, method = c("imputation", "ess"),
                                     p_saf = NULL, p_tox = NULL, cutoff_eli = 0.95,
                                     max_pending_ratio = NULL, min_completed = NULL) {

  method <- match.arg(method)
  if (is.null(p_saf)) p_saf <- 0.6 * target
  if (is.null(p_tox)) p_tox <- 1.4 * target
  rules <- tite_rules(method, max_pending_ratio, min_completed)
  max_pending_ratio <- rules$max_pending_ratio
  min_completed <- rules$min_completed

  bound <- boin_boundary(target, max_n, p_saf = p_saf, p_tox = p_tox,
                         cutoff_eli = cutoff_eli)
  max_n <- as.integer(max_n)

  # Every attainable state: n_tox + n_pending cannot exceed n.
  grid <- do.call(rbind, lapply(seq_len(max_n), function(k) {
    states <- expand.grid(n_pending = 0:k, n_tox = 0:k)
    states <- states[states$n_tox + states$n_pending <= k, , drop = FALSE]
    data.frame(n = rep(k, nrow(states)), n_tox = states$n_tox,
               n_pending = states$n_pending)
  }))

  cells <- tite_boin_bounds(grid$n, grid$n_tox, grid$n_pending, bound,
                            method = method,
                            max_pending_ratio = max_pending_ratio,
                            min_completed = as.integer(min_completed))
  rownames(cells) <- NULL

  structure(
    cells,
    class = c("tite_boin_decision_table", "data.frame"),
    method = method,
    statistic = if (method == "imputation") "STFT" else "ESS",
    target = target,
    p_saf = p_saf,
    p_tox = p_tox,
    lambda_e = bound$lambda_e,
    lambda_d = bound$lambda_d,
    cutoff_eli = cutoff_eli,
    max_pending_ratio = max_pending_ratio,
    min_completed = as.integer(min_completed)
  )
}
