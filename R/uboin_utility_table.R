#' Utility Table of the U-BOIN Design
#'
#' @description
#'   Tabulate, for a dose with a given number of patients, the posterior mean
#'   utility and the admissibility of the dose for every possible number of
#'   responses, of DLTs and of responses without a DLT. These are the decision
#'   tables of the U-BOIN design of Zhou, Lee and Yuan (2019), which can be
#'   written into the protocol before the trial starts.
#'
#' @param n
#'   Integer scalar. Number of patients treated at the dose.
#'
#' @param tox_upper
#'   Numeric scalar. Highest acceptable DLT probability.
#'
#' @param resp_lower
#'   Numeric scalar. Lowest acceptable response probability.
#'
#' @param utility
#'   Numeric vector of length four. Utility of the outcomes (no response, DLT),
#'   (no response, no DLT), (response, DLT) and (response, no DLT). Defaults to
#'   \code{c(0, 30, 50, 100)}.
#'
#' @param cutoff_eli
#'   Numeric scalar. Cutoff of the toxicity rule. Defaults to 0.95.
#'
#' @param cutoff_eli_resp
#'   Numeric scalar. Cutoff of the futility rule. Defaults to 0.90.
#'
#' @param prior_tox
#'   Numeric vector of length two. Beta prior of the toxicity rule.
#'   Defaults to \code{c(1, 1)}.
#'
#' @param prior_resp
#'   Numeric vector of length two. Beta prior of the futility rule.
#'   Defaults to \code{c(1, 1)}.
#'
#' @return
#'   A data frame with one row per possible outcome of the \code{n} patients
#'   and columns
#'   \item{n_resp}{Number of responses.}
#'   \item{n_tox}{Number of DLTs.}
#'   \item{n_resp_no_tox}{Number of responses without a DLT.}
#'   \item{utility}{Posterior mean utility of the dose.}
#'   \item{toxic}{Whether the dose is toxic.}
#'   \item{futile}{Whether the dose is futile.}
#'   \item{admissible}{Whether the dose is neither toxic nor futile.}
#'
#' @details
#'   The posterior mean utility is \code{sum(utility * (1/4 + n_k)) / (1 + n)},
#'   with \code{n_k} the patients with each outcome, as in
#'   \code{\link{uboin_simulate}}. A dose is toxic when the posterior
#'   probability that its DLT probability exceeds \code{tox_upper} is above
#'   \code{cutoff_eli}, and futile when the posterior probability that its
#'   response probability is below \code{resp_lower} is above
#'   \code{cutoff_eli_resp}.
#'
#'   Tables S4 to S6 of the supporting information of the article, for 3, 6 and
#'   9 patients and the utilities \code{c(0, 15, 25, 100)}, are reproduced
#'   with \code{prior_resp = c(0.05, 0.05)}; with the default Beta(1, 1) prior a
#'   dose without any response among 3, 6 or 9 patients is admissible, whereas
#'   those tables mark it inadmissible.
#'
#' @references
#'   Zhou, Y., Lee, J. J. and Yuan, Y. (2019). A utility-based Bayesian optimal
#'   interval (U-BOIN) phase I/II design to identify the optimal biological
#'   dose for targeted and immune therapies. Statistics in Medicine, 38(28),
#'   S5299-S5316. \doi{10.1002/sim.8361}
#'
#' @examples
#' # Table S4 of the supporting information of the article
#' tab <- uboin_utility_table(3, tox_upper = 0.30, resp_lower = 0.20,
#'                            utility = c(0, 15, 25, 100),
#'                            prior_resp = c(0.05, 0.05))
#' tab[tab$admissible, ]
#'
#' @seealso \code{\link{uboin_simulate}}, \code{\link{uboin_select_obd}}
#'
#' @export
uboin_utility_table <- function(n, tox_upper, resp_lower,
                                utility = c(0, 30, 50, 100),
                                cutoff_eli = 0.95, cutoff_eli_resp = 0.90,
                                prior_tox = c(1, 1), prior_resp = c(1, 1)) {

  check_count(n, "n", 1L)
  des <- uboin_check_design(tox_upper, resp_lower, NULL, NULL, NULL, utility,
                            cutoff_eli, cutoff_eli_resp, prior_tox,
                            prior_resp, stage1 = FALSE)
  n <- as.integer(n)

  grid <- expand.grid(n_resp_no_tox = 0:n, n_tox = 0:n, n_resp = 0:n)
  grid <- grid[, c("n_resp", "n_tox", "n_resp_no_tox")]
  resp_tox <- grid$n_resp - grid$n_resp_no_tox
  no_resp_tox <- grid$n_tox - resp_tox
  no_resp_no_tox <- n - grid$n_resp - no_resp_tox
  keep <- grid$n_resp_no_tox <= grid$n_resp & resp_tox <= grid$n_tox &
    no_resp_no_tox >= 0
  grid <- grid[keep, ]
  counts <- cbind(no_resp_tox[keep], no_resp_no_tox[keep], resp_tox[keep],
                  grid$n_resp_no_tox)

  utility_est <- as.numeric((counts + 0.25) %*% des$utility) / (1 + n)
  above <- 1 - stats::pbeta(des$tox_upper, des$prior_tox[1L] + grid$n_tox,
                            des$prior_tox[2L] + n - grid$n_tox)
  below <- stats::pbeta(des$resp_lower, des$prior_resp[1L] + grid$n_resp,
                        des$prior_resp[2L] + n - grid$n_resp)
  toxic <- above > des$cutoff_eli
  futile <- below > des$cutoff_eli_resp

  out <- data.frame(
    n_resp = grid$n_resp,
    n_tox = grid$n_tox,
    n_resp_no_tox = grid$n_resp_no_tox,
    utility = utility_est,
    toxic = toxic,
    futile = futile,
    admissible = !toxic & !futile
  )
  rownames(out) <- NULL
  out
}
