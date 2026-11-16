#' Desirability Score Table of the BOIN12 Design
#'
#' @description
#'   Tabulate the desirability of a dose for every possible outcome of the
#'   patients treated there, as in Table 3 of Lin et al. (2020): the
#'   posterior probability that the utility of the dose exceeds the
#'   benchmark, its rank-based desirability score (RDS) and whether the dose
#'   is admissible.
#'
#' @inheritParams boin12_simulate
#'
#' @param n
#'   Integer vector. Numbers of patients at a dose to tabulate. Defaults to
#'   \code{c(0, 3, 6, 9)}.
#'
#' @return
#'   A data frame of class \code{boin12_decision_table} with one row per
#'   outcome and columns \code{n}, \code{n_tox}, \code{n_resp},
#'   \code{n_tox_resp} (patients with both a DLT and a response; only when
#'   \code{psi00 + psi11} is not 100, since the utility otherwise depends on
#'   the margins alone), \code{admissible}, \code{desirability} and
#'   \code{score}. The design parameters are kept in the attribute
#'   \code{"design"}, a list with \code{lambda_e}, \code{lambda_d},
#'   \code{u_low} and \code{u_benchmark} among others.
#'
#' @details
#'   The utility of a patient is 100 for a response without DLT,
#'   \code{psi00} for neither, \code{psi11} for both and 0 for a DLT without
#'   response. The utility of a dose, divided by 100, has a beta posterior
#'   under the quasi-binomial model with the prior \code{utility_prior}, and
#'   the desirability is its posterior probability of exceeding
#'   \code{u_benchmark / 100}. The score ranks the desirability among the
#'   admissible outcomes of the table, from 1 for the lowest, equal
#'   desirabilities sharing the lowest rank. A dose is inadmissible
#'   (\code{"E"} in the printed table) when the posterior probability, under
#'   a uniform prior, that its DLT rate exceeds \code{target} is at least
#'   \code{cutoff_eli}, or that its response rate is below
#'   \code{target_resp} is at least \code{cutoff_eli_resp}. During a trial
#'   only the order of the desirabilities of the doses compared matters, so
#'   the score gives the same decisions.
#'
#' @references
#'   Lin, R., Zhou, Y., Yan, F., Li, D. and Yuan, Y. (2020). BOIN12: Bayesian
#'   optimal interval phase I/II trial design for utility-based dose finding
#'   in immunotherapy and targeted therapies. JCO Precision Oncology, 4,
#'   1393-1402. \doi{10.1200/PO.20.00257}
#'
#' @examples
#' # Table 3 of Lin et al. (2020)
#' tab <- boin12_decision_table(target = 0.35, target_resp = 0.25)
#' tab
#'
#' # A utility whose scores psi00 and psi11 do not add up to 100
#' boin12_decision_table(target = 0.35, target_resp = 0.25, n = c(0, 3),
#'                       psi00 = 30, psi11 = 50)
#'
#' @seealso \code{\link{boin12_simulate}}, \code{\link{sim_boin12}}
#'
#' @export
boin12_decision_table <- function(target, target_resp, n = c(0, 3, 6, 9),
                                  p_saf = NULL, p_tox = NULL,
                                  cutoff_eli = 0.95, cutoff_eli_resp = 0.90,
                                  psi00 = 40, psi11 = 60,
                                  benchmark_weight = 0.5,
                                  utility_prior = c(1, 1)) {

  des <- boin12_check_design(target, target_resp, p_saf = p_saf,
                             p_tox = p_tox, psi00 = psi00, psi11 = psi11,
                             benchmark_weight = benchmark_weight,
                             utility_prior = utility_prior)
  check_scalar_prob(cutoff_eli, "cutoff_eli")
  check_scalar_prob(cutoff_eli_resp, "cutoff_eli_resp")
  if (!is.numeric(n) || length(n) < 1L || any(!is.finite(n)) ||
      any(n < 0) || any(n != round(n))) {
    stop("'n' must hold whole numbers of patients, zero or more",
         call. = FALSE)
  }
  n <- sort(unique(as.integer(n)))
  joint <- abs(psi00 + psi11 - 100) > 1e-9

  rows <- lapply(n, function(m) {
    grid <- expand.grid(n_tox = 0:m, n_resp = 0:m,
                        n_tox_resp = if (joint) 0:m else 0L)
    if (joint) {
      keep <- grid$n_tox_resp <= pmin(grid$n_tox, grid$n_resp) &
        grid$n_resp - grid$n_tox_resp <= m - grid$n_tox
      grid <- grid[keep, , drop = FALSE]
    }
    grid <- grid[order(grid$n_tox, grid$n_resp, grid$n_tox_resp), ,
                 drop = FALSE]
    data.frame(n = rep(m, nrow(grid)), grid, row.names = NULL)
  })
  tab <- do.call(rbind, rows)

  y <- tab$n_tox
  x <- tab$n_resp
  z <- if (joint) tab$n_tox_resp else
    pmax(0, x + y - tab$n)
  q <- (100 * (x - z) + psi00 * (tab$n - y - x + z) + psi11 * z) / 100
  a <- des$utility_prior[1L]
  b <- des$utility_prior[2L]
  tab$admissible <-
    stats::pbeta(target, 1 + y, 1 + tab$n - y, lower.tail = FALSE) <
      cutoff_eli &
    stats::pbeta(target_resp, 1 + x, 1 + tab$n - x) < cutoff_eli_resp
  tab$desirability <- stats::pbeta(des$u_benchmark / 100, a + q,
                                   b + tab$n - q, lower.tail = FALSE)
  tab$score <- NA_integer_
  if (any(tab$admissible)) {
    tab$score[tab$admissible] <- as.integer(
      rank(round(tab$desirability[tab$admissible], 12), ties.method = "min")
    )
  }
  if (!joint) tab$n_tox_resp <- NULL
  rownames(tab) <- NULL

  des$cutoff_eli <- cutoff_eli
  des$cutoff_eli_resp <- cutoff_eli_resp
  attr(tab, "design") <- des
  class(tab) <- c("boin12_decision_table", "data.frame")
  tab
}
