#' Boundaries of the BOIN-ET Design
#'
#' @description
#'   Find the lower and upper DLT boundaries and the lower response boundary
#'   of the BOIN-ET design of Takeda, Taguri and Morita (2018). The three
#'   boundaries minimize the probability of an incorrect dose assignment under
#'   six hypotheses about the DLT and response probabilities at the current
#'   dose, and are searched on a grid of step 0.01, as in the 'boinet' package.
#'
#' @param target
#'   Numeric scalar. Target DLT probability, for example 0.30.
#'
#' @param target_resp
#'   Numeric scalar. Target response probability, for example 0.60.
#'
#' @param p_saf
#'   Numeric scalar. Highest DLT probability deemed subtherapeutic. Defaults
#'   to \code{0.1 * target}, the value recommended by Takeda et al. (2018),
#'   which is far below the \code{0.6 * target} of the BOIN design.
#'
#' @param p_tox
#'   Numeric scalar. Lowest DLT probability deemed overly toxic. Defaults to
#'   \code{1.4 * target}.
#'
#' @param p_saf_resp
#'   Numeric scalar. Highest response probability deemed subtherapeutic.
#'   Defaults to \code{0.6 * target_resp}.
#'
#' @param prior
#'   Numeric vector of length 6. Prior probabilities of the six hypotheses.
#'   Defaults to \code{rep(1 / 6, 6)}, the noninformative choice of the
#'   article. The third hypothesis, under which every decision is correct,
#'   does not enter the criterion.
#'
#' @return
#'   An object of class \code{boin_et_boundary}, which is a list with
#'   components
#'   \item{lambda1}{Lower DLT boundary. The dose is escalated when the observed
#'     DLT rate is at most \code{lambda1} and the observed response rate is at
#'     most \code{eta1}.}
#'   \item{lambda2}{Upper DLT boundary. The dose is de-escalated when the
#'     observed DLT rate is at least \code{lambda2}.}
#'   \item{eta1}{Lower response boundary.}
#'   \item{prob_incorrect}{The minimized probability of an incorrect
#'     decision.}
#'   \item{target, target_resp, p_saf, p_tox, p_saf_resp, prior}{The design
#'     parameters.}
#'
#' @details
#'   The grids run from \code{p_saf} to \code{target} for \code{lambda1}, from
#'   \code{target} to \code{p_tox} for \code{lambda2} and from
#'   \code{p_saf_resp} to \code{target_resp} for \code{eta1}, all by 0.01. At
#'   each grid point the probability of an incorrect decision is evaluated
#'   with binomial probabilities for 100 patients, and the first minimum in
#'   the order of the grids is returned. This is the search of
#'   \code{boinet::gridoptim()}, including its expression for the fourth
#'   hypothesis, and with the default design parameters and a target of 0.30
#'   for the DLT and 0.60 for the response it gives the boundaries 0.14, 0.35
#'   and 0.48 reported by Takeda et al. (2018).
#'
#' @references
#'   Takeda, K., Taguri, M. and Morita, S. (2018). BOIN-ET: Bayesian optimal
#'   interval design for dose finding based on both efficacy and toxicity
#'   outcomes. Pharmaceutical Statistics, 17(4), 383-395.
#'   \doi{10.1002/pst.1864}
#'
#'   Yamaguchi, Y., Takeda, K. and Maruo, K. (2025). boinet: An R Package for
#'   Bayesian Optimal Interval Design for Dose Finding Based on Both Efficacy
#'   and Toxicity Outcomes. The R Journal, 17(2), 15-34.
#'   \doi{10.32614/RJ-2025-012}
#'
#' @examples
#' bd <- boin_et_boundary(target = 0.30, target_resp = 0.60)
#' bd
#'
#' c(bd$lambda1, bd$lambda2, bd$eta1)
#'
#' @seealso \code{\link{boin_et_decision_table}}, \code{\link{sim_boin_et}}
#'
#' @export
boin_et_boundary <- function(target, target_resp, p_saf = NULL, p_tox = NULL,
                             p_saf_resp = NULL, prior = rep(1 / 6, 6)) {

  par <- boin_et_check_design(target, target_resp, p_saf, p_tox, p_saf_resp)
  if (!is.numeric(prior) || length(prior) != 6L || any(!is.finite(prior)) ||
      any(prior < 0)) {
    stop("'prior' must be six non-negative finite numbers", call. = FALSE)
  }

  phi <- par$target
  phi1 <- par$p_saf
  phi2 <- par$p_tox
  delta <- par$target_resp
  delta1 <- par$p_saf_resp
  n <- 100

  l1 <- seq(phi1, phi, by = 0.01)
  l2 <- seq(phi, phi2, by = 0.01)
  e1 <- seq(delta1, delta, by = 0.01)

  # Binomial probabilities at every grid value, as in the reference.
  l1p <- stats::pbinom(n * l1, n, phi)
  l1p1 <- stats::pbinom(n * l1, n, phi1)
  l1p2 <- stats::pbinom(n * l1, n, phi2)
  l2p <- stats::pbinom(n * l2 - 1, n, phi)
  l2p1 <- stats::pbinom(n * l2 - 1, n, phi1)
  l2p2 <- stats::pbinom(n * l2 - 1, n, phi2)
  e1d <- stats::pbinom(n * e1, n, delta)
  e1d1 <- stats::pbinom(n * e1, n, delta1)

  # Every combination of the three grids, the first grid varying fastest, so
  # that the first minimum is the one found by the reference.
  s1 <- rep(seq_along(l1), times = length(l2) * length(e1))
  s2 <- rep(rep(seq_along(l2), each = length(l1)), times = length(e1))
  s3 <- rep(seq_along(e1), each = length(l1) * length(l2))

  a <- l1p[s1]
  a1 <- l1p1[s1]
  a2 <- l1p2[s1]
  b <- l2p[s2]
  b1 <- l2p1[s2]
  b2 <- l2p2[s2]
  e <- e1d[s3]
  f <- e1d1[s3]

  value <- (prior[1] * (a1 * (1 - f) + 2 / 3 * (b1 - a1) * f + (b1 - a1) * (1 - f) + (1 - b1))
            + prior[2] * (a1 * e + 2 / 3 * (b1 - a1) * e + (1 - b1))
            + prior[4] * (a * e + 2 / 3 * (b - a1) * e + (1 - b))
            + prior[5] * (a2 * f + a2 * (1 - f) + 2 / 3 * (b2 - a2) * f + (b2 - a2) * (1 - f))
            + prior[6] * (a2 * e + a2 * (1 - e) + 2 / 3 * (b2 - a2) * e + (b2 - a2) * (1 - e)))

  best <- which(value == min(value))[1L]

  structure(
    list(
      lambda1 = l1[s1[best]],
      lambda2 = l2[s2[best]],
      eta1 = e1[s3[best]],
      prob_incorrect = value[best],
      target = phi,
      target_resp = delta,
      p_saf = phi1,
      p_tox = phi2,
      p_saf_resp = delta1,
      prior = prior
    ),
    class = "boin_et_boundary"
  )
}
