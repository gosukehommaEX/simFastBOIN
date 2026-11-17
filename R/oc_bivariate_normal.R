#' Joint Probability of a DLT and a Response under the Gaussian Copula
#'
#' @description
#'   Internal helper of \code{\link{oc_scenarios}}. The probability that a
#'   patient has both a DLT and a response when the two outcomes are
#'   dichotomized latent standard normal variables with correlation
#'   \code{rho}, as in the simulation functions with
#'   \code{resp_cor_type = "copula"}: Pr(Z1 < qnorm(p1), Z2 < qnorm(p2)).
#'
#' @param p1 Numeric scalar, the DLT probability.
#' @param p2 Numeric scalar, the response probability.
#' @param rho Numeric scalar strictly between -1 and 1, the correlation.
#'
#' @return Numeric scalar.
#'
#' @noRd
oc_bivariate_normal <- function(p1, p2, rho) {
  if (p1 <= 0 || p2 <= 0) return(0)
  if (p1 >= 1) return(p2)
  if (p2 >= 1) return(p1)
  if (abs(rho) < .Machine$double.eps) return(p1 * p2)
  a <- stats::qnorm(p1)
  b <- stats::qnorm(p2)
  s <- sqrt(1 - rho^2)
  # Pr(Z1 < a, Z2 < b) = integral over z < a of phi(z) Phi((b - rho z) / s)
  stats::integrate(
    function(z) stats::dnorm(z) * stats::pnorm((b - rho * z) / s),
    lower = -Inf, upper = a, rel.tol = 1e-10, abs.tol = 0
  )$value
}
