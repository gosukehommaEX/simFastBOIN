#' Points of the Importance Sampling of Comb-BOIN12
#'
#' @param n_draws Integer scalar. Number of points.
#' @param df Numeric scalar. Degrees of freedom of the t proposal.
#'
#' @return A list with \code{z}, an \code{n_draws} by 5 matrix of standard
#'   normal quantiles, and \code{w}, a vector of chi-square quantiles with
#'   \code{df} degrees of freedom, of the Halton points 1 to \code{n_draws}
#'   in the bases 2, 3, 5, 7, 11 and 13. No random number is used.
#'
#' @noRd
comb_boin12_points <- function(n_draws, df = 4) {
  bases <- c(2, 3, 5, 7, 11, 13)
  u <- vapply(bases, function(b) {
    i <- seq_len(n_draws)
    out <- numeric(n_draws)
    f <- 1
    while (any(i > 0)) {
      f <- f / b
      out <- out + f * (i %% b)
      i <- i %/% b
    }
    out
  }, numeric(n_draws))
  u <- matrix(u, nrow = n_draws)
  list(z = stats::qnorm(u[, 1:5, drop = FALSE]),
       w = stats::qchisq(u[, 6L], df))
}
