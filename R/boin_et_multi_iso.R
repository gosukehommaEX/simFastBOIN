#' Model-Averaged Unimodal Isotonic Estimate of the Response Rates
#'
#' @description
#'   Estimate the response rates of the doses tried in a BOIN-ET trial by
#'   averaging the unimodal isotonic regressions with every possible mode,
#'   weighted by their AIC, after Lin and Yin (2017). This follows
#'   \code{boinet::multi.iso()} (boinet 1.6.0, MIT license, by Yusuke Yamaguchi
#'   and Kentaro Takeda), including its unweighted fit to the observed rates.
#'   Internal.
#'
#' @param obs
#'   Integer vector. Responses at each dose tried, in dose order.
#'
#' @param n
#'   Integer vector. Patients at each dose tried.
#'
#' @return
#'   Numeric vector of the estimated response rates.
#'
#' @keywords internal
#' @noRd
boin_et_multi_iso <- function(obs, n) {
  prob <- obs / n
  ld <- length(prob)
  unimod <- array(0, dim = c(ld, ld))
  aic <- numeric(ld)
  for (k in seq_len(ld)) {
    unimod[k, ] <- (Iso::ufit(y = prob, x = seq_len(ld), lmode = k))$y
    lk <- prod(stats::dbinom(obs, n, unimod[k, ]))
    aic[k] <- (-2 * log(lk) + 2 * ld)
  }
  pik <- exp(-aic / 2) / sum(exp(-aic / 2))
  apply(as.matrix(seq_len(ld)), 1, function(x) sum(pik * unimod[, x]))
}
