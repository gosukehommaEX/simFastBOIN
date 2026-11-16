#' Admissibility Tables of BOIN12
#'
#' @param target,target_resp Numeric scalars. DLT and response limits.
#' @param cutoff_eli,cutoff_eli_resp Numeric scalars. Posterior probability
#'   cutoffs.
#' @param max_n Integer scalar. Largest number of patients at one dose.
#'
#' @return A list with the integer vectors \code{tox_ok} and \code{resp_ok}:
#'   element \code{n * (max_n + 1) + k + 1} is 1 when \code{k} DLTs (or
#'   responses) among \code{n} patients do not make the dose inadmissible,
#'   that is when the posterior probability, under a uniform prior, that the
#'   DLT rate exceeds \code{target} is below \code{cutoff_eli} (and that the
#'   response rate is below \code{target_resp} is below
#'   \code{cutoff_eli_resp}).
#'
#' @noRd
boin12_admissible_table <- function(target, target_resp, cutoff_eli,
                                    cutoff_eli_resp, max_n) {
  w <- max_n + 1L
  n <- rep(0:max_n, each = w)
  k <- rep(0:max_n, times = w)
  ok <- k <= n
  tox_ok <- integer(w * w)
  resp_ok <- integer(w * w)
  tox_ok[ok] <- as.integer(
    stats::pbeta(target, 1 + k[ok], 1 + n[ok] - k[ok], lower.tail = FALSE) <
      cutoff_eli
  )
  resp_ok[ok] <- as.integer(
    stats::pbeta(target_resp, 1 + k[ok], 1 + n[ok] - k[ok]) < cutoff_eli_resp
  )
  list(tox_ok = tox_ok, resp_ok = resp_ok)
}
