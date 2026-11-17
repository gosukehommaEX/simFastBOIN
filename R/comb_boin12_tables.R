#' Admissibility Tables of Comb-BOIN12
#'
#' @param target,target_resp Numeric scalars. DLT and response limits.
#' @param cutoff_eli,cutoff_eli_resp Numeric scalars. Posterior probability
#'   cutoffs.
#' @param extrasafe Logical scalar. Use \code{cutoff_eli - offset} at the
#'   lowest combination.
#' @param offset Numeric scalar.
#' @param max_n Integer scalar. Largest number of patients at one
#'   combination.
#'
#' @return A list with the integer vectors \code{tox_ok}, \code{resp_ok}
#'   (as in \code{boin12_admissible_table()}) and \code{tox_ok_lowest}, the
#'   DLT table of the lowest combination.
#'
#' @noRd
comb_boin12_tables <- function(target, target_resp, cutoff_eli,
                               cutoff_eli_resp, extrasafe, offset, max_n) {
  adm <- boin12_admissible_table(target, target_resp, cutoff_eli,
                                 cutoff_eli_resp, max_n)
  adm$tox_ok_lowest <- if (extrasafe) {
    boin12_admissible_table(target, target_resp, cutoff_eli - offset,
                            cutoff_eli_resp, max_n)$tox_ok
  } else {
    adm$tox_ok
  }
  adm
}
