#' DLT Exclusion Table of BOIN-ETC
#'
#' @description
#'   Whether the DLT data of a combination leave it in the trial, for every
#'   count of DLTs among 0 to \code{max_n} patients: a combination is
#'   excluded when Pr(DLT rate > \code{p_eli}) under the Beta(1, 1) prior is
#'   above \code{cutoff_eli}, as in \code{boin_et_admissible_table()}.
#'   Internal.
#'
#' @param p_eli,cutoff_eli
#'   Numeric scalars.
#'
#' @param max_n
#'   Integer scalar. Largest number of patients at one combination.
#'
#' @return An integer vector in the layout read by the C++ engine
#'   (\code{n * (max_n + 1) + y}): 1 when the combination stays.
#'
#' @keywords internal
#' @noRd
boin_etc_tox_table <- function(p_eli, cutoff_eli, max_n) {
  as.integer(boin_et_admissible_table(p_eli, 0.5, cutoff_eli, 0.5,
                                      max_n)$tox_ok)
}
