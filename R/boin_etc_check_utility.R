#' Check the Utility Settings of BOIN-ETC
#'
#' @description
#'   Check \code{psi00}, \code{psi11} and \code{utility_prior} of the final
#'   selection of BOIN-ETC. Internal.
#'
#' @param psi00,psi11
#'   Numeric scalars from 0 to 100.
#'
#' @param utility_prior
#'   Numeric vector of two positive numbers.
#'
#' @return \code{TRUE}, invisibly.
#'
#' @keywords internal
#' @noRd
boin_etc_check_utility <- function(psi00, psi11, utility_prior) {
  for (nm in c("psi00", "psi11")) {
    value <- get(nm)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0 || value > 100) {
      stop("'", nm, "' must be a single number from 0 to 100", call. = FALSE)
    }
  }
  if (!is.numeric(utility_prior) || length(utility_prior) != 2L ||
      any(!is.finite(utility_prior)) || any(utility_prior <= 0)) {
    stop("'utility_prior' must hold two positive numbers", call. = FALSE)
  }
  invisible(TRUE)
}
