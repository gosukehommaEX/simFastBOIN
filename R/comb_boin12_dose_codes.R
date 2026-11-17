#' Standardized Dose Levels of the Logistic Model of Comb-BOIN12
#'
#' @param dose Numeric vector of dose values or \code{NULL} for the level
#'   numbers.
#' @param n_levels Integer scalar. Number of levels of the drug.
#' @param name Character scalar. Argument name for messages.
#'
#' @return Numeric vector of the doses shifted to mean 0 and scaled to
#'   standard deviation 0.5 (\code{sd()}, denominator \code{n_levels - 1});
#'   zero for a single level.
#'
#' @noRd
comb_boin12_dose_codes <- function(dose, n_levels, name) {
  if (is.null(dose)) dose <- seq_len(n_levels)
  if (!is.numeric(dose) || length(dose) != n_levels ||
      any(!is.finite(dose)) || any(diff(dose) <= 0)) {
    stop("'", name, "' must hold one increasing dose value per level of ",
         "the drug", call. = FALSE)
  }
  if (n_levels == 1L) return(0)
  (dose - mean(dose)) / stats::sd(dose) * 0.5
}
