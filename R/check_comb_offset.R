#' Check the Offset of the Extra Safety Cutoff
#'
#' @param offset Numeric scalar, the amount by which \code{cutoff_eli} is
#'   lowered at the lowest combination.
#' @param cutoff_eli Numeric scalar, already checked.
#'
#' @return \code{TRUE}, invisibly, or an error.
#'
#' @noRd
check_comb_offset <- function(offset, cutoff_eli) {
  if (!is.numeric(offset) || length(offset) != 1L || !is.finite(offset) ||
      offset < 0 || offset >= cutoff_eli) {
    stop("'offset' must be a single number from 0 to below 'cutoff_eli'",
         call. = FALSE)
  }
  invisible(TRUE)
}
