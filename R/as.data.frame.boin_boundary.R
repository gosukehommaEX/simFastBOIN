#' Coerce BOIN Decision Boundaries to a Data Frame
#'
#' @description
#'   Return the integer decision boundaries produced by \code{\link{boin_boundary}}
#'   as one row per sample size.
#'
#' @param x
#'   An object of class \code{boin_boundary}.
#'
#' @param row.names
#'   Passed on for compatibility with the generic, normally unused.
#'
#' @param optional
#'   Passed on for compatibility with the generic, normally unused.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   A data frame with columns \code{n_pts}, \code{escalate_if_dlt_leq},
#'   \code{deescalate_if_dlt_geq} and \code{eliminate_if_dlt_geq}, plus
#'   \code{stop_if_dlt_geq} when the boundaries were built with
#'   \code{extrasafe = TRUE}.
#'
#' @examples
#' bd <- boin_boundary(target = 0.30, max_n = 18)
#' as.data.frame(bd)[seq(3, 18, by = 3), ]
#'
#' @export
as.data.frame.boin_boundary <- function(x, row.names = NULL, optional = FALSE, ...) {
  out <- data.frame(
    n_pts = x$n,
    escalate_if_dlt_leq = x$b_esc,
    deescalate_if_dlt_geq = x$b_deesc,
    eliminate_if_dlt_geq = x$b_elim,
    stringsAsFactors = FALSE
  )
  if (x$extrasafe) out$stop_if_dlt_geq <- x$b_stop
  out
}
