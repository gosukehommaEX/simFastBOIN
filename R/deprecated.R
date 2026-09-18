#' Deprecated Functions in simFastBOIN
#'
#' @description
#'   These functions were renamed in version 2.0.0. They still work but issue a
#'   deprecation warning and will be removed in a future release. Their return
#'   values are those of the replacement functions, which differ from the return
#'   values of version 1.3.2, so calls should be updated rather than left in place.
#'
#' @param ...
#'   Arguments passed to the replacement function. Note that the replacements use
#'   different argument names and a different argument order, so only calls that
#'   name every argument will carry over unchanged.
#'
#' @return
#'   The value of the corresponding replacement function.
#'
#' @section Replacements:
#'   \describe{
#'     \item{\code{get_boin_boundary()}}{use \code{\link{boin_lambda}} for the
#'       interval boundaries, or \code{\link{boin_boundary}} for the integer
#'       decision boundaries.}
#'     \item{\code{get_boin_decision()}}{use \code{\link{boin_decision_table}}.}
#'     \item{\code{get_boin_stopping_boundaries()}}{use \code{\link{boin_boundary}}
#'       with \code{extrasafe = TRUE} and read the \code{b_stop} component.}
#'     \item{\code{get_pts_and_tox()}}{use \code{\link{boin_simulate}}.}
#'     \item{\code{isotonic_regression()}}{use \code{\link{boin_isotonic}}.}
#'     \item{\code{select_mtd()}}{use \code{\link{boin_select_mtd}}, which now
#'       takes the patient and DLT counts instead of precomputed isotonic
#'       estimates.}
#'   }
#'
#' @name simFastBOIN-deprecated
#' @keywords internal
NULL

#' @rdname simFastBOIN-deprecated
#' @export
get_boin_boundary <- function(...) {
  .Deprecated("boin_lambda", package = "simFastBOIN")
  boin_lambda(...)
}

#' @rdname simFastBOIN-deprecated
#' @export
get_boin_decision <- function(...) {
  .Deprecated("boin_decision_table", package = "simFastBOIN")
  boin_decision_table(...)
}

#' @rdname simFastBOIN-deprecated
#' @export
get_boin_stopping_boundaries <- function(...) {
  .Deprecated("boin_boundary", package = "simFastBOIN")
  boin_boundary(...)
}

#' @rdname simFastBOIN-deprecated
#' @export
get_pts_and_tox <- function(...) {
  .Deprecated("boin_simulate", package = "simFastBOIN")
  boin_simulate(...)
}

#' @rdname simFastBOIN-deprecated
#' @export
isotonic_regression <- function(...) {
  .Deprecated("boin_isotonic", package = "simFastBOIN")
  boin_isotonic(...)
}

#' @rdname simFastBOIN-deprecated
#' @export
select_mtd <- function(...) {
  .Deprecated("boin_select_mtd", package = "simFastBOIN")
  boin_select_mtd(...)
}
