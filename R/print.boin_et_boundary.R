#' Print the Boundaries of the BOIN-ET Design
#'
#' @description
#'   Display the boundaries found by \code{\link{boin_et_boundary}} together
#'   with the design parameters they were found from.
#'
#' @param x
#'   An object of class \code{boin_et_boundary}.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   The object \code{x}, invisibly.
#'
#' @examples
#' print(boin_et_boundary(target = 0.30, target_resp = 0.60))
#'
#' @seealso \code{\link{boin_et_boundary}}
#'
#' @export
print.boin_et_boundary <- function(x, ...) {
  cat("BOIN-ET boundaries\n")
  cat("  target DLT rate      : ", format(x$target), "\n", sep = "")
  cat("  target response rate : ", format(x$target_resp), "\n", sep = "")
  cat("  p_saf / p_tox        : ", format(x$p_saf), " / ", format(x$p_tox),
      "\n", sep = "")
  cat("  p_saf_resp           : ", format(x$p_saf_resp), "\n", sep = "")
  cat("\n")
  cat("  lambda1 (escalate)   : ", format(round(x$lambda1, 4)), "\n", sep = "")
  cat("  lambda2 (de-escalate): ", format(round(x$lambda2, 4)), "\n", sep = "")
  cat("  eta1 (response)      : ", format(round(x$eta1, 4)), "\n", sep = "")
  cat("  probability of an incorrect decision: ",
      format(signif(x$prob_incorrect, 4)), "\n", sep = "")
  invisible(x)
}
