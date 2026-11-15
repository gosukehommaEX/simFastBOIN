#' Print a BOIN-ET Decision Table
#'
#' @description
#'   Display the decision table produced by \code{\link{boin_et_decision_table}}
#'   as one matrix per number of patients, with the DLTs in the rows and the
#'   responses in the columns.
#'
#' @param x
#'   An object of class \code{boin_et_decision_table}.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   The object \code{x}, invisibly.
#'
#' @examples
#' print(boin_et_decision_table(target = 0.30, target_resp = 0.60, n = 3))
#'
#' @seealso \code{\link{boin_et_decision_table}}
#'
#' @export
print.boin_et_decision_table <- function(x, ...) {
  df <- x
  class(df) <- "data.frame"
  bd <- attr(x, "boundary")

  cat("BOIN-ET decision table\n")
  cat("  lambda1 = ", format(round(bd$lambda1, 4)),
      ", lambda2 = ", format(round(bd$lambda2, 4)),
      ", eta1 = ", format(round(bd$eta1, 4)), "\n", sep = "")

  for (k in unique(df$n)) {
    part <- df[df$n == k, , drop = FALSE]
    entry <- ifelse(part$tox_excluded | part$resp_excluded,
                    paste0(part$decision, "*"), part$decision)
    out <- matrix("", nrow = k + 1L, ncol = k + 1L,
                  dimnames = list(DLTs = 0:k, Responses = 0:k))
    out[cbind(part$n_tox + 1L, part$n_resp + 1L)] <- entry
    cat("\n", k, " patients\n", sep = "")
    print(out, quote = FALSE)
  }

  cat("\nE = escalate, S = stay, D = de-escalate, E/S/D = choose among the",
      "current and neighboring doses;\n* = the dose is no longer admissible\n")
  invisible(x)
}
