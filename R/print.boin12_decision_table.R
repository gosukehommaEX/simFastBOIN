#' Print the Desirability Score Table of the BOIN12 Design
#'
#' @description
#'   Display the table of \code{\link{boin12_decision_table}} with the score
#'   of every outcome, \code{"E"} marking the outcomes that make the dose
#'   inadmissible.
#'
#' @param x
#'   An object of class \code{boin12_decision_table}.
#'
#' @param desirability
#'   Logical scalar. Show the posterior probabilities as well. Defaults to
#'   \code{FALSE}.
#'
#' @param ...
#'   Not used.
#'
#' @return
#'   \code{x}, invisibly.
#'
#' @examples
#' tab <- boin12_decision_table(target = 0.35, target_resp = 0.25, n = c(0, 3))
#' print(tab, desirability = TRUE)
#'
#' @export
print.boin12_decision_table <- function(x, desirability = FALSE, ...) {

  check_flag(desirability, "desirability")
  des <- attr(x, "design")

  out <- data.frame(n = x$n, n_tox = x$n_tox, n_resp = x$n_resp)
  if (!is.null(x$n_tox_resp)) out$n_tox_resp <- x$n_tox_resp
  out$score <- ifelse(x$admissible, format(x$score), "E")
  if (desirability) out$desirability <- round(x$desirability, 4)

  cat("BOIN12 desirability score table\n")
  cat("  DLT limit / response limit : ", format(des$target), " / ",
      format(des$target_resp), "\n", sep = "")
  cat("  lambda_e / lambda_d        : ", format(round(des$lambda_e, 3)), " / ",
      format(round(des$lambda_d, 3)), "\n", sep = "")
  cat("  utility scores             : 100 (response only), ",
      format(des$psi00), " (neither), ", format(des$psi11), " (both), ",
      "0 (DLT only)\n", sep = "")
  cat("  utility benchmark          : ", format(round(des$u_benchmark, 2)),
      "\n", sep = "")
  cat("  E: inadmissible (DLT cutoff ", format(des$cutoff_eli),
      ", response cutoff ", format(des$cutoff_eli_resp), ")\n\n", sep = "")
  print(out, row.names = FALSE, right = TRUE)

  invisible(x)
}
