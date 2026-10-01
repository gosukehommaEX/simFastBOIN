#' Print a 3+3 Decision Table
#'
#' @description
#'   Display the decision table produced by \code{\link{decision_table_3p3}}:
#'   one row per number of patients and range of DLTs sharing a decision, for
#'   dose escalation and, with \code{mtd_rule = "expand"}, for the search for
#'   the MTD, followed by what happens once escalation stops.
#'
#' @param x
#'   An object of class \code{decision_table_3p3}.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   The object \code{x}, invisibly.
#'
#' @details
#'   A data frame that no longer has the columns or attributes of a decision
#'   table, for example after selecting columns, is printed as a plain data
#'   frame.
#'
#' @examples
#' print(decision_table_3p3())
#'
#' print(decision_table_3p3(mtd_rule = "expand"))
#'
#' @seealso \code{\link{decision_table_3p3}},
#'   \code{\link{plot.decision_table_3p3}}
#'
#' @export
print.decision_table_3p3 <- function(x, ...) {

  required <- c("stage", "n", "n_tox", "decision")
  if (is.null(attr(x, "mtd_rule")) || !all(required %in% names(x))) {
    return(NextMethod())
  }

  rule <- attr(x, "mtd_rule")
  rule_label <- if (rule == "previous") {
    "previous (the dose below the one declared too toxic)"
  } else {
    "expand (the highest dose with at most 1 DLT in 6)"
  }
  cat("3+3 decision table, MTD rule: ", rule_label, "\n", sep = "")

  rows <- decision_table_3p3_rows(x)
  for (stage in unique(rows$Stage)) {
    cat("\n", stage, "\n", sep = "")
    print(rows[rows$Stage == stage, c("Patients", "DLTs", "Decision")],
          row.names = FALSE, right = FALSE)
  }

  cat("\n")
  cat(strwrap(paste(decision_table_3p3_notes(rule), collapse = " ")), sep = "\n")

  invisible(x)
}
