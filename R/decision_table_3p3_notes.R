#' What Follows the Rows of a 3+3 Decision Table
#'
#' @description
#'   Internal helper of \code{\link{print.decision_table_3p3}}. Sentences that
#'   state what happens once escalation stops, which depends on the definition
#'   of the MTD.
#'
#' @param mtd_rule
#'   Character scalar, \code{"previous"} or \code{"expand"}.
#'
#' @return
#'   A character vector of sentences.
#'
#' @noRd
decision_table_3p3_notes <- function(mtd_rule) {

  stop_note <- "Stop escalation means that this dose and all higher doses are too toxic."
  if (mtd_rule == "previous") {
    c(stop_note,
      paste("The MTD is the dose below the one at which escalation stopped, and",
            "no MTD is selected when escalation stops at the starting dose."),
      "When the highest dose is escalated from, it is the MTD.")
  } else {
    c(stop_note,
      paste("The search for the MTD starts at the dose below the one at which",
            "escalation stopped, or at the highest dose when it is escalated from."),
      "No MTD is selected when the search moves below the starting dose.")
  }
}
