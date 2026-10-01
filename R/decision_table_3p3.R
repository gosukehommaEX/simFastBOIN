#' Decision Table for the 3+3 Design
#'
#' @description
#'   Tabulate the dose assignment rule of the traditional 3+3 design used by
#'   \code{\link{oc_3p3}} and \code{\link{sim_3p3}}: for every number of DLTs
#'   among the three or six patients treated at the current dose, the decision
#'   taken during dose escalation and, for \code{mtd_rule = "expand"}, during the
#'   search for the MTD that follows it.
#'
#' @param mtd_rule
#'   Character scalar, either \code{"previous"} or \code{"expand"}. The
#'   definition of the MTD, as in \code{\link{oc_3p3}}. Defaults to
#'   \code{"previous"}.
#'
#' @return
#'   An object of class \code{decision_table_3p3}, which is a data frame with one
#'   row per state of the current dose and the columns
#'   \item{stage}{\code{"escalation"} for dose escalation, or \code{"search"} for
#'     the search for the MTD, which only \code{mtd_rule = "expand"} has.}
#'   \item{n}{Number of patients treated at the current dose, 3 or 6.}
#'   \item{n_tox}{Number of patients with a DLT, from 0 to \code{n}.}
#'   \item{decision}{Decision code. During escalation, \code{"E"} escalates to
#'     the next higher dose, \code{"S"} treats three more patients at the current
#'     dose and \code{"STOP"} stops escalation, this dose and all higher doses
#'     being too toxic. During the search, \code{"S"} treats three more patients
#'     at the current dose, \code{"MTD"} selects it as the MTD and \code{"D"}
#'     moves the search to the next lower dose.}
#'   The definition of the MTD is stored in the attribute \code{mtd_rule}. The
#'   object has \code{print} and \code{plot} methods.
#'
#' @details
#'   During escalation three patients are treated at the current dose. No DLT
#'   escalates, one DLT adds three more patients at the same dose, and two or
#'   more stop escalation. With six patients, at most one DLT escalates and two or
#'   more stop escalation. Escalating from the highest dose also ends escalation.
#'
#'   With \code{mtd_rule = "previous"} the MTD is the dose below the one at which
#'   escalation stopped, or the highest dose when it was escalated from, and no
#'   MTD is selected when escalation stops at the starting dose. With
#'   \code{mtd_rule = "expand"} the search starts at the dose below the one at
#'   which escalation stopped, or at the highest dose when it was escalated from.
#'   A dose with three patients receives three more; a dose with at most one DLT
#'   out of six is selected as the MTD, and one with two or more passes the search
#'   to the next lower dose. No MTD is selected when the search moves below the
#'   starting dose.
#'
#'   Every number of DLTs from 0 to \code{n} is listed for six patients, although
#'   in a trial the six patients at a dose during escalation always include one
#'   DLT among the first three. The search examines only doses that were
#'   escalated from, so a dose with three patients has no DLT there, and only that
#'   state is listed.
#'
#'   The rules are those of \code{\link{oc_3p3}}: the tests follow this table
#'   through every possible course of a trial and obtain the operating
#'   characteristics of \code{oc_3p3()}.
#'
#' @examples
#' decision_table_3p3()
#'
#' # The other definition of the MTD adds the search for the MTD
#' decision_table_3p3(mtd_rule = "expand")
#'
#' @seealso \code{\link{print.decision_table_3p3}},
#'   \code{\link{plot.decision_table_3p3}}, \code{\link{oc_3p3}},
#'   \code{\link{boin_decision_table}}
#'
#' @export
decision_table_3p3 <- function(mtd_rule = c("previous", "expand")) {

  mtd_rule <- match.arg(mtd_rule)

  escalation <- data.frame(
    stage = "escalation",
    n = c(rep(3L, 4L), rep(6L, 7L)),
    n_tox = c(0:3, 0:6),
    stringsAsFactors = FALSE
  )
  escalation$decision <- ifelse(
    escalation$n_tox >= 2L, "STOP",
    ifelse(escalation$n == 3L & escalation$n_tox == 1L, "S", "E")
  )

  out <- escalation
  if (mtd_rule == "expand") {
    search <- data.frame(
      stage = "search",
      n = c(3L, rep(6L, 7L)),
      n_tox = c(0L, 0:6),
      stringsAsFactors = FALSE
    )
    search$decision <- ifelse(
      search$n == 3L, "S", ifelse(search$n_tox <= 1L, "MTD", "D")
    )
    out <- rbind(escalation, search)
  }
  rownames(out) <- NULL

  structure(out, class = c("decision_table_3p3", "data.frame"),
            mtd_rule = mtd_rule)
}
