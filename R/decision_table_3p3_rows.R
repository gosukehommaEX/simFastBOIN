#' Rows of a 3+3 Decision Table as Printed
#'
#' @description
#'   Internal helper of \code{\link{print.decision_table_3p3}}. Collapse the
#'   consecutive numbers of DLTs that share a decision into one row, within each
#'   stage and number of patients.
#'
#' @param x
#'   An object of class \code{decision_table_3p3}.
#'
#' @return
#'   A data frame of character columns \code{Stage}, \code{Patients},
#'   \code{DLTs} and \code{Decision}.
#'
#' @noRd
decision_table_3p3_rows <- function(x) {

  stage_labels <- c(escalation = "Dose escalation", search = "Search for the MTD")
  decision_labels <- c(
    E = "Escalate",
    S = "Treat 3 more patients at this dose",
    STOP = "Stop escalation",
    MTD = "Select this dose as the MTD",
    D = "Move to the next lower dose"
  )

  rows <- list()
  for (stage in intersect(names(stage_labels), unique(x$stage))) {
    for (k in sort(unique(x$n[x$stage == stage]))) {
      idx <- which(x$stage == stage & x$n == k)
      idx <- idx[order(x$n_tox[idx])]
      tox <- x$n_tox[idx]
      dec <- x$decision[idx]
      new_run <- c(TRUE, dec[-1L] != dec[-length(dec)] | diff(tox) != 1L)
      run_id <- cumsum(new_run)
      for (r in unique(run_id)) {
        lo <- min(tox[run_id == r])
        hi <- max(tox[run_id == r])
        tox_label <- if (lo == hi) {
          as.character(lo)
        } else if (lo == 0L && hi == k) {
          paste0(lo, "-", hi)
        } else if (hi == k) {
          paste(">=", lo)
        } else if (lo == 0L) {
          paste("<=", hi)
        } else {
          paste0(lo, "-", hi)
        }
        rows[[length(rows) + 1L]] <- c(stage_labels[[stage]], as.character(k),
                                       tox_label,
                                       decision_labels[[dec[run_id == r][1L]]])
      }
    }
  }

  columns <- c("Stage", "Patients", "DLTs", "Decision")
  if (length(rows) == 0L) {
    out <- matrix(character(0), nrow = 0L, ncol = length(columns))
  } else {
    out <- do.call(rbind, rows)
  }
  colnames(out) <- columns
  as.data.frame(out, stringsAsFactors = FALSE)
}
