#' Rows of a TITE-BOIN Decision Table as Printed
#'
#' @description
#'   Internal helper of \code{\link{print.tite_boin_decision_table}}. Collapse
#'   the states of a decision table into the rows of the published tables:
#'   consecutive numbers of pending patients with the same decision share a row,
#'   and consecutive numbers of DLTs with the same decision throughout share a
#'   row as well.
#'
#' @param x
#'   An object of class \code{tite_boin_decision_table}.
#'
#' @param cohort_size
#'   Integer scalar or \code{NULL}. When supplied, only the sample sizes that are
#'   multiples of \code{cohort_size} are kept.
#'
#' @param digits
#'   Integer scalar. Number of decimal places of the boundaries.
#'
#' @return
#'   A data frame of character columns \code{Patients}, \code{DLTs},
#'   \code{Pending}, \code{Escalate}, \code{Stay} and \code{De-escalate}.
#'
#' @noRd
tite_boin_display_rows <- function(x, cohort_size = NULL, digits = 2L) {

  n_values <- sort(unique(x$n))
  if (!is.null(cohort_size)) {
    n_values <- n_values[n_values %% cohort_size == 0L]
  }

  fmt <- function(v) formatC(v, format = "f", digits = digits)
  esc_text <- ifelse(is.na(x$esc_bound), "", fmt(x$esc_bound))
  deesc_text <- ifelse(is.na(x$deesc_bound), "", fmt(x$deesc_bound))
  signature <- paste(x$decision, esc_text, deesc_text, sep = "|")

  # Text of the Escalate, Stay and De-escalate columns for one state.
  decision_text <- function(i) {
    decision <- x$decision[i]
    e <- esc_text[i]
    d <- deesc_text[i]
    parts <- strsplit(decision, "/", fixed = TRUE)[[1L]]
    if (length(parts) == 1L) {
      return(switch(EXPR = decision,
                    E = c("Y", "", ""),
                    S = c("", "Y", ""),
                    D = c("", "", "Y"),
                    DE = c("", "", "Y & Elim"),
                    SUS = c("Suspend", "Suspend", "Suspend")))
    }
    has_esc <- parts[1L] %in% c("E", "SUS")
    has_deesc <- parts[length(parts)] == "D"
    esc <- if (!has_esc) {
      ""
    } else if (parts[1L] == "E") {
      paste(">=", e)
    } else {
      paste("Suspend if >=", e)
    }
    stay <- if (has_esc && has_deesc) {
      paste0("> ", d, ", < ", e)
    } else if (has_esc) {
      paste("<", e)
    } else {
      paste(">", d)
    }
    deesc <- if (has_deesc) paste("<=", d) else ""
    c(esc, stay, deesc)
  }

  rows <- list()
  for (k in n_values) {
    in_n <- which(x$n == k)
    tox_values <- sort(unique(x$n_tox[in_n]))

    # For each number of DLTs, runs of consecutive pending counts that share a
    # decision.
    groups <- lapply(tox_values, function(s) {
      idx <- in_n[x$n_tox[in_n] == s]
      idx <- idx[order(x$n_pending[idx])]
      pend <- x$n_pending[idx]
      sig <- signature[idx]
      new_run <- c(TRUE, sig[-1L] != sig[-length(sig)] | diff(pend) != 1L)
      run_id <- cumsum(new_run)
      list(tox = s,
           first = idx[!duplicated(run_id)],
           last = idx[!duplicated(run_id, fromLast = TRUE)])
    })

    # A number of DLTs whose single run covers every attainable pending count
    # can share a row with its neighbors.
    full <- vapply(groups, function(g) {
      length(g$first) == 1L && x$n_pending[g$first] == 0L &&
        x$n_pending[g$last] == k - g$tox
    }, logical(1))

    i <- 1L
    while (i <= length(groups)) {
      g <- groups[[i]]
      j <- i
      if (full[i]) {
        while (j < length(groups) && full[j + 1L] &&
               groups[[j + 1L]]$tox == groups[[j]]$tox + 1L &&
               signature[groups[[j + 1L]]$first] == signature[g$first]) {
          j <- j + 1L
        }
      }

      if (j > i) {
        lo <- g$tox
        hi <- groups[[j]]$tox
        tox_label <- if (hi == k) {
          paste(">=", lo)
        } else if (hi == lo + 1L) {
          paste0(lo, ", ", hi)
        } else {
          paste0(lo, "-", hi)
        }
        pending_max <- k - lo
        pending_label <- if (pending_max > 0L) paste("<=", pending_max) else "0"
        rows[[length(rows) + 1L]] <- c(as.character(k), tox_label, pending_label,
                                       decision_text(g$first))
      } else {
        pending_max <- k - g$tox
        for (m in seq_along(g$first)) {
          c0 <- x$n_pending[g$first[m]]
          c1 <- x$n_pending[g$last[m]]
          pending_label <- if (c0 == 0L && c1 > 0L) {
            paste("<=", c1)
          } else if (c1 == pending_max && c0 > 0L) {
            paste(">=", c0)
          } else if (c0 == c1) {
            as.character(c0)
          } else {
            paste0(c0, "-", c1)
          }
          rows[[length(rows) + 1L]] <- c(as.character(k), as.character(g$tox),
                                         pending_label, decision_text(g$first[m]))
        }
      }
      i <- j + 1L
    }
  }

  columns <- c("Patients", "DLTs", "Pending", "Escalate", "Stay", "De-escalate")
  if (length(rows) == 0L) {
    out <- matrix(character(0), nrow = 0L, ncol = length(columns))
  } else {
    out <- do.call(rbind, rows)
  }
  colnames(out) <- columns
  as.data.frame(out, stringsAsFactors = FALSE, optional = TRUE)
}
