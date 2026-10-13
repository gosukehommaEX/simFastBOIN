#' Print the Next Dose Combination of a BOIN Drug-Combination Trial
#'
#' @description
#'   Display the decision returned by \code{\link{comb_boin_next}}: the data at
#'   the current combination, the boundaries that apply to it, the decision,
#'   the candidates that were compared and the combinations eliminated so far.
#'
#' @param x
#'   An object of class \code{comb_boin_next}.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   The object \code{x}, invisibly.
#'
#' @examples
#' n_pts <- matrix(c(3, 0, 0,
#'                   0, 0, 0), nrow = 2, byrow = TRUE)
#' n_tox <- matrix(0, nrow = 2, ncol = 3)
#'
#' print(comb_boin_next(target = 0.30, n_pts = n_pts, n_tox = n_tox,
#'                      current = c(1, 1)))
#'
#' @export
print.comb_boin_next <- function(x, ...) {

  name <- function(dose_a, dose_b) paste0("A", dose_a, "B", dose_b)
  a <- x$current[["dose_a"]]
  b <- x$current[["dose_b"]]
  n_cur <- x$n_pts[a, b]
  y_cur <- x$n_tox[a, b]
  bd <- x$boundary

  cat("BOIN drug-combination design: dose for the next cohort\n")
  cat("  target DLT rate : ", format(x$target), "\n", sep = "")
  cat("  current         : ", name(a, b), ", ", y_cur, " DLT",
      if (y_cur == 1L) "" else "s", " out of ", n_cur, " patient",
      if (n_cur == 1L) "" else "s", " (", format(round(y_cur / n_cur, 3)),
      ")\n", sep = "")
  cat("  boundaries      : escalate if DLTs <= ", bd$escalate_at_most,
      ", de-escalate if DLTs >= ", bd$deescalate_at_least,
      if (is.na(bd$eliminate_at_least)) {
        ", no elimination yet"
      } else {
        paste0(", eliminate if DLTs >= ", bd$eliminate_at_least)
      },
      " (", bd$n, " patients)\n", sep = "")
  cat("\n")

  targets <- name(x$next_dose[, "dose_a"], x$next_dose[, "dose_b"])
  sentence <- switch(
    x$decision,
    escalate = paste("escalate to", paste(targets, collapse = " or ")),
    `de-escalate` = paste("de-escalate to", paste(targets, collapse = " or ")),
    stay = paste("stay at", name(a, b)),
    stop = switch(
      x$stop_reason,
      lowest_dose_eliminated = paste("stop the trial: the lowest combination",
                                     "is eliminated, and no combination is",
                                     "selected"),
      lowest_dose_too_toxic = paste("stop the trial: the lowest combination",
                                    "is too toxic, and no combination is",
                                    "selected"),
      n_earlystop = paste0("stop the trial: ", n_cur, " patients have been ",
                           "treated at ", name(a, b), "; select the MTD ",
                           "with comb_boin_select_mtd()")
    )
  )
  cat("Decision: ", sentence, "\n", sep = "")
  if (x$tie) {
    cat("  The two candidates tie; the design chooses one of them at random.\n")
  }
  if (!is.na(x$note)) cat("  ", x$note, "\n", sep = "")

  if (!is.null(x$candidates)) {
    cand <- x$candidates
    shown <- data.frame(
      Direction = cand$direction,
      Combination = ifelse(is.na(cand$dose_a), "-",
                           name(cand$dose_a, cand$dose_b)),
      Patients = ifelse(is.na(cand$n_pts), "-", cand$n_pts),
      DLTs = ifelse(is.na(cand$n_tox), "-", cand$n_tox),
      "Pr(in interval)" = ifelse(is.na(cand$prob_in_interval), "-",
                                 format(round(cand$prob_in_interval, 3),
                                        nsmall = 3)),
      Score = ifelse(is.na(cand$score), "-",
                     format(round(cand$score, 4), nsmall = 4)),
      Available = ifelse(cand$available, "yes", "no"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    cat("\nCandidates (interval ", format(round(x$lambda_e, 4)), " to ",
        format(round(x$lambda_d, 4)), "):\n", sep = "")
    print(shown, row.names = FALSE, right = TRUE)
  }

  gone <- which(x$eliminated, arr.ind = TRUE)
  cat("\nEliminated combinations: ",
      if (nrow(gone) == 0L) {
        "none"
      } else {
        gone <- gone[order(gone[, 1L], gone[, 2L]), , drop = FALSE]
        paste(name(gone[, 1L], gone[, 2L]), collapse = ", ")
      },
      "\n", sep = "")

  invisible(x)
}
