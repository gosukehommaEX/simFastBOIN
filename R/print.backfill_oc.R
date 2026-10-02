#' Print Operating Characteristics of a BOIN Design with Backfilling
#'
#' @description
#'   Display the result of \code{\link{sim_bf_boin}} or
#'   \code{\link{sim_be_boin}}: the design, a table by dose of the true DLT and
#'   response rates, the MTD selection, the patients treated and backfilled and
#'   the DLTs, then the exposure above the overdose cutoff and the timing.
#'
#' @param x
#'   An object of class \code{backfill_oc}.
#'
#' @param digits
#'   Integer scalar. Number of decimal places. Defaults to 1.
#'
#' @param percent
#'   Logical scalar. Show patients, backfilled patients and DLTs as percentages
#'   of their totals rather than as averages per trial. Defaults to
#'   \code{FALSE}.
#'
#' @param kable
#'   Logical scalar. Print only the table, formatted by \code{knitr::kable()}.
#'   Defaults to \code{FALSE}.
#'
#' @param kable_format
#'   Character scalar passed to \code{knitr::kable()} as \code{format}.
#'   Defaults to \code{"pipe"}.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   The object \code{x}, invisibly.
#'
#' @examples
#' oc <- sim_bf_boin(
#'   target = 0.25,
#'   p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
#'   p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58),
#'   n_cohort = 10, cohort_size = 3, window = 1, accrual_rate = 3,
#'   n_trials = 200, seed = 1
#' )
#' print(oc)
#' print(oc, percent = TRUE)
#'
#' @seealso \code{\link{sim_bf_boin}}, \code{\link{sim_be_boin}}
#'
#' @export
print.backfill_oc <- function(x, digits = 1, percent = FALSE, kable = FALSE,
                              kable_format = "pipe", ...) {

  if (!is.logical(percent) || length(percent) != 1L || is.na(percent)) {
    stop("'percent' must be TRUE or FALSE", call. = FALSE)
  }

  dose_names <- paste0("DL", seq_len(x$n_doses))

  pts <- x$n_pts_dose
  bf <- x$n_bf_dose
  tox <- x$n_tox_dose
  bf_total <- x$total_n_bf
  pts_label <- "Patients treated"
  bf_label <- "  of which backfilled"
  tox_label <- "Patients with DLT"
  if (percent) {
    pts <- pts / x$total_n_pts * 100
    bf <- bf / x$total_n_pts * 100
    bf_total <- bf_total / x$total_n_pts * 100
    tox <- tox / x$total_n_tox * 100
    pts_label <- "Patients treated (%)"
    bf_label <- "  of which backfilled (%)"
    tox_label <- "Patients with DLT (%)"
  }

  tab <- rbind(
    c(x$p_true * 100, NA_real_),
    c(x$p_resp * 100, NA_real_),
    c(x$sel_percent, x$percent_no_mtd),
    c(pts, x$total_n_pts),
    c(bf, bf_total),
    c(tox, x$total_n_tox)
  )
  rownames(tab) <- c("True DLT rate (%)", "True response rate (%)",
                     "MTD selected (%)", pts_label, bf_label, tox_label)
  colnames(tab) <- c(dose_names, "Total / No MTD")
  tab <- round(tab, digits)

  if (kable) {
    if (!requireNamespace("knitr", quietly = TRUE)) {
      stop("Package 'knitr' is needed when 'kable' is TRUE", call. = FALSE)
    }
    print(knitr::kable(tab, format = kable_format, row.names = TRUE))
    return(invisible(x))
  }

  set <- x$settings
  if (identical(set$design, "be_boin")) {
    cat("BE-BOIN operating characteristics (Chen et al., 2026)\n")
  } else {
    cat("BF-BOIN operating characteristics (Zhao et al., 2024)\n")
  }
  cat("  target DLT rate : ", format(x$target), "\n", sep = "")
  cat("  trials          : ", x$n_trials, "\n", sep = "")
  cat("  cohorts         : ", set$n_cohort, " of size ",
      paste(unique(set$cohort_size), collapse = "/"), "\n", sep = "")
  cat("  escalation size : ", set$max_total_pts, "\n", sep = "")
  cat("  backfill cap    : ", set$n_cap, " patients per dose\n", sep = "")
  cat("  lambda_e / _d   : ", format(round(x$lambda_e, 4)), " / ",
      format(round(x$lambda_d, 4)), "\n", sep = "")
  cat("  window          : ", format(set$window), " (DLT), ",
      format(set$resp_window), " (response)\n", sep = "")
  cat("  accrual         : ", format(set$accrual_rate),
      " patients per unit of time (", set$accrual, ")\n", sep = "")
  cat("\n")

  print(tab, na.print = "")

  overdose <- x$overdose
  cat("\nDoses above a true DLT rate of ", format(overdose$cutoff), ": ",
      if (length(overdose$doses) > 0L) {
        paste0("DL", overdose$doses, collapse = ", ")
      } else {
        "none"
      },
      "\n", sep = "")
  cat("  Patients treated there              : ",
      format(round(overdose$pct_patients, digits)), "%\n", sep = "")
  cat("  Trials treating any patient there   : ",
      format(round(overdose$pct_trials_any, digits)), "%\n", sep = "")
  cat("  Trials selecting an MTD there       : ",
      format(round(overdose$pct_trials_mtd_above, digits)), "%\n", sep = "")

  cat("\nBackfilling and timing\n")
  cat("  Trials with backfilled patients     : ",
      format(round(x$pct_trials_backfill, digits)), "%\n", sep = "")
  cat("  Average trial duration              : ",
      format(round(x$duration_mean, digits)), "\n", sep = "")
  cat("  Trials in which the escalation waited: ",
      format(round(x$pct_trials_suspended, digits)), "%\n", sep = "")
  cat("  Average time the escalation waited  : ",
      format(round(x$avg_time_suspended, digits)), "\n", sep = "")
  if (identical(set$no_slot, "leave")) {
    cat("  Average patients turned away        : ",
        format(round(x$avg_n_turned_away, digits)), "\n", sep = "")
  }

  invisible(x)
}
