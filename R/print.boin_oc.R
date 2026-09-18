#' Print Operating Characteristics of a BOIN Design
#'
#' @description
#'   Display the summary table produced by \code{\link{sim_boin}}.
#'
#' @param x
#'   An object of class \code{boin_oc}.
#'
#' @param digits
#'   Integer scalar. Number of decimal places. Defaults to 1.
#'
#' @param percent
#'   Logical scalar. Express the average number of patients and DLTs at each dose
#'   as a percentage of the trial total instead of a count. Defaults to
#'   \code{FALSE}.
#'
#' @param kable
#'   Logical scalar. Return the table through \code{knitr::kable()} rather than
#'   printing it directly, which is convenient inside a report. Requires the
#'   \pkg{knitr} package. Defaults to \code{FALSE}.
#'
#' @param kable_format
#'   Character scalar passed to \code{knitr::kable()}, for example \code{"pipe"},
#'   \code{"html"} or \code{"latex"}. Defaults to \code{"pipe"}.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   The object \code{x}, invisibly.
#'
#' @examples
#' oc <- sim_boin(
#'   target = 0.30,
#'   p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   n_trials = 200,
#'   seed = 123
#' )
#'
#' print(oc)
#' print(oc, percent = TRUE)
#'
#' @export
print.boin_oc <- function(x, digits = 1, percent = FALSE, kable = FALSE,
                          kable_format = "pipe", ...) {

  if (!is.logical(percent) || length(percent) != 1L || is.na(percent)) {
    stop("'percent' must be TRUE or FALSE", call. = FALSE)
  }

  dose_names <- paste0("DL", seq_len(x$n_doses))

  pts <- x$n_pts_dose
  tox <- x$n_tox_dose
  pts_label <- "Patients treated"
  tox_label <- "Patients with DLT"
  if (percent) {
    pts <- pts / x$total_n_pts * 100
    tox <- tox / x$total_n_tox * 100
    pts_label <- "Patients treated (%)"
    tox_label <- "Patients with DLT (%)"
  }

  tab <- rbind(
    c(x$p_true * 100, NA_real_),
    c(x$sel_percent, x$percent_no_mtd),
    c(pts, x$total_n_pts),
    c(tox, x$total_n_tox)
  )
  rownames(tab) <- c("True DLT rate (%)", "MTD selected (%)", pts_label, tox_label)
  colnames(tab) <- c(dose_names, "Total / No MTD")
  tab <- round(tab, digits)

  if (kable) {
    if (!requireNamespace("knitr", quietly = TRUE)) {
      stop("Package 'knitr' is needed when 'kable' is TRUE", call. = FALSE)
    }
    out <- knitr::kable(tab, format = kable_format, row.names = TRUE)
    print(out)
    return(invisible(x))
  }

  cat("BOIN operating characteristics\n")
  cat("  target DLT rate : ", format(x$target), "\n", sep = "")
  cat("  trials          : ", x$n_trials, "\n", sep = "")
  cat("  cohorts         : ", x$settings$n_cohort, " of size ",
      paste(unique(x$settings$cohort_size), collapse = "/"), "\n", sep = "")
  cat("  max sample size : ", x$settings$max_total_pts, "\n", sep = "")
  cat("  lambda_e / _d   : ", format(round(x$lambda_e, 4)), " / ",
      format(round(x$lambda_d, 4)), "\n", sep = "")
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
  cat("  Trials treating over 60% there      : ",
      format(round(overdose$pct_trials_over_60, digits)), "%\n", sep = "")
  cat("  Trials treating over 80% there      : ",
      format(round(overdose$pct_trials_over_80, digits)), "%\n", sep = "")
  cat("  Trials selecting an MTD there       : ",
      format(round(overdose$pct_trials_mtd_above, digits)), "%\n", sep = "")

  invisible(x)
}
