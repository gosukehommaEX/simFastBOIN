#' Print Operating Characteristics of the 3+3 Design
#'
#' @description
#'   Display the summary table produced by \code{\link{oc_3p3}} or
#'   \code{\link{sim_3p3}}. The layout matches
#'   \code{\link{print.boin_oc}}, so a 3+3 table can be read beside a BOIN one.
#'
#' @param x
#'   An object of class \code{oc_3p3}.
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
#' oc <- oc_3p3(p_true = c(0.30, 0.48, 0.67))
#'
#' print(oc)
#' print(oc, percent = TRUE)
#'
#' @seealso \code{\link{oc_3p3}}, \code{\link{sim_3p3}}
#'
#' @export
print.oc_3p3 <- function(x, digits = 1, percent = FALSE, kable = FALSE,
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
    tox <- if (x$total_n_tox > 0) tox / x$total_n_tox * 100 else tox
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
    print(knitr::kable(tab, format = kable_format, row.names = TRUE))
    return(invisible(x))
  }

  rule_text <- if (x$mtd_rule == "previous") {
    "dose below the toxic one"
  } else {
    "at most one DLT in six"
  }

  cat("3+3 operating characteristics\n")
  cat("  MTD rule        : ", x$mtd_rule, " (", rule_text, ")\n", sep = "")
  cat("  obtained by     : ",
      if (identical(x$method, "exact")) {
        "exact enumeration"
      } else {
        paste0(x$n_trials, " simulated trials")
      },
      "\n", sep = "")
  cat("  start dose      : DL", x$start_dose, "\n", sep = "")
  cat("  max sample size : ", 6L * (x$n_doses - x$start_dose + 1L),
      " (six per dose)\n", sep = "")
  cat("\n")

  print(tab, na.print = "")

  overdose <- x$overdose
  cat("\nDoses above a true DLT rate of ", format(round(overdose$cutoff, 4)), ": ",
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

  invisible(x)
}
