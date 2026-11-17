#' Print the Operating Characteristics of the BOIN12 Design
#'
#' @description
#'   Display the summary produced by \code{\link{sim_boin12}} as a table with
#'   one column per dose.
#'
#' @param x
#'   An object of class \code{boin12_oc}.
#'
#' @param digits
#'   Integer scalar. Number of decimal places. Defaults to 1.
#'
#' @param percent
#'   Logical scalar. Express the average numbers of patients, DLTs and
#'   responses at each dose as percentages of the trial totals. Defaults to
#'   \code{FALSE}.
#'
#' @param kable
#'   Logical scalar. Print the table with \code{knitr::kable()}, for use in R
#'   Markdown. Defaults to \code{FALSE}.
#'
#' @param kable_format
#'   Character scalar. Format passed to \code{knitr::kable()}. Defaults to
#'   \code{"pipe"}.
#'
#' @param ...
#'   Not used.
#'
#' @return
#'   \code{x}, invisibly.
#'
#' @examples
#' oc <- sim_boin12(
#'   target = 0.35, target_resp = 0.25,
#'   p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
#'   p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38),
#'   n_cohort = 12, cohort_size = 3, n_trials = 200, seed = 1
#' )
#' print(oc)
#' print(oc, percent = TRUE)
#'
#' @export
print.boin12_oc <- function(x, digits = 1, percent = FALSE, kable = FALSE,
                            kable_format = "pipe", ...) {

  check_flag(percent, "percent")
  dose_names <- paste0("DL", seq_len(x$n_doses))

  pts <- x$n_pts_dose
  tox <- x$n_tox_dose
  resp <- x$n_resp_dose
  labels <- c("Patients treated", "Patients with DLT", "Patients with response")
  if (percent) {
    pts <- pts / x$total_n_pts * 100
    tox <- tox / x$total_n_tox * 100
    resp <- resp / x$total_n_resp * 100
    labels <- paste(labels, "(%)")
  }

  tab <- rbind(
    c(x$p_true * 100, NA_real_),
    c(x$p_resp * 100, NA_real_),
    c(x$sel_percent, x$percent_no_obd),
    c(pts, x$total_n_pts),
    c(tox, x$total_n_tox),
    c(resp, x$total_n_resp)
  )
  rownames(tab) <- c("True DLT rate (%)", "True response rate (%)",
                     "OBD selected (%)", labels)
  colnames(tab) <- c(dose_names, "Total / No OBD")
  tab <- round(tab, digits)

  if (kable) {
    if (!requireNamespace("knitr", quietly = TRUE)) {
      stop("Package 'knitr' is needed when 'kable' is TRUE", call. = FALSE)
    }
    print(knitr::kable(tab, format = kable_format, row.names = TRUE))
    return(invisible(x))
  }

  set <- x$settings
  design <- if (is.null(x$design)) "BOIN12" else x$design
  cat(design, " operating characteristics\n", sep = "")
  cat("  DLT limit            : ", format(x$target), "\n", sep = "")
  cat("  response limit       : ", format(x$target_resp), "\n", sep = "")
  cat("  trials               : ", x$n_trials, "\n", sep = "")
  cat("  cohorts              : ", set$n_cohort, " of size ",
      paste(unique(set$cohort_size), collapse = "/"), "\n", sep = "")
  cat("  max sample size      : ", set$max_total_pts, "\n", sep = "")
  cat("  lambda_e / lambda_d  : ", format(round(x$lambda_e, 4)), " / ",
      format(round(x$lambda_d, 4)), "\n", sep = "")
  cat("  utility scores       : 100, ", format(set$psi00), ", ",
      format(set$psi11), ", 0; benchmark ", format(round(x$u_benchmark, 2)),
      "\n", sep = "")
  cat("  ties / no candidate  : ", set$tie_break, " / ", set$no_candidate, "\n",
      sep = "")
  cat("  exploration from     : ", set$explore_from, "\n", sep = "")
  cat("  OBD selection        : utility by ", x$est_method, "\n", sep = "")
  cat("\n")

  print(tab, na.print = "")

  cat("\nTrials stopped early: ", format(round(x$percent_early_stop, digits)),
      "%\n", sep = "")
  if (!is.null(x$duration)) {
    cat("Average trial duration: ", format(round(x$duration, digits)), "\n",
        sep = "")
  }

  invisible(x)
}
