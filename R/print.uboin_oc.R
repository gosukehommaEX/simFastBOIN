#' Print Operating Characteristics of the U-BOIN Design
#'
#' @description
#'   Print the design, a table of the true rates and utilities, the selection
#'   percentages and the average numbers of patients, DLTs and responses of
#'   every dose, and the reasons the trials ended, for a \code{uboin_oc}
#'   object.
#'
#' @param x
#'   An object of class \code{uboin_oc}, as returned by \code{\link{sim_uboin}}.
#'
#' @param digits
#'   Integer scalar. Number of decimal places. Defaults to 1.
#'
#' @param percent
#'   Logical scalar. Show the numbers of patients, DLTs and responses as
#'   percentages of their totals. Defaults to \code{FALSE}.
#'
#' @param kable
#'   Logical scalar. Print the table with \code{knitr::kable()} instead of the
#'   plain console layout. Defaults to \code{FALSE}.
#'
#' @param kable_format
#'   Character scalar. Format passed to \code{knitr::kable()}.
#'   Defaults to \code{"pipe"}.
#'
#' @param ...
#'   Not used.
#'
#' @return
#'   \code{x}, invisibly.
#'
#' @examples
#' oc <- sim_uboin(
#'   tox_upper = 0.30,
#'   resp_lower = 0.20,
#'   p_true = c(0.02, 0.15, 0.30, 0.45, 0.60),
#'   p_resp = c(0.20, 0.65, 0.65, 0.65, 0.65),
#'   n_cohort = 18,
#'   cohort_size = 3,
#'   resp_assoc = 0.2,
#'   n_trials = 200,
#'   seed = 123
#' )
#'
#' print(oc)
#' print(oc, percent = TRUE)
#'
#' @export
print.uboin_oc <- function(x, digits = 1, percent = FALSE, kable = FALSE,
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
    c(x$true_utility, NA_real_),
    c(x$sel_percent, x$percent_no_obd),
    c(pts, x$total_n_pts),
    c(tox, x$total_n_tox),
    c(resp, x$total_n_resp)
  )
  rownames(tab) <- c("True DLT rate (%)", "True response rate (%)",
                     "True utility", "OBD selected (%)", labels)
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
  cat("U-BOIN operating characteristics\n")
  cat("  DLT upper limit      : ", format(x$tox_upper), "\n", sep = "")
  cat("  response lower limit : ", format(x$resp_lower), "\n", sep = "")
  cat("  stage I target       : ", format(x$target), "\n", sep = "")
  cat("  trials               : ", x$n_trials, "\n", sep = "")
  cat("  cohorts              : ", set$n_cohort, " of size ", set$cohort_size,
      "\n", sep = "")
  cat("  max sample size      : ", set$max_total_pts, "\n", sep = "")
  cat("  stage II from / stop : ", set$n_stage1, " / ", set$n_stop,
      " patients at a dose\n", sep = "")
  cat("  lambda_e / lambda_d  : ", format(round(x$lambda_e, 4)), " / ",
      format(round(x$lambda_d, 4)), "\n", sep = "")
  utilities <- paste(format(x$utility, trim = TRUE), collapse = ", ")
  cat("  utilities            : ", utilities, "\n", sep = "")
  cat("  cutoffs (DLT / resp) : ", format(set$cutoff_eli), " / ",
      format(set$cutoff_eli_resp), "\n", sep = "")
  cat("\n")

  print(tab, na.print = "")

  cat("\nTrials reaching stage II: ", format(round(x$percent_stage2, digits)),
      "%\n", sep = "")
  cat("Reason for ending (%):\n")
  print(round(x$stop_reason_percent, digits))

  invisible(x)
}
