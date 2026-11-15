#' Print the Operating Characteristics of the BOIN-ET Design
#'
#' @description
#'   Display the summary produced by \code{\link{sim_boin_et}} as a table with
#'   one column per dose.
#'
#' @param x
#'   An object of class \code{boin_et_oc}.
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
#'   Further arguments, currently ignored.
#'
#' @return
#'   The object \code{x}, invisibly.
#'
#' @examples
#' oc <- sim_boin_et(
#'   target = 0.30, target_resp = 0.60,
#'   p_true = c(0.05, 0.15, 0.25, 0.30, 0.35, 0.40),
#'   p_resp = c(0.05, 0.35, 0.55, 0.65, 0.75, 0.70),
#'   n_cohort = 12, cohort_size = 3, n_trials = 500, seed = 1
#' )
#' print(oc)
#' print(oc, percent = TRUE)
#'
#' @seealso \code{\link{sim_boin_et}}
#'
#' @export
print.boin_et_oc <- function(x, digits = 1, percent = FALSE, kable = FALSE,
                             kable_format = "pipe", ...) {

  if (!is.logical(percent) || length(percent) != 1L || is.na(percent)) {
    stop("'percent' must be TRUE or FALSE", call. = FALSE)
  }

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
  design <- if (is.null(x$design)) "BOIN-ET" else x$design
  cat(design, " operating characteristics\n", sep = "")
  cat("  target DLT rate      : ", format(x$target), "\n", sep = "")
  cat("  target response rate : ", format(x$target_resp), "\n", sep = "")
  cat("  trials               : ", x$n_trials, "\n", sep = "")
  cat("  cohorts              : ", set$n_cohort, " of size ",
      paste(unique(set$cohort_size), collapse = "/"), "\n", sep = "")
  cat("  max sample size      : ", set$max_total_pts, "\n", sep = "")
  cat("  lambda1 / lambda2    : ", format(round(x$lambda1, 4)), " / ",
      format(round(x$lambda2, 4)), "\n", sep = "")
  cat("  eta1                 : ", format(round(x$eta1, 4)), "\n", sep = "")
  cat("  OBD selection        : ", x$obd_method, " with ", x$est_method,
      " response rates\n", sep = "")
  if (isTRUE(set$boinet_compatible)) {
    cat("  rules                : those of boinet 1.6.0\n")
  }
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
