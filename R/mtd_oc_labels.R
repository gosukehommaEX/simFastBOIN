#' Row and Line Labels of the Printed Results of MTD-Finding Designs
#'
#' @description
#'   Internal helper of the print methods of \code{sim_boin()} and of the
#'   designs that share its result layout. The gBOIN design reports a
#'   toxicity score in place of the DLT, so its tables are labelled
#'   differently; every other design keeps the labels of the BOIN design.
#'
#' @param design Name of the design stored in the result, or \code{NULL}.
#'
#' @return A list of labels: \code{target} (the line of the target),
#'   \code{rate} (the row of the true values), \code{tox} and \code{tox_pct}
#'   (the row of the toxicity, in numbers and in percent), \code{trial_tox}
#'   (the row of the toxicity of simulated trials) and \code{above} (the
#'   start of the line of the overdose summary).
#'
#' @noRd
mtd_oc_labels <- function(design = NULL) {
  if (identical(design, "gBOIN")) {
    list(
      target = "  target score    : ",
      rate = "True score (%)",
      tox = "Toxicity score",
      tox_pct = "Toxicity score (%)",
      trial_tox = "Toxicity score",
      above = "Doses above a true score of "
    )
  } else {
    list(
      target = "  target DLT rate : ",
      rate = "True DLT rate (%)",
      tox = "Patients with DLT",
      tox_pct = "Patients with DLT (%)",
      trial_tox = "DLTs observed",
      above = "Doses above a true DLT rate of "
    )
  }
}
