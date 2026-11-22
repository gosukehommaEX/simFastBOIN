#' Labels of the Print Methods of the BOIN-ET Family
#'
#' @description
#'   Texts that the print methods of the BOIN-ET results show, by design: for
#'   the gBOIN-ET design the outcomes are scores, not DLTs and responses.
#'   Internal.
#'
#' @param design
#'   Character scalar or \code{NULL}. Name of the design.
#'
#' @return
#'   A list of character scalars.
#'
#' @keywords internal
#' @noRd
et_oc_labels <- function(design = NULL) {
  if (identical(design, "gBOIN-ET")) {
    list(
      target = "  target toxicity score: ",
      target_resp = "  target efficacy score: ",
      rate = "True toxicity score (%)",
      rate_resp = "True efficacy score (%)",
      tox = "Toxicity score",
      resp = "Efficacy score",
      trial_tox = "Toxicity score",
      trial_resp = "Efficacy score",
      estimate = " efficacy scores"
    )
  } else {
    list(
      target = "  target DLT rate      : ",
      target_resp = "  target response rate : ",
      rate = "True DLT rate (%)",
      rate_resp = "True response rate (%)",
      tox = "Patients with DLT",
      resp = "Patients with response",
      trial_tox = "DLTs observed",
      trial_resp = "Responses observed",
      estimate = " response rates"
    )
  }
}
