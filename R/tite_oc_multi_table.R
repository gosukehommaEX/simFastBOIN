#' Assemble the Cross-Scenario Summary Table of TITE-BOIN Designs
#'
#' @description
#'   Build the data frame shown by \code{\link{print.tite_boin_oc_multi}}. Each
#'   scenario has the four rows of the BOIN table, built by the internal
#'   \code{oc_multi_table()}, followed by the average trial duration and the
#'   percentage of trials with accrual suspended, both in the last column.
#'   Internal.
#'
#' @param results
#'   List of \code{tite_boin_oc} objects, one per scenario.
#'
#' @param scenario_names
#'   Character vector of scenario names.
#'
#' @param n_doses
#'   Integer scalar. Number of dose levels.
#'
#' @param percent
#'   Logical scalar. Express the average number of patients and DLTs at each dose
#'   as a percentage of the trial total instead of a count. Defaults to
#'   \code{FALSE}.
#'
#' @param digits
#'   Integer scalar. Number of decimal places. Defaults to 1.
#'
#' @return
#'   A data frame with six rows per scenario.
#'
#' @keywords internal
#' @noRd
tite_oc_multi_table <- function(results, scenario_names, n_doses,
                                percent = FALSE, digits = 1) {

  base <- oc_multi_table(results, scenario_names, n_doses, percent = percent,
                         digits = digits)
  dose_names <- paste0("DL", seq_len(n_doses))

  blocks <- lapply(seq_along(results), function(i) {
    res <- results[[i]]
    extra <- data.frame(
      Scenario = c("", ""),
      Item = c("Trial duration (mean)", "Accrual suspended (% trials)"),
      stringsAsFactors = FALSE
    )
    for (dose in dose_names) extra[[dose]] <- NA_real_
    extra[["Total / No MTD"]] <- c(round(res$duration_mean, digits),
                                   round(res$pct_trials_suspended, digits))
    rbind(base[(4L * i - 3L):(4L * i), , drop = FALSE], extra)
  })

  out <- do.call(rbind, blocks)
  rownames(out) <- NULL
  out
}
