#' Summary Table of Backfill Designs Across Scenarios
#'
#' @description
#'   Internal helper of \code{\link{sim_bf_boin_multi}},
#'   \code{\link{sim_be_boin_multi}} and their print method. Extend the table of
#'   \code{\link{sim_boin_multi}} with the true response rates, the backfilled
#'   patients and the average trial duration of every scenario.
#'
#' @param results
#'   List of \code{backfill_oc} objects, one per scenario.
#'
#' @param scenario_names
#'   Character vector of scenario names.
#'
#' @param n_doses
#'   Integer scalar. Number of dose levels.
#'
#' @param percent
#'   Logical scalar. Show patients, backfilled patients and DLTs as percentages
#'   of their totals.
#'
#' @param digits
#'   Integer scalar. Number of decimal places.
#'
#' @return
#'   A data frame with columns \code{Scenario}, \code{Item}, one per dose and
#'   \code{Total / No MTD}, and seven rows per scenario.
#'
#' @noRd
backfill_oc_multi_table <- function(results, scenario_names, n_doses,
                                    percent = FALSE, digits = 1) {

  base <- oc_multi_table(results, scenario_names, n_doses, percent = percent,
                         digits = digits)
  dose_names <- paste0("DL", seq_len(n_doses))

  blocks <- lapply(seq_along(results), function(i) {
    res <- results[[i]]
    rows <- base[(4L * i - 3L):(4L * i), , drop = FALSE]

    bf <- res$n_bf_dose
    bf_total <- res$total_n_bf
    bf_label <- "Backfilled patients"
    if (percent) {
      bf <- bf / res$total_n_pts * 100
      bf_total <- bf_total / res$total_n_pts * 100
      bf_label <- "Backfilled patients (%)"
    }

    extra <- data.frame(
      Scenario = c("", "", ""),
      Item = c("True response rate (%)", bf_label, "Trial duration (mean)"),
      stringsAsFactors = FALSE
    )
    values <- rbind(res$p_resp * 100, bf, rep(NA_real_, n_doses))
    for (k in seq_len(n_doses)) {
      extra[[dose_names[k]]] <- round(values[, k], digits)
    }
    extra[["Total / No MTD"]] <- c(NA_real_, round(bf_total, digits),
                                   round(res$duration_mean, digits))

    # True DLT rate, true response rate, MTD selected, patients treated,
    # backfilled patients, patients with DLT, trial duration.
    rbind(rows[1L, ], extra[1L, ], rows[2:3, ], extra[2L, ], rows[4L, ],
          extra[3L, ])
  })

  out <- do.call(rbind, blocks)
  rownames(out) <- NULL
  out
}
