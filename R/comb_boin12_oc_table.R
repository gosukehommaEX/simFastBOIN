#' Table of a Comb-BOIN12 Result by Combination
#'
#' @param x A \code{comb_boin12_oc} object.
#' @param digits Integer scalar. Decimal places.
#' @param percent Logical scalar. Patients, DLTs and responses as
#'   percentages of their totals.
#'
#' @return A data frame with one block of rows (one per level of drug A) per
#'   item and one column per level of drug B.
#'
#' @noRd
comb_boin12_oc_table <- function(x, digits, percent = FALSE) {
  n_rows <- x$n_rows
  n_cols <- x$n_cols
  pts <- x$n_pts_dose
  tox <- x$n_tox_dose
  resp <- x$n_resp_dose
  labels <- c("Patients treated", "Patients with DLT",
              "Patients with response")
  if (percent) {
    pts <- pts / x$total_n_pts * 100
    tox <- tox / x$total_n_tox * 100
    resp <- resp / x$total_n_resp * 100
    labels <- paste(labels, "(%)")
  }
  blocks <- list(x$p_true * 100, x$p_resp * 100, x$sel_percent, pts, tox,
                 resp)
  names(blocks) <- c("True DLT rate (%)", "True response rate (%)",
                     "OBD selected (%)", labels)
  rows <- lapply(names(blocks), function(item) {
    values <- matrix(round(as.numeric(blocks[[item]]), digits), n_rows, n_cols)
    colnames(values) <- paste0("B", seq_len(n_cols))
    cbind(
      data.frame(Item = c(item, rep("", n_rows - 1L)),
                 "Drug A" = paste0("A", seq_len(n_rows)),
                 check.names = FALSE, stringsAsFactors = FALSE),
      as.data.frame(values)
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
