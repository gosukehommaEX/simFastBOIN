#' Assemble the Table of a Drug-Combination Simulation
#'
#' @description
#'   Stack the true DLT rates, the selection percentages and the average numbers
#'   of patients and DLTs of a \code{comb_boin_oc} object into one data frame,
#'   with one row per level of drug A in each block. Internal.
#'
#' @param x
#'   An object of class \code{comb_boin_oc}.
#'
#' @param digits
#'   Integer scalar. Number of decimal places.
#'
#' @return
#'   A data frame.
#'
#' @keywords internal
#' @noRd
comb_oc_table <- function(x, digits) {
  n_rows <- x$n_rows
  n_cols <- x$n_cols
  blocks <- list(
    "True DLT rate (%)" = x$p_true * 100,
    "MTD selected (%)" = x$sel_percent,
    "Patients treated" = x$n_pts_dose,
    "Patients with DLT" = x$n_tox_dose
  )
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
