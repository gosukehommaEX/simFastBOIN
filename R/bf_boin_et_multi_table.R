#' Table of the BF-BOIN-ET Operating Characteristics Across Scenarios
#'
#' @description
#'   Stack the tables of \code{bf_boin_et_table()} of several scenarios, with
#'   the scenario name on the first row of each block.
#'
#' @param results List of \code{bf_boin_et_oc} objects.
#' @param scenario_names Character vector of scenario names.
#' @param percent Logical scalar, passed to \code{bf_boin_et_table()}.
#' @param digits Integer scalar. Number of decimal places.
#'
#' @return A data frame with columns \code{Scenario}, \code{Item}, one column
#'   per dose and \code{"Total / No OBD"}.
#'
#' @noRd
bf_boin_et_multi_table <- function(results, scenario_names, percent = FALSE,
                                   digits = 1) {
  blocks <- lapply(seq_along(results), function(i) {
    tab <- bf_boin_et_table(results[[i]], percent)
    num <- vapply(tab, is.numeric, logical(1))
    tab[num] <- lapply(tab[num], round, digits = digits)
    cbind(Scenario = c(scenario_names[i], rep("", nrow(tab) - 1L)), tab,
          stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, blocks)
  rownames(out) <- NULL
  out
}
