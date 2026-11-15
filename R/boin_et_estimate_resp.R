#' Response Estimates for the OBD Selection of BOIN-ET
#'
#' @description
#'   Estimate the response rate of every dose tried in each trial with the
#'   chosen method. The model-based estimates are computed once for each
#'   distinct set of final data, so that trials ending with the same counts
#'   share the computation. Internal.
#'
#' @param n_pts,n_resp
#'   Integer matrices, one row per trial and one column per dose.
#'
#' @param est_method
#'   Character scalar: \code{"observed"}, \code{"multi_iso"} or
#'   \code{"fp_logistic"}.
#'
#' @param need
#'   Logical vector, one entry per trial. The model-based estimates are only
#'   computed for these trials; the others get the observed rates.
#'
#' @return
#'   Numeric matrix of the same shape, \code{NA} for the doses not tried.
#'
#' @keywords internal
#' @noRd
boin_et_estimate_resp <- function(n_pts, n_resp, est_method, need) {
  probe <- n_resp / n_pts
  probe[n_pts == 0L] <- NA_real_
  if (identical(est_method, "observed")) return(probe)

  rows <- which(need & rowSums(n_pts > 0L) >= 2L)
  if (length(rows) == 0L) return(probe)

  key <- apply(cbind(n_pts[rows, , drop = FALSE], n_resp[rows, , drop = FALSE]),
               1, paste, collapse = "_")
  first <- !duplicated(key)
  owner <- match(key, key[first])
  first_rows <- rows[first]

  estimates <- lapply(first_rows, function(r) {
    eva <- which(n_pts[r, ] > 0L)
    if (identical(est_method, "multi_iso")) {
      est <- boin_et_multi_iso(obs = n_resp[r, eva], n = n_pts[r, eva])
    } else {
      est <- boin_et_fp_logit(obs = n_resp[r, eva], n = n_pts[r, eva],
                              dose = eva)
    }
    list(eva = eva, est = est)
  })

  for (i in seq_along(rows)) {
    e <- estimates[[owner[i]]]
    probe[rows[i], e$eva] <- e$est
  }
  probe
}
