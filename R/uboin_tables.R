#' Elimination and Admissibility Tables of U-BOIN
#'
#' @description
#'   Tabulate, for every number of patients at a dose and every number of
#'   DLTs or responses, whether the dose is eliminated in stage I, toxic or
#'   futile, and derive the boundaries of stage I. Internal.
#'
#' @param des
#'   List returned by \code{uboin_check_design()}.
#'
#' @param max_n
#'   Integer scalar. Largest number of patients at one dose.
#'
#' @return
#'   A list with the integer vectors \code{elim}, \code{toxic} and
#'   \code{futile}, laid out as \code{n * (max_n + 1) + count} (n and count
#'   from 0 to \code{max_n}), and \code{b_esc}, \code{b_deesc} and
#'   \code{b_elim}, indexed by the number of patients (1 to \code{max_n});
#'   \code{b_esc} and \code{b_deesc} only when \code{des} holds the
#'   boundaries of stage I.
#'
#' @details
#'   A dose with \code{y} DLTs among \code{n} patients is toxic when
#'   \code{Pr(p > tox_upper) > cutoff_eli} under the Beta posterior with the
#'   prior \code{prior_tox}; in stage I it is then eliminated if \code{n >= 3}.
#'   It is futile when \code{Pr(q < resp_lower) > cutoff_eli_resp} under the
#'   Beta posterior of its responses with the prior \code{prior_resp}. The
#'   escalation boundary is \code{floor(lambda_e * n)} and the de-escalation
#'   boundary that of \code{\link{boin_boundary}}, capped at the elimination
#'   boundary.
#'
#' @keywords internal
#' @noRd
uboin_tables <- function(des, max_n) {
  size <- max_n + 1L
  toxic <- matrix(0L, size, size)
  futile <- matrix(0L, size, size)
  for (n in 0:max_n) {
    count <- 0:n
    above <- 1 - stats::pbeta(des$tox_upper, des$prior_tox[1L] + count,
                              des$prior_tox[2L] + n - count)
    below <- stats::pbeta(des$resp_lower, des$prior_resp[1L] + count,
                          des$prior_resp[2L] + n - count)
    toxic[n + 1L, count + 1L] <- as.integer(above > des$cutoff_eli)
    futile[n + 1L, count + 1L] <- as.integer(below > des$cutoff_eli_resp)
  }
  elim <- toxic
  elim[seq_len(min(3L, size)), ] <- 0L

  n <- seq_len(max_n)
  b_elim <- vapply(n, function(k) {
    hit <- which(elim[k + 1L, seq_len(k + 1L)] == 1L)
    if (length(hit) > 0L) as.integer(hit[1L] - 1L) else NA_integer_
  }, integer(1))
  out <- list(elim = as.integer(t(elim)), toxic = as.integer(t(toxic)),
              futile = as.integer(t(futile)), b_elim = b_elim)
  if (is.na(des$lambda_e)) return(out)

  b_esc <- as.integer(floor(des$lambda_e * n))
  scaled <- des$lambda_d * n
  is_whole <- abs(round(scaled) - scaled) < 1e-12
  b_deesc <- as.integer(ifelse(is_whole, round(scaled) + 1, ceiling(scaled)))
  capped <- !is.na(b_elim) & b_deesc > b_elim
  b_deesc[capped] <- b_elim[capped]

  c(out, list(b_esc = b_esc, b_deesc = b_deesc))
}
