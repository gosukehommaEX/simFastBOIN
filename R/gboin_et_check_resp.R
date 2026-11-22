#' Check the Efficacy Categories of gBOIN-ET
#'
#' @description
#'   Validate the true probabilities of the efficacy categories and their
#'   scores, normalize the scores by their maximum, and compute the true mean
#'   normalized score of every dose and the probabilities of a category or a
#'   higher one. The toxicity categories are checked by
#'   \code{gboin_check_grades()}. Internal.
#'
#' @param p_resp_grade
#'   Numeric matrix. True probability of each efficacy category (rows) at
#'   each dose (columns).
#'
#' @param score_resp
#'   Numeric vector. Score of each efficacy category.
#'
#' @param n_doses
#'   Integer scalar. Number of doses of the toxicity scenario.
#'
#' @return
#'   A list with the components \code{score_norm}, \code{p_resp} and
#'   \code{cum_prob}.
#'
#' @keywords internal
#' @noRd
gboin_et_check_resp <- function(p_resp_grade, score_resp, n_doses) {
  if (!is.matrix(p_resp_grade) || !is.numeric(p_resp_grade) ||
      nrow(p_resp_grade) < 2L || any(!is.finite(p_resp_grade)) ||
      any(p_resp_grade < 0)) {
    stop("'p_resp_grade' must be a numeric matrix of probabilities with one ",
         "row for each efficacy category (at least two) and one column for ",
         "each dose", call. = FALSE)
  }
  if (ncol(p_resp_grade) != n_doses) {
    stop("'p_resp_grade' must have one column per dose in 'p_grade'",
         call. = FALSE)
  }
  if (any(abs(colSums(p_resp_grade) - 1) > 1e-8)) {
    stop("every column of 'p_resp_grade' must sum to 1", call. = FALSE)
  }
  if (!is.numeric(score_resp) || length(score_resp) != nrow(p_resp_grade) ||
      any(!is.finite(score_resp)) || any(score_resp < 0) ||
      max(score_resp) <= 0) {
    stop("'score_resp' must hold one non-negative score for each row of ",
         "'p_resp_grade', not all zero", call. = FALSE)
  }
  score_norm <- as.numeric(score_resp) / max(score_resp)
  p_resp_grade <- unname(p_resp_grade)
  cum_prob <- apply(p_resp_grade, 2L, function(p) rev(cumsum(rev(p))))
  list(
    score_norm = score_norm,
    p_resp = as.numeric(colSums(p_resp_grade * score_norm)),
    cum_prob = cum_prob
  )
}
