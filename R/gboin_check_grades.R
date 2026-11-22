#' Check the Grade Probabilities and Scores of the gBOIN Design
#'
#' @description
#'   Internal helper. Check the matrix of true grade probabilities and the
#'   scores of the grades, and derive the quantities the simulation needs.
#'
#' @param p_grade Matrix of true probabilities, one row per toxicity
#'   category and one column per dose; every column sums to 1.
#' @param score Score of each category, non-negative, not all zero.
#'
#' @return A list with \code{score_norm} (the scores divided by their
#'   maximum), \code{p_true} (the true mean of the normalized score at each
#'   dose) and \code{cum_prob} (matrix of the shape of \code{p_grade}: the
#'   probability of the category of the row or a higher one).
#'
#' @noRd
gboin_check_grades <- function(p_grade, score) {
  if (!is.matrix(p_grade) || !is.numeric(p_grade) || nrow(p_grade) < 2L ||
      any(!is.finite(p_grade)) || any(p_grade < 0)) {
    stop("'p_grade' must be a numeric matrix of probabilities with one row ",
         "for each toxicity category (at least two) and one column for ",
         "each dose", call. = FALSE)
  }
  if (any(abs(colSums(p_grade) - 1) > 1e-8)) {
    stop("every column of 'p_grade' must sum to 1", call. = FALSE)
  }
  if (!is.numeric(score) || length(score) != nrow(p_grade) ||
      any(!is.finite(score)) || any(score < 0) || max(score) <= 0) {
    stop("'score' must hold one non-negative score for each row of ",
         "'p_grade', not all zero", call. = FALSE)
  }
  score_norm <- as.numeric(score) / max(score)
  p_grade <- unname(p_grade)
  cum_prob <- apply(p_grade, 2L, function(p) rev(cumsum(rev(p))))
  list(
    score_norm = score_norm,
    p_true = as.numeric(colSums(p_grade * score_norm)),
    cum_prob = cum_prob
  )
}
