#' Plot a TITE-BOIN Decision Table
#'
#' @description
#'   Draw the decision table produced by \code{\link{tite_boin_decision_table}}
#'   as one panel per number of patients treated. Within a panel, each cell is a
#'   pair of pending patients and DLTs, colored by its decision and labeled with
#'   the decision code and, where the decision depends on the follow-up of the
#'   pending patients, the boundaries on the follow-up statistic.
#'
#' @param x
#'   An object of class \code{tite_boin_decision_table}.
#'
#' @param n
#'   Integer vector or \code{NULL}. Numbers of patients treated to draw, one
#'   panel each. Defaults to every number in the table.
#'
#' @param cohort_size
#'   Integer scalar or \code{NULL}. When supplied, only the numbers of patients
#'   that are multiples of \code{cohort_size} are drawn.
#'
#' @param text_size
#'   Numeric scalar. Multiplier applied to every text size in the figure.
#'   Defaults to 1.
#'
#' @param colors
#'   Named character vector of six colors, with names \code{"E"}, \code{"S"},
#'   \code{"D"}, \code{"DE"}, \code{"SUS"} and \code{"Depends"}, the last for the
#'   states whose decision depends on the follow-up statistic. Defaults to the
#'   palette of \code{\link{plot.boin_decision_table}} extended with gray for
#'   suspension and a pale violet for the follow-up dependent states.
#'
#' @param digits
#'   Integer scalar. Number of decimal places of the boundaries in the labels.
#'   Defaults to 2.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   A \pkg{ggplot2} object, which can be printed, saved with
#'   \code{ggplot2::ggsave()} or extended with further layers.
#'
#' @details
#'   A cell whose decision depends on the follow-up statistic shows its codes
#'   over the boundaries, kept short so that they fit in the cell. A label such
#'   as \code{"E/S"} over \code{">=2.15"} means that the dose is escalated when
#'   the statistic is at least 2.15 and kept otherwise, and \code{"S/D"} over
#'   \code{"<=0.52"} that it is de-escalated when the statistic is at most 0.52.
#'   The boundary marked \code{">="} always belongs to the first code, which is
#'   \code{"E"}, or \code{"SUS"} when escalation is blocked, and the one marked
#'   \code{"<="} to \code{"D"}. The statistic is STFT for
#'   \code{method = "imputation"} and ESS for \code{method = "ess"}.
#'
#'   When the table was built with a positive \code{min_follow_up}, the code of
#'   a cell whose escalation also requires the shortest follow-up of the pending
#'   patients to reach that fraction of the window is marked with an asterisk,
#'   as in \code{"E*"} or \code{"E/S*"}, and a caption explains the mark.
#'
#'   Every cell carries its decision code, so the figure does not rely on color
#'   alone. With many panels, or many patients per panel, the cells become small;
#'   draw fewer panels through \code{n}, enlarge the figure, or lower
#'   \code{text_size}. Requires the \pkg{ggplot2} package, which is only
#'   suggested by \pkg{simFastBOIN} rather than required.
#'
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' decisions <- tite_boin_decision_table(target = 0.2, max_n = 15)
#'
#' plot(decisions, n = c(9, 12, 15))
#'
#' # Every cohort end, with smaller text
#' plot(decisions, cohort_size = 3, text_size = 0.8)
#'
#' @seealso \code{\link{tite_boin_decision_table}},
#'   \code{\link{print.tite_boin_decision_table}}
#'
#' @export
plot.tite_boin_decision_table <- function(x, n = NULL, cohort_size = NULL,
                                          text_size = 1, colors = NULL,
                                          digits = 2, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed to plot a decision table", call. = FALSE)
  }
  required <- c("n", "n_tox", "n_pending", "decision", "esc_bound", "deesc_bound")
  if (is.null(attr(x, "statistic")) || !all(required %in% names(x))) {
    stop("'x' must be a decision table returned by tite_boin_decision_table()",
         call. = FALSE)
  }
  if (!is.numeric(text_size) || length(text_size) != 1L ||
      !is.finite(text_size) || text_size <= 0) {
    stop("'text_size' must be a single positive number", call. = FALSE)
  }
  check_count(digits, "digits", 0L)

  available <- sort(unique(x$n))
  if (is.null(n)) {
    n <- available
  } else {
    if (!is.numeric(n) || length(n) < 1L || any(!is.finite(n)) ||
        any(n != round(n))) {
      stop("'n' must contain whole numbers", call. = FALSE)
    }
    absent <- setdiff(n, available)
    if (length(absent) > 0L) {
      stop("'n' contains numbers of patients that are not in the table: ",
           paste(absent, collapse = ", "), call. = FALSE)
    }
    n <- sort(unique(as.integer(n)))
  }
  if (!is.null(cohort_size)) {
    check_count(cohort_size, "cohort_size", 1L)
    n <- n[n %% as.integer(cohort_size) == 0L]
  }
  if (length(n) == 0L) {
    stop("no number of patients left to plot; check 'n' and 'cohort_size'",
         call. = FALSE)
  }

  levels_fill <- c("E", "S", "D", "DE", "SUS", "Depends")
  if (is.null(colors)) {
    colors <- c(E = "#0CA30C", S = "#2A78D6", D = "#EDA100", DE = "#D03B3B",
                SUS = "#A6A6A6", Depends = "#E4DAF2")
  }
  if (!is.character(colors) || !all(levels_fill %in% names(colors))) {
    stop("'colors' must be a character vector named 'E', 'S', 'D', 'DE', ",
         "'SUS' and 'Depends'", call. = FALSE)
  }
  colors <- colors[levels_fill]

  statistic <- attr(x, "statistic")
  keep <- x$n %in% n
  decision <- x$decision[keep]
  esc <- x$esc_bound[keep]
  deesc <- x$deesc_bound[keep]
  fmt <- function(v) formatC(v, format = "f", digits = as.integer(digits))

  # Short labels that fit in a cell: the codes, then the boundary at or above
  # which the first code ("E" or "SUS") applies, then the boundary at or below
  # which "D" applies.
  label <- decision

  # Escalations that also require a minimum follow-up of the pending patients.
  min_follow_up <- attr(x, "min_follow_up")
  if (is.null(min_follow_up)) min_follow_up <- 0
  needs_mf <- min_follow_up > 0 & x$n_pending[keep] > 0L &
    sub("/.*$", "", decision) == "E"
  label[needs_mf] <- paste0(label[needs_mf], "*")
  caption <- if (any(needs_mf)) {
    paste0("* Escalation also requires every pending patient to have been ",
           "followed for at least ", format(100 * min_follow_up),
           "% of the assessment window; accrual is suspended otherwise")
  } else {
    NULL
  }

  has_esc <- !is.na(esc)
  label[has_esc] <- paste0(label[has_esc], "\n>=", fmt(esc[has_esc]))
  has_deesc <- !is.na(deesc)
  label[has_deesc] <- paste0(label[has_deesc], "\n<=", fmt(deesc[has_deesc]))

  panel_levels <- paste(n, "patients treated")
  cells <- data.frame(
    panel = factor(paste(x$n[keep], "patients treated"), levels = panel_levels),
    n_pending = x$n_pending[keep],
    n_tox = x$n_tox[keep],
    fill_group = factor(ifelse(grepl("/", decision, fixed = TRUE), "Depends",
                               decision), levels = levels_fill),
    label = label,
    stringsAsFactors = FALSE
  )

  legend_labels <- c(
    E = "E = escalate to the next higher dose",
    S = "S = stay at the current dose",
    D = "D = de-escalate to the next lower dose",
    DE = "DE = de-escalate and eliminate this dose and all higher doses",
    SUS = "SUS = suspend accrual until more data are available",
    Depends = paste0("Depends on ", statistic, ": the first code applies at or ",
                     "above the value after >=, D at or below the value after <=")
  )

  integer_breaks <- function(limits) seq(ceiling(limits[1L]), floor(limits[2L]))

  ggplot2::ggplot(
    cells,
    ggplot2::aes(x = n_pending, y = n_tox, fill = fill_group)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      colour = "black", size = 2.6 * text_size, lineheight = 0.9
    ) +
    ggplot2::scale_fill_manual(
      values = colors, labels = legend_labels, drop = FALSE, name = NULL
    ) +
    ggplot2::scale_x_continuous(breaks = integer_breaks, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = integer_breaks, expand = c(0, 0)) +
    ggplot2::facet_wrap(~ panel, scales = "free") +
    ggplot2::labs(
      x = "Number of patients with pending DLT data",
      y = "Number of DLTs observed",
      caption = caption
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 9 * text_size),
      axis.title = ggplot2::element_text(size = 11 * text_size, face = "bold"),
      strip.text = ggplot2::element_text(size = 10 * text_size, face = "bold"),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 9 * text_size),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA,
                                           linewidth = 0.6)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1, byrow = TRUE))
}
