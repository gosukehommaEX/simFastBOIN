#' Plot a 3+3 Decision Table
#'
#' @description
#'   Draw the decision table produced by \code{\link{decision_table_3p3}} as a
#'   grid of colored cells, one per pair of patient and DLT counts, with the
#'   decision code in each cell. With \code{mtd_rule = "expand"} the search for
#'   the MTD is drawn in a second panel.
#'
#' @param x
#'   An object of class \code{decision_table_3p3}.
#'
#' @param text_size
#'   Numeric scalar. Multiplier applied to every text size in the figure.
#'   Defaults to 1.
#'
#' @param colors
#'   Named character vector of colors with the names \code{"E"}, \code{"S"} and
#'   \code{"STOP"}, and also \code{"MTD"} and \code{"D"} for
#'   \code{mtd_rule = "expand"}. Defaults to the palette of
#'   \code{\link{plot.boin_decision_table}}, with its color for elimination used
#'   for \code{"STOP"} and a violet for \code{"MTD"}.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   A \pkg{ggplot2} object, which can be printed, saved with
#'   \code{ggplot2::ggsave()} or extended with further layers.
#'
#' @details
#'   Every cell carries its decision code, so the figure does not rely on color
#'   alone. Requires the \pkg{ggplot2} package, which is only suggested by
#'   \pkg{simFastBOIN} rather than required.
#'
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' plot(decision_table_3p3())
#'
#' plot(decision_table_3p3(mtd_rule = "expand"))
#'
#' @seealso \code{\link{decision_table_3p3}},
#'   \code{\link{print.decision_table_3p3}}
#'
#' @export
plot.decision_table_3p3 <- function(x, text_size = 1, colors = NULL, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed to plot a decision table", call. = FALSE)
  }
  required <- c("stage", "n", "n_tox", "decision")
  if (is.null(attr(x, "mtd_rule")) || !all(required %in% names(x))) {
    stop("'x' must be a decision table returned by decision_table_3p3()",
         call. = FALSE)
  }
  if (!is.numeric(text_size) || length(text_size) != 1L ||
      !is.finite(text_size) || text_size <= 0) {
    stop("'text_size' must be a single positive number", call. = FALSE)
  }

  levels_all <- c("E", "S", "STOP", "MTD", "D")
  levels_used <- levels_all[levels_all %in% x$decision]
  if (is.null(colors)) {
    colors <- c(E = "#0CA30C", S = "#2A78D6", STOP = "#D03B3B",
                MTD = "#A48AD3", D = "#EDA100")
  }
  if (!is.character(colors) || !all(levels_used %in% names(colors))) {
    stop("'colors' must be a character vector named ",
         paste0("'", levels_used, "'", collapse = ", "), call. = FALSE)
  }
  colors <- colors[levels_used]

  stage_labels <- c(escalation = "Dose escalation", search = "Search for the MTD")
  stages <- intersect(names(stage_labels), unique(x$stage))
  cells <- data.frame(
    panel = factor(unname(stage_labels[x$stage]),
                   levels = unname(stage_labels[stages])),
    n_pts = factor(x$n, levels = sort(unique(x$n))),
    n_tox = x$n_tox,
    decision = factor(x$decision, levels = levels_used),
    stringsAsFactors = FALSE
  )

  legend_labels <- c(
    E = "E = escalate to the next higher dose",
    S = "S = treat three more patients at the current dose",
    STOP = "STOP = stop escalation; this dose and all higher doses are too toxic",
    MTD = "MTD = select the current dose as the MTD",
    D = "D = move the search to the next lower dose"
  )[levels_used]

  ggplot2::ggplot(
    cells,
    ggplot2::aes(x = n_pts, y = n_tox, fill = decision)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = decision),
      colour = "black", fontface = "bold", size = 4 * text_size
    ) +
    ggplot2::scale_fill_manual(
      values = colors, labels = legend_labels, drop = FALSE, name = NULL
    ) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = seq(0L, max(x$n_tox)),
                                expand = c(0, 0)) +
    ggplot2::facet_wrap(~ panel) +
    ggplot2::labs(
      x = "Number of evaluable patients treated at the current dose",
      y = "Number of patients with a DLT"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 11 * text_size),
      axis.title = ggplot2::element_text(size = 12 * text_size, face = "bold"),
      strip.text = ggplot2::element_text(size = 11 * text_size, face = "bold"),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 10 * text_size),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA,
                                           linewidth = 0.8)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1, byrow = TRUE))
}
