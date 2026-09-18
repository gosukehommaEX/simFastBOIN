#' Plot a BOIN Decision Table
#'
#' @description
#'   Draw the decision table produced by \code{\link{boin_decision_table}} as a
#'   grid of colored cells, one per pair of DLT and patient counts, with the
#'   decision letter in each cell.
#'
#' @param x
#'   An object of class \code{boin_decision_table}.
#'
#' @param text_size
#'   Numeric scalar. Multiplier applied to every text size in the figure. Use a
#'   larger value when the figure is exported at high resolution, for example
#'   into a report. Defaults to 1.
#'
#' @param colors
#'   Named character vector of four colors, with names \code{"E"}, \code{"S"},
#'   \code{"D"} and \code{"DE"}. Defaults to a palette whose adjacent categories
#'   remain distinguishable under the common forms of color vision deficiency.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   A \pkg{ggplot2} object, which can be printed, saved with
#'   \code{ggplot2::ggsave()} or extended with further layers.
#'
#' @details
#'   Every cell carries its decision letter, so the figure does not rely on
#'   color alone to convey the decision and remains readable in grayscale.
#'
#'   Requires the \pkg{ggplot2} package, which is only suggested by
#'   \pkg{simFastBOIN} rather than required.
#'
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' decisions <- boin_decision_table(target = 0.30, max_n = 12)
#'
#' plot(decisions)
#'
#' # Larger text for a high resolution export
#' plot(decisions, text_size = 2)
#'
#' # A palette of your own
#' plot(decisions, colors = c(E = "#4DAF4A", S = "#377EB8",
#'                            D = "#FF7F00", DE = "#E41A1C"))
#'
#' @seealso \code{\link{boin_decision_table}}, \code{\link{print.boin_decision_table}}
#'
#' @export
plot.boin_decision_table <- function(x, text_size = 1, colors = NULL, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed to plot a decision table", call. = FALSE)
  }
  if (!is.numeric(text_size) || length(text_size) != 1L ||
      !is.finite(text_size) || text_size <= 0) {
    stop("'text_size' must be a single positive number", call. = FALSE)
  }

  levels_decision <- c("E", "S", "D", "DE")

  if (is.null(colors)) {
    colors <- c(E = "#0CA30C", S = "#2A78D6", D = "#EDA100", DE = "#D03B3B")
  }
  if (!is.character(colors) || !all(levels_decision %in% names(colors))) {
    stop("'colors' must be a character vector named 'E', 'S', 'D' and 'DE'",
         call. = FALSE)
  }
  colors <- colors[levels_decision]

  n_tox_levels <- as.integer(rownames(x))
  n_pts_levels <- as.integer(colnames(x))

  cells <- data.frame(
    n_pts = rep(n_pts_levels, each = length(n_tox_levels)),
    n_tox = rep(n_tox_levels, times = length(n_pts_levels)),
    decision = as.vector(unclass(x)),
    stringsAsFactors = FALSE
  )
  cells <- cells[!is.na(cells$decision), , drop = FALSE]
  cells$decision <- factor(cells$decision, levels = levels_decision)

  legend_labels <- c(
    E = "E = escalate to the next higher dose",
    S = "S = stay at the current dose",
    D = "D = de-escalate to the next lower dose",
    DE = "DE = de-escalate and eliminate this dose and all higher doses"
  )

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
    ggplot2::scale_x_continuous(breaks = n_pts_levels, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = n_tox_levels, expand = c(0, 0)) +
    ggplot2::labs(
      x = "Number of evaluable patients treated at the current dose",
      y = "Number of patients with a DLT"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 11 * text_size),
      axis.title = ggplot2::element_text(size = 12 * text_size, face = "bold"),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 10 * text_size),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA,
                                           linewidth = 0.8)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1, byrow = TRUE))
}
