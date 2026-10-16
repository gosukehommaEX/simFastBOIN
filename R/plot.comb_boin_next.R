#' Plot the Next Dose Combination of a BOIN Drug-Combination Trial
#'
#' @description
#'   Draw the dose matrix of a drug-combination trial with the data observed at
#'   every combination, the current combination, the combination chosen for the
#'   next cohort and the combinations eliminated, as returned by
#'   \code{\link{comb_boin_next}}.
#'
#' @param x
#'   An object of class \code{comb_boin_next}.
#'
#' @param text_size
#'   Numeric scalar. Multiplier applied to every text size in the figure. Use a
#'   larger value when the figure is exported at high resolution. Defaults to 1.
#'
#' @param colors
#'   Named character vector of five colors, with names \code{"Next"},
#'   \code{"Current"}, \code{"Eliminated"}, \code{"Treated"} and
#'   \code{"Not treated"}. Defaults to the palette of
#'   \code{\link{plot.boin_decision_table}}, with gray and white for the
#'   combinations that are neither current, next nor eliminated.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   A \pkg{ggplot2} object, which can be printed, saved with
#'   \code{ggplot2::ggsave()} or extended with further layers.
#'
#' @details
#'   Each cell shows the number of patients with a DLT and the number treated,
#'   and the status of the current, next and eliminated combinations is also
#'   written in the cell, so the figure does not rely on color alone. The
#'   candidates of an escalation or a de-escalation also show the posterior
#'   probability that their DLT rate lies between \code{lambda_e} and
#'   \code{lambda_d}. The levels of drug A run up the vertical axis and those of
#'   drug B along the horizontal axis.
#'
#'   Requires the \pkg{ggplot2} package, which is only suggested by
#'   \pkg{simFastBOIN} rather than required.
#'
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' n_pts <- matrix(c(6, 6, 0, 0,
#'                   3, 3, 0, 0,
#'                   0, 0, 0, 0), nrow = 3, byrow = TRUE)
#' n_tox <- matrix(c(0, 1, 0, 0,
#'                   1, 2, 0, 0,
#'                   0, 0, 0, 0), nrow = 3, byrow = TRUE)
#'
#' nxt <- comb_boin_next(target = 0.25, n_pts = n_pts, n_tox = n_tox,
#'                       current = c(2, 2))
#' plot(nxt)
#'
#' @seealso \code{\link{comb_boin_next}}
#'
#' @export
plot.comb_boin_next <- function(x, text_size = 1, colors = NULL, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed to plot the decision", call. = FALSE)
  }
  if (!is.numeric(text_size) || length(text_size) != 1L ||
      !is.finite(text_size) || text_size <= 0) {
    stop("'text_size' must be a single positive number", call. = FALSE)
  }

  levels_status <- c("Next", "Current", "Eliminated", "Treated", "Not treated")
  if (is.null(colors)) {
    colors <- c(Next = "#0CA30C", Current = "#2A78D6", Eliminated = "#D03B3B",
                Treated = "#D9D9D9", "Not treated" = "#FFFFFF")
  }
  if (!is.character(colors) || !all(levels_status %in% names(colors))) {
    stop("'colors' must be a character vector named 'Next', 'Current', ",
         "'Eliminated', 'Treated' and 'Not treated'", call. = FALSE)
  }
  colors <- colors[levels_status]

  n_rows <- nrow(x$n_pts)
  n_cols <- ncol(x$n_pts)
  a <- x$current[["dose_a"]]
  b <- x$current[["dose_b"]]

  cells <- data.frame(
    dose_a = rep(seq_len(n_rows), times = n_cols),
    dose_b = rep(seq_len(n_cols), each = n_rows),
    n = as.vector(x$n_pts),
    y = as.vector(x$n_tox),
    eliminated = as.vector(x$eliminated),
    stringsAsFactors = FALSE
  )
  is_next <- rep(FALSE, nrow(cells))
  if (x$decision %in% c("escalate", "de-escalate")) {
    for (r in seq_len(nrow(x$next_dose))) {
      is_next <- is_next | (cells$dose_a == x$next_dose[r, "dose_a"] &
                              cells$dose_b == x$next_dose[r, "dose_b"])
    }
  }
  is_current <- cells$dose_a == a & cells$dose_b == b

  status <- ifelse(cells$n > 0, "Treated", "Not treated")
  status[cells$eliminated] <- "Eliminated"
  status[is_current] <- "Current"
  status[is_next] <- "Next"
  cells$status <- factor(status, levels = levels_status)

  word <- rep("", nrow(cells))
  word[cells$eliminated] <- "eliminated"
  word[is_current] <- paste0(
    "current",
    if (x$decision == "stay") ", stay" else "",
    ifelse(cells$eliminated[is_current], ", eliminated", "")
  )
  word[is_next] <- "next"

  prob <- rep("", nrow(cells))
  if (!is.null(x$candidates)) {
    for (r in seq_len(nrow(x$candidates))) {
      cand <- x$candidates[r, ]
      if (is.na(cand$dose_a) || !cand$available) next
      hit <- cells$dose_a == cand$dose_a & cells$dose_b == cand$dose_b
      prob[hit] <- paste0("Pr = ", format(round(cand$prob_in_interval, 2),
                                          nsmall = 2))
    }
  }

  counts <- ifelse(cells$n > 0, paste0(cells$y, "/", cells$n), "-")
  cells$label <- vapply(seq_len(nrow(cells)), function(i) {
    paste(c(counts[i], word[i][nzchar(word[i])], prob[i][nzchar(prob[i])]),
          collapse = "\n")
  }, character(1))
  cells$text_colour <- ifelse(cells$status %in% c("Next", "Current",
                                                  "Eliminated"),
                              "white", "black")

  name <- function(dose_a, dose_b) paste0("A", dose_a, "B", dose_b)
  title <- switch(
    x$decision,
    escalate = paste("Escalate to", paste(name(x$next_dose[, "dose_a"],
                                                x$next_dose[, "dose_b"]),
                                           collapse = " or ")),
    `de-escalate` = paste("De-escalate to",
                          paste(name(x$next_dose[, "dose_a"],
                                     x$next_dose[, "dose_b"]),
                                collapse = " or ")),
    stay = paste("Stay at", name(a, b)),
    stop = "Stop the trial"
  )

  ggplot2::ggplot(
    cells,
    ggplot2::aes(x = dose_b, y = dose_a, fill = status)
  ) +
    # show.legend = TRUE keeps the legend keys of levels absent from the
    # figure, which ggplot2 3.5.0 and later otherwise leave blank.
    ggplot2::geom_tile(colour = "black", linewidth = 0.4, show.legend = TRUE) +
    ggplot2::geom_text(
      ggplot2::aes(label = label, colour = text_colour),
      size = 3.6 * text_size, lineheight = 0.9
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_manual(values = colors, limits = levels_status,
                               drop = FALSE, name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq_len(n_cols),
                                labels = paste0("B", seq_len(n_cols)),
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = seq_len(n_rows),
                                labels = paste0("A", seq_len(n_rows)),
                                expand = c(0, 0)) +
    ggplot2::labs(
      title = title,
      subtitle = "Cells show DLTs / patients treated",
      x = "Dose level of drug B",
      y = "Dose level of drug A"
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13 * text_size, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10 * text_size),
      axis.text = ggplot2::element_text(size = 11 * text_size),
      axis.title = ggplot2::element_text(size = 12 * text_size, face = "bold"),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 10 * text_size),
      panel.grid = ggplot2::element_blank()
    )
}
