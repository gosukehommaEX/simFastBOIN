# Run the R code of the cheat sheet and draw its figures
# -----------------------------------------------------------------------------
#
# The cheat sheet, pkgdown/assets/cheatsheet/simFastBOIN_cheatsheet.html, is
# a single HTML page that pkgdown copies to the site. This script evaluates
# every block of R code of the sheet, in order and in one environment, so
# that the code shown on the sheet is known to run, and draws the figures of
# the sheet from the objects that the code creates:
#
#   figures/decision.png    plot(decisions)
#   figures/comb_next.png   plot(nxt)
#   figures/results.png     plot(results)
#   figures/compare.png     plot(comparison)
#
# Usage, from the package root, with the package and 'ggplot2' installed:
#
#   source("pkgdown/cheatsheet_figures.R")
#   res <- run_cheatsheet()
#
# The printed output of the code is discarded; the warnings are collected
# and returned with the files written. Run it before every update of the
# sheet.
# -----------------------------------------------------------------------------

# The blocks of R code of the sheet, as text.
cheatsheet_code <- function(sheet) {
  page <- paste(readLines(sheet, warn = FALSE, encoding = "UTF-8"),
                collapse = "\n")
  m <- gregexpr('<pre><code class="language-r">[\\s\\S]*?</code></pre>',
                page, perl = TRUE)
  blocks <- regmatches(page, m)[[1]]
  blocks <- sub('^<pre><code class="language-r">', "", blocks)
  blocks <- sub("</code></pre>$", "", blocks)
  unescape <- function(x) {
    x <- gsub("&lt;", "<", x, fixed = TRUE)
    x <- gsub("&gt;", ">", x, fixed = TRUE)
    x <- gsub("&quot;", "\"", x, fixed = TRUE)
    gsub("&amp;", "&", x, fixed = TRUE)
  }
  vapply(blocks, unescape, character(1), USE.NAMES = FALSE)
}

# Lines of code that install packages or open a vignette are not run.
cheatsheet_skip <- function(code) {
  lines <- strsplit(code, "\n", fixed = TRUE)[[1]]
  skip <- grepl("install.packages(", lines, fixed = TRUE) |
    grepl("pak::pak(", lines, fixed = TRUE) |
    grepl("vignette(", lines, fixed = TRUE)
  paste(lines[!skip], collapse = "\n")
}

run_cheatsheet <- function(
    sheet = file.path("pkgdown", "assets", "cheatsheet",
                      "simFastBOIN_cheatsheet.html"),
    fig_dir = file.path("pkgdown", "assets", "cheatsheet", "figures")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("install the 'ggplot2' package first", call. = FALSE)
  }
  blocks <- cheatsheet_code(sheet)
  env <- new.env(parent = globalenv())
  warnings <- character(0)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (i in seq_along(blocks)) {
    message(sprintf("block %d of %d", i, length(blocks)))
    exprs <- parse(text = cheatsheet_skip(blocks[i]))
    for (e in exprs) {
      withCallingHandlers(
        utils::capture.output({
          r <- withVisible(eval(e, envir = env))
          if (r$visible) print(r$value)
        }),
        warning = function(w) {
          warnings <<- c(warnings, sprintf("block %d: %s", i,
                                           conditionMessage(w)))
          invokeRestart("muffleWarning")
        }
      )
    }
  }
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  figures <- list(
    "decision.png" = list(plot = plot(env$decisions), height = 3.6),
    "comb_next.png" = list(plot = plot(env$nxt), height = 3.6),
    "results.png" = list(plot = plot(env$results), height = 3.6),
    "compare.png" = list(plot = plot(env$comparison), height = 4.8)
  )
  files <- character(0)
  for (f in names(figures)) {
    path <- file.path(fig_dir, f)
    ggplot2::ggsave(path, figures[[f]]$plot, width = 7,
                    height = figures[[f]]$height, dpi = 150, bg = "white")
    files <- c(files, path)
  }
  list(blocks = length(blocks), files = files, warnings = warnings)
}
