# Report whether the optional packages used for cross-implementation checks were
# actually available, so that a silent skip is visible in the test output. The
# report is emitted once per package rather than once per test.
reported_packages <- new.env(parent = emptyenv())

have_package <- function(package) {
  available <- requireNamespace(package, quietly = TRUE)
  if (!exists(package, envir = reported_packages, inherits = FALSE)) {
    assign(package, TRUE, envir = reported_packages)
    message("Package '", package, "' available for cross-checks: ", available)
  }
  available
}

# Call a print method with its console output captured, and return what the
# method returned together with whether it returned visibly. Using this instead
# of expect_invisible() keeps the test output free of printed tables.
quiet_print <- function(x, ...) {
  result <- NULL
  invisible(capture.output(result <- withVisible(print(x, ...))))
  result
}
