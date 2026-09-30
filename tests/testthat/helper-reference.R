# Whether an optional package used for cross-implementation checks can be
# loaded. It is called after skip_if_not_installed(), so a missing package is
# reported by testthat as a skip, with its reason, in the S column and in the
# summary of the test results; no separate message is needed.
have_package <- function(package) {
  requireNamespace(package, quietly = TRUE)
}

# Call a print method with its console output captured, and return what the
# method returned together with whether it returned visibly. Using this instead
# of expect_invisible() keeps the test output free of printed tables.
quiet_print <- function(x, ...) {
  result <- NULL
  invisible(capture.output(result <- withVisible(print(x, ...))))
  result
}
