# Internal argument validation helpers. None of these are exported.

check_scalar_prob <- function(x, name, lower = 0, upper = 1) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop("'", name, "' must be a single finite number", call. = FALSE)
  }
  if (x <= lower || x >= upper) {
    stop("'", name, "' must lie strictly between ", lower, " and ", upper,
         call. = FALSE)
  }
  invisible(TRUE)
}

check_thresholds <- function(target, p_saf, p_tox) {
  check_scalar_prob(target, "target")
  check_scalar_prob(p_saf, "p_saf")
  check_scalar_prob(p_tox, "p_tox")
  if (target < 0.05) stop("'target' is too low; use a value of at least 0.05", call. = FALSE)
  if (target > 0.6) stop("'target' is too high; use a value of at most 0.6", call. = FALSE)
  if ((target - p_saf) < 0.1 * target) {
    stop("'p_saf' must be clearly below 'target'; ",
         "the difference has to exceed 0.1 * target", call. = FALSE)
  }
  if ((p_tox - target) < 0.1 * target) {
    stop("'p_tox' must be clearly above 'target'; ",
         "the difference has to exceed 0.1 * target", call. = FALSE)
  }
  invisible(TRUE)
}

check_p_true <- function(p_true) {
  if (!is.numeric(p_true) || length(p_true) < 1L || any(!is.finite(p_true))) {
    stop("'p_true' must be a numeric vector of finite DLT probabilities", call. = FALSE)
  }
  if (any(p_true < 0) || any(p_true > 1)) {
    stop("'p_true' must contain values between 0 and 1", call. = FALSE)
  }
  if (is.unsorted(p_true)) {
    warning("'p_true' is not non-decreasing across dose levels; ",
            "the BOIN design assumes that toxicity increases with dose",
            call. = FALSE)
  }
  invisible(TRUE)
}

check_count <- function(x, name, min_value = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != as.integer(x)) {
    stop("'", name, "' must be a single whole number", call. = FALSE)
  }
  if (x < min_value) {
    stop("'", name, "' must be at least ", min_value, call. = FALSE)
  }
  invisible(TRUE)
}

# Expand cohort_size to one entry per cohort. A scalar is recycled, a shorter
# vector is padded with its last element, and a longer one is truncated.
expand_cohort_size <- function(cohort_size, n_cohort) {
  if (!is.numeric(cohort_size) || length(cohort_size) < 1L ||
      any(!is.finite(cohort_size)) || any(cohort_size != as.integer(cohort_size)) ||
      any(cohort_size < 1)) {
    stop("'cohort_size' must contain whole numbers of at least 1", call. = FALSE)
  }
  cohort_size <- as.integer(cohort_size)
  n_given <- length(cohort_size)
  if (n_given == n_cohort) {
    out <- cohort_size
  } else if (n_given < n_cohort) {
    out <- c(cohort_size, rep(cohort_size[n_given], n_cohort - n_given))
  } else {
    out <- cohort_size[seq_len(n_cohort)]
  }
  out
}

# Capture and restore the state of the random number generator, so that setting
# a seed for a simulation does not change the state of the user's session.
get_random_seed <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}

restore_random_seed <- function(old_seed) {
  if (!is.null(old_seed)) {
    assign(".Random.seed", old_seed, envir = globalenv())
  }
  invisible(NULL)
}

# Coerce a vector or matrix of counts to an integer matrix with one row per trial.
as_count_matrix <- function(x, name) {
  if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
  if (!is.numeric(x)) stop("'", name, "' must be numeric", call. = FALSE)
  if (any(is.na(x))) stop("'", name, "' must not contain missing values", call. = FALSE)
  if (any(x < 0) || any(x != round(x))) {
    stop("'", name, "' must contain non-negative whole numbers", call. = FALSE)
  }
  storage.mode(x) <- "integer"
  x
}

# Accept either a list of list(name, p_true) or a named list of numeric vectors.
normalise_scenarios <- function(scenarios) {
  if (!is.list(scenarios) || length(scenarios) == 0L) {
    stop("'scenarios' must be a non-empty list", call. = FALSE)
  }
  nms <- names(scenarios)
  out <- vector("list", length(scenarios))
  for (i in seq_along(scenarios)) {
    element <- scenarios[[i]]
    if (is.numeric(element)) {
      name <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else paste("Scenario", i)
      out[[i]] <- list(name = name, p_true = element)
    } else if (is.list(element) && !is.null(element$p_true)) {
      name <- if (!is.null(element$name)) {
        as.character(element$name)
      } else if (!is.null(nms) && nzchar(nms[i])) {
        nms[i]
      } else {
        paste("Scenario", i)
      }
      out[[i]] <- list(name = name, p_true = element$p_true)
    } else {
      stop("each scenario must be a numeric vector or a list with a 'p_true' element",
           call. = FALSE)
    }
    check_p_true(out[[i]]$p_true)
  }
  n_doses <- length(out[[1L]]$p_true)
  if (any(vapply(out, function(z) length(z$p_true), integer(1)) != n_doses)) {
    stop("all scenarios must have the same number of doses", call. = FALSE)
  }
  if (anyDuplicated(vapply(out, function(z) z$name, character(1)))) {
    stop("scenario names must be unique", call. = FALSE)
  }
  out
}
