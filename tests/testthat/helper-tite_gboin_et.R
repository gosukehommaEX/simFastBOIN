# Scenarios 1, 3, 4 and 7 of Table S1 of Takeda et al. (2023): true
# probabilities of toxicity grades 0 or 1, 2, 3 and 4 and of progressive
# disease, stable disease, partial response and complete response (rows) at
# six doses (columns).
takeda2023_scenario <- function(s) {
  tox <- switch(as.character(s),
    "1" = rbind(
      c(0.94, 0.87, 0.79, 0.68, 0.62, 0.50),
      c(0.05, 0.10, 0.15, 0.20, 0.20, 0.20),
      c(0.01, 0.03, 0.05, 0.10, 0.15, 0.25),
      c(0.00, 0.00, 0.01, 0.02, 0.03, 0.05)
    ),
    "3" = rbind(
      c(0.83, 0.75, 0.62, 0.51, 0.34, 0.19),
      c(0.12, 0.15, 0.18, 0.19, 0.16, 0.11),
      c(0.04, 0.07, 0.11, 0.14, 0.15, 0.11),
      c(0.01, 0.03, 0.09, 0.16, 0.35, 0.59)
    ),
    "4" = rbind(
      c(0.82, 0.65, 0.41, 0.42, 0.34, 0.26),
      c(0.10, 0.20, 0.34, 0.28, 0.31, 0.34),
      c(0.05, 0.10, 0.15, 0.18, 0.21, 0.24),
      c(0.03, 0.05, 0.10, 0.12, 0.14, 0.16)
    ),
    "7" = rbind(
      c(0.50, 0.38, 0.29, 0.19, 0.13, 0.08),
      c(0.20, 0.24, 0.21, 0.16, 0.11, 0.07),
      c(0.25, 0.12, 0.12, 0.10, 0.08, 0.05),
      c(0.05, 0.26, 0.38, 0.55, 0.68, 0.80)
    ),
    stop("unknown scenario")
  )
  eff <- switch(as.character(s),
    "1" = rbind(
      c(0.64, 0.52, 0.45, 0.35, 0.20, 0.05),
      c(0.30, 0.40, 0.40, 0.40, 0.40, 0.15),
      c(0.05, 0.05, 0.10, 0.15, 0.20, 0.35),
      c(0.01, 0.03, 0.05, 0.10, 0.20, 0.45)
    ),
    "3" = rbind(
      c(0.60, 0.50, 0.30, 0.15, 0.10, 0.05),
      c(0.30, 0.20, 0.20, 0.15, 0.15, 0.15),
      c(0.05, 0.20, 0.35, 0.20, 0.20, 0.20),
      c(0.05, 0.10, 0.15, 0.50, 0.55, 0.60)
    ),
    "4" = rbind(
      c(0.30, 0.20, 0.05, 0.05, 0.05, 0.05),
      c(0.35, 0.30, 0.25, 0.20, 0.15, 0.10),
      c(0.30, 0.40, 0.20, 0.25, 0.30, 0.30),
      c(0.05, 0.10, 0.50, 0.50, 0.50, 0.55)
    ),
    "7" = rbind(
      c(0.15, 0.15, 0.10, 0.05, 0.05, 0.05),
      c(0.10, 0.10, 0.10, 0.10, 0.05, 0.02),
      c(0.25, 0.20, 0.20, 0.20, 0.20, 0.20),
      c(0.50, 0.55, 0.60, 0.65, 0.70, 0.73)
    )
  )
  list(p_grade = tox, p_resp_grade = eff)
}

# The design of the simulation study of the article, with the boundaries of
# boin_et_boundary() for the targets 0.47 / 1.5 and 1.75 / 3 unless given.
tite_takeda_args <- function(s = 1, ...) {
  sc <- takeda2023_scenario(s)
  args <- list(
    target = 0.47 / 1.5, target_resp = 1.75 / 3, p_grade = sc$p_grade,
    p_resp_grade = sc$p_resp_grade, score = c(0, 0.5, 1, 1.5),
    score_resp = c(0, 0.25, 1, 3), n_cohort = 12, cohort_size = 3,
    window = 30, resp_window = 60, accrual_rate = 1 / 10
  )
  utils::modifyList(args, list(...))
}

tite_takeda_trials <- function(s = 1, ...) {
  do.call(tite_gboin_et_simulate, tite_takeda_args(s, ...))
}

tite_takeda_oc <- function(s = 1, ...) {
  do.call(sim_tite_gboin_et, tite_takeda_args(s, ...))
}
