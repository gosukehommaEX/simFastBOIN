# Scenarios 1, 5 and 6 of Table S1 of Takeda et al. (2022): true
# probabilities of toxicity grades 0 or 1, 2, 3 and 4 and of stable or
# progressive disease, partial response and complete response (rows) at six
# doses (columns).
takeda_scenario <- function(s) {
  tox <- switch(as.character(s),
    "1" = rbind(
      c(0.83, 0.75, 0.62, 0.51, 0.34, 0.19),
      c(0.12, 0.15, 0.18, 0.19, 0.16, 0.11),
      c(0.04, 0.07, 0.11, 0.14, 0.15, 0.11),
      c(0.01, 0.03, 0.09, 0.16, 0.35, 0.59)
    ),
    "5" = rbind(
      c(1.00, 0.91, 0.88, 0.86, 0.80, 0.65),
      c(0.00, 0.06, 0.07, 0.08, 0.10, 0.13),
      c(0.00, 0.03, 0.04, 0.05, 0.08, 0.14),
      c(0.00, 0.00, 0.01, 0.01, 0.02, 0.08)
    ),
    "6" = rbind(
      c(0.50, 0.38, 0.29, 0.19, 0.13, 0.08),
      c(0.25, 0.24, 0.21, 0.16, 0.11, 0.07),
      c(0.11, 0.12, 0.12, 0.10, 0.08, 0.05),
      c(0.14, 0.26, 0.38, 0.55, 0.68, 0.80)
    ),
    stop("unknown scenario")
  )
  cr <- switch(as.character(s),
    "1" = c(0.05, 0.10, 0.15, 0.50, 0.55, 0.60),
    "5" = c(0.10, 0.20, 0.30, 0.55, 0.40, 0.30),
    "6" = c(0.50, 0.55, 0.60, 0.65, 0.70, 0.75)
  )
  pr <- switch(as.character(s),
    "1" = c(0.05, 0.20, 0.35, 0.20, 0.20, 0.20),
    "5" = c(0.30, 0.30, 0.30, 0.20, 0.20, 0.20),
    "6" = c(0.20, 0.20, 0.20, 0.20, 0.20, 0.20)
  )
  list(p_grade = tox, p_resp_grade = rbind(1 - cr - pr, pr, cr))
}

ees_score <- function() c(0, 1, 1.5)
ees_target <- function() 0.95 / 1.5

takeda_trials <- function(s = 1, ...) {
  sc <- takeda_scenario(s)
  gboin_et_simulate(
    target = 0.47 / 1.5, target_resp = ees_target(), p_grade = sc$p_grade,
    p_resp_grade = sc$p_resp_grade, score = c(0, 0.5, 1, 1.5),
    score_resp = ees_score(), n_cohort = 10, cohort_size = 3, ...
  )
}

takeda_oc <- function(s = 1, ...) {
  sc <- takeda_scenario(s)
  sim_gboin_et(
    target = 0.47 / 1.5, target_resp = ees_target(), p_grade = sc$p_grade,
    p_resp_grade = sc$p_resp_grade, score = c(0, 0.5, 1, 1.5),
    score_resp = ees_score(), n_cohort = 10, cohort_size = 3, ...
  )
}
