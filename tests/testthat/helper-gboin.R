# Scenarios of Table 4 of Mu et al. (2019): true probabilities of grades 0 or
# 1, 2, 3 and 4 (rows) at six doses (columns).
mu_scenario <- function(s) {
  switch(as.character(s),
    "1" = rbind(
      c(0.83, 0.75, 0.62, 0.51, 0.34, 0.19),
      c(0.12, 0.15, 0.18, 0.19, 0.16, 0.11),
      c(0.04, 0.07, 0.11, 0.14, 0.15, 0.11),
      c(0.01, 0.03, 0.09, 0.16, 0.35, 0.59)
    ),
    "6" = rbind(
      c(0.50, 0.38, 0.29, 0.19, 0.13, 0.08),
      c(0.25, 0.24, 0.21, 0.16, 0.11, 0.07),
      c(0.11, 0.12, 0.12, 0.10, 0.08, 0.05),
      c(0.14, 0.26, 0.38, 0.55, 0.68, 0.80)
    ),
    stop("unknown scenario")
  )
}

ets_score <- function() c(0, 0.5, 1, 1.5)
ets_target <- function() 0.47 / 1.5
