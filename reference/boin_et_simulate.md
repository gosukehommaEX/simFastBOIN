# Simulate BOIN-ET Trials

Run the BOIN-ET dose-finding algorithm of Takeda, Taguri and Morita
(2018) over many simulated trials and return the patient, DLT and
response counts at every dose. This is the simulation engine behind
[`sim_boin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_et.md),
exposed for users who want the raw trial data.

## Usage

``` r
boin_et_simulate(
  target,
  target_resp,
  p_true,
  p_resp,
  n_cohort,
  cohort_size,
  n_trials = 10000,
  start_dose = 1,
  n_earlystop = NULL,
  p_saf = NULL,
  p_tox = NULL,
  p_saf_resp = NULL,
  cutoff_eli = 0.95,
  cutoff_eli_resp = 0.99,
  p_eli = NULL,
  p_eli_resp = NULL,
  resp_cor = 0.2,
  resp_cor_type = c("copula", "binary"),
  boinet_compatible = FALSE,
  prior = rep(1/6, 6),
  seed = 123
)
```

## Arguments

- target:

  Numeric scalar. Target DLT probability, for example 0.30.

- target_resp:

  Numeric scalar. Target response probability, for example 0.60.

- p_true:

  Numeric vector. True DLT probability at each dose level, in increasing
  dose order.

- p_resp:

  Numeric vector of the same length as `p_true`. True response
  probability at each dose level.

- n_cohort:

  Integer scalar. Number of cohorts in a trial.

- cohort_size:

  Integer scalar or vector. Number of patients per cohort. A scalar is
  used for every cohort. A vector shorter than `n_cohort` is padded with
  its last element and a longer one is truncated.

- n_trials:

  Integer scalar. Number of trials to simulate. Defaults to 10000.

- start_dose:

  Integer scalar. Dose level for the first cohort. Defaults to 1.

- n_earlystop:

  Integer scalar or `NULL`. The trial stops once this many patients have
  been treated at any dose. Defaults to `NULL`, the maximum sample size,
  which is the default of the 'boinet' package.

- p_saf:

  Numeric scalar. Highest DLT probability deemed subtherapeutic.
  Defaults to `0.1 * target`, the value recommended by Takeda et al.
  (2018), which is far below the `0.6 * target` of the BOIN design.

- p_tox:

  Numeric scalar. Lowest DLT probability deemed overly toxic. Defaults
  to `1.4 * target`.

- p_saf_resp:

  Numeric scalar. Highest response probability deemed subtherapeutic.
  Defaults to `0.6 * target_resp`.

- cutoff_eli:

  Numeric scalar. A dose is no longer admissible when the posterior
  probability that its DLT rate is below `p_eli` falls below
  `1 - cutoff_eli`. Defaults to 0.95.

- cutoff_eli_resp:

  Numeric scalar. A dose is no longer admissible when the posterior
  probability that its response rate is above `p_eli_resp` falls below
  `1 - cutoff_eli_resp`. Defaults to 0.99.

- p_eli:

  Numeric scalar or `NULL`. DLT probability of the rule above. Defaults
  to `NULL`, which uses `target` as the 'boinet' package does. Takeda et
  al. (2018) used 0.40 with a cutoff of 0.90.

- p_eli_resp:

  Numeric scalar or `NULL`. Response probability of the rule above.
  Defaults to `NULL`, which uses `p_saf_resp` as the 'boinet' package
  does. Takeda et al. (2018) used 0.20 with a cutoff of 0.90.

- resp_cor:

  Numeric scalar strictly between -1 and 1. Correlation between the DLT
  and the response of a patient. Defaults to 0.2, the default of the
  'boinet' package.

- resp_cor_type:

  Character scalar. How `resp_cor` links the two outcomes: `"copula"`
  (the default) uses it as the correlation of a Gaussian copula, which
  is how the 'boinet' package generates its data; `"binary"` uses it as
  the correlation of the two binary outcomes, as in the simulations of
  Takeda et al. (2018).

- boinet_compatible:

  Logical scalar. Apply the rules of the 'boinet' package (version
  1.6.0) exactly. Defaults to `FALSE`, which keeps to Takeda et al.
  (2018), who do not allow dose skipping during escalation; see Details.

- prior:

  Numeric vector of length 6. Prior probabilities of the six hypotheses.
  Defaults to `rep(1 / 6, 6)`, the noninformative choice of the article.
  The third hypothesis, under which every decision is correct, does not
  enter the criterion.

- seed:

  Integer scalar or `NULL`. Random seed. When `NULL` the current state
  of the random number generator is used. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `boin_et_trials`, which is a list with components

- n_pts:

  Integer matrix, `n_trials` by number of doses, of patients treated.

- n_tox:

  Integer matrix of the same shape, of DLTs observed.

- n_resp:

  Integer matrix of the same shape, of responses observed.

- cohorts_used:

  Integer vector of cohorts enrolled in each trial.

- stop_reason:

  Character vector giving why each trial ended.

- early_stop:

  Logical vector, `TRUE` for the trials stopped early, which select no
  OBD.

- design:

  `"BOIN-ET"`.

- boundary:

  The
  [`boin_et_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_et_boundary.md)
  object used.

- settings:

  List of the design parameters.

## Details

Every cohort is treated at the current dose and its outcomes are
observed before the next decision, which follows
[`boin_et_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_et_decision_table.md).
After each cohort the admissibility of every dose is re-derived from its
data, and a dose excluded by its DLT data also excludes every higher
dose. The trial stops early when no dose is admissible, or when the dose
should be lowered and no admissible dose is found at or below the
intended one. When the intended dose is not admissible, an escalation
moves to the lowest admissible dose above it, which can only skip a dose
excluded by its response data (the next cohort stays when there is
none), and a stay or a de-escalation moves to the highest admissible
dose below it. An escalation from the highest dose stays there, or ends
the trial when that dose is not admissible. The trial also ends once a
dose has `n_earlystop` patients. These two ends are not early stops: the
OBD is still selected.

These are the rules of
[`boinet::boinet()`](https://rdrr.io/pkg/boinet/man/boinet.html) with
two changes, which follow Takeda et al. (2018), who do not allow dose
skipping during escalation, and the later rules of the same group
(Takeda et al., 2025, and the 'bfboinet' package). In 'boinet', a dose
excluded by its DLT data does not exclude the doses above it, so that an
escalation can skip it and reach an untried dose, and an escalation from
the highest dose stays there even when it is excluded.
`boinet_compatible = TRUE` applies these rules of 'boinet'.

Each patient draws two uniform variates in enrollment order, the first
for the DLT and the second for the response, and a tie under the rule
for the neighboring doses draws one more. The 'boinet' package draws its
random numbers differently (event times from a copula, and arrival
times), so the trials agree with it in distribution, not one by one.

Possible values of `stop_reason` are `"no_admissible_dose"`,
`"no_lower_admissible_dose"`, `"n_earlystop"`, `"max_cohorts"` and
`"top_dose_excluded"`; the first two are early stops.

## References

Takeda, K., Taguri, M. and Morita, S. (2018). BOIN-ET: Bayesian optimal
interval design for dose finding based on both efficacy and toxicity
outcomes. Pharmaceutical Statistics, 17(4), 383-395.
[doi:10.1002/pst.1864](https://doi.org/10.1002/pst.1864)

Takeda, K., Zhu, J. and Hirakawa, A. (2025). BF-BOIN-ET: A backfill
Bayesian optimal interval design using efficacy and toxicity outcomes
for dose optimization. Pharmaceutical Statistics, 24, e2470.
[doi:10.1002/pst.2470](https://doi.org/10.1002/pst.2470)

Yamaguchi, Y., Takeda, K. and Maruo, K. (2025). boinet: An R Package for
Bayesian Optimal Interval Design for Dose Finding Based on Both Efficacy
and Toxicity Outcomes. The R Journal, 17(2), 15-34.
[doi:10.32614/RJ-2025-012](https://doi.org/10.32614/RJ-2025-012)

## See also

[`sim_boin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_et.md),
[`boin_et_select_obd`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_et_select_obd.md)

## Examples

``` r
trials <- boin_et_simulate(
  target = 0.30, target_resp = 0.60,
  p_true = c(0.05, 0.15, 0.25, 0.30, 0.35, 0.40),
  p_resp = c(0.05, 0.35, 0.55, 0.65, 0.75, 0.70),
  n_cohort = 12, cohort_size = 3, n_trials = 200, seed = 123
)
trials
#> Simulated BOIN-ET trials
#>   trials               : 200
#>   doses                : 6
#>   target DLT rate      : 0.3
#>   target response rate : 0.6
#>   cohorts              : 12 of size 3
#>   max sample size      : 36
#> 
#> Average per trial:
#>                     DL1  DL2   DL3  DL4  DL5  DL6
#> Patients treated   3.32 9.42 14.22 7.56 1.41 0.07
#> DLTs observed      0.20 1.36  3.62 2.35 0.52 0.04
#> Responses observed 0.17 3.19  7.71 5.10 1.08 0.06
#> 
#> Stopping reason (%):
#> 
#> max_cohorts 
#>         100 

head(trials$n_resp)
#>      DL1 DL2 DL3 DL4 DL5 DL6
#> [1,]   0   0   8   6   3   0
#> [2,]   1   1  15   0   0   0
#> [3,]   0   7  10   0   0   0
#> [4,]   0   1   3   1   7   6
#> [5,]   0   1  17   0   0   0
#> [6,]   0   7   7   0   0   0
table(trials$stop_reason)
#> 
#> max_cohorts 
#>         200 
```
