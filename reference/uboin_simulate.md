# Simulate U-BOIN Trials

Simulate trials of the utility-based BOIN (U-BOIN) phase I/II design of
Zhou, Lee and Yuan (2019), which looks for the optimal biological dose
(OBD) from a binary toxicity and a binary efficacy outcome observed
before the next cohort is treated. This is the simulation engine behind
[`sim_uboin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_uboin.md),
exposed for users who want the raw trial data.

## Usage

``` r
uboin_simulate(
  tox_upper,
  resp_lower,
  p_true = NULL,
  p_resp = NULL,
  n_cohort,
  cohort_size,
  n_trials = 10000,
  target = NULL,
  resp_assoc = 0,
  p_outcome = NULL,
  utility = c(0, 30, 50, 100),
  n_stage1 = 12,
  n_stop = NULL,
  start_dose = 1,
  p_saf = NULL,
  p_tox = NULL,
  cutoff_eli = 0.95,
  cutoff_eli_resp = 0.9,
  prior_tox = c(1, 1),
  prior_resp = c(1, 1),
  seed = 123
)
```

## Arguments

- tox_upper:

  Numeric scalar. Highest acceptable DLT probability, for example 0.30.

- resp_lower:

  Numeric scalar. Lowest acceptable response probability, for example
  0.20.

- p_true:

  Numeric vector or `NULL`. True DLT probability of each dose.

- p_resp:

  Numeric vector or `NULL`. True response probability of each dose.

- n_cohort:

  Integer scalar. Number of cohorts in a trial; the maximum sample size
  is `n_cohort * cohort_size`.

- cohort_size:

  Integer scalar. Number of patients per cohort.

- n_trials:

  Integer scalar. Number of trials to simulate. Defaults to 10000.

- target:

  Numeric scalar or `NULL`. Target DLT probability of the BOIN design in
  stage I. Defaults to `NULL`, which uses `tox_upper - 0.05`.

- resp_assoc:

  Numeric scalar. Association parameter of the Gumbel model that joins
  the DLT and the response of a patient; 0 makes them independent.
  Defaults to 0.

- p_outcome:

  Numeric matrix or `NULL`. Probabilities of the four outcomes of a
  patient at each dose, one row per dose, in the columns (no response,
  DLT), (no response, no DLT), (response, DLT) and (response, no DLT).
  When given, `p_true`, `p_resp` and `resp_assoc` must be `NULL` or are
  not used. Defaults to `NULL`.

- utility:

  Numeric vector of length four. Utility of the four outcomes, in the
  order of `p_outcome`. Defaults to `c(0, 30, 50, 100)`, example 1 of
  the article.

- n_stage1:

  Integer scalar. Stage II starts once a dose has this many patients (s1
  of the article). Defaults to 12.

- n_stop:

  Integer scalar or `NULL`. The trial ends once a dose has this many
  patients (s2 of the article). Defaults to `NULL`, the maximum sample
  size.

- start_dose:

  Integer scalar. Dose level of the first cohort. Defaults to 1.

- p_saf:

  Numeric scalar. Highest DLT probability deemed subtherapeutic, for the
  boundaries of stage I. Defaults to `0.6 * target`.

- p_tox:

  Numeric scalar. Lowest DLT probability deemed overly toxic, for the
  boundaries of stage I. Defaults to `1.4 * target`.

- cutoff_eli:

  Numeric scalar. A dose is toxic when the posterior probability that
  its DLT probability exceeds `tox_upper` is above this value. Defaults
  to 0.95.

- cutoff_eli_resp:

  Numeric scalar. A dose is futile when the posterior probability that
  its response probability is below `resp_lower` is above this value.
  Defaults to 0.90.

- prior_tox:

  Numeric vector of length two. Beta prior of the DLT probability in the
  toxicity rule, which eliminates doses in stage I as well as in
  stage II. Defaults to `c(1, 1)`, the uniform prior of the BOIN design
  that stage I of the article uses.

- prior_resp:

  Numeric vector of length two. Beta prior of the response probability
  in the futility rule. Defaults to `c(1, 1)`.

- seed:

  Integer scalar or `NULL`. Random seed. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `uboin_trials`, which is a list with components

- n_outcome:

  Integer array with dimensions trial by dose by outcome: number of
  patients with each outcome.

- n_pts, n_tox, n_resp:

  Integer matrices with one row per trial and one column per dose:
  patients treated, DLTs and responses.

- eliminated:

  Logical matrix: doses eliminated in stage I.

- cohorts_used:

  Integer vector: cohorts treated in each trial.

- stage2:

  Logical vector: whether each trial reached stage II.

- stop_reason:

  Character vector: why each trial ended.

- boundary:

  List of the boundaries of stage I.

- design:

  `"U-BOIN"`.

- settings:

  List of the design parameters.

## Details

Stage I is the BOIN design with the target `target`. After each cohort,
a dose that is toxic and has at least three patients is eliminated
together with every higher dose, and the trial stops when the lowest
dose is eliminated; the next cohort escalates (unless the next dose is
eliminated), stays or de-escalates by the boundaries. Once a dose has
`n_stage1` patients, stage II starts with the next cohort:

- (B1) If the DLT rate observed at the highest dose tried is at most the
  escalation boundary, that dose is not the highest dose and the dose
  above it is not eliminated, the next cohort is treated at the dose
  above it.

- (B2) Otherwise the next cohort is treated at the admissible dose with
  the largest posterior mean utility, the lowest among ties. A dose is
  admissible when it has been tried, is not eliminated, and is neither
  toxic nor futile. The trial stops, without an OBD, when no dose is
  admissible.

The trial ends at the maximum sample size or when a dose has `n_stop`
patients. The posterior mean utility of a dose is
`sum(utility * (1/4 + n_k)) / (1 + n)`, the utilities weighted by the
posterior means of the outcome probabilities under the Dirichlet(1/4,
1/4, 1/4, 1/4) prior of the article. Utilities are compared after
rounding to ten decimals.

The article does not give the priors of the toxicity and futility rules.
With Beta(1, 1) priors the simulations of the article are reproduced
(`inst/reproduce/reproduce_Zhou2019.R`); the decision tables of its
supporting information are reproduced with `prior_resp = c(0.05, 0.05)`.
The article does not say either whether the admissible doses of rule B2
include doses not yet tried; they do not here, which the simulations of
the article support. The stage I elimination uses the toxicity rule, as
the article describes. A dose that becomes toxic in stage II does not
exclude the doses above it.

One uniform random variate is drawn per patient, in enrollment order,
and the outcome is the first whose cumulative probability exceeds it.

Only the pick-the-winner allocation of the article is provided; the
adaptive and equal randomization and the version for efficacy observed
late (Section 2.5 of the article) are not.

Values of `stop_reason` are `"max_sample_size"`, `"n_stop"`,
`"lowest_dose_eliminated"` and `"no_admissible_dose"`.

## References

Zhou, Y., Lee, J. J. and Yuan, Y. (2019). A utility-based Bayesian
optimal interval (U-BOIN) phase I/II design to identify the optimal
biological dose for targeted and immune therapies. Statistics in
Medicine, 38(28), S5299-S5316.
[doi:10.1002/sim.8361](https://doi.org/10.1002/sim.8361)

## See also

[`sim_uboin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_uboin.md),
[`uboin_select_obd`](https://gosukehommaex.github.io/simFastBOIN/reference/uboin_select_obd.md),
[`uboin_utility_table`](https://gosukehommaex.github.io/simFastBOIN/reference/uboin_utility_table.md)

## Examples

``` r
trials <- uboin_simulate(
  tox_upper = 0.30,
  resp_lower = 0.20,
  p_true = c(0.02, 0.15, 0.30, 0.45, 0.60),
  p_resp = c(0.20, 0.65, 0.65, 0.65, 0.65),
  n_cohort = 18,
  cohort_size = 3,
  resp_assoc = 0.2,
  n_trials = 200,
  seed = 123
)
trials
#> Simulated U-BOIN trials
#>   trials               : 200
#>   doses                : 5
#>   DLT upper limit      : 0.3
#>   response lower limit : 0.2
#>   cohorts              : 18 of size 3
#>   max sample size      : 54
#> 
#> Averages per trial:
#>                     DL1   DL2   DL3  DL4  DL5
#> Patients treated   6.00 31.22 12.70 3.57 0.51
#> DLTs observed      0.14  4.74  3.92 1.59 0.30
#> Responses observed 1.18 20.48  8.18 2.35 0.30
#> 
#> Trials reaching stage II: 100%
#> Reason for ending (%):
#> 
#>        max_sample_size                 n_stop lowest_dose_eliminated 
#>                    100                      0                      0 
#>     no_admissible_dose 
#>                      0 
```
