# Operating Characteristics of the U-BOIN Design Across Scenarios

Run
[`sim_uboin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_uboin.md)
under several scenarios of true DLT and response probabilities with the
same design, and collect the main operating characteristics of every
scenario in one table.

## Usage

``` r
sim_uboin_multi(
  tox_upper,
  resp_lower,
  scenarios,
  n_cohort,
  cohort_size,
  n_trials = 10000,
  target = NULL,
  resp_assoc = 0,
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
  keep_trials = FALSE,
  verbose = FALSE,
  seed = 123
)
```

## Arguments

- tox_upper:

  Numeric scalar. Highest acceptable DLT probability, for example 0.30.

- resp_lower:

  Numeric scalar. Lowest acceptable response probability, for example
  0.20.

- scenarios:

  List of scenarios, each a list with elements `p_true` and `p_resp` and
  optionally `name`. Names of the list are used as scenario names. All
  scenarios must have the same number of doses.

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
  toxicity rule. Defaults to `c(1, 1)`.

- prior_resp:

  Numeric vector of length two. Beta prior of the response probability
  in the futility rule. Defaults to `c(1, 1)`.

- keep_trials:

  Logical scalar. Keep the trial by trial data in the result of every
  scenario. Defaults to `FALSE`.

- verbose:

  Logical scalar. Report progress while the simulation runs. Defaults to
  `FALSE`.

- seed:

  Integer scalar or `NULL`. Random seed. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `uboin_oc_multi`, which is a list with components

- results:

  Named list of `uboin_oc` objects, one per scenario.

- summary_table:

  Data frame with the true rates, the selection percentages and the
  average numbers of patients, DLTs and responses of every scenario and
  dose.

- scenario_names:

  Character vector of scenario names.

- n_doses:

  Number of doses.

- call:

  The matched call.

## Details

Every scenario is simulated with the same seed, so that the scenarios
differ by their probabilities and not by the random numbers. The outcome
probabilities of a scenario are given by its marginal rates and
`resp_assoc`; to give them directly, run
[`sim_uboin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_uboin.md)
with `p_outcome` for each scenario.

## References

Zhou, Y., Lee, J. J. and Yuan, Y. (2019). A utility-based Bayesian
optimal interval (U-BOIN) phase I/II design to identify the optimal
biological dose for targeted and immune therapies. Statistics in
Medicine, 38(28), S5299-S5316.
[doi:10.1002/sim.8361](https://doi.org/10.1002/sim.8361)

## See also

[`sim_uboin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_uboin.md)

## Examples

``` r
scenarios <- list(
  "Scenario 1" = list(p_true = c(0.02, 0.15, 0.30, 0.45, 0.60),
                      p_resp = c(0.20, 0.65, 0.65, 0.65, 0.65)),
  "Scenario 3" = list(p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
                      p_resp = c(0.08, 0.46, 0.25, 0.20, 0.10))
)

oc <- sim_uboin_multi(
  tox_upper = 0.30,
  resp_lower = 0.20,
  scenarios = scenarios,
  n_cohort = 18,
  cohort_size = 3,
  resp_assoc = 0.2,
  n_trials = 200,
  seed = 123
)
oc
#> U-BOIN operating characteristics across 2 scenarios
#>   DLT upper limit      : 0.3
#>   response lower limit : 0.2
#>   trials each          : 200
#>   max sample size      : 54
#>   utilities            : 0, 30, 50, 100
#> 
#>    Scenario                   Item  DL1  DL2  DL3  DL4  DL5 Total / No OBD
#>  Scenario 1      True DLT rate (%)  2.0 15.0 30.0 45.0 60.0               
#>             True response rate (%) 20.0 65.0 65.0 65.0 65.0               
#>                   OBD selected (%)  1.0 82.0 15.0  2.0  0.0            0.0
#>                   Patients treated  6.0 31.2 12.7  3.6  0.5           54.0
#>                  Patients with DLT  0.1  4.7  3.9  1.6  0.3           10.7
#>             Patients with response  1.2 20.5  8.2  2.3  0.3           32.5
#>  Scenario 3      True DLT rate (%)  5.0 15.0 30.0 45.0 60.0               
#>             True response rate (%)  8.0 46.0 25.0 20.0 10.0               
#>                   OBD selected (%)  0.5 91.5  5.0  0.5  0.0            2.5
#>                   Patients treated  6.9 33.8  9.1  2.9  0.5           53.2
#>                  Patients with DLT  0.3  5.0  2.7  1.3  0.3            9.5
#>             Patients with response  0.5 15.4  2.3  0.5  0.1           18.7
```
