# Simulate BOIN Drug-Combination Trials

Run the BOIN drug-combination design of Lin and Yin (2017) over many
simulated trials of two drugs given together, and return the patient and
DLT counts at every dose combination. This is the simulation engine
behind
[`sim_comb_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin.md),
exposed for users who want the raw trial data.

## Usage

``` r
comb_boin_simulate(
  target,
  p_true,
  n_cohort,
  cohort_size,
  n_trials = 10000,
  start_dose = c(1, 1),
  n_earlystop = 18,
  p_saf = NULL,
  p_tox = NULL,
  cutoff_eli = 0.95,
  extrasafe = FALSE,
  offset = 0.05,
  titration = FALSE,
  seed = 123
)
```

## Arguments

- target:

  Numeric scalar. Target DLT probability, for example 0.30.

- p_true:

  Numeric matrix. True DLT probability of each dose combination, with
  one row per dose level of drug A and one column per dose level of drug
  B, both in increasing dose order.

- n_cohort:

  Integer scalar. Number of cohorts in a trial.

- cohort_size:

  Integer scalar or vector. Number of patients per cohort. A scalar is
  used for every cohort. A vector shorter than `n_cohort` is padded with
  its last element and a longer one is truncated.

- n_trials:

  Integer scalar. Number of trials to simulate. Defaults to 10000.

- start_dose:

  Integer vector of length two. Dose levels of drug A and drug B for the
  first cohort. Defaults to `c(1, 1)`. Titration, when requested, also
  starts here.

- n_earlystop:

  Integer scalar. The trial stops once this many patients have been
  treated at the current combination and the design would neither
  escalate nor de-escalate, or cannot. Defaults to 18. Set it to a value
  above the maximum sample size to switch this rule off.

- p_saf:

  Numeric scalar. Highest DLT probability deemed subtherapeutic.
  Defaults to `0.6 * target`.

- p_tox:

  Numeric scalar. Lowest DLT probability deemed overly toxic. Defaults
  to `1.4 * target`.

- cutoff_eli:

  Numeric scalar. Posterior probability cutoff for eliminating a
  combination. Defaults to 0.95.

- extrasafe:

  Logical scalar. Apply the stricter safety stopping rule at the lowest
  combination. Defaults to `FALSE`.

- offset:

  Numeric scalar between 0 and 0.5. Amount by which `cutoff_eli` is
  relaxed for the safety stopping rule. Defaults to 0.05.

- titration:

  Logical scalar. Treat one patient per combination, moving up one level
  of either drug, until the first DLT, before switching to full cohorts.
  Ignored when the first cohort size is one. Defaults to `FALSE`.

- seed:

  Integer scalar or `NULL`. Random seed. When `NULL` the current state
  of the random number generator is used. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `comb_boin_trials`, which is a list with components

- n_pts:

  Integer array of patients treated, with dimensions trial by level of
  drug A by level of drug B.

- n_tox:

  Integer array of the same shape, of DLTs observed.

- eliminated:

  Logical array of the same shape, marking combinations eliminated
  during the trial.

- cohorts_used:

  Integer vector of cohorts completed in each trial.

- stop_reason:

  Character vector giving why each trial ended.

- boundary:

  The `boin_boundary` object used.

- settings:

  List of the design parameters.

## Details

The decision boundaries are those of the single-agent design, from
[`boin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_boundary.md).
After each cohort the observed DLT rate at the current combination is
compared with them. When the design escalates, the next cohort goes one
level up of drug A or of drug B, whichever has the larger posterior
probability, under a Beta(0.5, 0.5) prior, that its DLT probability lies
between `lambda_e` and `lambda_d`; de-escalation chooses one level down
in the same way. Eliminated combinations are not escalated to. When a
combination is eliminated, so is every combination with at least the
same level of both drugs, and the trial stops when the lowest
combination is eliminated.

The rules, and the order in which they are applied, are those of
[`BOIN::get.oc.comb()`](https://rdrr.io/pkg/BOIN/man/get.oc.comb.html)
with `mtd.contour = FALSE`. That implementation adds two things to the
article: a candidate's posterior probability is increased by 0.0005 for
every patient already treated there, so that ties go to the better
informed combination, and the trial can stop at `n_earlystop`. Random
variates are drawn in the same order as there: one per patient, one for
the direction of every titration step that can go either way, and one
for every escalation or de-escalation that has a candidate, even when
there is only one. With the same seed and matching arguments a run
therefore reproduces the reference implementation trial by trial.

During titration the patients treated before the first DLT are added to
the planned cohorts, as in the reference implementation: the first full
cohort only tops up the combination reached. A trial can therefore treat
a few more patients than `n_cohort` cohorts hold.

Possible values of `stop_reason` are `"lowest_dose_eliminated"`,
`"lowest_dose_too_toxic"`, `"n_earlystop"` and `"max_cohorts"`.

## References

Lin, R. and Yin, G. (2017). Bayesian Optimal Interval Design for Dose
Finding in Drug-Combination Trials. Statistical Methods in Medical
Research, 26(5), 2155-2167.

Yan, F., Zhang, L., Zhou, Y., Pan, H., Liu, S. and Yuan, Y. (2020).
BOIN: An R Package for Designing Single-Agent and Drug-Combination
Dose-Finding Trials Using Bayesian Optimal Interval Designs. Journal of
Statistical Software, 94(13), 1-32.

## See also

[`sim_comb_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin.md),
[`comb_boin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin_select_mtd.md)

## Examples

``` r
p_true <- matrix(c(0.02, 0.04, 0.08, 0.14,
                   0.08, 0.25, 0.42, 0.48,
                   0.25, 0.45, 0.50, 0.60), nrow = 3, byrow = TRUE)

trials <- comb_boin_simulate(
  target = 0.25,
  p_true = p_true,
  n_cohort = 16,
  cohort_size = 3,
  n_trials = 200,
  seed = 123
)
trials
#> Simulated BOIN drug-combination trials
#>   trials          : 200
#>   dose levels     : 3 of drug A, 4 of drug B
#>   target DLT rate : 0.25
#>   cohorts         : 16 of size 3
#>   max sample size : 48
#> 
#> Average patients treated per trial:
#>      B1   B2   B3   B4
#> A1 3.98 3.70 2.71 2.84
#> A2 5.80 9.03 3.93 1.50
#> A3 4.74 4.04 1.03 0.30
#> 
#> Average DLTs observed per trial:
#>      B1   B2   B3   B4
#> A1 0.07 0.18 0.18 0.41
#> A2 0.40 2.12 1.65 0.72
#> A3 1.22 1.73 0.50 0.20
#> 
#> Stopping reason (%):
#> 
#> max_cohorts n_earlystop 
#>        49.5        50.5 

# Patients treated at each combination in the first trial
trials$n_pts[1, , ]
#>    B1 B2 B3 B4
#> A1  3  3 15  0
#> A2  0  0 12  9
#> A3  0  0  0  6
```
