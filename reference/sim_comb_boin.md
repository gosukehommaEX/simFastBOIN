# Operating Characteristics of the BOIN Drug-Combination Design

Simulate a trial of two drugs given together, run with the BOIN
drug-combination design of Lin and Yin (2017), many times under one
matrix of true DLT probabilities, and summarize how often each dose
combination is selected as the MTD, how many patients are treated at
each combination and how many DLTs are observed.

## Usage

``` r
sim_comb_boin(
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
  bound_mtd = FALSE,
  overdose_cutoff = NULL,
  keep_trials = FALSE,
  verbose = FALSE,
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
  first cohort. Defaults to `c(1, 1)`.

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
  of either drug, until the first DLT. Ignored when the first cohort
  size is one. Defaults to `FALSE`.

- bound_mtd:

  Logical scalar. Require the isotonic estimate at the selected
  combination to be below the de-escalation boundary. Defaults to
  `FALSE`.

- overdose_cutoff:

  Numeric scalar or `NULL`. Combinations whose true DLT probability
  exceeds this value count as overdoses in the `overdose` component of
  the result. Defaults to `NULL`, which uses `target`.

- keep_trials:

  Logical scalar. Keep the full trial by trial data in the result.
  Defaults to `FALSE`.

- verbose:

  Logical scalar. Report progress while the simulation runs. Defaults to
  `FALSE`.

- seed:

  Integer scalar or `NULL`. Random seed. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `comb_boin_oc`, which is a list with components

- sel_percent:

  Matrix of the percentage of trials selecting each combination as the
  MTD.

- percent_no_mtd:

  Percentage of trials ending without an MTD.

- n_pts_dose:

  Matrix of the average number of patients treated at each combination.

- n_tox_dose:

  Matrix of the average number of DLTs observed at each combination.

- total_n_pts:

  Average total number of patients per trial.

- total_n_tox:

  Average total number of DLTs per trial.

- true_mtd:

  Two-column matrix of the dose levels of drug A and drug B of the
  combinations whose true DLT probability is closest to `target`.

- pct_true_mtd_selected:

  Percentage of trials selecting one of these combinations.

- pct_patients_true_mtd:

  Percentage of all patients treated at one of them.

- overdose:

  List describing the exposure to combinations above `overdose_cutoff`,
  with the components of the same list returned by
  [`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md);
  `combinations` gives the dose levels of the combinations above the
  cutoff.

- stop_reason_percent:

  Percentage of trials by reason for stopping.

- trials:

  Trial level data when `keep_trials` is `TRUE`, otherwise `NULL`.

together with the design parameters and the call.

## Details

The trials are simulated with
[`comb_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin_simulate.md)
and the MTD is selected with
[`comb_boin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin_select_mtd.md).
A trial stopped for safety at the lowest combination selects none. Both
follow
[`BOIN::get.oc.comb()`](https://rdrr.io/pkg/BOIN/man/get.oc.comb.html)
with `mtd.contour = FALSE`: with the same seed and matching arguments
the two implementations agree trial by trial. `pct_true_mtd_selected`
and `pct_patients_true_mtd` are the quantities that the reference
implementation reports as the selection percentage of the MTD and the
percentage of patients treated at the MTD.

Note that `n_earlystop` defaults to 18 here, whereas the reference
implementation defaults to 100, which in practice switches the rule off.

Only the design that searches for a single MTD combination is provided;
the waterfall design of Zhang and Yuan (2016), which searches for the
MTD contour, is not.

## References

Lin, R. and Yin, G. (2017). Bayesian Optimal Interval Design for Dose
Finding in Drug-Combination Trials. Statistical Methods in Medical
Research, 26(5), 2155-2167.

Yan, F., Zhang, L., Zhou, Y., Pan, H., Liu, S. and Yuan, Y. (2020).
BOIN: An R Package for Designing Single-Agent and Drug-Combination
Dose-Finding Trials Using Bayesian Optimal Interval Designs. Journal of
Statistical Software, 94(13), 1-32.

Zhang, L. and Yuan, Y. (2016). A Practical Bayesian Design to Identify
the Maximum Tolerated Dose Contour for Drug Combination Trials.
Statistics in Medicine, 35(27), 4924-4936.

## See also

[`sim_comb_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin_multi.md),
[`comb_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin_simulate.md),
[`comb_boin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin_select_mtd.md)

## Examples

``` r
# The example of Yan and colleagues (2020): three levels of drug A, four of
# drug B and a target DLT rate of 0.25
p_true <- matrix(c(0.02, 0.04, 0.08, 0.14,
                   0.08, 0.25, 0.42, 0.48,
                   0.25, 0.45, 0.50, 0.60), nrow = 3, byrow = TRUE)

oc <- sim_comb_boin(
  target = 0.25,
  p_true = p_true,
  n_cohort = 16,
  cohort_size = 3,
  n_trials = 500,
  seed = 123
)
oc
#> BOIN drug-combination operating characteristics
#>   target DLT rate : 0.25
#>   dose levels     : 3 of drug A, 4 of drug B
#>   trials          : 500
#>   cohorts         : 16 of size 3
#>   max sample size : 48
#>   lambda_e / _d   : 0.1968 / 0.2984
#> 
#>               Item Drug A   B1   B2   B3   B4
#>  True DLT rate (%)     A1  2.0  4.0  8.0 14.0
#>                        A2  8.0 25.0 42.0 48.0
#>                        A3 25.0 45.0 50.0 60.0
#>   MTD selected (%)     A1  0.0  1.0  5.0 16.0
#>                        A2  4.2 35.4  8.6  1.0
#>                        A3 23.8  4.8  0.2  0.0
#>   Patients treated     A1  4.0  3.6  3.0  2.9
#>                        A2  5.8  8.4  4.0  1.5
#>                        A3  5.0  3.9  1.0  0.3
#>  Patients with DLT     A1  0.1  0.2  0.2  0.4
#>                        A2  0.4  2.1  1.6  0.7
#>                        A3  1.2  1.8  0.5  0.2
#> 
#> No MTD selected                      : 0%
#> Average number of patients           : 43.4
#> Average number of DLTs               : 9.4
#> 
#> True MTD (closest to the target)     : A3B1, A2B2
#>   Selected as the MTD                : 59.2%
#>   Patients treated there             : 30.9%
#> 
#> Combinations above a true DLT rate of 0.25: 5 of 12
#>   Patients treated there             : 24.6%
#>   Trials treating any patient there  : 82%
#>   Trials treating over 60% there     : 3%
#>   Trials treating over 80% there     : 0.2%
#>   Trials selecting an MTD there      : 14.6%

# \donttest{
# The settings of the article, which turn early stopping off
oc_jss <- sim_comb_boin(
  target = 0.25,
  p_true = p_true,
  n_cohort = 16,
  cohort_size = 3,
  n_earlystop = 100,
  n_trials = 1000,
  seed = 6
)
oc_jss$pct_true_mtd_selected
#> [1] 59.4
# }
```
