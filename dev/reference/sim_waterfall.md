# Operating Characteristics of the Waterfall Drug-Combination Design

Simulate a trial of two drugs given together, run with the waterfall
design of Zhang and Yuan (2016), many times under one matrix of true DLT
probabilities, and summarize how often each dose combination is selected
as an MTD of the contour, how often the whole contour is selected, how
many patients are treated at each combination and how many DLTs are
observed.

## Usage

``` r
sim_waterfall(
  target,
  p_true,
  cohort_size,
  n_cohort = NULL,
  n_trials = 10000,
  n_earlystop = 12,
  p_saf = NULL,
  p_tox = NULL,
  cutoff_eli = 0.95,
  extrasafe = FALSE,
  offset = 0.05,
  titration = FALSE,
  bound_mtd = FALSE,
  boin_compatible = FALSE,
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
  B, both in increasing dose order. It needs at least two rows and no
  more rows than columns.

- cohort_size:

  Integer scalar. Number of patients per cohort.

- n_cohort:

  Integer vector with one element per dose level of drug A, or `NULL`.
  Maximum number of cohorts of the first, second, ... subtrial. Defaults
  to `NULL`, which gives every subtrial about four patients per dose of
  the subtrial (see Details).

- n_trials:

  Integer scalar. Number of trials to simulate. Defaults to 10000.

- n_earlystop:

  Integer scalar. A subtrial ends once this many patients have been
  treated at the dose chosen for its next cohort. Defaults to 12.

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
  dose of every subtrial. Defaults to `FALSE`.

- offset:

  Numeric scalar between 0 and 0.5. Amount by which `cutoff_eli` is
  relaxed for the safety stopping rule. Defaults to 0.05.

- titration:

  Logical scalar. Start the first subtrial with one patient per dose, up
  to the first DLT. Ignored when `cohort_size` is one. Defaults to
  `FALSE`.

- bound_mtd:

  Logical scalar. Require the estimate at the candidate MTD of a
  subtrial to be below the de-escalation boundary. Defaults to `FALSE`.

- boin_compatible:

  Logical scalar. Drop the data of a subtrial that ends without a
  candidate MTD, as
  [`BOIN::get.oc.comb()`](https://rdrr.io/pkg/BOIN/man/get.oc.comb.html)
  does (see Details). Defaults to `FALSE`.

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

An object of class `waterfall_oc`, which is a list with components

- sel_percent:

  Matrix of the percentage of trials selecting each combination as the
  MTD of its level of drug A. A trial can select one combination per
  level of drug A, so the percentages can add up to more than 100.

- percent_no_mtd:

  Percentage of trials selecting no MTD at any level of drug A.

- n_pts_dose:

  Matrix of the average number of patients treated at each combination.

- n_tox_dose:

  Matrix of the average number of DLTs observed at each combination.

- total_n_pts:

  Average total number of patients per trial.

- total_n_tox:

  Average total number of DLTs per trial.

- total_n_enrolled:

  Average number of patients treated per trial, including those of
  subtrials whose data were dropped; it differs from `total_n_pts` only
  with `boin_compatible = TRUE`.

- true_mtd:

  Two-column matrix of the dose levels of drug A and drug B of the true
  MTD contour.

- pct_contour_selected:

  Percentage of trials selecting exactly the true MTD contour; `NA` when
  no level of drug A has a true MTD.

- pct_patients_contour:

  Percentage of all patients treated on the true MTD contour.

- pct_patients_above:

  Percentage of all patients treated above it.

- pct_patients_below:

  Percentage of all patients treated below it.

- avg_n_subtrials:

  Average number of subtrials per trial.

- stop_reason_percent:

  Percentage of trials by reason for ending.

- trials:

  Trial level data when `keep_trials` is `TRUE`, otherwise `NULL`; `mtd`
  is a matrix with one row per trial and one column per level of drug A,
  holding the level of drug B of the MTD.

together with the design parameters and the call.

## Details

The trials are simulated with
[`waterfall_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/waterfall_simulate.md)
and the MTD contour is selected with
[`waterfall_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/waterfall_select_mtd.md);
the rules, and the points in which they differ from the article and from
the reference implementation, are described there. With
`boin_compatible = TRUE`, the same seed and matching arguments, the
results equal those of
[`BOIN::get.oc.comb()`](https://rdrr.io/pkg/BOIN/man/get.oc.comb.html)
with `mtd.contour = TRUE`; with a seed of 6, 1000 trials and the default
`n_cohort`, they are the results of the article.

The true MTD of a level of drug A is the combination of that row whose
DLT probability is closest to `target`, provided it does not exceed
`target + 0.05`. A combination is below the contour when both of its
dose levels are at most those of some true MTD, and above it when both
are at least those of some true MTD; every other combination, including
the true MTDs, counts as on the contour. These are the definitions of
the reference implementation. The three percentages of patients are `NA`
when no level of drug A has a true MTD.

The phase I/II extension of the article is not provided.

## References

Zhang, L. and Yuan, Y. (2016). A Practical Bayesian Design to Identify
the Maximum Tolerated Dose Contour for Drug Combination Trials.
Statistics in Medicine, 35(27), 4924-4936.

Yan, F., Zhang, L., Zhou, Y., Pan, H., Liu, S. and Yuan, Y. (2020).
BOIN: An R Package for Designing Single-Agent and Drug-Combination
Dose-Finding Trials Using Bayesian Optimal Interval Designs. Journal of
Statistical Software, 94(13), 1-32.

## See also

[`sim_waterfall_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_waterfall_multi.md),
[`waterfall_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/waterfall_simulate.md),
[`waterfall_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/waterfall_select_mtd.md),
[`sim_comb_boin`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_comb_boin.md)

## Examples

``` r
# Scenario 1 of Zhang and Yuan (2016): the MTD contour is A1B3 and A2B2
p_true <- matrix(c(0.03, 0.10, 0.28,
                   0.10, 0.30, 0.50), nrow = 2, byrow = TRUE)

oc <- sim_waterfall(
  target = 0.30,
  p_true = p_true,
  cohort_size = 3,
  n_trials = 500,
  seed = 123
)
oc
#> Waterfall drug-combination operating characteristics
#>   target DLT rate : 0.3
#>   dose levels     : 2 of drug A, 3 of drug B
#>   trials          : 500
#>   cohorts         : at most 6/3 per subtrial, of size 3
#>   max sample size : 27
#>   lambda_e / _d   : 0.2365 / 0.3585
#> 
#>               Item Drug A   B1   B2   B3
#>  True DLT rate (%)     A1  3.0 10.0 28.0
#>                        A2 10.0 30.0 50.0
#>   MTD selected (%)     A1  0.0 16.2 83.8
#>                        A2 23.8 58.2 17.6
#>   Patients treated     A1  3.4  1.9  7.1
#>                        A2  5.4  6.6  2.5
#>  Patients with DLT     A1  0.1  0.2  2.0
#>                        A2  0.5  2.0  1.3
#> 
#> No MTD selected at any level of drug A : 0%
#> Average number of patients             : 26.9
#> Average number of DLTs                 : 6.1
#> Average number of subtrials            : 2
#> 
#> True MTD contour                       : A1B3, A2B2
#>   Trials selecting exactly the contour : 48.6%
#>   Patients treated on the contour      : 51%
#>   Patients treated above the contour   : 9.3%
#>   Patients treated below the contour   : 39.7%

# \donttest{
# The settings of the article
oc_article <- sim_waterfall(
  target = 0.30,
  p_true = p_true,
  cohort_size = 3,
  n_trials = 1000,
  boin_compatible = TRUE,
  seed = 6
)
oc_article$sel_percent
#>      B1   B2   B3
#> A1  0.3 15.5 84.2
#> A2 21.1 59.8 18.5
# }
```
