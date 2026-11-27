# Simulate Waterfall Drug-Combination Trials

Simulate trials of two drugs given together, run with the waterfall
design of Zhang and Yuan (2016), which looks for the MTD contour, that
is, for one MTD at every dose level of drug A. This is the simulation
engine behind
[`sim_waterfall`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_waterfall.md),
exposed for users who want the raw trial data.

## Usage

``` r
waterfall_simulate(
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

- seed:

  Integer scalar or `NULL`. Random seed. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `waterfall_trials`, which is a list with components

- n_pts:

  Integer array with dimensions trial by level of drug A by level of
  drug B: number of patients treated at each combination.

- n_tox:

  Integer array of the same dimensions: number of DLTs observed.

- eliminated:

  Logical array of the same dimensions: whether the combination is
  excluded from the selection of the MTD contour.

- n_subtrials:

  Integer vector: number of subtrials run in each trial.

- n_enrolled:

  Integer vector: number of patients treated in each trial, including
  those of a subtrial whose data were dropped.

- stop_reason:

  Character vector: why each trial ended.

- boundary:

  The `boin_boundary` object used for the decisions.

- settings:

  List of the design parameters.

## Details

The dose matrix is divided into one-dimensional subtrials, each run with
the BOIN design for a single agent. The first subtrial starts at the
lowest combination, moves up the first column (the lowest level of drug
B) to the highest level of drug A and then along that row. When a
subtrial ends, the dose whose isotonic estimate of the DLT probability
is closest to `target` is its candidate MTD. The next subtrial runs
along the row below the row of the candidate, over the second to the
last level of drug B, and starts one level of drug B above the candidate
(at the last level when the candidate is there). The trial ends after
the subtrial of the lowest level of drug A.

Within a subtrial, the dose escalates, stays or de-escalates by the
boundaries of
[`boin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_boundary.md);
a dose with too many DLTs is eliminated together with the doses above
it, and the subtrial ends without a candidate when its lowest dose is
eliminated. Otherwise a subtrial ends when `n_earlystop` patients have
been treated at the dose chosen for the next cohort, or when its maximum
number of cohorts, `n_cohort`, is reached. The default of `n_cohort` is
`ceiling(round(4 / cohort_size, 2) * m)` for a subtrial of `m` doses,
the default of the reference implementation, which follows the article's
suggestion of four patients per dose.

The rules are those of
[`BOIN::get.oc.comb()`](https://rdrr.io/pkg/BOIN/man/get.oc.comb.html)
with `mtd.contour = TRUE` in version 2.7.2 of the 'BOIN' package, the
implementation by the authors of the article. Random variates are drawn
one per patient in enrollment order, as that function draws them, so
with `boin_compatible = TRUE`, the same seed and matching arguments the
two give the same trials. Three rules of that function are not in the
article:

- When the candidate MTD of the first subtrial lies in the first column
  below the highest level of drug A, the rows above it are excluded, and
  when the data at the candidate would lead to an escalation, an
  additional subtrial explores the rest of its row.

- The combinations to the right of a candidate MTD in its row are
  excluded from the selection of the MTD contour.

- The trial ends when a subtrial ends without a candidate MTD; the rows
  below are then not explored.

With `extrasafe = TRUE`, that function evaluates the safety rule with a
`Beta(y + 0.05, n - y + 0.1)` distribution during a subtrial and with a
`Beta(y + 1, n - y + 1)` distribution when it chooses the candidate; the
package does the same.

The package differs from that function in two points. When a subtrial
other than the additional one ends without a candidate, that function
drops its data, so that its patients and DLTs are missing from the
counts and a trial stopped in the first subtrial has no patients. By
default the package keeps these data and marks every dose of that
subtrial as eliminated; `boin_compatible = TRUE` drops them. And that
function takes the escalation and de-escalation boundaries from the
default values of `p_saf` and `p_tox` whatever values are given; the
package uses the values given.

Values of `stop_reason` are `"completed"`, `"first_subtrial_no_mtd"`,
`"later_subtrial_no_mtd"` and `"max_sample_size"`.

## References

Zhang, L. and Yuan, Y. (2016). A Practical Bayesian Design to Identify
the Maximum Tolerated Dose Contour for Drug Combination Trials.
Statistics in Medicine, 35(27), 4924-4936.

Yan, F., Zhang, L., Zhou, Y., Pan, H., Liu, S. and Yuan, Y. (2020).
BOIN: An R Package for Designing Single-Agent and Drug-Combination
Dose-Finding Trials Using Bayesian Optimal Interval Designs. Journal of
Statistical Software, 94(13), 1-32.

## See also

[`sim_waterfall`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_waterfall.md),
[`waterfall_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/reference/waterfall_select_mtd.md)

## Examples

``` r
p_true <- matrix(c(0.03, 0.10, 0.28,
                   0.10, 0.30, 0.50), nrow = 2, byrow = TRUE)

trials <- waterfall_simulate(
  target = 0.30,
  p_true = p_true,
  cohort_size = 3,
  n_trials = 200,
  seed = 123
)
trials
#> Simulated waterfall drug-combination trials
#>   trials          : 200
#>   dose levels     : 2 of drug A, 3 of drug B
#>   target DLT rate : 0.3
#>   cohorts         : at most 6/3 per subtrial, of size 3
#>   max sample size : 27
#> 
#> Average patients treated per trial:
#>      B1   B2   B3
#> A1 3.35 2.01 6.99
#> A2 5.61 6.62 2.37
#> 
#> Average DLTs observed per trial:
#>      B1   B2   B3
#> A1 0.09 0.17 1.97
#> A2 0.50 2.08 1.20
#> 
#> Reason for ending (%):
#> 
#> completed 
#>       100 
```
