# Simulate BOIN-ETC Trials

Simulate trials of the BOIN-ETC design of Kakizume et al. (2024), which
finds the optimal dose combination of two drugs given together by
dividing the dose matrix into subtrials along which the DLT probability
increases (the waterfall approach) and running the BOIN-ET design of
Takeda et al. (2018) within each subtrial. Return the patients, DLTs and
responses at every combination. This is the simulation engine behind
[`sim_boin_etc`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_etc.md).

## Usage

``` r
boin_etc_simulate(
  target,
  target_resp,
  p_true,
  p_resp,
  n_cohort,
  cohort_size,
  n_trials = 10000,
  subtrials = c("etc1", "etc2", "etc3"),
  n_earlystop = 6,
  etc1_start = c("diagonal", "below"),
  rule4 = c("paper", "authors_code"),
  p_saf = NULL,
  p_tox = NULL,
  p_saf_resp = NULL,
  prior = rep(1/6, 6),
  cutoff_eli = 0.95,
  p_eli = NULL,
  resp_cor = 0,
  resp_cor_type = c("copula", "binary"),
  seed = 123
)
```

## Arguments

- target:

  Numeric scalar. Target DLT probability, phi.

- target_resp:

  Numeric scalar. Target response probability, delta, as in
  [`sim_boin_et`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_et.md).

- p_true:

  Numeric matrix. True DLT probability of every combination, with one
  row per level of drug A and one column per level of drug B. At least
  two levels of each drug.

- p_resp:

  Numeric matrix of the same size. True response probability of every
  combination.

- n_cohort:

  Integer scalar. Number of cohorts in a trial; with `cohort_size` it
  gives the maximum sample size.

- cohort_size:

  Integer scalar or vector. Number of patients per cohort. A scalar is
  used for every cohort; a shorter vector is padded with its last
  element.

- n_trials:

  Integer scalar. Number of trials to simulate. Defaults to 10000.

- subtrials:

  Character scalar. The division of the dose matrix into subtrials and
  their starting combinations: `"etc1"` (the default), `"etc2"` or
  `"etc3"`, BOIN-ETC1, 2 and 3 of the figures and tables of Kakizume et
  al. (2024). See Details.

- n_earlystop:

  Integer scalar. A subtrial ends when the combination chosen for the
  next cohort has this many patients (n_stop of Kakizume et al., 2024).
  Defaults to 6.

- etc1_start:

  Character scalar. Starting combination of a row subtrial of `"etc1"`:
  `"diagonal"` (the default) one level of drug B below the candidate of
  the previous subtrial, as the text and the figures of Kakizume et
  al. (2024) state; `"below"` the same level of drug B, the rule of the
  code of the authors that produced their published results. See
  Details.

- rule4:

  Character scalar. Choice among the neighboring combinations when the
  DLT rate is between `lambda1` and `lambda2` and the response rate at
  most `eta1`: `"paper"` (the default) the rule of the text;
  `"authors_code"` the rule of the code of the authors, needed to
  reproduce their results for BOIN-ETC1. See Details.

- p_saf, p_tox, p_saf_resp, prior:

  Design parameters of the BOIN-ET boundaries, as in
  [`boin_et_boundary`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_boundary.md).
  The defaults give `(lambda1, lambda2, eta1) = (0.14, 0.35, 0.48)` for
  a target DLT probability of 0.30 and a target response probability of
  0.60.

- cutoff_eli:

  Numeric scalar. A combination is excluded, with the higher
  combinations of its subtrial, when the posterior probability that its
  DLT rate exceeds `p_eli` is above this cutoff. Defaults to 0.95.

- p_eli:

  Numeric scalar or `NULL`. DLT probability of the exclusion rule.
  Defaults to `NULL`, `target`.

- resp_cor:

  Numeric scalar between -1 and 1. Correlation between the DLT and the
  response of a patient, see `resp_cor_type`. Defaults to 0, as in the
  simulations of Kakizume et al. (2024); the BOIN-ET functions default
  to 0.2.

- resp_cor_type:

  Character scalar. `"copula"` (the default) or `"binary"`, as in
  [`sim_boin_et`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_et.md).

- seed:

  Integer scalar or `NULL`. Random seed. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `boin_etc_trials`, which is a list with components

- n_pts, n_tox, n_resp, n_tox_resp:

  Integer arrays with dimensions trial by level of drug A by level of
  drug B, of the patients treated, the DLTs, the responses and the
  patients with both a DLT and a response.

- excluded:

  Logical array of the same dimensions, `TRUE` for the combinations
  excluded by the DLT rule during the trial.

- cohorts_used, n_subtrials:

  Integer vectors of the cohorts enrolled and of the subtrials run.

- stop_reason:

  Character vector: `"lowest_excluded"` (the lowest combination was
  excluded), `"subtrials_completed"` or `"max_cohorts"` (the maximum
  sample size was reached).

- early_stop:

  Logical vector, `TRUE` for the trials stopped because the lowest
  combination was excluded. These trials select no optimal dose
  combination.

- design:

  `"BOIN-ETC1"`, `"BOIN-ETC2"` or `"BOIN-ETC3"`.

- boundary:

  List with `lambda1`, `lambda2` and `eta1`.

- settings:

  List of the design parameters.

## Details

With `L` levels of drug A and `M` of drug B (`AaBb` is the combination
of level `a` of drug A and level `b` of drug B), the first subtrial is
`A1B1, A2B1, ..., ALB1, ALB2, ..., ALBM`, started at `A1B1`. Then:

- `"etc1"`, `"etc2"`:

  one subtrial per row `l = L - 1, ..., 1`, `AlB2, ..., AlBM`. When the
  candidate of the first subtrial is `AlB1`, the rows start at row
  `l - 1` (at row 1 when `l` is 1) and the rows above are skipped.
  `"etc2"` starts every row at its lowest combination. `"etc1"` starts
  the row below the candidate `AlBm` of the previous subtrial at level
  `m - 1` of drug B (`etc1_start = "diagonal"`) or at level `m`
  (`"below"`), and at level 2 when that level is 1, after a candidate in
  the first column and after a subtrial that ended without one.

- `"etc3"`:

  L-shaped subtrials `j = 1, ..., min(L, M) - 1`: column `j + 1` of drug
  B from level 1 of drug A up to level `L - j`, then row `L - j` from
  column `j + 2` to `M`, each started at its lowest combination. With
  `L > M` the shapes run out of columns after `M - 1` of them, which
  still covers the matrix.

Within a subtrial, after each cohort at a combination with DLT rate `p`
and response rate `q`: the combination and the higher combinations of
the subtrial are excluded when its DLT data fail the rule of
`cutoff_eli`; the subtrial ends without a candidate when its lowest
combination is excluded, and the trial ends when that combination is
`A1B1`. Otherwise the next cohort escalates when `p <= lambda1` and
`q <= eta1` (it stays at the highest combination of the subtrial and
below an excluded one), de-escalates when `p >= lambda2` (it stays at
the lowest), stays when `q > eta1`, and otherwise (rule 4) goes to the
next combination up if it has not been treated, or else to the treated,
not excluded combination of the current one and its two neighbors in the
subtrial with the highest response rate, a tie being broken at random.
The subtrial ends when the combination chosen for the next cohort has
`n_earlystop` patients (that combination is the candidate of the
subtrial), and the trial ends when all cohorts have been enrolled or all
subtrials have been run.

Kakizume et al. (2024) do not state the rules in every case; the
following were checked against their code and their published results
(`inst/reproduce/reproduce_Kakizume2024.R`):

- The text of Section 2.2 describes BOIN-ETC2 with the L-shaped
  subtrials and BOIN-ETC3 with the rows started at their lowest
  combination. Their Figure 2 and Supplementary Figures S1 and S2, their
  tables and their code have it the other way round, and so has this
  function.

- The text and the figures start a row of BOIN-ETC1 one level of drug B
  below the previous candidate (`"diagonal"`); the code that produced
  the tables starts it at the same level (`"below"`), and only that rule
  reproduces the published allocation.

- In rule 4, the code of the authors drops the response rate of an
  untreated combination below the current one but not the combination,
  so that the choice moves one combination down (from the current
  combination to the untreated one, or from the one above to the current
  one); it also breaks a tie among the first eligible combinations
  rather than among the tied ones. `rule4 = "authors_code"` applies
  these rules. An untreated combination below the current one only
  happens in BOIN-ETC1, and its published allocation is reproduced only
  with these rules.

- The code excludes only the current combination and may treat an
  excluded lowest combination again; here the text is followed (the
  published results do not depend on the difference). The text does not
  say what happens when `A1B1` is excluded.

The response of a patient does not exclude a combination, as in Kakizume
et al. (2024). Each patient draws two uniform variates in enrollment
order, the first for the DLT and the second for the response; a tie
under rule 4 draws one more after the cohort.

## References

Kakizume, T., Takeda, K., Taguri, M. and Morita, S. (2024). BOIN-ETC: A
Bayesian optimal interval design considering efficacy and toxicity to
identify the optimal dose combinations. Statistical Methods in Medical
Research, 33(4), 716-727.
[doi:10.1177/09622802241236936](https://doi.org/10.1177/09622802241236936)

Takeda, K., Taguri, M. and Morita, S. (2018). BOIN-ET: Bayesian optimal
interval design for dose finding based on both efficacy and toxicity
outcomes. Pharmaceutical Statistics, 17(4), 383-395.
[doi:10.1002/pst.1864](https://doi.org/10.1002/pst.1864)

## See also

[`sim_boin_etc`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_etc.md),
[`boin_etc_select_obd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_etc_select_obd.md),
[`boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_simulate.md)

## Examples

``` r
# Scenario 1 of Kakizume et al. (2024)
p_true <- matrix(c(0.05, 0.10, 0.10, 0.30, 0.30, 0.50), nrow = 2)
p_resp <- matrix(c(0.10, 0.25, 0.25, 0.60, 0.60, 0.65), nrow = 2)
trials <- boin_etc_simulate(
  target = 0.30, target_resp = 0.60, p_true = p_true, p_resp = p_resp,
  n_cohort = 8, cohort_size = 3, n_trials = 200, seed = 123
)
trials
#> Simulated BOIN-ETC1 trials
#>   trials               : 200
#>   dose levels          : 2 of drug A, 3 of drug B
#>   target DLT rate      : 0.3
#>   target response rate : 0.6
#>   cohorts              : 8 of size 3
#>   max sample size      : 24
#>   subtrial end (n_stop): 6
#> 
#> Patients treated, average per trial:
#>      B1   B2   B3
#> A1 3.22 4.16 4.36
#> A2 4.22 5.24 1.32
#> 
#> DLTs observed, average per trial:
#>      B1   B2   B3
#> A1 0.18 0.38 1.33
#> A2 0.44 1.52 0.69
#> 
#> Responses observed, average per trial:
#>      B1   B2   B3
#> A1 0.27 1.01 2.70
#> A2 1.00 3.21 0.81
#> 
#> Average number of subtrials: 2
#> 
#> Stopping reason (%):
#> 
#>         max_cohorts subtrials_completed 
#>                  61                  39 
```
