# Operating Characteristics of the iBOIN Design

Simulate trials of the iBOIN design of Zhou et al. (2021), which brings
prior information on the DLT probability of every dose into the BOIN
design, and summarize how often each dose is selected as the MTD and how
the patients are distributed over the doses.

## Usage

``` r
sim_iboin(
  target,
  p_true,
  skeleton,
  prior_ess,
  n_cohort,
  cohort_size,
  n_trials = 10000,
  start_dose = 1,
  n_earlystop = 18,
  p_saf = NULL,
  p_tox = NULL,
  cutoff_eli = 0.95,
  extrasafe = FALSE,
  offset = 0.05,
  robust = FALSE,
  select_with_prior = TRUE,
  min_mtd_sample = 1,
  overdose_cutoff = NULL,
  n_earlystop_rule = c("with_stay", "simple"),
  keep_trials = FALSE,
  verbose = FALSE,
  seed = 123
)
```

## Arguments

- target:

  Numeric scalar. Target DLT probability, for example 0.30.

- p_true:

  Numeric vector. True DLT probability at each dose level, in increasing
  dose order.

- skeleton:

  Numeric vector, one value per dose. Prior DLT probability at each
  dose, strictly between 0 and 1.

- prior_ess:

  Whole number or vector of whole numbers, one per dose. Prior effective
  sample size. See
  [`iboin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/iboin_boundary.md).

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

  Integer scalar. The trial stops once this many patients have been
  treated at the current dose and the design would stay there. Defaults
  to 18. The simulations of the article do not stop early; a value above
  the maximum sample size switches the rule off.

- p_saf:

  Numeric scalar. Highest DLT probability deemed subtherapeutic.
  Defaults to `0.6 * target`.

- p_tox:

  Numeric scalar. Lowest DLT probability deemed overly toxic. Defaults
  to `1.4 * target`.

- cutoff_eli:

  Numeric scalar. Posterior probability cutoff for dose elimination.
  Defaults to 0.95.

- extrasafe:

  Logical scalar. Apply the stricter safety stopping rule at the lowest
  dose. Defaults to `FALSE`.

- offset:

  Numeric scalar between 0 and 0.5. Amount by which `cutoff_eli` is
  relaxed for the safety stopping rule. Defaults to 0.05.

- robust:

  Logical scalar. Use the robust prior, which drops the prior
  information above the prior MTD when the prior MTD lies in the upper
  half of the doses. Defaults to `FALSE`. See
  [`iboin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/iboin_boundary.md).

- select_with_prior:

  Logical scalar. Add the prior information to the estimates of the DLT
  probabilities from which the MTD is selected, as
  [`iboin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/reference/iboin_select_mtd.md)
  does. Defaults to `TRUE`, which reproduces the operating
  characteristics of the article. With `FALSE` the MTD is selected from
  the data alone, by
  [`boin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_select_mtd.md).

- min_mtd_sample:

  Integer scalar. Smallest number of patients a dose must have been
  given to be selected as the MTD. Defaults to 1.

- overdose_cutoff:

  Numeric scalar or `NULL`. Doses with a true DLT probability above this
  value count as overdoses in the summary. Defaults to `NULL`, the
  target.

- n_earlystop_rule:

  Character scalar, either `"with_stay"` or `"simple"`. Under
  `"with_stay"` the trial stops at `n_earlystop` only when the design
  would remain at the current dose. Under `"simple"` it stops as soon as
  the threshold is reached.

- keep_trials:

  Logical scalar. Keep the trial level data in the result. Defaults to
  `FALSE`.

- verbose:

  Logical scalar. Report the progress in messages. Defaults to `FALSE`.

- seed:

  Integer scalar or `NULL`. Random seed. When `NULL` the current state
  of the random number generator is used. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `c("iboin_oc", "boin_oc")`, a list with the
components of the result of
[`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md)
and in addition

- design:

  The name of the design, `"iBOIN"`.

- skeleton:

  The skeleton.

- prior_ess:

  The prior effective sample sizes used, after the robust prior has been
  applied.

- select_with_prior:

  Whether the selection of the MTD used the prior information.

- boundary:

  The `iboin_boundary` object used.

The components `lambda_e` and `lambda_d` are the boundaries of the BOIN
design, which the iBOIN boundaries approach as patients are treated at a
dose.

## Details

The trials are simulated by
[`iboin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/iboin_simulate.md)
with the boundaries of
[`iboin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/iboin_boundary.md).
The article gives two ways to select the MTD at the end of the trial,
with and without the prior information in the estimates of the DLT
probabilities, and does not say which one its simulation study uses. Its
operating characteristics are reproduced with the prior information (see
`inst/reproduce/reproduce_Zhou2021.R`), which is therefore the default.
The selection from the data alone (`select_with_prior = FALSE`) keeps a
misspecified prior from biasing the selection, at the price of a lower
probability of selecting the true MTD when the prior is right.

With a prior effective sample size of 0 at every dose the result is, for
the same seed, that of
[`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md).

The result can be passed to
[`oc_results`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_results.md)
and
[`oc_compare`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_compare.md),
where it appears under the name iBOIN.

## References

Zhou, Y., Lee, J. J., Wang, S., Bailey, S. and Yuan, Y. (2021).
Incorporating historical information to improve phase I clinical trials.
Pharmaceutical Statistics, 20, 1017-1034.
[doi:10.1002/pst.2121](https://doi.org/10.1002/pst.2121)

## See also

[`iboin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/iboin_boundary.md),
[`iboin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/iboin_simulate.md),
[`sim_iboin_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_iboin_multi.md),
[`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md)

## Examples

``` r
# Scenario 3 of Zhou et al. (2021): the prior MTD is the true MTD
oc <- sim_iboin(
  target = 0.30,
  p_true = c(0.08, 0.15, 0.31, 0.45, 0.55),
  skeleton = c(0.10, 0.19, 0.30, 0.42, 0.54),
  prior_ess = 3,
  n_cohort = 10,
  cohort_size = 3,
  n_earlystop = 100,
  n_trials = 500,
  seed = 123
)
oc
#> iBOIN operating characteristics
#>   target DLT rate : 0.3
#>   trials          : 500
#>   cohorts         : 10 of size 3
#>   max sample size : 30
#>   skeleton        : 0.10 0.19 0.30 0.42 0.54
#>   prior ESS       : 3 3 3 3 3
#> 
#>                   DL1  DL2  DL3  DL4  DL5 Total / No MTD
#> True DLT rate (%) 8.0 15.0 31.0 45.0 55.0               
#> MTD selected (%)  1.6 23.8 62.2 12.0  0.2            0.2
#> Patients treated  3.8  9.4 12.8  3.6  0.4           29.9
#> Patients with DLT 0.3  1.5  4.0  1.7  0.2            7.6
#> 
#> Doses above a true DLT rate of 0.3: DL3, DL4, DL5
#>   Patients treated there              : 55.7%
#>   Trials treating any patient there   : 95.6%
#>   Trials treating over 60% there      : 47.2%
#>   Trials treating over 80% there      : 0%
#>   Trials selecting an MTD there       : 74.4%
```
