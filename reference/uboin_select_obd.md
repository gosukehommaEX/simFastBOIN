# Select the OBD from Completed U-BOIN Trials

Choose the optimal biological dose (OBD) of the U-BOIN design for each
of a set of completed trials from the numbers of patients with each
outcome at each dose.

## Usage

``` r
uboin_select_obd(
  n_outcome,
  tox_upper,
  resp_lower,
  eliminated = NULL,
  utility = c(0, 30, 50, 100),
  cutoff_eli = 0.95,
  cutoff_eli_resp = 0.9,
  prior_tox = c(1, 1),
  prior_resp = c(1, 1)
)
```

## Arguments

- n_outcome:

  Integer matrix with one row per dose and four columns, the numbers of
  patients with the outcomes (no response, DLT), (no response, no DLT),
  (response, DLT) and (response, no DLT) in one trial, or an integer
  array with dimensions trial by dose by outcome for several trials, as
  returned by
  [`uboin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/uboin_simulate.md).

- tox_upper:

  Numeric scalar. Highest acceptable DLT probability.

- resp_lower:

  Numeric scalar. Lowest acceptable response probability.

- eliminated:

  Logical vector with one element per dose for one trial, a logical
  matrix with one row per trial for several, or `NULL`. Doses eliminated
  in stage I, which cannot be selected. Defaults to `NULL`, none.

- utility:

  Numeric vector of length four. Utility of the four outcomes, in the
  order of the columns of `n_outcome`. Defaults to `c(0, 30, 50, 100)`.

- cutoff_eli:

  Numeric scalar. Cutoff of the toxicity rule. Defaults to 0.95.

- cutoff_eli_resp:

  Numeric scalar. Cutoff of the futility rule. Defaults to 0.90.

- prior_tox:

  Numeric vector of length two. Beta prior of the toxicity rule.
  Defaults to `c(1, 1)`.

- prior_resp:

  Numeric vector of length two. Beta prior of the futility rule.
  Defaults to `c(1, 1)`.

## Value

A data frame with one row per trial and columns

- trial:

  Trial index.

- obd:

  Selected dose level, or `NA` when no dose is admissible.

- reason:

  `"selected"` or `"no_admissible_dose"`.

- utility_DL1, ...:

  Posterior mean utility of each dose.

- admissible_DL1, ...:

  Whether each dose is admissible: tried, not eliminated, neither toxic
  nor futile.

## Details

The OBD is the admissible dose with the largest posterior mean utility,
the lowest among ties, with the rules of rule B2 described in
[`uboin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/uboin_simulate.md).

## References

Zhou, Y., Lee, J. J. and Yuan, Y. (2019). A utility-based Bayesian
optimal interval (U-BOIN) phase I/II design to identify the optimal
biological dose for targeted and immune therapies. Statistics in
Medicine, 38(28), S5299-S5316.
[doi:10.1002/sim.8361](https://doi.org/10.1002/sim.8361)

## See also

[`uboin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/uboin_simulate.md),
[`sim_uboin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_uboin.md)

## Examples

``` r
# Table S3 of the supporting information of the article, with its utilities
n_outcome <- rbind(c(0, 3, 0, 0),
                   c(0, 2, 0, 1),
                   c(0, 3, 1, 1),
                   c(1, 3, 1, 4),
                   c(1, 2, 1, 2))
uboin_select_obd(n_outcome, tox_upper = 0.30, resp_lower = 0.20,
                 utility = c(0, 15, 25, 100))
#>   trial obd   reason utility_DL1 utility_DL2 utility_DL3 utility_DL4
#> 1     1   4 selected          20       41.25    34.16667        50.5
#>   utility_DL5 admissible_DL1 admissible_DL2 admissible_DL3 admissible_DL4
#> 1    41.42857           TRUE           TRUE           TRUE           TRUE
#>   admissible_DL5
#> 1           TRUE
```
