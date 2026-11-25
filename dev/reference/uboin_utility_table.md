# Utility Table of the U-BOIN Design

Tabulate, for a dose with a given number of patients, the posterior mean
utility and the admissibility of the dose for every possible number of
responses, of DLTs and of responses without a DLT. These are the
decision tables of the U-BOIN design of Zhou, Lee and Yuan (2019), which
can be written into the protocol before the trial starts.

## Usage

``` r
uboin_utility_table(
  n,
  tox_upper,
  resp_lower,
  utility = c(0, 30, 50, 100),
  cutoff_eli = 0.95,
  cutoff_eli_resp = 0.9,
  prior_tox = c(1, 1),
  prior_resp = c(1, 1)
)
```

## Arguments

- n:

  Integer scalar. Number of patients treated at the dose.

- tox_upper:

  Numeric scalar. Highest acceptable DLT probability.

- resp_lower:

  Numeric scalar. Lowest acceptable response probability.

- utility:

  Numeric vector of length four. Utility of the outcomes (no response,
  DLT), (no response, no DLT), (response, DLT) and (response, no DLT).
  Defaults to `c(0, 30, 50, 100)`.

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

A data frame with one row per possible outcome of the `n` patients and
columns

- n_resp:

  Number of responses.

- n_tox:

  Number of DLTs.

- n_resp_no_tox:

  Number of responses without a DLT.

- utility:

  Posterior mean utility of the dose.

- toxic:

  Whether the dose is toxic.

- futile:

  Whether the dose is futile.

- admissible:

  Whether the dose is neither toxic nor futile.

## Details

The posterior mean utility is `sum(utility * (1/4 + n_k)) / (1 + n)`,
with `n_k` the patients with each outcome, as in
[`uboin_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/uboin_simulate.md).
A dose is toxic when the posterior probability that its DLT probability
exceeds `tox_upper` is above `cutoff_eli`, and futile when the posterior
probability that its response probability is below `resp_lower` is above
`cutoff_eli_resp`.

Tables S4 to S6 of the supporting information of the article, for 3, 6
and 9 patients and the utilities `c(0, 15, 25, 100)`, are reproduced
with `prior_resp = c(0.05, 0.05)`; with the default Beta(1, 1) prior a
dose without any response among 3, 6 or 9 patients is admissible,
whereas those tables mark it inadmissible.

## References

Zhou, Y., Lee, J. J. and Yuan, Y. (2019). A utility-based Bayesian
optimal interval (U-BOIN) phase I/II design to identify the optimal
biological dose for targeted and immune therapies. Statistics in
Medicine, 38(28), S5299-S5316.
[doi:10.1002/sim.8361](https://doi.org/10.1002/sim.8361)

## See also

[`uboin_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/uboin_simulate.md),
[`uboin_select_obd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/uboin_select_obd.md)

## Examples

``` r
# Table S4 of the supporting information of the article
tab <- uboin_utility_table(3, tox_upper = 0.30, resp_lower = 0.20,
                           utility = c(0, 15, 25, 100),
                           prior_resp = c(0.05, 0.05))
tab[tab$admissible, ]
#>    n_resp n_tox n_resp_no_tox utility toxic futile admissible
#> 5       1     0             1   41.25 FALSE  FALSE       TRUE
#> 6       1     1             0   22.50 FALSE  FALSE       TRUE
#> 7       1     1             1   37.50 FALSE  FALSE       TRUE
#> 8       1     2             0   18.75 FALSE  FALSE       TRUE
#> 9       1     2             1   33.75 FALSE  FALSE       TRUE
#> 11      2     0             2   62.50 FALSE  FALSE       TRUE
#> 12      2     1             1   43.75 FALSE  FALSE       TRUE
#> 13      2     1             2   58.75 FALSE  FALSE       TRUE
#> 14      2     2             0   25.00 FALSE  FALSE       TRUE
#> 15      2     2             1   40.00 FALSE  FALSE       TRUE
#> 17      3     0             3   83.75 FALSE  FALSE       TRUE
#> 18      3     1             2   65.00 FALSE  FALSE       TRUE
#> 19      3     2             1   46.25 FALSE  FALSE       TRUE
```
