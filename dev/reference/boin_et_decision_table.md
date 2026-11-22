# Decision Table for the BOIN-ET Design

Tabulate the decision of the BOIN-ET design for every number of DLTs and
responses at the current dose, for the given numbers of patients, in the
form of Table 2 of Takeda, Taguri and Morita (2018). The table also
marks the counts for which the dose is no longer admissible.

## Usage

``` r
boin_et_decision_table(
  target,
  target_resp,
  n = c(3, 6, 9, 12),
  p_saf = NULL,
  p_tox = NULL,
  p_saf_resp = NULL,
  cutoff_eli = 0.95,
  cutoff_eli_resp = 0.99,
  p_eli = NULL,
  p_eli_resp = NULL,
  prior = rep(1/6, 6)
)
```

## Arguments

- target:

  Numeric scalar. Target DLT probability, for example 0.30.

- target_resp:

  Numeric scalar. Target response probability, for example 0.60.

- n:

  Integer vector. Numbers of patients at the current dose to tabulate.
  Defaults to `c(3, 6, 9, 12)`.

- p_saf:

  Numeric scalar. Highest DLT probability deemed subtherapeutic.
  Defaults to `0.1 * target`, the value recommended by Takeda et al.
  (2018), which is far below the `0.6 * target` of the BOIN design.

- p_tox:

  Numeric scalar. Lowest DLT probability deemed overly toxic. Defaults
  to `1.4 * target`.

- p_saf_resp:

  Numeric scalar. Highest response probability deemed subtherapeutic.
  Defaults to `0.6 * target_resp`.

- cutoff_eli:

  Numeric scalar. A dose is no longer admissible when the posterior
  probability that its DLT rate is below `p_eli` falls below
  `1 - cutoff_eli`. Defaults to 0.95.

- cutoff_eli_resp:

  Numeric scalar. A dose is no longer admissible when the posterior
  probability that its response rate is above `p_eli_resp` falls below
  `1 - cutoff_eli_resp`. Defaults to 0.99.

- p_eli:

  Numeric scalar or `NULL`. DLT probability of the rule above. Defaults
  to `NULL`, which uses `target` as the 'boinet' package does. Takeda et
  al. (2018) used 0.40 with a cutoff of 0.90.

- p_eli_resp:

  Numeric scalar or `NULL`. Response probability of the rule above.
  Defaults to `NULL`, which uses `p_saf_resp` as the 'boinet' package
  does. Takeda et al. (2018) used 0.20 with a cutoff of 0.90.

- prior:

  Numeric vector of length 6. Prior probabilities of the six hypotheses.
  Defaults to `rep(1 / 6, 6)`, the noninformative choice of the article.
  The third hypothesis, under which every decision is correct, does not
  enter the criterion.

## Value

An object of class `boin_et_decision_table`, which is a data frame with
one row per combination of `n` patients, `n_tox` DLTs and `n_resp`
responses, and the columns

- n, n_tox, n_resp:

  The counts.

- decision:

  `"E"` to escalate, `"S"` to stay, `"D"` to de-escalate and `"E/S/D"`
  when the next dose is chosen among the current dose and its neighbors
  by their observed response rates.

- tox_excluded, resp_excluded:

  Whether the DLT data or the response data make the dose inadmissible.

The boundaries are attached as the attribute `"boundary"`, an object of
class
[`boin_et_boundary`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_boundary.md).

## Details

When the observed DLT rate is at most `lambda1` and the observed
response rate at most `eta1`, the dose is escalated; when the DLT rate
is below `lambda2` and the response rate above `eta1`, the next cohort
stays; when the DLT rate is at least `lambda2`, the dose is
de-escalated. Otherwise the dose is escalated when the next dose has not
been tried, and the next cohort goes to the dose with the highest
observed response rate among the current dose and its neighbors, with
ties broken at random. An inadmissible dose is skipped by
[`boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_simulate.md).

## References

Takeda, K., Taguri, M. and Morita, S. (2018). BOIN-ET: Bayesian optimal
interval design for dose finding based on both efficacy and toxicity
outcomes. Pharmaceutical Statistics, 17(4), 383-395.
[doi:10.1002/pst.1864](https://doi.org/10.1002/pst.1864)

## See also

[`boin_et_boundary`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_boundary.md),
[`boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_simulate.md)

## Examples

``` r
# Table 2 of Takeda et al. (2018)
tab <- boin_et_decision_table(target = 0.30, target_resp = 0.60, n = c(3, 6))
tab
#> BOIN-ET decision table
#>   lambda1 = 0.14, lambda2 = 0.35, eta1 = 0.48
#> 
#> 3 patients
#>     Responses
#> DLTs 0     1     2  3 
#>    0 E     E     S  S 
#>    1 E/S/D E/S/D S  S 
#>    2 D     D     D  D 
#>    3 D*    D*    D* D*
#> 
#> 6 patients
#>     Responses
#> DLTs 0     1     2     3  4  5  6 
#>    0 E     E     E     S  S  S  S 
#>    1 E/S/D E/S/D E/S/D S  S  S  S 
#>    2 E/S/D E/S/D E/S/D S  S  S  S 
#>    3 D     D     D     D  D  D  D 
#>    4 D*    D*    D*    D* D* D* D*
#>    5 D*    D*    D*    D* D* D* D*
#>    6 D*    D*    D*    D* D* D* D*
#> 
#> E = escalate, S = stay, D = de-escalate, E/S/D = choose among the current and neighboring doses;
#> * = the dose is no longer admissible

head(as.data.frame(tab))
#>   n n_tox n_resp decision tox_excluded resp_excluded
#> 1 3     0      0        E        FALSE         FALSE
#> 2 3     0      1        E        FALSE         FALSE
#> 3 3     0      2        S        FALSE         FALSE
#> 4 3     0      3        S        FALSE         FALSE
#> 5 3     1      0    E/S/D        FALSE         FALSE
#> 6 3     1      1    E/S/D        FALSE         FALSE
```
