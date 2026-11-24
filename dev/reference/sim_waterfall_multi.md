# Operating Characteristics of the Waterfall Design Across Scenarios

Run
[`sim_waterfall`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_waterfall.md)
under several matrices of true DLT probabilities with the same design,
and collect the main operating characteristics of every scenario in one
table.

## Usage

``` r
sim_waterfall_multi(
  target,
  scenarios,
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

- scenarios:

  List of scenarios. Each element is either a matrix of true DLT
  probabilities, with one row per dose level of drug A and one column
  per dose level of drug B, or a list with elements `name` and `p_true`.
  Names of the list are used as scenario names. All matrices must have
  the same dimensions.

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

An object of class `waterfall_oc_multi`, which is a list with components

- results:

  Named list of `waterfall_oc` objects, one per scenario.

- summary_table:

  Data frame with one row per scenario: the true MTD contour, the
  percentage of trials selecting exactly that contour, the percentages
  of patients treated on and above it, the percentage of trials without
  any MTD, and the average numbers of patients and DLTs.

- scenario_names:

  Character vector of scenario names.

- n_rows, n_cols:

  Numbers of dose levels of drug A and drug B.

- call:

  The matched call.

## Details

Every scenario is simulated with the same seed, so that the scenarios
differ by their DLT probabilities and not by the random numbers.

## References

Zhang, L. and Yuan, Y. (2016). A Practical Bayesian Design to Identify
the Maximum Tolerated Dose Contour for Drug Combination Trials.
Statistics in Medicine, 35(27), 4924-4936.

## See also

[`sim_waterfall`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_waterfall.md)

## Examples

``` r
scenarios <- list(
  "Contour A1B3, A2B2" = matrix(c(0.03, 0.10, 0.28,
                                  0.10, 0.30, 0.50), nrow = 2, byrow = TRUE),
  "Contour A1B1" = matrix(c(0.30, 0.40, 0.50,
                            0.42, 0.49, 0.55), nrow = 2, byrow = TRUE)
)

oc <- sim_waterfall_multi(
  target = 0.30,
  scenarios = scenarios,
  cohort_size = 3,
  n_trials = 200,
  seed = 123
)
oc
#> Waterfall drug-combination operating characteristics across 2 scenarios
#>   target DLT rate : 0.3
#>   dose levels     : 2 of drug A, 3 of drug B
#>   trials each     : 200
#>   max sample size : 27
#> 
#>            Scenario MTD contour Contour selected (%) Patients on contour (%)
#>  Contour A1B3, A2B2  A1B3, A2B2                 50.0                    50.5
#>        Contour A1B1        A1B1                 55.5                    52.0
#>  Patients above contour (%) No MTD (%) Patients DLTs
#>                         8.8          0     26.9  6.0
#>                        48.0         11     18.2  6.6

# The full results of one scenario
oc$results[["Contour A1B1"]]
#> Waterfall drug-combination operating characteristics
#>   target DLT rate : 0.3
#>   dose levels     : 2 of drug A, 3 of drug B
#>   trials          : 200
#>   cohorts         : at most 6/3 per subtrial, of size 3
#>   max sample size : 27
#>   lambda_e / _d   : 0.2365 / 0.3585
#> 
#>               Item Drug A   B1   B2   B3
#>  True DLT rate (%)     A1 30.0 40.0 50.0
#>                        A2 42.0 49.0 55.0
#>   MTD selected (%)     A1 62.0 21.5  5.5
#>                        A2 20.5  4.5  0.5
#>   Patients treated     A1  9.5  2.8  0.8
#>                        A2  4.0  1.0  0.1
#>  Patients with DLT     A1  2.7  1.1  0.4
#>                        A2  1.8  0.5  0.1
#> 
#> No MTD selected at any level of drug A : 11%
#> Average number of patients             : 18.2
#> Average number of DLTs                 : 6.6
#> Average number of subtrials            : 1.4
#> 
#> True MTD contour                       : A1B1
#>   Trials selecting exactly the contour : 55.5%
#>   Patients treated on the contour      : 52%
#>   Patients treated above the contour   : 48%
#>   Patients treated below the contour   : 0%
```
