# Print the Operating Characteristics of Several Designs

Print the table built by
[`oc_results`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_results.md):
for every scenario, the chosen per-dose quantities with one row per
design and one column per dose (for a combination, one matrix per design
with rows for drug A and columns for drug B), followed by the quantities
of the whole trial with one row per design.

## Usage

``` r
# S3 method for class 'oc_results'
print(x, metric = c("sel_percent", "n_pts"), digits = 1, ...)
```

## Arguments

- x:

  An object of class `oc_results`.

- metric:

  Character vector. Per-dose quantities to print, among `"sel_percent"`,
  `"n_pts"`, `"n_tox"` and `"n_resp"`. Defaults to
  `c("sel_percent", "n_pts")`.

- digits:

  Integer scalar. Number of decimals. Defaults to 1.

- ...:

  Further arguments, currently ignored.

## Value

`x`, invisibly.

## See also

[`oc_results`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_results.md),
[`plot.oc_results`](https://gosukehommaex.github.io/simFastBOIN/reference/plot.oc_results.md)

## Examples

``` r
boin <- sim_boin(target = 0.30, p_true = c(0.10, 0.25, 0.40, 0.55),
                 n_cohort = 10, cohort_size = 3, n_trials = 200, seed = 1)
print(oc_results(BOIN = boin, "3+3" = oc_3p3(c(0.10, 0.25, 0.40, 0.55))))
#> Operating characteristics of 2 design(s) in 1 scenario(s)
#> 
#> Scenario 1
#> Selection (%) by dose
#>         1    2    3   4
#> BOIN 11.5 53.0 31.0 4.0
#> 3+3  36.3 37.5 14.8 2.0
#> Patients by dose
#>        1    2   3   4
#> BOIN 7.0 11.8 7.6 1.9
#> 3+3  3.7  3.9 2.3 0.7
#> Per trial (averages)
#>      No selection (%) Patients DLTs
#> BOIN              0.5     28.4  7.9
#> 3+3               9.4     10.6  2.6
```
