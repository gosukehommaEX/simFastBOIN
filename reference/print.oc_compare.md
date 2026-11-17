# Print a Comparison of Designs

Print the comparison built by
[`oc_compare`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_compare.md):
one block per scenario with one row per design.

## Usage

``` r
# S3 method for class 'oc_compare'
print(x, digits = 1, ...)
```

## Arguments

- x:

  An object of class `oc_compare`.

- digits:

  Integer scalar. Number of decimals. Defaults to 1.

- ...:

  Further arguments, currently ignored.

## Value

`x`, invisibly.

## See also

[`oc_compare`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_compare.md),
[`plot.oc_compare`](https://gosukehommaex.github.io/simFastBOIN/reference/plot.oc_compare.md)

## Examples

``` r
truth <- oc_scenarios(list(c(0.10, 0.25, 0.40, 0.55)), target = 0.30)
boin <- sim_boin(target = 0.30, p_true = c(0.10, 0.25, 0.40, 0.55),
                 n_cohort = 10, cohort_size = 3, n_trials = 200, seed = 1)
print(oc_compare(BOIN = boin, "3+3" = oc_3p3(c(0.10, 0.25, 0.40, 0.55)),
                 truth = truth))
#> Comparison of 2 design(s) with the true doses
#> Correct: the true MTD for designs selecting an MTD, the true OBD for designs selecting an OBD
#> Overdose: true DLT probability above 0.3 
#> 
#> Scenario 1
#>            Correct (%) Pts at correct (%) Overdose sel. (%) Pts at overdose
#> BOIN (MTD)        53.0               41.7              35.0             9.5
#> 3+3 (MTD)         37.5               36.5              16.8             3.0
#>            No selection (%) Patients
#> BOIN (MTD)              0.5     28.4
#> 3+3 (MTD)               9.4     10.6
```
