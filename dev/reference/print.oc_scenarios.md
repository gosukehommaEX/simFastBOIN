# Print True Scenarios

Print the scenarios tabulated by
[`oc_scenarios`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_scenarios.md):
for a single agent, one block per scenario with the DLT and response
probabilities and the utility of every dose and marks under the true
MTDs and OBDs; for a combination, the same quantities as matrices (rows
for drug A, columns for drug B) with the true MTDs and OBDs marked by an
asterisk.

## Usage

``` r
# S3 method for class 'oc_scenarios'
print(x, digits = 2, ...)
```

## Arguments

- x:

  An object of class `oc_scenarios`.

- digits:

  Integer scalar. Number of decimals of the probabilities; the utility
  is printed with one decimal. Defaults to 2.

- ...:

  Further arguments, currently ignored.

## Value

`x`, invisibly.

## See also

[`oc_scenarios`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_scenarios.md),
[`plot.oc_scenarios`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/plot.oc_scenarios.md)

## Examples

``` r
truth <- oc_scenarios(
  list(list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
            p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38))),
  target = 0.35, target_resp = 0.25
)
print(truth)
#> True scenarios: 1 scenario(s), 5 doses 
#> True MTD: DLT probability closest to 0.35 
#> True OBD: highest utility among the doses with DLT probability <= 0.35 and response probability >= 0.25 
#> Utilities: 100 (response, no DLT), 40 (neither), 60 (both), 0 (DLT, no response)
#> 
#> Scenario 1
#>             1    2    3    4    5
#> p_true   0.05 0.12 0.27 0.35 0.50
#> p_resp   0.20 0.35 0.36 0.37 0.38
#> utility  50.0 56.2 50.8 48.2 42.8
#> true MTD                   *     
#> true OBD         *               
```
