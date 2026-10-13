# Print BOIN Decision Boundaries

Display the interval boundaries and the integer decision boundaries
produced by
[`boin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_boundary.md).

## Usage

``` r
# S3 method for class 'boin_boundary'
print(x, cohort_size = NULL, ...)
```

## Arguments

- x:

  An object of class `boin_boundary`.

- cohort_size:

  Integer scalar or `NULL`. When supplied, only the sample sizes that
  are multiples of `cohort_size` are shown, which is how the boundaries
  are consulted during a trial with equally sized cohorts.

- ...:

  Further arguments, currently ignored.

## Value

The object `x`, invisibly.

## Examples

``` r
bd <- boin_boundary(target = 0.30, max_n = 18, extrasafe = TRUE)
print(bd, cohort_size = 3)
#> BOIN decision boundaries
#>   target DLT rate : 0.3
#>   p_saf / p_tox   : 0.18 / 0.42
#>   lambda_e        : 0.2365
#>   lambda_d        : 0.3585
#>   cutoff_eli      : 0.95
#>   safety cutoff   : 0.9
#> 
#>                                                  
#> Number of patients treated         3 6 9 12 15 18
#> Escalate if # of DLT <=            0 1 2  2  3  4
#> Deescalate if # of DLT >=          2 3 4  5  6  7
#> Eliminate if # of DLT >=           3 4 5  7  8  9
#> Stop at lowest dose if # of DLT >= 2 4 5  6  7  8
```
