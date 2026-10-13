# Print the Next Dose Combination of a BOIN Drug-Combination Trial

Display the decision returned by
[`comb_boin_next`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin_next.md):
the data at the current combination, the boundaries that apply to it,
the decision, the candidates that were compared and the combinations
eliminated so far.

## Usage

``` r
# S3 method for class 'comb_boin_next'
print(x, ...)
```

## Arguments

- x:

  An object of class `comb_boin_next`.

- ...:

  Further arguments, currently ignored.

## Value

The object `x`, invisibly.

## Examples

``` r
n_pts <- matrix(c(3, 0, 0,
                  0, 0, 0), nrow = 2, byrow = TRUE)
n_tox <- matrix(0, nrow = 2, ncol = 3)

print(comb_boin_next(target = 0.30, n_pts = n_pts, n_tox = n_tox,
                     current = c(1, 1)))
#> BOIN drug-combination design: dose for the next cohort
#>   target DLT rate : 0.3
#>   current         : A1B1, 0 DLTs out of 3 patients (0)
#>   boundaries      : escalate if DLTs <= 0, de-escalate if DLTs >= 2, eliminate if DLTs >= 3 (3 patients)
#> 
#> Decision: escalate to A2B1 or A1B2
#>   The two candidates tie; the design chooses one of them at random.
#> 
#> Candidates (interval 0.2365 to 0.3585):
#>  Direction Combination Patients DLTs Pr(in interval)  Score Available
#>     drug A        A2B1        0    0           0.085 0.0854       yes
#>     drug B        A1B2        0    0           0.085 0.0854       yes
#> 
#> Eliminated combinations: none
```
