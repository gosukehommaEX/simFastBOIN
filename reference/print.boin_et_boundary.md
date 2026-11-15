# Print the Boundaries of the BOIN-ET Design

Display the boundaries found by
[`boin_et_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_et_boundary.md)
together with the design parameters they were found from.

## Usage

``` r
# S3 method for class 'boin_et_boundary'
print(x, ...)
```

## Arguments

- x:

  An object of class `boin_et_boundary`.

- ...:

  Further arguments, currently ignored.

## Value

The object `x`, invisibly.

## See also

[`boin_et_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_et_boundary.md)

## Examples

``` r
print(boin_et_boundary(target = 0.30, target_resp = 0.60))
#> BOIN-ET boundaries
#>   target DLT rate      : 0.3
#>   target response rate : 0.6
#>   p_saf / p_tox        : 0.03 / 0.42
#>   p_saf_resp           : 0.36
#> 
#>   lambda1 (escalate)   : 0.14
#>   lambda2 (de-escalate): 0.35
#>   eta1 (response)      : 0.48
#>   probability of an incorrect decision: 0.04703
```
