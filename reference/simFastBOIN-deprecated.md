# Deprecated Functions in simFastBOIN

These functions were renamed in version 2.0.0. They still work but issue
a deprecation warning and will be removed in a future release. Their
return values are those of the replacement functions, which differ from
the return values of version 1.3.2, so calls should be updated rather
than left in place.

## Usage

``` r
get_boin_boundary(...)

get_boin_decision(...)

get_boin_stopping_boundaries(...)

get_pts_and_tox(...)

isotonic_regression(...)

select_mtd(...)
```

## Arguments

- ...:

  Arguments passed to the replacement function. Note that the
  replacements use different argument names and a different argument
  order, so only calls that name every argument will carry over
  unchanged.

## Value

The value of the corresponding replacement function.

## Replacements

- `get_boin_boundary()`:

  use
  [`boin_lambda`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_lambda.md)
  for the interval boundaries, or
  [`boin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_boundary.md)
  for the integer decision boundaries.

- `get_boin_decision()`:

  use
  [`boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_decision_table.md).

- `get_boin_stopping_boundaries()`:

  use
  [`boin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_boundary.md)
  with `extrasafe = TRUE` and read the `b_stop` component.

- `get_pts_and_tox()`:

  use
  [`boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_simulate.md).

- `isotonic_regression()`:

  use
  [`boin_isotonic`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_isotonic.md).

- `select_mtd()`:

  use
  [`boin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_select_mtd.md),
  which now takes the patient and DLT counts instead of precomputed
  isotonic estimates.
