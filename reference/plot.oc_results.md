# Plot the Operating Characteristics of Several Designs

Draw one per-dose quantity of the table built by
[`oc_results`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_results.md).
For a single agent, one panel per scenario shows bars for every dose,
colored by design and labeled with their values. For a combination, one
panel per design and scenario shows the value of every combination as a
tile.

## Usage

``` r
# S3 method for class 'oc_results'
plot(x, metric = "sel_percent", text_size = 1, ...)
```

## Arguments

- x:

  An object of class `oc_results`.

- metric:

  Character scalar. The per-dose quantity: `"sel_percent"` (the
  default), `"n_pts"`, `"n_tox"` or `"n_resp"`.

- text_size:

  Numeric scalar. Multiplier applied to every text size in the figure.
  Defaults to 1.

- ...:

  Further arguments, currently ignored.

## Value

A ggplot2 object.

## Details

Requires the ggplot2 package, which is only suggested by simFastBOIN.

## See also

[`oc_results`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_results.md),
[`print.oc_results`](https://gosukehommaex.github.io/simFastBOIN/reference/print.oc_results.md)

## Examples

``` r
boin <- sim_boin(target = 0.30, p_true = c(0.10, 0.25, 0.40, 0.55),
                 n_cohort = 10, cohort_size = 3, n_trials = 200, seed = 1)
res <- oc_results(BOIN = boin, "3+3" = oc_3p3(c(0.10, 0.25, 0.40, 0.55)))
plot(res)

plot(res, metric = "n_pts")
```
