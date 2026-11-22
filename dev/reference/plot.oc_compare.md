# Plot a Comparison of Designs

Draw the comparison built by
[`oc_compare`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_compare.md):
one panel per quantity, with a point for every design and scenario, the
points of a design joined by a line.

## Usage

``` r
# S3 method for class 'oc_compare'
plot(
  x,
  metrics = c("pct_correct", "pct_patients_correct", "pct_overdose_selected",
    "n_pts_overdose"),
  text_size = 1,
  ...
)
```

## Arguments

- x:

  An object of class `oc_compare`.

- metrics:

  Character vector. Columns of `x` to draw. Defaults to
  `c("pct_correct", "pct_patients_correct", "pct_overdose_selected", "n_pts_overdose")`.

- text_size:

  Numeric scalar. Multiplier applied to every text size in the figure.
  Defaults to 1.

- ...:

  Further arguments, currently ignored.

## Value

A ggplot2 object.

## Details

Requires the ggplot2 package, which is only suggested by simFastBOIN.
Missing values (for example the correct selection in a scenario without
a true OBD) are left out.

## See also

[`oc_compare`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_compare.md),
[`print.oc_compare`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/print.oc_compare.md)

## Examples

``` r
truth <- oc_scenarios(list(c(0.10, 0.25, 0.40, 0.55)), target = 0.30)
boin <- sim_boin(target = 0.30, p_true = c(0.10, 0.25, 0.40, 0.55),
                 n_cohort = 10, cohort_size = 3, n_trials = 200, seed = 1)
plot(oc_compare(BOIN = boin, "3+3" = oc_3p3(c(0.10, 0.25, 0.40, 0.55)),
                truth = truth))
```
