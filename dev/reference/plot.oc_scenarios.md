# Plot True Scenarios

Draw the scenarios tabulated by
[`oc_scenarios`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_scenarios.md).
For a single agent, one panel per scenario shows the DLT probability,
the response probability and the utility divided by 100 against the
dose, with the true MTDs labeled "MTD" on the DLT curve and the true
OBDs labeled "OBD" on the utility curve. For a combination, one panel
per scenario shows the DLT probability of every combination as a tile
carrying its value, the true MTDs outlined with a solid line and the
true OBDs with a dashed line.

## Usage

``` r
# S3 method for class 'oc_scenarios'
plot(x, text_size = 1, ...)
```

## Arguments

- x:

  An object of class `oc_scenarios`.

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

[`oc_scenarios`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_scenarios.md),
[`print.oc_scenarios`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/print.oc_scenarios.md)

## Examples

``` r
truth <- oc_scenarios(
  list(list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
            p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38))),
  target = 0.35, target_resp = 0.25
)
plot(truth)
```
