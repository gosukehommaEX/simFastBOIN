# Plot a BOIN Decision Table

Draw the decision table produced by
[`boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_decision_table.md)
as a grid of colored cells, one per pair of DLT and patient counts, with
the decision letter in each cell.

## Usage

``` r
# S3 method for class 'boin_decision_table'
plot(x, text_size = 1, colors = NULL, ...)
```

## Arguments

- x:

  An object of class `boin_decision_table`.

- text_size:

  Numeric scalar. Multiplier applied to every text size in the figure.
  Use a larger value when the figure is exported at high resolution, for
  example into a report. Defaults to 1.

- colors:

  Named character vector of four colors, with names `"E"`, `"S"`, `"D"`
  and `"DE"`. Defaults to a palette whose adjacent categories remain
  distinguishable under the common forms of color vision deficiency.

- ...:

  Further arguments, currently ignored.

## Value

A ggplot2 object, which can be printed, saved with
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
or extended with further layers.

## Details

Every cell carries its decision letter, so the figure does not rely on
color alone to convey the decision and remains readable in grayscale.

Requires the ggplot2 package, which is only suggested by simFastBOIN
rather than required.

## See also

[`boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_decision_table.md),
[`print.boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/print.boin_decision_table.md)

## Examples

``` r
decisions <- boin_decision_table(target = 0.30, max_n = 12)

plot(decisions)


# Larger text for a high resolution export
plot(decisions, text_size = 2)


# A palette of your own
plot(decisions, colors = c(E = "#4DAF4A", S = "#377EB8",
                           D = "#FF7F00", DE = "#E41A1C"))
```
