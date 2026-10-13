# Simulate the 3+3 Design

Simulate the traditional 3+3 dose-finding design and summarize the same
quantities that
[`oc_3p3`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_3p3.md)
computes exactly. The two are meant to be used together:
[`oc_3p3()`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_3p3.md)
gives the answer without Monte Carlo error, and this function confirms
it from independent trials.

## Usage

``` r
sim_3p3(
  p_true,
  n_trials = 10000,
  mtd_rule = c("previous", "expand"),
  start_dose = 1,
  overdose_cutoff = NULL,
  seed = 123
)
```

## Arguments

- p_true:

  Numeric vector. True DLT probability at each dose level, in increasing
  dose order.

- n_trials:

  Integer scalar. Number of trials to simulate. Defaults to 10000.

- mtd_rule:

  Character scalar, either `"previous"` or `"expand"`. See
  [`oc_3p3`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_3p3.md).
  Defaults to `"previous"`.

- start_dose:

  Integer scalar. Dose level for the first cohort. Defaults to 1.

- overdose_cutoff:

  Numeric scalar or `NULL`. Doses whose true DLT probability exceeds
  this value count as overdoses in the `overdose` component of the
  result. Defaults to `NULL`, which uses one third.

- seed:

  Integer scalar or `NULL`. Random seed. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `oc_3p3`, with the components described in
[`oc_3p3`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_3p3.md).
The `method` component is the number of trials rather than `"exact"`.

## Details

The design is described under
[`oc_3p3`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_3p3.md).
Since the operating characteristics of the 3+3 design can be written in
closed form, this function is a check on that closed form rather than
the way to obtain the numbers. Simulation is the only option when a
quantity is wanted that the closed form does not provide, such as the
distribution of the sample size rather than its mean.

## See also

[`oc_3p3`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_3p3.md),
[`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md)

## Examples

``` r
exact <- oc_3p3(p_true = c(0.30, 0.48, 0.67))
simulated <- sim_3p3(p_true = c(0.30, 0.48, 0.67), n_trials = 2000, seed = 1)

round(exact$sel_percent, 2)
#> [1] 39.77  9.23  0.42
round(simulated$sel_percent, 2)
#> [1] 40.30  9.60  0.55
```
