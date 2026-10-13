# Exact Operating Characteristics of the 3+3 Design

Compute the operating characteristics of the traditional 3+3
dose-finding design exactly. The probability of every path the design
can take is evaluated in closed form, so the result carries no Monte
Carlo error and no seed is involved. Use it as a comparator for the BOIN
results returned by
[`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md).

## Usage

``` r
oc_3p3(
  p_true,
  mtd_rule = c("previous", "expand"),
  start_dose = 1,
  overdose_cutoff = NULL
)
```

## Arguments

- p_true:

  Numeric vector. True DLT probability at each dose level, in increasing
  dose order.

- mtd_rule:

  Character scalar, either `"previous"` or `"expand"`. The two
  definitions of the MTD in common use. See Details. Defaults to
  `"previous"`.

- start_dose:

  Integer scalar. Dose level for the first cohort. Defaults to 1.

- overdose_cutoff:

  Numeric scalar or `NULL`. Doses whose true DLT probability exceeds
  this value count as overdoses in the `overdose` component of the
  result. Defaults to `NULL`, which uses one third.

## Value

An object of class `oc_3p3`, which is a list with components

- sel_percent:

  Percentage of trials selecting each dose as the MTD.

- percent_no_mtd:

  Percentage of trials ending without an MTD.

- n_pts_dose:

  Average number of patients treated at each dose.

- n_tox_dose:

  Average number of DLTs observed at each dose.

- total_n_pts:

  Average total number of patients per trial.

- total_n_tox:

  Average total number of DLTs per trial.

- overdose:

  List describing the design's relationship with doses above
  `overdose_cutoff`: `pct_patients` (the percentage of all patients
  treated at such a dose, that is the probability that a patient is
  dosed above the cutoff), `avg_n_patients`, `pct_trials_any` (the
  percentage of trials treating at least one patient there),
  `pct_trials_mtd_above` (the percentage of all trials whose selected
  MTD is above the cutoff) and `pct_mtd_above_when_selected` (the same
  among the trials that selected an MTD). The `cutoff` and the dose
  levels `doses` that exceed it are also returned.

together with `p_true`, `n_doses`, `mtd_rule`, `start_dose`, `method`
and the call.

## Details

The escalation rule is the traditional one. Three patients are treated
at the current dose. No DLT leads to the next dose up. One DLT out of
three leads to three more patients at the same dose, and then to the
next dose up if no further DLT is seen and to stopping otherwise. Two or
more DLTs out of three stop the trial at once. A trial that escalates
past the highest dose ends there.

The two MTD rules differ only in what happens once escalation has
stopped. `"previous"` names the dose below the one declared too toxic,
and the highest dose when escalation ran out of doses. `"expand"`
instead takes the MTD to be the highest dose at which at most one of six
patients had a DLT: starting from the highest dose still eligible, a
dose that has only three patients is expanded to six and assessed, and
the search moves down a level whenever the expanded dose fails. Neither
rule is universally the standard, so the comparator being matched should
be checked before a rule is chosen. The two agree on which doses are
visited during escalation, and so on the probability that a patient is
exposed to an overly toxic dose; they differ in the sample size and in
the selection percentages.

The 3+3 design has no fixed maximum sample size, so the within-trial
exposure percentage and the 60 percent and 80 percent thresholds
reported for BOIN in
[`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md)
have no counterpart here.

## See also

[`sim_3p3`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_3p3.md),
[`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md)

## Examples

``` r
oc_3p3(p_true = c(0.30, 0.48, 0.67))
#> 3+3 operating characteristics
#>   MTD rule        : previous (dose below the toxic one)
#>   obtained by     : exact enumeration
#>   start dose      : DL1
#>   max sample size : 18 (six per dose)
#> 
#>                    DL1  DL2  DL3 Total / No MTD
#> True DLT rate (%) 30.0 48.0 67.0               
#> MTD selected (%)  39.8  9.2  0.4           50.6
#> Patients treated   4.3  2.1  0.4            6.7
#> Patients with DLT  1.3  1.0  0.2            2.5
#> 
#> Doses above a true DLT rate of 0.3333: DL2, DL3
#>   Patients treated there              : 35.8%
#>   Trials treating any patient there   : 49.4%
#>   Trials selecting an MTD there       : 9.7%

# The other definition of the MTD
oc_3p3(p_true = c(0.30, 0.48, 0.67), mtd_rule = "expand")
#> 3+3 operating characteristics
#>   MTD rule        : expand (at most one DLT in six)
#>   obtained by     : exact enumeration
#>   start dose      : DL1
#>   max sample size : 18 (six per dose)
#> 
#>                    DL1  DL2  DL3 Total / No MTD
#> True DLT rate (%) 30.0 48.0 67.0               
#> MTD selected (%)  36.5  6.3  0.2           57.0
#> Patients treated   5.2  2.3  0.4            7.8
#> Patients with DLT  1.6  1.1  0.2            2.9
#> 
#> Doses above a true DLT rate of 0.3333: DL2, DL3
#>   Patients treated there              : 33.5%
#>   Trials treating any patient there   : 49.4%
#>   Trials selecting an MTD there       : 6.4%

# Exposure to doses above a stated DLT rate
oc_3p3(p_true = c(0.20, 0.40, 0.60), overdose_cutoff = 1 / 3)$overdose
#> $cutoff
#> [1] 0.3333333
#> 
#> $doses
#> [1] 2 3
#> 
#> $pct_patients
#> [1] 48.37808
#> 
#> $avg_n_patients
#> [1] 3.891095
#> 
#> $pct_trials_any
#> [1] 70.8608
#> 
#> $pct_trials_mtd_above
#> [1] 21.9181
#> 
#> $pct_mtd_above_when_selected
#> [1] 30.9312
#> 
```
