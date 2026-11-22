# Optimal dose finding with the BOIN-ET family

## Why the optimal biological dose

For chemotherapy the most effective dose is usually taken to be the most
toxic one that is still acceptable, the maximum tolerated dose (MTD).
For targeted and immune therapies the response may level off, or even
fall, below the MTD, and the aim is then the optimal biological dose
(OBD), which weighs the response against the toxicity. The BOIN-ET
design of Takeda, Taguri and Morita (2018) decides the next dose from
both the DLT rate and the response rate at the current dose, and selects
the OBD at the end of the trial.

The package provides four designs of this family:

- **BOIN-ET** (Takeda et al., 2018), in which every cohort is assessed
  before the next one is dosed:
  [`sim_boin_et()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_et.md)
  and
  [`sim_boin_et_multi()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_et_multi.md);
- **TITE-BOIN-ET** (Takeda, Morita and Taguri, 2020), which doses new
  patients while earlier ones are still assessed:
  [`sim_tite_boin_et()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_tite_boin_et.md)
  and
  [`sim_tite_boin_et_multi()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_tite_boin_et_multi.md);
- **BF-BOIN-ET** (Takeda, Zhu and Hirakawa, 2025), which treats at lower
  doses the patients who arrive while a cohort is assessed:
  [`sim_bf_boin_et()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_bf_boin_et.md)
  and
  [`sim_bf_boin_et_multi()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_bf_boin_et_multi.md);
- **BOIN-ETC** (Kakizume, Takeda, Taguri and Morita, 2024), for two
  drugs given together, which runs BOIN-ET in subtrials of the dose
  matrix:
  [`sim_boin_etc()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_etc.md)
  and
  [`sim_boin_etc_multi()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_etc_multi.md).

The rules follow the R packages of the authors, ‘boinet’ (Yamaguchi,
Takeda and Maruo, 2025) for the first two and ‘bfboinet’ for the third,
with the changes described below, and the article for BOIN-ETC, with the
rules of the code of the authors as options. The trials and the
selection of the OBD run in C++.

## The boundaries and the decision table

The DLT rate at the current dose is compared with two boundaries,
`lambda1` and `lambda2`, and the response rate with one, `eta1`. They
minimize the probability of an incorrect decision, as found by the grid
search of
[`boinet::gridoptim()`](https://rdrr.io/pkg/boinet/man/gridoptim.html).

``` r

bd <- boin_et_boundary(target = 0.30, target_resp = 0.60)
bd
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

With a target DLT rate of 0.30, a target response rate of 0.60 and the
design parameters of Takeda et al. (2018), the boundaries are 0.14 for
`lambda1`, 0.35 for `lambda2` and 0.48 for `eta1`, the values of the
article. The decision at the current dose is then

- **E**, escalate, when the DLT rate is at most `lambda1` and the
  response rate at most `eta1`;
- **S**, stay, when the DLT rate is below `lambda2` and the response
  rate above `eta1`;
- **D**, de-escalate, when the DLT rate is at least `lambda2`;
- **E/S/D** otherwise, which moves to the dose with the highest observed
  response rate among the current dose and its neighbors, or escalates
  when the next higher dose has not been tried.

``` r

boin_et_decision_table(target = 0.30, target_resp = 0.60, n = c(3, 6))
#> BOIN-ET decision table
#>   lambda1 = 0.14, lambda2 = 0.35, eta1 = 0.48
#> 
#> 3 patients
#>     Responses
#> DLTs 0     1     2  3 
#>    0 E     E     S  S 
#>    1 E/S/D E/S/D S  S 
#>    2 D     D     D  D 
#>    3 D*    D*    D* D*
#> 
#> 6 patients
#>     Responses
#> DLTs 0     1     2     3  4  5  6 
#>    0 E     E     E     S  S  S  S 
#>    1 E/S/D E/S/D E/S/D S  S  S  S 
#>    2 E/S/D E/S/D E/S/D S  S  S  S 
#>    3 D     D     D     D  D  D  D 
#>    4 D*    D*    D*    D* D* D* D*
#>    5 D*    D*    D*    D* D* D* D*
#>    6 D*    D*    D*    D* D* D* D*
#> 
#> E = escalate, S = stay, D = de-escalate, E/S/D = choose among the current and neighboring doses;
#> * = the dose is no longer admissible
```

These are the decisions of Table 2 of Takeda et al. (2018). A star marks
the counts that make a dose inadmissible: by default, as in ‘boinet’, a
dose is excluded when the posterior probability that its DLT rate
exceeds `target` is above `cutoff_eli = 0.95`, or when the posterior
probability that its response rate is below `p_saf_resp` is above
`cutoff_eli_resp = 0.99`. `p_eli` and `p_eli_resp` change the two rates;
the simulations of the article exclude a dose when the posterior
probability that its DLT rate is below 0.40 is less than 0.10, or that
its response rate is above 0.20 is less than 0.10 (`p_eli = 0.40`,
`p_eli_resp = 0.20`, `cutoff_eli = 0.90` and `cutoff_eli_resp = 0.90`).

## Operating characteristics for one scenario

Scenario 1 of Takeda et al. (2018) has six doses, twelve cohorts of
three patients and a trial that starts at the lowest dose.

``` r

p_true <- c(0.05, 0.15, 0.25, 0.30, 0.35, 0.40)
p_resp <- c(0.05, 0.35, 0.55, 0.65, 0.75, 0.70)

oc <- sim_boin_et(target = 0.30, target_resp = 0.60, p_true = p_true,
                  p_resp = p_resp, n_cohort = 12, cohort_size = 3,
                  n_trials = 2000, seed = 123)
oc
#> BOIN-ET operating characteristics
#>   target DLT rate      : 0.3
#>   target response rate : 0.6
#>   trials               : 2000
#>   cohorts              : 12 of size 3
#>   max sample size      : 36
#>   lambda1 / lambda2    : 0.14 / 0.35
#>   eta1                 : 0.48
#>   OBD selection        : max_resp with observed response rates
#> 
#>                        DL1  DL2  DL3  DL4  DL5  DL6 Total / No OBD
#> True DLT rate (%)      5.0 15.0 25.0 30.0 35.0 40.0               
#> True response rate (%) 5.0 35.0 55.0 65.0 75.0 70.0               
#> OBD selected (%)       1.5 16.1 47.1 26.7  7.0  0.5            1.1
#> Patients treated       3.5  8.9 14.3  7.1  1.8  0.2           35.8
#> Patients with DLT      0.2  1.3  3.6  2.2  0.6  0.1            8.0
#> Patients with response 0.2  3.1  7.9  4.7  1.4  0.1           17.3
#> 
#> Trials stopped early: 1.1%
```

The output gives the percentage of trials selecting each dose as the
OBD, the average numbers of patients, DLTs and responses at each dose,
and the percentage of trials stopped early. Here dose 3 is selected most
often, in 47 percent of the trials, and 1.1 percent of the trials select
no OBD.

Several scenarios are run with
[`sim_boin_et_multi()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_et_multi.md).
In scenario 7 of the article the response rate peaks at dose 2 while the
DLT rate goes on rising.

``` r

scenarios <- list(
  list(name = "Scenario 1", p_true = p_true, p_resp = p_resp),
  list(name = "Scenario 7", p_true = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55),
       p_resp = c(0.30, 0.70, 0.60, 0.50, 0.40, 0.30))
)
sim_boin_et_multi(target = 0.30, target_resp = 0.60, scenarios = scenarios,
                  n_cohort = 12, cohort_size = 3, n_trials = 1000,
                  seed = 123)
#> BOIN-ET operating characteristics across 2 scenarios
#>   target DLT rate      : 0.3
#>   target response rate : 0.6
#>   trials each          : 1000
#>   max sample size      : 36
#>   OBD selection        : max_resp with observed response rates
#> 
#>    Scenario                   Item  DL1  DL2  DL3  DL4  DL5  DL6 Total / No OBD
#>  Scenario 1      True DLT rate (%)  5.0 15.0 25.0 30.0 35.0 40.0               
#>             True response rate (%)  5.0 35.0 55.0 65.0 75.0 70.0               
#>                   OBD selected (%)  1.7 17.0 47.4 26.5  6.1  0.3            1.0
#>                   Patients treated  3.5  9.1 14.4  7.0  1.7  0.2           35.9
#>                  Patients with DLT  0.2  1.4  3.6  2.2  0.6  0.1            8.0
#>             Patients with response  0.2  3.1  7.9  4.6  1.2  0.1           17.3
#>  Scenario 7      True DLT rate (%)  5.0 15.0 25.0 35.0 45.0 55.0               
#>             True response rate (%) 30.0 70.0 60.0 50.0 40.0 30.0               
#>                   OBD selected (%)  3.0 80.8 12.1  3.3  0.7  0.0            0.1
#>                   Patients treated  5.6 25.5  3.4  1.1  0.4  0.1           36.0
#>                  Patients with DLT  0.3  3.9  0.8  0.4  0.1  0.1            5.6
#>             Patients with response  1.7 17.8  2.1  0.6  0.1  0.0           22.3
```

## Selecting the OBD

At the end of a trial the DLT rates of the doses tried are made monotone
by isotonic regression, and the response rates are estimated with one of
three estimators (`est_method`): the observed rates (the default), the
model-averaged unimodal isotonic regression of Lin and Yin (2017)
(`"multi_iso"`), or the fractional polynomial logistic regression of
Takeda et al. (2018) (`"fp_logistic"`, which needs the ‘mfp’ package).
The OBD is then chosen among the admissible doses with one of four rules
(`obd_method`): the highest estimated response rate up to the estimated
MTD (`"max_resp"`, the default), or the highest value of one of the
three utilities compared by Yamaguchi et al. (2024)
(`"utility_weighted"`, `"utility_truncated_linear"` and
`"utility_scoring"`).

[`boin_et_select_obd()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_select_obd.md)
selects the OBD from the final data of trials, so the same simulated
trials can be judged under every rule.

``` r

oc_kept <- sim_boin_et(target = 0.30, target_resp = 0.60, p_true = p_true,
                       p_resp = p_resp, n_cohort = 12, cohort_size = 3,
                       n_trials = 2000, keep_trials = TRUE, seed = 123)
tr <- oc_kept$trials

select_percent <- function(est_method, obd_method) {
  obd <- boin_et_select_obd(n_pts = tr$n_pts, n_tox = tr$n_tox,
                            n_resp = tr$n_resp, target = 0.30,
                            target_resp = 0.60, est_method = est_method,
                            obd_method = obd_method,
                            early_stop = tr$early_stop)$obd
  as.numeric(table(factor(obd, levels = 1:6))) / length(obd) * 100
}

rules <- c("max_resp", "utility_weighted", "utility_truncated_linear",
           "utility_scoring")
selection <- rbind(
  t(sapply(rules, function(m) select_percent("observed", m))),
  "max_resp, multi_iso" = select_percent("multi_iso", "max_resp")
)
colnames(selection) <- paste0("DL", 1:6)
round(selection, 1)
#>                          DL1  DL2  DL3  DL4  DL5 DL6
#> max_resp                 1.5 16.1 47.1 26.7  7.0 0.5
#> utility_weighted         1.5 16.4 44.0 28.5  7.6 0.8
#> utility_truncated_linear 5.9 34.3 39.1 15.7  3.3 0.4
#> utility_scoring          2.1 16.5 39.4 29.8 10.1 1.0
#> max_resp, multi_iso      0.2 14.4 48.9 27.7  7.0 0.5

all.equal(unname(selection["max_resp", ]), unname(oc_kept$sel_percent))
#> [1] TRUE
```

The first row is the selection of
[`sim_boin_et()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_et.md),
which calls
[`boin_et_select_obd()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_select_obd.md)
with the same defaults. The fractional polynomial is fitted once for
every distinct set of final data and takes longer:

``` r

sim_boin_et(target = 0.30, target_resp = 0.60, p_true = p_true,
            p_resp = p_resp, n_cohort = 12, cohort_size = 3,
            n_trials = 2000, est_method = "fp_logistic", seed = 123)
```

## Dose moves that differ from ‘boinet’

After each cohort the admissibility of every dose is derived again from
its data, and the trial stops early when no dose is admissible, or when
the dose should be lowered and no admissible dose is found at or below
the intended one. When the intended dose is not admissible, a stay or a
de-escalation moves to the highest admissible dose below it, and an
escalation to the lowest admissible dose above it, or stays when there
is none.

Two rules differ from those of ‘boinet’ 1.6.0. Takeda et al. (2018) do
not allow dose skipping during escalation, so a dose excluded by its DLT
data also excludes every higher dose, as in the later BF-BOIN-ET design
and the ‘bfboinet’ package, and an escalation from the highest dose ends
the trial when that dose is excluded. In ‘boinet’, an escalation can
skip a dose excluded by its DLT data and reach an untried higher dose,
and an escalation from the highest dose stays there even when it is
excluded. `boinet_compatible = TRUE` applies these rules of ‘boinet’ in
every BOIN-ET and TITE-BOIN-ET function.

``` r

skip_design <- list(target = 0.30, target_resp = 0.60,
                    p_true = c(0.10, 0.20, 0.45, 0.50, 0.55),
                    p_resp = c(0.10, 0.20, 0.30, 0.50, 0.60),
                    n_cohort = 12, cohort_size = 3, n_trials = 2000,
                    seed = 123)
oc_default <- do.call(sim_boin_et, skip_design)
oc_boinet <- do.call(sim_boin_et,
                     c(skip_design, list(boinet_compatible = TRUE)))
round(rbind(default = oc_default$n_pts_dose,
            boinet = oc_boinet$n_pts_dose), 2)
#>          DL1   DL2  DL3  DL4  DL5
#> default 5.20 14.65 9.73 4.77 1.13
#> boinet  4.84 11.57 9.47 6.89 3.05
```

Doses 4 and 5, with DLT probabilities of 0.50 and 0.55, receive 5.9
patients per trial on average under the default rules and 9.9 under
those of ‘boinet’, which let an escalation from dose 2 skip dose 3 once
it is excluded.

## TITE-BOIN-ET

When the DLT and the response take long to assess compared with the time
between arrivals, TITE-BOIN-ET doses new patients while earlier ones are
still assessed. A patient still under assessment counts in the DLT and
response rates of a dose by the fraction of the window followed so far
(Takeda et al., 2020), and the decision rule of BOIN-ET is applied to
these rates. The admissibility of the doses counts such a patient as a
patient without the event, as in ‘boinet’.

The first patient of a cohort arrives at the time of the previous
decision, and the next decision is due when the next patient arrives.
While fewer than a share `min_completed_share` of the patients at the
current dose have finished both assessments, the decision waits for the
next assessment to end. The default, 0.5, is the rule of
[`boinet::tite.boinet()`](https://rdrr.io/pkg/boinet/man/tite.boinet.html).
The simulations of Takeda et al. (2020) did not wait
(`min_completed_share = 0`), and the article recommends such a rule when
the accrual is fast. By default the times to DLT and to response are
Weibull, with `late_fraction` and `resp_late_fraction` of the events in
the second half of each window, and they are linked by a Gaussian copula
with correlation `resp_cor`.

The setting of Takeda et al. (2020) is a DLT window of 28 days, a
response window of 56 days and one patient every 14 days; the unit of
time is the day.

``` r

tite_design <- list(target = 0.30, target_resp = 0.60, p_true = p_true,
                    p_resp = p_resp, n_cohort = 12, cohort_size = 3,
                    window = 28, resp_window = 56, accrual_rate = 1 / 14,
                    n_trials = 1000, seed = 123)
oc_tite <- do.call(sim_tite_boin_et, tite_design)
oc_tite
#> TITE-BOIN-ET operating characteristics
#>   target DLT rate      : 0.3
#>   target response rate : 0.6
#>   trials               : 1000
#>   cohorts              : 12 of size 3
#>   max sample size      : 36
#>   lambda1 / lambda2    : 0.14 / 0.35
#>   eta1                 : 0.48
#>   OBD selection        : max_resp with observed response rates
#> 
#>                        DL1  DL2  DL3  DL4  DL5  DL6 Total / No OBD
#> True DLT rate (%)      5.0 15.0 25.0 30.0 35.0 40.0               
#> True response rate (%) 5.0 35.0 55.0 65.0 75.0 70.0               
#> OBD selected (%)       1.1 18.0 44.8 26.9  6.1  0.4            2.7
#> Patients treated       4.0  9.3 13.5  6.9  1.8  0.2           35.6
#> Patients with DLT      0.2  1.5  3.4  2.1  0.7  0.1            7.9
#> Patients with response 0.2  3.2  7.4  4.5  1.4  0.2           16.9
#> 
#> Trials stopped early: 2.7%
#> Average trial duration: 616.5

oc_tite_0 <- do.call(sim_tite_boin_et,
                     c(tite_design, list(min_completed_share = 0)))
round(rbind(
  "wait for half" = c(oc_tite$sel_percent, weeks = oc_tite$duration / 7),
  "no wait" = c(oc_tite_0$sel_percent, weeks = oc_tite_0$duration / 7)
), 1)
#>               DL1  DL2  DL3  DL4 DL5 DL6 weeks
#> wait for half 1.1 18.0 44.8 26.9 6.1 0.4  88.1
#> no wait       0.9 20.7 45.6 23.0 7.2 0.3  77.1
```

The trials last 88 weeks on average with the default rule and 77 weeks
without waiting. The duration runs from the arrival of the first patient
to the last decision, which is the end of both assessments of every
patient when all cohorts are enrolled. The durations reported by Takeda
et al. (2020) correspond to counting also the interval between arrivals
before the first patient.

## BF-BOIN-ET

In BF-BOIN-ET a cohort is treated at the current dose, and the next
cohort starts with the first patient who arrives once the DLT outcomes
of the cohort are known. The patients who arrive in between are
backfilled one by one, when the current dose is not the lowest, to the
lower dose with the highest utility among those with a promising
response rate. The utility is defined by the scores `psi00`,
`psi_dlt_no_resp`, `psi_resp_no_dlt` and `psi11` of a patient without
DLT or response, with a DLT only, with a response only and with both.
The dose escalation uses the DLT rate at the current dose with the
boundaries of BOIN, and stays after 1 DLT in 3 patients. A dose with
`n_elimination` patients or more is excluded for good when its DLT or
response data are poor, a dose excluded by its DLT data taking every
higher dose with it. The OBD is the admissible dose up to the estimated
MTD with the highest utility. The rules are those of
[`bfboinet::get.oc.backboinet()`](https://rdrr.io/pkg/bfboinet/man/get.oc.backboinet.html).

The example of ‘bfboinet’ has five doses, ten cohorts of three, windows
of one month and three patients a month; the unit of time is the month.

``` r

bf_design <- list(target = 0.30, target_resp = 0.25,
                  p_true = c(0.03, 0.05, 0.20, 0.22, 0.45),
                  p_resp = c(0.05, 0.10, 0.50, 0.68, 0.70),
                  n_cohort = 10, cohort_size = 3, window = 1,
                  accrual_rate = 3, n_trials = 1000, seed = 123)
oc_bf <- do.call(sim_bf_boin_et, bf_design)
oc_bf
#> BF-BOIN-ET operating characteristics
#>   target DLT rate      : 0.3
#>   target response rate : 0.25
#>   trials               : 1000
#>   cohorts              : 10 of size 3
#>   patients per dose    : 12 at most with backfill
#>   lambda1 / lambda2    : 0.2365 / 0.3585
#>   eta1                 : 0.1968
#> 
#>                        DL1  DL2  DL3  DL4  DL5 Total / No OBD
#> True DLT rate (%)      3.0  5.0 20.0 22.0 45.0               
#> True response rate (%) 5.0 10.0 50.0 68.0 70.0               
#> OBD selected (%)       2.2  3.5 40.7 48.0  5.5            0.1
#> Patients treated       4.1  5.9 11.3  9.8  5.3           36.3
#> Patients backfilled    0.7  1.4  4.5  2.8  0.2            9.5
#> Patients with DLT      0.1  0.3  2.2  2.1  2.4            7.2
#> Patients with response 0.2  0.6  5.6  6.6  3.7           16.7
#> 
#> Trials stopped early: 0%
#> Average trial duration: 16.8
```

With `suspend = 0.33` the next cohort may start while a third of the
patients at the current dose are still assessed for DLT, which is the
BF-BOIN-ET-A design of Takeda et al. (2025); the default, `suspend = 0`,
is their BF-BOIN-ET-C design.

``` r

oc_bf_a <- do.call(sim_bf_boin_et, c(bf_design, list(suspend = 0.33)))
round(rbind(
  "BF-BOIN-ET-C" = c(oc_bf$sel_percent, patients = oc_bf$total_n_pts,
                     backfilled = oc_bf$total_n_bf,
                     months = oc_bf$duration),
  "BF-BOIN-ET-A" = c(oc_bf_a$sel_percent, patients = oc_bf_a$total_n_pts,
                     backfilled = oc_bf_a$total_n_bf,
                     months = oc_bf_a$duration)
), 1)
#>              DL1 DL2  DL3  DL4 DL5 patients backfilled months
#> BF-BOIN-ET-C 2.2 3.5 40.7 48.0 5.5     36.3        9.5   16.8
#> BF-BOIN-ET-A 3.2 3.5 40.7 45.6 6.2     33.6        6.1   13.3
```

The patients include the backfilled ones. The trials last 16.8 months on
average under BF-BOIN-ET-C and 13.3 months under BF-BOIN-ET-A.

## BOIN-ETC

When two drugs are given together, the doses form a matrix, with the `L`
levels of drug A in the rows and the `M` levels of drug B in the
columns. BOIN-ETC (Kakizume et al., 2024) divides the matrix into
subtrials along which the DLT rate increases and runs the rules of
BOIN-ET in each of them, from its lowest combination. A subtrial ends
when the combination chosen for the next cohort already has
`n_earlystop` patients (6 by default, `n_stop` of the article), and that
combination is the candidate of the subtrial. The first subtrial goes up
drug A at the lowest level of drug B and then along the highest level of
drug A, `A1B1, A2B1, ..., ALB1, ALB2, ..., ALBM`. The three versions
then differ (`subtrials`):

- **BOIN-ETC1** (`"etc1"`, the default) and **BOIN-ETC2** (`"etc2"`) run
  one subtrial per level of drug A, from level `L - 1` down to level 1,
  over levels 2 to `M` of drug B; when the candidate of the first
  subtrial is at level `l` of drug A and level 1 of drug B, they start
  at level `l - 1` of drug A. BOIN-ETC2 starts every subtrial at level 2
  of drug B. BOIN-ETC1 starts it next to the candidate of the previous
  subtrial, and the article and its code differ there: the text and the
  figures start one level of drug B lower, diagonally below the
  candidate (`etc1_start = "diagonal"`, the default), and the code of
  the authors at the same level of drug B, directly below it
  (`etc1_start = "below"`), in both cases at level 2 or above.
- **BOIN-ETC3** (`"etc3"`) runs L-shaped subtrials: the subtrial `j`
  goes up level `j + 1` of drug B to level `L - j` of drug A and then
  along that level of drug A to level `M` of drug B.

The names follow the figures, the tables and the code of the article;
the text of its Section 2.2 describes BOIN-ETC2 and BOIN-ETC3 the other
way round. A combination is excluded when the posterior probability that
its DLT rate exceeds `target` is above `cutoff_eli = 0.95`, and with it
the higher combinations of its subtrial; the subtrial ends when its
lowest combination is excluded, and the trial ends when that combination
is A1B1. In the code of the authors, the choice among the neighboring
combinations when the decision of BOIN-ET is E/S/D moves one combination
down when the combination below the current one is untreated;
`rule4 = "authors_code"` applies this rule, and `rule4 = "paper"`, the
default, the rule of the text. Only `etc1_start = "below"` with
`rule4 = "authors_code"` reproduces the allocation published for
BOIN-ETC1.

At the end of the trial the DLT rates are made monotone along both drugs
by bivariate isotonic regression, the MTD combination of every level of
drug A bounds the candidates of that level, and the optimal dose
combination (ODC) is the treated, not excluded candidate with the
highest posterior mean of the utility of BOIN12, with the scores `psi00`
and `psi11` (40 and 60 by default).
[`boin_etc_select_obd()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_etc_select_obd.md)
applies this selection to completed trials. When candidates tie,
`odc_tie = "lowest"` (the default) selects the lowest one, and
`odc_tie = "all"` counts every tied candidate, as the percentages of the
article do.

Scenario 1 of the article has two levels of drug A, three of drug B, a
target DLT rate of 0.30, a target response rate of 0.60 and at most 24
patients.

``` r

p_true_ab <- matrix(c(0.05, 0.10, 0.10, 0.30, 0.30, 0.50), nrow = 2)
p_resp_ab <- matrix(c(0.10, 0.25, 0.25, 0.60, 0.60, 0.65), nrow = 2)

oc_etc <- sim_boin_etc(target = 0.30, target_resp = 0.60, p_true = p_true_ab,
                       p_resp = p_resp_ab, n_cohort = 8, cohort_size = 3,
                       n_trials = 2000, seed = 123)
oc_etc
#> BOIN-ETC1 operating characteristics
#>   target DLT rate      : 0.3
#>   target response rate : 0.6
#>   dose levels          : 2 of drug A, 3 of drug B
#>   trials               : 2000
#>   cohorts              : 8 of size 3
#>   max sample size      : 24
#>   subtrial end (n_stop): 6
#>   lambda1/lambda2/eta1 : 0.14 / 0.35 / 0.48
#>   start of a row       : diagonal
#>   rule 4               : paper
#>   utility scores       : 100, 40, 60, 0; ties lowest
#> 
#>                    Item Drug A   B1   B2   B3
#>       True DLT rate (%)     A1  5.0 10.0 30.0
#>                             A2 10.0 30.0 50.0
#>  True response rate (%)     A1 10.0 25.0 60.0
#>                             A2 25.0 60.0 65.0
#>        OBD selected (%)     A1  3.4 11.3 31.3
#>                             A2 13.7 36.4  3.9
#>        Patients treated     A1  3.2  4.2  4.5
#>                             A2  4.3  4.9  1.3
#>       Patients with DLT     A1  0.2  0.4  1.3
#>                             A2  0.4  1.5  0.6
#>  Patients with response     A1  0.3  1.0  2.7
#>                             A2  1.1  3.0  0.8
#> 
#> No OBD selected                      : 0%
#> Trials stopped early                 : 0%
#> Average number of patients           : 22.3
#> Average number of subtrials          : 2
#> 
#> True OBD (highest true utility)      : A2B2, A1B3
#>   Selected as the OBD                : 67.8%
#>   Patients treated there             : 42.1%
#> Combinations above the target DLT    : 1 of 6
#>   Selected as the OBD                : 3.9%
#>   Patients treated there             : 5.7%
```

The true ODCs, the combinations with the highest utility among those
with a DLT rate at most 0.30, are A2B2, A1B3. The three versions, and
BOIN-ETC1 with the rules of the code of the authors, compare as follows,
with tied ODCs counted as in the article:

``` r

etc_versions <- list(
  "BOIN-ETC1" = list(subtrials = "etc1"),
  "BOIN-ETC1, code" = list(subtrials = "etc1", etc1_start = "below",
                           rule4 = "authors_code"),
  "BOIN-ETC2" = list(subtrials = "etc2"),
  "BOIN-ETC3" = list(subtrials = "etc3")
)
etc_runs <- lapply(etc_versions, function(v) {
  do.call(sim_boin_etc,
          c(list(target = 0.30, target_resp = 0.60, p_true = p_true_ab,
                 p_resp = p_resp_ab, n_cohort = 8, cohort_size = 3,
                 n_trials = 2000, odc_tie = "all", seed = 123), v))
})
round(t(sapply(etc_runs, function(o) {
  c("ODC selected (%)" = o$pct_true_obd_selected,
    "Patients at ODC (%)" = o$pct_patients_true_obd,
    "Patients" = o$total_n_pts, "Subtrials" = o$mean_subtrials)
})), 1)
#>                 ODC selected (%) Patients at ODC (%) Patients Subtrials
#> BOIN-ETC1                   73.4                42.1     22.3         2
#> BOIN-ETC1, code             75.1                42.7     22.3         2
#> BOIN-ETC2                   73.4                42.1     22.3         2
#> BOIN-ETC3                   73.4                42.1     22.3         2
```

## How the implementation is checked

The scripts in `inst/reproduce/` and `inst/validation/` compare the
package with the articles and with the packages and the code of the
authors, and write the tables of judgments in their `results/` folders,
which are summarized here. A simulated cell passes when it lies within
four standard errors of the difference between the two simulations plus
the rounding of the published value.

- **Takeda et al. (2018).** 68 of the 68 entries of the boundaries and
  of the decisions of Table 2 agree exactly. Of the 156 cells of Table
  3, 142 pass; 10 of the 14 others are in scenario 11, where the printed
  response rate of 0.55 at the lowest dose does not fit the published
  results, and ‘boinet’ too selects that dose far more often than the
  article reports.
- **Takeda et al. (2020).** Of the 324 cells of Table 2, 295 pass. The
  BOIN-ET rows fail in 12 cells, in scenario 11 and in the cells that
  also fail for Table 3 of the 2018 article. Of the 17 cells of
  TITE-BOIN-ET that fail, 8 are in scenario 11 and 6 are early
  terminations, durations and allocations in scenarios 8, 9, 10 and 12,
  which agree better with an exclusion of doses that uses only the
  patients whose outcomes are known; the article does not say which
  patients this rule uses, and the package keeps the rule of ‘boinet’.
  The remaining 3 are not explained. The durations are compared after
  adding the mean interval between arrivals before the first patient.
- **Takeda et al. (2025).** 336 of the 336 cells of Tables 1 and 2 pass,
  with and without the suspension rule.
- **‘boinet’.** The operating characteristics of
  [`sim_boin_et()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_et.md)
  and
  [`sim_tite_boin_et()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_tite_boin_et.md)
  with `boinet_compatible = TRUE` are compared with those of
  [`boinet::boinet()`](https://rdrr.io/pkg/boinet/man/boinet.html) and
  [`boinet::tite.boinet()`](https://rdrr.io/pkg/boinet/man/tite.boinet.html):
  312 of 312 cells pass for BOIN-ET and 196 of 196 for TITE-BOIN-ET. The
  packages draw their random numbers differently, so they agree in
  distribution, not trial by trial. The tests also compare the
  boundaries with
  [`boinet::gridoptim()`](https://rdrr.io/pkg/boinet/man/gridoptim.html)
  and the OBD of simulated trials with
  [`boinet::obd.select()`](https://rdrr.io/pkg/boinet/man/obd.select.html)
  exactly.
- **‘bfboinet’.** 85 of the 85 cells comparing
  [`sim_bf_boin_et()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_bf_boin_et.md)
  with
  [`bfboinet::get.oc.backboinet()`](https://rdrr.io/pkg/bfboinet/man/get.oc.backboinet.html)
  pass.
- **Kakizume et al. (2024).** With `etc1_start = "below"` and
  `rule4 = "authors_code"` for BOIN-ETC1, 1168 of the 1176 values of the
  allocation (the mean numbers of patients at every combination of
  Supplementary Tables S4 and S5, and %nODC and \#N of Table 3 and
  Supplementary Tables S1 and S6 to S9) pass, and the 8 others are
  misprints that disagree with the other tables of the article. The 996
  values of the selection (%ODC and Supplementary Tables S2 and S3) are
  listed for information, because the code of the authors estimates the
  DLT rates with
  [`Iso::biviso()`](https://rdrr.io/pkg/Iso/man/biviso.html), whose
  convergence errors decide between the combinations that tie for the
  MTD combination of a level of drug A in the exact isotonic regression
  of the package; 982 of them lie within the tolerance.
- **The code of the authors.** `inst/validation/compare-with-BOINETC.R`
  runs the R code of the authors. With the random numbers of R after
  `set.seed(1234)` it gives 1657 of the 1692 published values of the
  allocation and the selection of Scenarios 1 to 5 and 7 to 9 to the
  printed digits, and the 35 others are misprints. On the same trials,
  [`boin_etc_select_obd()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_etc_select_obd.md)
  selects other combinations than the code in 4.6 percent of the trials
  on average, because combinations tie for the MTD combination of a
  level of drug A in the exact isotonic regression or because the code
  selects an untreated or excluded combination; no trial differs for
  another reason.
- **Separate implementations.** The trials of the four engines agree
  exactly, trial by trial, with implementations written in Python, from
  the R code of ‘boinet’ and ‘bfboinet’ and, for BOIN-ETC, from the
  rules described above, and fed with the same random numbers; so does
  the selection of BOIN-ETC.

## References

Takeda, K., Taguri, M. and Morita, S. (2018). BOIN-ET: Bayesian optimal
interval design for dose finding based on both efficacy and toxicity
outcomes. *Pharmaceutical Statistics*, 17(4), 383-395.

Takeda, K., Morita, S. and Taguri, M. (2020). TITE-BOIN-ET:
Time-to-event Bayesian optimal interval design to accelerate
dose-finding based on both efficacy and toxicity outcomes.
*Pharmaceutical Statistics*, 19(3), 335-349.

Takeda, K., Zhu, J. and Hirakawa, A. (2025). BF-BOIN-ET: A backfill
Bayesian optimal interval design using efficacy and toxicity outcomes
for dose optimization. *Pharmaceutical Statistics*, 24(2), e2470.

Yamaguchi, Y., Takeda, K. and Maruo, K. (2025). boinet: An R package for
Bayesian optimal interval design for dose finding based on both efficacy
and toxicity outcomes. *The R Journal*, 17(2), 15-34.

Yamaguchi, Y., Takeda, K., Yoshida, S. and Maruo, K. (2024). Optimal
biological dose selection in dose-finding trials with model-assisted
designs based on efficacy and toxicity: a simulation study. *Journal of
Biopharmaceutical Statistics*, 34(3), 379-393.

Kakizume, T., Takeda, K., Taguri, M. and Morita, S. (2024). BOIN-ETC: A
Bayesian optimal interval design considering efficacy and toxicity to
identify the optimal dose combinations. *Statistical Methods in Medical
Research*, 33(4), 716-727.

Lin, R. and Yin, G. (2017). STEIN: A simple toxicity and efficacy
interval design for seamless phase I/II clinical trials. *Statistics in
Medicine*, 36(26), 4106-4120.
