# Backfilling with the BF-BOIN and BE-BOIN designs

## Why backfill

In the BOIN design the next cohort is dosed only once the current one
has completed the DLT assessment, and patients who become eligible in
the meantime are usually turned away. Backfilling treats them at a lower
dose that has already been cleared, so that more patients contribute
data on safety and activity at doses that may later be carried forward,
while the dose escalation goes on as before.

Two designs are provided. BF-BOIN (Zhao et al., 2024) keeps the BOIN
decision table and fully staggered cohorts. BE-BOIN (Chen et al., 2026)
combines backfilling with the time-to-event BOIN design, so that the
escalation itself decides while some patients are still pending. In
both, the package selects the MTD from the data of all patients,
backfilled ones included. The second stage of Chen et al. (2026), which
randomizes patients between doses to choose the optimal biological dose,
is not part of the package.

## The rules of BF-BOIN

The trial has two parts that run side by side.

- **Dose escalation.** Cohorts are fully staggered. When every patient
  of the current cohort has completed the assessment, the next dose
  follows the BOIN decision table applied to the completed patients at
  the current dose, backfilled ones included. Elimination is checked at
  the same moments, at every dose.
- **Backfilling.** A patient who arrives while the escalation waits is
  treated at the highest dose below the current one that has not been
  eliminated, has fewer than `n_cap` patients (12 by default) and has a
  response observed at or below it. A dose is closed when its own DLT
  rate exceeds the de-escalation boundary and so does the rate pooled
  with the next higher dose; the doses above it close with it.
- **Conflicts.** When the data of a backfilled dose call for a more
  cautious decision than those of the current dose (Table 2 of Zhao et
  al., 2024), the decision is taken from the data pooled from that dose
  up to the current one.

The trial ends when the escalation has treated `n_cohort` cohorts, when
`n_earlystop` patients, backfilled ones included, have been treated at
the current dose and the decision is to stay, or when the lowest dose is
eliminated.

## Operating characteristics of BF-BOIN

Scenario 3 of Zhao et al. (2024) has the MTD at dose 3 and responses
that rise with the dose. The design is that of the article: a target of
0.25, ten cohorts of three, `n_earlystop = 9`, one DLT out of three
treated as a stay, three patients a month and an assessment window of
one month.

``` r

zhao_design <- list(target = 0.25, n_cohort = 10, cohort_size = 3,
                    window = 1, accrual_rate = 3, n_earlystop = 9,
                    stay_on_1_of_3 = TRUE, accrual = "uniform",
                    no_slot = "leave", n_trials = 2000, seed = 123)
p_true <- c(0.04, 0.12, 0.25, 0.43, 0.63)
p_resp <- c(0.10, 0.20, 0.30, 0.45, 0.58)

oc_bf <- do.call(sim_bf_boin,
                 c(zhao_design, list(p_true = p_true, p_resp = p_resp)))
oc_bf
#> BF-BOIN operating characteristics (Zhao et al., 2024)
#>   target DLT rate : 0.25
#>   trials          : 2000
#>   cohorts         : 10 of size 3
#>   escalation size : 30
#>   backfill cap    : 12 patients per dose
#>   lambda_e / _d   : 0.1968 / 0.2984
#>   window          : 1 (DLT), 1 (response)
#>   accrual         : 3 patients per unit of time (uniform)
#> 
#>                         DL1  DL2  DL3  DL4  DL5 Total / No MTD
#> True DLT rate (%)       4.0 12.0 25.0 43.0 63.0               
#> True response rate (%) 10.0 20.0 30.0 45.0 58.0               
#> MTD selected (%)        2.3 30.3 57.8  9.6  0.0            0.0
#> Patients treated        6.0 10.8  9.8  3.7  0.5           30.8
#>   of which backfilled   1.8  3.4  1.9  0.3  0.0            7.4
#> Patients with DLT       0.2  1.3  2.4  1.6  0.3            5.9
#> 
#> Doses above a true DLT rate of 0.25: DL4, DL5
#>   Patients treated there              : 13.6%
#>   Trials treating any patient there   : 49.5%
#>   Trials selecting an MTD there       : 9.6%
#> 
#> Backfilling and timing
#>   Trials with backfilled patients     : 85.2%
#>   Average trial duration              : 14.2
#>   Trials in which the escalation waited: 100%
#>   Average time the escalation waited  : 6.2
#>   Average patients turned away        : 11.2
```

Without any response no dose is opened for backfilling and the design is
the BOIN design. The table below sets the two runs next to the values
that Zhao et al. (2024) report in their Table 4: the percentage of
trials selecting the MTD and a dose above it, the average number of
patients below, at and above the MTD and in all, and the average
duration in months.

``` r

oc_boin <- do.call(sim_bf_boin,
                   c(zhao_design, list(p_true = p_true, p_resp = rep(0, 5))))

zhao_metrics <- function(oc) {
  c(oc$sel_percent[[3]], sum(oc$sel_percent[4:5]),
    sum(oc$n_pts_dose[1:2]), oc$n_pts_dose[[3]], sum(oc$n_pts_dose[4:5]),
    oc$total_n_pts, oc$duration_mean)
}
comparison <- rbind(
  "BF-BOIN, Table 4" = c(57.6, 9.7, 16.9, 9.8, 4.2, 30.9, 14.3),
  "BF-BOIN, simulated" = zhao_metrics(oc_bf),
  "BOIN, Table 4" = c(56.7, 12.3, 11.4, 8.0, 4.6, 23.9, 14.3),
  "BOIN, simulated" = zhao_metrics(oc_boin)
)
colnames(comparison) <- c("MTD (%)", "Above (%)", "n below", "n at",
                          "n above", "n total", "Months")
round(comparison, 1)
#>                    MTD (%) Above (%) n below n at n above n total Months
#> BF-BOIN, Table 4      57.6       9.7    16.9  9.8     4.2    30.9   14.3
#> BF-BOIN, simulated    57.8       9.6    16.8  9.8     4.2    30.8   14.2
#> BOIN, Table 4         56.7      12.3    11.4  8.0     4.6    23.9   14.3
#> BOIN, simulated       57.1      11.7    11.5  8.1     4.6    24.2   14.5
```

With 2,000 trials the simulated values carry a simulation error of about
one point for the percentages. The tests of the package compare all
twelve scenarios of Table 4 with 10,000 trials each; over 20 seeds the
largest gaps were 2.1 points for the selection, 0.4 patients and 0.17
months. The script `inst/validation/compare-with-Zhao2024.R` repeats the
comparison.

The response profile decides which doses are backfilled. Scenarios 2 and
6 of the article share the DLT probabilities but not the response
probabilities.

``` r

scenarios <- list(
  list(name = "Scenario 2", p_true = c(0.12, 0.25, 0.42, 0.49, 0.55),
       p_resp = c(0.20, 0.30, 0.40, 0.50, 0.60)),
  list(name = "Scenario 6", p_true = c(0.12, 0.25, 0.42, 0.49, 0.55),
       p_resp = c(0.30, 0.35, 0.36, 0.36, 0.36))
)
zhao_multi <- zhao_design[names(zhao_design) != "n_trials"]
do.call(sim_bf_boin_multi,
        c(zhao_multi, list(scenarios = scenarios, n_trials = 1000)))
#> BF-BOIN operating characteristics across 2 scenarios
#>   target DLT rate : 0.25
#>   trials each     : 1000
#>   escalation size : 30
#>   backfill cap    : 12 patients per dose
#>   window          : 1 (DLT), 1 (response)
#>   accrual         : 3 patients per unit of time (uniform)
#> 
#>    Scenario                   Item  DL1  DL2  DL3  DL4  DL5 Total / No MTD
#>  Scenario 2      True DLT rate (%) 12.0 25.0 42.0 49.0 55.0               
#>             True response rate (%) 20.0 30.0 40.0 50.0 60.0               
#>                   MTD selected (%) 32.4 57.5  9.2  0.6  0.0            0.3
#>                   Patients treated 10.4 10.0  3.9  0.6  0.0           24.8
#>                Backfilled patients  3.0  1.8  0.3  0.0  0.0            5.2
#>                  Patients with DLT  1.2  2.5  1.6  0.3  0.0            5.7
#>              Trial duration (mean)                                    11.8
#>  Scenario 6      True DLT rate (%) 12.0 25.0 42.0 49.0 55.0               
#>             True response rate (%) 30.0 35.0 36.0 36.0 36.0               
#>                   MTD selected (%) 31.5 57.3 10.4  0.5  0.0            0.3
#>                   Patients treated 11.3 10.2  3.7  0.6  0.0           25.8
#>                Backfilled patients  3.6  1.8  0.4  0.0  0.0            5.9
#>                  Patients with DLT  1.4  2.6  1.5  0.3  0.0            5.8
#>              Trial duration (mean)                                    12.0
```

## Points the articles leave open

The articles do not settle every detail of a simulated trial. The
package makes the following choices, and offers an argument where more
than one reading is reasonable.

- **A patient with nowhere to go.** When the escalation waits and no
  dose is open for backfilling, `no_slot = "wait"` (the default) lets
  the patient wait until one can be treated, as in
  [`sim_tite_boin()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin.md),
  and `no_slot = "leave"` turns the patient away. Zhao et al. (2024)
  state a Poisson process of arrivals, but their durations, for BOIN as
  well as for BF-BOIN, agree with uniform times between arrivals
  (`accrual = "uniform"`) and patients turned away, which is why the
  comparison above uses these settings.
- **The end of backfilling.** Backfilling stops when the last patient of
  the escalation is enrolled, and every enrolled patient is then
  followed to the end of the assessment.
- **Pending backfill patients.** The escalation of BF-BOIN uses only the
  patients who have completed the assessment.
- **Responses.** Neither article states how the time to response is
  distributed or how a response is related to a DLT. The time to
  response is Weibull, with `resp_late_fraction` of the responses in the
  second half of `resp_window`, and the DLT and the response are linked
  by a Gaussian copula with correlation `resp_cor`, as in Takeda et
  al. (2025), who set it to 0.2. The default here is 0, independence.
- **Which conflicting dose anchors the pooled estimate.** With several
  conflicting backfilled doses, Zhao et al. (2024) pool from the highest
  of them and Chen et al. (2026) from the lowest. `conflict_dose`
  chooses, with the default of each article.

## BE-BOIN

BE-BOIN decides the dose for a cohort when its first patient arrives,
with the time-to-event BOIN rule, so the escalation waits only when the
rule asks it to. Two rules of Chen et al. (2025) decide when that is:
rule 1 suspends accrual when fewer than 51 percent of the patients at
the current dose have completed the assessment
(`max_pending_ratio = 0.49`), and rule 2 allows escalation only when
every pending patient at the current dose has been followed for at least
a quarter of the window (`min_follow_up = 0.25`). The decision table
with both rules is Table A1 of the supplementary materials of Chen et
al. (2026), which
[`tite_boin_decision_table()`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_decision_table.md)
reproduces entry by entry; an escalation that needs rule 2 is marked
`MF >= 0.25`.

``` r

be_table <- tite_boin_decision_table(target = 0.25, max_n = 9,
                                     max_pending_ratio = 0.49,
                                     min_follow_up = 0.25)
print(be_table, cohort_size = 3)
#> TITE-BOIN decision table with target 0.25
#> Method: single mean imputation (Yuan et al., 2018); boundaries refer to STFT
#> 
#>  Patients DLTs Pending Escalate             Stay    De-escalate
#>  3        0    0       Y                                       
#>  3        0    1       MF >= 0.25                              
#>  3        0    >= 2    Suspend              Suspend Suspend    
#>  3        1, 2 <= 2                                 Y          
#>  3        3    0                                    Y & Elim   
#>  6        0    0       Y                                       
#>  6        0    1-2     MF >= 0.25                              
#>  6        0    >= 3    Suspend              Suspend Suspend    
#>  6        1    0       Y                                       
#>  6        1    1       >= 0.22 & MF >= 0.25 < 0.22             
#>  6        1    2       >= 1.38 & MF >= 0.25 < 1.38             
#>  6        1    >= 3    Suspend              Suspend Suspend    
#>  6        2, 3 <= 4                                 Y          
#>  6        >= 4 <= 2                                 Y & Elim   
#>  9        0    0       Y                                       
#>  9        0    1-4     MF >= 0.25                              
#>  9        0    >= 5    Suspend              Suspend Suspend    
#>  9        1    0       Y                                       
#>  9        1    1-3     MF >= 0.25                              
#>  9        1    4       >= 0.66 & MF >= 0.25 < 0.66             
#>  9        1    >= 5    Suspend              Suspend Suspend    
#>  9        2    <= 4                         Y                  
#>  9        2    >= 5    Suspend              Suspend Suspend    
#>  9        3, 4 <= 6                                 Y          
#>  9        >= 5 <= 4                                 Y & Elim   
#> 
#> Patients = treated at the current dose, DLTs = DLTs observed, Pending = assessment not yet completed
#> STFT = total follow-up time of the pending patients divided by the length of the assessment window
#> Accrual is suspended when more than 49% of the patients are pending
#> MF = shortest follow-up time of the pending patients divided by the length of the assessment window
#> While patients are pending, escalation requires MF >= 0.25; otherwise accrual is suspended
#> Y & Elim = de-escalate and eliminate this dose and all higher doses
```

Every estimate of BE-BOIN, for the escalation, for opening and closing
backfill doses and for conflicts, imputes the pending patients with the
single mean imputation of TITE-BOIN. Patients who arrive while the
escalation is suspended are backfilled with the rules of BF-BOIN. With a
window of three months and two patients a month, the two designs can be
compared on the same scenario; Chen et al. (2026) report shorter trials
for BE-BOIN.

``` r

slow_design <- list(target = 0.25, p_true = p_true, p_resp = p_resp,
                    n_cohort = 10, cohort_size = 3, window = 3,
                    accrual_rate = 2, n_earlystop = 9, n_trials = 2000,
                    seed = 123)
oc_be <- do.call(sim_be_boin, slow_design)
oc_be
#> BE-BOIN operating characteristics (Chen et al., 2026)
#>   target DLT rate : 0.25
#>   trials          : 2000
#>   cohorts         : 10 of size 3
#>   escalation size : 30
#>   backfill cap    : 12 patients per dose
#>   lambda_e / _d   : 0.1968 / 0.2984
#>   window          : 3 (DLT), 3 (response)
#>   accrual         : 2 patients per unit of time (exponential)
#> 
#>                         DL1  DL2  DL3  DL4  DL5 Total / No MTD
#> True DLT rate (%)       4.0 12.0 25.0 43.0 63.0               
#> True response rate (%) 10.0 20.0 30.0 45.0 58.0               
#> MTD selected (%)        3.5 31.6 54.1 10.4  0.4            0.0
#> Patients treated        7.0 10.6  8.5  3.2  0.4           29.8
#>   of which backfilled   1.9  2.8  1.2  0.2  0.0            6.1
#> Patients with DLT       0.3  1.3  2.1  1.4  0.2            5.4
#> 
#> Doses above a true DLT rate of 0.25: DL4, DL5
#>   Patients treated there              : 12.2%
#>   Trials treating any patient there   : 49.3%
#>   Trials selecting an MTD there       : 10.8%
#> 
#> Backfilling and timing
#>   Trials with backfilled patients     : 70.7%
#>   Average trial duration              : 22.3
#>   Trials in which the escalation waited: 99.7%
#>   Average time the escalation waited  : 8.1

oc_bf_slow <- do.call(sim_bf_boin, slow_design)
round(rbind(
  "BF-BOIN" = c(MTD = oc_bf_slow$sel_percent[[3]],
                patients = oc_bf_slow$total_n_pts,
                backfilled = oc_bf_slow$total_n_bf,
                months = oc_bf_slow$duration_mean),
  "BE-BOIN" = c(MTD = oc_be$sel_percent[[3]],
                patients = oc_be$total_n_pts,
                backfilled = oc_be$total_n_bf,
                months = oc_be$duration_mean)
), 1)
#>          MTD patients backfilled months
#> BF-BOIN 54.1     33.2        9.7   31.4
#> BE-BOIN 54.1     29.8        6.1   22.3
```

## Agreement with the BOIN and TITE-BOIN designs

The DLT variates of the escalation are drawn from R’s stream exactly as
in
[`sim_boin()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md),
and arrivals, backfill patients and responses have generators of their
own. When no dose is ever opened for backfilling, the trials are
therefore those of
[`sim_boin()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md)
with the same seed, and those of
[`sim_tite_boin()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin.md)
for BE-BOIN with `no_slot = "wait"`.

``` r

no_backfill <- sim_bf_boin(target = 0.25, p_true = p_true,
                           p_resp = rep(0, 5), n_cohort = 10,
                           cohort_size = 3, window = 1, accrual_rate = 3,
                           n_trials = 1000, seed = 7)
plain_boin <- sim_boin(target = 0.25, p_true = p_true, n_cohort = 10,
                       cohort_size = 3, n_trials = 1000, seed = 7)

all.equal(no_backfill$sel_percent, plain_boin$sel_percent)
#> [1] TRUE
all.equal(no_backfill$n_pts_dose, plain_boin$n_pts_dose)
#> [1] TRUE
```

## How the implementation is checked

- The operating characteristics of Table 4 of Zhao et al. (2024), twelve
  scenarios of BF-BOIN and five of BOIN, are reproduced within
  simulation error, durations included.
- The trial example of Figure 1 of Zhao et al. (2024) and that of
  Section S2 of the supplementary materials of Chen et al. (2026) are
  replayed patient by patient, with the random numbers scripted to give
  the outcomes of the articles.
- The decision tables of Chen et al. (2025) and of the supplementary
  materials of Chen et al. (2026) are reproduced entry by entry.
- The simulated trials agree exactly, trial by trial, with a separate
  implementation written in Python and run on a replica of R’s random
  number stream, over 2,700 trials of both designs and several settings.
- Without backfilling the trials are identical to those of the BOIN and
  TITE-BOIN engines, as shown above.

Chen et al. (2026) give their scenarios only as figures, so their
operating characteristics cannot be reproduced number for number.

## References

Zhao, Y., Yuan, Y., Korn, E. L. and Freidlin, B. (2024). Backfilling
Patients in Phase I Dose-Escalation Trials Using Bayesian Optimal
Interval Design (BOIN). *Clinical Cancer Research*, 30(4), 673-679.

Chen, K., Zhao, Y., Takeda, K. and Yuan, Y. (2026). BE-BOIN: A Dose
Optimization Design Accommodating Backfill and Late-Onset Toxicity.
*Therapeutic Innovation and Regulatory Science*.

Chen, K., Chen, T.-Y., Zhang, Y., Lin, R. and Yuan, Y. (2025). Practical
Considerations for Using the TITE-BOIN Design to Handle Late-Onset
Toxicity or Fast Accrual in Phase I Trials. *Clinical Cancer Research*,
31(13), 2573-2580.

Takeda, K., Zhu, J. and Hirakawa, A. (2025). BF-BOIN-ET: A Backfill
Bayesian Optimal Interval Design Using Efficacy and Toxicity Outcomes
for Dose Optimization. *Pharmaceutical Statistics*, 24(2), e2470.

Liu, S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase
I Clinical Trials. *Journal of the Royal Statistical Society: Series C*,
64, 507-523.
