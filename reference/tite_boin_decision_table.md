# Decision Table for the TITE-BOIN Design

Tabulate the dose assignment rule of the time-to-event Bayesian optimal
interval (TITE-BOIN) design, which allows new patients to be treated
while the DLT assessment of earlier patients is still pending. For every
attainable combination of patients treated, DLTs observed and patients
pending at the current dose, the table gives the decision or, when the
decision depends on how long the pending patients have been followed,
the values of the follow-up statistic at which it changes.

## Usage

``` r
tite_boin_decision_table(
  target,
  max_n,
  method = c("imputation", "ess"),
  p_saf = NULL,
  p_tox = NULL,
  cutoff_eli = 0.95,
  max_pending_ratio = NULL,
  min_completed = NULL,
  min_follow_up = 0
)
```

## Arguments

- target:

  Numeric scalar. Target DLT probability, for example 0.30.

- max_n:

  Integer scalar. Largest number of patients treated at a single dose.

- method:

  Character scalar. How the pending patients enter the estimate of the
  DLT probability: `"imputation"` (the default) for the single mean
  imputation of Yuan et al. (2018), or `"ess"` for the effective sample
  size of Lin and Yuan (2020). See the details.

- p_saf:

  Numeric scalar. Highest DLT probability deemed subtherapeutic.
  Defaults to `0.6 * target`.

- p_tox:

  Numeric scalar. Lowest DLT probability deemed overly toxic. Defaults
  to `1.4 * target`.

- cutoff_eli:

  Numeric scalar. Posterior probability cutoff for dose elimination.
  Defaults to 0.95.

- max_pending_ratio:

  Numeric scalar greater than 0 and at most 1. Accrual is suspended when
  the proportion of the patients at the current dose whose assessment is
  pending exceeds this value, unless the dose is de-escalated whatever
  their outcomes. Defaults to 0.5 for `"imputation"`, following Yuan et
  al. (2018), and to 1, which never suspends on this ground, for
  `"ess"`.

- min_completed:

  Integer scalar. While some patients are pending, escalation requires
  at least this many patients at the current dose to have completed the
  assessment, and accrual is suspended instead when fewer have. Defaults
  to 0 for `"imputation"` and to 2 for `"ess"`, following Lin and Yuan
  (2020).

- min_follow_up:

  Numeric scalar between 0 and 1. While some patients are pending,
  escalation also requires the shortest follow-up among the pending
  patients at the current dose, as a fraction of the assessment window
  (MF), to be at least this value, and accrual is suspended instead when
  it is not. Defaults to 0, which never suspends on this ground. Chen et
  al. (2025) recommend 0.25, their rule 2 with B = 25. The rule leaves
  the boundaries on the follow-up statistic unchanged, so it is recorded
  as an attribute rather than in the rows of the table.

## Value

An object of class `tite_boin_decision_table`, which is a data frame
with one row per attainable state of the current dose, ordered by `n`,
`n_tox` and `n_pending`, and with the columns

- n:

  Number of patients treated at the current dose, from 1 to `max_n`.

- n_tox:

  Number of DLTs observed, which can only come from patients who
  completed the assessment.

- n_pending:

  Number of patients whose assessment is pending.

- decision:

  Decision code. `"E"` escalates, `"S"` stays, `"D"` de-escalates,
  `"DE"` de-escalates and eliminates the dose together with all higher
  doses, and `"SUS"` suspends accrual until more data are available.
  These do not depend on the follow-up of the pending patients. The
  codes `"E/S"`, `"S/D"` and `"E/S/D"` mark states whose decision
  depends on the follow-up statistic through `esc_bound` and
  `deesc_bound`; in them `"SUS"` takes the place of `"E"` when
  escalation is blocked by `min_completed`.

- esc_bound:

  Escalate, or suspend when escalation is blocked, when the follow-up
  statistic is at least this value. `NA` when no such boundary applies.

- deesc_bound:

  De-escalate when the follow-up statistic is at most this value. `NA`
  when no such boundary applies.

The design parameters are stored as attributes: `method`, `statistic`
(`"STFT"` or `"ESS"`), `target`, `p_saf`, `p_tox`, `lambda_e`,
`lambda_d`, `cutoff_eli`, `max_pending_ratio`, `min_completed` and
`min_follow_up`. When `min_follow_up` is positive, an escalation (code
`"E"`, or `"E"` as the first part of a code with a slash) in a state
with pending patients takes place only when MF is at least
`min_follow_up`, and accrual is suspended otherwise. The object has
`print` and `plot` methods.

## Details

Suppose that n patients have been treated at the current dose, that r of
them have completed the DLT assessment with `n_tox` DLTs, and that c =
n - r are pending. The standardized total follow-up time (STFT) is the
sum of the follow-up times of the pending patients divided by the length
of the assessment window, so that STFT lies between 0 and c. With no
pending patient the decision is that of the BOIN design, taken from
[`boin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_boundary.md).

`method = "imputation"` replaces each pending outcome by its expected
value under a uniform time to DLT. The estimate of the DLT probability
is `(n_tox + q * (c - STFT)) / n`, where `q = p / (1 - p)` and
`p = (n_tox + a) / (r + 1)` is the posterior mean under a Beta(a, 1 - a)
prior with `a = target / 2` (Yuan et al., 2018, Supplementary Appendix
A). The estimate decreases as STFT grows, so the rule becomes a pair of
boundaries on STFT. As in the BOIN design, escalation is possible only
while `n_tox / n` is below the target and de-escalation only once it
reaches the target.

`method = "ess"` uses the approximated likelihood of Lin and Yuan
(2020), under which the estimate is `n_tox / ESS` with the effective
sample size ESS = r + STFT, which lies between r and n. The dose is
escalated when ESS is at least `n_tox / lambda_e` and de-escalated when
it is at most `n_tox / lambda_d`. The boundaries in the table refer to
ESS.

The rules are applied in this order: elimination; de-escalation that
holds whatever the outcomes of the pending patients, that is when
`n_tox / n` is at least `lambda_d`; suspension for too many pending
patients; and finally the boundaries on the follow-up statistic, in
which an escalation blocked by `min_completed` becomes a suspension. An
escalation that remains is then subject to `min_follow_up`. Elimination
requires `Pr(p > target | data) > cutoff_eli` computed with all n
treated patients, the pending ones counted as without DLT, which is how
both articles define it. It is not evaluated before three patients have
been treated.

Neither the length of the assessment window nor the weights given to the
follow-up times change the table, because both enter only through STFT.
A weighted STFT, such as that of a piecewise uniform time to DLT, is
compared with the same boundaries.

## Agreement with the published tables

With the defaults, `method = "imputation"` reproduces every entry of
Table 1 (target 0.2) of Yuan et al. (2018) and of Table S1 (target 0.3)
of its supplementary appendix. These tables allow de-escalation when
`n_tox / n` equals the target (Table 1, 15 patients with 3 DLTs), which
equation (5) of the appendix, read literally, would not. The
implementation follows the tables, which agree with the verbal
description of the rule in the appendix.

With `max_pending_ratio = 0.49` and `min_follow_up = 0.25`,
`method = "imputation"` reproduces every entry of Table 1 (target 0.3,
up to six patients) of Chen et al. (2025) and of Table A1 (target 0.25,
up to nine patients) of the supplementary materials of Chen et al.
(2026). Their rule 1 suspends accrual when fewer than 51 percent of the
patients have completed the assessment, which for whole numbers of
patients is the same as more than 49 percent pending, unless the
observed DLT rate already exceeds `lambda_d`; this exception is the
de-escalation that holds whatever the outcomes of the pending patients.

For `method = "ess"` no table for the BOIN design has been published;
the tables of Lin and Yuan (2020) are for the keyboard and mTPI designs.

## References

Yuan, Y., Lin, R., Li, D., Nie, L. and Warren, K. E. (2018).
Time-to-Event Bayesian Optimal Interval Design to Accelerate Phase I
Trials. Clinical Cancer Research, 24(20), 4921-4930.

Lin, R. and Yuan, Y. (2020). Time-to-Event Model-Assisted Designs for
Dose-Finding Trials with Delayed Toxicity. Biostatistics, 21(4),
807-824.

Chen, K., Chen, T.-Y., Zhang, Y., Lin, R. and Yuan, Y. (2025). Practical
Considerations for Using the TITE-BOIN Design to Handle Late-Onset
Toxicity or Fast Accrual in Phase I Trials. Clinical Cancer Research,
31(13), 2573-2580.

Chen, K., Zhao, Y., Takeda, K. and Yuan, Y. (2026). BE-BOIN: A Dose
Optimization Design Accommodating Backfill and Late-Onset Toxicity.
Therapeutic Innovation and Regulatory Science.
[doi:10.1007/s43441-026-00994-0](https://doi.org/10.1007/s43441-026-00994-0)

## See also

[`print.tite_boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/print.tite_boin_decision_table.md),
[`plot.tite_boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/plot.tite_boin_decision_table.md),
[`boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_decision_table.md)

## Examples

``` r
# Table 1 of Yuan et al. (2018)
decisions <- tite_boin_decision_table(target = 0.2, max_n = 15)
print(decisions, cohort_size = 3)
#> TITE-BOIN decision table with target 0.2
#> Method: single mean imputation (Yuan et al., 2018); boundaries refer to STFT
#> 
#>  Patients DLTs Pending Escalate Stay    De-escalate
#>  3        0    <= 1    Y                           
#>  3        0    >= 2    Suspend  Suspend Suspend    
#>  3        1    <= 2                     Y          
#>  3        >= 2 <= 1                     Y & Elim   
#>  6        0    <= 3    Y                           
#>  6        0    >= 4    Suspend  Suspend Suspend    
#>  6        1    <= 3             Y                  
#>  6        1    >= 4    Suspend  Suspend Suspend    
#>  6        2    <= 4                     Y          
#>  6        >= 3 <= 3                     Y & Elim   
#>  9        0    <= 4    Y                           
#>  9        0    >= 5    Suspend  Suspend Suspend    
#>  9        1    <= 2    Y                           
#>  9        1    3       >= 0.77  < 0.77             
#>  9        1    4       >= 2.15  < 2.15             
#>  9        1    >= 5    Suspend  Suspend Suspend    
#>  9        2    0                Y                  
#>  9        2    1                > 0.52  <= 0.52    
#>  9        2    2                > 1.59  <= 1.59    
#>  9        2    3                > 2.66  <= 2.66    
#>  9        2    4                > 3.73  <= 3.73    
#>  9        2    >= 5    Suspend  Suspend Suspend    
#>  9        3    <= 6                     Y          
#>  9        >= 4 <= 5                     Y & Elim   
#>  12       0    <= 6    Y                           
#>  12       0    >= 7    Suspend  Suspend Suspend    
#>  12       1    <= 5    Y                           
#>  12       1    6       >= 1.24  < 1.24             
#>  12       1    >= 7    Suspend  Suspend Suspend    
#>  12       2    <= 6             Y                  
#>  12       2    >= 7    Suspend  Suspend Suspend    
#>  12       3, 4 <= 9                     Y          
#>  12       >= 5 <= 7                     Y & Elim   
#>  15       0    <= 7    Y                           
#>  15       0    >= 8    Suspend  Suspend Suspend    
#>  15       1    <= 7    Y                           
#>  15       1    >= 8    Suspend  Suspend Suspend    
#>  15       2    <= 2    Y                           
#>  15       2    3       >= 1.14  < 1.14             
#>  15       2    4       >= 2.31  < 2.31             
#>  15       2    5       >= 3.48  < 3.48             
#>  15       2    6       >= 4.65  < 4.65             
#>  15       2    7       >= 5.82  < 5.82             
#>  15       2    >= 8    Suspend  Suspend Suspend    
#>  15       3    <= 2             Y                  
#>  15       3    3                > 1.16  <= 1.16    
#>  15       3    4                > 2.34  <= 2.34    
#>  15       3    5                > 3.53  <= 3.53    
#>  15       3    6                > 4.72  <= 4.72    
#>  15       3    7                > 5.90  <= 5.90    
#>  15       3    >= 8    Suspend  Suspend Suspend    
#>  15       4, 5 <= 11                    Y          
#>  15       >= 6 <= 9                     Y & Elim   
#> 
#> Patients = treated at the current dose, DLTs = DLTs observed, Pending = assessment not yet completed
#> STFT = total follow-up time of the pending patients divided by the length of the assessment window
#> Accrual is suspended when more than 50% of the patients are pending
#> Y & Elim = de-escalate and eliminate this dose and all higher doses

# Nine patients, one DLT and four pending: escalate once STFT reaches 2.15
decisions[decisions$n == 9 & decisions$n_tox == 1 & decisions$n_pending == 4, ]
#> TITE-BOIN decision table with target 0.2
#> Method: single mean imputation (Yuan et al., 2018); boundaries refer to STFT
#> 
#>  Patients DLTs Pending Escalate Stay   De-escalate
#>  9        1    4       >= 2.15  < 2.15            
#> 
#> Patients = treated at the current dose, DLTs = DLTs observed, Pending = assessment not yet completed
#> STFT = total follow-up time of the pending patients divided by the length of the assessment window
#> Accrual is suspended when more than 50% of the patients are pending
#> Y & Elim = de-escalate and eliminate this dose and all higher doses

# Table 1 of Chen et al. (2025): at least 51 percent of the patients with a
# completed assessment, and escalation only once every pending patient has
# been followed for a quarter of the window
decisions_chen <- tite_boin_decision_table(target = 0.3, max_n = 6,
                                           max_pending_ratio = 0.49,
                                           min_follow_up = 0.25)
print(decisions_chen, cohort_size = 3)
#> TITE-BOIN decision table with target 0.3
#> Method: single mean imputation (Yuan et al., 2018); boundaries refer to STFT
#> 
#>  Patients DLTs Pending Escalate             Stay    De-escalate
#>  3        0    0       Y                                       
#>  3        0    1       MF >= 0.25                              
#>  3        0    >= 2    Suspend              Suspend Suspend    
#>  3        1    0                            Y                  
#>  3        1    1                            > 0.88  <= 0.88    
#>  3        1    >= 2    Suspend              Suspend Suspend    
#>  3        2    <= 1                                 Y          
#>  3        3    0                                    Y & Elim   
#>  6        0    0       Y                                       
#>  6        0    1-2     MF >= 0.25                              
#>  6        0    >= 3    Suspend              Suspend Suspend    
#>  6        1    0       Y                                       
#>  6        1    1       MF >= 0.25                              
#>  6        1    2       >= 0.60 & MF >= 0.25 < 0.60             
#>  6        1    >= 3    Suspend              Suspend Suspend    
#>  6        2    0                            Y                  
#>  6        2    1                            > 0.73  <= 0.73    
#>  6        2    2                            > 1.80  <= 1.80    
#>  6        2    >= 3    Suspend              Suspend Suspend    
#>  6        3    <= 3                                 Y          
#>  6        >= 4 <= 2                                 Y & Elim   
#> 
#> Patients = treated at the current dose, DLTs = DLTs observed, Pending = assessment not yet completed
#> STFT = total follow-up time of the pending patients divided by the length of the assessment window
#> Accrual is suspended when more than 49% of the patients are pending
#> MF = shortest follow-up time of the pending patients divided by the length of the assessment window
#> While patients are pending, escalation requires MF >= 0.25; otherwise accrual is suspended
#> Y & Elim = de-escalate and eliminate this dose and all higher doses

# The effective sample size of Lin and Yuan (2020)
decisions_ess <- tite_boin_decision_table(target = 0.3, max_n = 12, method = "ess")
print(decisions_ess, cohort_size = 3)
#> TITE-BOIN decision table with target 0.3
#> Method: effective sample size (Lin and Yuan, 2020); boundaries refer to ESS
#> 
#>  Patients DLTs Pending Escalate           Stay           De-escalate
#>  3        0    <= 1    Y                                            
#>  3        0    >= 2    Suspend            Suspend        Suspend    
#>  3        1    0                          Y                         
#>  3        1    >= 1                       > 2.79         <= 2.79    
#>  3        2    <= 1                                      Y          
#>  3        3    0                                         Y & Elim   
#>  6        0    <= 4    Y                                            
#>  6        0    >= 5    Suspend            Suspend        Suspend    
#>  6        1    <= 1    Y                                            
#>  6        1    2-3     >= 4.23            < 4.23                    
#>  6        1    4       >= 4.23            > 2.79, < 4.23 <= 2.79    
#>  6        1    >= 5    Suspend if >= 4.23 > 2.79, < 4.23 <= 2.79    
#>  6        2    0                          Y                         
#>  6        2    >= 1                       > 5.58         <= 5.58    
#>  6        3    <= 3                                      Y          
#>  6        >= 4 <= 2                                      Y & Elim   
#>  9        0    <= 7    Y                                            
#>  9        0    >= 8    Suspend            Suspend        Suspend    
#>  9        1    <= 4    Y                                            
#>  9        1    5-6     >= 4.23            < 4.23                    
#>  9        1    7       >= 4.23            > 2.79, < 4.23 <= 2.79    
#>  9        1    >= 8    Suspend if >= 4.23 > 2.79, < 4.23 <= 2.79    
#>  9        2    0       Y                                            
#>  9        2    1-3     >= 8.46            < 8.46                    
#>  9        2    >= 4    >= 8.46            > 5.58, < 8.46 <= 5.58    
#>  9        3    0                          Y                         
#>  9        3    >= 1                       > 8.37         <= 8.37    
#>  9        4    <= 5                                      Y          
#>  9        >= 5 <= 4                                      Y & Elim   
#>  12       0    <= 10   Y                                            
#>  12       0    >= 11   Suspend            Suspend        Suspend    
#>  12       1    <= 7    Y                                            
#>  12       1    8-9     >= 4.23            < 4.23                    
#>  12       1    10      >= 4.23            > 2.79, < 4.23 <= 2.79    
#>  12       1    >= 11   Suspend if >= 4.23 > 2.79, < 4.23 <= 2.79    
#>  12       2    <= 3    Y                                            
#>  12       2    4-6     >= 8.46            < 8.46                    
#>  12       2    >= 7    >= 8.46            > 5.58, < 8.46 <= 5.58    
#>  12       3    <= 3                       Y                         
#>  12       3    >= 4                       > 8.37         <= 8.37    
#>  12       4    0                          Y                         
#>  12       4    >= 1                       > 11.16        <= 11.16   
#>  12       5, 6 <= 7                                      Y          
#>  12       >= 7 <= 5                                      Y & Elim   
#> 
#> Patients = treated at the current dose, DLTs = DLTs observed, Pending = assessment not yet completed
#> STFT = total follow-up time of the pending patients divided by the length of the assessment window
#> ESS = number of patients with a completed assessment plus STFT
#> While patients are pending, escalation requires at least 2 patients with a completed assessment
#> Y & Elim = de-escalate and eliminate this dose and all higher doses
```
