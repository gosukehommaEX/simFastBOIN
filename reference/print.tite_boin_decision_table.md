# Print a TITE-BOIN Decision Table

Display the decision table produced by
[`tite_boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_decision_table.md)
in the layout of the published tables: one row per number of patients
treated, number of DLTs and range of pending patients sharing a
decision, with the boundaries on the follow-up statistic in the
Escalate, Stay and De-escalate columns.

## Usage

``` r
# S3 method for class 'tite_boin_decision_table'
print(x, cohort_size = NULL, digits = 2, ...)
```

## Arguments

- x:

  An object of class `tite_boin_decision_table`.

- cohort_size:

  Integer scalar or `NULL`. When supplied, only the sample sizes that
  are multiples of `cohort_size` are shown, which is how the published
  tables are laid out.

- digits:

  Integer scalar. Number of decimal places of the boundaries. Defaults
  to 2, as in the published tables.

- ...:

  Further arguments, currently ignored.

## Value

The object `x`, invisibly.

## Details

A boundary shown as `">= 2.15"` under Escalate means that the dose is
escalated when the follow-up statistic is at least 2.15; the statistic
is STFT for `method = "imputation"` and ESS for `method = "ess"`.
"Suspend if \>= 4.23" under Escalate marks an escalation blocked by
`min_completed`. When the table was built with a positive
`min_follow_up`, an escalation with pending patients carries the
condition on MF, for example `">= 0.22 & MF >= 0.25"`; accrual is
suspended when the condition on the follow-up statistic holds but the
one on MF does not.

A data frame that no longer has the columns or attributes of a decision
table, for example after selecting columns, is printed as a plain data
frame.

## See also

[`tite_boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_decision_table.md),
[`plot.tite_boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/plot.tite_boin_decision_table.md)

## Examples

``` r
decisions <- tite_boin_decision_table(target = 0.3, max_n = 15)

# Table S1 of Yuan et al. (2018)
print(decisions, cohort_size = 3)
#> TITE-BOIN decision table with target 0.3
#> Method: single mean imputation (Yuan et al., 2018); boundaries refer to STFT
#> 
#>  Patients DLTs Pending Escalate Stay    De-escalate
#>  3        0    <= 1    Y                           
#>  3        0    >= 2    Suspend  Suspend Suspend    
#>  3        1    0                Y                  
#>  3        1    1                > 0.88  <= 0.88    
#>  3        1    >= 2    Suspend  Suspend Suspend    
#>  3        2    <= 1                     Y          
#>  3        3    0                        Y & Elim   
#>  6        0    <= 3    Y                           
#>  6        0    >= 4    Suspend  Suspend Suspend    
#>  6        1    <= 1    Y                           
#>  6        1    2       >= 0.60  < 0.60             
#>  6        1    3       >= 1.96  < 1.96             
#>  6        1    >= 4    Suspend  Suspend Suspend    
#>  6        2    0                Y                  
#>  6        2    1                > 0.73  <= 0.73    
#>  6        2    2                > 1.80  <= 1.80    
#>  6        2    3                > 2.87  <= 2.87    
#>  6        2    >= 4    Suspend  Suspend Suspend    
#>  6        3    <= 3                     Y          
#>  6        >= 4 <= 2                     Y & Elim   
#>  9        0    <= 4    Y                           
#>  9        0    >= 5    Suspend  Suspend Suspend    
#>  9        1    <= 4    Y                           
#>  9        1    >= 5    Suspend  Suspend Suspend    
#>  9        2    0       Y                           
#>  9        2    1       >= 0.59  < 0.59             
#>  9        2    2       >= 1.65  < 1.65             
#>  9        2    3       >= 2.71  < 2.71             
#>  9        2    4       >= 3.77  < 3.77             
#>  9        2    >= 5    Suspend  Suspend Suspend    
#>  9        3    0                Y                  
#>  9        3    1                > 0.58  <= 0.58    
#>  9        3    2                > 1.65  <= 1.65    
#>  9        3    3                > 2.72  <= 2.72    
#>  9        3    4                > 3.79  <= 3.79    
#>  9        3    >= 5    Suspend  Suspend Suspend    
#>  9        4    <= 5                     Y          
#>  9        >= 5 <= 4                     Y & Elim   
#>  12       0    <= 6    Y                           
#>  12       0    >= 7    Suspend  Suspend Suspend    
#>  12       1    <= 6    Y                           
#>  12       1    >= 7    Suspend  Suspend Suspend    
#>  12       2    <= 3    Y                           
#>  12       2    4       >= 1.33  < 1.33             
#>  12       2    5       >= 2.72  < 2.72             
#>  12       2    6       >= 4.11  < 4.11             
#>  12       2    >= 7    Suspend  Suspend Suspend    
#>  12       3    <= 6             Y                  
#>  12       3    >= 7    Suspend  Suspend Suspend    
#>  12       4    0                Y                  
#>  12       4    1                > 0.43  <= 0.43    
#>  12       4    2                > 1.50  <= 1.50    
#>  12       4    3                > 2.57  <= 2.57    
#>  12       4    4                > 3.65  <= 3.65    
#>  12       4    5                > 4.72  <= 4.72    
#>  12       4    6                > 5.79  <= 5.79    
#>  12       4    >= 7    Suspend  Suspend Suspend    
#>  12       5, 6 <= 7                     Y          
#>  12       >= 7 <= 5                     Y & Elim   
#>  15       0    <= 7    Y                           
#>  15       0    >= 8    Suspend  Suspend Suspend    
#>  15       1    <= 7    Y                           
#>  15       1    >= 8    Suspend  Suspend Suspend    
#>  15       2    <= 5    Y                           
#>  15       2    6       >= 0.35  < 0.35             
#>  15       2    7       >= 2.07  < 2.07             
#>  15       2    >= 8    Suspend  Suspend Suspend    
#>  15       3    <= 1    Y                           
#>  15       3    2       >= 0.11  < 0.11             
#>  15       3    3       >= 1.29  < 1.29             
#>  15       3    4       >= 2.46  < 2.46             
#>  15       3    5       >= 3.64  < 3.64             
#>  15       3    6       >= 4.81  < 4.81             
#>  15       3    7       >= 5.98  < 5.98             
#>  15       3    >= 8    Suspend  Suspend Suspend    
#>  15       4    <= 7             Y                  
#>  15       4    >= 8    Suspend  Suspend Suspend    
#>  15       5    0                Y                  
#>  15       5    1                > 0.28  <= 0.28    
#>  15       5    2                > 1.35  <= 1.35    
#>  15       5    3                > 2.42  <= 2.42    
#>  15       5    4                > 3.50  <= 3.50    
#>  15       5    5                > 4.57  <= 4.57    
#>  15       5    6                > 5.64  <= 5.64    
#>  15       5    7                > 6.72  <= 6.72    
#>  15       5    >= 8    Suspend  Suspend Suspend    
#>  15       6, 7 <= 9                     Y          
#>  15       >= 8 <= 7                     Y & Elim   
#> 
#> Patients = treated at the current dose, DLTs = DLTs observed, Pending = assessment not yet completed
#> STFT = total follow-up time of the pending patients divided by the length of the assessment window
#> Accrual is suspended when more than 50% of the patients are pending
#> Y & Elim = de-escalate and eliminate this dose and all higher doses

# Boundaries with three decimal places
print(decisions[decisions$n == 9, ], digits = 3)
#> TITE-BOIN decision table with target 0.3
#> Method: single mean imputation (Yuan et al., 2018); boundaries refer to STFT
#> 
#>  Patients DLTs Pending Escalate Stay    De-escalate
#>  9        0    <= 4    Y                           
#>  9        0    >= 5    Suspend  Suspend Suspend    
#>  9        1    <= 4    Y                           
#>  9        1    >= 5    Suspend  Suspend Suspend    
#>  9        2    0       Y                           
#>  9        2    1       >= 0.591 < 0.591            
#>  9        2    2       >= 1.651 < 1.651            
#>  9        2    3       >= 2.710 < 2.710            
#>  9        2    4       >= 3.770 < 3.770            
#>  9        2    >= 5    Suspend  Suspend Suspend    
#>  9        3    0                Y                  
#>  9        3    1                > 0.579 <= 0.579   
#>  9        3    2                > 1.651 <= 1.651   
#>  9        3    3                > 2.723 <= 2.723   
#>  9        3    4                > 3.795 <= 3.795   
#>  9        3    >= 5    Suspend  Suspend Suspend    
#>  9        4    <= 5                     Y          
#>  9        >= 5 <= 4                     Y & Elim   
#> 
#> Patients = treated at the current dose, DLTs = DLTs observed, Pending = assessment not yet completed
#> STFT = total follow-up time of the pending patients divided by the length of the assessment window
#> Accrual is suspended when more than 50% of the patients are pending
#> Y & Elim = de-escalate and eliminate this dose and all higher doses
```
