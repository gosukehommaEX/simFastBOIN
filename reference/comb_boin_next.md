# Next Dose Combination in a BOIN Drug-Combination Trial

Decide the dose combination for the next cohort of a trial run with the
BOIN drug-combination design of Lin and Yin (2017), from the numbers of
patients and DLTs observed so far at every combination and the
combination that treated the last cohort. The rules are those of
[`comb_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin_simulate.md),
so the decisions taken during the trial are the ones whose operating
characteristics
[`sim_comb_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin.md)
reports.

## Usage

``` r
comb_boin_next(
  target,
  n_pts,
  n_tox,
  current,
  n_earlystop = 18,
  p_saf = NULL,
  p_tox = NULL,
  cutoff_eli = 0.95,
  extrasafe = FALSE,
  offset = 0.05
)
```

## Arguments

- target:

  Numeric scalar. Target DLT probability, for example 0.30.

- n_pts:

  Integer matrix. Number of patients treated so far at each dose
  combination, with one row per dose level of drug A and one column per
  dose level of drug B.

- n_tox:

  Integer matrix of the same shape as `n_pts`. Number of patients with a
  DLT at each combination.

- current:

  Integer vector of length two. Dose levels of drug A and drug B of the
  combination that treated the last cohort.

- n_earlystop:

  Integer scalar. The trial stops once this many patients have been
  treated at the current combination and the design would neither
  escalate nor de-escalate, or cannot. Defaults to 18. Use the value of
  the simulations that support the protocol.

- p_saf:

  Numeric scalar. Highest DLT probability deemed subtherapeutic.
  Defaults to `0.6 * target`.

- p_tox:

  Numeric scalar. Lowest DLT probability deemed overly toxic. Defaults
  to `1.4 * target`.

- cutoff_eli:

  Numeric scalar. Posterior probability cutoff for eliminating a
  combination. Defaults to 0.95.

- extrasafe:

  Logical scalar. Apply the stricter safety stopping rule at the lowest
  combination. Defaults to `FALSE`.

- offset:

  Numeric scalar between 0 and 0.5. Amount by which `cutoff_eli` is
  relaxed for the safety stopping rule. Defaults to 0.05.

## Value

An object of class `comb_boin_next`, which is a list with components

- decision:

  One of `"escalate"`, `"de-escalate"`, `"stay"` and `"stop"`.

- next_dose:

  Integer matrix with columns `dose_a` and `dose_b`: one row for the
  next combination, two rows when two candidates tie, and no row when
  the trial stops.

- tie:

  Whether two candidates tie, in which case the design picks one of them
  at random.

- stop_reason:

  `NA`, or why the trial stops: `"lowest_dose_eliminated"`,
  `"lowest_dose_too_toxic"` or `"n_earlystop"`.

- note:

  A sentence explaining a decision to stay when the rule called for a
  move that is not available, otherwise `NA`.

- candidates:

  Data frame of the two combinations considered for an escalation or a
  de-escalation, one level of drug A and one level of drug B away, with
  their data, the posterior probability that their DLT rate lies between
  `lambda_e` and `lambda_d`, the score used to choose between them and
  whether they are available; `NULL` when no move was considered.

- eliminated:

  Logical matrix of the combinations eliminated by the data.

- boundary:

  The escalation, de-escalation and elimination boundaries at the number
  of patients treated at the current combination.

together with the data, the current combination and the design
parameters.

## Details

The decision follows the same steps as a simulated trial after each
cohort. A combination whose DLT count reaches its elimination boundary
is eliminated together with every combination with at least the same
level of both drugs, and the trial stops if the lowest combination is
eliminated or, with `extrasafe`, judged too toxic. The trial also stops
when `n_earlystop` patients have been treated at the current combination
and the design would neither escalate nor de-escalate. Otherwise the
observed DLT rate at the current combination is compared with the
boundaries of
[`boin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_boundary.md),
and an escalation or a de-escalation moves one level of drug A or of
drug B to the candidate with the higher score: the posterior
probability, under a Beta(0.5, 0.5) prior, that its DLT rate lies
between `lambda_e` and `lambda_d`, plus 0.0005 per patient already
treated there. An eliminated combination is not a candidate for
escalation. When the two scores are equal the design chooses at random,
and the function returns both.

During a simulated trial a combination is eliminated as soon as its data
reach the boundary and is never treated again, so the combinations
eliminated by the data are exactly those eliminated in the course of the
trial.

[`BOIN::next.comb()`](https://rdrr.io/pkg/BOIN/man/next.comb.html)
applies somewhat different rules from the simulations of
[`BOIN::get.oc.comb()`](https://rdrr.io/pkg/BOIN/man/get.oc.comb.html):
it stops as soon as `n.earlystop` patients have been treated at the
current combination, it does not escalate drug A when a combination with
the next level of drug A and at most the current level of drug B has an
observed DLT rate at or above `lambda_d` (and likewise for drug B), and
it draws a random number to break ties. This function follows the
simulations instead, so that a trial is run with the rules whose
operating characteristics were evaluated. Where the rules coincide, the
two functions agree.

## References

Lin, R. and Yin, G. (2017). Bayesian Optimal Interval Design for Dose
Finding in Drug-Combination Trials. Statistical Methods in Medical
Research, 26(5), 2155-2167.

## See also

[`plot.comb_boin_next`](https://gosukehommaex.github.io/simFastBOIN/reference/plot.comb_boin_next.md),
[`boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_decision_table.md),
[`sim_comb_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin.md)

## Examples

``` r
# Three levels of drug A and four of drug B; the last cohort was treated at
# A2B2, where two of three patients had a DLT
n_pts <- matrix(c(6, 6, 0, 0,
                  3, 3, 0, 0,
                  0, 0, 0, 0), nrow = 3, byrow = TRUE)
n_tox <- matrix(c(0, 1, 0, 0,
                  1, 2, 0, 0,
                  0, 0, 0, 0), nrow = 3, byrow = TRUE)

nxt <- comb_boin_next(target = 0.25, n_pts = n_pts, n_tox = n_tox,
                      current = c(2, 2))
nxt
#> BOIN drug-combination design: dose for the next cohort
#>   target DLT rate : 0.25
#>   current         : A2B2, 2 DLTs out of 3 patients (0.667)
#>   boundaries      : escalate if DLTs <= 0, de-escalate if DLTs >= 1, eliminate if DLTs >= 3 (3 patients)
#> 
#> Decision: de-escalate to A1B2
#> 
#> Candidates (interval 0.1968 to 0.2984):
#>  Direction Combination Patients DLTs Pr(in interval)  Score Available
#>     drug A        A1B2        6    1           0.218 0.2211       yes
#>     drug B        A2B1        3    1           0.167 0.1687       yes
#> 
#> Eliminated combinations: none
nxt$next_dose
#>      dose_a dose_b
#> [1,]      1      2
```
