# simFastBOIN: Fast Simulation of Bayesian Optimal Interval Designs

Tools for designing and evaluating phase I dose-finding trials that use
the Bayesian optimal interval (BOIN) design of Liu and Yuan (2015). The
simulation engine is written in C++ and reproduces the reference
implementation in the BOIN package exactly, including the order in which
random variates are consumed, so that results obtained with the same
seed agree to the last trial.

## Design tools

- [`boin_lambda`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_lambda.md):

  Escalation and de-escalation interval boundaries.

- [`boin_p_tox`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_p_tox.md):

  The toxic threshold giving a required de-escalation boundary.

- [`boin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_boundary.md):

  Integer decision boundaries as a function of sample size.

- [`boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_decision_table.md):

  Decision table indexed by DLTs and patients, with print and plot
  methods.

- [`boin_stopping_table`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_stopping_table.md):

  Safety stopping boundary as a two-row table.

- [`tite_boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_decision_table.md):

  Decision table of the time-to-event BOIN design, which treats new
  patients while earlier DLT data are pending, with print and plot
  methods.

## Simulation tools

- [`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md):

  Operating characteristics for one scenario.

- [`sim_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_multi.md):

  Operating characteristics for several scenarios.

- [`boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_simulate.md):

  Raw trial data from the simulation engine.

- [`boin_isotonic`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_isotonic.md):

  Isotonic estimates of the DLT probability.

- [`boin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_select_mtd.md):

  MTD selection from completed trials.

- [`sim_tite_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin.md):

  Operating characteristics of the time-to-event BOIN design for one
  scenario, with trial duration and suspensions of accrual.

- [`sim_tite_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin_multi.md):

  The same for several scenarios.

- [`tite_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_simulate.md):

  Raw trial data from the time-to-event BOIN engine, with patients
  arriving over time and DLTs observed late.

## Backfilling

- [`sim_bf_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin.md):

  Operating characteristics of the BOIN design with backfilling
  (BF-BOIN) for one scenario of DLT and response probabilities.

- [`sim_bf_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin_multi.md):

  The same for several scenarios.

- [`bf_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/bf_boin_simulate.md):

  Raw trial data from the BF-BOIN engine, with backfilled patients and
  responses.

- [`sim_be_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_be_boin.md):

  Operating characteristics of BE-BOIN, which combines backfilling with
  the time-to-event BOIN design, for one scenario.

- [`sim_be_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_be_boin_multi.md):

  The same for several scenarios.

- [`be_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/be_boin_simulate.md):

  Raw trial data from the BE-BOIN engine.

## Drug combinations

- [`sim_comb_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin.md):

  Operating characteristics of the BOIN drug-combination design of Lin
  and Yin (2017) for one matrix of DLT probabilities.

- [`sim_comb_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin_multi.md):

  The same for several matrices.

- [`comb_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin_simulate.md):

  Raw trial data from the drug-combination engine.

- [`comb_boin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin_select_mtd.md):

  Selection of the MTD combination from completed trials.

- [`comb_boin_next`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin_next.md):

  The combination for the next cohort of a trial in progress, with print
  and plot methods.

## The 3+3 design

- [`oc_3p3`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_3p3.md):

  Operating characteristics of the traditional 3+3 design, in closed
  form.

- [`sim_3p3`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_3p3.md):

  The same design by simulation, as a check on the closed form.

- [`decision_table_3p3`](https://gosukehommaex.github.io/simFastBOIN/reference/decision_table_3p3.md):

  Decision table of the 3+3 design, with print and plot methods.

## References

Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase
I Clinical Trials. Journal of the Royal Statistical Society: Series C,
64, 507-523.

Yan, F., Zhang, L., Zhou, Y., Pan, H., Liu, S. and Yuan, Y. (2020).
BOIN: An R Package for Designing Single-Agent and Drug-Combination
Dose-Finding Trials Using Bayesian Optimal Interval Designs. Journal of
Statistical Software, 94(13), 1-32.

Yuan, Y., Lin, R., Li, D., Nie, L. and Warren, K. E. (2018).
Time-to-Event Bayesian Optimal Interval Design to Accelerate Phase I
Trials. Clinical Cancer Research, 24(20), 4921-4930.

Lin, R. and Yuan, Y. (2020). Time-to-Event Model-Assisted Designs for
Dose-Finding Trials with Delayed Toxicity. Biostatistics, 21(4),
807-824.

Zhao, Y., Yuan, Y., Korn, E. L. and Freidlin, B. (2024). Backfilling
Patients in Phase I Dose-Escalation Trials Using Bayesian Optimal
Interval Design (BOIN). Clinical Cancer Research, 30(4), 673-679.

Chen, K., Zhao, Y., Takeda, K. and Yuan, Y. (2026). BE-BOIN: A Dose
Optimization Design Accommodating Backfill and Late-Onset Toxicity.
Therapeutic Innovation and Regulatory Science.
[doi:10.1007/s43441-026-00994-0](https://doi.org/10.1007/s43441-026-00994-0)

Lin, R. and Yin, G. (2017). Bayesian Optimal Interval Design for Dose
Finding in Drug-Combination Trials. Statistical Methods in Medical
Research, 26(5), 2155-2167.

## See also

Useful links:

- <https://gosukehommaex.github.io/simFastBOIN/>

- <https://github.com/gosukehommaEX/simFastBOIN>

- Report bugs at <https://github.com/gosukehommaEX/simFastBOIN/issues>

## Author

**Maintainer**: Gosuke Homma <my.name.is.gosuke@gmail.com>
([ORCID](https://orcid.org/0000-0002-6854-5627))

Authors:

- Gosuke Homma <my.name.is.gosuke@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-6854-5627))
