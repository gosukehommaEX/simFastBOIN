# simFastBOIN: Fast Simulation of Bayesian Optimal Interval Designs

Tools for designing and evaluating phase I dose-finding trials that use
the Bayesian optimal interval (BOIN) design of Liu and Yuan (2015). The
simulation engine is written in C++ and reproduces the reference
implementation in the BOIN package exactly, including the order in which
random variates are consumed, so that results obtained with the same
seed agree to the last trial.

## Design tools

- [`boin_lambda`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_lambda.md):

  Escalation and de-escalation interval boundaries.

- [`boin_p_tox`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_p_tox.md):

  The toxic threshold giving a required de-escalation boundary.

- [`boin_boundary`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_boundary.md):

  Integer decision boundaries as a function of sample size.

- [`boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_decision_table.md):

  Decision table indexed by DLTs and patients, with print and plot
  methods.

- [`boin_stopping_table`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_stopping_table.md):

  Safety stopping boundary as a two-row table.

- [`tite_boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/tite_boin_decision_table.md):

  Decision table of the time-to-event BOIN design, which treats new
  patients while earlier DLT data are pending, with print and plot
  methods.

## Simulation tools

- [`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin.md):

  Operating characteristics for one scenario.

- [`sim_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_multi.md):

  Operating characteristics for several scenarios.

- [`boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_simulate.md):

  Raw trial data from the simulation engine.

- [`boin_isotonic`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_isotonic.md):

  Isotonic estimates of the DLT probability.

- [`boin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_select_mtd.md):

  MTD selection from completed trials.

- [`sim_tite_boin`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_tite_boin.md):

  Operating characteristics of the time-to-event BOIN design for one
  scenario, with trial duration and suspensions of accrual.

- [`sim_tite_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_tite_boin_multi.md):

  The same for several scenarios.

- [`tite_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/tite_boin_simulate.md):

  Raw trial data from the time-to-event BOIN engine, with patients
  arriving over time and DLTs observed late.

## Backfilling

- [`sim_bf_boin`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_bf_boin.md):

  Operating characteristics of the BOIN design with backfilling
  (BF-BOIN) for one scenario of DLT and response probabilities.

- [`sim_bf_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_bf_boin_multi.md):

  The same for several scenarios.

- [`bf_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/bf_boin_simulate.md):

  Raw trial data from the BF-BOIN engine, with backfilled patients and
  responses.

- [`sim_be_boin`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_be_boin.md):

  Operating characteristics of BE-BOIN, which combines backfilling with
  the time-to-event BOIN design, for one scenario.

- [`sim_be_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_be_boin_multi.md):

  The same for several scenarios.

- [`be_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/be_boin_simulate.md):

  Raw trial data from the BE-BOIN engine.

## Drug combinations

- [`sim_comb_boin`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_comb_boin.md):

  Operating characteristics of the BOIN drug-combination design of Lin
  and Yin (2017) for one matrix of DLT probabilities.

- [`sim_comb_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_comb_boin_multi.md):

  The same for several matrices.

- [`comb_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/comb_boin_simulate.md):

  Raw trial data from the drug-combination engine.

- [`comb_boin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/comb_boin_select_mtd.md):

  Selection of the MTD combination from completed trials.

- [`comb_boin_next`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/comb_boin_next.md):

  The combination for the next cohort of a trial in progress, with print
  and plot methods.

- [`sim_waterfall`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_waterfall.md):

  Operating characteristics of the waterfall design of Zhang and Yuan
  (2016), which finds the MTD contour, for one matrix of DLT
  probabilities.

- [`sim_waterfall_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_waterfall_multi.md):

  The same for several matrices.

- [`waterfall_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/waterfall_simulate.md):

  Raw trial data from the waterfall engine.

- [`waterfall_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/waterfall_select_mtd.md):

  Selection of the MTD contour from completed trials.

## Optimal biological dose (BOIN-ET)

- [`boin_et_boundary`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_boundary.md):

  Boundaries of the BOIN-ET design of Takeda, Taguri and Morita (2018),
  which decides from both the DLT and the response.

- [`boin_et_decision_table`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_decision_table.md):

  Decision table indexed by DLTs and responses, with a print method.

- [`sim_boin_et`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_et.md):

  Operating characteristics of BOIN-ET for one scenario of DLT and
  response probabilities, with the selection of the optimal biological
  dose (OBD).

- [`sim_boin_et_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_et_multi.md):

  The same for several scenarios.

- [`boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_simulate.md):

  Raw trial data from the BOIN-ET engine.

- [`boin_et_select_obd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_select_obd.md):

  Selection of the OBD from completed trials, with the estimators and
  selection rules of the 'boinet' package.

- [`sim_tite_boin_et`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_tite_boin_et.md):

  Operating characteristics of the time-to-event TITE-BOIN-ET design of
  Takeda, Morita and Taguri (2020).

- [`sim_tite_boin_et_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_tite_boin_et_multi.md):

  The same for several scenarios.

- [`tite_boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/tite_boin_et_simulate.md):

  Raw trial data from the TITE-BOIN-ET engine.

- [`sim_bf_boin_et`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_bf_boin_et.md):

  Operating characteristics of the backfill BF-BOIN-ET design of Takeda,
  Zhu and Hirakawa (2025).

- [`sim_bf_boin_et_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_bf_boin_et_multi.md):

  The same for several scenarios.

## Optimal biological dose (BOIN12 and TITE-BOIN12)

- [`boin12_decision_table`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin12_decision_table.md):

  Desirability scores of the BOIN12 design of Lin et al. (2020), indexed
  by patients, DLTs and responses.

- [`sim_boin12`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin12.md):

  Operating characteristics of BOIN12 for one scenario of DLT and
  response probabilities, with the selection of the OBD.

- [`sim_boin12_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin12_multi.md):

  The same for several scenarios.

- [`boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin12_simulate.md):

  Raw trial data from the BOIN12 engine.

- [`boin12_select_obd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin12_select_obd.md):

  OBD selection from completed trials.

- [`sim_tite_boin12`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_tite_boin12.md):

  Operating characteristics of the time-to-event TITE-BOIN12 design of
  Zhou et al. (2022).

- [`sim_tite_boin12_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_tite_boin12_multi.md):

  The same for several scenarios.

- [`tite_boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/tite_boin12_simulate.md):

  Raw trial data from the TITE-BOIN12 engine.

## Optimal biological dose (U-BOIN)

- [`uboin_utility_table`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/uboin_utility_table.md):

  Posterior mean utility and admissibility of a dose for every outcome
  of its patients, the decision tables of the U-BOIN design of Zhou, Lee
  and Yuan (2019).

- [`sim_uboin`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_uboin.md):

  Operating characteristics of U-BOIN for one scenario of DLT and
  response probabilities, with the selection of the OBD.

- [`sim_uboin_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_uboin_multi.md):

  The same for several scenarios.

- [`uboin_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/uboin_simulate.md):

  Raw trial data from the U-BOIN engine.

- [`uboin_select_obd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/uboin_select_obd.md):

  OBD selection from completed trials.

## The 3+3 design

- [`oc_3p3`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_3p3.md):

  Operating characteristics of the traditional 3+3 design, in closed
  form.

- [`sim_3p3`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_3p3.md):

  The same design by simulation, as a check on the closed form.

- [`decision_table_3p3`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/decision_table_3p3.md):

  Decision table of the 3+3 design, with print and plot methods.

## Tables and figures

- [`oc_scenarios`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_scenarios.md):

  True scenarios with the utilities and the true MTDs and OBDs, with
  print and plot methods.

- [`oc_results`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_results.md):

  Results of several designs in one long table, with print and plot
  methods.

- [`oc_compare`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_compare.md):

  Several designs judged against the true doses, with print and plot
  methods.

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

Takeda, K., Taguri, M. and Morita, S. (2018). BOIN-ET: Bayesian optimal
interval design for dose finding based on both efficacy and toxicity
outcomes. Pharmaceutical Statistics, 17(4), 383-395.
[doi:10.1002/pst.1864](https://doi.org/10.1002/pst.1864)

Takeda, K., Morita, S. and Taguri, M. (2020). TITE-BOIN-ET:
Time-to-event Bayesian optimal interval design to accelerate
dose-finding based on both efficacy and toxicity outcomes.
Pharmaceutical Statistics, 19, 335-349.
[doi:10.1002/pst.1995](https://doi.org/10.1002/pst.1995)

Takeda, K., Zhu, J. and Hirakawa, A. (2025). BF-BOIN-ET: A backfill
Bayesian optimal interval design using efficacy and toxicity outcomes
for dose optimization. Pharmaceutical Statistics, 24, e2470.
[doi:10.1002/pst.2470](https://doi.org/10.1002/pst.2470)

Lin, R., Zhou, Y., Yan, F., Li, D. and Yuan, Y. (2020). BOIN12: Bayesian
optimal interval phase I/II trial design for utility-based dose finding
in immunotherapy and targeted therapies. JCO Precision Oncology, 4,
1393-1402.
[doi:10.1200/PO.20.00257](https://doi.org/10.1200/PO.20.00257)

Zhou, Y., Lin, R., Lee, J. J., Li, D., Wang, L., Li, R. and Yuan, Y.
(2022). TITE-BOIN12: A Bayesian phase I/II trial design to find the
optimal biological dose with late-onset toxicity and efficacy.
Statistics in Medicine, 41(11), 1918-1931.
[doi:10.1002/sim.9337](https://doi.org/10.1002/sim.9337)

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
