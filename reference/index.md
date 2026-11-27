# Package index

## Decision rules

Boundaries and decision tables of the BOIN design.

- [`boin_lambda()`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_lambda.md)
  : Escalation and De-escalation Interval Boundaries
- [`boin_p_tox()`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_p_tox.md)
  : Toxic Threshold Giving a Required De-escalation Boundary
- [`boin_boundary()`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_boundary.md)
  : Integer Decision Boundaries for the BOIN Design
- [`print(`*`<boin_boundary>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.boin_boundary.md)
  : Print BOIN Decision Boundaries
- [`as.data.frame(`*`<boin_boundary>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/as.data.frame.boin_boundary.md)
  : Coerce BOIN Decision Boundaries to a Data Frame
- [`boin_decision_table()`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_decision_table.md)
  : Decision Table for the BOIN Design
- [`print(`*`<boin_decision_table>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.boin_decision_table.md)
  : Print a BOIN Decision Table
- [`plot(`*`<boin_decision_table>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/plot.boin_decision_table.md)
  : Plot a BOIN Decision Table
- [`boin_stopping_table()`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_stopping_table.md)
  : Safety Stopping Boundary as a Two-Row Table

## Single-agent BOIN design

Operating characteristics, simulated trials and selection of the MTD.

- [`sim_boin()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md)
  : Operating Characteristics of a BOIN Design
- [`sim_boin_multi()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_multi.md)
  : Operating Characteristics Across Several Scenarios
- [`boin_simulate()`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_simulate.md)
  : Simulate BOIN Trials
- [`boin_isotonic()`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_isotonic.md)
  : Isotonic Estimates of the DLT Probability
- [`boin_select_mtd()`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_select_mtd.md)
  : Select the MTD from Completed Trials
- [`print(`*`<boin_oc>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.boin_oc.md)
  : Print Operating Characteristics of a BOIN Design
- [`print(`*`<boin_oc_multi>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.boin_oc_multi.md)
  : Print Operating Characteristics Across Scenarios
- [`print(`*`<boin_trials>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.boin_trials.md)
  : Print Simulated BOIN Trials

## Prior information (iBOIN)

Boundaries, simulated trials and operating characteristics of the iBOIN
design, which brings a skeleton and a prior effective sample size into
the BOIN design.

- [`iboin_boundary()`](https://gosukehommaex.github.io/simFastBOIN/reference/iboin_boundary.md)
  : Decision Boundaries of the iBOIN Design
- [`print(`*`<iboin_boundary>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.iboin_boundary.md)
  : Print the Decision Boundaries of the iBOIN Design
- [`sim_iboin()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_iboin.md)
  : Operating Characteristics of the iBOIN Design
- [`sim_iboin_multi()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_iboin_multi.md)
  : Operating Characteristics of the iBOIN Design in Several Scenarios
- [`iboin_simulate()`](https://gosukehommaex.github.io/simFastBOIN/reference/iboin_simulate.md)
  : Simulate iBOIN Trials
- [`iboin_select_mtd()`](https://gosukehommaex.github.io/simFastBOIN/reference/iboin_select_mtd.md)
  : Select the MTD of an iBOIN Trial

## Toxicity grades (gBOIN)

Simulated trials, selection of the MTD and operating characteristics of
the gBOIN design for toxicity categories that carry scores.

- [`sim_gboin()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_gboin.md)
  : Operating Characteristics of the gBOIN Design
- [`sim_gboin_multi()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_gboin_multi.md)
  : Operating Characteristics of the gBOIN Design in Several Scenarios
- [`gboin_simulate()`](https://gosukehommaex.github.io/simFastBOIN/reference/gboin_simulate.md)
  : Simulate gBOIN Trials
- [`gboin_select_mtd()`](https://gosukehommaex.github.io/simFastBOIN/reference/gboin_select_mtd.md)
  : Select the MTD of a gBOIN Trial

## Late-onset toxicity (TITE-BOIN)

- [`tite_boin_decision_table()`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_decision_table.md)
  : Decision Table for the TITE-BOIN Design
- [`print(`*`<tite_boin_decision_table>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.tite_boin_decision_table.md)
  : Print a TITE-BOIN Decision Table
- [`plot(`*`<tite_boin_decision_table>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/plot.tite_boin_decision_table.md)
  : Plot a TITE-BOIN Decision Table
- [`sim_tite_boin()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin.md)
  : Operating Characteristics of a TITE-BOIN Design
- [`sim_tite_boin_multi()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin_multi.md)
  : Operating Characteristics of a TITE-BOIN Design Across Several
  Scenarios
- [`tite_boin_simulate()`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_simulate.md)
  : Simulate TITE-BOIN Trials
- [`print(`*`<tite_boin_oc>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.tite_boin_oc.md)
  : Print Operating Characteristics of a TITE-BOIN Design
- [`print(`*`<tite_boin_oc_multi>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.tite_boin_oc_multi.md)
  : Print Operating Characteristics of a TITE-BOIN Design Across
  Scenarios
- [`print(`*`<tite_boin_trials>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.tite_boin_trials.md)
  : Print Simulated TITE-BOIN Trials

## Backfilling (BF-BOIN and BE-BOIN)

- [`sim_bf_boin()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin.md)
  : Operating Characteristics of a BF-BOIN Design
- [`sim_bf_boin_multi()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin_multi.md)
  : Operating Characteristics of a BF-BOIN Design Across Several
  Scenarios
- [`bf_boin_simulate()`](https://gosukehommaex.github.io/simFastBOIN/reference/bf_boin_simulate.md)
  : Simulate BF-BOIN Trials
- [`sim_be_boin()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_be_boin.md)
  : Operating Characteristics of a BE-BOIN Design
- [`sim_be_boin_multi()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_be_boin_multi.md)
  : Operating Characteristics of a BE-BOIN Design Across Several
  Scenarios
- [`be_boin_simulate()`](https://gosukehommaex.github.io/simFastBOIN/reference/be_boin_simulate.md)
  : Simulate BE-BOIN Trials
- [`print(`*`<backfill_oc>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.backfill_oc.md)
  : Print Operating Characteristics of a BOIN Design with Backfilling
- [`print(`*`<backfill_oc_multi>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.backfill_oc_multi.md)
  : Print Operating Characteristics of a Backfill Design Across
  Scenarios
- [`print(`*`<backfill_trials>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.backfill_trials.md)
  : Print Simulated Trials of a BOIN Design with Backfilling

## Drug combinations

- [`sim_comb_boin()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin.md)
  : Operating Characteristics of the BOIN Drug-Combination Design
- [`sim_comb_boin_multi()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin_multi.md)
  : Operating Characteristics of the BOIN Drug-Combination Design Across
  Scenarios
- [`comb_boin_simulate()`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin_simulate.md)
  : Simulate BOIN Drug-Combination Trials
- [`comb_boin_select_mtd()`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin_select_mtd.md)
  : Select the MTD Combination from Completed Drug-Combination Trials
- [`comb_boin_next()`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin_next.md)
  : Next Dose Combination in a BOIN Drug-Combination Trial
- [`print(`*`<comb_boin_oc>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.comb_boin_oc.md)
  : Print Operating Characteristics of the BOIN Drug-Combination Design
- [`print(`*`<comb_boin_oc_multi>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.comb_boin_oc_multi.md)
  : Print Operating Characteristics of the BOIN Drug-Combination Design
  Across Scenarios
- [`print(`*`<comb_boin_trials>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.comb_boin_trials.md)
  : Print Simulated BOIN Drug-Combination Trials
- [`print(`*`<comb_boin_next>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.comb_boin_next.md)
  : Print the Next Dose Combination of a BOIN Drug-Combination Trial
- [`plot(`*`<comb_boin_next>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/plot.comb_boin_next.md)
  : Plot the Next Dose Combination of a BOIN Drug-Combination Trial

## MTD contour of a drug combination (waterfall)

Simulated trials, selection of the MTD contour and operating
characteristics of the waterfall design.

- [`sim_waterfall()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_waterfall.md)
  : Operating Characteristics of the Waterfall Drug-Combination Design
- [`sim_waterfall_multi()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_waterfall_multi.md)
  : Operating Characteristics of the Waterfall Design Across Scenarios
- [`waterfall_simulate()`](https://gosukehommaex.github.io/simFastBOIN/reference/waterfall_simulate.md)
  : Simulate Waterfall Drug-Combination Trials
- [`waterfall_select_mtd()`](https://gosukehommaex.github.io/simFastBOIN/reference/waterfall_select_mtd.md)
  : Select the MTD Contour from Completed Waterfall Trials
- [`print(`*`<waterfall_oc>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.waterfall_oc.md)
  : Print Operating Characteristics of the Waterfall Design
- [`print(`*`<waterfall_oc_multi>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.waterfall_oc_multi.md)
  : Print Operating Characteristics of the Waterfall Design Across
  Scenarios
- [`print(`*`<waterfall_trials>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.waterfall_trials.md)
  : Print Simulated Waterfall Trials

## Optimal biological dose (BOIN-ET)

Decision rules, simulated trials and selection of the OBD, for BOIN-ET
and its time-to-event (TITE-BOIN-ET) and backfill (BF-BOIN-ET) versions.

- [`boin_et_boundary()`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_et_boundary.md)
  : Boundaries of the BOIN-ET Design
- [`print(`*`<boin_et_boundary>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.boin_et_boundary.md)
  : Print the Boundaries of the BOIN-ET Design
- [`boin_et_decision_table()`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_et_decision_table.md)
  : Decision Table for the BOIN-ET Design
- [`print(`*`<boin_et_decision_table>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.boin_et_decision_table.md)
  : Print a BOIN-ET Decision Table
- [`sim_boin_et()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_et.md)
  : Operating Characteristics of the BOIN-ET Design
- [`sim_boin_et_multi()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_et_multi.md)
  : Operating Characteristics of the BOIN-ET Design Across Scenarios
- [`boin_et_simulate()`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_et_simulate.md)
  : Simulate BOIN-ET Trials
- [`boin_et_select_obd()`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_et_select_obd.md)
  : Select the Optimal Biological Dose of BOIN-ET Trials
- [`print(`*`<boin_et_oc>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.boin_et_oc.md)
  : Print the Operating Characteristics of the BOIN-ET Design
- [`print(`*`<boin_et_oc_multi>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.boin_et_oc_multi.md)
  : Print the Operating Characteristics of BOIN-ET Across Scenarios
- [`print(`*`<boin_et_trials>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.boin_et_trials.md)
  : Print Simulated BOIN-ET Trials
- [`sim_tite_boin_et()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin_et.md)
  : Operating Characteristics of the TITE-BOIN-ET Design
- [`sim_tite_boin_et_multi()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin_et_multi.md)
  : Operating Characteristics of the TITE-BOIN-ET Design Across
  Scenarios
- [`tite_boin_et_simulate()`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_et_simulate.md)
  : Simulate TITE-BOIN-ET Trials
- [`sim_bf_boin_et()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin_et.md)
  : Operating Characteristics of the BF-BOIN-ET Design
- [`sim_bf_boin_et_multi()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin_et_multi.md)
  : Operating Characteristics of the BF-BOIN-ET Design Across Scenarios
- [`print(`*`<bf_boin_et_oc>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.bf_boin_et_oc.md)
  : Print the Operating Characteristics of the BF-BOIN-ET Design
- [`print(`*`<bf_boin_et_oc_multi>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.bf_boin_et_oc_multi.md)
  : Print the Operating Characteristics of BF-BOIN-ET Across Scenarios

## Optimal biological dose from grades (gBOIN-ET and TITE-gBOIN-ET)

Simulated trials, selection of the OBD and operating characteristics of
the gBOIN-ET design for toxicity and efficacy categories that carry
scores, and of its time-to-event version TITE-gBOIN-ET.

- [`sim_gboin_et()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_gboin_et.md)
  : Operating Characteristics of the gBOIN-ET Design
- [`sim_gboin_et_multi()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_gboin_et_multi.md)
  : Operating Characteristics of the gBOIN-ET Design in Several
  Scenarios
- [`gboin_et_simulate()`](https://gosukehommaex.github.io/simFastBOIN/reference/gboin_et_simulate.md)
  : Simulate gBOIN-ET Trials
- [`gboin_et_select_obd()`](https://gosukehommaex.github.io/simFastBOIN/reference/gboin_et_select_obd.md)
  : Select the Optimal Biological Dose of gBOIN-ET Trials
- [`sim_tite_gboin_et()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_gboin_et.md)
  : Operating Characteristics of the TITE-gBOIN-ET Design
- [`sim_tite_gboin_et_multi()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_gboin_et_multi.md)
  : Operating Characteristics of the TITE-gBOIN-ET Design in Several
  Scenarios
- [`tite_gboin_et_simulate()`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_gboin_et_simulate.md)
  : Simulate TITE-gBOIN-ET Trials

## Optimal biological dose (BOIN12 and TITE-BOIN12)

Desirability scores, simulated trials and selection of the OBD of the
BOIN12 design and of its time-to-event version, TITE-BOIN12.

- [`boin12_decision_table()`](https://gosukehommaex.github.io/simFastBOIN/reference/boin12_decision_table.md)
  : Desirability Score Table of the BOIN12 Design
- [`print(`*`<boin12_decision_table>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.boin12_decision_table.md)
  : Print the Desirability Score Table of the BOIN12 Design
- [`sim_boin12()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin12.md)
  : Operating Characteristics of the BOIN12 Design
- [`sim_boin12_multi()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin12_multi.md)
  : Operating Characteristics of the BOIN12 Design Across Scenarios
- [`boin12_simulate()`](https://gosukehommaex.github.io/simFastBOIN/reference/boin12_simulate.md)
  : Simulate BOIN12 Trials
- [`boin12_select_obd()`](https://gosukehommaex.github.io/simFastBOIN/reference/boin12_select_obd.md)
  : Select the Optimal Biological Dose of BOIN12 Trials
- [`print(`*`<boin12_oc>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.boin12_oc.md)
  : Print the Operating Characteristics of the BOIN12 Design
- [`print(`*`<boin12_oc_multi>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.boin12_oc_multi.md)
  : Print the Operating Characteristics of the BOIN12 Design Across
  Scenarios
- [`print(`*`<boin12_trials>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.boin12_trials.md)
  : Print Simulated BOIN12 Trials
- [`sim_tite_boin12()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin12.md)
  : Operating Characteristics of the TITE-BOIN12 Design
- [`sim_tite_boin12_multi()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin12_multi.md)
  : Operating Characteristics of the TITE-BOIN12 Design Across Scenarios
- [`tite_boin12_simulate()`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin12_simulate.md)
  : Simulate TITE-BOIN12 Trials

## Optimal biological dose (U-BOIN)

Utility tables, simulated trials, selection of the OBD and operating
characteristics of the U-BOIN design.

- [`uboin_utility_table()`](https://gosukehommaex.github.io/simFastBOIN/reference/uboin_utility_table.md)
  : Utility Table of the U-BOIN Design
- [`sim_uboin()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_uboin.md)
  : Operating Characteristics of the U-BOIN Design
- [`sim_uboin_multi()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_uboin_multi.md)
  : Operating Characteristics of the U-BOIN Design Across Scenarios
- [`uboin_simulate()`](https://gosukehommaex.github.io/simFastBOIN/reference/uboin_simulate.md)
  : Simulate U-BOIN Trials
- [`uboin_select_obd()`](https://gosukehommaex.github.io/simFastBOIN/reference/uboin_select_obd.md)
  : Select the OBD from Completed U-BOIN Trials
- [`print(`*`<uboin_oc>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.uboin_oc.md)
  : Print Operating Characteristics of the U-BOIN Design
- [`print(`*`<uboin_oc_multi>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.uboin_oc_multi.md)
  : Print Operating Characteristics of the U-BOIN Design Across
  Scenarios
- [`print(`*`<uboin_trials>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.uboin_trials.md)
  : Print Simulated U-BOIN Trials

## Optimal biological dose combinations (Comb-BOIN12)

Simulated trials and selection of the optimal biological dose
combination of the Comb-BOIN12 and TITE-Comb-BOIN12 designs for two
drugs given together.

- [`sim_comb_boin12()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin12.md)
  : Operating Characteristics of the Comb-BOIN12 Design
- [`sim_comb_boin12_multi()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin12_multi.md)
  : Operating Characteristics of the Comb-BOIN12 Design Across Scenarios
- [`comb_boin12_simulate()`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin12_simulate.md)
  : Simulate Comb-BOIN12 Trials
- [`comb_boin12_select_obd()`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin12_select_obd.md)
  : Select the OBD Combination of Comb-BOIN12 Trials
- [`sim_tite_comb_boin12()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_comb_boin12.md)
  : Operating Characteristics of the TITE-Comb-BOIN12 Design
- [`sim_tite_comb_boin12_multi()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_comb_boin12_multi.md)
  : Operating Characteristics of the TITE-Comb-BOIN12 Design Across
  Scenarios
- [`tite_comb_boin12_simulate()`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_comb_boin12_simulate.md)
  : Simulate Trials of the TITE-Comb-BOIN12 Design
- [`print(`*`<comb_boin12_oc>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.comb_boin12_oc.md)
  : Print the Operating Characteristics of the Comb-BOIN12 Design
- [`print(`*`<comb_boin12_oc_multi>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.comb_boin12_oc_multi.md)
  : Print the Operating Characteristics of Comb-BOIN12 Across Scenarios
- [`print(`*`<comb_boin12_trials>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.comb_boin12_trials.md)
  : Print Simulated Comb-BOIN12 Trials

## Optimal dose combinations (BOIN-ETC)

Simulated trials and selection of the optimal dose combination of the
BOIN-ETC designs, which run BOIN-ET in subtrials of the dose matrix.

- [`sim_boin_etc()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_etc.md)
  : Operating Characteristics of the BOIN-ETC Design
- [`sim_boin_etc_multi()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_etc_multi.md)
  : Operating Characteristics of the BOIN-ETC Design Across Scenarios
- [`boin_etc_simulate()`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_etc_simulate.md)
  : Simulate BOIN-ETC Trials
- [`boin_etc_select_obd()`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_etc_select_obd.md)
  : Select the Optimal Dose Combination of BOIN-ETC Trials
- [`print(`*`<boin_etc_oc>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.boin_etc_oc.md)
  : Print the Operating Characteristics of the BOIN-ETC Design
- [`print(`*`<boin_etc_oc_multi>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.boin_etc_oc_multi.md)
  : Print the Operating Characteristics of the BOIN-ETC Design Across
  Scenarios
- [`print(`*`<boin_etc_trials>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.boin_etc_trials.md)
  : Print Simulated BOIN-ETC Trials

## The 3+3 design

- [`oc_3p3()`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_3p3.md)
  : Exact Operating Characteristics of the 3+3 Design
- [`print(`*`<oc_3p3>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.oc_3p3.md)
  : Print Operating Characteristics of the 3+3 Design
- [`sim_3p3()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_3p3.md)
  : Simulate the 3+3 Design
- [`decision_table_3p3()`](https://gosukehommaex.github.io/simFastBOIN/reference/decision_table_3p3.md)
  : Decision Table for the 3+3 Design
- [`print(`*`<decision_table_3p3>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.decision_table_3p3.md)
  : Print a 3+3 Decision Table
- [`plot(`*`<decision_table_3p3>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/plot.decision_table_3p3.md)
  : Plot a 3+3 Decision Table

## Tables and figures

True scenarios, and the results of several designs set side by side and
judged against the true doses.

- [`oc_scenarios()`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_scenarios.md)
  : True Dose-Toxicity and Dose-Response Scenarios
- [`print(`*`<oc_scenarios>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.oc_scenarios.md)
  : Print True Scenarios
- [`plot(`*`<oc_scenarios>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/plot.oc_scenarios.md)
  : Plot True Scenarios
- [`oc_results()`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_results.md)
  : Operating Characteristics of Several Designs in One Table
- [`print(`*`<oc_results>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.oc_results.md)
  : Print the Operating Characteristics of Several Designs
- [`plot(`*`<oc_results>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/plot.oc_results.md)
  : Plot the Operating Characteristics of Several Designs
- [`oc_compare()`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_compare.md)
  : Compare Several Designs Against the True Doses
- [`print(`*`<oc_compare>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/print.oc_compare.md)
  : Print a Comparison of Designs
- [`plot(`*`<oc_compare>`*`)`](https://gosukehommaex.github.io/simFastBOIN/reference/plot.oc_compare.md)
  : Plot a Comparison of Designs
