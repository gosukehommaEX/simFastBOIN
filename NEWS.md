# simFastBOIN 2.3.0

This version adds designs to the single-agent BOIN design: the
iBOIN design of Zhou et al. (2021), which brings prior information on the DLT
probability of every dose into the decision boundaries, the gBOIN design of
Mu et al. (2019), which decides from toxicity grades instead of from the DLT
alone, the gBOIN-ET design of Takeda et al. (2022), which finds the
optimal biological dose from graded toxicity and graded efficacy, and its
time-to-event version, the TITE-gBOIN-ET design of Takeda et al. (2023). For
two drugs given together, it adds the waterfall design of Zhang and Yuan
(2016), which finds the MTD contour, and for the optimal biological dose the
U-BOIN design of Zhou, Lee and Yuan (2019), which weighs the response
against the toxicity by the utility of the joint outcome of a patient.
Existing functions keep their results, except that BOIN-ET, TITE-BOIN-ET
and BF-BOIN-ET now compare the utilities of the doses after rounding, so
that a tie between two doses no longer depends on rounding error (see "Doses
with the same utility" below). The pkgdown site has a cheat sheet and
articles on how the package is checked and on its computation time.

## iBOIN

* `iboin_boundary()` tabulates the escalation and de-escalation boundaries of
  the iBOIN design from a skeleton, the prior DLT probability of every dose,
  and a prior effective sample size, by equations 3 and 5 of the article. The
  boundaries depend on the dose and on the number of patients treated there,
  and reproduce Table 1 of the article. Where the closed form puts the
  escalation boundary above the de-escalation boundary, the integer boundaries
  with the smallest probability of an incorrect decision are searched. The
  robust prior of the article is available (`robust = TRUE`). It has a
  `print()` method.

* `iboin_simulate()` simulates iBOIN trials in C++ (`src/iboin_core.h`). The
  trials run as BOIN trials with the boundaries of `iboin_boundary()`; dose
  elimination and the safety stopping rule keep the uniform prior of the BOIN
  design. With a prior effective sample size of 0 the trials are, for the same
  seed, those of `boin_simulate()`.

* `iboin_select_mtd()` selects the MTD with the prior information in the
  estimates of the DLT probabilities, `(n_tox + prior_ess * skeleton) /
  (n_pts + prior_ess)`, or from the data alone. Only doses given to patients
  can be selected, and the elimination uses the data alone.

* `sim_iboin()` and `sim_iboin_multi()` give the operating characteristics in
  one and in several scenarios. Their results print like those of `sim_boin()`
  and enter `oc_results()` and `oc_compare()` under the name iBOIN. The MTD is
  selected with the prior information by default (`select_with_prior`): the
  article gives both selections without saying which one its simulation study
  uses, and its operating characteristics are reproduced with the prior
  information and not without it.

* `inst/reproduce/reproduce_Zhou2021.R` reproduces the 100 boundaries of
  Table 1 of the article exactly and compares the rows BOIN, iBOIN and iBOINR
  of Table S.1 of its supporting information (10 scenarios, 5 quantities
  each) with simulations: 144 of the 150 values lie within simulation error.
  The other six are the percentage of patients above the MTD in Scenarios 2
  and 10, where the published values agree with the percentage of patients
  treated at doses whose true DLT probability exceeds 0.40; they are
  registered in `inst/validation/explained.csv`.

* The C++ core was compared trial by trial with an implementation of the
  rules of the article in Python on the same random numbers (18 settings of
  3,000 trials, no difference), and without prior information with the BOIN
  core.

## gBOIN

* `gboin_simulate()` simulates gBOIN trials in C++ (`src/gboin_core.h`) for a
  toxicity outcome with several categories that each carry a score, such as
  the equivalent toxicity score of Yuan et al. (2007). The scores are divided
  by their maximum and their total at a dose takes the place of the number of
  DLTs: the design escalates when the mean score is at most `lambda_e` and
  de-escalates when it is at least `lambda_d`, the BOIN boundaries for the
  target. A dose is eliminated when the posterior probability that its mean
  score exceeds the target, evaluated at the total score, is above the cutoff
  (`elim_rule = "posterior"`, the rule of the article), or when the total
  score reaches the elimination boundary of the number of DLTs
  (`elim_rule = "count"`, the rule of the 'UnifiedDoseFinding' package). One
  random number is drawn per patient, so that two categories of scores 0 and
  1 give, for the same seed, the trials of `boin_simulate()`. The total
  normalized score of a dose is the total of the scores as given, divided
  once by the largest score, as in `gboin_et_simulate()`: with scores that
  are binary fractions, such as 0, 0.5, 1 and 1.5, the total is exact and
  the decisions are the same on every platform.

* `gboin_select_mtd()` selects the MTD from the total scores by the isotonic
  regression of the BOIN design.

* `sim_gboin()` and `sim_gboin_multi()` give the operating characteristics in
  one and in several scenarios, with the average number of patients of each
  toxicity category at each dose. Their results print with the labels of a
  score and enter `oc_results()` and `oc_compare()` under the name gBOIN.

* The continuous end points of the article (normal distribution, total
  toxicity burden) are not implemented.

* The C++ core was compared trial by trial with an implementation of the
  rules in Python on the same random numbers (17 settings of 3,000 trials
  under each elimination rule, no difference), and with two categories with
  the BOIN core.

* `inst/reproduce/reproduce_Mu2019.R` reproduces the 14 boundaries of Table 2
  of the article and compares the operating characteristics of its Table 5
  (10 scenarios, the selection percentage and the number of patients of six
  doses) with simulations under `elim_rule = "posterior"`: 101 of the 120
  cells lie within simulation error and three hold misprints, registered in
  `inst/validation/explained.csv`. The other 16 cells, seven of them in
  Scenario 10, do not agree, and the reason was not found; they are left as
  they are. Of these, 13 lie outside the allowance also for the values of
  the 'UnifiedDoseFinding' package, which the first author of the article
  co-wrote (with 4,000 trials on both sides).

* `inst/validation/compare-with-UnifiedDoseFinding.R` compares the selection
  percentages and the numbers of patients of the same ten scenarios with the
  'UnifiedDoseFinding' package: under `elim_rule = "count"`, the rule of that
  package, all 130 values agree within simulation error.

## gBOIN-ET

* `gboin_et_simulate()` simulates gBOIN-ET trials in C++
  (`src/gboinet_core.h`) for a toxicity outcome and an efficacy outcome with
  several categories that each carry a score, such as the equivalent toxicity
  score and the equivalent efficacy score of Takeda et al. (2022). The scores
  are divided by their maximum, and the mean scores of the patients at a dose
  take the place of the observed DLT and response rates in the rules of
  `boin_et_simulate()`, with posterior probabilities from the total scores.
  The toxicity event and the efficacy event of a patient are linked by a
  Gaussian copula and the category of a patient with an event is drawn
  independently, the joint distribution of `boinet::gboinet()`;
  `boinet_compatible = TRUE` gives the rules of that function. With two
  categories of scores 0 and 1 for both outcomes the trials are, for the
  same seed, those of `boin_et_simulate()`.

* The boundaries are those of `boin_et_boundary()` for the normalized
  targets, or the three values given through `lambda1`, `lambda2` and
  `eta1`, for example those printed in an article: the mean scores are not
  multiples of `1 / n`, so that boundaries rounded to two decimals can give
  other decisions.

* By default a dose is excluded by the posterior probability rules of
  'boinet': for its toxicity as in Section 2.2 of Takeda et al. (2022), and
  for its efficacy when the posterior probability that its mean efficacy
  score exceeds `0.6 * target_resp` is below 0.01. The efficacy rule of the
  article, a posterior probability above 0.50 that the mean efficacy score
  is below half the target, is given by `p_eli_resp = target_resp / 2` and
  `cutoff_eli_resp = 0.50`.

* `gboin_et_select_obd()` selects the OBD from the total scores with the
  rules of `boin_et_select_obd()`. The mean efficacy scores are the observed
  ones, those of the linear fractional polynomial of Takeda et al. (2022)
  (`est_method = "fp_linear"`, fitted in C++), or those of the logistic
  fractional polynomial of 'boinet'. `min_obd_sample` restricts the selection
  to doses with a minimum number of patients, six in the simulations of the
  article. The estimated efficacy scores and the utilities are compared
  after rounding to ten decimals, so that doses that are equal as fractions
  tie on every platform; the mean scores during the trial are computed
  exactly for scores that are binary fractions.

* `sim_gboin_et()` and `sim_gboin_et_multi()` give the operating
  characteristics in one and in several scenarios, with the average number
  of patients of each category at each dose. Their results print with the
  labels of scores and enter `oc_results()` and `oc_compare()` under the
  name gBOIN-ET.

* The C++ core was compared trial by trial with an implementation of the
  rules in Python on the same random numbers (18 settings of 3,000 trials,
  no difference), and with two categories with the BOIN-ET core.

* `inst/reproduce/reproduce_Takeda2022.R` reproduces the 12 boundaries
  printed in the article and compares the gBOIN-ET rows of its Table 1 (10
  scenarios, the percentage of early terminations, the selection percentage
  and the number of patients of six doses). The SAS program of the authors
  does not draw the categories from the probabilities of their Table S1 but
  one Bernoulli variate per category, the highest category that comes up
  being taken; the script simulates the probabilities that the program
  draws, and lists the other points in which the program differs from the
  article. Of the 130 cells, 119 lie within simulation error. Four are
  explained by rules of the program that the package does not have and are
  registered in `inst/validation/explained.csv`: a transcription of the
  program reproduces them. The other 7 cells, six of them numbers of
  patients, do not agree, also not with the transcription, and the reason
  was not found; they are left as they are.

* `inst/validation/compare-with-gboinet.R` compares the selection
  percentages and the numbers of patients of six of these scenarios with
  `boinet::gboinet()`, under four selection rules: with
  `boinet_compatible = TRUE` all 156 values agree within simulation error.

## TITE-gBOIN-ET

* `tite_gboin_et_simulate()` simulates TITE-gBOIN-ET trials in C++
  (`src/tite_gboinet_core.h`): the gBOIN-ET design of `gboin_et_simulate()`
  for outcomes that are observed late, in which the next cohort is enrolled
  while earlier patients are still assessed (Takeda et al., 2023). The
  patients arrive over time as in `tite_boin_et_simulate()`, with the same
  arguments for the windows, the arrivals, the event times and the share of
  patients that must have finished before a decision. The mean score of a
  dose is the total of the scores observed so far divided by the effective
  number of patients, a patient still under assessment counting by the
  fraction of the window followed so far, and the rules of gBOIN-ET are
  applied to these means; the admissibility of a dose uses the scores
  observed so far and the number of patients treated. These are the rules
  of `boinet::tite.gboinet()`; `boinet_compatible = TRUE` also drops the two
  changes of the package described in `boin_et_simulate()`. The boundaries
  can be given through `lambda1`, `lambda2` and `eta1`.

* With two categories of scores 0 and 1 for both outcomes the trials are,
  for the same seed, those of `tite_boin_et_simulate()`, and with arrivals a
  fixed time apart that is longer than both windows they are those of
  `gboin_et_simulate()`.

* `sim_tite_gboin_et()` and `sim_tite_gboin_et_multi()` give the operating
  characteristics in one and in several scenarios, with the average trial
  duration. The OBD is selected by `gboin_et_select_obd()` from the data at
  the last decision; the default, the highest observed mean efficacy score
  among the doses up to the MTD, is the selection of the article. The
  results print with the labels of scores and enter `oc_results()` and
  `oc_compare()` under the name TITE-gBOIN-ET.

* The C++ core was compared trial by trial with an implementation of the
  rules in Python on the same random numbers (22 settings of 3,000 trials,
  no difference), with two categories with the TITE-BOIN-ET core (4
  settings) and with arrivals far apart with the gBOIN-ET core (3 settings).
  The counts, the stopping reasons and the selected doses are the same
  whether or not the compiler fuses a multiplication and an addition into
  one operation.

* `inst/reproduce/reproduce_Takeda2023.R` compares the TITE-gBOIN-ET cells
  of Table 1 and of Tables S2 to S6 of the article and its supporting
  information (10 scenarios under three settings of the efficacy window and
  of the accrual, 570 cells) with simulations. The R program of the authors
  was written again in C++ with R's random numbers, and the transcription
  reproduces all 4,560 summaries that the supporting information holds for
  the program exactly. The program differs from the article and from
  `boinet::tite.gboinet()` in two points: where a single dose has the
  highest mean efficacy score among the current dose and its neighbors, it
  passes that dose to `sample()`, which then draws from all doses up to it,
  and it excludes doses with the effective numbers of patients. The first
  point moves the published tables substantially (in Scenario 1 the
  percentage of correct selections is 64.3 as published and about 75
  without it). The package has neither. The cells that these two points
  move out of the simulation error are registered in
  `inst/validation/explained.csv`, and the script also compares the
  simulations with the transcription of the program with the two points
  removed (`inst/reproduce/takeda2023-tables.csv`). Of the 570 cells, 273
  lie within simulation error of the published tables and 297 are explained
  by the two points; against the program without them, all 570 lie within
  simulation error.

* `inst/validation/compare-with-tite-gboinet.R` compares the selection
  percentages, the numbers of patients and the trial duration of six of
  these scenarios with `boinet::tite.gboinet()`, under four selection rules.
  With `boinet_compatible = TRUE`, all 168 values agree with that function
  within simulation error.

## Waterfall

* `waterfall_simulate()` simulates trials of the waterfall design of Zhang
  and Yuan (2016) in C++ (`src/waterfall_core.h`). The design looks for the
  MTD contour of two drugs given together, one MTD at every dose level of
  drug A, by dividing the dose matrix into subtrials that are run with the
  single-agent BOIN design: the first subtrial goes up the lowest level of
  drug B and then along the highest level of drug A, and each following
  subtrial runs along the row below the last candidate MTD, starting one
  level of drug B above it. The rules are those of `BOIN::get.oc.comb()`
  with `mtd.contour = TRUE` in 'BOIN' 2.7.2, the implementation by the
  authors of the article, with the same order of the random numbers. Four
  rules of that function are not in the article and are described in the
  help page: a candidate in the first column below the highest level of drug
  A removes the rows above it and can add a subtrial on its own row, the
  combinations to the right of a candidate in its row are excluded from the
  selection, the trial ends when a subtrial ends without a candidate, and
  the MTD of a row that is not to the right of the MTD of the row above is
  moved to the column of that MTD.

* When a subtrial ends without a candidate MTD, `BOIN::get.oc.comb()` drops
  its data, so that its patients are missing from the counts and a trial
  stopped in the first subtrial has no patients. The package keeps these
  data by default; `boin_compatible = TRUE` drops them, and then the trials
  are, for the same seed, those of `BOIN::get.oc.comb()`. That function also
  takes the escalation and de-escalation boundaries from the default values
  of `p.saf` and `p.tox` whatever values are given; the package uses the
  values given.

* `waterfall_select_mtd()` selects the MTD of every level of drug A from the
  data of a trial and the combinations excluded during it, with the
  bivariate isotonic regression and the rules of `BOIN::get.oc.comb()`.

* `sim_waterfall()` and `sim_waterfall_multi()` give the operating
  characteristics in one and in several scenarios: the percentage of trials
  selecting each combination as an MTD, the percentage selecting exactly
  the true MTD contour, the percentages of patients treated on, above and
  below the contour, and the average numbers of patients and DLTs. They
  have `print()` methods. The phase I/II extension of the article is not
  provided.

* `BOIN::get.oc.comb(mtd.contour = TRUE)` was written again in Python with
  R's random numbers. With the default seed of that function and 1,000
  trials, the transcription gives all 398 cells that the article and its
  supporting information print for the waterfall design (Tables IV, V and
  S1, 14 scenarios) exactly, so the article was produced with these rules.
  The C++ core was compared trial by trial with the transcription on the
  same random numbers (41 settings with and without `boin_compatible`,
  164,000 trials, no difference). The trials and the selections are the same
  whether or not the compiler fuses a multiplication and an addition into
  one operation.

* `inst/reproduce/reproduce_Zhang2016.R` compares the 398 cells of the
  article with `sim_waterfall()` under `boin_compatible = TRUE`, the seed 6
  and 1,000 trials, and judges them by equality with the printed digits,
  without an allowance for simulation error
  (`inst/reproduce/zhang2016-tables.csv` holds the cells). All 398 cells
  are equal. The sample sizes of Table V, which the article matched between
  its two designs, are in all 14 scenarios the multiples of 3 closest to the
  average number of patients treated, those of the dropped subtrials
  included.

* `inst/validation/compare-waterfall-with-BOIN.R` compares every cell of the
  output of `BOIN::get.oc.comb(mtd.contour = TRUE)` with `sim_waterfall()`
  in 20 settings (seven dose matrices, cohorts of 1, 2 and 3, titration,
  the extra safety rule and the bound on the MTD), again by equality. All
  870 values are equal.

## U-BOIN

* `uboin_simulate()` simulates trials of the U-BOIN design of Zhou, Lee and
  Yuan (2019) in C++ (`src/uboin_core.h`), for a binary toxicity and a
  binary efficacy outcome observed before the next cohort and the
  pick-the-winner allocation of the article. Stage I runs as the BOIN design
  with the target `tox_upper - 0.05`; once a dose has `n_stage1` patients,
  the next cohorts escalate above the highest dose tried while its DLT rate
  is at most the escalation boundary, and otherwise go to the admissible
  dose with the largest posterior mean utility. The outcomes of a patient
  are drawn from the Gumbel model of the article (`resp_assoc`) or from
  given probabilities (`p_outcome`).

* The article does not give the priors of its toxicity and futility rules,
  nor whether doses not yet tried can be assigned in stage II. Before the
  code was written, the design was written in Python and its simulations
  compared with Table 3 of the article and Tables S8 and S9 of its
  supporting information: Beta(1, 1) priors (`prior_tox`, `prior_resp`) and
  the tried doses only agree with them, and are the defaults. The decision
  tables of the supporting information (Tables S4 to S6) are reproduced with
  `prior_resp = c(0.05, 0.05)`. The results of Scenario 7 of Table 3 agree
  with the rates of Scenario 7 of Table S8 rather than with those printed,
  and Table S8 agrees with s1 = 12 rather than with the s1 = 9 of its text.

* `uboin_utility_table()` tabulates the posterior mean utility and the
  admissibility of a dose for every outcome of its patients, the decision
  tables of the article. `uboin_select_obd()` selects the OBD from the data
  of completed trials. `sim_uboin()` and `sim_uboin_multi()` give the
  operating characteristics in one and in several scenarios, with `print()`
  methods, and their results enter `oc_results()` and `oc_compare()` under
  the name U-BOIN. The adaptive and equal randomization, the version for
  efficacy observed late (Section 2.5 of the article) and outcomes with more
  than two levels, such as the three levels of toxicity and of efficacy of
  Example 3 of the article, are not provided.

* The C++ core was compared trial by trial with the implementation in
  Python on the same random numbers (36 settings of 3,000 trials, no
  difference), and the trials and the selections are the same whether or
  not the compiler fuses a multiplication and an addition into one
  operation.

* `inst/reproduce/reproduce_Zhou2019.R` compares the boundaries of Table 2
  and Table S2, the utility tables S4 to S6, the true utilities and the 264
  simulated cells of Table 3, Table S9 and Table S8 with the package
  (`inst/reproduce/zhou2019-tables.csv` and
  `inst/reproduce/zhou2019-utility-tables.csv` hold the published values).
  All 264 simulated cells lie within simulation error (10,000 trials; the
  early stopping of Scenario 3 of Table S8, 2.7 printed and 1.33 simulated,
  lies close to its allowance of 1.43). The 72 boundaries of Table S2, the
  163 utilities of Tables S4 to S6 and the 120 true utilities are equal to
  the printed values, and so are 11 of the 12 boundaries of Table 2: the
  article prints 0.358 for lambda_d at the target 0.30, which is 0.35852.
  Tables S4 to S6 judge the admissibility of all 323 counts as
  `prior_resp = c(0.05, 0.05)` does; 12 of them differ under the default
  prior. With the rates printed for Scenario 7 of Table 3, 7 of its 11
  cells lie outside simulation error, and with s1 = 9 in Table S8, 16 of
  its 88 cells.

## Doses with the same utility in BOIN-ET, TITE-BOIN-ET and BF-BOIN-ET

* `boin_et_select_obd()`, and with it `sim_boin_et()`, `sim_tite_boin_et()`
  and their versions for several scenarios, compare the response estimates
  and the utilities of the doses after rounding to ten decimals; ties go to
  the lowest dose, as before. Two doses can have the same utility as
  fractions: with the scoring utility and its default scores, 1 DLT and 7
  responses in 15 patients and 4 DLTs and 9 responses in 15 patients both
  give 196 / 3. Version 2.2.0 compared the utilities as computed, as
  'boinet' does, so that rounding error decided between such doses, and the
  selected dose could differ between platforms (for example where the
  compiler fuses a multiplication and an addition into one operation). In
  these trials the OBD can now differ from that of 2.2.0 and from that of
  `boinet::obd.select()`, also with `boinet_compatible = TRUE`; the simulated
  trials themselves are unchanged. In 96,000 simulated trials (four
  scenarios, the four selection rules, with and without
  `boinet_compatible`) the selection of 2.2.0 differed from a selection in
  exact fractions in fewer than one trial in a thousand, and in some trials
  between a compilation with and one without the fused operation; the
  present selection agrees with the exact fractions in every trial and is
  the same for both compilations.

* `sim_bf_boin_et()` rounds its utility in the same way, for the backfill
  dose and for the OBD. When the backfill dose changes, the rest of the
  trial changes, and the trials of one simulation share one stream of random
  numbers, so that for a given seed the results differ from those of 2.2.0
  from the first such trial on; they agree within simulation error. In
  40,000 trials simulated one by one, a few trials in a thousand changed,
  and the present core gave the same trials with and without the fused
  operation.

## Other changes

* The vignettes describe the new designs: iBOIN and gBOIN in
  `vignette("simFastBOIN-introduction")`, gBOIN-ET and TITE-gBOIN-ET in
  `vignette("simFastBOIN-boin-et")`, the waterfall design in
  `vignette("simFastBOIN-combination")` and U-BOIN in
  `vignette("simFastBOIN-boin12")`, each with an example, the points in which
  the rules differ from the article or from the reference implementation,
  and the results of its reproduction and comparison scripts, which the
  vignettes read from the result files. The README has a section for each
  of them.

* A new vignette, `vignette("simFastBOIN-comparison")`, sets designs side by
  side with `oc_scenarios()`, `oc_results()` and `oc_compare()`: designs
  that select the MTD (BOIN, iBOIN, TITE-BOIN, BF-BOIN and 3+3), the two
  elimination rules of gBOIN and gBOIN-ET with TITE-gBOIN-ET on their mean
  scores, designs that select the OBD (BOIN-ET, BOIN12 and U-BOIN) and
  drug-combination designs (BOIN combination, Comb-BOIN12 and BOIN-ETC).

* `inst/reproduce/reproduce_Takeda2018.R` and
  `inst/reproduce/reproduce_Takeda2020.R` now simulate scenario 11 with a
  response rate of 0.05 at the lowest dose, as Figure 2 of both articles
  plots and labels it, in place of the 0.55 that Table 3 (2018) and Table 2
  (2020) print: the published results agree with 0.05, and all the cells of
  the scenario lie within simulation error. The results with 0.55 are kept
  as INFO rows. The cells that fail go from 14 to 4 (2018) and from 29 to 12
  (2020).

* The pkgdown site has an article, "How the package is checked", that
  summarizes the judgments of every script in `inst/reproduce/` and
  `inst/validation/`, the reasons registered in
  `inst/validation/explained.csv` and the values that fail, all read from
  the results files of the scripts. The article is not part of the package.

* The pkgdown site has an article, "Computation time", with the time that
  every design takes on one machine, measured by `inst/benchmark/benchmark.R`
  (`run_benchmark()`): the designs that were compared with a reference
  implementation, in another R package or in the code of the authors, are
  run by both with the same settings and the same number of trials. The
  results are in `inst/benchmark/results/benchmark.csv`. The section "Speed"
  of the README summarizes them.

* The pkgdown site has a cheat sheet, a single page of cards in two columns
  that prints on A4 in landscape, reached by "Cheat Sheet" in the
  navigation bar: the designs and their functions, the decision rules, the
  common arguments, examples for every group of designs, the comparison of
  designs and the conduct of a trial. `pkgdown/cheatsheet_figures.R` runs
  its code and draws its figures.

* The print methods of `sim_boin()`, `sim_boin_multi()` and `boin_simulate()`
  take the name of the design and the labels of the rows from the result, so
  that the results of iBOIN and gBOIN print through them. The output for the
  BOIN design is unchanged.

* The print methods of `sim_boin_et()`, `sim_boin_et_multi()` and
  `boin_et_simulate()` do the same for the results of gBOIN-ET and of
  TITE-gBOIN-ET. The output for the BOIN-ET design is unchanged.

# simFastBOIN 2.2.0

This version adds designs that decide from both the DLT and the response and
select the optimal biological dose (OBD): the BOIN-ET design of Takeda, Taguri
and Morita (2018) with its time-to-event (TITE-BOIN-ET) and backfill
(BF-BOIN-ET) versions, the utility-based BOIN12 design of Lin et al. (2020)
with its time-to-event version TITE-BOIN12 (Zhou et al., 2022), and, for two
drugs given together, Comb-BOIN12 and TITE-Comb-BOIN12 (Lu et al., 2025) and
BOIN-ETC (Kakizume et al., 2024). It also adds the accelerated BOIN design of
Kojima et al. (2024) and functions that tabulate and plot the results of
different designs side by side. Existing functions only gain arguments whose
defaults keep their results, so no existing result changes.

## BOIN-ET

* `boin_et_boundary()` finds the boundaries `lambda1`, `lambda2` and `eta1` by
  the grid search of `boinet::gridoptim()`, and reproduces the values 0.14,
  0.35 and 0.48 of the article. It has a `print()` method.

* `boin_et_decision_table()` tabulates the decision for every number of DLTs
  and responses, marking the counts that make a dose inadmissible, and
  reproduces Table 2 of the article. It has a `print()` method.

* `boin_et_simulate()` simulates BOIN-ET trials in C++ (`src/boinet_core.h`).
  The rules follow `boinet::boinet()`: the admissibility of every dose is
  re-derived after each cohort, an inadmissible dose is skipped, and the trial
  stops once a dose reaches `n_earlystop` patients (by default the maximum
  sample size). Two rules differ, to keep to the article, which does not allow
  dose skipping during escalation: a dose excluded by its DLT data also
  excludes every higher dose, as in Takeda et al. (2025) and the 'bfboinet'
  package, so that an escalation cannot skip it, and an escalation from the
  highest dose ends the trial when that dose is excluded (stop reason
  `"top_dose_excluded"`). `boinet_compatible = TRUE` applies the rules of
  'boinet' 1.6.0 exactly, in this function and in all the BOIN-ET functions
  below. The DLT and the response of a patient are linked by a Gaussian
  copula, as in 'boinet', or by the correlation of the binary outcomes, as in
  the simulations of the article (`resp_cor_type`). `p_eli` and `p_eli_resp`
  set the rates of the rules that make a dose inadmissible, so that the rules
  of the article can be used as well as those of 'boinet'.

* `boin_et_select_obd()` selects the OBD from completed trials with the
  estimators (`est_method`: observed rates, model-averaged unimodal isotonic
  regression, or fractional polynomial logistic regression with the 'mfp'
  package) and selection rules (`obd_method`: the highest response up to the
  estimated MTD, or a weighted, truncated linear or scoring utility) of
  `boinet::obd.select()`. The selection runs in C++; the model-based estimates
  are computed once for every distinct set of final data.

* `sim_boin_et()` and `sim_boin_et_multi()` report the percentage of trials
  selecting each dose as the OBD and the average numbers of patients, DLTs
  and responses at each dose, for one or several scenarios, with `print()`
  methods.

* The C++ engine and the OBD selection agree trial by trial with a separate
  implementation written in Python, fed with the same random numbers. The
  tests compare the boundaries with `boinet::gridoptim()` and the OBD of
  simulated trials with `boinet::obd.select()` exactly.
  `inst/reproduce/reproduce_Takeda2018.R` reproduces the boundaries, Table 2
  and Table 3 of the article, and `inst/validation/compare-with-boinet.R`
  compares the operating characteristics with `boinet::boinet()`, which draws
  its random numbers differently.

* 'boinet' and 'mfp' are suggested packages, used by the tests, the
  comparison and `est_method = "fp_logistic"`.

## TITE-BOIN-ET

* `tite_boin_et_simulate()` simulates TITE-BOIN-ET trials (Takeda, Morita and
  Taguri, 2020) in C++ (`src/tite_boinet_core.h`), following
  `boinet::tite.boinet()`: patients arrive over time, a pending patient counts
  by the fraction of the DLT or response window followed, and the decision
  waits while fewer than half of the patients at the current dose have
  finished both assessments. `min_completed_share` changes this share; the
  simulations of the article did not wait (`min_completed_share = 0`). The
  times to DLT and to response come from Weibull or uniform distributions
  linked by a Gaussian copula. The dose rules are those of
  `boin_et_simulate()`.

* `sim_tite_boin_et()` and `sim_tite_boin_et_multi()` select the OBD with
  `boin_et_select_obd()` and report the operating characteristics and the
  average trial duration, with the `print()` methods of BOIN-ET.

## BF-BOIN-ET

* `sim_bf_boin_et()` and `sim_bf_boin_et_multi()` simulate the backfill
  BOIN-ET design (Takeda, Zhu and Hirakawa, 2025) in C++
  (`src/bf_boinet_core.h`), following `bfboinet::get.oc.backboinet()`: the
  patients who arrive while a cohort is assessed are backfilled one by one to
  the dose with the highest utility among the doses with a promising response
  rate below the current dose, the dose escalation uses the DLT rate with the
  boundaries of BOIN, a dose with enough patients is excluded for good, and
  the OBD is the admissible dose up to the estimated MTD with the highest
  utility. They report the OBD selection, the patients treated and
  backfilled, the DLTs and responses, and the trial duration, with `print()`
  methods.

* `suspend` sets the share of the patients at the current dose whose DLT
  assessment may still be running when the next cohort starts. The default,
  0, waits until the DLT outcomes of the cohort are known, the BF-BOIN-ET-C
  design of Takeda et al. (2025), and 0.33 gives their BF-BOIN-ET-A design.
  `backfill_extension` sets how long backfilling continues after the last
  cohort has been enrolled, in the unit of `window`; the default, 100, is the
  value fixed in 'bfboinet'.

* The utility is set by the scores `psi00` (no DLT and no response), `psi11`
  (both), `psi_resp_no_dlt` and `psi_dlt_no_resp`, which 'bfboinet' calls
  `ppsi10` and `ppsi01`. The last two names spell out the outcome because the
  articles index the scores in opposite orders.

## BOIN12

* `boin12_decision_table()` tabulates the desirability of a dose, the
  posterior probability that its utility exceeds the benchmark under the
  quasi-binomial model, and the rank-based score of Lin et al. (2020), and
  reproduces all 166 cells of their Table 3. It has a `print()` method.

* `boin12_simulate()` simulates BOIN12 trials in C++ (`src/boin12_core.h`):
  the BOIN boundaries rule out an escalation from a toxic dose, the next
  cohort goes to the most desirable admissible dose among the current dose and
  its neighbors, ties to the higher dose, and the dose exploration rule
  escalates after 9 patients when the dose above is untried. A dose made
  inadmissible by its DLT data also makes every higher dose inadmissible, and
  the trial stops when no dose is admissible. `n_earlystop` (off by default)
  stops a trial when the current dose has that many patients and the next
  cohort would stay. The engine agrees trial by trial with a separate Python
  implementation fed with the same random numbers.

* Three arguments of `boin12_simulate()`, `sim_boin12()` and
  `sim_boin12_multi()` select the rules of the simulations of Lin et al.
  (2020) where they differ from the BOIN12 software. `tie_break = "random"`
  breaks ties of desirability at random instead of choosing the higher dose
  (`"higher"`, the default). `no_candidate = "stop"` ends the trial without an
  OBD, with the stopping reason `"no_admissible_candidate"`, when no dose
  allowed for the next cohort is admissible, instead of going to the nearest
  admissible dose (`"nearest"`, the default). `explore_from = "admissible"`
  applies the dose exploration rule only at an admissible dose instead of at
  any dose (`"any"`, the default). The last two rules are inferred from the
  published results, which `inst/reproduce/reproduce_Lin2020.R` compares with
  the three arguments set.

* `boin12_select_obd()` selects the OBD as the admissible dose up to the MTD,
  estimated by isotonic regression, with the highest estimated utility:
  the posterior mean (`est_method = "posterior_mean"`, the default), the
  observed mean or the desirability.

* `sim_boin12()` and `sim_boin12_multi()` report the OBD selection and the
  patients, DLTs and responses at each dose, with `print()` methods.

* `inst/reproduce/reproduce_Lin2020.R` compares the table, the boundaries,
  the example of the Supplement of Lin et al. (2020) and the BOIN12
  desirabilities of Table S.3 of Zhou et al. (2022), sets 10,000 simulated
  trials of each of the 16 scenarios beside the BOIN12 results of Figure 3
  of Lin et al. and Figures S3 and S4 of their Supplement, and
  `inst/validation/compare-with-escalation.R` compares the operating
  characteristics with `escalation::get_boin12()`.

## TITE-BOIN12

* `tite_boin12_simulate()` simulates TITE-BOIN12 trials (Zhou et al., 2022)
  in C++ (`src/tite_boin12_core.h`) with the approximated likelihood: patients
  arrive over time, the DLT and the response are assessed over their windows,
  and the decision waits while fewer than half of the patients at the current
  dose have finished both assessments (`min_completed_share`). A pending
  outcome is replaced by its conditional expectation given the time followed,
  as in Equation (11) of the article (`pending_estimate = "conditional"`, the
  default), or by the estimated probability itself (`"marginal"`), which
  reproduces the desirabilities of Table S.3 of the Supplement. The
  admissibility of a dose uses the same approximated likelihood, so that a
  response that is merely late does not exclude a dose. After the last cohort,
  or when `n_earlystop` is reached, the trial waits for every outcome and the
  OBD is selected from the complete data. With `min_completed_share = 1` and
  fixed arrivals the trials are exactly those of `boin12_simulate()`, also
  with the arguments `tie_break`, `no_candidate` and `explore_from` of BOIN12,
  which `tite_boin12_simulate()`, `sim_tite_boin12()` and
  `sim_tite_boin12_multi()` share.

* `sim_tite_boin12()` and `sim_tite_boin12_multi()` report the OBD selection,
  the patients, DLTs and responses at each dose and the average trial
  duration, with the `print()` methods of BOIN12.

* The engine agrees trial by trial with a separate Python implementation fed
  with the same random numbers. `inst/reproduce/reproduce_Zhou2022.R` compares
  the desirabilities of Table S.3 of the Supplement and the operating
  characteristics of BOIN12 and TITE-BOIN12 in Figure 3 of the article.

## Comb-BOIN12

* `comb_boin12_simulate()` simulates trials of the Comb-BOIN12 design of Lu,
  Zhang, Yuan and Lin (2025) for two drugs given together in C++
  (`src/comb_boin12_core.h`), with the utility, the desirability and the
  admissibility of BOIN12 at every combination. The utility `pE - wT * pT` of
  the article is `psi00 = 100 * wT / (1 + wT)` and `psi11 = 100 / (1 + wT)`,
  and its benchmark weight `w_u` is `benchmark_weight`. A combination made
  inadmissible by its DLT data also makes every combination with at least
  the same levels of both drugs inadmissible; `extrasafe` and `offset` lower
  the cutoff at the lowest combination, as the article does with 0.85. The
  arguments `tie_break` and `no_candidate` of BOIN12 apply along both drugs:
  `"higher"` keeps the tied combinations with the highest sum of the levels
  and draws one of them when several remain. There is no exploration rule;
  with one row or one column the trials are those of `boin12_simulate()`
  without it.

* `comb_boin12_select_obd()` selects the OBD combination: the DLT rates are
  made monotone by bivariate isotonic regression, computed in C++ to a
  tolerance of `1e-13`, the MTD combination bounds the candidates, and the
  admissible candidate with the highest estimated utility is selected.
  Distances to the target within `1e-8` and utilities within `1e-6` count as
  ties, so that combinations of one block of the isotonic regression, which
  tie exactly, are ordered by the tie rule and not by the rounding of an
  iterative algorithm. The response probabilities are the
  posterior means of the quadratic logistic model of the article with Cauchy
  priors (`est_method = "logistic"`, the default), computed by adaptive
  importance sampling with Halton points so that no random number is used, or
  the observed rates. `obd_bound = "rectangle"` (the default) takes the
  candidates up to the levels of the MTD combination, as the article states;
  `obd_bound = "contour"` bounds each level of drug B by its own MTD along
  drug A, the rule that reproduces the published selection percentages.

* `tite_comb_boin12_simulate()`, `sim_tite_comb_boin12()` and
  `sim_tite_comb_boin12_multi()` simulate the TITE-Comb-BOIN12 design of the
  same article, in which the DLT and the response may be observed late: the
  rules of Comb-BOIN12 at the decision times and with the interim estimates,
  admissibility and quasi-events of TITE-BOIN12
  (`src/tite_comb_boin12_core.h`), and the trial durations.
  `pending_use = "utility"` counts the pending outcomes only in the utility,
  as Sections 2.2 and 3.1 of the article state the rules, the decision and
  the admissibility using the observed outcomes. With `accrual = "fixed"` and
  `min_completed_share = 1` the trials are those of `comb_boin12_simulate()`
  with the same seed, and with one row of combinations those of
  `tite_boin12_simulate()` without exploration.

* `sim_comb_boin12()` and `sim_comb_boin12_multi()` report the selection,
  the patients, DLTs and responses at every combination and the selection of
  and the patients at the true OBD combinations and at the combinations above
  the DLT limit, with `print()` methods. `oc_results()` and `oc_compare()`
  accept their results.

* The engines and the selection agree trial by trial with a separate Python
  implementation fed with the same random numbers, and the logistic estimates
  agree with long Markov chain Monte Carlo runs within 0.01. When
  `psi00 + psi11` is 100, combinations with the same numbers of patients,
  DLTs and responses tie exactly, whatever the rounding of the two scores;
  BOIN12 and TITE-BOIN12 now compute their quasi-events in the same way (no
  change with whole numbers `psi00` and `psi11`, such as the defaults).

* The documentation of `comb_boin_select_mtd()` now states the effect of the
  small amount added to the isotonic estimates correctly: among equal
  estimates it favors the higher combination below the target and the lower
  one above it, as `BOIN::select.mtd.comb()` does. The selection itself is
  unchanged.

## BOIN-ETC

* `boin_etc_simulate()` simulates trials of the BOIN-ETC design of Kakizume,
  Takeda, Taguri and Morita (2024) for two drugs given together in C++
  (`src/boin_etc_core.h`): the dose matrix is divided into subtrials along
  which the DLT probability increases, and BOIN-ET (with the boundaries of
  `boin_et_boundary()`) is run in each of them until the next combination
  has `n_earlystop` patients. `subtrials = "etc1"`, `"etc2"` or `"etc3"`
  are BOIN-ETC1, 2 and 3 of the figures, tables and code of the article (its
  Section 2.2 describes BOIN-ETC2 and BOIN-ETC3 the other way round). A
  combination excluded by its DLT data also excludes the higher combinations
  of its subtrial, the subtrial ends when its lowest combination is
  excluded, and the trial ends when that combination is A1B1.
  `etc1_start` chooses the starting combination of a row of BOIN-ETC1, one
  level of drug B below the previous candidate (`"diagonal"`, the text and
  figures of the article, the default) or at the same level (`"below"`, the
  code of the authors), and `rule4 = "authors_code"` applies the choice among
  the neighboring combinations of the code of the authors, which moves one
  combination down when the combination below the current one is untreated.
  Only `etc1_start = "below"` with `rule4 = "authors_code"` reproduces the
  published allocation of BOIN-ETC1.

* `boin_etc_select_obd()` selects the optimal dose combination: the DLT
  rates are made monotone by unweighted bivariate isotonic regression,
  computed in C++ to a tolerance of `1e-13`, the MTD combination of every
  row (the highest of the tied ones, within `1e-8`) bounds the candidates of
  that row, and the treated, not excluded candidate with the highest
  posterior mean of the utility is selected. `odc_tie = "all"` counts every
  tied candidate, the definition of the percentages of the article.

* `sim_boin_etc()` and `sim_boin_etc_multi()` report the selection, the
  patients, DLTs and responses at every combination and the selection of and
  the patients at the true optimal combinations and at the combinations above
  the target, with `print()` methods. `oc_results()` and `oc_compare()`
  accept their results.

* The engine and the selection agree trial by trial with a separate Python
  implementation fed with the same random numbers.
  `inst/reproduce/reproduce_Kakizume2024.R` compares the allocation and the
  selection with Table 3 and Supplementary Tables S1 to S9 of the article,
  and `inst/validation/compare-with-BOINETC.R` runs the code of the authors
  in R. With the random numbers of R after `set.seed(1234)`, the code gives
  the published allocation and selection of 48 of the 60 settings to the
  printed digits, but for values that appear to be misprints. On the same
  trials, `boin_etc_select_obd()` selects other combinations only when
  combinations tie for the MTD combination of a row in the exact isotonic
  regression, of which the code keeps those that the convergence errors of
  `Iso::biviso()` leave closest to the target, or when the code selects an
  untreated or excluded combination.

## Accelerated BOIN

* `boin_simulate()`, `sim_boin()` and `sim_boin_multi()` gain the argument
  `titration_doses`, which ends the titration at a given dose: when no DLT is
  seen up to it, full cohorts start at the next dose. With
  `n_earlystop = 6` this is the accelerated BOIN design of Kojima et al.
  (2024), and the warning on a small `n_earlystop` is then not given. The
  default, `NULL`, keeps the titration of the BOIN package, and the results
  of earlier versions do not change. `inst/reproduce/reproduce_Kojima2024.R`
  compares Tables 1 and 2 of the article.

## Tables and figures

* `oc_scenarios()` tabulates scenarios of DLT and response probabilities,
  for a single agent or a drug combination, with the true utility of every
  dose and the true MTDs and OBDs. The utility takes the correlation of the
  DLT and the response as the simulation functions do; the true OBDs can
  also be given directly when an article defines them otherwise.

* `oc_results()` collects the results of any of the simulation functions,
  of their `_multi` versions and of `oc_3p3()` into one long data frame, and
  `oc_compare()` judges them against the true doses of `oc_scenarios()`:
  selection of the true MTD or OBD, patients treated there, selection of and
  patients at overdoses, trials without a selection, sample size and
  duration. All three have `print()` and `plot()` methods.

## Checks of the time-to-event and backfill designs

* `inst/reproduce/reproduce_Takeda2020.R` sets the simulations beside Table 2
  of Takeda et al. (2020) (TITE-BOIN-ET and BOIN-ET), and
  `inst/reproduce/reproduce_Takeda2025.R` beside Tables 1 and 2 of Takeda et
  al. (2025) (BF-BOIN-ET with and without the suspension rule).

* The C++ engines of TITE-BOIN-ET and BF-BOIN-ET agree trial by trial with
  separate implementations written in Python from the R code of 'boinet' and
  'bfboinet', fed with the same random numbers.
  `inst/validation/compare-with-boinet.R` now also compares
  `sim_tite_boin_et()` with `boinet::tite.boinet()`, and
  `inst/validation/compare-with-bfboinet.R` compares `sim_bf_boin_et()` with
  `bfboinet::get.oc.backboinet()`.

## Documentation

* The new vignette `vignette("simFastBOIN-boin-et")` describes BOIN-ET,
  TITE-BOIN-ET, BF-BOIN-ET and BOIN-ETC, the selection of the OBD and the dose
  moves that differ from 'boinet', and summarizes the reproductions and
  comparisons from the tables of judgments that the scripts write.

* The new vignette `vignette("simFastBOIN-boin12")` describes BOIN12,
  TITE-BOIN12, Comb-BOIN12 and TITE-Comb-BOIN12: the desirability score
  table, the rules of the simulations of the articles that the defaults do
  not follow, a comparison of BOIN, BOIN-ET and BOIN12 on the same scenarios
  with `oc_compare()`, and the reproductions and comparisons.

* `vignette("simFastBOIN-introduction")` describes the accelerated BOIN design
  and the functions `oc_scenarios()`, `oc_results()` and `oc_compare()`, and
  `vignette("simFastBOIN-combination")` points to the designs that select an
  optimal dose combination.

# simFastBOIN 2.1.0

This version adds the time-to-event BOIN (TITE-BOIN) design, the BF-BOIN and
BE-BOIN designs, which backfill patients to lower doses while the dose
escalation waits, the BOIN design for trials of two drugs given together, and a
decision table for the 3+3 design. No existing function or result changes.

## TITE-BOIN decision table

* `tite_boin_decision_table()` tabulates the decision rule of the TITE-BOIN
  design, which lets new patients be treated while the DLT assessment of earlier
  patients is still pending. Two methods are offered: the single mean
  imputation of Yuan et al. (2018), `method = "imputation"`, with boundaries on
  the standardized total follow-up time (STFT), and the effective sample size
  of Lin and Yuan (2020), `method = "ess"`, with boundaries on the effective
  sample size (ESS). The table has `print()` and `plot()` methods, and the
  printed layout follows the published tables.

* The suspension rules of the two articles are available for both methods:
  `max_pending_ratio` (accrual is suspended when more than half of the patients
  are pending, the default of the imputation method) and `min_completed`
  (escalation requires two patients with a completed assessment, the default of
  the effective sample size method).

* `min_follow_up` adds rule 2 of Chen et al. (2025): while patients are
  pending, escalation requires every pending patient at the current dose to
  have been followed for at least this fraction of the assessment window, and
  accrual is suspended until then. It is available in
  `tite_boin_decision_table()`, `tite_boin_simulate()`, `sim_tite_boin()` and
  `sim_tite_boin_multi()`, and defaults to 0, which changes nothing. With
  `max_pending_ratio = 0.49` (their rule 1 with A = 51) and
  `min_follow_up = 0.25`, the imputation method reproduces every entry of
  Table 1 (target 0.3) of Chen et al. (2025) and of Table A1 (target 0.25) of
  the supplementary materials of Chen et al. (2026), 131 states in all. The
  printed table shows the condition, for example `>= 0.22 & MF >= 0.25`.

* With the defaults, the imputation method reproduces every entry of Table 1
  (target 0.2) of Yuan et al. (2018) and Table S1 (target 0.3) of its
  supplementary appendix, 640 states in all, and the printed table has the same
  rows. The published tables allow de-escalation when the observed DLT rate
  equals the target, which equation (5) of the appendix, read literally, does
  not; the implementation follows the tables. The dose decisions of the trial
  example in Figure 1 of that article are reproduced as well.

* Elimination counts the pending patients as treated without DLT, as both
  articles define it. States without pending patients agree with
  `boin_decision_table()`.

## TITE-BOIN simulation engine

* `tite_boin_simulate()` simulates TITE-BOIN trials in which patients arrive
  over time and DLTs can occur late in the assessment window. It returns the
  patient and DLT counts of `boin_simulate()` together with the duration of each
  trial, the number of suspensions of accrual and the time spent suspended. The
  engine is written in C++ (`src/tite_core.h`), like the BOIN engine.

* The dose for a cohort is decided when its first patient arrives, with the
  rules of `tite_boin_decision_table()`; the engine and the table are checked
  against each other at every state. A suspension lasts until the next pending
  patient at the current dose completes the assessment or, with a positive
  `min_follow_up`, until the minimum follow-up is reached, whichever comes
  first. Because DLTs can be
  observed after a dose has been left, the elimination rule is checked at every
  dose, not only the current one.

* Times between arrivals can be exponential, uniform or fixed, and the time to
  DLT Weibull, with the share of DLTs in the second half of the window set by
  `late_fraction`, or uniform.

* `prior_weights` weights the follow-up of the pending patients by a piecewise
  uniform prior for the time to DLT over the three thirds of the window, as in
  Appendix D of Yuan et al. (2018) and Section S1 of Lin and Yuan (2020). The
  decision table is unchanged; only the follow-up compared with it is. Equal
  weights, the default, use the unweighted computation, so the default results
  do not change. The weighted trials agree exactly with the separate Python
  implementation, which uses the formula of Lin and Yuan's reference code.

* One uniform variate per patient, drawn from R's stream in the order used by
  `boin_simulate()`, decides both whether a DLT occurs and when; arrivals come
  from a separate generator. Whenever no patient is pending at a decision the
  trials are therefore identical to those of `boin_simulate()` with the same
  seed, which the tests check for both methods and several designs.

* The trials agree exactly with a separate implementation written in Python,
  run on a replica of R's random number stream.

## TITE-BOIN operating characteristics

* `sim_tite_boin()` and `sim_tite_boin_multi()` summarize simulated TITE-BOIN
  trials in the same way as `sim_boin()` and `sim_boin_multi()`, with the same
  options for selecting the MTD (`bound_mtd`, `mtd_max_estimate`,
  `min_mtd_sample`) and the same `overdose` component. They add the average
  and standard deviation of the trial duration, the percentage of trials with
  accrual suspended, and the average number and total time of suspensions.

* The results inherit from the BOIN classes (`c("tite_boin_oc", "boin_oc")` and
  `c("tite_boin_oc_multi", "boin_oc_multi")`), so code that reads the result of
  `sim_boin()` also reads them. They have their own `print()` methods, which
  also show the design and the timing; the table of `sim_tite_boin_multi()`
  has two rows more per scenario, for the duration and the suspensions.

* Whenever no patient is pending at a decision, every summary agrees exactly
  with that of `sim_boin()` under the same seed, including the options of MTD
  selection and the safety rules; the tests check this for both methods.

## Backfilling: BF-BOIN and BE-BOIN

* `bf_boin_simulate()` simulates the BOIN design with backfilling (BF-BOIN) of
  Zhao et al. (2024). The dose escalation proceeds in fully staggered cohorts
  with the BOIN decision table, and a patient who arrives while it waits is
  backfilled to a lower dose that has not been eliminated, has fewer than
  `n_cap` patients and has a response observed at or below it. Doses are
  closed for backfilling, and conflicts between the backfilled doses and the
  current dose are resolved with pooled data, as in Table 2 of the article.

* `be_boin_simulate()` simulates BE-BOIN of Chen et al. (2026), which combines
  backfilling with TITE-BOIN. The escalation decides while patients are
  pending, with rules 1 and 2 of the article (`max_pending_ratio = 0.49` and
  `min_follow_up = 0.25`), and every estimate imputes the pending patients,
  over several doses with equation (2) of the article.

* `sim_bf_boin()`, `sim_be_boin()`, `sim_bf_boin_multi()` and
  `sim_be_boin_multi()` summarize the simulated trials. The results inherit
  from the BOIN classes, count backfilled patients among the patients treated,
  and add the backfilled patients and the responses at each dose, the trial
  duration, the waits of the escalation and the patients turned away. The MTD
  is selected from the data of all patients, as in `sim_boin()`. They have
  their own `print()` methods, and the table of the `_multi` functions has
  seven rows per scenario.

* Points that the articles leave open are arguments: `no_slot` decides
  whether a patient with nowhere to go waits or is turned away,
  `conflict_dose` whether the highest (Zhao et al., 2024) or the lowest (Chen
  et al., 2026) conflicting dose anchors the pooled estimate, and
  `backfill_dose` which open dose receives a backfill patient. The time to
  response is Weibull (`resp_window`, `resp_late_fraction`), and the DLT and
  the response are linked by a Gaussian copula (`resp_cor`, 0 by default), as
  in Takeda et al. (2025).

* The engine is written in C++ (`src/backfill_core.h`). The DLT variates of
  the escalation are drawn from R's stream as in `boin_simulate()`; arrivals,
  backfill patients and responses have generators of their own. When no dose
  is ever opened for backfilling the trials are identical to those of
  `boin_simulate()` (BF-BOIN) and of `tite_boin_simulate()` (BE-BOIN with
  `no_slot = "wait"`) with the same seed, which the tests check.

* With `stay_on_1_of_3 = TRUE`, `n_earlystop = 9`, `accrual = "uniform"` and
  `no_slot = "leave"`, the operating characteristics of Table 4 of Zhao et al.
  (2024), twelve scenarios of BF-BOIN and five of BOIN, are reproduced within
  simulation error, durations included; over 20 seeds of 10,000 trials the
  largest gaps were 2.1 points for the selection, 0.4 patients and 0.17
  months. The article states a Poisson process of arrivals, but its durations
  agree with uniform times between arrivals and patients turned away while the
  escalation waits. `inst/validation/compare-with-Zhao2024.R` repeats the
  comparison.

* The trial examples of Figure 1 of Zhao et al. (2024) and of Section S2 of
  the supplementary materials of Chen et al. (2026) are replayed patient by
  patient, and the simulated trials agree exactly with a separate
  implementation written in Python, over 2,700 trials of both designs.

## BOIN design for drug combinations

* `comb_boin_simulate()` simulates the BOIN drug-combination design of Lin and
  Yin (2017), which looks for a single MTD in a matrix of dose combinations of
  two drugs. After each cohort the decision of the single-agent design is
  taken at the current combination; an escalation goes one level up of either
  drug, and a de-escalation one level down, to the neighbor with the larger
  posterior probability that its DLT rate lies between the boundaries. The
  engine is written in C++ (`src/comb_core.h`) and follows
  `BOIN::get.oc.comb()` with `mtd.contour = FALSE`: the same rules, applied in
  the same order, with the random variates drawn in the same order, so that
  with the same seed the two produce the same trials.

* `comb_boin_select_mtd()` selects the MTD combination from completed trials,
  with the bivariate isotonic regression of `Iso::biviso()`, the routine the
  reference implementation calls. Combinations that the regression pools are
  equal in exact arithmetic but can differ in their last digits, and these
  digits decide between combinations on the same anti-diagonal, so another
  algorithm would select differently in a few trials out of a thousand. The
  'Iso' package, which version 1.4.0 moved to Suggests, is therefore imported
  again.

* `sim_comb_boin()` and `sim_comb_boin_multi()` summarize the simulated trials:
  the selection percentage and the average numbers of patients and DLTs at
  every combination, the percentages of trials selecting a true MTD
  combination and of patients treated at one, and the exposure to
  combinations above `overdose_cutoff`. They have their own `print()`
  methods, and the table of `sim_comb_boin_multi()` has one row per scenario.

* `comb_boin_next()` gives the combination for the next cohort of a trial in
  progress, from the numbers of patients and DLTs observed at every combination
  and the combination of the last cohort. It reports the boundaries that apply,
  the two candidates with their posterior probabilities and scores, both
  candidates when they tie, the combinations eliminated, and why the trial
  stops when it does. It has `print()` and `plot()` methods; the plot draws the
  dose matrix with the data, the current and next combinations and the
  eliminated ones. The decision is taken by the same C++ code as in the
  simulations (`comb_decide()` in `src/comb_core.h`), so a trial is run with the
  rules whose operating characteristics `sim_comb_boin()` reports. These differ
  from `BOIN::next.comb()` in three respects: `BOIN::next.comb()` stops as soon
  as `n.earlystop` patients have been treated at the current combination,
  blocks an escalation of one drug when a combination at the next level of that
  drug has an observed DLT rate at or above `lambda_d`, and breaks ties at
  random. The tests check that `comb_boin_next()` takes the decision of
  simulated trials after every cohort, and that it agrees with
  `BOIN::next.comb()` on random data wherever the rules coincide. The C++ code
  also takes the same decision as the Python implementation at each of the
  40,664 states reached by 3,300 simulated trials in 11 settings.

* The tests reproduce the example of Yan et al. (2020): every number printed
  by `summary()` in the article agrees to the printed precision. The one
  average exactly half-way between two printed values, 1.905 patients at A2B4,
  is printed 1.90 in the article and 1.91 by R 4.6.0; how such a value is
  rounded depends on the version of R. The tests also check exact agreement
  with `BOIN::get.oc.comb()` in nine settings whenever the BOIN package is
  installed. The simulated trials also agree exactly with a separate
  implementation written in Python and run on a replica of R's random number
  stream, over 5,500 trials in 11 settings.
  `inst/validation/compare-comb-with-BOIN.R` extends the comparison with BOIN
  to 200 configurations, and `inst/validation/compare-with-LinYin2017.R` sets
  the results beside Tables 2 and 3 of Lin and Yin (2017), whose rules differ
  in some details.

* The waterfall design of Zhang and Yuan (2016), which looks for the MTD
  contour, and a time-to-event version of the drug-combination design are not
  provided.

## 3+3 decision table

* `decision_table_3p3()` tabulates the rule of the 3+3 design used by
  `oc_3p3()` and `sim_3p3()`: the decision at three and at six patients during
  dose escalation and, for `mtd_rule = "expand"`, during the search for the MTD
  that follows. It has `print()` and `plot()` methods in the style of the BOIN
  decision table. The tests follow the table through every possible course of a
  trial and obtain the operating characteristics of `oc_3p3()`.

## Documentation

* A second vignette, `vignette("simFastBOIN-tite")`, describes the TITE-BOIN
  design. The introductory vignette gains a pointer to it and a paragraph on
  `decision_table_3p3()`.
  The README describes the design, how it is checked, and its speed: 10,000
  trials of five doses and 20 cohorts of three take 0.09 s with the imputation
  method and 0.11 s with the effective sample size, against 0.05 s for
  `sim_boin()` on the same machine.

* The labels of `plot.tite_boin_decision_table()` are shorter, the codes over
  the boundaries (for example `E/S` over `>=2.15`), so that they fit in the
  cells.

* A third vignette, `vignette("simFastBOIN-backfill")`, describes the BF-BOIN
  and BE-BOIN designs, and the README and the introductory vignette point to
  it.

* The DESCRIPTION cites the two TITE-BOIN articles and the three articles on
  backfilling and its suspension rules.

* A fourth vignette, `vignette("simFastBOIN-combination")`, describes the
  drug-combination design, and the README and the introductory vignette point
  to it. The DESCRIPTION cites Lin and Yin (2017). The vignette also shows the
  decision table for a protocol and how to run a trial with
  `comb_boin_next()`.

* A pkgdown configuration (`_pkgdown.yml`) groups the reference pages by design
  and lists the four vignettes as articles.

* The `plot()` methods of `boin_decision_table()`,
  `tite_boin_decision_table()` and `comb_boin_next()` draw the legend key of
  every decision or status, including those that do not occur in the figure,
  which `ggplot2` 3.5.0 and later otherwise leave blank.

# simFastBOIN 2.0.0

The first release since 1.3.2, and a breaking one. Versions 1.4.0 and 1.5.0 were
development versions on GitHub and were never released; their entries are kept
below because the changes they describe are all part of this release. Anyone
upgrading from 1.3.2 should read the 1.4.0 entry, which lists the renamed
functions and the corrected defects.

This section covers three design options requested for a quality control review.
Every new option is off by default, so a simulation run with version 1.4.0
arguments produces exactly the same trials, selections and summary figures. Two
names in the returned objects did change, and are listed under their respective
headings below: `overdose60` and `overdose80` moved into the new `overdose`
component, and one reason code was renamed.

## Staying on one DLT out of three

* `stay_on_1_of_3` makes one DLT out of three patients a stay rather than a
  de-escalation, which is the modification offered by the BOIN web application.
  It is available on `boin_boundary()`, `boin_decision_table()`,
  `boin_simulate()`, `sim_boin()` and `sim_boin_multi()`.

  The option raises the de-escalation boundary at three patients from one DLT to
  two, so exactly one cell of the decision table changes and nothing else moves.
  It is applied only where one DLT out of three currently de-escalates: it never
  overrides an escalation, and it never overrides an elimination, which is a
  safety rule. `boin_boundary()` reports whether it took effect in
  `stay_on_1_of_3_applied`.

  With the default thresholds it takes effect for target rates from about 0.098
  to 0.279. Note that the web application describes the range as 0.25 to 0.279;
  the upper end agrees, but the lower end is 0.098, below which one DLT out of
  three already eliminates the dose.

## Bounding the MTD by a stated DLT rate

* `mtd_max_estimate` caps the isotonic estimate a dose may have and still be
  selected as the MTD, on `boin_select_mtd()`, `sim_boin()` and
  `sim_boin_multi()`. It changes only the selection, never the dose-finding, so
  the decision table is untouched.

  `bound_mtd` caps the estimate at the de-escalation boundary, which always lies
  above the target rate. A cap at or below the target therefore cannot be
  expressed with `bound_mtd` at all, and needs `mtd_max_estimate`.

* The reason code `"no_dose_below_lambda_d"` returned by `boin_select_mtd()` is
  now `"no_dose_below_bound"`, since the bound need not be the de-escalation
  boundary.

## Exposure to overly toxic doses

* `sim_boin()` and `sim_boin_multi()` gain `overdose_cutoff`, and report an
  `overdose` component describing how far patients were exposed to doses whose
  true DLT probability exceeds it:

  - `pct_patients`, the percentage of all simulated patients treated at such a
    dose, that is the probability that a patient is dosed above the cutoff;
  - `pct_patients_by_trial`, the same percentage computed within each trial and
    then averaged;
  - `avg_n_patients`, `pct_trials_any`, `pct_trials_over_60` and
    `pct_trials_over_80`.

  and how often the design ends up recommending such a dose, which is a separate
  question from exposure during the trial:

  - `pct_trials_mtd_above`, the percentage of all trials whose selected MTD is
    above the cutoff;
  - `pct_mtd_above_when_selected`, the same among the trials that selected an
    MTD at all.

  The cutoff defaults to `target`, which reproduces the earlier behavior.

* This replaces the `overdose60` and `overdose80` components of the result, which
  are now `overdose$pct_trials_over_60` and `overdose$pct_trials_over_80`. They
  are reported as zero rather than `NA` when no dose exceeds the cutoff, and the
  `doses` component names the dose levels that do.

## Choosing a de-escalation boundary directly

* `boin_p_tox()` returns the value of `p_tox` for which the de-escalation
  boundary equals a stated value, inverting `boin_lambda()`. A design team that
  wants to tighten the boundary usually states the boundary rather than the
  threshold behind it, and this turns that statement into the argument the other
  functions take.

  A boundary close to the target needs a threshold close to the target, and the
  rest of the package refuses a `p_tox` within ten percent of `target`. The
  inverted value is still returned, but a warning names the smallest boundary
  that leaves a usable threshold: about 0.315 for a target of 0.30.

## Printing

* `print()` on the result of `sim_boin_multi()` gains `percent` and `digits`,
  which the method for a single scenario already had. Passing `percent = TRUE`
  to the cross-scenario method previously fell into `...` and was ignored, so
  the counts were printed as though nothing had been asked for.

## Documentation

* `start_dose`, added in version 1.4.0, is now documented as being ignored when
  `titration = TRUE`, because the titration phase always begins at the lowest
  dose. This matches the reference implementation.

* The spelling throughout the documentation is now the United States spelling
  declared by the `Language` field, and `tests/spelling.R` keeps it that way.
  The `colours` argument of `plot.boin_decision_table()` is renamed `colors`
  accordingly. The argument was introduced in version 1.4.0, which was never
  released, so no released code is affected.

## The 3+3 design as a comparator

* `oc_3p3()` returns the operating characteristics of the traditional 3+3
  design. The probability of every path the design can take is evaluated in
  closed form, so the result carries no Monte Carlo error and takes no seed.
  `sim_3p3()` simulates the same design and is there to confirm the closed form.
  Both return an object of class `oc_3p3` whose components carry the same names
  as those of `sim_boin()`, so the two designs can be tabulated side by side, and
  `print()` lays them out the same way.

* `mtd_rule` selects between the two definitions of the MTD in common use.
  `"previous"` names the dose below the one declared too toxic. `"expand"` takes
  the MTD to be the highest dose at which at most one of six patients had a DLT,
  expanding a dose that has only three patients before assessing it. Neither is
  universally the standard, so the comparator being matched should be checked
  before a rule is chosen. The two agree on which doses are visited during
  escalation, and therefore on the probability that a patient is exposed to an
  overly toxic dose; they differ in the sample size and in the selection
  percentages.

* `overdose_cutoff` works as it does for `sim_boin()` and defaults to one third
  here, which is the rate the regulatory question is asked about. The 3+3 design
  has no fixed maximum sample size, so the within-trial exposure percentage and
  the 60 percent and 80 percent thresholds reported for BOIN have no counterpart
  and are not returned.

## Validation

* `inst/validation/compare-with-BOIN.R` compares the package against `BOIN` over
  200 configurations, five dose-toxicity scenarios by five designs by eight
  option combinations, and covers `extrasafe`, `titration` and `boundMTD`. All
  200 agree to within floating point on the selection percentages, the patients
  and DLTs at each dose, the percentage of trials ending without an MTD and the
  totals. This is a defect test rather than a Monte Carlo comparison, because the
  two implementations draw the same variates in the same order for a given seed.
  The options that `BOIN` has no counterpart for are checked separately, to
  confirm that each leaves the simulated trials untouched when it is off.

# simFastBOIN 1.4.0

This is a substantial revision. The simulation engine has been rewritten in C++,
the user-facing functions have been renamed, and several defects that affected
numerical results have been corrected. Results obtained with version 1.3.2 are
not reproduced by this version, even with the same seed.

## Exact agreement with the BOIN package

The engine now draws one uniform variate per patient, in enrollment order, and
applies the decision rules in the same order as
`BOIN::get.oc()`. With the same seed and matching arguments the two
implementations agree trial by trial rather than only on average. The test suite
checks this directly against `BOIN::get.oc()` and `BOIN::get.boundary()` when the
BOIN package is installed.

Version 1.3.2 claimed in its README that random number generation was the same as
that of the BOIN package when the seed was fixed. That claim was incorrect,
because trials were generated across trials rather than one at a time, and it has
been removed.

## Corrected defects

* `max_total_pts` was computed as `n_cohort * cohort_size[1]`, so a trial with a
  vector `cohort_size` used the wrong maximum sample size. It is now the sum of
  the cohort sizes.

* With `titration = TRUE` and a vector `cohort_size`, the titration phase tested
  `if (cohort_size > 1)` on a vector, which stops with an error on R 4.2 and
  later.

* The early stopping rule at `n_earlystop` was evaluated after the dose had
  already been changed, and therefore against the wrong dose. It is now evaluated
  at the dose that was just treated, before the transition, as in the reference
  implementation. This changes operating characteristics, in particular when a
  de-escalation lands on a dose that has already accrued `n_earlystop` patients.

* MTD selection now re-derives dose elimination from the final data of the trial
  instead of carrying the elimination status over from the dose-finding stage.
  The two differ when a trial ends at its maximum sample size, because the last
  cohort was then never checked against the elimination boundary.

* The isotonic fit used for MTD selection is now computed over the admissible
  doses only. Previously every treated dose entered the fit, so an eliminated
  dose could change the estimates at the doses below it and hence the selected
  MTD.

* The tie-breaking perturbation was applied twice, once inside
  `isotonic_regression()` and once inside `select_mtd()`, and was indexed by dose
  level rather than by position in the admissible set. It is now applied once,
  during MTD selection, and `boin_isotonic()` returns the unperturbed estimates.

* `get_pts_and_tox()` called `set.seed()` without restoring the state of the
  random number generator. The new functions restore it on exit.

* Random number generation was not consistent within a trial: the titration
  phase, the ordinary cohorts and the final partial cohort used three different
  schemes, the last of which drew a binomial variate. All patients now use a
  single scheme.

## Renamed functions

The old names remain available and issue a deprecation warning, but they return
the value of the replacement function, which is not the value returned by version
1.3.2. See `?"simFastBOIN-deprecated"`.

| Version 1.3.2 | Version 1.4.0 |
|---|---|
| `get_boin_boundary()` | `boin_lambda()`, and `boin_boundary()` for the integer boundaries |
| `get_boin_decision()` | `boin_decision_table()` |
| `get_boin_stopping_boundaries()` | `boin_boundary(extrasafe = TRUE)` |
| `get_pts_and_tox()` | `boin_simulate()` |
| `isotonic_regression()` | `boin_isotonic()` |
| `select_mtd()` | `boin_select_mtd()` |

`sim_boin()` and `sim_boin_multi()` keep their names.

## Changed interfaces

* The design arguments come first: `sim_boin(target, p_true, n_cohort,
  cohort_size, ...)`, with `n_trials` now an optional argument. Calls that relied
  on positional matching must be updated.

* `boundMTD` is now `bound_mtd`, and `return_details` is now `keep_trials`.

* `sim_boin()` returns an object of class `boin_oc` directly, rather than a list
  with a `summary` element. Selection percentages, average patient counts and
  average DLT counts are separate components, and the percentage of trials
  selecting no MTD is no longer appended to the selection vector.

* `sim_boin()` and `sim_boin_multi()` now return their result visibly, so calling
  them at the prompt prints the summary and assigning the result is silent.
  Previously the result was invisible unless `verbose = TRUE`.

* `sim_boin_multi()` uses the same seed for every scenario. Version 1.3.2 used
  `seed + i` for the i-th scenario, so a scenario simulated in a multi-scenario
  run did not match the same scenario simulated on its own.

* `sim_boin_multi()` accepts a plain named list of `p_true` vectors in addition
  to the list of `list(name, p_true)` used previously.

* `boin_select_mtd()` takes the patient and DLT counts, not precomputed isotonic
  estimates, since the fit depends on which doses are admissible.

* Progress messages are emitted with `message()` rather than `cat()`, so they can
  be suppressed with `suppressMessages()`.

## New

* `start_dose`, matching the reference implementation, which did not exist in
  version 1.3.2.
* `boin_boundary()` returns the escalation, de-escalation, elimination and safety
  stopping boundaries as integer DLT counts indexed by sample size, with a print
  method that can restrict the table to the end of each cohort.
* `boin_decision_table()` returns a classed object with `print` and `plot`
  methods. The print method blanks the impossible combinations and can restrict
  the table to the sample sizes reached at the end of a cohort; the plot method
  draws the table as a grid of colored cells, each carrying its decision letter
  so that the figure does not rely on color alone. It needs `ggplot2`, which is
  suggested rather than required.
* `boin_stopping_table()` lays the safety stopping boundary out as a two-row
  table, one column per sample size.
* `boin_simulate()` returns the raw trial data with a `stop_reason` for every
  trial, and has a print method.
* Operating characteristics now include the risk of overdosing 60 and 80 percent
  of patients, and the distribution of the reason for stopping.
* Arguments are validated. Out-of-range probabilities, thresholds too close to
  the target, non-monotone `p_true` and inconsistent counts are reported instead
  of propagating silently.

## Dependencies

* Added `Rcpp`, and `utils` for one call to `globalVariables()`.
* Added `ggplot2` to `Suggests`, used only by `plot.boin_decision_table()`.
* Removed `Iso`. The pool adjacent violators algorithm is implemented in C++;
  `Iso` is now only suggested, and used in a test that checks the two agree.
* `knitr` and `kableExtra` moved out of `Imports`. `knitr` is suggested and used
  for the vignette and for the optional table output; `kableExtra` is no longer
  used.
* `utils` was imported in `NAMESPACE` but absent from `DESCRIPTION`. Neither
  imports it now.

## Performance

The C++ engine replaces the vectorised R implementation, and MTD selection no
longer calls out to R once per trial. For five doses, 20 cohorts of three,
`n_earlystop = 18` and 10000 simulated trials, one measurement gave 9.98 seconds
elapsed for `BOIN::get.oc()` against 0.07 seconds for `sim_boin()`. The ratio
depends on the design and on the machine, and the README carries the code that
reproduces the measurement.

# simFastBOIN 1.3.2

## New Features

### Customizable Safety and Toxicity Thresholds

* **Added `p_saf` and `p_tox` parameters to core functions**
  - `sim_boin()`, `sim_boin_multi()`, and `get_pts_and_tox()` now accept custom safety and toxicity thresholds
  - `p_saf`: Highest toxicity probability deemed acceptable for safety (default: 0.6 * target)
  - `p_tox`: Lowest toxicity probability deemed unacceptable for toxicity (default: 1.4 * target)
  - These parameters are passed to `get_boin_boundary()` for boundary calculation
  - Stored in summary objects for reference and reproducibility

* **Usage examples**
  ```r
  # Using default thresholds (0.6 * target and 1.4 * target)
  result_default <- sim_boin(
    n_trials = 10000,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )
  
  # Using custom thresholds
  result_custom <- sim_boin(
    n_trials = 10000,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    p_saf = 0.15,  # Custom safety threshold
    p_tox = 0.45,  # Custom toxicity threshold
    seed = 123
  )
  ```

## Documentation Improvements

* **Enhanced parameter documentation**
  - Added comprehensive roxygen2 documentation for `p_saf` and `p_tox` parameters
  - All @param entries now use consistent multi-line format for better readability
  - Updated examples to demonstrate custom threshold usage

* **Expanded test coverage**
  - Added tests for `p_saf` and `p_tox` parameter handling in `sim_boin()`
  - Added tests for `p_saf` and `p_tox` parameter handling in `sim_boin_multi()`
  - Added tests for `p_saf` and `p_tox` parameter handling in `get_pts_and_tox()`
  - Verified that custom values are correctly stored in summary objects
  - Verified that default values are correctly calculated when not specified

## Breaking Changes

None. All existing code continues to work as before. The new `p_saf` and `p_tox` parameters have sensible defaults (0.6 * target and 1.4 * target respectively) that match the standard BOIN methodology, ensuring full backward compatibility.

## Internal Changes

* Updated function signatures to include `p_saf` and `p_tox` parameters with NULL defaults
* Modified internal logic to calculate default values when parameters are not specified
* Updated summary objects to store `p_saf` and `p_tox` values for reproducibility

---

# simFastBOIN 1.3.1

## Bug Fixes

* **Fixed DESCRIPTION Date field**
  - Corrected Date from 2024-12-06 to 2025-12-06
  - This was a typo in the year field that was caught during CRAN submission process

No functional changes. All code and features remain identical to version 1.3.0.

---

# simFastBOIN 1.3.0

## New Features

### Progress Message Control

* **Added `verbose` parameter to `sim_boin()` and `sim_boin_multi()`**
  - Control whether progress messages are printed to console
  - `verbose = FALSE` (default): Run silently without progress messages
  - `verbose = TRUE`: Display progress messages as in previous versions
  - Particularly useful for R Markdown documents and vignettes where clean output is desired
  - Results are identical regardless of verbose setting

* **Usage examples**
  ```r
  # Silent mode (default) - ideal for vignettes and reports
  result <- sim_boin(
    n_trials = 10000,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )
  
  # With progress messages
  result <- sim_boin(
    n_trials = 10000,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    verbose = TRUE,
    seed = 123
  )
  ```

## Documentation Improvements

* **Standardized roxygen2 documentation format**
  - All parameter descriptions now use consistent multi-line format
  - Improved readability and maintenance
  - Enhanced consistency across all functions

* **Updated vignettes**
  - Clean output without progress messages
  - Better integration with R Markdown workflows
  - Improved presentation quality

## Breaking Changes

None. All existing code continues to work as before. The new `verbose` parameter defaults to `FALSE`, which changes the default behavior to silent mode, but all functionality remains identical.

## Migration Guide

For users who prefer the previous behavior with progress messages:

```r
# Add verbose = TRUE to see progress messages
result <- sim_boin(
  n_trials = 10000,
  target = 0.30,
  p_true = p_true,
  n_cohort = 48,
  cohort_size = 3,
  verbose = TRUE,  # Add this line
  seed = 123
)
```

---

# simFastBOIN 1.2.1

## New Features: Multi-Scenario Simulation and HTML Output

### New Multi-Scenario Simulation Function

* **Added `sim_boin_multi()` function**
  - Run BOIN simulations across multiple dose-toxicity scenarios simultaneously
  - Automatically orchestrates simulations for each scenario using `sim_boin()`
  - Aggregates results into a unified comparison table
  - Ideal for protocol development evaluating multiple dose-toxicity relationships
  - Returns results organized by scenario for easy comparison
  - Progress messages track simulation status for each scenario

* **Usage example**
  ```r
  scenarios <- list(
    list(name = "Scenario 1: MTD at DL4",
         p_true = c(0.05, 0.10, 0.20, 0.30, 0.45)),
    list(name = "Scenario 2: MTD at DL3",
         p_true = c(0.10, 0.15, 0.30, 0.45, 0.60))
  )
  
  result <- sim_boin_multi(
    scenarios = scenarios,
    target = 0.30,
    n_trials = 10000,
    n_cohort = 48,
    cohort_size = 3,
    seed = 123
  )
  ```

### Enhanced Print Methods: HTML Table Output

* **Added `html` format option to `print.boin_summary()` and `print.boin_multi_summary()`**
  - `kable_format = "html"`: Generate HTML tables with enhanced styling
  - Includes striped rows, hover effects, and responsive formatting via kableExtra
  - Automatically applies visual formatting including bold headers and borders
  - Useful for web-based reports and interactive documents
  - Full support for embedded HTML display in R Markdown documents

* **Updated kable_format parameter documentation**
  - `"pipe"` (default): Markdown pipe table format
  - `"simple"`: Minimal text table format
  - `"latex"`: LaTeX table format
  - `"html"`: HTML table format with enhanced styling (NEW)

* **Enhanced print output examples**
  ```r
  # HTML table for web display
  print(result$summary, kable = TRUE, kable_format = "html")
  
  # Multi-scenario results as HTML
  print(multi_result, kable = TRUE, kable_format = "html")
  ```

## Documentation and Output Formatting Enhancements

### Comprehensive roxygen2 Documentation

* **Added detailed roxygen2-formatted documentation to all functions**
  - `get_boin_boundary()`: Escalation and de-escalation boundary calculation
  - `get_boin_decision()`: Decision table generation with decision rules
  - `get_boin_stopping_boundaries()`: Safety stopping rule table generation
  - `get_pts_and_tox()`: Patient enrollment and toxicity simulation
  - `isotonic_regression()`: Isotonic regression with PAVA algorithm
  - `select_mtd()`: MTD selection with optional boundMTD constraint
  - `sim_boin()`: Main simulation workflow
  - `sim_boin_multi()`: Multi-scenario simulation orchestration
  - `print.boin_summary()`: Summary table printing and formatting
  - `print.boin_multi_summary()`: Multi-scenario summary printing and formatting

* **Enhanced function comments**
  - Added section headers for major processing blocks
  - Detailed explanations of algorithm steps and decision logic
  - Clear descriptions of parameter usage and constraints
  - Clarified relationships between input parameters and outputs

### Print Method Enhancements: print.boin_summary()

* **Added `percent` parameter**
  - `percent = FALSE` (default): Display Avg Pts and Avg DLTs as absolute numbers
  - `percent = TRUE`: Display as percentages of total

All changes are backward compatible. New parameters have sensible defaults.

---

# simFastBOIN 1.2.0

## Bug Fixes and Compatibility

### BOIN Package Compatibility

* **Restored exact compatibility with BOIN package**
  - Changed DLT generation from `runif()` back to `rbinom()`
  - `rbinom()` ensures identical random number sequence as BOIN with same seed
  - Verified: MTD selection results now match BOIN within <0.5% across all scenarios

### Default Parameter Corrections

* **Updated default values to match BOIN standard**
  - `min_mtd_sample`: Changed from 6 to 1
    - Doses with ≥1 patient can now be considered for MTD selection
    - Matches BOIN package default behavior
  - `n_earlystop_rule`: Changed default from "simple" to "with_stay"
    - Trial now stops when n ≥ n_earlystop AND next decision = "Stay"
    - Ensures algorithm convergence before stopping
    - Follows BOIN standard implementation

## Performance Optimizations

### Vectorized Implementation

* **Batch processing for improved performance**
  - Uses vectorized operations for fast simulation
  - Pre-allocated vectors and vectorized computations
  - Early exit for trials with no valid doses

### DLT Generation Optimization

* **Optimized random number generation in sim_boin()**
  - Uses `rbinom()` for accurate DLT generation matching BOIN package
  - Vectorized operations for efficient computation
  - Performance optimization particularly notable for large-scale simulations

### MTD Selection Enhancement

* **Optimized select_mtd() function**
  - Returns NA immediately when no valid MTD candidates exist
  - Avoids unnecessary computations for trials without viable doses
  - Improves overall simulation efficiency

## Internal Improvements

* Enhanced code documentation with detailed optimization rationale
* Improved memory efficiency through pre-allocation
* Maintained backward compatibility with existing APIs

## Documentation Updates

* Updated roxygen2 documentation in sim_boin.R
* Updated parameter descriptions
* Updated README.md with corrected default values

## Breaking Changes

None. The changes restore compatibility and fix defaults to match standard BOIN behavior.
Users running with explicit parameters should see no change.

## Migration Guide

For users upgrading from simFastBOIN 1.1.0:

If you were using defaults:
```r
# Old code (simFastBOIN 1.0.0)
result <- sim_boin(n_trials = 10000, target = 0.30, p_true = p_true, ...)

# New code (simFastBOIN 1.2.0) - No change needed!
# Default behavior now matches BOIN package
result <- sim_boin(n_trials = 10000, target = 0.30, p_true = p_true, ...)
```

If you were using custom parameters:
```r
# Explicitly setting these will ensure consistent behavior across versions
result <- sim_boin(
  n_trials = 10000,
  target = 0.30,
  p_true = p_true,
  min_mtd_sample = 1,        # Now the default
  n_earlystop_rule = "with_stay",  # Now the default
  ...
)
```

---

# simFastBOIN 1.0.0

## Initial Release

### Core Features

* **High-Performance BOIN Simulation**
  - Vectorized implementation for 2-5x speedup over traditional approaches
  - Batch processing of all trials simultaneously at each cohort
  - Efficient matrix operations for state management

* **Automatic Decision Table Generation**
  - `sim_boin()` automatically generates BOIN boundaries and decision tables
  - No need for manual pre-computation
  - User-friendly API with minimal required parameters

* **Safety Features**
  - Optional `extrasafe` parameter for safety stopping at lowest dose
  - Dose elimination rules based on posterior probability
  - Configurable safety thresholds

* **Professional Output**
  - Publication-ready summary tables
  - Optional knitr::kable format for RMarkdown
  - Detailed trial-level results available

### Main Functions

* `sim_boin()`: Run BOIN trial simulations
* `get_boin_boundary()`: Calculate BOIN interval boundaries
* `get_boin_decision()`: Generate decision table
* `get_boin_stopping_boundaries()`: Generate safety stopping boundaries
* `isotonic_regression()`: Apply isotonic regression for dose-toxicity estimation
* `summarize_simulation_boin()`: Aggregate simulation results

### Performance

* 10,000 trials with 9 doses and 48 cohorts: ~1-2 seconds
* Memory efficient: <100MB for typical simulations
* Scales with number of cohorts, not number of trials

### Documentation

* Comprehensive function documentation with examples
* Detailed README with quick start guide
* Performance benchmarks and comparisons
