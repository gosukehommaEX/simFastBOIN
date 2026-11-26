# How the package is checked

## Four kinds of checks

Every design of the package is checked in four ways.

- **Tests.** The test suite in `tests/testthat/` checks every function
  on cases whose results are known: decision tables printed in the
  articles, limiting cases in which a design reduces to a simpler one
  with the same random numbers, and values computed beforehand by a
  separate implementation.
- **Separate implementations.** The trials of every simulation engine
  were compared, trial by trial, with an implementation of the same
  rules written in Python and fed with the random numbers of R, whose
  default generator was written again in Python for that purpose. These
  implementations are development files and are not part of the package.
- **Reference implementations.** The scripts in `inst/validation/`
  compare the package with other R packages and with the code of the
  authors of a design, exactly where the rules and the random numbers
  are the same and within simulation error otherwise.
- **Published results.** The scripts in `inst/reproduce/` compare the
  package with the tables and figures of the articles that introduced
  the designs.

This article summarizes the last two. The scripts write their judgments
to the `results/` folders of `inst/reproduce/` and `inst/validation/`,
and every number below is read from those files.

## Judgments

Every compared value receives one of four judgments.

- **PASS.** The value agrees: exactly, or to the printed digits, for a
  value that involves no random numbers, and within the tolerance stated
  in the script for a simulated value. Where the package and the
  reference draw the same trials, the values must be equal; otherwise
  the tolerance is a number of standard errors of the difference between
  two independent simulations, three or four as the script states, plus
  the rounding of the printed value.
- **EXPLAINED.** The value lies outside the tolerance for a reason that
  was found, such as a misprint in the article, or a rule of the program
  of the authors that the article does not describe. The reason is
  registered in `inst/validation/explained.csv`, which has 399 rows.
- **FAIL.** The value lies outside the tolerance and no reason was
  found. These values are kept and listed below.
- **INFO.** The value is listed for information and not judged: for
  example the results of the default rules of the package where an
  article used other rules, or a quantity that an article does not
  describe well enough to be compared.

A script is run from the root of the source package after
`devtools::load_all()`. It defines its functions, and the function named
in the usage at the top of the script runs the comparison and writes the
results file, for example

``` r

devtools::load_all()
source("inst/reproduce/reproduce_Zhou2021.R")
res <- run_reproduce_zhou2021()
table(res$design, res$judgement)
```

Some scripts take a long time, up to a few hours, and those that compare
with other packages need these packages.

## Summary

| Design | Script | Compared with | PASS | EXPLAINED | FAIL | INFO |
|:---|:---|:---|---:|---:|---:|---:|
| iBOIN | reproduce_Zhou2021.R | Zhou et al. (2021) | 244 | 6 | 0 | 100 |
| gBOIN | reproduce_Mu2019.R | Mu et al. (2019) | 115 | 2 | 17 | 120 |
| gBOIN | compare-with-UnifiedDoseFinding.R | ‘UnifiedDoseFinding’ | 130 | 0 | 0 | 130 |
| Accelerated BOIN | reproduce_Kojima2024.R | Kojima et al. (2024) | 165 | 0 | 3 | 168 |
| Waterfall | reproduce_Zhang2016.R | Zhang and Yuan (2016) | 398 | 0 | 0 | 426 |
| Waterfall | compare-waterfall-with-BOIN.R | ‘BOIN’ | 870 | 0 | 0 | 870 |
| BOIN-ET | reproduce_Takeda2018.R | Takeda et al. (2018) | 210 | 0 | 14 | 0 |
| TITE-BOIN-ET | reproduce_Takeda2020.R | Takeda et al. (2020) | 295 | 0 | 29 | 0 |
| BOIN-ET, TITE-BOIN-ET | compare-with-boinet.R | ‘boinet’ | 508 | 0 | 0 | 0 |
| BF-BOIN-ET | reproduce_Takeda2025.R | Takeda et al. (2025) | 336 | 0 | 0 | 0 |
| BF-BOIN-ET | compare-with-bfboinet.R | ‘bfboinet’ | 85 | 0 | 0 | 0 |
| BOIN-ETC | reproduce_Kakizume2024.R | Kakizume et al. (2024) | 1168 | 8 | 0 | 996 |
| BOIN-ETC | compare-with-BOINETC.R | Code of Kakizume et al. (2024) | 1801 | 35 | 0 | 144 |
| gBOIN-ET | reproduce_Takeda2022.R | Takeda et al. (2022) | 131 | 4 | 7 | 260 |
| gBOIN-ET | compare-with-gboinet.R | ‘boinet’ | 156 | 0 | 0 | 156 |
| TITE-gBOIN-ET | reproduce_Takeda2023.R | Takeda et al. (2023) | 846 | 297 | 0 | 190 |
| TITE-gBOIN-ET | compare-with-tite-gboinet.R | ‘boinet’ | 168 | 0 | 0 | 168 |
| BOIN12 | reproduce_Lin2020.R | Lin et al. (2020) | 242 | 1 | 0 | 66 |
| BOIN12 | compare-with-escalation.R | ‘escalation’ | 93 | 6 | 0 | 0 |
| TITE-BOIN12 | reproduce_Zhou2022.R | Zhou et al. (2022) | 121 | 4 | 1 | 162 |
| Comb-BOIN12, TITE-Comb-BOIN12 | reproduce_Lu2025.R | Lu et al. (2025) | 1508 | 0 | 56 | 190 |
| U-BOIN | reproduce_Zhou2019.R | Zhou, Lee and Yuan (2019) | 1340 | 36 | 0 | 0 |

Of the 11,456 values judged by the 22 scripts, 10,930 pass, 399 are
explained and 127 fail, and 4,146 more are listed for information. The
sections below describe what each script compares.

## Designs that select the MTD

**iBOIN.** `reproduce_Zhou2021.R` compares
[`iboin_boundary()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/iboin_boundary.md)
with the boundaries of Table 1 of Zhou et al. (2021), exactly, and the
simulations with the rows BOIN, iBOIN and iBOINR (the robust prior) of
Table S.1 of its supporting information, five operating characteristics
in ten scenarios, within four standard errors plus 0.05. The article
does not say whether its simulations select the MTD with or without the
prior information; they are judged with it, and the selection from the
data alone is INFO.

| Part                                        | PASS | EXPLAINED | FAIL | INFO |
|:--------------------------------------------|-----:|----------:|-----:|-----:|
| Table 1                                     |  100 |         0 |    0 |    0 |
| Table S.1, BOIN                             |   48 |         2 |    0 |    0 |
| Table S.1, iBOIN                            |   48 |         2 |    0 |    0 |
| Table S.1, iBOINR                           |   48 |         2 |    0 |    0 |
| Table S.1, iBOIN (MTD from the data alone)  |    0 |         0 |    0 |   50 |
| Table S.1, iBOINR (MTD from the data alone) |    0 |         0 |    0 |   50 |

**gBOIN.** `reproduce_Mu2019.R` compares the boundaries of Table 2 of Mu
et al. (2019) and the selection percentages and numbers of patients of
Table 5 under the elimination rule of the article
(`elim_rule = "posterior"`), within four standard errors plus the
rounding; the rule of the ‘UnifiedDoseFinding’ package
(`elim_rule = "count"`) is INFO. `compare-with-UnifiedDoseFinding.R`
compares the same scenarios with that package, judging the rows with its
elimination rule.

| Part                                  | PASS | EXPLAINED | FAIL | INFO |
|:--------------------------------------|-----:|----------:|-----:|-----:|
| Table 2                               |   14 |         0 |    0 |    0 |
| Table 5, gBOIN                        |  101 |         2 |   17 |    0 |
| Table 5, gBOIN (elimination by count) |    0 |         0 |    0 |  120 |

| Part                  | PASS | EXPLAINED | FAIL | INFO |
|:----------------------|-----:|----------:|-----:|-----:|
| elim_rule = count     |  130 |         0 |    0 |    0 |
| elim_rule = posterior |    0 |         0 |    0 |  130 |

**Accelerated BOIN.** `reproduce_Kojima2024.R` compares the accelerated
BOIN design (`ABOIN` in the table, `titration_doses = 2` and
`n_earlystop = 6`) and the BOIN design with Tables 1 and 2 of Kojima et
al. (2024), within four standard errors plus the rounding. The BOIN rows
of the article agree with no early stop (`n_earlystop = 100`), although
the text states 12 patients; the rows with 12 patients and the 3+3 rows,
computed in closed form by
[`oc_3p3()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_3p3.md),
are INFO.

| Part                                 | PASS | EXPLAINED | FAIL | INFO |
|:-------------------------------------|-----:|----------:|-----:|-----:|
| Target 0.30, ABOIN                   |   42 |         0 |    0 |    0 |
| Target 0.30, BOIN                    |   42 |         0 |    0 |    0 |
| Target 0.30, BOIN (n_earlystop = 12) |    0 |         0 |    0 |   42 |
| Target 0.30, 3+3                     |    0 |         0 |    0 |   42 |
| Target 0.25, ABOIN                   |   41 |         0 |    1 |    0 |
| Target 0.25, BOIN                    |   40 |         0 |    2 |    0 |
| Target 0.25, BOIN (n_earlystop = 12) |    0 |         0 |    0 |   42 |
| Target 0.25, 3+3                     |    0 |         0 |    0 |   42 |

**BOIN, TITE-BOIN, BF-BOIN and BE-BOIN.** `compare-with-BOIN.R` compares
[`sim_boin()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin.md)
with [`BOIN::get.oc()`](https://rdrr.io/pkg/BOIN/man/get.oc.html) over a
grid of configurations, by equality, since the two draw the same trials
for the same seed, and `compare-with-Zhao2024.R` compares BF-BOIN with
Table 4 of Zhao et al. (2024). These scripts print their results and do
not keep a results file in the package; the tests check the agreement
with the ‘BOIN’ package on smaller settings. The decision tables of
TITE-BOIN are reproduced by the tests entry by entry.

## Drug combinations

**BOIN drug-combination design.** `compare-comb-with-BOIN.R` compares
[`sim_comb_boin()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_comb_boin.md)
with
[`BOIN::get.oc.comb()`](https://rdrr.io/pkg/BOIN/man/get.oc.comb.html)
by equality, and `compare-with-LinYin2017.R` sets the results beside
Tables 2 and 3 of Lin and Yin (2017), as INFO, since the package follows
the ‘BOIN’ package where it differs from the article. Neither keeps a
results file.

**Waterfall.** `reproduce_Zhang2016.R` compares
[`sim_waterfall()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_waterfall.md)
with Tables IV, V and S1 of Zhang and Yuan (2016) with
`boin_compatible = TRUE`, the seed 6 and 1,000 trials, by equality with
the printed digits: the article was produced with these settings. The
default rules of the package, which keep the data of a subtrial that
ends without a candidate MTD, and the sample sizes of Table V are INFO.
`compare-waterfall-with-BOIN.R` compares every value of the output of
`BOIN::get.oc.comb(mtd.contour = TRUE)` in 20 settings, again by
equality.

| Part                                       | PASS | EXPLAINED | FAIL | INFO |
|:-------------------------------------------|-----:|----------:|-----:|-----:|
| Table IV, Waterfall                        |  178 |         0 |    0 |    0 |
| Table S1, Waterfall                        |  178 |         0 |    0 |    0 |
| Table V, Waterfall                         |   42 |         0 |    0 |   14 |
| Table IV, Waterfall (rules of the package) |    0 |         0 |    0 |  178 |
| Table S1, Waterfall (rules of the package) |    0 |         0 |    0 |  178 |
| Table V, Waterfall (rules of the package)  |    0 |         0 |    0 |   56 |

| Part                             | PASS | EXPLAINED | FAIL | INFO |
|:---------------------------------|-----:|----------:|-----:|-----:|
| Waterfall                        |  870 |         0 |    0 |    0 |
| Waterfall (rules of the package) |    0 |         0 |    0 |  870 |

## The BOIN-ET family

**BOIN-ET and TITE-BOIN-ET.** `reproduce_Takeda2018.R` compares the
boundaries and the decision tables of Table 2 of Takeda et al. (2018),
exactly, and the operating characteristics of Table 3, within four
standard errors plus 0.5 for the rounding. `reproduce_Takeda2020.R` does
the same for Table 2 of Takeda et al. (2020), whose simulations did not
wait for finished assessments (`min_completed_share = 0`).
`compare-with-boinet.R` compares
[`sim_boin_et()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_et.md)
and
[`sim_tite_boin_et()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_tite_boin_et.md)
with `boinet_compatible = TRUE` with
[`boinet::boinet()`](https://rdrr.io/pkg/boinet/man/boinet.html) and
[`boinet::tite.boinet()`](https://rdrr.io/pkg/boinet/man/tite.boinet.html),
within four standard errors plus 0.05.

| Part       | PASS | EXPLAINED | FAIL | INFO |
|:-----------|-----:|----------:|-----:|-----:|
| Boundaries |    3 |         0 |    0 |    0 |
| Table 2    |   65 |         0 |    0 |    0 |
| Table 3    |  142 |         0 |   14 |    0 |

| Part                  | PASS | EXPLAINED | FAIL | INFO |
|:----------------------|-----:|----------:|-----:|-----:|
| Table 2, TITE-BOIN-ET |  151 |         0 |   17 |    0 |
| Table 2, BOIN-ET      |  144 |         0 |   12 |    0 |

| Part         | PASS | EXPLAINED | FAIL | INFO |
|:-------------|-----:|----------:|-----:|-----:|
| BOIN-ET      |  312 |         0 |    0 |    0 |
| TITE-BOIN-ET |  196 |         0 |    0 |    0 |

Of the values that fail, listed at the end of this article, 10 of the 14
of Takeda et al. (2018) and 17 of the 29 of Takeda et al. (2020) are in
scenario 11, whose printed response rate of 0.55 at the lowest dose does
not fit the published results; the ‘boinet’ package too selects that
dose far more often than the article reports. In Takeda et al. (2020)
the early terminations, durations and allocations of TITE-BOIN-ET in
scenarios 8, 9, 10 and 12 agree better with an exclusion of doses that
uses only the patients whose outcomes are known, which the article does
not state. The other values that fail are not explained.

**BF-BOIN-ET.** `reproduce_Takeda2025.R` compares BF-BOIN-ET-C and
BF-BOIN-ET-A with Tables 1 and 2 of Takeda et al. (2025), within four
standard errors plus 0.05, and `compare-with-bfboinet.R` compares
[`sim_bf_boin_et()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_bf_boin_et.md)
with
[`bfboinet::get.oc.backboinet()`](https://rdrr.io/pkg/bfboinet/man/get.oc.backboinet.html),
within four standard errors.

| Part                            | PASS | EXPLAINED | FAIL | INFO |
|:--------------------------------|-----:|----------:|-----:|-----:|
| Table 1 (30 days), BF-BOIN-ET-C |   84 |         0 |    0 |    0 |
| Table 1 (30 days), BF-BOIN-ET-A |   84 |         0 |    0 |    0 |
| Table 2 (60 days), BF-BOIN-ET-C |   84 |         0 |    0 |    0 |
| Table 2 (60 days), BF-BOIN-ET-A |   84 |         0 |    0 |    0 |

| Part                                         | PASS | EXPLAINED | FAIL | INFO |
|:---------------------------------------------|-----:|----------:|-----:|-----:|
| Example of bfboinet                          |   17 |         0 |    0 |    0 |
| Plateau                                      |   17 |         0 |    0 |    0 |
| Toxic                                        |   14 |         0 |    0 |    0 |
| Futile high doses                            |   20 |         0 |    0 |    0 |
| Plateau (uniform, exponential, suspend 0.33) |   17 |         0 |    0 |    0 |

**BOIN-ETC.** `reproduce_Kakizume2024.R` compares the allocation of the
three versions of BOIN-ETC with Table 3 and Tables S1 and S4 to S9 of
Kakizume et al. (2024), within four standard errors plus the rounding,
with the rules of the code of the authors for BOIN-ETC1
(`etc1_start = "below"`, `rule4 = "authors_code"`). The selection
percentages are INFO: the code estimates the DLT rates with
[`Iso::biviso()`](https://rdrr.io/pkg/Iso/man/biviso.html), whose
convergence decides between combinations that tie in the exact isotonic
regression of the package. `compare-with-BOINETC.R` runs the R code of
the authors with the random numbers of R after `set.seed(1234)`, which
gives the published values to the printed digits, and judges the
selection of the package on the same trials as INFO.

| Part                      | PASS | EXPLAINED | FAIL | INFO |
|:--------------------------|-----:|----------:|-----:|-----:|
| Kakizume 2024 SI Table S2 |    0 |         0 |    0 |  408 |
| Kakizume 2024 SI Table S4 |  406 |         2 |    0 |    0 |
| Kakizume 2024 Table 3     |   60 |         0 |    0 |   30 |
| Kakizume 2024 SI Table S6 |   58 |         2 |    0 |   30 |
| Kakizume 2024 SI Table S7 |   60 |         0 |    0 |   30 |
| Kakizume 2024 SI Table S3 |    0 |         0 |    0 |  408 |
| Kakizume 2024 SI Table S5 |  408 |         0 |    0 |    0 |
| Kakizume 2024 SI Table S1 |   60 |         0 |    0 |   30 |
| Kakizume 2024 SI Table S8 |   59 |         1 |    0 |   30 |
| Kakizume 2024 SI Table S9 |   57 |         3 |    0 |   30 |

| Part                                  | PASS | EXPLAINED | FAIL | INFO |
|:--------------------------------------|-----:|----------:|-----:|-----:|
| Kakizume 2024 SI Table S2             |  315 |         0 |    0 |    0 |
| Kakizume 2024 SI Table S4             |  314 |         1 |    0 |    0 |
| Kakizume 2024 Table 3                 |   72 |         0 |    0 |    0 |
| Copy of the selection of CALC.SUM()   |  144 |         0 |    0 |    0 |
| simFastBOIN on the trials of the code |    0 |         0 |    0 |  144 |
| Kakizume 2024 SI Table S6             |   68 |         4 |    0 |    0 |
| Kakizume 2024 SI Table S7             |   72 |         0 |    0 |    0 |
| Kakizume 2024 SI Table S3             |  295 |        20 |    0 |    0 |
| Kakizume 2024 SI Table S5             |  311 |         4 |    0 |    0 |
| Kakizume 2024 SI Table S1             |   72 |         0 |    0 |    0 |
| Kakizume 2024 SI Table S8             |   69 |         3 |    0 |    0 |
| Kakizume 2024 SI Table S9             |   69 |         3 |    0 |    0 |

**gBOIN-ET.** `reproduce_Takeda2022.R` compares the boundaries printed
in Takeda et al. (2022) and Table 1. The SAS program of the authors
draws one Bernoulli variate for every category and takes the highest
that comes up, so that the categories are not drawn with the
probabilities of Table S1; the judged rows simulate the probabilities
that the program draws, with the printed boundaries, the linear
fractional polynomial, `min_obd_sample = 6`, independent outcomes and
`boinet_compatible = TRUE`, within four standard errors plus half a unit
of the printed digit. The rows with the probabilities of Table S1 and
with the default rules are INFO. `compare-with-gboinet.R` compares
[`sim_gboin_et()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_gboin_et.md)
with [`boinet::gboinet()`](https://rdrr.io/pkg/boinet/man/gboinet.html)
under four selection rules, within four standard errors plus 0.05.

| Part                                          | PASS | EXPLAINED | FAIL | INFO |
|:----------------------------------------------|-----:|----------:|-----:|-----:|
| Boundaries                                    |   12 |         0 |    0 |    0 |
| Table 1, gBOIN-ET                             |  119 |         4 |    7 |    0 |
| Table 1, gBOIN-ET (rules of the package)      |    0 |         0 |    0 |  130 |
| Table 1, gBOIN-ET (probabilities of Table S1) |    0 |         0 |    0 |  130 |

| Part                            | PASS | EXPLAINED | FAIL | INFO |
|:--------------------------------|-----:|----------:|-----:|-----:|
| gBOIN-ET                        |  156 |         0 |    0 |    0 |
| gBOIN-ET (rules of the package) |    0 |         0 |    0 |  156 |

**TITE-gBOIN-ET.** `reproduce_Takeda2023.R` compares the boundaries
printed in Takeda et al. (2023) and Table 1 and Tables S2 to S6. The R
program of its supporting information, written again in C++ with the
random numbers of R, reproduces the summaries that the supporting
information holds for the program exactly. It differs from the article
and from ‘boinet’ in two points: where one dose has the highest mean
efficacy score among the current dose and its neighbors, it passes that
dose to [`sample()`](https://rdrr.io/r/base/sample.html), which then
draws from all doses up to it, and it excludes doses with the effective
numbers of patients. The published cells are judged with the rules of
‘boinet’, and those that these two points move out of the tolerance are
explained; the same simulations are also judged against the program
without the two points. `compare-with-tite-gboinet.R` compares
[`sim_tite_gboin_et()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_tite_gboin_et.md)
with
[`boinet::tite.gboinet()`](https://rdrr.io/pkg/boinet/man/tite.gboinet.html),
within four standard errors plus 0.05.

| Part | PASS | EXPLAINED | FAIL | INFO |
|:---|---:|---:|---:|---:|
| Boundaries | 3 | 0 | 0 | 0 |
| Table 1, TITE-gBOIN-ET | 28 | 32 | 0 | 0 |
| Tables S2 to S6, TITE-gBOIN-ET | 245 | 265 | 0 | 0 |
| Table 1, TITE-gBOIN-ET against the corrected program | 60 | 0 | 0 | 0 |
| Tables S2 to S6, TITE-gBOIN-ET against the corrected program | 510 | 0 | 0 | 0 |
| Table 1, TITE-gBOIN-ET (rules of the package) | 0 | 0 | 0 | 60 |
| Tables S2 to S6, TITE-gBOIN-ET (rules of the package) | 0 | 0 | 0 | 130 |

| Part                                 | PASS | EXPLAINED | FAIL | INFO |
|:-------------------------------------|-----:|----------:|-----:|-----:|
| TITE-gBOIN-ET                        |  168 |         0 |    0 |    0 |
| TITE-gBOIN-ET (rules of the package) |    0 |         0 |    0 |  168 |

## Utility-based designs

**BOIN12.** `reproduce_Lin2020.R` compares the desirability scores of
Table 3 of Lin et al. (2020), the boundaries of Figure 1, the example of
Section 6 of its supplement and a few values of Zhou et al. (2022), to
the printed digits, and the simulations of Figure 3 and Figures S3 and
S4 of the supplement, within three standard errors plus 0.05, with the
rules of the simulations of the article (`tie_break = "random"`,
`no_candidate = "stop"`, `explore_from = "admissible"`); the default
rules are INFO. `compare-with-escalation.R` compares
[`sim_boin12()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin12.md)
with `escalation::get_boin12()` in nine scenarios, within four standard
errors; the values that differ are explained by the rules of
‘escalation’ when a dose becomes inadmissible.

| Part                    | PASS | EXPLAINED | FAIL | INFO |
|:------------------------|-----:|----------:|-----:|-----:|
| Lin 2020 Table 3        |  166 |         0 |    0 |    0 |
| Lin 2020 Figure 1       |    2 |         0 |    0 |    0 |
| Zhou 2022 Section 2.2.4 |    1 |         0 |    0 |    0 |
| Lin 2020 SI Section 6   |    6 |         0 |    0 |    0 |
| Zhou 2022 Table S.3     |    5 |         1 |    0 |    0 |
| Lin 2020 Figure 3       |   36 |         0 |    0 |   36 |
| Lin 2020 SI Figure S3   |   14 |         0 |    0 |   14 |
| Lin 2020 SI Figure S4   |   12 |         0 |    0 |   16 |

| Part           | PASS | EXPLAINED | FAIL | INFO |
|:---------------|-----:|----------:|-----:|-----:|
| Nine scenarios |   93 |         6 |    0 |    0 |

**TITE-BOIN12.** `reproduce_Zhou2022.R` compares the desirabilities and
estimates of Table S.3 of the supplement of Zhou et al. (2022) with
`pending_estimate = "marginal"`, to the printed digits, and the
simulations of Figure 3, within three standard errors plus half a unit
of the printed digit. The values of Equation (11) of the article and the
trial durations, whose accrual the article does not describe, are INFO.

| Part                              | PASS | EXPLAINED | FAIL | INFO |
|:----------------------------------|-----:|----------:|-----:|-----:|
| Zhou 2022 Table S.3               |   50 |         4 |    0 |    0 |
| Zhou 2022 Table S.3 (Equation 11) |    0 |         0 |    0 |   54 |
| Zhou 2022 Figure 3                |   71 |         0 |    1 |  108 |

**Comb-BOIN12 and TITE-Comb-BOIN12.** `reproduce_Lu2025.R` compares the
utilities of Tables 1 and 2 and Tables S2 and S3 of Lu et al. (2025) and
the tail probabilities of Figure S1, to the printed digits, and the
simulations of Figure 2, Figure S2 and Tables 3 and 4, within three
standard errors plus 0.05, with `obd_bound = "contour"` and, for
TITE-Comb-BOIN12, `pending_use = "utility"`. The selection percentages
of Tables 3 and 4, which differ from those of Figure 2 for the same
trials, and the durations are INFO. The values that fail are bars of
TITE-Comb-BOIN12, whose accrual and waiting rules the article does not
describe, and selection percentages of Comb-BOIN12, for a reason that
was not found.

| Part                 | PASS | EXPLAINED | FAIL | INFO |
|:---------------------|-----:|----------:|-----:|-----:|
| Lu 2025 Table 1      |  153 |         0 |    0 |    0 |
| Lu 2025 SI Table S2  |  288 |         0 |    0 |    0 |
| Lu 2025 Table 2      |  130 |         0 |    0 |    0 |
| Lu 2025 SI Table S3  |  240 |         0 |    0 |    0 |
| Lu 2025 SI Figure S1 |  164 |         0 |    0 |    0 |
| Lu 2025 Table 3      |   27 |         0 |    0 |   54 |
| Lu 2025 SI Figure S2 |  217 |         0 |   11 |    0 |
| Lu 2025 Figure 2     |  259 |         0 |   45 |   76 |
| Lu 2025 Table 4      |   30 |         0 |    0 |   60 |

**U-BOIN.** There is no reference implementation. `reproduce_Zhou2019.R`
compares the boundaries of Table 2 and Table S2 of Zhou, Lee and Yuan
(2019), the utility tables S4 to S6 of its supporting information and
the true utilities, to the printed digits, and the simulations of Table
3 and Tables S8 and S9, within four standard errors plus half a unit of
the printed digit, with the rates of Scenario 7 of Table S8 for Scenario
7 of Table 3 and with s1 = 12 for Table S8, with which the published
results agree. The admissibility of Tables S4 to S6 is judged with the
default prior and with `prior_resp = c(0.05, 0.05)`, which reproduces
it, and the simulations with the printed rates and with the s1 = 9 of
the text are judged as well.

| Part                                              | PASS | EXPLAINED | FAIL | INFO |
|:--------------------------------------------------|-----:|----------:|-----:|-----:|
| Table 2, BOIN                                     |   11 |         1 |    0 |    0 |
| Table S2, U-BOIN                                  |   72 |         0 |    0 |    0 |
| Table S4, U-BOIN, efficacy prior Beta(0.05, 0.05) |   20 |         0 |    0 |    0 |
| Table S4, U-BOIN                                  |   30 |         3 |    0 |    0 |
| Table S5, U-BOIN, efficacy prior Beta(0.05, 0.05) |   84 |         0 |    0 |    0 |
| Table S5, U-BOIN                                  |  126 |         4 |    0 |    0 |
| Table S6, U-BOIN, efficacy prior Beta(0.05, 0.05) |  219 |         0 |    0 |    0 |
| Table S6, U-BOIN                                  |  318 |         5 |    0 |    0 |
| Table 3, U-BOIN                                   |  128 |         0 |    0 |    0 |
| Table 3, U-BOIN, rates of Table 3 as printed      |    4 |         7 |    0 |    0 |
| Table S9, U-BOIN                                  |  128 |         0 |    0 |    0 |
| Table S8, U-BOIN                                  |  128 |         0 |    0 |    0 |
| Table S8, U-BOIN, s1 = 9 as in the text           |   72 |        16 |    0 |    0 |

## Values that are explained

The rows of `inst/validation/explained.csv`, grouped by script and
reason.

| Script | Values | Reason |
|:---|---:|:---|
| reproduce_Zhou2021.R | 6 | The risk of overdosing of the same row, which also counts the patients above the MTD, is reproduced, and the published value agrees with the percentage of patients treated at doses whose true DLT probability exceeds 0.40 (column alternative of reproduce_Zhou2021.csv: about 6.4 and 5.9 for BOIN, 3.8 and 0.4 for iBOIN and iBOINR in Scenarios 2 and 10, against published 6.0, 5.4, 3.6 and 0.6). In the other eight scenarios the doses above the MTD are the doses above 0.40, and the values are reproduced |
| reproduce_Mu2019.R | 1 | The published average of 2.7 patients at dose 1 is impossible: every trial starts at dose 1 with a cohort of three, so that the average is at least 3 (reproduced: about 3.65). The published row sums to 30.0, so that another cell of the row must be off as well; doses 1 and 4 together hold 11.2 patients in the article and about 11.1 here |
| reproduce_Mu2019.R | 1 | The published row of selection proportions sums to 0.91 (0.91 at dose 2 and 0.00 elsewhere) although 4.1 patients are treated at dose 3 on average and no trial is reported as stopped early; the reproduced value is about 9.9, which completes the row to 100. The printed 0.00 appears to be a misprint of 0.09 |
| reproduce_Kakizume2024.R | 2 | Tables S8 and S9 apply other utility scores to the trials of Table S1, so that their %nODC and \#N are those of Table S1: Table S1 and Table S8 give a %nODC of 39.8 and a \#N of 54.6 for this setting. Table S9 prints 54.6 as %nODC and, as \#N, 38.6, the \#N of the same scenario and design with n_stop = 6 (Table 3). |
| reproduce_Kakizume2024.R | 1 | Tables S6 and S7 apply other utility scores to the trials of Table 3, so that their %nODC and \#N are those of Table 3: Table 3 and Table S7 give 17.8 for this setting. The value 26.1 is the \#N of the same scenario and design with n_stop = 9 (Table S1). |
| reproduce_Kakizume2024.R | 1 | Tables S6 and S7 apply other utility scores to the trials of Table 3, so that their %nODC and \#N are those of Table 3: Table 3 and Table S7 give 22.3 for this setting. The value 32.1 is the \#N of the same scenario and design with n_stop = 9 (Table S1). |
| reproduce_Kakizume2024.R | 1 | Tables S8 and S9 apply other utility scores to the trials of Table S1, so that their %nODC and \#N are those of Table S1: Table S1 and Table S8 give 22.9 for this setting, and 222.9 is not a percentage. |
| reproduce_Kakizume2024.R | 1 | Tables S8 and S9 apply other utility scores to the trials of Table S1, so that their %nODC and \#N are those of Table S1: Table S1 and Table S9 give 60.3 for this setting. The row of Table S8 (60.3 / 58.3 / 26.4) appears shifted: 58.3 is the %nODC of the same scenario and design with n_stop = 6 (Table 3). |
| reproduce_Kakizume2024.R | 1 | The mean numbers of patients of this row of Table S4 add up to 40.24, whereas \#N of Table 3 is 39.2 for the same trials; with 2.81 at A3B3 they add up to 39.24. |
| reproduce_Kakizume2024.R | 1 | The mean numbers of patients of this row of Table S4 add up to 45.61, whereas \#N of Table 3 is 45.3 for the same trials; with 1.14 at A2B5 they add up to 45.31. A translation of the code of the authors into Python, fed with the random numbers of R after set.seed(1234), gives 1.137 at A2B5 and every other value of the run to the printed digits. |
| compare-with-BOINETC.R | 10 | In Table S3, levels 1 and 2 of drug A of BOIN-ETC2 for this scenario are the selection percentages of BOIN-ETC3 of the code of the authors (all 10 agree to the printed digits); level 3 and the %ODC of Table S1 (80.6) agree with BOIN-ETC2 of the code. |
| compare-with-BOINETC.R | 10 | In Table S3, levels 1 and 2 of drug A of BOIN-ETC3 for this scenario agree neither with BOIN-ETC2 nor with BOIN-ETC3 of the code of the authors, whose BOIN-ETC3 values for them are printed in the row of BOIN-ETC2; level 3 and the %ODC of Table S1 (77.9) agree with BOIN-ETC3 of the code. |
| compare-with-BOINETC.R | 2 | Tables S8 and S9 apply other utility scores to the trials of Table S1, so that their %nODC and \#N are those of Table S1: Table S1 and Table S8 give a %nODC of 39.8 and a \#N of 54.6 for this setting. Table S9 prints 54.6 as %nODC and, as \#N, 38.6, the \#N of the same scenario and design with n_stop = 6 (Table 3). |
| compare-with-BOINETC.R | 1 | 80.3 is the %ODC of Table S1 (utility 100-40-60-0) for the same trials. The code of the authors gives 78.8 for this utility, and its %ODC for the other two utilities agrees with Tables S1 and S9. |
| compare-with-BOINETC.R | 1 | A translation of the code of the authors into Python, fed with the random numbers of R after set.seed(1234), gives 3.036 at A2B2 and every other value of the run to the printed digits except A2B1, whose published value is a misprint (the row then adds up to \#N of Table S1); the published 3.05 appears to be a misprint as well. |
| compare-with-BOINETC.R | 1 | A translation of the code of the authors into Python, fed with the random numbers of R after set.seed(1234), gives 6.762 at A1B5 and every other value of the run to the printed digits, including \#N and %nODC of Table S1; the published 6.79 appears to be a misprint. |
| compare-with-BOINETC.R | 1 | A translation of the code of the authors into Python, fed with the random numbers of R after set.seed(1234), gives 7.317 at A1B4 and every other value of the run to the printed digits, including \#N and %nODC of Table S1; the published 7.31 appears to be a misprint. |
| compare-with-BOINETC.R | 1 | Tables S6 and S7 apply other utility scores to the trials of Table 3, so that their %nODC and \#N are those of Table 3: Table 3 and Table S7 give 17.8 for this setting. The value 26.1 is the \#N of the same scenario and design with n_stop = 9 (Table S1). |
| compare-with-BOINETC.R | 1 | Tables S6 and S7 apply other utility scores to the trials of Table 3, so that their %nODC and \#N are those of Table 3: Table 3 and Table S7 give 22.3 for this setting. The value 32.1 is the \#N of the same scenario and design with n_stop = 9 (Table S1). |
| compare-with-BOINETC.R | 1 | Tables S8 and S9 apply other utility scores to the trials of Table S1, so that their %nODC and \#N are those of Table S1: Table S1 and Table S8 give 22.9 for this setting, and 222.9 is not a percentage. |
| compare-with-BOINETC.R | 1 | Tables S8 and S9 apply other utility scores to the trials of Table S1, so that their %nODC and \#N are those of Table S1: Table S1 and Table S9 give 60.3 for this setting. The row of Table S8 (60.3 / 58.3 / 26.4) appears shifted: 58.3 is the %nODC of the same scenario and design with n_stop = 6 (Table 3). |
| compare-with-BOINETC.R | 1 | The mean numbers of patients of this row of Table S4 add up to 45.61, whereas \#N of Table 3 is 45.3 for the same trials; with 1.14 at A2B5 they add up to 45.31. A translation of the code of the authors into Python, fed with the random numbers of R after set.seed(1234), gives 1.137 at A2B5 and every other value of the run to the printed digits. |
| compare-with-BOINETC.R | 1 | The mean numbers of patients of this row of Table S5 add up to 26.21, whereas \#N of Table S1 is 26.1 for the same trials; with the values of the code at A2B1 (7.66) and A2B2 (3.04) they add up to 26.10. A translation of the code of the authors into Python, fed with the random numbers of R after set.seed(1234), gives 7.656 at A2B1 and every other value of the run to the printed digits except A2B2. |
| compare-with-BOINETC.R | 1 | The row of Table S6 (17.8 / 58.3 / 18.0) appears shifted: 17.8 is the \#N of the same scenario with BOIN-ETC2 and n_stop = 6 (Table 3), whereas the code of the authors gives a %ODC of 81.4 for these trials and this utility, and its %ODC for the other two utilities agrees with Tables 3 and S7 (81.4 and 90.8). |
| compare-with-BOINETC.R | 1 | The row of Table S6 (80.3 / 59.6 / 26.1) takes 80.3 and 26.1 from Table S1 (n_stop = 9, utility 100-40-60-0) for the same scenario and design. The code of the authors gives a %ODC of 78.7 for these trials and this utility, and its %ODC for the other two utilities agrees with Tables 3 and S7. |
| compare-with-BOINETC.R | 1 | The row of Table S8 (60.3 / 58.3 / 26.4) appears shifted: 60.3 is the %nODC of these trials (Tables S1 and S9), whereas the code of the authors gives a %ODC of 79.3 for these trials and this utility, and its %ODC for the other two utilities agrees with Tables S1 and S9. |
| reproduce_Takeda2022.R | 2 | The SAS program of the authors sorts the 36 fractional polynomials by the first value of the fit table of PROC GENMOD, the deviance, in descending order and takes the first, the model with the largest deviance; the article and the package take the best-fitting model. A transcription of the program in Python gives about 5.6 for dose 2 and 13.2 for dose 6 with the largest deviance, within the allowance (published 4 and 14), and about 6.2 and 10.7 with the smallest, as the package (about 6.5 and 10.3) |
| reproduce_Takeda2022.R | 1 | The difference is at the edge of the allowance (about 2.5 against a published 1, allowance about 1.5). A transcription of the SAS program of the authors in Python gives about 1.9 (10,000 trials), within the allowance. The program differs from the package in the points listed in inst/reproduce/reproduce_Takeda2022.R (no posterior probability rules during the trial or at the selection, the fractional polynomial with the largest deviance); which of them accounts for this cell was not isolated |
| reproduce_Takeda2022.R | 1 | The SAS program of the authors has no posterior probability rules. It ends a trial at the lowest dose when the design would stay or de-escalate there and the mean efficacy score is at most half the target or the toxicity score reaches a count boundary, and it recommends no dose when no dose up to the MTD has six patients; the package ends a trial when the posterior probability rules of the article exclude every dose, which happens less often here (about 0.6 against a published 2). A transcription of the program in Python gives about 1.7 (10,000 trials), within the allowance |
| reproduce_Takeda2023.R | 297 | The program of the supporting information of the article passes the dose with the highest mean efficacy score to sample(), which draws from all doses up to it when it is a single dose, and excludes doses with the effective numbers of patients; neither is in the article or in boinet::tite.gboinet(). A transcription of the program reproduces all its published summaries exactly, and with these two points removed it gives the value of the package within simulation error (row ‘TITE-gBOIN-ET against the corrected program’) |
| reproduce_Lin2020.R | 1 | Six patients with two DLTs and three responses give 0.191. The same row of Table S.3 reports a response rate of 0.50 for dose 4, whose three patients had no response, so the row appears to hold a misprint |
| compare-with-escalation.R | 6 | Doses 1 and 2 have low response rates and become inadmissible. The rules of ‘escalation’ (admissibility of the neighbouring doses only, de-escalation to the dose below whatever its data, end of the trial when no neighbouring dose is admissible, OBD selected without the admissibility check) were written again in Python and give the values of ‘escalation’ within simulation error, and the rules of simFastBOIN written in Python give those of sim_boin12() (8000 trials each; Scenario 2: DL1 20.8 and 14.8, DL3 19.1 and 27.3, no OBD 13.0 and 4.7; Scenario 4: no OBD 1.3 and 0.4) |
| reproduce_Zhou2022.R | 3 | Table 2 gives 0.118 when a pending outcome is replaced by pi_hat and 0.117 by Equation (11), with pi_hat = (0.43, 0.44). The table reports pi_hat = (0.37, 0.55), but at day 315 two responses are observed and two patients have completed the response window without response at dose 3, so that pi_E is at most 0.50 whatever the convention; the row was computed from data other than Table 2 |
| reproduce_Zhou2022.R | 1 | The three patients at dose 4 had no response, so pi_E is 0. The value 0.56 repeats that of dose 3, as in the BOIN12 row of the same analysis, and appears to be a misprint |
| reproduce_Zhou2019.R | 16 | The supporting information gives s1 = 9 for Table S8, but its U-BOIN-CD rows agree with s1 = 12, which the rows of the design U-BOIN use. |
| reproduce_Zhou2019.R | 12 | Tables S4 to S6 mark a dose without any response inadmissible, which the efficacy prior Beta(0.05, 0.05) reproduces (rows U-BOIN, efficacy prior Beta(0.05, 0.05)). The simulations of the article agree with the default Beta(1, 1) prior, under which these counts are admissible. |
| reproduce_Zhou2019.R | 7 | The results of Scenario 7 of Table 3 agree with the rates of Scenario 7 of Table S8 (DLT 0.03, 0.08, 0.25, 0.40, 0.55; efficacy 0.05, 0.25, 0.35, 0.40, 0.45), not with the rates printed in Table 3; the rows of the design U-BOIN use the rates of Table S8. |
| reproduce_Zhou2019.R | 1 | The article prints 0.358, the value common in the BOIN literature; lambda_d for the target 0.30 is 0.35852, which rounds to 0.359. |

## Values that fail

The 127 values that fail, by script. For Kojima et al. (2024), Mu et
al. (2019), Takeda et al. (2022) and Zhou et al. (2022) no reason was
found; for Takeda et al. (2018) and (2020) and Lu et al. (2025) the
sections above describe them.

| Script | Value | Published | Reproduced | Tolerance |
|:---|:---|:---|:---|:---|
| reproduce_Mu2019.R | gBOIN, Scenario 1, Patients DL3 | 8.4 | 9.358 | 0.495 |
| reproduce_Mu2019.R | gBOIN, Scenario 1, Patients DL5 | 3.7 | 2.872 | 0.345 |
| reproduce_Mu2019.R | gBOIN, Scenario 2, Selection (%) DL2 | 0 | 1.26 | 1.092 |
| reproduce_Mu2019.R | gBOIN, Scenario 3, Selection (%) DL1 | 2 | 4.19 | 1.796 |
| reproduce_Mu2019.R | gBOIN, Scenario 3, Selection (%) DL5 | 0 | 1.36 | 1.115 |
| reproduce_Mu2019.R | gBOIN, Scenario 4, Patients DL2 | 9 | 9.961 | 0.571 |
| reproduce_Mu2019.R | gBOIN, Scenario 4, Patients DL3 | 10.8 | 9.492 | 0.543 |
| reproduce_Mu2019.R | gBOIN, Scenario 6, Selection (%) DL2 | 29 | 23.13 | 3.785 |
| reproduce_Mu2019.R | gBOIN, Scenario 6, Patients DL3 | 1.7 | 1.402 | 0.29 |
| reproduce_Mu2019.R | gBOIN, Scenario 8, Patients DL4 | 8.5 | 7.47 | 0.453 |
| reproduce_Mu2019.R | gBOIN, Scenario 10, Selection (%) DL3 | 3 | 5.79 | 2.034 |
| reproduce_Mu2019.R | gBOIN, Scenario 10, Selection (%) DL4 | 71 | 76.32 | 3.796 |
| reproduce_Mu2019.R | gBOIN, Scenario 10, Selection (%) DL5 | 23 | 15.51 | 3.451 |
| reproduce_Mu2019.R | gBOIN, Scenario 10, Patients DL3 | 5 | 5.638 | 0.47 |
| reproduce_Mu2019.R | gBOIN, Scenario 10, Patients DL4 | 13.3 | 13.91 | 0.602 |
| reproduce_Mu2019.R | gBOIN, Scenario 10, Patients DL5 | 3.8 | 2.777 | 0.359 |
| reproduce_Mu2019.R | gBOIN, Scenario 10, Patients DL6 | 1.3 | 0.989 | 0.236 |
| reproduce_Kojima2024.R | target 0.25, ABOIN, Scenario 1, SS | 14.75 | 14.95 | 0.13 |
| reproduce_Kojima2024.R | target 0.25, BOIN, Scenario 1, MTD_DL4 | 0.49 | 1.01 | 0.49 |
| reproduce_Kojima2024.R | target 0.25, BOIN, Scenario 1, MTD_DL5 | 99.39 | 98.81 | 0.54 |
| reproduce_Takeda2018.R | Scenario 2, allocation, 2 | 5 | 3.87 | 0.98 |
| reproduce_Takeda2018.R | Scenario 4, allocation, 2 | 5 | 3.98 | 0.99 |
| reproduce_Takeda2018.R | Scenario 5, allocation, 1 | 4 | 3.17 | 0.64 |
| reproduce_Takeda2018.R | Scenario 11, recommendation (%), 1 | 0 | 42.05 | 8.15 |
| reproduce_Takeda2018.R | Scenario 11, recommendation (%), 2 | 2 | 0.35 | 1.41 |
| reproduce_Takeda2018.R | Scenario 11, recommendation (%), 3 | 33 | 19.05 | 6.58 |
| reproduce_Takeda2018.R | Scenario 11, recommendation (%), 4 | 42 | 23.45 | 7.06 |
| reproduce_Takeda2018.R | Scenario 11, recommendation (%), 5 | 19 | 12.2 | 5.57 |
| reproduce_Takeda2018.R | Scenario 11, allocation, 1 | 3 | 17.76 | 2.84 |
| reproduce_Takeda2018.R | Scenario 11, allocation, 2 | 5 | 2.83 | 1.06 |
| reproduce_Takeda2018.R | Scenario 11, allocation, 3 | 13 | 7.16 | 1.99 |
| reproduce_Takeda2018.R | Scenario 11, allocation, 4 | 10 | 5.5 | 1.83 |
| reproduce_Takeda2018.R | Scenario 11, allocation, 5 | 4 | 2.27 | 1.38 |
| reproduce_Takeda2018.R | Scenario 12, allocation, 1 | 4 | 6.83 | 1.25 |
| reproduce_Takeda2020.R | TITE-BOIN-ET, Scenario 5, allocation, 1 | 5 | 4.07 | 0.83 |
| reproduce_Takeda2020.R | TITE-BOIN-ET, Scenario 7, recommendation (%), 1 | 1 | 6.3 | 4.26 |
| reproduce_Takeda2020.R | TITE-BOIN-ET, Scenario 8, ET (%), NA | 1 | 0 | 0.5 |
| reproduce_Takeda2020.R | TITE-BOIN-ET, Scenario 9, ET (%), NA | 4 | 1.4 | 2.32 |
| reproduce_Takeda2020.R | TITE-BOIN-ET, Scenario 9, recommendation (%), 1 | 67 | 74.65 | 7.24 |
| reproduce_Takeda2020.R | TITE-BOIN-ET, Scenario 9, duration (weeks), NA | 77 | 79.11 | 2.11 |
| reproduce_Takeda2020.R | TITE-BOIN-ET, Scenario 10, duration (weeks), NA | 76 | 78.94 | 2.16 |
| reproduce_Takeda2020.R | TITE-BOIN-ET, Scenario 11, ET (%), NA | 1 | 0 | 0.5 |
| reproduce_Takeda2020.R | TITE-BOIN-ET, Scenario 11, recommendation (%), 1 | 2 | 41.75 | 8.14 |
| reproduce_Takeda2020.R | TITE-BOIN-ET, Scenario 11, recommendation (%), 3 | 34 | 14.4 | 5.94 |
| reproduce_Takeda2020.R | TITE-BOIN-ET, Scenario 11, recommendation (%), 4 | 37 | 21.25 | 6.84 |
| reproduce_Takeda2020.R | TITE-BOIN-ET, Scenario 11, allocation, 1 | 4 | 14.21 | 2.62 |
| reproduce_Takeda2020.R | TITE-BOIN-ET, Scenario 11, allocation, 2 | 6 | 4.49 | 1.16 |
| reproduce_Takeda2020.R | TITE-BOIN-ET, Scenario 11, allocation, 3 | 12 | 7.86 | 1.83 |
| reproduce_Takeda2020.R | TITE-BOIN-ET, Scenario 11, allocation, 4 | 8 | 6 | 1.72 |
| reproduce_Takeda2020.R | TITE-BOIN-ET, Scenario 12, allocation, 1 | 6 | 8.71 | 1.39 |
| reproduce_Takeda2020.R | TITE-BOIN-ET, Scenario 12, duration (weeks), NA | 14 | 20.72 | 2.61 |
| reproduce_Takeda2020.R | BOIN-ET, Scenario 2, allocation, 2 | 5 | 3.76 | 0.96 |
| reproduce_Takeda2020.R | BOIN-ET, Scenario 4, allocation, 2 | 5 | 3.94 | 0.94 |
| reproduce_Takeda2020.R | BOIN-ET, Scenario 5, allocation, 1 | 4 | 3.21 | 0.67 |
| reproduce_Takeda2020.R | BOIN-ET, Scenario 11, ET (%), NA | 1 | 0.05 | 0.85 |
| reproduce_Takeda2020.R | BOIN-ET, Scenario 11, recommendation (%), 1 | 0 | 38.55 | 8.04 |
| reproduce_Takeda2020.R | BOIN-ET, Scenario 11, recommendation (%), 3 | 34 | 19.8 | 6.67 |
| reproduce_Takeda2020.R | BOIN-ET, Scenario 11, recommendation (%), 4 | 42 | 25.25 | 7.23 |
| reproduce_Takeda2020.R | BOIN-ET, Scenario 11, allocation, 1 | 3 | 16.7 | 2.81 |
| reproduce_Takeda2020.R | BOIN-ET, Scenario 11, allocation, 2 | 5 | 3.03 | 1.07 |
| reproduce_Takeda2020.R | BOIN-ET, Scenario 11, allocation, 3 | 13 | 7.65 | 2 |
| reproduce_Takeda2020.R | BOIN-ET, Scenario 11, allocation, 4 | 10 | 5.8 | 1.84 |
| reproduce_Takeda2020.R | BOIN-ET, Scenario 11, allocation, 5 | 4 | 2.35 | 1.38 |
| reproduce_Takeda2022.R | gBOIN-ET, Scenario 2, Selection (%) DL2 | 9 | 12.27 | 2.807 |
| reproduce_Takeda2022.R | gBOIN-ET, Scenario 3, Patients DL3 | 6.6 | 5.962 | 0.606 |
| reproduce_Takeda2022.R | gBOIN-ET, Scenario 5, Patients DL6 | 2.4 | 2.009 | 0.342 |
| reproduce_Takeda2022.R | gBOIN-ET, Scenario 6, Patients DL1 | 17.7 | 20.599 | 0.787 |
| reproduce_Takeda2022.R | gBOIN-ET, Scenario 7, Patients DL1 | 18.7 | 20.075 | 0.887 |
| reproduce_Takeda2022.R | gBOIN-ET, Scenario 8, Patients DL5 | 1.3 | 0.815 | 0.244 |
| reproduce_Takeda2022.R | gBOIN-ET, Scenario 8, Patients DL6 | 0.2 | 0.06 | 0.091 |
| reproduce_Zhou2022.R | Zhou 2022 Figure 3, Scenario 7, TITE-BOIN12AL, n_pts_obd | 14.3 | 15.44 | 0.94 |
| reproduce_Lu2025.R | Lu 2025 SI Figure S2, Scenario 3, sel_overtox, wT = 0.05 | 6.54 | 4.66 | 1.53 |
| reproduce_Lu2025.R | Lu 2025 SI Figure S2, Scenario 4, sel_overtox, wT = 0.11 | 23.64 | 20.22 | 2.6 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 4, sel_overtox, Comb-BOIN12 (AT=1, AE=2) | 23.64 | 20.22 | 2.6 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 4, sel_overtox, Comb-BOIN12 (AT=2, AE=3) | 23.64 | 20.22 | 2.6 |
| reproduce_Lu2025.R | Lu 2025 SI Figure S2, Scenario 4, sel_overtox, wT = 0.05 | 26.02 | 22.10 | 2.68 |
| reproduce_Lu2025.R | Lu 2025 SI Figure S2, Scenario 4, sel_overtox, wT = 0.17 | 21.48 | 18.92 | 2.51 |
| reproduce_Lu2025.R | Lu 2025 SI Figure S2, Scenario 7, sel_overtox, wT = 0.11 | 23.56 | 20.06 | 2.6 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 7, sel_overtox, Comb-BOIN12 (AT=1, AE=2) | 23.56 | 20.06 | 2.6 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 7, sel_overtox, Comb-BOIN12 (AT=2, AE=3) | 23.56 | 20.06 | 2.6 |
| reproduce_Lu2025.R | Lu 2025 SI Figure S2, Scenario 7, sel_overtox, wT = 0.05 | 26.06 | 20.44 | 2.68 |
| reproduce_Lu2025.R | Lu 2025 SI Figure S2, Scenario 7, sel_overtox, wT = 0.17 | 22.4 | 18.74 | 2.55 |
| reproduce_Lu2025.R | Lu 2025 SI Figure S2, Scenario 13, sel_obdc, wT = 0.17 | 79.6 | 76.90 | 2.47 |
| reproduce_Lu2025.R | Lu 2025 SI Figure S2, Scenario 14, sel_obdc, wT = 0.11 | 42.76 | 37.08 | 3.02 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 14, sel_obdc, Comb-BOIN12 (AT=1, AE=2) | 42.76 | 37.08 | 3.02 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 14, sel_obdc, Comb-BOIN12 (AT=2, AE=3) | 42.76 | 37.08 | 3.02 |
| reproduce_Lu2025.R | Lu 2025 SI Figure S2, Scenario 14, sel_obdc, wT = 0.05 | 43.4 | 37.52 | 3.02 |
| reproduce_Lu2025.R | Lu 2025 SI Figure S2, Scenario 14, sel_obdc, wT = 0.17 | 41.74 | 35.88 | 3.01 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 1, sel_obdc, TITE-Comb-BOIN12 (AT=2, AE=3) | 46.16 | 50.48 | 3.05 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 4, sel_obdc, TITE-Comb-BOIN12 (AT=1, AE=2) | 28.11 | 32.96 | 2.87 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 4, pts_obdc, TITE-Comb-BOIN12 (AT=1, AE=2) | 17.81 | 19.35 | 1.44 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 4, sel_overtox, TITE-Comb-BOIN12 (AT=1, AE=2) | 28.12 | 22.36 | 2.75 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 4, sel_obdc, TITE-Comb-BOIN12 (AT=2, AE=3) | 28.74 | 34.08 | 2.89 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 4, sel_overtox, TITE-Comb-BOIN12 (AT=2, AE=3) | 29.72 | 23.30 | 2.79 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 5, sel_obdc, TITE-Comb-BOIN12 (AT=1, AE=2) | 37.86 | 43.94 | 3.03 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 5, pts_obdc, TITE-Comb-BOIN12 (AT=1, AE=2) | 30.33 | 32.24 | 1.52 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 5, sel_overtox, TITE-Comb-BOIN12 (AT=1, AE=2) | 45.12 | 39.24 | 3.04 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 5, sel_obdc, TITE-Comb-BOIN12 (AT=2, AE=3) | 37.1 | 43.18 | 3.02 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 5, sel_overtox, TITE-Comb-BOIN12 (AT=2, AE=3) | 45.36 | 40.76 | 3.04 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 7, sel_obdc, TITE-Comb-BOIN12 (AT=1, AE=2) | 20.42 | 24.86 | 2.64 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 7, sel_overtox, TITE-Comb-BOIN12 (AT=1, AE=2) | 27.76 | 21.20 | 2.74 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 7, sel_obdc, TITE-Comb-BOIN12 (AT=2, AE=3) | 20.84 | 25.08 | 2.65 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 7, sel_overtox, TITE-Comb-BOIN12 (AT=2, AE=3) | 28.18 | 21.82 | 2.75 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 8, sel_obdc, TITE-Comb-BOIN12 (AT=2, AE=3) | 27.02 | 30.74 | 2.82 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 8, pts_obdc, TITE-Comb-BOIN12 (AT=2, AE=3) | 13.47 | 15.17 | 1.55 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 12, pts_overtox, TITE-Comb-BOIN12 (AT=2, AE=3) | 21.76 | 20.63 | 1.13 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 13, pts_obdc, TITE-Comb-BOIN12 (AT=1, AE=2) | 52.91 | 58.40 | 2.5 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 13, pts_overtox, TITE-Comb-BOIN12 (AT=1, AE=2) | 18.14 | 15.05 | 1.36 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 13, pts_obdc, TITE-Comb-BOIN12 (AT=2, AE=3) | 53.27 | 59.44 | 2.52 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 13, pts_overtox, TITE-Comb-BOIN12 (AT=2, AE=3) | 17.54 | 14.89 | 1.37 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 14, pts_obdc, TITE-Comb-BOIN12 (AT=1, AE=2) | 16.14 | 17.50 | 1.25 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 14, sel_obdc, TITE-Comb-BOIN12 (AT=2, AE=3) | 41.62 | 37.02 | 3.01 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 14, pts_obdc, TITE-Comb-BOIN12 (AT=2, AE=3) | 15.83 | 17.48 | 1.27 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 15, pts_overtox, TITE-Comb-BOIN12 (AT=1, AE=2) | 21.73 | 19.90 | 1.14 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 15, sel_overtox, TITE-Comb-BOIN12 (AT=2, AE=3) | 6.22 | 4.56 | 1.5 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 15, pts_overtox, TITE-Comb-BOIN12 (AT=2, AE=3) | 20.91 | 19.73 | 1.15 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 16, sel_overtox, TITE-Comb-BOIN12 (AT=2, AE=3) | 24.28 | 21.54 | 2.62 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 16, pts_overtox, TITE-Comb-BOIN12 (AT=2, AE=3) | 20.4 | 22.07 | 1.28 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 17, sel_obdc, TITE-Comb-BOIN12 (AT=2, AE=3) | 28.78 | 31.94 | 2.85 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 17, pts_obdc, TITE-Comb-BOIN12 (AT=2, AE=3) | 22.21 | 23.53 | 1.24 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 17, sel_overtox, TITE-Comb-BOIN12 (AT=2, AE=3) | 53.2 | 49.90 | 3.04 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 17, pts_overtox, TITE-Comb-BOIN12 (AT=2, AE=3) | 34.96 | 36.95 | 1.42 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 18, pts_obdc, TITE-Comb-BOIN12 (AT=1, AE=2) | 30.0 | 32.63 | 1.5 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 18, pts_obdc, TITE-Comb-BOIN12 (AT=2, AE=3) | 30.43 | 33.13 | 1.53 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 19, pts_obdc, TITE-Comb-BOIN12 (AT=1, AE=2) | 35.27 | 37.48 | 1.58 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 19, sel_obdc, TITE-Comb-BOIN12 (AT=2, AE=3) | 61.8 | 64.92 | 2.91 |
| reproduce_Lu2025.R | Lu 2025 Figure 2, Scenario 19, pts_obdc, TITE-Comb-BOIN12 (AT=2, AE=3) | 34.88 | 37.54 | 1.58 |

## References

Kakizume, T., Takeda, K., Taguri, M. and Morita, S. (2024). BOIN-ETC: A
Bayesian optimal interval design considering efficacy and toxicity to
identify the optimal dose combinations. *Statistical Methods in Medical
Research*, 33(4), 716-727.

Kojima, M., Wende, W. and Zhao, H. (2024). Tips for Accelerating BOIN
Design. *Therapeutic Innovation & Regulatory Science*, 58, 1129-1137.

Lin, R. and Yin, G. (2017). Bayesian Optimal Interval Design for Dose
Finding in Drug-Combination Trials. *Statistical Methods in Medical
Research*, 26(5), 2155-2167.

Lin, R., Zhou, Y., Yan, F., Li, D. and Yuan, Y. (2020). BOIN12: Bayesian
optimal interval phase I/II trial design for utility-based dose finding
in immunotherapy and targeted therapies. *JCO Precision Oncology*, 4,
1393-1402.

Lu, M., Zhang, J., Yuan, Y. and Lin, R. (2025). Comb-BOIN12: A
utility-based Bayesian optimal interval design for dose optimization in
cancer drug-combination trials. *Statistics in Biopharmaceutical
Research*, 17(2), 266-276.

Mu, R., Yuan, Y., Xu, J., Mandrekar, S. J. and Yin, J. (2019). gBOIN: a
unified model-assisted phase I trial design accounting for toxicity
grades, and binary or continuous end points. *Journal of the Royal
Statistical Society: Series C*, 68, 289-308.

Takeda, K., Taguri, M. and Morita, S. (2018). BOIN-ET: Bayesian optimal
interval design for dose finding based on both efficacy and toxicity
outcomes. *Pharmaceutical Statistics*, 17(4), 383-395.

Takeda, K., Morita, S. and Taguri, M. (2020). TITE-BOIN-ET:
Time-to-event Bayesian optimal interval design to accelerate
dose-finding based on both efficacy and toxicity outcomes.
*Pharmaceutical Statistics*, 19(3), 335-349.

Takeda, K., Morita, S. and Taguri, M. (2022). gBOIN-ET: The generalized
Bayesian optimal interval design for optimal dose-finding accounting for
ordinal graded efficacy and toxicity in early clinical trials.
*Biometrical Journal*, 64(7), 1178-1191.

Takeda, K., Yamaguchi, Y., Taguri, M. and Morita, S. (2023).
TITE-gBOIN-ET: Time-to-event generalized Bayesian optimal interval
design to accelerate dose-finding accounting for ordinal graded efficacy
and toxicity outcomes. *Biometrical Journal*, 65(7), 2200265.

Takeda, K., Zhu, J. and Hirakawa, A. (2025). BF-BOIN-ET: A backfill
Bayesian optimal interval design using efficacy and toxicity outcomes
for dose optimization. *Pharmaceutical Statistics*, 24(2), e2470.

Zhang, L. and Yuan, Y. (2016). A Practical Bayesian Design to Identify
the Maximum Tolerated Dose Contour for Drug Combination Trials.
*Statistics in Medicine*, 35(27), 4924-4936.

Zhao, Y., Yuan, Y., Korn, E. L. and Freidlin, B. (2024). Backfilling
patients in phase I dose-escalation trials using Bayesian optimal
interval design (BOIN). *Clinical Cancer Research*, 30(4), 673-679.

Zhou, Y., Lee, J. J. and Yuan, Y. (2019). A utility-based Bayesian
optimal interval (U-BOIN) phase I/II design to identify the optimal
biological dose for targeted and immune therapies. *Statistics in
Medicine*, 38(28), S5299-S5316.

Zhou, Y., Lee, J. J., Wang, S., Bailey, S. and Yuan, Y. (2021).
Incorporating historical information to improve phase I clinical trials.
*Pharmaceutical Statistics*, 20, 1017-1034.

Zhou, Y., Lin, R., Lee, J. J., Li, D., Wang, L., Li, R. and Yuan, Y.
(2022). TITE-BOIN12: A Bayesian phase I/II trial design to find the
optimal biological dose with late-onset toxicity and efficacy.
*Statistics in Medicine*, 41(11), 1918-1931.
