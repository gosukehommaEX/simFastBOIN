# simFastBOIN 1.5.0

Three design options requested for a quality control review. Every new option is
off by default, so a simulation run with version 1.4.0 arguments produces exactly
the same trials, selections and summary figures. Two names in the returned
objects did change, and are listed under their respective headings below:
`overdose60` and `overdose80` moved into the new `overdose` component, and one
reason code was renamed.

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

  The cutoff defaults to `target`, which reproduces the earlier behaviour.

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

## Documentation

* `start_dose`, added in version 1.4.0, is now documented as being ignored when
  `titration = TRUE`, because the titration phase always begins at the lowest
  dose. This matches the reference implementation.

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
  draws the table as a grid of coloured cells, each carrying its decision letter
  so that the figure does not rely on colour alone. It needs `ggplot2`, which is
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
