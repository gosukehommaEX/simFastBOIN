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
