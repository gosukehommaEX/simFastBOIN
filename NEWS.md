# simFastBOIN 1.2.1

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
  - `print.boin_summary()`: Summary table printing and formatting

* **Enhanced function comments**
  - Added section headers for major processing blocks
  - Detailed explanations of algorithm steps and decision logic
  - Clear descriptions of parameter usage and constraints
  - Clarified relationships between input parameters and outputs

### Print Method Enhancements: print.boin_summary()

* **Added `percent` parameter**
  - `percent = FALSE` (default): Display Avg Pts and Avg DLTs as absolute numbers
  - `percent = TRUE`: Display as percentages of trial totals
  - Enables flexible presentation of operating characteristics
  - Row labels update automatically based on selection

* **Added `kable_format` parameter**
  - `kable_format = "pipe"` (default): Markdown pipe table format
  - `kable_format = "simple"`: Minimal text table format
  - `kable_format = "latex"`: LaTeX table format
  - Enables output for RMarkdown documents and reports
  - Provides publication-ready table formatting

* **Enhanced print output**
  - `scenario_name` parameter for identifying multiple result sets
  - Improved table formatting with flexible options
  - Better integration with RMarkdown workflows

## Documentation Files

* Updated DESCRIPTION file with version 1.2.1
* All roxygen2 documentation includes @param, @return, @details, @examples, @importFrom, and @export tags
* Consistent formatting across all function headers

## Breaking Changes

None. All changes are backward compatible. New parameters have sensible defaults.

## Migration Guide

For users upgrading from simFastBOIN 1.2.0:

```r
# Enhanced printing with percentages
result <- sim_boin(n_trials = 1000, target = 0.30, p_true = p_true, ...)
print(result$summary, percent = TRUE)

# Kable format for RMarkdown
print(result$summary, kable = TRUE, kable_format = "pipe")

# Combine options
print(result$summary, scenario_name = "My Analysis", 
      percent = TRUE, kable = TRUE, kable_format = "latex")
```

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

### Isotonic Regression Acceleration

* **C-based PAVA implementation** 
  - Uses `Iso::pava()` (C implementation) for fast isotonic regression
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

For users upgrading from simFastBOIN 1.1.x:

If you were using defaults:
```r
# Old code (simFastBOIN 1.1.x)
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
