# simFastBOIN 1.2.0

## Performance Optimizations

### Isotonic Regression Acceleration

* **Adopted C-based PAVA implementation** 
  - Replaced pure R implementation with `Iso::pava()` (C implementation)
  - Significant speedup for isotonic regression computation
  - Pre-allocated all vectors to avoid repeated memory allocations
  - Vectorized pseudocount and variance calculations
  - Early exit for trials with no valid doses

### DLT Generation Optimization

* **Optimized random number generation in sim_boin()**
  - Replaced `rbinom()` with `runif()` for DLT generation
  - Vectorized threshold comparison: `as.integer(runif(n_active) < p_true_current)`
  - Significantly faster computation due to:
    1. `runif()` is simpler and faster than `rbinom()`
    2. Enables vectorized threshold comparison
    3. Reduces function call overhead
  - Performance improvement particularly notable for large-scale simulations

### MTD Selection Enhancement

* **Optimized select_mtd() function**
  - Returns NA immediately when no valid MTD candidates exist
  - Avoids unnecessary computations for trials without viable doses
  - Improves overall simulation efficiency

## Internal Improvements

* Enhanced code documentation with detailed optimization rationale
* Improved memory efficiency through pre-allocation
* Maintained backward compatibility with existing APIs

---

# simFastBOIN 1.1.0

## New Features

### Conservative MTD Selection: boundMTD

* Added `boundMTD` parameter to `sim_boin()` for more conservative MTD selection
  - When `boundMTD = TRUE`, selected MTD must have isotonic-estimated toxicity rate below de-escalation boundary (lambda_d)
  - Prevents selection of doses too close to overly toxic doses
  - Returns `"no_dose_below_lambda_d"` stopping reason when no dose satisfies constraint
  - Maintains backward compatibility with `boundMTD = FALSE` (default)

### Flexible Early Stopping Rules: n_earlystop_rule

* Added `n_earlystop_rule` parameter to `sim_boin()` with two options:
  - `"simple"` (default): Stop when sample size at current dose reaches `n_earlystop`
    - Backward compatible with previous versions
    - Faster trials with fewer patients
  - `"with_stay"`: Stop when sample size reaches `n_earlystop` AND next decision is "Stay"
    - Follows BOIN standard implementation
    - Ensures algorithm convergence before stopping
    - Results in more thorough dose evaluation

### Enhanced Stopping Reason Tracking

* Improved tracking of trial termination reasons:
  - `"n_earlystop_reached"`: Simple early stopping rule triggered
  - `"n_earlystop_with_stay"`: Convergence-based early stopping rule triggered
  - `"no_dose_below_lambda_d"`: boundMTD constraint violation
  - Existing reasons maintained for backward compatibility

## Bug Fixes

### Critical Escalation Logic Fix

* **Fixed index out-of-bounds error in escalation logic** (#issue)
  - Problem: When trials at maximum dose received "Escalate" decision, range check was performed after index access
  - Impact: Caused simulation crashes with error "subscript out of bounds"
  - Solution: Implemented two-stage check:
    1. Filter trials not at maximum dose
    2. Check if next dose is not eliminated (only for eligible trials)
  - This fix ensures stable operation for all dose configurations

### MTD Selection Logic Correction

* **Fixed MTD selection for early-stopped trials** (#issue)
  - Problem: Trials stopped due to `n_earlystop` rules were incorrectly excluded from MTD selection
  - Impact: All trials showed "No MTD Selected = 100%" in some scenarios
  - Solution: Distinguish between terminal stopping reasons (prevent MTD selection) and completion reasons (allow MTD selection)
  - Early stopping due to sample size limits now correctly proceeds to MTD selection

## Internal Improvements

### Enhanced .select_mtd_batch()

* Refactored internal MTD selection function for clarity and correctness:
  - Separated stopping reasons into two categories:
    - Terminal reasons: `"lowest_dose_too_toxic"`, `"lowest_dose_eliminated"`
    - Completion reasons: `"n_earlystop_reached"`, `"n_earlystop_with_stay"`
  - Trials with completion reasons now undergo normal MTD selection
  - Added support for `boundMTD` constraint checking
  - Improved code readability with explicit reason categorization

### Vectorized Escalation Processing

* Optimized escalation logic in main simulation loop:
  - Safer index operations with explicit bounds checking
  - Maintained vectorization for performance
  - Added clear comments for maintainability

## Documentation

### Updated Documentation

* Enhanced `sim_boin()` documentation:
  - Added detailed descriptions of `boundMTD` and `n_earlystop_rule` parameters
  - Included practical examples comparing different configurations
  - Added references to BOIN standard implementation
  - Documented all possible stopping reasons

### New README Sections

* Added comprehensive feature descriptions:
  - Detailed explanation of `boundMTD` with examples
  - Comparison of `"simple"` vs `"with_stay"` stopping rules
  - Design comparison table showing trade-offs
  - Extended complete workflow examples

### New Files

* Created NEWS.md for version history tracking

## Testing

* Validated new features with comprehensive test scenarios:
  - Baseline configuration (backward compatibility)
  - boundMTD only
  - n_earlystop_rule = "with_stay" only
  - Combined BOIN standard (boundMTD + with_stay)
  - Maximum conservatism (all safety features)
  - Edge cases (boundary violations, all-toxic scenarios)

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
