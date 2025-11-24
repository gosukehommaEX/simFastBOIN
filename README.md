# simFastBOIN <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/simFastBOIN)](https://CRAN.R-project.org/package=simFastBOIN)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

Fast and efficient simulation tools for Bayesian Optimal Interval (BOIN) designs in Phase I clinical trials.

## Overview

simFastBOIN provides comprehensive functions for simulating Phase I dose-finding trials using the BOIN design methodology. The package enables researchers to evaluate operating characteristics and design performance across different dose-toxicity scenarios, essential for protocol development and regulatory submissions.

**Key Features:**
- **High Performance**: Optimized vectorized implementation (2-5x faster than traditional approaches)
- **User-Friendly**: Simplified API with automatic decision table generation
- **Flexible**: Multiple safety stopping rules and customizable design parameters
- **Conservative MTD Selection**: Optional boundMTD constraint for enhanced safety
- **BOIN Standard Compliance**: Supports both backward-compatible and BOIN-standard implementations
- **Professional Output**: Publication-ready tables for reports and manuscripts

## Installation

You can install the development version from GitHub:

```r
# Install if not already installed
# install.packages("devtools")
devtools::install_github("gosukehommaEX/simFastBOIN")
```

## Quick Start

### Basic Usage

```r
library(simFastBOIN)

# Define design parameters
target <- 0.30  # Target DLT rate (30%)
p_true <- seq(0.05, 0.45, by = 0.05)  # True toxicity probabilities

# Run simulation (all decision tables generated automatically!)
result <- sim_boin(
  n_trials = 10000,
  target = target,
  p_true = p_true,
  n_doses = 9,
  n_cohort = 48,
  cohort_size = 3,
  seed = 123
)

# Display results
print(result$summary)
```

### Output

The code above produces a unified, easy-to-read summary table:

```
BOIN Simulation Summary
| Metric                           |    DL1 |    DL2 |    DL3 |    DL4 |    DL5 |    DL6 |    DL7 |    DL8 |    DL9 | Total/No MTD |
| -------------------------------- |--------|--------|--------|--------|--------|--------|--------|--------|--------|--------------|
| True Toxicity (%)                |    5.0 |   10.0 |   15.0 |   20.0 |   25.0 |   30.0 |   35.0 |   40.0 |   45.0 |            - |
| MTD Selected (%)                 |    0.3 |    2.5 |    8.1 |   18.9 |   25.8 |   23.5 |   14.2 |    5.4 |    1.3 |          0.0 |
| Participants Treated (mean)      |    3.7 |    4.9 |    6.5 |    8.1 |    8.6 |    7.0 |    4.5 |    2.1 |    0.7 |         46.1 |
| Participants w/ DLTs (mean)      |    0.2 |    0.5 |    1.0 |    1.6 |    2.2 |    2.1 |    1.6 |    0.8 |    0.3 |         10.2 |
```

**Interpretation:**

1. **True Toxicity (%)**: The actual DLT rate at each dose (simulation ground truth)
2. **MTD Selected (%)**: Percentage of trials selecting each dose as MTD
   - Last column shows "No MTD" percentage (trials that failed to identify MTD)
3. **Participants Treated (mean)**: Average patient enrollment at each dose
   - Last column shows total across all doses
4. **Participants w/ DLTs (mean)**: Average DLT counts at each dose
   - Last column shows total DLTs across all doses

### Advanced Usage: BOIN Standard Implementation

```r
# BOIN standard with conservative MTD selection and convergence-based stopping
result_standard <- sim_boin(
  n_trials = 10000,
  target = 0.30,
  p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
  n_doses = 5,
  n_cohort = 10,
  cohort_size = 3,
  boundMTD = TRUE,              # Conservative MTD selection
  n_earlystop_rule = "with_stay",  # Stop only when converged
  seed = 123
)

print(result_standard$summary)
```

### Maximum Conservatism

```r
# All safety features enabled
result_conservative <- sim_boin(
  n_trials = 10000,
  target = 0.30,
  p_true = seq(0.05, 0.45, by = 0.05),
  n_doses = 9,
  n_cohort = 48,
  cohort_size = 3,
  extrasafe = TRUE,             # Safety stopping at lowest dose
  boundMTD = TRUE,              # Conservative MTD selection
  n_earlystop_rule = "with_stay",  # Convergence-based stopping
  seed = 123
)

print(result_conservative$summary)
```

## Key Features

### 1. boundMTD: Conservative MTD Selection

The `boundMTD` option provides an additional safety constraint during MTD selection:

- Selected MTD must have isotonic-estimated toxicity rate **below** the de-escalation boundary (lambda_d)
- Prevents selection of doses that are too close to overly toxic doses
- Returns `"no_dose_below_lambda_d"` if no dose satisfies the constraint

**Example:**
```r
# Without boundMTD: may select doses near toxicity boundary
result1 <- sim_boin(
  n_trials = 1000,
  target = 0.30,
  p_true = c(0.10, 0.20, 0.28, 0.36, 0.50, 0.65),
  n_doses = 6,
  n_cohort = 20,
  cohort_size = 3,
  boundMTD = FALSE,
  seed = 123
)

# With boundMTD: more conservative selection
result2 <- sim_boin(
  n_trials = 1000,
  target = 0.30,
  p_true = c(0.10, 0.20, 0.28, 0.36, 0.50, 0.65),
  n_doses = 6,
  n_cohort = 20,
  cohort_size = 3,
  boundMTD = TRUE,
  seed = 123
)

# Compare MTD selection rates
result1$summary$mtd_selection_percent
result2$summary$mtd_selection_percent
```

### 2. n_earlystop_rule: Stopping Criteria

Two stopping rules are available:

#### "simple" (Default - Backward Compatible)
- Stop when sample size at current dose reaches `n_earlystop`
- Faster trials, fewer patients

#### "with_stay" (BOIN Standard)
- Stop when sample size reaches `n_earlystop` **AND** next decision is "Stay"
- Ensures algorithm convergence before stopping
- More patients, better MTD identification

**Comparison:**
```r
# Simple rule: Stop immediately at n_earlystop
result_simple <- sim_boin(
  n_trials = 1000,
  target = 0.30,
  p_true = c(0.05, 0.10, 0.20, 0.30, 0.45, 0.60),
  n_doses = 6,
  n_cohort = 20,
  cohort_size = 3,
  n_earlystop_rule = "simple",
  seed = 123
)

# With stay: Wait for convergence
result_with_stay <- sim_boin(
  n_trials = 1000,
  target = 0.30,
  p_true = c(0.05, 0.10, 0.20, 0.30, 0.45, 0.60),
  n_doses = 6,
  n_cohort = 20,
  cohort_size = 3,
  n_earlystop_rule = "with_stay",
  seed = 123
)

# Compare average total patients
cat("Simple rule:", result_simple$summary$avg_total_n_pts, "patients\n")
cat("With stay:", result_with_stay$summary$avg_total_n_pts, "patients\n")
```

### 3. extrasafe: Safety Stopping at Lowest Dose

Additional safety rule to stop the entire trial if the lowest dose is overly toxic:

```r
result <- sim_boin(
  n_trials = 10000,
  target = 0.30,
  p_true = seq(0.05, 0.45, by = 0.05),
  n_doses = 9,
  n_cohort = 48,
  cohort_size = 3,
  extrasafe = TRUE,   # Enable safety stopping
  offset = 0.05,      # Cutoff = cutoff_eli - offset
  seed = 123
)
```

## Main Functions

### Simulation

- **`sim_boin()`**: Run BOIN trial simulations with automatic decision table generation
  - Automatically generates BOIN boundaries and decision tables
  - Optional safety stopping rules (`extrasafe`, `boundMTD`)
  - Flexible stopping criteria (`n_earlystop_rule`)
  - Returns detailed trial-level results and summary statistics

### Design Setup

- **`get_boin_boundary()`**: Calculate BOIN interval boundaries (lambda_e, lambda_d)
- **`get_boin_decision()`**: Generate decision table for dose escalation/de-escalation
- **`get_boin_stopping_boundaries()`**: Generate safety stopping boundaries

### MTD Selection

- **`isotonic_regression()`**: Apply isotonic regression to estimate dose-toxicity curve
- Internal functions for batch processing in simulations

### Utilities

- **`summarize_simulation_boin()`**: Aggregate simulation results into operating characteristics
- **`print.boin_simulation_summary()`**: Print formatted summary tables

## Performance

simFastBOIN uses optimized vectorized operations and efficient algorithms for exceptional performance:

**Benchmark** (9 doses, 48 cohorts, all safety features enabled):

| Number of Trials | Computation Time |
|-----------------|------------------|
| 1,000           | ~0.03 seconds    |
| 10,000          | ~0.25 seconds    |
| 100,000         | ~2.7  seconds     |

Performance scales primarily with the number of cohorts (not trials), making it ideal for large-scale simulations. The package is **2-5x faster** than traditional sequential implementations.

### Key Optimizations

simFastBOIN achieves superior performance through several optimization techniques:

1. **C-based Isotonic Regression**
   - Uses `Iso::pava()` (C implementation) for isotonic regression
   - Significantly faster than pure R implementations
   - Pre-allocated vectors and vectorized computations
   - Early exit for trials with no valid doses

2. **Optimized Random Number Generation**
   - DLT generation uses `runif()` instead of `rbinom()`
   - Vectorized threshold comparison: `as.integer(runif(n) < p_true)`
   - Reduces function call overhead and improves memory access patterns
   - Particularly beneficial for large-scale simulations

3. **Batch Processing Architecture**
   - Processes all trials simultaneously at each cohort
   - Matrix-based operations instead of nested loops
   - Reduces iterations from `n_trials × n_cohorts` to just `n_cohorts`

4. **Efficient MTD Selection**
   - Batch processing of MTD selection across all trials
   - Early termination when no valid candidates exist
   - Vectorized distance calculations and candidate identification

These optimizations maintain full compatibility with BOIN methodology while delivering substantial performance gains.

## Design Comparison

Compare different safety configurations using the same scenario:

```r
library(simFastBOIN)

# Common parameters
target <- 0.30
p_true <- c(0.05, 0.10, 0.20, 0.30, 0.45, 0.60)
n_doses <- 6
n_trials <- 1000
seed_val <- 123

# Test 1: Baseline (backward compatible)
result_baseline <- sim_boin(
  n_trials = n_trials,
  target = target,
  p_true = p_true,
  n_doses = n_doses,
  n_cohort = 20,
  cohort_size = 3,
  seed = seed_val
)

# Test 2: boundMTD only
result_boundMTD <- sim_boin(
  n_trials = n_trials,
  target = target,
  p_true = p_true,
  n_doses = n_doses,
  n_cohort = 20,
  cohort_size = 3,
  boundMTD = TRUE,
  seed = seed_val
)

# Test 3: with_stay only
result_with_stay <- sim_boin(
  n_trials = n_trials,
  target = target,
  p_true = p_true,
  n_doses = n_doses,
  n_cohort = 20,
  cohort_size = 3,
  n_earlystop_rule = "with_stay",
  seed = seed_val
)

# Test 4: BOIN standard (boundMTD + with_stay)
result_standard <- sim_boin(
  n_trials = n_trials,
  target = target,
  p_true = p_true,
  n_doses = n_doses,
  n_cohort = 20,
  cohort_size = 3,
  boundMTD = TRUE,
  n_earlystop_rule = "with_stay",
  seed = seed_val
)

# Test 5: Maximum conservatism (all features)
result_conservative <- sim_boin(
  n_trials = n_trials,
  target = target,
  p_true = p_true,
  n_doses = n_doses,
  n_cohort = 20,
  cohort_size = 3,
  extrasafe = TRUE,
  boundMTD = TRUE,
  n_earlystop_rule = "with_stay",
  seed = seed_val
)

# Create comparison table
comparison <- data.frame(
  Setting = c("Baseline", "boundMTD", "with_stay", "Standard", "Conservative"),
  boundMTD = c("✗", "✓", "✗", "✓", "✓"),
  n_earlystop_rule = c("simple", "simple", "with_stay", "with_stay", "with_stay"),
  extrasafe = c("✗", "✗", "✗", "✗", "✓"),
  Avg_Patients = c(
    result_baseline$summary$avg_total_n_pts,
    result_boundMTD$summary$avg_total_n_pts,
    result_with_stay$summary$avg_total_n_pts,
    result_standard$summary$avg_total_n_pts,
    result_conservative$summary$avg_total_n_pts
  ),
  MTD_Selection = c(
    result_baseline$summary$mtd_selection_percent[4],  # DL4 = true MTD
    result_boundMTD$summary$mtd_selection_percent[4],
    result_with_stay$summary$mtd_selection_percent[4],
    result_standard$summary$mtd_selection_percent[4],
    result_conservative$summary$mtd_selection_percent[4]
  )
)

print(comparison)
```

**Expected Output:**

| Setting | boundMTD | n_earlystop_rule | extrasafe | Avg_Patients | MTD_Selection* |
|---------|----------|------------------|-----------|--------------|----------------|
| Baseline | ✗ | simple | ✗ | 38.6 | 53.0% |
| boundMTD | ✓ | simple | ✗ | 38.6 | 46.9% |
| with_stay | ✗ | with_stay | ✗ | 40.7 | 50.0% |
| Standard | ✓ | with_stay | ✗ | 40.7 | 49.2% |
| Conservative | ✓ | with_stay | ✓ | 40.9 | 50.3% |

*Selection rate at true MTD (DL4, p=0.30)

**Key Observations:**
- `boundMTD` reduces selection of doses near the toxicity boundary (more conservative)
- `with_stay` increases average patient enrollment (ensures convergence)
- `extrasafe` has minimal impact when lower doses are safe

## Complete Workflow Example

```r
library(simFastBOIN)

# Step 1: Define scenario
target <- 0.30
p_true <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45)

# Step 2: Run simulation with BOIN standard settings
result <- sim_boin(
  n_trials = 10000,
  target = target,
  p_true = p_true,
  n_doses = 9,
  n_cohort = 48,
  cohort_size = 3,
  boundMTD = TRUE,
  n_earlystop_rule = "with_stay",
  extrasafe = TRUE,
  seed = 123
)

# Step 3: Display results
print(result$summary, scenario_name = "Linear Toxicity Scenario")

# Step 4: Export for publication (if using RMarkdown)
print(result$summary, kable_output = TRUE)

# Step 5: Access detailed trial-level results
# Each trial contains: n_pts, n_tox, mtd, iso_est, reason, cohorts_completed
head(result$detailed_results)

# Step 6: Check stopping reasons
table(sapply(result$detailed_results, function(x) x$reason))
```

## Stopping Reasons

The package tracks why each trial terminated:

- `"trial_completed"`: Normal completion
- `"n_earlystop_reached"`: Reached sample size limit (simple rule)
- `"n_earlystop_with_stay"`: Converged at sample size limit (with_stay rule)
- `"lowest_dose_too_toxic"`: Safety stopping at lowest dose (extrasafe)
- `"lowest_dose_eliminated"`: Lowest dose eliminated during trial
- `"no_dose_below_lambda_d"`: No dose satisfies boundMTD constraint

## References

Liu, S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical Trials. *Journal of the Royal Statistical Society: Series C*, 64, 507–523.

## Version History

See [NEWS.md](NEWS.md) for detailed version history.

## License

MIT © Gosuke Homma
