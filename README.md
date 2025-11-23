# simFastBOIN

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
result$summary
```

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
result1 <- sim_boin(..., boundMTD = FALSE)

# With boundMTD: more conservative selection
result2 <- sim_boin(..., boundMTD = TRUE)
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
result_simple <- sim_boin(..., n_earlystop_rule = "simple")
# Average: 38.4 patients

# With stay: Wait for convergence
result_with_stay <- sim_boin(..., n_earlystop_rule = "with_stay")
# Average: 40.8 patients (more thorough evaluation)
```

### 3. extrasafe: Safety Stopping at Lowest Dose

Additional safety rule to stop the entire trial if the lowest dose is overly toxic:

```r
result <- sim_boin(
  ...,
  extrasafe = TRUE,   # Enable safety stopping
  offset = 0.05       # Cutoff = cutoff_eli - offset
)
```

## Output

The simulation produces a unified, easy-to-read table:

```
BOIN Simulation Summary
| Metric                           |    DL1 |    DL2 |    DL3 |    DL4 |    DL5 |    DL6 | Total/No MTD |
| -------------------------------- |--------|--------|--------|--------|--------|--------|--------------|
| True Toxicity (%)                |    5.0 |   10.0 |   20.0 |   30.0 |   45.0 |   60.0 |            - |
| MTD Selected (%)                 |    0.2 |    5.2 |   29.6 |   50.9 |   13.6 |    0.4 |          0.0 |
| Participants Treated (mean)      |    3.7 |    5.6 |   10.1 |   12.2 |    5.9 |    1.0 |         38.4 |
| Participants w/ DLTs (mean)      |    0.2 |    0.6 |    2.0 |    3.7 |    2.6 |    0.6 |          9.6 |
```

### Interpretation

1. **True Toxicity (%)**: The actual DLT rate at each dose (simulation ground truth)
2. **MTD Selected (%)**: Percentage of trials selecting each dose as MTD
   - Last column shows "No MTD" percentage (trials that failed to identify MTD)
3. **Participants Treated (mean)**: Average patient enrollment at each dose
   - Last column shows total across all doses
4. **Participants w/ DLTs (mean)**: Average DLT counts at each dose
   - Last column shows total DLTs across all doses

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

simFastBOIN uses optimized vectorized operations for high performance:

- **10,000 trials, 5 doses, 10 cohorts**: ~0.1-0.3 seconds
- **10,000 trials, 9 doses, 48 cohorts**: ~1-2 seconds
- **100,000 trials, 9 doses, 48 cohorts**: ~10-20 seconds
- **2-5x faster** than traditional sequential implementations

Performance scales with the number of cohorts (not trials), making it ideal for large-scale simulations.

## Design Comparison

| Feature | Baseline | boundMTD | with_stay | Standard | Conservative |
|---------|----------|----------|-----------|----------|--------------|
| boundMTD | ✗ | ✓ | ✗ | ✓ | ✓ |
| n_earlystop_rule | simple | simple | with_stay | with_stay | with_stay |
| extrasafe | ✗ | ✗ | ✗ | ✗ | ✓ |
| Avg Patients | 38.4 | 38.4 | 40.8 | 40.8 | 40.5 |
| MTD Selection* | 52.9% | 48.6% | 50.3% | 48.8% | 48.2% |

*Selection rate at true MTD (p=0.30)

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
head(result$detailed_results)
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

## Support

For bug reports and feature requests, please open an issue on [GitHub](https://github.com/gosukehommaEX/simFastBOIN/issues).
