# simFastBOIN

Fast and efficient simulation tools for Bayesian Optimal Interval (BOIN) designs in Phase I clinical trials.

## Overview

simFastBOIN provides comprehensive functions for simulating Phase I dose-finding trials using the BOIN design methodology. The package enables researchers to evaluate operating characteristics and design performance across different dose-toxicity scenarios, essential for protocol development and regulatory submissions.

**Key Features:**
- **High Performance**: Optimized vectorized implementation (2-5x faster than traditional approaches)
- **User-Friendly**: Simplified API with automatic decision table generation
- **Flexible**: Optional safety stopping rules and customizable design parameters
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

### Output

The simulation produces a unified, easy-to-read table:

```
BOIN Simulation Summary
| Metric                           |    DL1 |    DL2 |    DL3 |    DL4 |    DL5 |    DL6 |    DL7 |    DL8 |    DL9 | Total/No MTD |
| -------------------------------- |--------|--------|--------|--------|--------|--------|--------|--------|--------|--------------|
| True Toxicity (%)                |    5.0 |   10.0 |   15.0 |   20.0 |   25.0 |   30.0 |   35.0 |   40.0 |   45.0 |            - |
| MTD Selected (%)                 |    0.2 |    2.2 |    8.3 |   19.1 |   25.2 |   24.1 |   13.9 |    4.9 |    1.4 |          0.7 |
| Participants Treated (mean)      |    3.6 |    4.8 |    6.5 |    8.1 |    8.3 |    7.1 |    4.5 |    2.0 |    0.7 |         45.7 |
| Participants w/ DLTs (mean)      |    0.2 |    0.5 |    1.0 |    1.6 |    2.1 |    2.1 |    1.6 |    0.8 |    0.3 |         10.2 |
```

## Main Functions

### Simulation

- **`sim_boin()`**: Run BOIN trial simulations with automatic decision table generation
  - Automatically generates BOIN boundaries and decision tables
  - Optional safety stopping rule (`extrasafe = TRUE`)
  - Returns detailed trial-level results and summary statistics

### Design Setup (Advanced Users)

- **`get_boin_boundary()`**: Calculate BOIN design boundaries (λ_e, λ_d)
- **`get_boin_decision()`**: Generate decision table for dose escalation/de-escalation
- **`get_boin_stopping_boundaries()`**: Generate safety stopping boundaries

### Analysis

- **`summarize_simulation_boin()`**: Aggregate simulation results
- **`print.boin_summary()`**: Display formatted summary (S3 method)

## Key Parameters

### Required Parameters

- **`n_trials`**: Number of trials to simulate (default: 10000)
- **`target`**: Target toxicity probability (e.g., 0.30 for 30%)
- **`p_true`**: True toxicity probabilities for each dose
- **`n_doses`**: Number of doses to evaluate
- **`n_cohort`**: Maximum number of cohorts per trial
- **`cohort_size`**: Patients per cohort (scalar or vector)

### Optional Parameters

- **`n_earlystop`**: Sample size triggering early stopping (default: 18)
- **`cutoff_eli`**: Cutoff for dose elimination (default: 0.95)
- **`extrasafe`**: Enable additional safety stopping at lowest dose (default: FALSE)
- **`offset`**: Safety stopping offset when `extrasafe = TRUE` (default: 0.05)
- **`min_mtd_sample`**: Minimum sample size for MTD consideration (default: 6)
- **`seed`**: Random seed for reproducibility (default: 123)

## Advanced Examples

### With Safety Stopping Rule

```r
# Enable extra safety stopping at the lowest dose
result_safe <- sim_boin(
  n_trials = 10000,
  target = 0.30,
  p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
  n_doses = 5,
  n_cohort = 10,
  cohort_size = 3,
  extrasafe = TRUE,  # Enable safety stopping
  offset = 0.05,     # Safety cutoff = 0.95 - 0.05 = 0.90
  seed = 123
)

result_safe$summary
```

### Custom Design Parameters

```r
# Customize elimination threshold and early stopping
result_custom <- sim_boin(
  n_trials = 10000,
  target = 0.25,
  p_true = seq(0.05, 0.45, by = 0.05),
  n_doses = 9,
  n_cohort = 48,
  cohort_size = 3,
  n_earlystop = 21,      # Higher early stopping threshold
  cutoff_eli = 0.90,     # Less conservative elimination
  extrasafe = TRUE,
  offset = 0.10,         # More conservative safety stopping
  min_mtd_sample = 9,    # Require more patients for MTD
  seed = 123
)
```

### Output Formats

```r
# Standard console output
print(result$summary)

# For RMarkdown/knitr documents
print(result$summary, kable_output = TRUE)

# With scenario name
print(result$summary, scenario_name = "Scenario 1: Linear Increase")
```

### Accessing Detailed Results

```r
# Get results from individual trials
first_trial <- result$detailed_results[[1]]

# Check patient allocation in first trial
first_trial$n_pts

# Check selected MTD in first trial
first_trial$mtd

# Check termination reason
first_trial$reason
```

## Output Interpretation

The unified summary table provides four key metrics:

1. **True Toxicity (%)**: The actual DLT rate at each dose (simulation ground truth)
2. **MTD Selected (%)**: Percentage of trials selecting each dose as MTD
   - Last column shows "No MTD" percentage (trials that failed to identify MTD)
3. **Participants Treated (mean)**: Average patient enrollment at each dose
   - Last column shows total across all doses
4. **Participants w/ DLTs (mean)**: Average DLT counts at each dose
   - Last column shows total DLTs across all doses

These operating characteristics are essential for:
- Protocol development and optimization
- Regulatory submissions
- Sample size justification
- Risk-benefit assessment

## Performance

simFastBOIN uses optimized vectorized operations for high performance:

- **10,000 trials, 5 doses, 10 cohorts**: ~0.1-0.3 seconds
- **100,000 trials, 9 doses, 48 cohorts**: ~2-5 seconds
- **2-5x faster** than traditional sequential implementations

Performance scales with the number of cohorts (not trials), making it ideal for large-scale simulations.

## Complete Workflow Example

```r
library(simFastBOIN)

# Step 1: Define scenario
target <- 0.30
p_true <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45)

# Step 2: Run simulation
result <- sim_boin(
  n_trials = 10000,
  target = target,
  p_true = p_true,
  n_doses = 9,
  n_cohort = 48,
  cohort_size = 3,
  extrasafe = TRUE,
  seed = 123
)

# Step 3: Display results
print(result$summary, scenario_name = "Linear Toxicity Scenario")

# Step 4: Export for publication (if using RMarkdown)
print(result$summary, kable_output = TRUE)
```

## References

Liu, S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical Trials. *Journal of the Royal Statistical Society: Series C*, 64, 507–523.

## License

MIT © Gosuke Homma

## Support

For bug reports and feature requests, please open an issue on [GitHub](https://github.com/gosukehommaEX/simFastBOIN/issues).
