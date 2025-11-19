# simFastBOIN

Fast and efficient simulation tools for Bayesian Optimal Interval (BOIN) designs in Phase I clinical trials.

## Overview

simFastBOIN provides comprehensive functions for simulating Phase I dose-finding trials using the BOIN design methodology. The package enables researchers to evaluate operating characteristics and design performance across different dose-toxicity scenarios, essential for protocol development and regulatory submissions.

## Installation

You can install the development version from GitHub:

```r
# Install if not already installed
# install.packages("devtools")
devtools::install_github("gosukehommaEX/simFastBOIN")
```

## Quick Start

### 1. Load the package

```r
library(simFastBOIN)
```

### 2. Define design parameters

```r
# Target DLT (Dose-Limiting Toxicity) rate
target <- 0.30

# True toxicity probabilities for each dose level
p_true <- seq(0.05, 0.45, by = 0.05)

# Number of doses
n_doses <- length(p_true)  # 9 doses
```

### 3. Generate BOIN boundaries and decision table

```r
# Get BOIN boundaries for the target DLT rate
boin_bound <- get_boin_boundary(target)

# Generate decision table
decision_table <- get_boin_decision(
  target = target,
  lambda_e = boin_bound$lambda_e,
  lambda_d = boin_bound$lambda_d,
  n_earlystop = 18,
  cutoff_eli = 0.95
)

# Generate stopping boundaries
stopping_boundaries <- get_boin_stopping_boundaries(
  target = target,
  n_earlystop = 18,
  p0 = 0.90
)
```

### 4. Run simulation

```r
result <- sim_boin(
  n_trials = 10000,
  target = target,
  p_true = p_true,
  n_doses = 9,
  n_cohort = 48,
  cohort_size = 3,
  decision_table = decision_table,
  stopping_boundaries = stopping_boundaries,
  cutoff_eli = 0.95, 
  min_mtd_sample = 6,
  seed = 123
)
```

### 5. Display results

```r
# Standard format
print_simulation_summary_boin(result$summary, scenario_name = "Base Case")

# knitr::kable format (recommended for reports)
print_simulation_summary_boin(result$summary, 
                              scenario_name = "Base Case",
                              kable_output = TRUE)
```

## Main Functions

### Core Simulation Functions

- **`sim_boin()`**: Run multiple BOIN trial simulations with specified design parameters
- **`sim_boin_one_trial()`**: Simulate a single BOIN trial
- **`summarize_simulation_boin()`**: Aggregate results from multiple simulations

### Design Setup Functions

- **`get_boin_boundary()`**: Calculate BOIN design boundaries (lambda_e, lambda_d) for a target DLT rate
- **`get_boin_decision()`**: Generate the BOIN decision table for dose escalation/de-escalation
- **`get_boin_stopping_boundaries()`**: Generate trial stopping boundaries

### Output Functions

- **`print_simulation_summary_boin()`**: Display formatted simulation results with operating characteristics
  - `scenario_name`: Optional scenario label
  - `kable_output`: If `TRUE`, displays results in knitr::kable format (suitable for RMarkdown and reports)

## Output Interpretation

The simulation summary provides four key operating characteristics:

1. **MTD Selected (%)**: Percentage of simulations selecting each dose as the Maximum Tolerated Dose
2. **Number of Participants Treated (mean)**: Average patient enrollment at each dose level
3. **Number of Participants w/ DLTs (mean)**: Average number of dose-limiting toxicity events at each dose
4. **% No MTD Selected (N/S)**: Percentage of trials that failed to identify an MTD

These metrics are essential for evaluating trial design performance and for inclusion in clinical trial protocols.

## Example: Complete Workflow

```r
library(simFastBOIN)

# Design parameters
target <- 0.30
p_true <- seq(0.05, 0.45, by = 0.05)

# Get boundaries and decision table
boin_bound <- get_boin_boundary(target)
decision_table <- get_boin_decision(
  target, 
  boin_bound$lambda_e, 
  boin_bound$lambda_d, 
  18, 
  0.95
)
stopping_boundaries <- get_boin_stopping_boundaries(target, 18, 0.90)

# Run simulation
result <- sim_boin(
  n_trials = 10000,
  target = target,
  p_true = p_true,
  n_doses = 9,
  n_cohort = 48,
  cohort_size = 3,
  decision_table = decision_table,
  stopping_boundaries = stopping_boundaries,
  cutoff_eli = 0.95, 
  min_mtd_sample = 6,
  seed = 123
)

# Display results in knitr::kable format
print_simulation_summary_boin(result$summary, kable_output = TRUE)
```

### Sample Output

```
MTD Selected (%)
|  DL1 |  DL2 |  DL3 |  DL4 |  DL5 |  DL6 |  DL7 |  DL8 |  DL9 |
|------|------|------|------|------|------|------|------|------|
|  0.3 |  2.1 |  8.8 | 18.1 | 25.9 | 23.8 | 13.9 |  5.3 |  1.3 |

Number of Participants Treated (mean)
|  DL1 |  DL2 |  DL3 |  DL4 |  DL5 |  DL6 |  DL7 |  DL8 |  DL9 | Total |
|------|------|------|------|------|------|------|------|------|-------|
|  3.7 |  4.9 |  6.5 |    8 |  8.5 |  7.1 |  4.5 |  2.1 |  0.7 |  45.9 |

Number of Participants w/ DLTs (mean)
|  DL1 |  DL2 |  DL3 |  DL4 |  DL5 |  DL6 |  DL7 |  DL8 |  DL9 | Total |
|------|------|------|------|------|------|------|------|------|-------|
|  0.2 |  0.5 |    1 |  1.6 |  2.1 |  2.2 |  1.6 |  0.8 |  0.3 |  10.2 |

% No MTD Selected (N/S): 0.6 %
```

## References

Liu, S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical Trials. *Journal of the Royal Statistical Society: Series C*, 64, 507–523.

## License

MIT © Gosuke Homma
