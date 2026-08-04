# CALM-VOI replication code

This repository reproduces the simulations, generated tables and
figures, and real-data applications for "Value-of-Information
Allocation for Stochastic Surrogate Covariates in Randomized Trials."

## Quick start

Install R 4.4 or later and the required packages:

```r
install.packages(c("BART", "dplyr", "foreign", "ggplot2", "rpart", "tidyr"))
```

Download the two large voter-study files and verify their checksums:

```sh
Rscript download_data.R
```

Run the full replication pipeline from the repository root:

```sh
Rscript run_all.R
```

The submission profile uses 5,000 replications for the main and
score-sensitivity grids; 2,000 for CALM-UQ, external-policy, and
query-noise experiments; 1,000 ACTG, turnout, and robocall query-bank
replays; 200 IST replays; and 200 synthetic implementation runs. The
serial run can take several hours.

For a smoke test, supply one common replication count:

```sh
Rscript run_all.R 100
```

On macOS or Linux, independent simulation cells can run in parallel:

```sh
CALM_VOI_WORKERS=5 Rscript run_all.R
```

## Pipeline

`run_all.R` runs the simulation grid, score sensitivity, CALM-UQ,
external-policy and query-noise experiments, the synthetic component
check, and four real-data applications. It then regenerates the six
figures, LaTeX table fragments, and simulation metadata before running
`validate_outputs.R`.

The validator checks reported summaries against the stored Monte Carlo
draws, confirms paired-comparison invariants, and verifies generated
artifact completeness. The archived submission outputs are included so
the validator can be run immediately after cloning and downloading the
two voter files:

```sh
Rscript validate_outputs.R
```

## Main files

| File | Purpose |
|---|---|
| `calm_voi.R` | Exact allocation, AIPW analysis, data-generating processes, and pilot rules |
| `run_main_grid.R` | Main simulation grid |
| `run_sensitivity.R` | Score-quality sensitivity analysis |
| `run_calmuq.R` | Fixed-gamma CALM-UQ experiment |
| `run_external_pilot.R` | Independently trained frozen policies |
| `run_noise_robustness.R` | Skewed and correlated query-noise checks |
| `worked_example.R` | Synthetic component check |
| `run_real_application.R` | ACTG 175 application |
| `run_turnout_screen.R` | Household turnout application |
| `run_robocall_screen.R` | Treatment-placebo robocall application |
| `run_ist_application.R` | International Stroke Trial application |
| `real_significance.R` | Wald and paired-bootstrap comparisons |
| `make_plots.R` | Generated figures |
| `make_tables.R` | Generated table fragments and metadata |
| `validate_outputs.R` | Reproducibility checks |
| `download_data.R` | Verified downloads for the two voter datasets |

## Software

The archived results were generated and checked with R 4.4.1, BART
2.9.10, ggplot2 4.0.2, tidyr 1.3.2, and dplyr 1.2.0. Base R is
sufficient for the synthetic simulations. The applications additionally
use BART, foreign, and rpart. All stochastic scripts set explicit seeds.
Newer compatible versions may produce small floating-point or MCMC
differences.

See `DATA_SOURCES.md` for data provenance, download locations, and
third-party licenses.
