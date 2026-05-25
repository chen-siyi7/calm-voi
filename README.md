# CALM-VOI

R implementation of

**Value-of-Information Allocation for Stochastic Surrogate Covariates in Randomized Trials**

This repository contains the R simulation pipeline that produced the
empirical results in the manuscript. The code implements:

- **Six estimators**: AIPW (no surrogate), CALM (single surrogate query per
  patient), CALM-Unif (uniform budget allocation), CALM-Neyman
  (uncertainty-weighted Neyman allocation), CALM-VOI (plug-in water-filling on
  $q_i = \hat B_i \hat\sigma_i$), and CALM-VOI-oracle (water-filling on the true $q_i$).
- **Six data-generating processes**: 1a, 1b, 1c (constant leverage, varying noise),
  1d (heterogeneous noise), 2 (heterogeneous leverage), 3 (arm-cancellation).
- **Cross-fitted AIPW** with stratified outcome regression, water-filling
  allocation, integer rounding, and sandwich variance estimation.

## Quick start

```sh
git clone https://github.com/[your-username]/calm-voi.git
cd calm-voi
Rscript scripts/run_all.R          # full B = 1000 grid, ~30 min
Rscript scripts/run_all.R 100      # quick B = 100, ~5 min
Rscript scripts/run_all.R 500      # B = 500, ~15 min
```

Output CSVs and RDS files are written to `outputs/`.

## Layout

```
calm-voi/
├── R/
│   └── calm_voi.R              core library
├── scripts/
│   ├── run_main_grid.R         36-cell main grid (6 DGPs × 3 n × 2 K)
│   ├── run_sensitivity.R       16-cell score-quality sensitivity
│   ├── worked_example.R        250-patient pipeline check
│   └── run_all.R               master driver
├── LICENSE                     MIT
└── README.md
```

## Dependencies

- **R** ≥ 4.0
- Base R packages only (`MASS`, `stats`)

No CRAN dependencies needed for the simulation pipeline.

## Library API

`R/calm_voi.R` exposes:

| Function | Purpose |
|---|---|
| `one_replication(n, K_budget, dgp, sigma_eta_spec)` | One Monte Carlo replication of all 6 estimators |
| `aipw_xfit(Y, A, X, M, ...)` | Cross-fitted AIPW with stratified outcome regression |
| `water_filling(q, total_budget, K_min)` | Solve water-filling for continuous allocation |
| `integer_allocate(K_float, total_budget, K_min)` | Largest-fractional-remainder rounding |
| `make_outcomes(dgp, X, Tstar, A)` | Generate outcomes for DGP $\in$ {"1","2","3"} |
| `make_sigma_eta(spec, n)` | Per-patient noise scale; spec is numeric or "bimodal" |
| `llm_queries(X, Tstar, sigma_eta, K_max)` | Draw the surrogate prediction matrix |

### Minimal example

```r
source("R/calm_voi.R")
set.seed(20260522)

out <- one_replication(
  n              = 200,
  K_budget       = 5,
  dgp            = "3",      # arm-cancellation
  sigma_eta_spec = "1.0"     # homogeneous noise
)

print(sapply(out, `[[`, "tau"))    # tau-hat per estimator
print(sapply(out, `[[`, "se"))     # sandwich SE per estimator
```

## Reproducibility

The pipeline is deterministic.

- `SEED_BASE = 20260522`
- Each cell uses `SEED_BASE + cell_index * 1e6` as its base seed.
- Within a cell, replication $b$ uses `base_seed + b`.

With identical `B` and on the same R version and platform, results are
reproducible to the bit.

## DGPs

All DGPs share a common skeleton:

- $X_i \sim \mathcal{N}(0, \Sigma)$ with $\Sigma_{jk} = 0.3^{|j-k|}$, $p = 5$
- $T_i^* \sim \mathcal{N}(0, 1)$, $A_i \sim \mathrm{Bernoulli}(0.5)$
- $m_i^{(k)} = 2(T_i^* + \eta_i^{(k)}) + X_i^\top\beta_X$ with $\beta_X = (1, 0.5, -0.3, 0, 0)^\top$
- $\eta_i^{(k)} \sim \mathcal{N}(0, \sigma_{\eta, i}^2)$
- True ATE $\tau = 0.4$

DGPs differ only in the arm-specific outcome models and noise specification:

| DGP key | underlying id | $\sigma_{\eta,i}$ | Leverage structure |
|---|---|---|---|
| 1a | 1 | 0.5 | Constant, $B_i \equiv 2$ |
| 1b | 1 | 1.0 | Constant, $B_i \equiv 2$ |
| 1c | 1 | 2.0 | Constant, $B_i \equiv 2$ |
| 1d | 1 | bimodal $\{0.3, 2.0\}$ | Constant, $B_i \equiv 2$ |
| 2  | 2 | 1.0 | Smoothly heterogeneous, $\mathrm{CV}(B_i) \approx 0.25$ |
| 3  | 3 | 1.0 | Arm-cancellation: $B_i \in \{0, 2\}$ |

DGP 3 is the regime where VOI is theoretically distinctive: patients with
$X_{i,2} < 0$ have $\beta_{1i} = -\beta_{0i}$, hence $B_i = 0$, hence VOI
assigns them the minimum $K_i = 1$ and reallocates the freed budget to the
$B_i > 0$ patients.

## Citation

```bibtex
@article{calm_voi_2026,
  title  = {Value-of-Information Allocation for Stochastic Surrogate Covariates in Randomized Trials},
  author = {[Author names blinded for peer review]},
  journal = {Submitted to Statistica Sinica},
  year   = {2026}
}
```

## License

MIT. See `LICENSE`.
