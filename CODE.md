# File guide

For the conceptual overview see [`README.md`](README.md).

## Estimation script (entry point)
| File | What it does |
|---|---|
| `estimate_double_heston_RAUX_new2.m` | Computes the empirical RAUX target and the bootstrap weighting matrix, then estimates the double Heston (8 parameters): multistart CMA-ES with adaptive relaunch, pattern search, out-of-sample check and noise floor; saves the results. |

## Setup
| File | What it does |
|---|---|
| `setup_data_new.m` | Loads `spy_data.mat`, builds RV and daily returns, estimates the empirical HAR/LHAR (diagnostic), generates the replication seeds and starts the parallel pool. |

## Objective function and simulator
| File | What it does |
|---|---|
| `objective_double_heston_RAUX.m` | Given θ: simulates the double Heston replications, computes RAUX on each, returns χ², mean simulated moments, per-moment contributions and t-stats. |
| `simulate_double_heston_new.m` | Simulates intraday returns with two independent CIR factors (exact-moment scheme, full truncation). |

## RAUX auxiliary criterion
| File | What it does |
|---|---|
| `rough_aux_estimate.m` | Computes the 11 RAUX statistics (HAR, variogram, leverage, skewness/kurtosis). The "rough" in the name is historical: the criterion is the same one used for the rough Heston. |
| `rough_aux_bootstrap_vcv.m` | Covariance of the RAUX moments via moving-block bootstrap (with ridge) → weighting matrix. |
| `HAR_estimate.m` | HAR-RV regression (used inside `rough_aux_estimate` and in the setup). |
| `LHAR_estimate.m` | LHAR regression (computed in the setup only as a reference). |
| `aggregateAvg.m` | Backward-looking rolling average (weekly/monthly regressors). |
| `nwest.m` | OLS with Newey-West covariance. |

## Optimizers
| File | What it does |
|---|---|
| `cmaes_bnd_new.m` | Box-constrained CMA-ES with per-generation relative stopping rule. |
| `pattern_search_bnd_new.m` | Bound-constrained pattern search with Latin Hypercube polls and plateau stopping. |

## Output and logging
| File | What it does |
|---|---|
| `print_estimation_results_new.m` | Prints and saves the results table (observed vs simulated moments, t-stats, per-moment χ²). |
| `dual_log_new.m` | Writes messages to screen and to the logfile. |

## Data and results
| File | Content |
|---|---|
| `spy_data.mat` | Cleaned SPY 5-minute returns (days × intervals); the only data file read by the code. |
| `results_double_heston_RAUX_new2.mat / .txt` | Estimation results. |
| `logfile_double_heston_RAUX_new2.txt` | Estimation log. |
