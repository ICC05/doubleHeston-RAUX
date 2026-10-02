# Double Heston — indirect inference estimation with the RAUX auxiliary criterion

MATLAB code for the **indirect inference** (EMSM) estimation of the two-factor **double Heston** model on 5-minute returns of the SPDR S&P 500 ETF (SPY). The auxiliary criterion is **RAUX**, a vector of 11 statistics summarizing the persistence, roughness, leverage and tails of realized volatility.

The code accompanies the research project *Indirect inference estimation of the rough Heston model* (I. Cardosi Carrara, Sant'Anna School of Advanced Studies, 2026), Section 7.4 of the project report. The project regards this as its **decisive** estimation.

---

## 1. Why this repository

Under the LHAR criterion (8 moments), the double Heston (8 parameters) is **just-identified** (q = p):
- the system b(θ) = β̂ generically has an exact solution, so a near-zero χ² (7.03, see [Heston-doubleHeston-roughHeston-LHAR](https://github.com/ICC05/Heston-doubleHeston-roughHeston-LHAR)) can arise "by construction";
- indeed the **out-of-sample χ² explodes** to 174.6, a sign that part of the fit depended on the Monte Carlo seeds.

To turn that result into a **falsifiable test**, the model is re-estimated with RAUX: 11 moments for 8 parameters, i.e. **3 over-identifying restrictions**. A low χ² can then no longer be an artifact of the design.

There is a second motivation: moments and weighting matrix are **bit-for-bit identical** to those of the rough Heston RAUX estimation ([roughHeston-RAUX](https://github.com/ICC05/roughHeston-RAUX)), so the two models are compared on the same metric.

## 2. The idea: indirect inference

1. Auxiliary statistics on the **observed data** → β̂.
2. For each candidate θ: **simulate** many paths, compute the same statistics, average over replications → simulated *binding function* b̂(θ).
3. Minimize χ²(θ) = [β̂ − b̂(θ)]ᵀ Ω [β̂ − b̂(θ)], where Ω is the inverse covariance of the empirical moments.

The random shocks are held fixed across candidate θ (common random numbers).

## 3. Data

SPY, 1990–2013, cleaned 5-minute returns: low-activity days and outliers removed, August–December 2008 excluded. The final sample has **T = 5,867 days** (`spy_data.mat`). Daily realized variance is the sum of squared intraday returns, excluding overnight.

## 4. Structural model: double Heston (Christoffersen et al., 2010)

dSₜ/Sₜ = μ dt + √v₁,ₜ dW₁,ₜ + √v₂,ₜ dW₂,ₜ

dvₙ,ₜ = κₙ(v̄ₙ − vₙ,ₜ) dt + σₙ √vₙ,ₜ dZₙ,ₜ,  dWₙ,ₜ dZₙ,ₜ = ρₙ dt,  n = 1, 2

The two factors are mutually independent. The parameters are **θ = (κ₁, κ₂, v̄₁, v̄₂, σ₁, σ₂, ρ₁, ρ₂)**, so p = 8. The idea is that a fast and a slow factor, added together, mimic the near-hyperbolic decay of volatility autocorrelations.

### Simulation
- Grid of **M = 80 intraday steps** (matching the 80 five-minute bars in the data) over **6,000 days**.
- Each CIR factor is propagated with the **exact conditional-moment-matching** scheme (Corsi and Renò, 2012) with full-truncation correction. Returns are the sum of the two contributions.
- Each factor starts at its long-run mean, so no burn-in is needed.
- **H = 494 effective replications**: 247 independent streams plus antithetics. Shocks are regenerated from 247 fixed Threefry seeds and replications are evaluated in parallel.
- Safeguards:
  - candidates with Feller ratio < 10⁻⁴ or σₙ/κₙ > 500 are penalized without simulating them;
  - invalid replications are discarded.

## 5. RAUX auxiliary criterion (q = 11)

A vector of statistics computed with the same routine on the data and on the simulations:

| Block | Statistics |
|---|---|
| Level and persistence (HAR-RV) | mean of log RV; β⁽ᵈ⁾, β⁽ʷ⁾, β⁽ᵐ⁾; residual variance |
| Roughness | slope and intercept of the log-RV variogram (lags 2, 5, 10, 22 days) |
| Leverage | corr(rₜ, log RVₜ₊₁); corr((RS⁻ − RS⁺)/RV, log RVₜ₊₁) from realized semivariances |
| Tail shape | skewness and excess kurtosis of daily returns |

For the double Heston, a smooth model with no roughness parameter, the variogram block acts as a **pure specification test**: it checks whether two smooth factors can mimic the apparent roughness of realized volatility.

**Weighting matrix**: covariance from a **moving-block bootstrap** (B = 1,000, 44-day blocks, joint resampling of the series, ridge for invertibility). It uses the same seeds and settings as the rough Heston pipeline.

## 6. Optimization and diagnostics

1. **CMA-ES**: population 20, up to 300 generations, 10 restarts with Latin Hypercube initialization. Restarts with χ² ≥ 5× the best are relaunched from an antithetic point (at most 3).
2. **Bound-constrained pattern search**, from 5 starting points.

Bounds: κₙ ∈ [0.001, 50], v̄ₙ ∈ [0.01, 3], σₙ ∈ [0.01, 10], ρₙ ∈ [−0.999, 0.999].

Diagnostics:
- **out-of-sample** χ² under independent seeds;
- **noise floor** over 10 seeds;
- **cross-restart dispersion of ρ₁ and ρ₂**;
- **per-moment** χ² decomposition.

## 7. Results (from the project report)

| Parameter | Estimate | | RAUX moment | Observed | Simulated |
|---|---|---|---|---|---|
| κ₁ | 2.7082 | | mean log RV | 4.650 | 4.666 |
| κ₂ | 0.0130 | | β⁽ᵈ⁾ | 0.340 | 0.341 |
| v̄₁ | 0.2177 | | β⁽ʷ⁾ | 0.389 | 0.413 |
| v̄₂ | 0.4175 | | β⁽ᵐ⁾ | 0.218 | 0.185 |
| σ₁ | 1.1845 | | variogram slope | 0.119 | 0.140 |
| σ₂ | 0.1685 | | leverage corr | −0.125 | −0.085 |
| ρ₁ | 0.1601 | | skewness | −0.185 | −0.138 |
| ρ₂ | −0.9990 | | excess kurtosis | 4.573 | 3.215 |

**In-sample χ² = 11.84, out-of-sample χ² = 12.13** (ratio 1.02). The noise floor is 12.3 ± 0.40.

- **Central result of the project**: on an over-identified metric the double Heston drives the distance down to the Monte Carlo noise level and keeps it there out of sample. The rough Heston, on the same 11 moments, attains **χ² = 483.87**, about **40 times worse**.
- **Mechanism**: a fast factor (half-life ≈ 0.26 days) and a very persistent one (half-life ≈ 53 days). Leverage is carried mainly by the slow factor, with ρ₂ ≈ −1.
- **Roughness**: the variogram slope is nearly reproduced (0.140 vs 0.119) **without any rough driver**. This is consistent with the literature arguing that measured roughness is partly an artifact of measurement noise.
- **Only sizeable misfit**: kurtosis is underestimated (3.2 vs 4.6). This points to adding **jumps** (jump-augmented multifactor model) as the natural extension.
- **Caveat**: the cross-restart dispersion of ρ₁ and ρ₂ is large (≈ 1.2), due to the factor-exchange symmetry of the model. The result concerns the **reproduction** of the stylized facts, not the sharp **identification** of the parameters, which is left for future work.

## 8. How to run

Requirements: base MATLAB plus Parallel Computing Toolbox (recommended).

From the repository folder, run in MATLAB:

```matlab
estimate_double_heston_RAUX_new2
```

Output in the current folder:
- `results_double_heston_RAUX_new2.mat`: struct `est` plus `coeffs_aux`, `VCV_aux` and more;
- `results_double_heston_RAUX_new2.txt`: results table;
- `logfile_double_heston_RAUX_new2.txt`: full log.

The number of replications (`cfg.n_repl = 247` in `setup_data_new.m`) is tuned for a 247-worker pool. On a standard PC it is worth lowering it, keeping in mind that Monte Carlo noise increases.

A description of every file is in [`CODE.md`](CODE.md).

## Main references
- Christoffersen, P., Jacobs, K. & Mimouni, K. (2010). *Volatility dynamics for the S&P 500*. RFS.
- Corsi, F. & Renò, R. (2012). JBES.
- Gouriéroux, C., Monfort, A. & Renault, E. (1993). *Indirect inference*. JAE.
- Patton, A. J. & Sheppard, K. (2015). *Good volatility, bad volatility*. REStat.
