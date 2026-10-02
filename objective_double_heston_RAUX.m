function [y, betamean, stdbeta, chi2_contrib, t_stats] = objective_double_heston_RAUX(coeffs, par, invVCV, W)
% OBJECTIVE_DOUBLE_HESTON_RAUX  SMM objective for Double Heston with the
% AUGMENTED auxiliary (RAUX, 11 moments) instead of the LHAR (8 moments).
% Drop-in replacement for objective_double_heston_LHAR_new: identical
% Monte-Carlo machinery (per-replication Threefry seeds, common random
% numbers + antithetics, parfor, feasibility guard, chi^2 decomposition);
% the ONLY change is that each simulated path is summarized by
% rough_aux_estimate (HAR slopes + log-RV variogram + realized-semivariance /
% lag-1 leverage + return skew/kurtosis) rather than LHAR_estimate.
%
% MOTIVATION: with the LHAR (8 moments) the Double Heston is JUST-IDENTIFIED
% (8 params = 8 moments), so chi^2 has 0 degrees of freedom and the apparent
% in-sample fit is illusory (it collapses out-of-sample). RAUX gives 11
% moments -> 3 over-identifying restrictions, a genuine J-test, and the same
% moment set as the rough Heston pipeline for a clean cross-model comparison.
% The vg_slope moment binds no structural parameter here (no roughness in a
% classical 2-factor model); it acts as a specification test of whether two
% smooth CIR factors can reproduce the empirical short-lag scaling of log-RV.
%
% INPUTS:
%   coeffs - (11 x 1) target RAUX moments (empirical, from rough_aux_estimate)
%   par    - 1 x 8 row vector [kappa1, kappa2, Vbar1, Vbar2, sigma1, sigma2,
%            rho1, rho2]
%   invVCV - (11 x 11) inverse moment covariance (block bootstrap), the
%            efficient weighting matrix
%   W      - struct with .seeds, .n_intra, .n_sim, .n_repl
%
% OUTPUTS (same contract as objective_double_heston_LHAR_new):
%   y            - SMM weighted quadratic loss betamean' * invVCV * betamean
%   betamean     - (11 x 1) mean (simulated - empirical) moment gap
%   stdbeta      - (11 x 1) MC standard error of the mean across replications
%   chi2_contrib - (11 x 1) per-moment chi^2 decomposition b .* (W*b)
%   t_stats      - (11 x 1) pseudo t-stat betamean ./ stdbeta

PENALTY = 1e10;

p.kappa1 = par(1);
p.kappa2 = par(2);
p.Vbar1  = par(3);
p.Vbar2  = par(4);
p.sigma1 = par(5);
p.sigma2 = par(6);
p.rho1   = par(7);
p.rho2   = par(8);
p.mu     = 0;

nintra = W.n_intra;
ndays  = W.n_sim;
nrepl  = W.n_repl;
seeds  = W.seeds;
nmom   = length(coeffs);

% ------------------ feasibility ------------------
% Same thresholds as objective_double_heston_LHAR_new (2026-05-22 alignment):
% f1, f2 >= 1e-4 keeps each variance well-defined; sigma_k/kappa_k <= 500
% rules out the pathological random-walk regime.
f1 = 2 * p.kappa1 * p.Vbar1 / (p.sigma1^2);
f2 = 2 * p.kappa2 * p.Vbar2 / (p.sigma2^2);
if f1 < 1e-4 || f2 < 1e-4 || p.sigma1/p.kappa1 > 500 || p.sigma2/p.kappa2 > 500
    y = PENALTY;
    betamean     = NaN(nmom, 1);
    stdbeta      = NaN(nmom, 1);
    chi2_contrib = NaN(nmom, 1);
    t_stats      = NaN(nmom, 1);
    return
end

n_repl_total = 2 * nrepl;
beta = NaN(nmom, n_repl_total);

parfor j = 1:n_repl_total
    beta_j = NaN(nmom, 1);
    try
        if j <= nrepl
            j_idx = j;       s_val = +1;
        else
            j_idx = j - nrepl; s_val = -1;
        end

        % Common random numbers + antithetics: same seed -> same path across
        % optimizer iterations; s_val flips the sign for the antithetic wave.
        rs = RandStream('Threefry', 'Seed', seeds(j_idx));
        w1j = s_val * randn(rs, nintra, ndays);
        w2j = s_val * randn(rs, nintra, ndays);
        w3j = s_val * randn(rs, nintra, ndays);
        w4j = s_val * randn(rs, nintra, ndays);

        r = simulate_double_heston_new(ndays, nintra, p, w1j, w2j, w3j, w4j);

        RV       = sum(r.^2, 1) * 252;
        dailyret = sum(r, 1)    * sqrt(252);

        if all(isfinite(RV)) && all(RV > 0) && all(isfinite(dailyret))
            logRV = log(RV);
            if all(isfinite(logRV))
                % r is (intervals x days); rough_aux_estimate wants
                % (days x intervals), so pass the transpose. logRV/dailyret
                % are rows -> pass as columns.
                result = rough_aux_estimate(logRV', dailyret', r');
                if all(isfinite(result.moments))
                    beta_j = result.moments - coeffs;
                end
            end
        end
    catch
        % beta_j stays NaN
    end
    beta(:, j) = beta_j;
end

valid_cols = all(isfinite(beta), 1);
n_valid = sum(valid_cols);

if n_valid < max(2, 0.25 * n_repl_total)
    y = PENALTY;
    betamean = NaN(nmom, 1);
    stdbeta  = NaN(nmom, 1);
    chi2_contrib = NaN(nmom, 1);
    t_stats      = NaN(nmom, 1);
    return
end

betamean = mean(beta(:, valid_cols), 2);
stdbeta  = std(beta(:, valid_cols), 0, 2) / sqrt(n_valid);

if any(~isfinite(betamean))
    y = PENALTY;
    chi2_contrib = NaN(nmom, 1);
    t_stats      = NaN(nmom, 1);
    return
end

y = betamean' * invVCV * betamean;
if ~isfinite(y)
    y = PENALTY;
end

% Per-moment chi^2 decomposition. With W non-diagonal individual entries can
% be negative; their magnitude still flags which moment is being sacrificed.
Wb = invVCV * betamean;
chi2_contrib = betamean .* Wb;
t_stats      = betamean ./ max(stdbeta, 1e-12);

end
