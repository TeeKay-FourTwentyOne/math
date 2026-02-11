"""
Phase 3: Correlation structure analysis of constraint violations.

The make-or-break question for a probabilistic proof: are the bad events
(B_d violates constraint at d) positively or negatively correlated under
a uniformly random D12?

For each prime p ≡ 3 mod 4:
1. Fix a known D11
2. Compute exact covariance matrix of B(d) under random D12
3. Estimate joint probability of all constraints being satisfied
4. Compare union bound, LLL bound, and actual probability

If negative correlation: proof via FKG/second moment
If positive correlation: need LLL dependency graph or abandon probabilistic approach

Key variable: For each d, X_d = Delta(D12, D12, d) = #{a in D12 : (a-d) in D12}.
The binding constraint is: A(d) + X_d <= (p-3)/2 for all d in D11.
Define bad event B_d = {X_d > (p-3)/2 - A(d)}.
"""

import sys
import os
import json
import time
import math
from collections import defaultdict

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import Delta

# Same known solutions as sample_d12_space.py
KNOWN_P3MOD4 = {
    7: {"n": 4, "D11": {1, 2, 5, 6}, "D12": {0, 3, 5}},
    11: {"n": 6, "D11": {1, 2, 4, 7, 9, 10}, "D12": {0, 4, 5, 7, 10}},
    19: {"n": 10, "D11": {1, 2, 3, 6, 8, 11, 13, 16, 17, 18},
         "D12": {0, 2, 6, 8, 9, 12, 13, 17, 18}},
    23: {"n": 12, "D11": {5, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17, 18},
         "D12": {0, 1, 2, 6, 10, 13, 14, 16, 18, 20, 21}},
    31: {"n": 16,
         "D11": {6, 7, 8, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 23, 24, 25},
         "D12": {0, 1, 2, 3, 8, 11, 12, 13, 15, 18, 20, 21, 24, 27, 29}},
    43: {"n": 22,
         "D11": {1, 2, 5, 10, 11, 13, 16, 17, 18, 19, 20, 23, 24, 25, 26, 27,
                 30, 32, 33, 38, 41, 42},
         "D12": {0, 2, 5, 6, 8, 11, 15, 16, 20, 24, 25, 27, 28, 31, 32, 34,
                 35, 36, 37, 39, 41}},
}


def load_registry_solutions():
    """Load additional solutions from registry."""
    registry_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "solutions_registry.json"
    )
    if not os.path.exists(registry_path):
        return
    with open(registry_path) as f:
        registry = json.load(f)
    for sol in registry["solutions"]:
        p = sol["m"]
        if p in {47, 59}:
            D11_raw = set(sol["D11"])
            D12_raw = set(sol["D12"])
            target = (p + 1) // 2
            if len(D11_raw) != target:
                D22 = set(range(1, p)) - D11_raw
                D12T = {(-x) % p for x in D12_raw}
                D11_raw, D12_raw = D22, D12T
            n = (p + 1) // 2
            KNOWN_P3MOD4[p] = {"n": n, "D11": D11_raw, "D12": D12_raw}


def autocorrelation_fft(indicator, m):
    """Compute Delta(S,S,d) for all d via FFT."""
    fft_val = np.fft.fft(indicator)
    autocorr = np.fft.ifft(np.abs(fft_val) ** 2).real
    return np.round(autocorr).astype(int)


def compute_exact_moments(p, D11, num_d, num_mc=100000):
    """Compute E[X_d] and Cov(X_d, X_{d'}) for d, d' in D11 via Monte Carlo.

    X_d = Delta(D12, D12, d) where D12 is a random subset of Z_p of size
    (p-1)/2 containing 0.
    """
    d12_size = (p - 1) // 2
    D11_list = sorted(D11)
    num_constraints = len(D11_list)
    d_to_idx = {d: i for i, d in enumerate(D11_list)}

    rng = np.random.default_rng(42)

    # Collect samples of X_d for all d in D11
    X_samples = np.zeros((num_mc, num_constraints), dtype=np.float64)

    for trial in range(num_mc):
        # Generate random D12 containing 0
        others = rng.choice(range(1, p), size=d12_size - 1, replace=False)
        d12_indicator = np.zeros(p, dtype=np.float64)
        d12_indicator[0] = 1.0
        for x in others:
            d12_indicator[x] = 1.0

        B = autocorrelation_fft(d12_indicator, p)

        for i, d in enumerate(D11_list):
            X_samples[trial, i] = B[d]

    # Compute statistics
    means = X_samples.mean(axis=0)
    cov_matrix = np.cov(X_samples, rowvar=False)

    return D11_list, means, cov_matrix, X_samples


def analyze_prime(p, D11, D12_known, num_mc=100000):
    """Full correlation analysis for a single prime."""
    n = (p + 1) // 2
    d11_size = (p + 1) // 2
    d12_size = (p - 1) // 2
    threshold = (p - 3) // 2  # binding constraint: A(d) + B(d) <= threshold

    # Compute A(d) for the known D11
    d11_indicator = np.zeros(p, dtype=np.float64)
    for j in D11:
        d11_indicator[j] = 1.0
    A = autocorrelation_fft(d11_indicator, p)

    D11_list = sorted(D11)
    num_constraints = len(D11_list)

    # Compute per-constraint slack: how much room does B(d) have?
    # Need A(d) + B(d) <= threshold, so B(d) <= threshold - A(d)
    B_thresholds = {d: threshold - int(A[d]) for d in D11_list}

    print(f"\n  A(d) values for d in D11:")
    A_vals = [int(A[d]) for d in D11_list]
    print(f"    min={min(A_vals)}, max={max(A_vals)}, mean={sum(A_vals) / len(A_vals):.2f}")
    print(f"    B(d) threshold: min={min(B_thresholds.values())}, "
          f"max={max(B_thresholds.values())}")

    # Monte Carlo: sample random D12 and measure correlations
    print(f"\n  Running {num_mc:,} Monte Carlo trials...")
    t0 = time.time()
    D11_list, means, cov_matrix, X_samples = compute_exact_moments(
        p, D11, num_constraints, num_mc=num_mc
    )
    elapsed = time.time() - t0
    print(f"    Completed in {elapsed:.1f}s")

    # Display mean and variance of X_d = B(d) = Delta(D12,D12,d)
    variances = np.diag(cov_matrix)
    print(f"\n  E[B(d)] for d in D11: min={means.min():.3f}, max={means.max():.3f}, "
          f"mean={means.mean():.3f}")
    print(f"  Var[B(d)]: min={variances.min():.3f}, max={variances.max():.3f}, "
          f"mean={variances.mean():.3f}")
    print(f"  Std[B(d)]: {np.sqrt(variances.mean()):.3f}")

    # Theoretical prediction: E[B(d)] ≈ (|D12|-1)*(|D12|-1)/(p-1) for d != 0
    # Exact: E[Delta(D12,D12,d)] where D12 is uniform random with |D12|=k, 0 in D12
    k = d12_size
    # For d != 0: expected number of a in D12 with a-d in D12.
    # If 0 ∈ D12 always, and D12 has k elements, the k-1 non-zero elements are
    # chosen uniformly from {1,...,p-1}.
    # For a = 0: count [(-d) in D12\{0}] = Pr[-d in D12\{0}] = (k-1)/(p-1)
    # For a != 0, a in D12: count [(a-d) in D12] = [(a-d)=0]*Pr[a=d] + ...
    # This is getting complicated. Just use MC means.
    print(f"  Theoretical E[B(d)] ≈ {(k - 1) * (k - 1) / (p - 1):.3f} "
          f"(crude estimate)")

    # Correlation matrix
    std_devs = np.sqrt(np.maximum(variances, 1e-10))
    corr_matrix = cov_matrix / np.outer(std_devs, std_devs)
    np.fill_diagonal(corr_matrix, 1.0)

    # Extract off-diagonal correlations
    upper_tri = corr_matrix[np.triu_indices(num_constraints, k=1)]
    print(f"\n  CORRELATION STRUCTURE ({num_constraints} constraints):")
    print(f"    Off-diagonal correlations:")
    print(f"      min={upper_tri.min():.4f}, max={upper_tri.max():.4f}, "
          f"mean={upper_tri.mean():.4f}")
    print(f"      Fraction positive: {(upper_tri > 0).sum()}/{len(upper_tri)} "
          f"({100 * (upper_tri > 0).mean():.1f}%)")
    print(f"      Fraction > 0.1: {(upper_tri > 0.1).sum()}/{len(upper_tri)}")
    print(f"      Fraction < -0.1: {(upper_tri < -0.1).sum()}/{len(upper_tri)}")

    mean_corr = upper_tri.mean()
    if mean_corr > 0.05:
        print(f"    >>> POSITIVE CORRELATION (mean={mean_corr:.4f})")
        print(f"    >>> LLL or dependency graph needed")
    elif mean_corr < -0.05:
        print(f"    >>> NEGATIVE CORRELATION (mean={mean_corr:.4f})")
        print(f"    >>> FKG/second moment method may work")
    else:
        print(f"    >>> NEAR ZERO CORRELATION (mean={mean_corr:.4f})")
        print(f"    >>> Events are approximately independent")

    # Joint probability estimation
    # Count fraction of MC trials where ALL constraints are satisfied
    all_satisfied = np.ones(num_mc, dtype=bool)
    per_constraint_satisfied = np.zeros(num_constraints, dtype=int)
    for i, d in enumerate(D11_list):
        bt = B_thresholds[d]
        satisfied = X_samples[:, i] <= bt
        per_constraint_satisfied[i] = satisfied.sum()
        all_satisfied &= satisfied

    joint_prob = all_satisfied.sum() / num_mc
    print(f"\n  JOINT PROBABILITY ANALYSIS:")
    print(f"    Per-constraint satisfaction rates:")
    rates = per_constraint_satisfied / num_mc
    print(f"      min={rates.min():.4f}, max={rates.max():.4f}, "
          f"mean={rates.mean():.4f}")

    # Independence prediction: product of per-constraint rates
    if rates.min() > 0:
        log_indep = sum(np.log(rates))
        indep_pred = np.exp(log_indep)
    else:
        indep_pred = 0.0

    print(f"    Joint satisfaction (MC): {all_satisfied.sum()}/{num_mc} "
          f"= {joint_prob:.6f}")
    print(f"    Independence prediction: {indep_pred:.6f}")
    if indep_pred > 0 and joint_prob > 0:
        ratio = joint_prob / indep_pred
        print(f"    Ratio (actual/independent): {ratio:.4f}")
        if ratio > 1.5:
            print(f"    >>> POSITIVE correlation HELPS (events tend to co-satisfy)")
        elif ratio < 0.67:
            print(f"    >>> POSITIVE correlation HURTS (events tend to co-fail)")
        else:
            print(f"    >>> Near independence")

    # Union bound comparison
    per_violation = 1.0 - rates
    union_bound = min(1.0, per_violation.sum())
    print(f"\n    Union bound on P[any violation]: {union_bound:.4f}")
    print(f"    Actual P[any violation]: {1 - joint_prob:.4f}")
    print(f"    Union bound gap: {union_bound - (1 - joint_prob):.4f}")

    # LLL check: is max individual probability < 1/(e*max_degree)?
    # Dependency graph: constraints at d and d' are dependent if they share
    # elements in their Delta computation
    max_individual_fail = per_violation.max()
    # Each constraint depends on all others (shared D12), so degree = num_constraints-1
    lll_threshold = 1.0 / (math.e * num_constraints)
    print(f"\n    LLL analysis:")
    print(f"      Max P[violation at single d]: {max_individual_fail:.4f}")
    print(f"      LLL threshold (1/eD): {lll_threshold:.4f}")
    print(f"      LLL applicable (naive)? "
          f"{'YES' if max_individual_fail < lll_threshold else 'NO'}")

    # MVN approximation
    print(f"\n  MULTIVARIATE NORMAL APPROXIMATION:")
    if cov_matrix.shape[0] <= 50:
        try:
            # Center the thresholds
            centered_thresholds = np.array([
                B_thresholds[d] - means[i] for i, d in enumerate(D11_list)
            ]) / std_devs

            # P[all X_i <= t_i] under MVN with correlation matrix
            # Use Monte Carlo with the estimated correlation structure
            mvn_samples = np.random.default_rng(99).multivariate_normal(
                np.zeros(num_constraints), corr_matrix, size=100000
            )
            mvn_all_below = np.all(
                mvn_samples <= centered_thresholds[np.newaxis, :], axis=1
            )
            mvn_prob = mvn_all_below.mean()
            print(f"    MVN joint satisfaction: {mvn_prob:.6f}")
            print(f"    MC joint satisfaction:  {joint_prob:.6f}")
            if mvn_prob > 0:
                print(f"    MVN/MC ratio: {mvn_prob / max(joint_prob, 1e-10):.4f}")
        except Exception as e:
            print(f"    MVN computation failed: {e}")
    else:
        print(f"    Skipped (too many constraints: {num_constraints})")

    return {
        "p": p,
        "n": n,
        "num_constraints": num_constraints,
        "mean_B": float(means.mean()),
        "std_B": float(np.sqrt(variances.mean())),
        "mean_corr": float(mean_corr),
        "joint_prob_mc": float(joint_prob),
        "indep_prediction": float(indep_pred),
        "union_bound": float(union_bound),
        "lll_applicable": bool(max_individual_fail < lll_threshold),
        "per_constraint_rates": rates.tolist(),
    }


def main():
    load_registry_solutions()

    # MC sample sizes per prime
    mc_sizes = {7: 500000, 11: 500000, 19: 200000, 23: 200000,
                31: 100000, 43: 50000, 47: 50000, 59: 30000}

    all_results = {}
    for p in sorted(KNOWN_P3MOD4.keys()):
        info = KNOWN_P3MOD4[p]
        D11 = info["D11"]
        D12 = info["D12"]
        n = info["n"]

        print(f"\n{'=' * 80}")
        print(f"CORRELATION ANALYSIS: p = {p}, n = {n}")
        print(f"{'=' * 80}")

        num_mc = mc_sizes.get(p, 30000)
        result = analyze_prime(p, D11, D12, num_mc=num_mc)
        all_results[str(p)] = result

    # Cross-prime summary
    print(f"\n{'=' * 80}")
    print("CORRELATION SUMMARY ACROSS PRIMES")
    print(f"{'=' * 80}")
    print(f"{'p':>4} {'n':>3} {'#constr':>7} {'mean_corr':>10} "
          f"{'joint_MC':>10} {'indep':>10} {'ratio':>8} {'LLL':>5}")
    print("-" * 70)
    for p_str in sorted(all_results.keys(), key=int):
        r = all_results[p_str]
        p = r["p"]
        ratio = (r["joint_prob_mc"] / r["indep_prediction"]
                 if r["indep_prediction"] > 0 else float("inf"))
        print(f"{p:4d} {r['n']:3d} {r['num_constraints']:7d} "
              f"{r['mean_corr']:10.4f} {r['joint_prob_mc']:10.6f} "
              f"{r['indep_prediction']:10.6f} {ratio:8.3f} "
              f"{'Y' if r['lll_applicable'] else 'N':>5}")

    # Key insight: how does the correlation structure evolve with p?
    primes_sorted = sorted(all_results.keys(), key=int)
    if len(primes_sorted) >= 3:
        corrs = [all_results[p]["mean_corr"] for p in primes_sorted]
        primes_vals = [int(p) for p in primes_sorted]
        print(f"\n  Correlation trend: {list(zip(primes_vals, [f'{c:.4f}' for c in corrs]))}")
        if all(c < 0.05 for c in corrs):
            print("  >>> Correlations are uniformly small → independence assumption viable")
        elif corrs[-1] < corrs[0]:
            print("  >>> Correlations DECREASE with p → probabilistic proof gets easier")
        else:
            print("  >>> Correlations INCREASE with p → need careful analysis")

    # Joint probability trend
    print(f"\n  Joint probability trend:")
    for p_str in primes_sorted:
        r = all_results[p_str]
        print(f"    p={r['p']:3d}: P[all ok] = {r['joint_prob_mc']:.6f}, "
              f"E[B]={r['mean_B']:.2f}, std={r['std_B']:.2f}")

    # Save
    output_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "correlation_results.json"
    )
    with open(output_path, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nResults saved to {output_path}")


if __name__ == "__main__":
    main()
