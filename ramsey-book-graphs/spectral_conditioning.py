#!/usr/bin/env python3
"""
Spectral Conditioning Analysis: Why are B-events positively associated?

HYPOTHESIS: The positive association arises because the "spectral quality" of D12
acts as a latent variable. Specifically:

B(d) = (1/p) Σ_k |D̂₁₂(k)|² e^{-2πikd/p}

If we define Q(k) = |D̂₁₂(k)|² for k ≠ 0, then:
- B(d) = (|D12|²/p) + (1/p) Σ_{k≠0} Q(k) e^{-2πikd/p}
- Σ_k Q(k) = p|D12| (Parseval)
- Q(k) = Q(p-k) (conjugate symmetry)

"Spectrally flat" D12 (all Q(k) ≈ |D12|) → all B(d) ≈ mean → easy to satisfy
"Spectrally peaked" D12 (some Q(k) large) → some B(d) large → violations

Key test: if we condition on spectral flatness (max Q(k)), is the positive
association explained (i.e., does it vanish or weaken)?

This script:
1. Computes the spectral profile Q(k) for random D12
2. Measures correlation between spectral flatness and constraint satisfaction
3. Tests whether conditioning on spectral flatness explains the positive association
4. Identifies the "spectral signature" of valid D12
"""

import numpy as np
import time
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def autocorrelation_fft(indicator, p):
    fft_val = np.fft.fft(indicator)
    autocorr = np.fft.ifft(np.abs(fft_val) ** 2).real
    return np.round(autocorr).astype(int)


KNOWN_SOLUTIONS = {
    11: {"D11": {1, 2, 4, 7, 9, 10}},
    19: {"D11": {1, 2, 3, 6, 8, 11, 13, 16, 17, 18}},
    23: {"D11": {5, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17, 18}},
    31: {"D11": {6, 7, 8, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 23, 24, 25}},
}


def spectral_analysis(p, num_mc=200000):
    """Analyze spectral profile of D12 and its relation to constraint satisfaction."""
    n = (p + 1) // 2
    d12_size = (p - 1) // 2
    threshold = (p - 3) // 2
    k_nonzero = d12_size - 1

    D11 = sorted(KNOWN_SOLUTIONS[p]["D11"])

    # Compute A(d)
    d11_ind = np.zeros(p)
    for d in D11:
        d11_ind[d] = 1.0
    A = autocorrelation_fft(d11_ind, p)
    T = {d: threshold - int(A[d]) for d in D11}

    print(f"\n{'='*80}")
    print(f"SPECTRAL CONDITIONING: p={p}, n={n}")
    print(f"{'='*80}")

    rng = np.random.default_rng(42)

    # Collect data
    max_Q_vals = []   # spectral peak
    total_Q2_vals = []  # Σ Q(k)² (spectral energy concentration)
    all_ok_flags = []
    marginal_ok_per_d = {d: [] for d in D11}

    t0 = time.time()
    for trial in range(num_mc):
        S = rng.choice(range(1, p), size=k_nonzero, replace=False)
        d12_ind = np.zeros(p, dtype=complex)
        d12_ind[0] = 1.0
        for x in S:
            d12_ind[x] = 1.0

        # Compute spectrum
        fft_val = np.fft.fft(d12_ind.real)
        Q = np.abs(fft_val) ** 2  # |D̂₁₂(k)|²

        # Spectral statistics (k ≠ 0)
        Q_nonzero = Q[1:]
        max_Q = Q_nonzero.max()
        total_Q2 = np.sum(Q_nonzero ** 2)

        max_Q_vals.append(max_Q)
        total_Q2_vals.append(total_Q2)

        # Compute B(d) from autocorrelation
        B_vals = np.fft.ifft(Q).real
        B = np.round(B_vals).astype(int)

        all_ok = True
        for d in D11:
            ok = (B[d] <= T[d])
            marginal_ok_per_d[d].append(ok)
            if not ok:
                all_ok = False
        all_ok_flags.append(all_ok)

    elapsed = time.time() - t0

    max_Q_arr = np.array(max_Q_vals)
    total_Q2_arr = np.array(total_Q2_vals)
    all_ok_arr = np.array(all_ok_flags)

    print(f"  {num_mc:,} trials in {elapsed:.1f}s")
    print(f"  Overall Pr[all ok] = {all_ok_arr.mean():.6f}")

    # E[Q(k)] for k ≠ 0 = p|D12| - |D12|² / ...
    # Actually Σ Q(k) = p|D12|, Q(0) = |D12|², so Σ_{k≠0} Q(k) = p|D12| - |D12|²
    E_Q = (p * d12_size - d12_size ** 2) / (p - 1)
    print(f"  E[Q(k)] for k≠0 = {E_Q:.2f}")
    print(f"  Mean max Q(k): {max_Q_arr.mean():.2f}")
    print(f"  Mean Σ Q²: {total_Q2_arr.mean():.2f}")

    # Correlation between spectral flatness and constraint satisfaction
    print(f"\n  SPECTRAL FLATNESS vs CONSTRAINT SATISFACTION:")

    # Split by max_Q quantiles
    quantiles = np.percentile(max_Q_arr, [25, 50, 75])
    labels = ["Q1 (flattest)", "Q2", "Q3", "Q4 (peakiest)"]

    for i, (lo, hi) in enumerate(zip(
        [-1, quantiles[0], quantiles[1], quantiles[2]],
        [quantiles[0], quantiles[1], quantiles[2], max_Q_arr.max() + 1]
    )):
        mask = (max_Q_arr > lo) & (max_Q_arr <= hi)
        if mask.sum() == 0:
            continue
        ok_rate = all_ok_arr[mask].mean()
        print(f"    {labels[i]:20s} (max_Q ∈ ({lo:.1f}, {hi:.1f}]): "
              f"Pr[all ok] = {ok_rate:.6f} (n={mask.sum():,})")

    # Conditional positive association test
    # For each spectral bin, compute ratio joint/indep
    print(f"\n  CONDITIONAL RATIO (joint/indep) by spectral quartile:")
    for i, (lo, hi) in enumerate(zip(
        [-1, quantiles[0], quantiles[1], quantiles[2]],
        [quantiles[0], quantiles[1], quantiles[2], max_Q_arr.max() + 1]
    )):
        mask = (max_Q_arr > lo) & (max_Q_arr <= hi)
        if mask.sum() == 0:
            continue

        ok_in_bin = all_ok_arr[mask]
        joint_rate = ok_in_bin.mean()

        # Per-constraint rates in this bin
        marginal_rates = []
        for d in D11:
            d_ok = np.array(marginal_ok_per_d[d])[mask]
            marginal_rates.append(d_ok.mean())

        if all(r > 0 for r in marginal_rates) and joint_rate > 0:
            log_indep = sum(np.log(r) for r in marginal_rates)
            indep_pred = np.exp(log_indep)
            ratio = joint_rate / indep_pred if indep_pred > 0 else float('inf')
            print(f"    {labels[i]:20s}: joint={joint_rate:.6f}, "
                  f"indep={indep_pred:.6f}, ratio={ratio:.4f}")
        else:
            print(f"    {labels[i]:20s}: insufficient data")

    # Correlation between max_Q and number of satisfied constraints
    num_ok_per_trial = np.zeros(num_mc)
    for d in D11:
        num_ok_per_trial += np.array(marginal_ok_per_d[d], dtype=float)

    corr = np.corrcoef(max_Q_arr, num_ok_per_trial)[0, 1]
    print(f"\n  Corr(max_Q, #constraints satisfied) = {corr:.4f}")

    corr2 = np.corrcoef(total_Q2_arr, num_ok_per_trial)[0, 1]
    print(f"  Corr(Σ Q², #constraints satisfied) = {corr2:.4f}")

    # Spectral profile of valid D12
    print(f"\n  SPECTRAL PROFILE OF VALID D12:")
    if all_ok_arr.sum() > 0:
        valid_max_Q = max_Q_arr[all_ok_arr]
        invalid_max_Q = max_Q_arr[~all_ok_arr]
        print(f"    Valid D12:   max_Q mean={valid_max_Q.mean():.2f}, "
              f"std={valid_max_Q.std():.2f}")
        print(f"    Invalid D12: max_Q mean={invalid_max_Q.mean():.2f}, "
              f"std={invalid_max_Q.std():.2f}")
        print(f"    Valid D12 are {(invalid_max_Q.mean() - valid_max_Q.mean()) / invalid_max_Q.std():.2f} "
              f"sigma flatter")

        valid_Q2 = total_Q2_arr[all_ok_arr]
        invalid_Q2 = total_Q2_arr[~all_ok_arr]
        print(f"    Valid D12:   Σ Q² mean={valid_Q2.mean():.2f}")
        print(f"    Invalid D12: Σ Q² mean={invalid_Q2.mean():.2f}")
    else:
        print(f"    No valid D12 found")

    return {
        "p": p,
        "overall_rate": float(all_ok_arr.mean()),
        "corr_maxQ": float(corr),
        "corr_Q2": float(corr2),
    }


def main():
    print("=" * 80)
    print("SPECTRAL CONDITIONING ANALYSIS")
    print("Why are B-events positively associated for fixed D11?")
    print("=" * 80)

    results = {}
    for p in [11, 19, 23, 31]:
        results[p] = spectral_analysis(p, num_mc=200000)

    # Summary
    print(f"\n\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    print(f"\n  {'p':>4s} {'Pr[all ok]':>12s} {'Corr(maxQ)':>12s} {'Corr(ΣQ²)':>12s}")
    print(f"  {'-'*45}")
    for p in sorted(results.keys()):
        r = results[p]
        print(f"  {p:4d} {r['overall_rate']:12.6f} {r['corr_maxQ']:12.4f} "
              f"{r['corr_Q2']:12.4f}")

    print(f"\n  INTERPRETATION:")
    print(f"  If Corr(maxQ, #ok) < -0.2: spectral flatness explains positive association")
    print(f"  If conditional ratios remain > 1: positive association is INTRINSIC to B(d)")
    print(f"  Both findings together → the mechanism is spectral quality as latent variable")


if __name__ == "__main__":
    main()
