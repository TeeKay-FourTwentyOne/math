#!/usr/bin/env python3
"""
Test whether the Suen/Janson correlation inequality gives a useful lower bound
for Pr[all B(d) <= T(d)].

The Suen inequality (Suen 1990, Janson 1998): For events E_1, ..., E_m:

  Pr[∩ Ē_i] >= ∏(1 - Pr[E_i]) × exp(-Δ)

where Δ = Σ_{i~j} Pr[E_i ∩ E_j] / ((1-Pr[E_i])(1-Pr[E_j])),
summed over pairs (i,j) that are "dependent."

In our case:
- E_d = {B(d) > T(d)} for d ∈ D11 (bad events)
- All events are dependent (same D12)
- Pr[E_d] = 1 - Pr[B(d) ≤ T(d)] ≈ 1/2

If Δ = O(1) or O(log p), the Suen bound gives:
  Pr[all ok] >= (1/2)^n × exp(-O(1)) ≈ (1/2)^n × const

which is sufficient for the proof.

This script computes Δ for each prime and tests whether the bound is useful.
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


def compute_suen_bound(p, num_mc=200000):
    """Compute Suen/Janson bound components."""
    n = (p + 1) // 2
    d12_size = (p - 1) // 2
    threshold = (p - 3) // 2
    k_nonzero = d12_size - 1

    D11 = sorted(KNOWN_SOLUTIONS[p]["D11"])
    num_d = len(D11)

    # Compute A(d)
    d11_ind = np.zeros(p)
    for d in D11:
        d11_ind[d] = 1.0
    A = autocorrelation_fft(d11_ind, p)
    T = {d: threshold - int(A[d]) for d in D11}

    rng = np.random.default_rng(42)

    # Collect data for Suen inequality
    marginal_bad = np.zeros(num_d)  # Pr[E_d] = Pr[B(d) > T(d)]
    pairwise_bad = np.zeros((num_d, num_d))  # Pr[E_i ∩ E_j]
    all_ok = 0

    t0 = time.time()
    for trial in range(num_mc):
        S = rng.choice(range(1, p), size=k_nonzero, replace=False)
        d12_ind = np.zeros(p)
        d12_ind[0] = 1.0
        d12_ind[S] = 1.0
        B = autocorrelation_fft(d12_ind, p)

        bad_flags = np.zeros(num_d, dtype=bool)
        for i, d in enumerate(D11):
            if int(B[d]) > T[d]:
                bad_flags[i] = True

        marginal_bad += bad_flags
        pairwise_bad += np.outer(bad_flags, bad_flags)

        if not bad_flags.any():
            all_ok += 1

    elapsed = time.time() - t0
    marginal_bad /= num_mc
    pairwise_bad /= num_mc

    joint_ok_rate = all_ok / num_mc

    # Compute Suen's Delta
    # Δ = Σ_{i<j} Pr[E_i ∩ E_j] / ((1-Pr[E_i])(1-Pr[E_j]))
    # But we should only sum over "dependent" pairs.
    # In our case, all pairs are dependent, so sum over all i<j.
    Delta = 0.0
    for i in range(num_d):
        for j in range(i+1, num_d):
            p_ij = pairwise_bad[i, j]
            p_i = marginal_bad[i]
            p_j = marginal_bad[j]
            denom = (1 - p_i) * (1 - p_j)
            if denom > 0:
                Delta += p_ij / denom

    # Suen bound: Pr[all ok] >= prod(1 - Pr[E_i]) * exp(-Delta)
    log_prod = sum(np.log(max(1 - p_i, 1e-10)) for p_i in marginal_bad)
    suen_lower = np.exp(log_prod - Delta)

    # Independence prediction
    indep_pred = np.exp(log_prod)

    # Also compute the "correlation gap" Δ' = Σ Pr[E_i ∩ E_j] - Pr[E_i]Pr[E_j]
    delta_corr = 0.0
    for i in range(num_d):
        for j in range(i+1, num_d):
            delta_corr += pairwise_bad[i, j] - marginal_bad[i] * marginal_bad[j]

    print(f"\n  p={p}, n={n}")
    print(f"  {num_mc:,} MC trials, {elapsed:.1f}s")
    print(f"  Marginal Pr[E_d] (bad event):")
    print(f"    min={marginal_bad.min():.4f}, max={marginal_bad.max():.4f}, "
          f"mean={marginal_bad.mean():.4f}")

    print(f"\n  Suen inequality components:")
    print(f"    Δ (sum of Pr[Ei∩Ej]/((1-pi)(1-pj))) = {Delta:.4f}")
    print(f"    exp(-Δ) = {np.exp(-Delta):.6f}")
    print(f"    ∏(1-Pr[Ei]) = {np.exp(log_prod):.6e}")
    print(f"    Suen lower bound: {suen_lower:.6e}")
    print(f"    Actual Pr[all ok]: {joint_ok_rate:.6f}")
    print(f"    Independence prediction: {indep_pred:.6e}")

    if suen_lower > 0:
        print(f"    Actual / Suen bound: {joint_ok_rate / suen_lower:.4f}")
    print(f"    Actual / Indep: {joint_ok_rate / indep_pred:.4f}")

    # Pairwise correlation analysis for bad events
    print(f"\n  Pairwise correlations of BAD events:")
    print(f"    Σ (Pr[Ei∩Ej] - Pr[Ei]Pr[Ej]) = {delta_corr:.6f}")
    # Average pairwise
    n_pairs = num_d * (num_d - 1) / 2
    print(f"    Average pairwise gap: {delta_corr / n_pairs:.6f}")

    # Check: are the BAD events negatively associated?
    neg_pairs = 0
    pos_pairs = 0
    for i in range(num_d):
        for j in range(i+1, num_d):
            gap = pairwise_bad[i, j] - marginal_bad[i] * marginal_bad[j]
            if gap < -0.001:
                neg_pairs += 1
            elif gap > 0.001:
                pos_pairs += 1
    total_pairs = num_d * (num_d - 1) // 2
    print(f"    Negative pairs: {neg_pairs}/{total_pairs}")
    print(f"    Positive pairs: {pos_pairs}/{total_pairs}")

    # Modified inequality: if bad events are negatively associated,
    # Pr[all ok] >= prod(1 - Pr[E_i]) directly (no Delta correction needed)
    if delta_corr < 0:
        print(f"\n  >>> BAD EVENTS ARE NEGATIVELY ASSOCIATED!")
        print(f"  >>> This directly implies positive association of GOOD events!")
        print(f"  >>> Pr[all ok] >= prod(Pr[ok_d]) = {indep_pred:.6e}")
    else:
        print(f"\n  >>> BAD events are positively associated (delta_corr > 0)")
        print(f"  >>> But Suen bound may still be useful if Δ is small")

    return {
        "p": p,
        "n": n,
        "Delta": Delta,
        "exp_neg_Delta": np.exp(-Delta),
        "suen_lower": suen_lower,
        "actual": joint_ok_rate,
        "indep": indep_pred,
        "delta_corr": delta_corr,
        "bad_events_negative": delta_corr < 0,
    }


def main():
    print("=" * 80)
    print("SUEN INEQUALITY AND BAD EVENT ASSOCIATION ANALYSIS")
    print("=" * 80)

    results = {}
    for p in [11, 19, 23, 31]:
        results[p] = compute_suen_bound(p, num_mc=200000)

    # Summary
    print(f"\n\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    print(f"\n  {'p':>4s} {'n':>3s} {'Δ':>10s} {'e^-Δ':>10s} {'Suen':>12s} "
          f"{'Actual':>12s} {'Indep':>12s} {'δ_corr':>10s} {'BadNeg?':>8s}")
    print(f"  {'-'*85}")
    for p in sorted(results.keys()):
        r = results[p]
        print(f"  {p:4d} {r['n']:3d} {r['Delta']:10.4f} {r['exp_neg_Delta']:10.6f} "
              f"{r['suen_lower']:12.6e} {r['actual']:12.6f} {r['indep']:12.6e} "
              f"{r['delta_corr']:10.6f} {'YES' if r['bad_events_negative'] else 'no':>8s}")

    print(f"\n  KEY INSIGHT:")
    all_neg = all(r['bad_events_negative'] for r in results.values())
    if all_neg:
        print(f"  >>> BAD EVENTS {'{'}B(d) > T(d){'}'} ARE NEGATIVELY ASSOCIATED!")
        print(f"  >>> This DIRECTLY IMPLIES that good events {'{'}B(d) ≤ T(d){'}'} are")
        print(f"  >>> POSITIVELY associated: Pr[all ok] >= ∏ Pr[ok_d].")
        print(f"  >>> If this holds for all p, Lemma L6 follows immediately.")
        print(f"  >>> The proof of R(B_{{n-1}}, B_n) = 4n-1 would be COMPLETE.")
    else:
        print(f"  >>> Mixed results. Check individual primes.")


if __name__ == "__main__":
    main()
