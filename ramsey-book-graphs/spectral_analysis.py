#!/usr/bin/env python3
"""
Full spectral analysis of known valid constructions at p=43, 47, 59.

Task 1: Full power spectrum |hat(1_{D11})(k)|^2 for all k=1..p-1
Task 2: Spectral complementarity between D11 and D12 indicator functions

For each valid (D11, D12) pair, compute:
- Power spectrum of D11 and D12
- Cross-spectrum analysis
- Spectral complementarity measure
"""

import json
import cmath
import math
import numpy as np


def power_spectrum(indicator_set, p):
    """Compute |hat(1_S)(k)|^2 for k=0,...,p-1."""
    spectrum = []
    for k in range(p):
        omega = cmath.exp(2j * cmath.pi * k / p)
        s = sum(omega ** d for d in indicator_set)
        spectrum.append(abs(s) ** 2)
    return spectrum


def cross_spectrum(set1, set2, p):
    """Compute hat(1_{S1})(k) * conj(hat(1_{S2})(k)) for k=0,...,p-1."""
    result = []
    for k in range(p):
        omega = cmath.exp(2j * cmath.pi * k / p)
        s1 = sum(omega ** d for d in set1)
        s2 = sum(omega ** d for d in set2)
        result.append(s1 * s2.conjugate())
    return result


def analyze_solution(sol, p):
    """Full spectral analysis of a (D11, D12) pair."""
    d11 = set(sol['D11'])
    d12 = set(sol['D12'])
    d22 = set(sol.get('D22', set(range(1, p)) - d11))

    n = (p + 1) // 2  # book parameter

    # Power spectra
    spec_d11 = power_spectrum(d11, p)
    spec_d12 = power_spectrum(d12, p)

    # Non-trivial spectra (k != 0)
    nt_d11 = spec_d11[1:]
    nt_d12 = spec_d12[1:]

    # Stats for D11
    d11_mean = sum(nt_d11) / len(nt_d11)
    d11_max = max(nt_d11)
    d11_min = min(nt_d11)
    d11_flat = d11_max / d11_mean if d11_mean > 0 else 0
    d11_var = sum((x - d11_mean)**2 for x in nt_d11) / len(nt_d11)
    d11_std = math.sqrt(d11_var)

    # Stats for D12
    d12_mean = sum(nt_d12) / len(nt_d12)
    d12_max = max(nt_d12)
    d12_min = min(nt_d12)
    d12_flat = d12_max / d12_mean if d12_mean > 0 else 0
    d12_var = sum((x - d12_mean)**2 for x in nt_d12) / len(nt_d12)
    d12_std = math.sqrt(d12_var)

    # Cross-spectrum analysis
    cross = cross_spectrum(d11, d12, p)
    cross_nt = cross[1:]  # non-trivial
    cross_mags = [abs(c) for c in cross_nt]
    cross_reals = [c.real for c in cross_nt]

    # Spectral complementarity: correlation between D11 and D12 power spectra
    # If they're complementary, when D11 power is high, D12 power should be low
    nt_d11_arr = np.array(nt_d11)
    nt_d12_arr = np.array(nt_d12)

    if np.std(nt_d11_arr) > 0 and np.std(nt_d12_arr) > 0:
        spectral_corr = np.corrcoef(nt_d11_arr, nt_d12_arr)[0, 1]
    else:
        spectral_corr = 0.0

    # Complementarity score: sum of |hat_D11|^2 * |hat_D12|^2 (low = complementary)
    overlap_score = sum(a * b for a, b in zip(nt_d11, nt_d12)) / (sum(nt_d11) * sum(nt_d12)) * len(nt_d11)
    # If independent/uncorrelated: overlap_score ≈ 1
    # If complementary: overlap_score < 1
    # If aligned: overlap_score > 1

    # Product spectrum: |hat_D11|^2 * |hat_D12|^2 controls B-value distribution
    product_spec = [a * b for a, b in zip(nt_d11, nt_d12)]
    product_mean = sum(product_spec) / len(product_spec)
    product_max = max(product_spec)
    product_flat = product_max / product_mean if product_mean > 0 else 0

    # Sorted spectra for detailed view
    d11_sorted = sorted(enumerate(nt_d11), key=lambda x: -x[1])
    d12_sorted = sorted(enumerate(nt_d12), key=lambda x: -x[1])

    # Top 5 spectral peaks for D11 and their D12 values
    top5_d11 = [(k+1, nt_d11[k], nt_d12[k]) for k, _ in d11_sorted[:5]]
    top5_d12 = [(k+1, nt_d12[k], nt_d11[k]) for k, _ in d12_sorted[:5]]

    return {
        'p': p,
        'n': n,
        'd11_size': len(d11),
        'd12_size': len(d12),
        'd11_spectrum': {
            'mean': round(d11_mean, 4),
            'max': round(d11_max, 4),
            'min': round(d11_min, 4),
            'flatness': round(d11_flat, 4),
            'std': round(d11_std, 4),
            'cv': round(d11_std / d11_mean, 4) if d11_mean > 0 else 0,
            'max_over_p': round(d11_max / p, 4),
        },
        'd12_spectrum': {
            'mean': round(d12_mean, 4),
            'max': round(d12_max, 4),
            'min': round(d12_min, 4),
            'flatness': round(d12_flat, 4),
            'std': round(d12_std, 4),
            'cv': round(d12_std / d12_mean, 4) if d12_mean > 0 else 0,
            'max_over_p': round(d12_max / p, 4),
        },
        'spectral_complementarity': {
            'correlation': round(float(spectral_corr), 4),
            'overlap_score': round(overlap_score, 4),
            'product_flatness': round(product_flat, 4),
        },
        'top5_d11_peaks': [(k, round(v11, 2), round(v12, 2))
                           for k, v11, v12 in top5_d11],
        'top5_d12_peaks': [(k, round(v12, 2), round(v11, 2))
                           for k, v12, v11 in top5_d12],
        'full_spectrum_d11': [round(x, 4) for x in nt_d11],
        'full_spectrum_d12': [round(x, 4) for x in nt_d12],
    }


def main():
    with open('solutions_registry.json') as f:
        reg = json.load(f)

    results = {}

    for sol in sorted(reg['solutions'], key=lambda s: s['m']):
        m = sol['m']
        is_prime = all(m % i != 0 for i in range(2, int(m**0.5)+1)) and m > 1
        if not is_prime or m % 4 != 3:
            continue

        p = m
        print(f"\n{'='*60}")
        print(f"p = {p} (n = {sol['n']})")
        print(f"{'='*60}")

        result = analyze_solution(sol, p)
        results[str(p)] = result

        # Print summary
        print(f"  |D11| = {result['d11_size']}, |D12| = {result['d12_size']}")
        print(f"\n  D11 power spectrum:")
        print(f"    mean = {result['d11_spectrum']['mean']:.4f}")
        print(f"    max  = {result['d11_spectrum']['max']:.4f} (max/p = {result['d11_spectrum']['max_over_p']:.4f})")
        print(f"    min  = {result['d11_spectrum']['min']:.4f}")
        print(f"    flatness = {result['d11_spectrum']['flatness']:.4f}")
        print(f"    CV = {result['d11_spectrum']['cv']:.4f}")

        print(f"\n  D12 power spectrum:")
        print(f"    mean = {result['d12_spectrum']['mean']:.4f}")
        print(f"    max  = {result['d12_spectrum']['max']:.4f} (max/p = {result['d12_spectrum']['max_over_p']:.4f})")
        print(f"    min  = {result['d12_spectrum']['min']:.4f}")
        print(f"    flatness = {result['d12_spectrum']['flatness']:.4f}")
        print(f"    CV = {result['d12_spectrum']['cv']:.4f}")

        print(f"\n  Spectral complementarity:")
        print(f"    Correlation(|hat_D11|^2, |hat_D12|^2) = {result['spectral_complementarity']['correlation']:.4f}")
        print(f"    Overlap score = {result['spectral_complementarity']['overlap_score']:.4f} (1.0 = independent)")
        print(f"    Product flatness = {result['spectral_complementarity']['product_flatness']:.4f}")

        print(f"\n  Top 5 D11 spectral peaks (freq, |hat_D11|^2, |hat_D12|^2):")
        for k, v11, v12 in result['top5_d11_peaks']:
            ratio = v12/v11 if v11 > 0 else 0
            print(f"    k={k:>3}: D11={v11:>8.2f}, D12={v12:>8.2f}, ratio={ratio:.3f}")

        print(f"\n  Top 5 D12 spectral peaks (freq, |hat_D12|^2, |hat_D11|^2):")
        for k, v12, v11 in result['top5_d12_peaks']:
            ratio = v11/v12 if v12 > 0 else 0
            print(f"    k={k:>3}: D12={v12:>8.2f}, D11={v11:>8.2f}, ratio={ratio:.3f}")

    # Cross-prime comparison
    print(f"\n{'='*60}")
    print("CROSS-PRIME SUMMARY")
    print(f"{'='*60}")
    print(f"\n{'p':>4} {'D11_flat':>9} {'D12_flat':>9} {'corr':>7} {'overlap':>8} {'D11max/p':>9} {'D12max/p':>9}")
    print('-' * 60)
    for pkey in sorted(results.keys(), key=int):
        r = results[pkey]
        print(f"{pkey:>4} {r['d11_spectrum']['flatness']:>9.4f} {r['d12_spectrum']['flatness']:>9.4f} "
              f"{r['spectral_complementarity']['correlation']:>7.4f} "
              f"{r['spectral_complementarity']['overlap_score']:>8.4f} "
              f"{r['d11_spectrum']['max_over_p']:>9.4f} {r['d12_spectrum']['max_over_p']:>9.4f}")

    with open('spectral_analysis_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to spectral_analysis_results.json")


if __name__ == '__main__':
    main()
