"""
Phase 4: Algebraic analysis and construction for primes p ≡ 3 mod 4.

Key findings from Phases 1-3:
- Every symmetric D11 is automatically QR-balanced (since -1 ∈ QNR for p ≡ 3 mod 4)
- The D22 "loose" constraints are actually binding for p ≥ 31
- All margins in known solutions are exactly 0 (extremally tight)
- Constraints couple D11 and D22: pushing A(d)+B(d) down at D11 pushes it up at D22

This script performs deep algebraic analysis of ALL known solutions to find
a constructive pattern for D12 given D11, and attempts to prove existence.
"""

import sys
import os
import json
import math
import numpy as np
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import BlockCirculantGraph, verify_construction, Delta


# Collect all known solutions for primes p ≡ 3 mod 4
KNOWN = {
    11: {"D11": {1, 2, 4, 7, 9, 10}, "D12": {0, 4, 5, 7, 10}},
    19: {"D11": {1, 2, 3, 6, 8, 11, 13, 16, 17, 18},
         "D12": {0, 2, 6, 8, 9, 12, 13, 17, 18}},
    23: {"D11": {5, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17, 18},
         "D12": {0, 1, 2, 6, 10, 13, 14, 16, 18, 20, 21}},
    31: {"D11": {6, 7, 8, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 23, 24, 25},
         "D12": {0, 1, 2, 3, 8, 11, 12, 13, 15, 18, 20, 21, 24, 27, 29}},
    43: {"D11": {1, 2, 5, 10, 11, 13, 16, 17, 18, 19, 20, 23, 24, 25, 26, 27,
                 30, 32, 33, 38, 41, 42},
         "D12": {0, 2, 5, 6, 8, 11, 15, 16, 20, 24, 25, 27, 28, 31, 32, 34,
                 35, 36, 37, 39, 41}},
    47: {"D11": {2, 6, 8, 10, 11, 13, 14, 16, 19, 20, 21, 26, 27, 28, 31, 33,
                 34, 36, 37, 39, 41, 45},
         "D12": {0, 1, 8, 9, 10, 12, 19, 23, 24, 25, 26, 28, 29, 32, 33, 35,
                 38, 39, 40, 42, 44, 45, 46}},
    59: {"D11": {1, 2, 3, 5, 6, 8, 9, 10, 11, 17, 19, 25, 27, 29, 30, 32, 34,
                 40, 42, 48, 49, 50, 51, 53, 54, 56, 57, 58},
         "D12": {0, 1, 4, 6, 7, 9, 10, 11, 13, 18, 20, 21, 23, 24, 25, 26, 32,
                 33, 36, 38, 39, 42, 43, 45, 49, 50, 51, 54, 55}},
}


def primitive_root(p):
    """Find smallest primitive root mod p."""
    for g in range(2, p):
        seen = set()
        x = 1
        for _ in range(p - 1):
            seen.add(x)
            x = (x * g) % p
        if len(seen) == p - 1:
            return g
    return None


def discrete_log(x, g, p):
    """Compute discrete log of x base g mod p. Returns i such that g^i = x mod p."""
    val = 1
    for i in range(p - 1):
        if val == x:
            return i
        val = (val * g) % p
    return None


def autocorrelation_fft(indicator, m):
    """Compute Delta(S,S,d) for all d via FFT."""
    fft_val = np.fft.fft(indicator)
    autocorr = np.fft.ifft(np.abs(fft_val) ** 2).real
    return np.round(autocorr).astype(int)


def power_spectrum_fft(indicator, m):
    """Compute |hat{S}(k)|^2 for all k via FFT."""
    fft_val = np.fft.fft(indicator)
    return np.abs(fft_val) ** 2


def check_all_constraints(p, D11, D12):
    """Check all 4 constraints. Returns (valid, max_violation, details)."""
    D11_set = set(D11)
    D22_set = set(range(1, p)) - D11_set
    n = (p + 1) // 2

    ind11 = np.zeros(p, dtype=np.float64)
    for j in D11:
        ind11[j] = 1.0
    ind12 = np.zeros(p, dtype=np.float64)
    for j in D12:
        ind12[j] = 1.0

    A = autocorrelation_fft(ind11, p)
    B = autocorrelation_fft(ind12, p)

    max_viol = 0
    details = {"C1": [], "C2": [], "C3": [], "C4": []}

    # C1: V1V1 red, d ∈ D11: A(d) + B(d) ≤ n-2 = (p-3)/2
    thresh_c1 = (p - 3) // 2
    for d in sorted(D11_set):
        val = int(A[d]) + int(B[d])
        viol = val - thresh_c1
        if viol > max_viol:
            max_viol = viol
        details["C1"].append((d, val, thresh_c1, viol))

    # C2: V2V2 blue, d ∈ D11: A(d) + B(p-d) ≤ (p-3)/2
    # From V2V2 theorem: V2V2(d) = A(d) + B(-d) + (m-2-2k) + 2*[d∈D11]
    # For blue (d ∈ D11): V2V2(d) ≤ n-1
    # V2V2(d) = A(d) + B(p-d) - 3 + 2 = A(d) + B(p-d) - 1 (with k=(p+1)/2)
    # A(d) + B(p-d) - 1 ≤ (p-1)/2 → A(d) + B(p-d) ≤ (p+1)/2
    thresh_c2 = (p + 1) // 2
    for d in sorted(D11_set):
        val = int(A[d]) + int(B[(p - d) % p])
        viol = val - thresh_c2
        if viol > max_viol:
            max_viol = viol
        details["C2"].append((d, val, thresh_c2, viol))

    # C3: V1V1 blue, d ∈ D22: A(d) + B(d) ≤ n-1 = (p-1)/2
    thresh_c3 = (p - 1) // 2
    for d in sorted(D22_set):
        val = int(A[d]) + int(B[d])
        viol = val - thresh_c3
        if viol > max_viol:
            max_viol = viol
        details["C3"].append((d, val, thresh_c3, viol))

    # C4: V2V2 red, d ∈ D22: A(d) + B(p-d) ≤ (p+3)/2
    # V2V2(d) = A(d) + B(p-d) - 3 (with 2*[d∈D11]=0)
    # Red: V2V2(d) ≤ n-2 = (p-3)/2
    # A(d) + B(p-d) - 3 ≤ (p-3)/2 → A(d) + B(p-d) ≤ (p+3)/2
    thresh_c4 = (p + 3) // 2
    for d in sorted(D22_set):
        val = int(A[d]) + int(B[(p - d) % p])
        viol = val - thresh_c4
        if viol > max_viol:
            max_viol = viol
        details["C4"].append((d, val, thresh_c4, viol))

    return max_viol <= 0, max_viol, details


def analyze_spectral(p, D11, D12):
    """Spectral analysis of (D11, D12) pair."""
    ind11 = np.zeros(p, dtype=np.float64)
    for j in D11:
        ind11[j] = 1.0
    ind12 = np.zeros(p, dtype=np.float64)
    for j in D12:
        ind12[j] = 1.0

    spec11 = power_spectrum_fft(ind11, p)  # |hat{D11}(k)|^2
    spec12 = power_spectrum_fft(ind12, p)  # |hat{D12}(k)|^2
    P = spec11 + spec12  # power spectrum

    fft11 = np.fft.fft(ind11)
    fft12 = np.fft.fft(ind12)

    return {
        "spec11": spec11,
        "spec12": spec12,
        "P": P,
        "fft11": fft11,
        "fft12": fft12,
        "P_mean": np.mean(P[1:]),
        "P_std": np.std(P[1:]),
        "P_max": np.max(P[1:]),
        "P_min": np.min(P[1:]),
        "rms_from_p": np.sqrt(np.mean((P[1:] - p) ** 2)),
    }


def analyze_index_structure(p, S, g):
    """Analyze set S in terms of discrete log structure."""
    indices = sorted(discrete_log(x, g, p) for x in S if x != 0)
    # Mod 2 (QR/QNR)
    qr_count = sum(1 for i in indices if i % 2 == 0)
    qnr_count = len(indices) - qr_count
    # Mod 3
    mod3 = [0, 0, 0]
    for i in indices:
        mod3[i % 3] += 1
    # Mod 4
    mod4 = [0, 0, 0, 0]
    for i in indices:
        mod4[i % 4] += 1
    # Gaps between consecutive indices
    gaps = [indices[i + 1] - indices[i] for i in range(len(indices) - 1)]
    if len(indices) > 1:
        gaps.append(p - 1 - indices[-1] + indices[0])

    return {
        "indices": indices,
        "qr": qr_count,
        "qnr": qnr_count,
        "mod3": mod3,
        "mod4": mod4,
        "gap_min": min(gaps) if gaps else 0,
        "gap_max": max(gaps) if gaps else 0,
        "gap_mean": sum(gaps) / len(gaps) if gaps else 0,
    }


def analyze_d12_as_translate(p, D11, D12, g):
    """Check if D12 \ {0} is a multiplicative translate of D11."""
    D12_nonzero = D12 - {0}
    best_shift = None
    best_overlap = 0
    for a in range(1, p - 1):
        multiplier = pow(g, a, p)
        shifted = frozenset((multiplier * d) % p for d in D11)
        overlap = len(shifted & D12_nonzero)
        if overlap > best_overlap:
            best_overlap = overlap
            best_shift = a
    return best_shift, best_overlap, len(D12_nonzero)


def analyze_d12_as_union(p, D11, D12, g):
    """Check if D12 \ {0} is a union of D11-halves (QR part and QNR part) with shifts."""
    D12_nonzero = D12 - {0}

    # Split D11 into QR and QNR parts
    QR = set()
    for x in range(1, p):
        QR.add(pow(x, 2, p))
    d11_qr = D11 & QR
    d11_qnr = D11 - QR

    best = {"desc": "none", "overlap": 0}

    # Try: D12 = {0} ∪ (g^a · d11_qr) ∪ (g^b · d11_qnr) for various a, b
    for a in range(p - 1):
        ma = pow(g, a, p)
        shifted_qr = frozenset((ma * d) % p for d in d11_qr)
        overlap_qr = len(shifted_qr & D12_nonzero)
        for b in range(p - 1):
            mb = pow(g, b, p)
            shifted_qnr = frozenset((mb * d) % p for d in d11_qnr)
            combined = shifted_qr | shifted_qnr
            overlap = len(combined & D12_nonzero)
            if overlap > best["overlap"]:
                best = {
                    "desc": f"g^{a}·QR(D11) ∪ g^{b}·QNR(D11)",
                    "overlap": overlap,
                    "a": a,
                    "b": b,
                    "size_combined": len(combined),
                }

    return best


def try_legendre_corrector(p, D11):
    """Try to construct D12 as a 'Legendre corrector' - minimizing max P(k) deviation from p."""
    ind11 = np.zeros(p, dtype=np.float64)
    for j in D11:
        ind11[j] = 1.0
    spec11 = power_spectrum_fft(ind11, p)

    # Ideal: |hat{D12}(k)|^2 = p - |hat{D11}(k)|^2 for k > 0
    ideal_spec12 = np.zeros(p, dtype=np.float64)
    ideal_spec12[0] = ((p - 1) / 2) ** 2  # |D12|^2
    for k in range(1, p):
        ideal_spec12[k] = max(0, p - spec11[k])

    # Check feasibility
    infeasible = sum(1 for k in range(1, p) if p - spec11[k] < 0)

    return {
        "ideal_spec12": ideal_spec12,
        "infeasible_count": infeasible,
        "ideal_total": np.sum(ideal_spec12[1:]),
        "needed_total": ((p - 1) / 2) * ((p - 1) / 2 - 1),  # |D12|(|D12|-1)
    }


def analyze_constraint_margins(p, D11, D12):
    """Detailed margin analysis across all 4 constraints."""
    valid, max_viol, details = check_all_constraints(p, D11, D12)

    margins = {"C1": [], "C2": [], "C3": [], "C4": []}
    for c in ["C1", "C2", "C3", "C4"]:
        for (d, val, thresh, viol) in details[c]:
            margins[c].append(thresh - val)

    summary = {}
    for c in ["C1", "C2", "C3", "C4"]:
        if margins[c]:
            summary[c] = {
                "min_margin": min(margins[c]),
                "max_margin": max(margins[c]),
                "mean_margin": sum(margins[c]) / len(margins[c]),
                "num_tight": sum(1 for m in margins[c] if m == 0),
                "num_violated": sum(1 for m in margins[c] if m < 0),
            }
    return summary, valid


def try_algebraic_constructions(p, D11, g):
    """Try several algebraic constructions for D12 and check validity."""
    results = []
    d12_size = (p - 1) // 2

    QR = set()
    for x in range(1, p):
        QR.add(pow(x, 2, p))
    QNR = set(range(1, p)) - QR

    # Construction 1: D12 = {0} ∪ QR
    D12_try = {0} | QR
    if len(D12_try) == d12_size:
        valid, max_viol, _ = check_all_constraints(p, D11, D12_try)
        results.append(("QR ∪ {0}", valid, max_viol))

    # Construction 2: D12 = {0} ∪ QNR
    D12_try = {0} | QNR
    if len(D12_try) == d12_size:
        valid, max_viol, _ = check_all_constraints(p, D11, D12_try)
        results.append(("QNR ∪ {0}", valid, max_viol))

    # Construction 3: D12 = {0} ∪ (g·D11_half)
    # Take the "first half" of D11 (indices 0..len/2-1 in sorted order)
    sorted_d11 = sorted(D11)
    half = len(sorted_d11) // 2
    d11_first = set(sorted_d11[:half])
    d11_second = set(sorted_d11[half:])

    for shift in range(p - 1):
        mult = pow(g, shift, p)
        D12_try = {0} | frozenset((mult * d) % p for d in d11_first)
        if len(D12_try) == d12_size:
            valid, max_viol, _ = check_all_constraints(p, D11, D12_try)
            if valid:
                results.append((f"{{0}} ∪ g^{shift}·D11_first_half", valid, max_viol))
                break
        D12_try = {0} | frozenset((mult * d) % p for d in d11_second)
        if len(D12_try) == d12_size:
            valid, max_viol, _ = check_all_constraints(p, D11, D12_try)
            if valid:
                results.append((f"{{0}} ∪ g^{shift}·D11_second_half", valid, max_viol))
                break

    # Construction 4: D12 = {0} ∪ {d : index(d) ∈ S} for various index sets S
    for r in [2, 3, 4, 6]:
        if (p - 1) % r != 0:
            continue
        for offset in range(r):
            idx_set = set(i for i in range(p - 1) if i % r == offset)
            D12_elements = {pow(g, i, p) for i in idx_set}
            D12_try = {0} | D12_elements
            if len(D12_try) == d12_size:
                valid, max_viol, _ = check_all_constraints(p, D11, D12_try)
                results.append((f"index ≡ {offset} mod {r} ∪ {{0}}", valid, max_viol))

    # Construction 5: D12 = {0} ∪ (D11 - best_element + worst_from_D22)
    # Greedy: start from D11 and swap elements to satisfy D12 constraints
    # This is more of a search than algebraic, skip for now

    return results


def analyze_b_symmetry(p, D12):
    """Analyze B(d) vs B(p-d) to understand asymmetry of D12."""
    ind12 = np.zeros(p, dtype=np.float64)
    for j in D12:
        ind12[j] = 1.0
    B = autocorrelation_fft(ind12, p)

    asym = []
    for d in range(1, (p + 1) // 2):
        asym.append((d, int(B[d]), int(B[p - d]), int(B[d]) - int(B[p - d])))
    return asym


def main():
    primes = sorted(KNOWN.keys())

    print("=" * 100)
    print("PHASE 4: ALGEBRAIC ANALYSIS OF p ≡ 3 mod 4 CONSTRUCTIONS")
    print("=" * 100)

    # Key mathematical observation
    print("""
KEY OBSERVATION: For p ≡ 3 mod 4, -1 is a quadratic non-residue.
Therefore, for ANY symmetric D11: each pair {d, p-d} contains one QR and one QNR.
This means ALL symmetric D11 of size (p+1)/2 are automatically QR-balanced.
The QR-balance seen in Phase 1 is not a distinguishing property — it's universal.
""")

    all_results = {}

    for p in primes:
        n = (p + 1) // 2
        g = primitive_root(p)
        D11 = KNOWN[p]["D11"]
        D12 = KNOWN[p]["D12"]

        print(f"\n{'=' * 100}")
        print(f"PRIME p = {p}, n = {n}, primitive root g = {g}")
        print(f"|D11| = {len(D11)}, |D12| = {len(D12)}, expected |D12| = {(p-1)//2}")
        print(f"{'=' * 100}")

        # 1. Index structure analysis
        print(f"\n  1. INDEX STRUCTURE (discrete log base {g}):")
        d11_idx = analyze_index_structure(p, D11, g)
        d12_idx = analyze_index_structure(p, D12 - {0}, g)

        print(f"    D11 indices: {d11_idx['indices']}")
        print(f"    D11 QR/QNR: {d11_idx['qr']}/{d11_idx['qnr']}")
        print(f"    D11 mod 3: {d11_idx['mod3']}")
        print(f"    D11 mod 4: {d11_idx['mod4']}")
        print(f"    D11 gap stats: min={d11_idx['gap_min']}, max={d11_idx['gap_max']}, "
              f"mean={d11_idx['gap_mean']:.1f}")

        print(f"\n    D12\\{{0}} indices: {d12_idx['indices']}")
        print(f"    D12\\{{0}} QR/QNR: {d12_idx['qr']}/{d12_idx['qnr']}")
        print(f"    D12\\{{0}} mod 3: {d12_idx['mod3']}")
        print(f"    D12\\{{0}} mod 4: {d12_idx['mod4']}")
        print(f"    D12\\{{0}} gap stats: min={d12_idx['gap_min']}, max={d12_idx['gap_max']}, "
              f"mean={d12_idx['gap_mean']:.1f}")

        # 2. Spectral analysis
        print(f"\n  2. SPECTRAL ANALYSIS:")
        spec = analyze_spectral(p, D11, D12)

        print(f"    P(k) = |D̂11(k)|² + |D̂12(k)|² for k > 0:")
        print(f"      Mean: {spec['P_mean']:.2f} (ideal = {p})")
        print(f"      Std:  {spec['P_std']:.2f}")
        print(f"      Range: [{spec['P_min']:.2f}, {spec['P_max']:.2f}]")
        print(f"      RMS from p: {spec['rms_from_p']:.2f}")

        # Detailed P(k) values
        P_vals = spec["P"]
        if p <= 31:
            print(f"      P(k) values: {[f'{v:.1f}' for v in P_vals[1:]]}")

        # Phase analysis: angle between D̂11(k) and D̂12(k)
        fft11 = spec["fft11"]
        fft12 = spec["fft12"]
        phases = []
        for k in range(1, p):
            if abs(fft11[k]) > 0.01 and abs(fft12[k]) > 0.01:
                angle = np.angle(fft11[k] * np.conj(fft12[k]))
                phases.append(angle)
        if phases:
            print(f"    Phase angles between D̂11 and D̂12:")
            print(f"      Mean: {np.mean(phases):.4f}")
            print(f"      Std:  {np.std(phases):.4f}")

        # 3. Multiplicative translate check
        print(f"\n  3. D12 AS MULTIPLICATIVE TRANSLATE OF D11:")
        shift, overlap, total = analyze_d12_as_translate(p, D11, D12, g)
        print(f"    Best shift: g^{shift}, overlap: {overlap}/{total} "
              f"({overlap/total:.1%})")

        # 4. Union of cosets check
        if p <= 43:
            print(f"\n  4. D12 AS SHIFTED QR/QNR PARTS OF D11:")
            union_result = analyze_d12_as_union(p, D11, D12, g)
            print(f"    Best: {union_result['desc']}")
            print(f"    Overlap: {union_result['overlap']}/{len(D12)-1}")

        # 5. Constraint margins
        print(f"\n  5. CONSTRAINT MARGIN ANALYSIS:")
        margin_summary, valid = analyze_constraint_margins(p, D11, D12)
        for c in ["C1", "C2", "C3", "C4"]:
            if c in margin_summary:
                s = margin_summary[c]
                name = {"C1": "V1V1 red  (d∈D11)",
                        "C2": "V2V2 blue (d∈D11)",
                        "C3": "V1V1 blue (d∈D22)",
                        "C4": "V2V2 red  (d∈D22)"}[c]
                print(f"    {name}: min_margin={s['min_margin']}, "
                      f"mean={s['mean_margin']:.2f}, "
                      f"tight={s['num_tight']}, violated={s['num_violated']}")

        # 6. B(d) asymmetry
        print(f"\n  6. B(d) ASYMMETRY [B(d) - B(p-d)]:")
        asym = analyze_b_symmetry(p, D12)
        diffs = [a[3] for a in asym]
        print(f"    Range: [{min(diffs)}, {max(diffs)}]")
        print(f"    Mean: {sum(diffs)/len(diffs):.2f}")
        print(f"    Num zero: {sum(1 for d in diffs if d == 0)}/{len(diffs)}")
        if p <= 23:
            print(f"    Values: {diffs}")

        # 7. Legendre corrector analysis
        print(f"\n  7. LEGENDRE CORRECTOR FEASIBILITY:")
        leg = try_legendre_corrector(p, D11)
        print(f"    Infeasible ideal frequencies: {leg['infeasible_count']}/{p-1}")
        print(f"    Ideal D12 energy (sum |D̂12(k)|² for k>0): {leg['ideal_total']:.1f}")
        print(f"    Needed energy (|D12|(|D12|-1)): {leg['needed_total']:.1f}")

        # 8. Try algebraic constructions
        print(f"\n  8. ALGEBRAIC CONSTRUCTION ATTEMPTS:")
        constructions = try_algebraic_constructions(p, D11, g)
        for desc, ok, mv in constructions:
            status = "✓ VALID" if ok else f"✗ max_viol={mv}"
            print(f"    {desc}: {status}")

        all_results[p] = {
            "p": p, "g": g,
            "d11_idx": d11_idx,
            "d12_idx": d12_idx,
            "P_mean": spec["P_mean"],
            "P_std": spec["P_std"],
            "rms_from_p": spec["rms_from_p"],
            "margin_summary": margin_summary,
        }

    # Global summary
    print(f"\n{'=' * 100}")
    print("GLOBAL PATTERN ANALYSIS")
    print(f"{'=' * 100}")

    print("\n  Power spectrum P(k) quality across primes:")
    print(f"  {'p':>5} {'P_mean':>8} {'P_std':>8} {'RMS':>8} {'P/p':>8} {'std/√p':>8}")
    for p in primes:
        r = all_results[p]
        print(f"  {p:5d} {r['P_mean']:8.2f} {r['P_std']:8.2f} "
              f"{r['rms_from_p']:8.2f} {r['P_mean']/p:8.4f} "
              f"{r['P_std']/math.sqrt(p):8.4f}")

    print(f"\n  Tightness analysis:")
    print(f"  {'p':>5} {'C1_tight':>10} {'C2_tight':>10} {'C3_tight':>10} {'C4_tight':>10}")
    for p in primes:
        ms = all_results[p]["margin_summary"]
        row = f"  {p:5d}"
        for c in ["C1", "C2", "C3", "C4"]:
            if c in ms:
                row += f" {ms[c]['num_tight']:10d}"
            else:
                row += f" {'N/A':>10}"
        print(row)

    # Save results
    output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "algebraic_results.json")
    save_data = {}
    for p in primes:
        r = all_results[p]
        save_data[str(p)] = {
            "p": r["p"], "g": r["g"],
            "P_mean": float(r["P_mean"]),
            "P_std": float(r["P_std"]),
            "rms_from_p": float(r["rms_from_p"]),
        }
    with open(output_path, "w") as f:
        json.dump(save_data, f, indent=2)
    print(f"\n  Results saved to {output_path}")

    # Now attempt construction for ALL primes up to ~100
    print(f"\n{'=' * 100}")
    print("CONSTRUCTION SEARCH: Exhaustive D12 search for small p")
    print(f"{'=' * 100}")

    # For p=7, try all D12
    p = 7
    g = primitive_root(p)
    d12_size = (p - 1) // 2  # 3
    print(f"\n  p={p}: Trying all symmetric D11 with all D12 of size {d12_size}")

    from itertools import combinations
    pairs_7 = [(d, (p - d) % p) for d in range(1, (p + 1) // 2)]
    num_pairs = len(pairs_7)
    d11_pairs_needed = (p + 1) // 4  # 2

    found_7 = 0
    for d11_combo in combinations(range(num_pairs), d11_pairs_needed):
        D11 = set()
        for i in d11_combo:
            D11.add(pairs_7[i][0])
            D11.add(pairs_7[i][1])

        # Try all D12 of correct size containing 0
        others = list(range(1, p))
        for d12_combo in combinations(others, d12_size - 1):
            D12 = {0} | set(d12_combo)
            valid, max_viol, _ = check_all_constraints(p, D11, D12)
            if valid:
                found_7 += 1
                if found_7 <= 3:
                    print(f"    VALID: D11={sorted(D11)}, D12={sorted(D12)}")

    print(f"    Total valid pairs: {found_7}")
    if found_7 == 0:
        # Try without 0 ∈ D12 constraint and with different D11/D12 sizes
        print(f"    No solutions with standard sizes. Trying alternative conventions...")

        # p=7, n=4: m=7, N=14. Need d1+d2=m-1=6.
        # Standard: d1=(p+1)/2=4, d2=(p-1)/2=3
        # Alternative: d1=(p-1)/2=3, d2=(p+1)/2=4 (swap V1/V2 roles)
        d11_size_alt = (p - 1) // 2  # 3
        d12_size_alt = (p + 1) // 2  # 4

        found_7_alt = 0
        for d11_combo in combinations(range(1, p), d11_size_alt):
            D11 = set(d11_combo)
            if not all((p - d) % p in D11 for d in D11):
                continue  # must be symmetric

            for d12_combo in combinations(range(1, p), d12_size_alt - 1):
                D12 = {0} | set(d12_combo)
                valid, max_viol, _ = check_all_constraints(p, D11, D12)
                if valid:
                    found_7_alt += 1
                    if found_7_alt <= 3:
                        print(f"    ALT VALID: D11={sorted(D11)}, D12={sorted(D12)}, "
                              f"|D11|={len(D11)}, |D12|={len(D12)}")
        print(f"    Alternative convention valid pairs: {found_7_alt}")


if __name__ == "__main__":
    main()
