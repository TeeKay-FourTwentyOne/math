#!/usr/bin/env python3
"""
Character-sum approach to the second moment for L6.

Key idea: Express N(D11) in terms of characters of Z_p, then compute
E_{D11}[N^2] via character orthogonality.

For a symmetric D11 of size n = (p+1)/2:
- D11 is determined by choosing which (p-1)/2 pairs {d, p-d} to include
- D12 = {0} union S, S a k-subset of {1,...,p-1}
- N(D11) = #{S : B(d) <= tau(d) for all d}

The constraint B(d) <= tau(d) depends on D11 through:
(a) tau(d) depends on A(d) and whether d is in D11
(b) A(d) depends on D11

Instead of the Gaussian proxy, try a DIRECT approach:

For the Paley graph P_p (p = 3 mod 4), the adjacency matrix has eigenvalues:
  lambda_0 = (p-1)/2
  lambda_chi = (-1 + eta(x) * sqrt(p)) / 2 for nontrivial character chi

where eta(x) = Legendre symbol.

Key question: Can we express N(D11) in terms of character sums and
use Weil-type bounds?

This script explores the structure computationally.
"""

import numpy as np
from itertools import combinations
from math import comb, sqrt
from collections import defaultdict


def legendre(a, p):
    """Legendre symbol (a/p)."""
    if a % p == 0:
        return 0
    return pow(a, (p - 1) // 2, p) * (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def legendre_correct(a, p):
    """Legendre symbol using Euler criterion."""
    if a % p == 0:
        return 0
    val = pow(a % p, (p - 1) // 2, p)
    if val == 1:
        return 1
    elif val == p - 1:
        return -1
    else:
        return 0


def get_qr(p):
    """Get quadratic residues mod p."""
    return set(pow(x, 2, p) for x in range(1, p))


def compute_A(D11_set, p):
    """Compute A(d) for all d = 1,...,p-1."""
    A = [0] * p
    for a in D11_set:
        for b in D11_set:
            if a != b:
                A[(a - b) % p] += 1
    return A


def compute_B(S_set, p):
    """Compute B(d) for all d = 1,...,p-1."""
    B = [0] * p
    for a in S_set:
        for b in S_set:
            if a != b:
                B[(a - b) % p] += 1
    return B


def get_all_symmetric_d11(p):
    """Generate all symmetric D11 of size n = (p+1)/2."""
    n = (p + 1) // 2
    pairs = []
    seen = set()
    for d in range(1, p):
        if d not in seen:
            pairs.append((d, p - d))
            seen.add(d)
            seen.add(p - d)

    reps_needed = n // 2
    all_d11 = []
    for combo in combinations(range(len(pairs)), reps_needed):
        D11 = []
        for idx in combo:
            d1, d2 = pairs[idx]
            D11.extend([d1, d2])
        if len(D11) == n:
            all_d11.append(frozenset(D11))
    return all_d11


def count_valid_d12(D11_set, p):
    """Count valid D12 for a given D11."""
    n = (p + 1) // 2
    k = (p - 3) // 2
    D11_set = set(D11_set)
    D22_set = set(range(1, p)) - D11_set
    A = compute_A(D11_set, p)

    # Thresholds
    tau = [0] * p
    for d in range(1, p):
        if d in D11_set:
            tau[d] = (p - 3) // 2 - A[d]
        else:
            tau[d] = (p + 3) // 2 - A[d]

    count = 0
    elements = list(range(1, p))
    for S in combinations(elements, k):
        S_full = frozenset({0} | set(S))
        B = compute_B(S_full, p)
        valid = True
        for d in range(1, p):
            if B[d] > tau[d]:
                valid = False
                break
        if valid:
            count += 1
    return count


def analyze_character_structure(p):
    """Analyze the character-sum structure of the second moment."""
    print(f"\n{'='*70}")
    print(f"p = {p}")
    print(f"{'='*70}")

    n = (p + 1) // 2
    k = (p - 3) // 2
    QR = get_qr(p)
    QNR = set(range(1, p)) - QR - {0}

    print(f"n = {n}, k = {k}")
    print(f"|QR| = {len(QR)}, |QNR| = {len(QNR)}")

    # Generate all symmetric D11
    all_d11 = get_all_symmetric_d11(p)
    print(f"Total symmetric D11: {len(all_d11)}")

    # For each D11, count valid D12 and compute character-theoretic data
    results = []
    for D11 in all_d11:
        D11_set = set(D11)
        N = count_valid_d12(D11_set, p)

        # Character-theoretic: how much does D11 overlap with QR?
        qr_overlap = len(D11_set & QR)

        # The "bias" of D11: sum of Legendre symbols
        chi_sum = sum(legendre_correct(d, p) for d in D11_set)

        # Also compute A-flatness
        A = compute_A(D11_set, p)
        max_A_d11 = max(A[d] for d in D11_set)
        max_A_d22 = max(A[d] for d in range(1, p) if d not in D11_set)

        results.append({
            'D11': sorted(D11),
            'N': N,
            'qr_overlap': qr_overlap,
            'chi_sum': chi_sum,
            'max_A_d11': max_A_d11,
            'max_A_d22': max_A_d22,
        })

    # Sort by N
    results.sort(key=lambda r: -r['N'])

    print(f"\nD11 sorted by N(D11):")
    print(f"{'N':>6} {'|D11∩QR|':>10} {'chi_sum':>10} {'max_A(D11)':>12} {'max_A(D22)':>12}")
    for r in results:
        print(f"{r['N']:6d} {r['qr_overlap']:10d} {r['chi_sum']:10d} {r['max_A_d11']:12d} {r['max_A_d22']:12d}")

    # Key question: does chi_sum predict N?
    # For p = 3 mod 4, the Paley graph adjacency is eta(a-b).
    # D11 = QR would give chi_sum = (p-1)/2 (all +1).
    # D11 = QNR would give chi_sum = -(p-1)/2 (all -1).

    # Correlation between chi_sum and N
    chi_sums = [r['chi_sum'] for r in results]
    N_vals = [r['N'] for r in results]

    if len(set(N_vals)) > 1 and len(set(chi_sums)) > 1:
        corr = np.corrcoef(chi_sums, N_vals)[0, 1]
        print(f"\nCorrelation(chi_sum, N) = {corr:.4f}")

    # Group by chi_sum
    by_chi = defaultdict(list)
    for r in results:
        by_chi[r['chi_sum']].append(r['N'])

    print(f"\nN by chi_sum:")
    for chi in sorted(by_chi.keys()):
        vals = by_chi[chi]
        print(f"  chi_sum = {chi:4d}: N values = {sorted(set(vals))}, count = {len(vals)}")

    # Analyze the Fourier-analytic relationship
    # For D11 with indicator function 1_{D11}, the DFT is:
    # hat{D11}(t) = sum_{d in D11} omega^{td}
    # For symmetric D11: hat{D11}(t) is real.
    #
    # The autocorrelation A(d) = (1/p) sum_t |hat{D11}(t)|^2 omega^{-td}
    #
    # For D11 = QR: hat{QR}(t) = sum_{d in QR} omega^{td} = gauss sum variant
    #   = (eta(t) * sqrt(p) - 1) / 2  for t != 0
    #   and hat{QR}(0) = (p-1)/2

    print(f"\n--- Fourier analysis ---")
    omega = np.exp(2j * np.pi / p)

    for i, r in enumerate(results[:min(5, len(results))]):
        D11 = r['D11']
        # Compute DFT
        dft = np.zeros(p, dtype=complex)
        for t in range(p):
            for d in D11:
                dft[t] += omega ** (t * d)

        # Power spectrum
        ps = np.abs(dft) ** 2

        # Is it real (symmetric D11)?
        is_real = np.allclose(dft.imag, 0, atol=1e-8)

        # Compare to Paley: |hat{QR}(t)|^2 = (p + eta(t)*sqrt(p) - ...)/4 for t != 0
        # Actually |hat{QR}(t)|^2 = (p - 2*eta(t)*sqrt(p) + 1)/4 for t != 0 when p=3 mod 4
        # Let's just compute the power spectrum statistics
        ps_nonzero = ps[1:]

        print(f"\n  D11 #{i}: N = {r['N']}, chi_sum = {r['chi_sum']}")
        print(f"    hat{{D11}}(0) = {dft[0].real:.1f}")
        print(f"    |hat{{D11}}|^2 stats (t>0): mean={np.mean(ps_nonzero):.2f}, "
              f"std={np.std(ps_nonzero):.2f}, min={np.min(ps_nonzero):.2f}, max={np.max(ps_nonzero):.2f}")
        print(f"    is_real: {is_real}")

        # Legendre-weighted sum: sum_t eta(t) * |hat{D11}(t)|^2
        leg_sum = sum(legendre_correct(t, p) * ps[t].real for t in range(1, p))
        print(f"    sum_t eta(t)*|hat(t)|^2 = {leg_sum:.2f}")

    return results


def analyze_switchings(p):
    """Analyze the effect of 'switching' one pair in D11.

    Key idea for the second moment: if switching one pair {d, p-d} in/out
    of D11 changes N(D11) by at most a constant factor, then the ratio
    E[N^2]/E[N]^2 is bounded.
    """
    print(f"\n{'='*70}")
    print(f"SWITCHING ANALYSIS: p = {p}")
    print(f"{'='*70}")

    all_d11 = get_all_symmetric_d11(p)

    # Compute N for all D11
    d11_N = {}
    for D11 in all_d11:
        N = count_valid_d12(set(D11), p)
        d11_N[D11] = N

    # Build the switching graph: two D11 are adjacent if they differ by one pair
    pairs = []
    seen = set()
    for d in range(1, p):
        if d not in seen:
            pairs.append(frozenset([d, p - d]))
            seen.add(d)
            seen.add(p - d)

    print(f"\nSwitching ratios (N(D11')/N(D11) when differing by one pair):")
    ratios = []
    for i, D11_1 in enumerate(all_d11):
        for j, D11_2 in enumerate(all_d11):
            if i >= j:
                continue
            sym_diff = D11_1 ^ D11_2  # symmetric difference
            if len(sym_diff) == 2:  # differ by exactly one pair
                N1 = d11_N[D11_1]
                N2 = d11_N[D11_2]
                if N1 > 0 and N2 > 0:
                    ratio = max(N1, N2) / min(N1, N2)
                    ratios.append(ratio)
                elif N1 > 0 or N2 > 0:
                    ratios.append(float('inf'))

    finite_ratios = [r for r in ratios if r < float('inf')]
    inf_count = sum(1 for r in ratios if r == float('inf'))

    if finite_ratios:
        print(f"  Finite ratios: min={min(finite_ratios):.4f}, max={max(finite_ratios):.4f}, "
              f"mean={np.mean(finite_ratios):.4f}, count={len(finite_ratios)}")
    print(f"  Infinite ratios (one side N=0): {inf_count}")
    print(f"  Total switching pairs: {len(ratios)}")

    # More detail: group by A-flatness
    n = (p + 1) // 2
    E_A = (p + 1) / 4

    aflat_ratios = []
    for i, D11_1 in enumerate(all_d11):
        for j, D11_2 in enumerate(all_d11):
            if i >= j:
                continue
            sym_diff = D11_1 ^ D11_2
            if len(sym_diff) == 2:
                # Check if both are A-flat
                A1 = compute_A(set(D11_1), p)
                A2 = compute_A(set(D11_2), p)
                max_A1 = max(A1[d] for d in D11_1)
                max_A2 = max(A2[d] for d in D11_2)
                if max_A1 <= int(E_A) and max_A2 <= int(E_A):
                    N1 = d11_N[D11_1]
                    N2 = d11_N[D11_2]
                    if N1 > 0 and N2 > 0:
                        ratio = max(N1, N2) / min(N1, N2)
                        aflat_ratios.append(ratio)
                    elif N1 > 0 or N2 > 0:
                        aflat_ratios.append(float('inf'))

    if aflat_ratios:
        finite = [r for r in aflat_ratios if r < float('inf')]
        inf_c = sum(1 for r in aflat_ratios if r == float('inf'))
        print(f"\n  Among A-flat D11:")
        if finite:
            print(f"    Finite switching ratios: min={min(finite):.4f}, max={max(finite):.4f}")
        print(f"    Infinite (one side N=0): {inf_c}")
        print(f"    Total A-flat switching pairs: {len(aflat_ratios)}")


def explore_orbit_structure(p):
    """Explore the multiplicative orbit structure.

    Z_p^* acts on D11 by multiplication: g * D11 = {g*d mod p : d in D11}.
    For D11 made of complementary pairs, multiplication by g preserves
    the symmetric structure iff g = ±1 mod p or g maps pairs to pairs.

    Actually: if D11 is symmetric, then -D11 = D11. Multiplication by g
    sends D11 to gD11, which is symmetric iff -gD11 = gD11 iff g(-D11) = gD11.
    Since -D11 = D11, this is always true. So gD11 is symmetric for all g.

    Key: the (p-1)/2 cosets of {+1,-1} in Z_p^* act on symmetric D11.
    If D11 is a union of cosets, then gD11 = D11 for g in those cosets.

    The QR set is exactly the set of cosets of {1,-1} that form a subgroup.
    """
    print(f"\n{'='*70}")
    print(f"ORBIT STRUCTURE: p = {p}")
    print(f"{'='*70}")

    n = (p + 1) // 2
    all_d11 = get_all_symmetric_d11(p)

    # Compute N for all D11
    d11_N = {}
    for D11 in all_d11:
        N = count_valid_d12(set(D11), p)
        d11_N[D11] = N

    # Orbits under multiplicative action
    # The group Z_p^* / {±1} has order (p-1)/2.
    # A generator of Z_p^* gives orbits.

    # Find primitive root
    def primitive_root(p):
        for g in range(2, p):
            s = set()
            val = 1
            for _ in range(p - 1):
                s.add(val)
                val = (val * g) % p
            if len(s) == p - 1:
                return g
        return None

    g = primitive_root(p)
    print(f"Primitive root: {g}")

    # Group elements of Z_p^* / {±1}
    # Cosets: {g^i, -g^i} = {g^i, g^i * g^{(p-1)/2}}
    # Since g^{(p-1)/2} = -1 (mod p) for prime p.

    # Action on D11: g * D11 = {g*d mod p : d in D11}
    # Since D11 is symmetric, g * D11 = {g*d mod p : d in D11} is also symmetric.

    # Compute orbits
    visited = set()
    orbits = []
    for D11 in all_d11:
        if D11 in visited:
            continue
        orbit = []
        current = D11
        for _ in range((p - 1) // 2):
            # Apply g (mod the {±1} quotient)
            new_d11 = frozenset((g * d) % p for d in current)
            if new_d11 not in d11_N:
                break
            orbit.append(new_d11)
            visited.add(new_d11)
            current = new_d11
        orbits.append(orbit)

    print(f"Number of orbits: {len(orbits)}")
    for i, orbit in enumerate(orbits):
        N_vals = [d11_N[D11] for D11 in orbit]
        print(f"  Orbit {i}: size {len(orbit)}, N values = {sorted(set(N_vals))}")

    # Key observation: N should be CONSTANT on orbits!
    # Because multiplying D11 by g also multiplies D12 by g, preserving the constraint.
    # If N(D11) = k, then N(gD11) = k because the valid D12 are exactly {g*D12 : D12 valid for D11}.

    all_constant = True
    for orbit in orbits:
        N_vals = [d11_N[D11] for D11 in orbit]
        if len(set(N_vals)) > 1:
            all_constant = False
            print(f"  WARNING: non-constant N on orbit! Values: {N_vals}")

    if all_constant:
        print(f"\nN(D11) is CONSTANT on all orbits (as expected by multiplicative symmetry)")

    # This means: E[N^2]/E[N]^2 depends only on the ORBIT structure
    # If there are K orbits with sizes s_1,...,s_K and N-values n_1,...,n_K:
    # E[N] = sum_i s_i * n_i / total
    # E[N^2] = sum_i s_i * n_i^2 / total
    # Ratio = E[N^2]/E[N]^2 = (total * sum s_i n_i^2) / (sum s_i n_i)^2

    total = len(all_d11)
    E_N = sum(len(orbit) * d11_N[orbit[0]] for orbit in orbits) / total
    E_N2 = sum(len(orbit) * d11_N[orbit[0]]**2 for orbit in orbits) / total
    ratio = E_N2 / E_N**2 if E_N > 0 else float('inf')

    print(f"\nE[N] = {E_N:.4f}")
    print(f"E[N^2] = {E_N2:.4f}")
    print(f"Ratio = {ratio:.4f}")

    # How many distinct N-values are there?
    distinct_N = sorted(set(d11_N[orbit[0]] for orbit in orbits))
    print(f"Distinct N values: {distinct_N}")
    print(f"Number of distinct N values: {len(distinct_N)}")

    return orbits, d11_N


def main():
    for p in [11, 19, 23]:
        if p % 4 != 3:
            continue
        analyze_character_structure(p)
        analyze_switchings(p)
        explore_orbit_structure(p)


if __name__ == '__main__':
    main()
