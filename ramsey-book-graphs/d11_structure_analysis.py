"""
Structural analysis of LP-feasible D11 sets.

Goal: Understand what makes a D11 LP-feasible, and identify algebraic
patterns that could lead to an explicit construction for all p ≡ 3 mod 4.

Key questions:
1. What is the Fourier spectrum of LP-feasible D11?
2. How do LP-feasible D11 relate to cyclotomic cosets?
3. Is there a pattern in which symmetric pairs are selected?
"""

import numpy as np
import random
from scipy.optimize import linprog
from collections import Counter


def quadratic_residues(p):
    return set((x * x) % p for x in range(1, p))


def primitive_root(p):
    """Find a primitive root mod p."""
    for g in range(2, p):
        seen = set()
        x = 1
        for _ in range(p - 1):
            x = (x * g) % p
            seen.add(x)
        if len(seen) == p - 1:
            return g
    return None


def discrete_log(x, g, p):
    """Compute log_g(x) mod (p-1)."""
    val = 1
    for i in range(p - 1):
        if val == x:
            return i
        val = (val * g) % p
    return None


def autocorrelation(S, m):
    indicator = np.zeros(m, dtype=int)
    for x in S:
        indicator[x % m] = 1
    result = np.zeros(m, dtype=int)
    for d in range(m):
        count = 0
        for x in range(m):
            if indicator[x] and indicator[(x + d) % m]:
                count += 1
        result[d] = count
    return result


def fourier_spectrum(S, m):
    indicator = np.zeros(m, dtype=complex)
    for x in S:
        indicator[x % m] = 1.0
    dft = np.fft.fft(indicator)
    return np.abs(dft) ** 2


def check_lp_feasibility(D11, p, n):
    """Check if a symmetric D11 admits an LP-feasible D12 spectrum."""
    k = len(D11)
    d12_size = n - 1
    d11_thresh = n - 2
    d22_thresh = 2 * k - n + 1

    A = autocorrelation(D11, p)

    omega = np.exp(2j * np.pi / p)
    nvar = p - 1
    F = np.zeros((nvar, nvar))
    for d in range(1, p):
        for j in range(1, p):
            F[d - 1, j - 1] = np.real(omega ** (j * d)) / p

    base = d12_size ** 2 / p
    total_Q = d12_size * (d12_size - 1)

    b_ub = np.zeros(nvar)
    for d in range(1, p):
        if d in D11:
            b_ub[d - 1] = d11_thresh - A[d] - base
        else:
            b_ub[d - 1] = d22_thresh - A[d] - base

    c = np.zeros(nvar)
    A_eq = np.ones((1, nvar))
    b_eq = np.array([total_Q])
    bounds = [(0, None)] * nvar

    result = linprog(c, A_ub=F, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                     bounds=bounds, method='highs')

    if result.success:
        # Also compute the margin: how much slack in the tightest constraint?
        Q_opt = result.x
        B_opt = np.array([base + F[d] @ Q_opt for d in range(nvar)])
        margins = []
        for d in range(1, p):
            ab = A[d] + B_opt[d - 1]
            thresh = d11_thresh if d in D11 else d22_thresh
            margins.append(thresh - ab)
        return True, min(margins), Q_opt
    return False, None, None


def analyze_known_d11_cyclotomic(p, D11):
    """Analyze D11 in terms of cyclotomic structure."""
    g = primitive_root(p)
    qr = quadratic_residues(p)

    # Compute discrete logs
    dlogs = {}
    for x in range(1, p):
        dlogs[x] = discrete_log(x, g, p)

    # Symmetric pairs: {d, -d mod p}
    pairs = []
    seen = set()
    for x in range(1, p):
        neg_x = (p - x) % p
        pair = (min(x, neg_x), max(x, neg_x))
        if pair not in seen:
            seen.add(pair)
            pairs.append(pair)

    # For each pair, which element is QR?
    pair_info = []
    for a, b in pairs:
        a_qr = a in qr
        b_qr = b in qr
        in_d11 = a in D11
        pair_info.append({
            'pair': (a, b),
            'a_qr': a_qr,
            'b_qr': b_qr,
            'in_D11': in_d11,
            'a_dlog': dlogs.get(a, None),
            'b_dlog': dlogs.get(b, None),
        })

    # Pattern: do LP-feasible D11 prefer certain cyclotomic classes?
    in_d11_count = sum(1 for pi in pair_info if pi['in_D11'])
    not_in_d11_count = sum(1 for pi in pair_info if not pi['in_D11'])

    # Count by quartile of discrete log
    qrt = (p - 1) // 4
    quartile_counts = Counter()
    for pi in pair_info:
        # Use the QR element's discrete log
        a, b = pi['pair']
        qr_elem = a if pi['a_qr'] else b
        dl = dlogs[qr_elem]
        q = dl // qrt if qrt > 0 else 0
        if pi['in_D11']:
            quartile_counts[('in', q)] += 1
        else:
            quartile_counts[('out', q)] += 1

    return pair_info, quartile_counts


def main():
    print("=" * 80)
    print("CYCLOTOMIC STRUCTURE OF KNOWN SOLUTIONS")
    print("=" * 80)

    # Known solutions for p ≡ 3 mod 4
    solutions = {
        43: set([1, 2, 5, 10, 11, 13, 16, 17, 18, 19, 20, 23, 24, 25, 26, 27, 30, 32, 33, 38, 41, 42]),
        47: set([2, 6, 8, 10, 11, 13, 14, 16, 19, 20, 21, 26, 27, 28, 31, 33, 34, 36, 37, 39, 41, 45]),
        59: set([1, 2, 3, 5, 6, 8, 9, 10, 11, 17, 19, 25, 27, 29, 30, 32, 34, 40, 42, 48, 49, 50, 51, 53, 54, 56, 57, 58]),
    }

    for p, D11 in sorted(solutions.items()):
        n = (p + 1) // 2
        g = primitive_root(p)
        qr = quadratic_residues(p)

        print(f"\np={p}, n={n}, g={g}, |D11|={len(D11)}")

        pair_info, qcounts = analyze_known_d11_cyclotomic(p, D11)

        # How many pairs with QR element "small" are in D11?
        selected_pairs = [pi for pi in pair_info if pi['in_D11']]
        not_selected = [pi for pi in pair_info if not pi['in_D11']]

        # For each selected pair, is the QR element the smaller one?
        qr_small_in = sum(1 for pi in selected_pairs if pi['a_qr'])
        qr_large_in = sum(1 for pi in selected_pairs if pi['b_qr'])
        qr_small_out = sum(1 for pi in not_selected if pi['a_qr'])
        qr_large_out = sum(1 for pi in not_selected if pi['b_qr'])

        print(f"  Selected pairs: {len(selected_pairs)}/{len(pair_info)}")
        print(f"  QR-element-is-smaller: in={qr_small_in}, out={qr_small_out}")
        print(f"  QR-element-is-larger:  in={qr_large_in}, out={qr_large_out}")

        # Discrete log analysis
        dlogs = {}
        for x in range(1, p):
            dlogs[x] = discrete_log(x, g, p)

        # For selected pairs, what are the discrete logs of the QR element?
        qr_dlogs_in = []
        qr_dlogs_out = []
        for pi in pair_info:
            a, b = pi['pair']
            qr_elem = a if pi['a_qr'] else b
            dl = dlogs[qr_elem]
            if pi['in_D11']:
                qr_dlogs_in.append(dl)
            else:
                qr_dlogs_out.append(dl)

        print(f"  DLog of QR element in selected pairs: {sorted(qr_dlogs_in)}")
        print(f"  DLog of QR element in NOT selected:   {sorted(qr_dlogs_out)}")

        # Are the selected dlogs a structured subset of even numbers?
        # QR elements have even dlogs (by definition: QR iff dlog is even)
        # So qr_dlogs should all be even
        even_in = [d for d in qr_dlogs_in if d % 4 == 0]
        odd_in = [d for d in qr_dlogs_in if d % 4 == 2]  # Even but ≡ 2 mod 4
        print(f"  QR dlogs ≡ 0 mod 4 (in): {len(even_in)}/{len(qr_dlogs_in)}")
        print(f"  QR dlogs ≡ 2 mod 4 (in): {len(odd_in)}/{len(qr_dlogs_in)}")

    # Now: exhaustive search for LP-feasible D11 at small primes
    print("\n" + "=" * 80)
    print("EXHAUSTIVE SEARCH FOR LP-FEASIBLE D11 (small primes)")
    print("=" * 80)

    for p in [23, 31]:
        n = (p + 1) // 2

        # Standard |D11| = n (or nearest even)
        k = n if n % 2 == 0 else n - 1

        # Build pairs
        pairs = []
        for x in range(1, (p + 1) // 2):
            pairs.append((x, p - x))
        total_pairs = len(pairs)
        num_pairs_needed = k // 2

        print(f"\np={p}, n={n}, |D11|={k}, choosing {num_pairs_needed}/{total_pairs} pairs")

        # For p=23, exhaustive is C(11,6)=462 — feasible
        from itertools import combinations

        total_combos = 1
        for i in range(num_pairs_needed):
            total_combos = total_combos * (total_pairs - i) // (i + 1)

        if total_combos > 50000:
            print(f"  Too many combinations ({total_combos}), sampling...")
            # Sample instead
            feasible_d11s = []
            for trial in range(2000):
                selected = random.sample(range(total_pairs), num_pairs_needed)
                D11 = set()
                for i in selected:
                    D11.add(pairs[i][0])
                    D11.add(pairs[i][1])

                is_feasible, margin, _ = check_lp_feasibility(D11, p, n)
                if is_feasible:
                    feasible_d11s.append((sorted(D11), margin))

            print(f"  Found {len(feasible_d11s)}/2000 LP-feasible D11")
        else:
            print(f"  Exhaustive search over {total_combos} combinations...")
            feasible_d11s = []
            for combo in combinations(range(total_pairs), num_pairs_needed):
                D11 = set()
                for i in combo:
                    D11.add(pairs[i][0])
                    D11.add(pairs[i][1])

                is_feasible, margin, _ = check_lp_feasibility(D11, p, n)
                if is_feasible:
                    feasible_d11s.append((sorted(D11), margin))

            print(f"  Found {len(feasible_d11s)}/{total_combos} LP-feasible D11 ({100*len(feasible_d11s)/total_combos:.1f}%)")

        if feasible_d11s:
            # Analyze the feasible D11s
            feasible_d11s.sort(key=lambda x: -x[1])  # Sort by margin, best first
            print(f"  Best margin: {feasible_d11s[0][1]:.4f}")
            print(f"  Worst margin: {feasible_d11s[-1][1]:.4f}")

            # Check: what's the autocorrelation spectrum of the best D11?
            best_D11 = set(feasible_d11s[0][0])
            A = autocorrelation(best_D11, p)
            S = fourier_spectrum(best_D11, p)

            print(f"\n  Best D11: {feasible_d11s[0][0]}")
            print(f"  A(d) at D11: {[A[d] for d in range(1,p) if d in best_D11]}")
            max_A_d11 = max(A[d] for d in range(1, p) if d in best_D11)
            print(f"  max A(d) at D11: {max_A_d11} (QR-constant would give {(p-3)//4})")
            print(f"  |D11_hat|^2 range: [{S[1:].min():.2f}, {S[1:].max():.2f}]")
            print(f"  |D11_hat|^2 std: {S[1:].std():.2f}")

            # Cyclotomic analysis of best D11
            g = primitive_root(p)
            qr = quadratic_residues(p)
            dlogs = {}
            for x in range(1, p):
                dlogs[x] = discrete_log(x, g, p)

            print(f"\n  Cyclotomic analysis (g={g}):")
            for d11_sorted, margin in feasible_d11s[:3]:
                d11_set = set(d11_sorted)
                # Which pairs are selected?
                qr_dlogs = []
                for x in range(1, (p + 1) // 2):
                    pair = (x, p - x)
                    qr_elem = x if x in qr else p - x
                    if x in d11_set:
                        qr_dlogs.append(dlogs[qr_elem])
                print(f"    D11 QR dlogs: {sorted(qr_dlogs)} (margin={margin:.4f})")


if __name__ == '__main__':
    main()
