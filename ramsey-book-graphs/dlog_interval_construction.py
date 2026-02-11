"""
Test explicit algebraic constructions for D11 based on discrete log intervals.

Idea: For p ≡ 3 mod 4, each symmetric pair {d, -d} has one QR and one QNR element.
Select pairs based on the QR element's discrete log position.

Constructions tested:
A) Contiguous interval: QR dlogs in {0, 2, ..., 2(k-1)} for some starting position
B) Biquadratic: QR dlogs ≡ 0 mod 4 (i.e., fourth powers)
C) Alternating: Every other even dlog
D) Half-period: dlogs in first half of the cycle
"""

import numpy as np
from scipy.optimize import linprog


def primitive_root(p):
    for g in range(2, p):
        seen = set()
        x = 1
        for _ in range(p - 1):
            x = (x * g) % p
            seen.add(x)
        if len(seen) == p - 1:
            return g
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


def check_lp_feasibility(D11_set, p, n):
    k = len(D11_set)
    d12_size = n - 1
    d11_thresh = n - 2
    d22_thresh = 2 * k - n + 1

    A = autocorrelation(D11_set, p)

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
        if d in D11_set:
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
        Q_opt = result.x
        B_opt = np.array([base + F[d] @ Q_opt for d in range(nvar)])
        margins = []
        for d in range(1, p):
            ab = A[d] + B_opt[d - 1]
            thresh = d11_thresh if d in D11_set else d22_thresh
            margins.append(thresh - ab)
        return True, min(margins)
    return False, None


def build_d11_from_dlog_set(dlog_set, g, p):
    """Build D11 from a set of even discrete logs of the QR elements in selected pairs."""
    D11 = set()
    for ell in dlog_set:
        d = pow(g, ell, p)  # QR element
        neg_d = (p - d) % p  # QNR element
        D11.add(d)
        D11.add(neg_d)
    return D11


def test_constructions():
    print("=" * 80)
    print("ALGEBRAIC D11 CONSTRUCTIONS: Feasibility Test")
    print("=" * 80)

    for p in [7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83]:
        n = (p + 1) // 2
        g = primitive_root(p)
        total_pairs = (p - 1) // 2
        target_pairs = (p + 1) // 4  # for |D11| = (p+1)/2

        # Even dlogs: 0, 2, 4, ..., p-3
        even_dlogs = list(range(0, p - 1, 2))
        assert len(even_dlogs) == total_pairs

        print(f"\np={p}, n={n}, g={g}, pairs={total_pairs}, target={target_pairs}")

        best_constructions = []

        # Construction A: Contiguous intervals of even dlogs
        # Try all starting positions
        for start in range(total_pairs):
            dlog_set = set()
            for i in range(target_pairs):
                idx = (start + i) % total_pairs
                dlog_set.add(even_dlogs[idx])

            D11 = build_d11_from_dlog_set(dlog_set, g, p)
            A = autocorrelation(D11, p)

            # Quick check: max A at D11 positions
            max_A_d11 = max(A[d] for d in range(1, p) if d in D11)
            qr_const = (p - 3) // 4

            # Full LP check only if max_A seems reasonable
            is_feas, margin = check_lp_feasibility(D11, p, n)

            if is_feas:
                best_constructions.append(('interval', start, margin, sorted(dlog_set)))

        # Construction B: Biquadratic residues (dlogs ≡ 0 mod 4)
        if (p - 1) % 4 == 0:
            bqr_dlogs = set(range(0, p - 1, 4))  # Fourth powers
            if len(bqr_dlogs) >= target_pairs:
                # Take first target_pairs of them
                dlog_set = set(list(bqr_dlogs)[:target_pairs])
                D11 = build_d11_from_dlog_set(dlog_set, g, p)
                is_feas, margin = check_lp_feasibility(D11, p, n)
                if is_feas:
                    best_constructions.append(('biquadratic', 0, margin, sorted(dlog_set)))

        # Construction C: Spread selection (every other even dlog, i.e., dlogs ≡ 0 mod 4)
        for offset in [0, 2]:
            dlog_set = set()
            for ell in range(offset, p - 1, 4):
                dlog_set.add(ell)
                if len(dlog_set) >= target_pairs:
                    break
            if len(dlog_set) == target_pairs:
                D11 = build_d11_from_dlog_set(dlog_set, g, p)
                is_feas, margin = check_lp_feasibility(D11, p, n)
                if is_feas:
                    best_constructions.append(('spread', offset, margin, sorted(dlog_set)))

        if best_constructions:
            print(f"  LP-FEASIBLE constructions found: {len(best_constructions)}")
            for ctype, param, margin, dlogs in best_constructions[:3]:
                print(f"    {ctype}(param={param}): margin={margin:.6f}, dlogs={dlogs}")
        else:
            print(f"  No LP-feasible algebraic construction found")

            # Report the best interval construction anyway
            best_max_a = None
            for start in range(total_pairs):
                dlog_set = set()
                for i in range(target_pairs):
                    idx = (start + i) % total_pairs
                    dlog_set.add(even_dlogs[idx])

                D11 = build_d11_from_dlog_set(dlog_set, g, p)
                A = autocorrelation(D11, p)
                max_A_d11 = max(A[d] for d in range(1, p) if d in D11)

                if best_max_a is None or max_A_d11 < best_max_a:
                    best_max_a = max_A_d11
                    best_start = start

            print(f"  Best interval: start={best_start}, max A at D11 = {best_max_a} "
                  f"(threshold {(p-3)//4})")


def test_non_standard_sizes():
    """Test different |D11| sizes — the known solutions sometimes use k ≠ n."""
    print("\n" + "=" * 80)
    print("NON-STANDARD D11 SIZES: Testing k = n-2, n, n+2")
    print("=" * 80)

    for p in [23, 31, 43, 47, 59, 67]:
        n = (p + 1) // 2
        g = primitive_root(p)
        total_pairs = (p - 1) // 2
        even_dlogs = list(range(0, p - 1, 2))

        print(f"\np={p}, n={n}")

        for k in [n - 2, n, n + 2]:
            if k <= 0 or k % 2 != 0 or k > p - 1:
                continue
            target_pairs = k // 2
            if target_pairs > total_pairs:
                continue

            d11_thresh = n - 2
            d22_thresh = 2 * k - n + 1

            # Try interval constructions
            feasible_count = 0
            for start in range(total_pairs):
                dlog_set = set()
                for i in range(target_pairs):
                    idx = (start + i) % total_pairs
                    dlog_set.add(even_dlogs[idx])

                D11 = build_d11_from_dlog_set(dlog_set, g, p)
                is_feas, margin = check_lp_feasibility(D11, p, n)
                if is_feas:
                    feasible_count += 1

            print(f"  |D11|={k:3d}: D11_thresh={d11_thresh}, D22_thresh={d22_thresh}, "
                  f"interval-feasible: {feasible_count}/{total_pairs}")


if __name__ == '__main__':
    test_constructions()
    test_non_standard_sizes()
