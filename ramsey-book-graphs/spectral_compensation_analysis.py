"""
Spectral compensation analysis for the D12 construction problem.

Key question: Given a symmetric D11 of size k, what spectral properties
must D12 have to satisfy all constraints?

In Fourier space:
  A(d) + B(d) = (1/m) * sum_k P(k) * omega^{kd}
  where P(k) = |D11_hat(k)|^2 + |D12_hat(k)|^2

For the constraints, we need the inverse DFT of P(k) to be bounded:
  P_check(d) <= n-2       for d in D11
  P_check(d) <= 2k-n+1    for d in D22

This means |D12_hat(k)|^2 must compensate for |D11_hat(k)|^2.

Specifically, if we write T(k) = target P(k), then:
  |D12_hat(k)|^2 = T(k) - |D11_hat(k)|^2

For this to be realizable, T(k) - |D11_hat(k)|^2 >= 0 for all k,
AND there must exist an actual set D12 of size n-1 with this spectrum.
"""

import numpy as np
from collections import defaultdict


def quadratic_residues(p):
    return sorted(set((x * x) % p for x in range(1, p)))


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
    """Compute |D_hat(k)|^2 for k = 0, ..., m-1."""
    indicator = np.zeros(m, dtype=complex)
    for x in S:
        indicator[x % m] = 1.0
    dft = np.fft.fft(indicator)
    return np.abs(dft) ** 2


def analyze_spectral_compensation(n, D11, D12, p):
    """Analyze the spectral compensation between D11 and D12."""
    S11 = fourier_spectrum(D11, p)
    S12 = fourier_spectrum(D12, p)
    P = S11 + S12  # Total power spectrum

    # Inverse DFT of P gives m * (A(d) + B(d))
    # A(d) + B(d) = (1/m) * sum_k P(k) * exp(2pi*i*k*d/m)

    A = autocorrelation(D11, p)
    B = autocorrelation(D12, p)

    # Spectral correlation
    # Exclude k=0 (it's just sizes squared)
    s11_nz = S11[1:]
    s12_nz = S12[1:]
    corr = np.corrcoef(s11_nz, s12_nz)[0, 1]

    return {
        'S11': S11,
        'S12': S12,
        'P': P,
        'A': A,
        'B': B,
        'spectral_correlation': corr,
        'S11_mean': np.mean(S11[1:]),
        'S12_mean': np.mean(S12[1:]),
        'S11_std': np.std(S11[1:]),
        'S12_std': np.std(S12[1:]),
    }


def analyze_all_known_solutions():
    """Analyze spectral compensation in known solutions."""
    solutions = {
        6: {
            'D11': [1, 3, 4, 7, 8, 10],
            'D12': [0, 2, 3, 4, 5],
        },
        10: {
            'D11': [1, 4, 5, 6, 7, 9, 11, 16, 17],
            'D12': [0, 1, 2, 4, 5, 6, 7, 9, 11],
        },
        12: {
            'D11': [1, 2, 3, 5, 8, 9, 14, 15, 18, 20, 21, 22],
            'D12': [0, 1, 2, 3, 4, 5, 9, 14, 18, 19, 20],
        },
        22: {
            'D11': [1, 2, 5, 10, 11, 13, 16, 17, 18, 19, 20, 23, 24, 25, 26, 27, 30, 32, 33, 38, 41, 42],
            'D12': [0, 2, 5, 6, 8, 11, 15, 16, 20, 24, 25, 27, 28, 31, 32, 34, 35, 36, 37, 39, 41],
        },
        24: {
            'D11': [2, 6, 8, 10, 11, 13, 14, 16, 19, 20, 21, 26, 27, 28, 31, 33, 34, 36, 37, 39, 41, 45],
            'D12': [0, 1, 8, 9, 10, 12, 19, 23, 24, 25, 26, 28, 29, 32, 33, 35, 38, 39, 40, 42, 44, 45, 46],
        },
        30: {
            'D11': [1, 2, 3, 5, 6, 8, 9, 10, 11, 17, 19, 25, 27, 29, 30, 32, 34, 40, 42, 48, 49, 50, 51, 53, 54, 56, 57, 58],
            'D12': [0, 1, 4, 6, 7, 9, 10, 11, 13, 18, 20, 21, 23, 24, 25, 26, 32, 33, 36, 38, 39, 42, 43, 45, 49, 50, 51, 54, 55],
        },
    }

    # Also include the verified small cases from prior work
    # n=6 (m=11), n=10 (m=19), n=12 (m=23)

    print("=" * 80)
    print("SPECTRAL COMPENSATION ANALYSIS OF KNOWN SOLUTIONS")
    print("=" * 80)

    for n, sol in sorted(solutions.items()):
        p = 2 * n - 1
        D11 = sol['D11']
        D12 = sol['D12']
        qr = quadratic_residues(p)

        result = analyze_spectral_compensation(n, D11, D12, p)

        k = len(D11)
        d11_thresh = n - 2
        d22_thresh = 2 * k - n + 1

        # Check constraint satisfaction
        max_ab_d11 = max(result['A'][d] + result['B'][d] for d in range(1, p) if d in set(D11))
        max_ab_d22 = max(result['A'][d] + result['B'][d] for d in range(1, p) if d not in set(D11))

        # Spectral analysis of actual D12 vs QR
        qr_spec = fourier_spectrum(qr, p)

        print(f"\nn={n}, p={p}, |D11|={k}, |D12|={len(D12)}")
        print(f"  D11 threshold: A+B ≤ {d11_thresh}, actual max = {max_ab_d11}")
        print(f"  D22 threshold: A+B ≤ {d22_thresh}, actual max = {max_ab_d22}")
        print(f"  Spectral correlation(|D11_hat|^2, |D12_hat|^2): {result['spectral_correlation']:.4f}")
        print(f"  |D11_hat|^2 stats (k>0): mean={result['S11_mean']:.2f}, std={result['S11_std']:.2f}")
        print(f"  |D12_hat|^2 stats (k>0): mean={result['S12_mean']:.2f}, std={result['S12_std']:.2f}")

        # How much does B(d) compensate for A(d)?
        A_vals = result['A']
        B_vals = result['B']
        d11_set = set(D11)

        # At D11 positions where A is high, is B low?
        A_at_d11 = [(d, A_vals[d], B_vals[d]) for d in range(1, p) if d in d11_set]
        A_at_d11.sort(key=lambda x: x[1], reverse=True)

        print(f"  Top-5 A(d) at D11 positions (d, A, B, A+B):")
        for d, a, b in A_at_d11[:5]:
            print(f"    d={d:3d}: A={a:2d}, B={b:2d}, A+B={a+b:2d} {'!' if a+b > d11_thresh else '✓'}")

        # Correlation between A and B at D11 positions
        a_d11 = np.array([A_vals[d] for d in range(1, p) if d in d11_set])
        b_d11 = np.array([B_vals[d] for d in range(1, p) if d in d11_set])
        if len(a_d11) > 2:
            corr_d11 = np.corrcoef(a_d11, b_d11)[0, 1]
            print(f"  Correlation of A(d), B(d) at D11 positions: {corr_d11:.4f}")


def investigate_required_b_spectrum():
    """
    Given known D11, compute what |D12_hat(k)|^2 MUST be for constraints to hold.

    If A(d) + B(d) <= T(d) for appropriate thresholds, then in Fourier space:
    sum_k P(k) omega^{kd} / m <= T(d)

    This means P(k) must satisfy a system of linear constraints on its inverse DFT.
    The required |D12_hat(k)|^2 = P(k) - |D11_hat(k)|^2.
    """
    from scipy.optimize import linprog

    solutions = {
        22: {
            'D11': [1, 2, 5, 10, 11, 13, 16, 17, 18, 19, 20, 23, 24, 25, 26, 27, 30, 32, 33, 38, 41, 42],
        },
        30: {
            'D11': [1, 2, 3, 5, 6, 8, 9, 10, 11, 17, 19, 25, 27, 29, 30, 32, 34, 40, 42, 48, 49, 50, 51, 53, 54, 56, 57, 58],
        },
    }

    print("\n" + "=" * 80)
    print("FOURIER LP: What |D12_hat(k)|^2 is needed?")
    print("=" * 80)

    for n, sol in sorted(solutions.items()):
        p = 2 * n - 1
        D11 = set(sol['D11'])
        k = len(D11)

        d11_thresh = n - 2
        d22_thresh = 2 * k - n + 1

        # Compute A(d) for known D11
        A = autocorrelation(D11, p)
        S11 = fourier_spectrum(D11, p)

        # The constraint: A(d) + B(d) <= thresh(d)
        # B(d) = (1/m) sum_{j=1}^{m-1} Q(j) * cos(2*pi*j*d/m)  + |D12|^2/m
        # where Q(j) = |D12_hat(j)|^2 for j > 0
        #
        # We want to find Q(1), ..., Q(m-1) >= 0 such that:
        # B(d) <= thresh(d) - A(d) for each relevant d
        #
        # Also: sum_j Q(j) = |D12|^2 - |D12| (from Parseval + B(0) = |D12|)
        # And: Q(j) >= 0 for all j

        d12_size = n - 1
        avg_Q = (d12_size * (d12_size - 1)) / (p - 1)  # average |D12_hat(k)|^2 for k > 0

        print(f"\nn={n}, p={p}, |D11|={k}, |D12|={d12_size}")
        print(f"  avg |D12_hat(k)|^2 for k>0 = {avg_Q:.2f}")
        print(f"  |D11_hat(k)|^2 for k>0: mean={np.mean(S11[1:]):.2f}, max={np.max(S11[1:]):.2f}")

        # Required B(d) bounds
        required_B = {}
        for d in range(1, p):
            if d in D11:
                required_B[d] = d11_thresh - A[d]
            else:
                required_B[d] = d22_thresh - A[d]

        print(f"  Required B(d) at D11 positions: [{min(required_B[d] for d in D11)}, "
              f"{max(required_B[d] for d in D11)}]")
        print(f"  Required B(d) at D22 positions: [{min(required_B[d] for d in range(1,p) if d not in D11)}, "
              f"{max(required_B[d] for d in range(1,p) if d not in D11)}]")

        # Solve LP: minimize max Q(j) subject to constraints
        # Variables: Q(1), ..., Q(p-1)  (these are |D12_hat(k)|^2)
        nvar = p - 1

        # B(d) = |D12|^2/p + (1/p) * sum_{j=1}^{p-1} Q(j) * Re(omega^{jd})
        # where omega = exp(2*pi*i/p)

        omega = np.exp(2j * np.pi / p)
        F = np.zeros((p - 1, nvar))  # F[d-1, j-1] = Re(omega^{j*d}) / p
        for d in range(1, p):
            for j in range(1, p):
                F[d - 1, j - 1] = np.real(omega ** (j * d)) / p

        # B(d) = |D12|^2/p + F[d-1] @ Q
        base = d12_size ** 2 / p

        # Constraint: B(d) <= required_B[d]
        # F[d-1] @ Q <= required_B[d] - base

        # Also: sum Q(j) = |D12|^2 - |D12| = d12_size * (d12_size - 1)
        total_Q = d12_size * (d12_size - 1)

        # Objective: minimize sum Q (feasibility) or minimize max Q
        # For feasibility, just check if LP is feasible
        c = np.zeros(nvar)  # Feasibility check

        A_ub = F  # Each row: F[d-1] @ Q <= rhs
        b_ub = np.array([required_B[d + 1] - base for d in range(p - 1)])

        # Equality: sum Q = total_Q
        A_eq = np.ones((1, nvar))
        b_eq = np.array([total_Q])

        # Bounds: Q >= 0
        bounds = [(0, None)] * nvar

        result = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                         bounds=bounds, method='highs')

        if result.success:
            Q_opt = result.x
            print(f"  LP FEASIBLE! Q range: [{Q_opt.min():.2f}, {Q_opt.max():.2f}]")
            print(f"  Q mean: {Q_opt.mean():.2f} (expected {avg_Q:.2f})")

            # Check: what B(d) does this give?
            B_lp = np.array([base + F[d] @ Q_opt for d in range(p - 1)])
            print(f"  LP B(d) range: [{B_lp.min():.2f}, {B_lp.max():.2f}]")

            # Required spectral compensation
            compensation = np.array([S11[j + 1] - avg_Q for j in range(nvar)])
            actual_compensation = np.array([Q_opt[j] - avg_Q for j in range(nvar)])
            corr_comp = np.corrcoef(S11[1:], Q_opt)[0, 1]
            print(f"  Correlation(|D11_hat|^2, Q_opt): {corr_comp:.4f}")

            # Compare Q_opt with QR spectrum
            qr = quadratic_residues(p)
            S_qr = fourier_spectrum(qr, p)
            print(f"  QR |D12_hat|^2 is constant at {S_qr[1]:.2f}")
            print(f"  Q_opt deviates from QR by: [{(Q_opt - S_qr[1:]).min():.2f}, "
                  f"{(Q_opt - S_qr[1:]).max():.2f}]")
        else:
            print(f"  LP INFEASIBLE for this D11!")


def investigate_d11_existence_via_lp():
    """
    The deepest question: for p ≡ 3 mod 4, does there always exist
    a symmetric D11 of size ~n such that the Fourier LP for D12 is feasible?

    We sample random symmetric D11 and check LP feasibility.
    """
    import random
    from scipy.optimize import linprog

    print("\n" + "=" * 80)
    print("D11 EXISTENCE: How often is the Fourier LP feasible for random D11?")
    print("=" * 80)

    for p in [11, 19, 23, 31, 43]:
        n = (p + 1) // 2
        d12_size = n - 1

        # Build pairs
        pairs = []
        for x in range(1, (p + 1) // 2):
            pairs.append((x, p - x))
        total_pairs = len(pairs)

        # Try different |D11| sizes
        for k in [n - 2, n, n + 2] if n + 2 <= p - 1 else [n - 2, n]:
            if k <= 0 or k % 2 != 0 or k // 2 > total_pairs:
                continue

            num_pairs_needed = k // 2
            d11_thresh = n - 2
            d22_thresh = 2 * k - n + 1

            trials = 200 if p <= 23 else 50
            feasible = 0

            omega = np.exp(2j * np.pi / p)
            F = np.zeros((p - 1, p - 1))
            for d in range(1, p):
                for j in range(1, p):
                    F[d - 1, j - 1] = np.real(omega ** (j * d)) / p

            base = d12_size ** 2 / p
            total_Q = d12_size * (d12_size - 1)

            for trial in range(trials):
                selected = random.sample(range(total_pairs), num_pairs_needed)
                D11 = set()
                for i in selected:
                    D11.add(pairs[i][0])
                    D11.add(pairs[i][1])

                A = autocorrelation(D11, p)

                # Build LP
                b_ub = np.zeros(p - 1)
                for d in range(1, p):
                    if d in D11:
                        b_ub[d - 1] = d11_thresh - A[d] - base
                    else:
                        b_ub[d - 1] = d22_thresh - A[d] - base

                c = np.zeros(p - 1)
                A_eq = np.ones((1, p - 1))
                b_eq = np.array([total_Q])
                bounds = [(0, None)] * (p - 1)

                result = linprog(c, A_ub=F, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                                 bounds=bounds, method='highs')

                if result.success:
                    feasible += 1

            rate = feasible / trials * 100
            print(f"  p={p:3d}, |D11|={k:3d}: {feasible}/{trials} LP-feasible ({rate:.1f}%)")


if __name__ == '__main__':
    analyze_all_known_solutions()
    investigate_required_b_spectrum()
    investigate_d11_existence_via_lp()
