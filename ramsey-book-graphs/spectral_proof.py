"""
Spectral analysis toward a proof of the key lemma.

KEY QUESTION: For a random symmetric D11 of size (p+1)/2, what is the distribution
of max |D̂11(k)|²? If we can show this is bounded by O(p), then the LP is feasible.

APPROACH: D11 is a union of (p+1)/4 negation pairs from (p-1)/2 available pairs.
D̂11(k) = 2 Σ_{j ∈ S} cos(2πkd_j/p) where S ⊂ [(p-1)/2] with |S| = (p+1)/4.

This is a sum of random cosines. We analyze:
1. E[|D̂11(k)|²] — expected spectral power
2. max_k |D̂11(k)|² — spectral peak
3. A(d) = (1/p) Σ_k |D̂11(k)|² ω^{kd} — autocorrelation from spectrum
4. The LP feasibility condition in terms of spectral properties
"""

import sys
import os
import math
import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def get_pairs(p):
    pairs = []
    seen = set()
    for d in range(1, p):
        if d not in seen:
            pairs.append((d, (-d) % p))
            seen.add(d)
            seen.add((-d) % p)
    return pairs


def d11_spectrum(p, D11):
    """Compute |D̂11(k)|² for all k."""
    ind = np.zeros(p, dtype=np.float64)
    for j in D11:
        ind[j] = 1.0
    dft = np.fft.fft(ind)
    return np.abs(dft) ** 2


def theoretical_spectrum_moments(p):
    """
    Compute theoretical E[|D̂11(k)|²] and Var[|D̂11(k)|²] for random D11.

    D11 = union of s = (p+1)/4 pairs from F = (p-1)/2 available pairs.
    D̂11(k) = 2 Σ_{j ∈ S} cos(2πkd_j/p) where d_j are the pair representatives.

    |D̂11(k)|² = 4 (Σ_{j ∈ S} cos(2πkd_j/p))² = 4 Σ_{j,j' ∈ S} c_j c_{j'}
    where c_j = cos(2πkd_j/p).

    E[|D̂11(k)|²] = 4 Σ_{j=1}^F Σ_{j'=1}^F E[1_S(j)·1_S(j')] c_j c_{j'}
    """
    F = (p - 1) // 2  # number of pairs
    s = (p + 1) // 4  # number of pairs selected

    # E[1_S(j)] = s/F
    q1 = s / F

    # E[1_S(j)·1_S(j')] = s(s-1)/(F(F-1)) for j ≠ j'
    q2 = s * (s - 1) / (F * (F - 1)) if F > 1 else 0

    pairs = get_pairs(p)
    # d_j = pairs[j][0]

    results = {}
    for k in range(1, (p + 1) // 2):
        c = [math.cos(2 * math.pi * k * pairs[j][0] / p) for j in range(F)]

        # E[|D̂11(k)|²] = 4 [q1 Σ c_j² + q2 Σ_{j≠j'} c_j c_{j'}]
        # = 4 [q1 Σ c_j² + q2 ((Σ c_j)² - Σ c_j²)]
        # = 4 [(q1 - q2) Σ c_j² + q2 (Σ c_j)²]

        sum_c = sum(c)
        sum_c2 = sum(ci ** 2 for ci in c)

        E_power = 4 * ((q1 - q2) * sum_c2 + q2 * sum_c ** 2)

        # Now: Σ c_j = Σ_{j=1}^F cos(2πkd_j/p)
        # Since {d_j} ∪ {-d_j} = {1,...,p-1}:
        # 2·Σ c_j = Σ_{d=1}^{p-1} cos(2πkd/p) = Re(Σ_{d=1}^{p-1} ω^{kd}) = -1
        # So Σ c_j = -1/2 for all k ≠ 0.

        # And: Σ c_j² = Σ cos²(2πkd_j/p) = (1/2) Σ (1 + cos(4πkd_j/p))
        # = F/2 + (1/2) Σ cos(4πkd_j/p)
        # 2·Σ cos(4πkd_j/p) = Σ_{d=1}^{p-1} cos(4πkd/p)
        # = Re(Σ_{d=1}^{p-1} ω^{2kd})
        # If 2k mod p ≠ 0: = -1, so Σ cos(4πkd_j/p) = -1/2
        # So Σ c_j² = F/2 - 1/4 = (p-1)/4 - 1/4 = (p-2)/4

        results[k] = {
            "E_power": E_power,
            "sum_c": sum_c,
            "sum_c2": sum_c2,
            "theoretical_sum_c": -0.5,
            "theoretical_sum_c2": (p - 2) / 4,
        }

    return results, q1, q2


def main():
    print("=" * 90)
    print("SPECTRAL ANALYSIS TOWARD PROOF")
    print("=" * 90)

    primes = [11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83]

    print(f"\n{'p':>5s} {'F':>4s} {'s':>4s} {'E[|D̂|²]':>10s} {'theory':>10s} "
          f"{'E[max_A]':>10s} {'LP_rate':>8s}")
    print("-" * 70)

    for p in primes:
        n = (p + 1) // 2
        F = (p - 1) // 2
        s = (p + 1) // 4
        pairs = get_pairs(p)

        # Theoretical expected spectrum
        theory, q1, q2 = theoretical_spectrum_moments(p)

        # E[|D̂11(k)|²] for generic k:
        # Using Σ c_j = -1/2 and Σ c_j² = (p-2)/4:
        # E = 4[(q1-q2)(p-2)/4 + q2·1/4]
        # = (q1-q2)(p-2) + q2
        # = q1(p-2) - q2(p-2) + q2
        # = q1(p-2) - q2(p-3)

        E_generic = q1 * (p - 2) - q2 * (p - 3)

        # Monte Carlo: sample random D11, compute statistics
        rng = np.random.default_rng(42)
        num_samples = min(5000, max(500, 500000 // (p * p)))

        max_A_vals = []
        max_power_vals = []
        mean_power_vals = []

        for _ in range(num_samples):
            chosen = rng.choice(len(pairs), size=s, replace=False)
            D11 = set()
            for i in chosen:
                D11.add(pairs[i][0])
                D11.add(pairs[i][1])

            ind = np.zeros(p, dtype=np.float64)
            for j in D11:
                ind[j] = 1.0

            dft = np.fft.fft(ind)
            power = np.abs(dft) ** 2
            autocorr = np.fft.ifft(power).real
            A = np.round(autocorr).astype(int)

            max_A_d11 = max(int(A[d]) for d in D11)
            max_power = float(max(power[1:F + 1]))
            mean_power = float(np.mean(power[1:F + 1]))

            max_A_vals.append(max_A_d11)
            max_power_vals.append(max_power)
            mean_power_vals.append(mean_power)

        E_max_A = sum(max_A_vals) / len(max_A_vals)

        # LP feasibility rate (subsample)
        lp_sub = min(200, num_samples)
        lp_feasible = 0
        for trial in range(lp_sub):
            chosen = rng.choice(len(pairs), size=s, replace=False)
            D11 = set()
            for i in chosen:
                D11.add(pairs[i][0])
                D11.add(pairs[i][1])

            # Quick LP check
            ind11 = np.zeros(p, dtype=np.float64)
            for j in D11:
                ind11[j] = 1.0
            dft11 = np.fft.fft(ind11)
            power11 = np.abs(dft11) ** 2

            d11_size = len(D11)
            d12_size = (p - 1) // 2
            P0 = d11_size ** 2 + d12_size ** 2
            S_val = p * (d11_size + d12_size) - P0

            c_obj = np.zeros(F + 1)
            c_obj[F] = 1.0

            A_ub_list = []
            b_ub_list = []
            D11_set = set(D11)

            for d in sorted(D11):
                row = np.zeros(F + 1)
                for k in range(1, F + 1):
                    row[k - 1] = (2.0 / p) * math.cos(2 * math.pi * k * d / p)
                row[F] = -1.0
                A_ub_list.append(row)
                b_ub_list.append(-P0 / p)

            row_t = np.zeros(F + 1)
            row_t[F] = 1.0
            A_ub_list.append(row_t)
            b_ub_list.append(n - 2)

            D22 = set(range(1, p)) - D11_set
            for d in sorted(D22):
                row = np.zeros(F + 1)
                for k in range(1, F + 1):
                    row[k - 1] = (2.0 / p) * math.cos(2 * math.pi * k * d / p)
                A_ub_list.append(row)
                b_ub_list.append(n + 1 - P0 / p)

            A_eq = np.zeros((1, F + 1))
            for k in range(F):
                A_eq[0, k] = 2.0
            b_eq = np.array([S_val])

            bounds = [(power11[k], None) for k in range(1, F + 1)] + [(None, None)]

            res = linprog(c_obj, A_ub=np.array(A_ub_list), b_ub=np.array(b_ub_list),
                          A_eq=A_eq, b_eq=b_eq, bounds=bounds, method='highs',
                          options={'presolve': True})
            if res.success and res.x[F] <= n - 2 + 1e-6:
                lp_feasible += 1

        lp_rate = lp_feasible / lp_sub

        print(f"{p:5d} {F:4d} {s:4d} {np.mean(mean_power_vals):10.2f} {E_generic:10.2f} "
              f"{E_max_A:10.2f} {lp_rate:8.3f}")

    # Detailed analysis for key primes
    print(f"\n{'='*90}")
    print("DETAILED SPECTRAL ANALYSIS")
    print(f"{'='*90}")

    for p in [11, 23, 43, 67, 83]:
        n = (p + 1) // 2
        F = (p - 1) // 2
        s = (p + 1) // 4
        pairs = get_pairs(p)

        print(f"\n  p={p}, n={n}, F={F}, s={s}")

        # Theoretical analysis
        q1 = s / F
        q2 = s * (s - 1) / (F * (F - 1))

        # E[|D̂11(k)|²]:
        E_power = q1 * (p - 2) - q2 * (p - 3)
        print(f"    q1 = {q1:.4f}, q2 = {q2:.4f}")
        print(f"    E[|D̂11(k)|²] = {E_power:.4f}")

        # For the LP to be feasible, we need max_k |D̂11(k)|² to be not too large.
        # The LP margin depends on how flat |D̂11(k)|² is.

        # If |D̂11(k)|² = c for all k>0, what value?
        # Parseval: Σ_{k=1}^{p-1} |D̂11(k)|² = p|D11| - |D11|²
        # = p(p+1)/2 - (p+1)²/4 = (p+1)(2p - (p+1))/4 = (p+1)(p-1)/4
        # So c = (p+1)(p-1)/(4(p-1)) = (p+1)/4
        ideal_power = (p + 1) / 4
        print(f"    Ideal flat |D̂11(k)|² = (p+1)/4 = {ideal_power:.4f}")

        # MC for max |D̂11(k)|²
        rng = np.random.default_rng(42)
        num_mc = min(2000, max(300, 200000 // (p * p)))

        max_powers = []
        power_stds = []

        for _ in range(num_mc):
            chosen = rng.choice(len(pairs), size=s, replace=False)
            D11 = set()
            for i in chosen:
                D11.add(pairs[i][0])
                D11.add(pairs[i][1])

            spec = d11_spectrum(p, D11)
            pk = spec[1:F + 1]
            max_powers.append(float(max(pk)))
            power_stds.append(float(np.std(pk)))

        print(f"    MC ({num_mc} samples):")
        print(f"      E[max |D̂11(k)|²] = {np.mean(max_powers):.2f}")
        print(f"      max |D̂11(k)|² / p = {np.mean(max_powers)/p:.4f}")
        print(f"      E[std |D̂11(k)|²] = {np.mean(power_stds):.2f}")

        # Key ratio: max_power / ideal_power
        ratio = np.mean(max_powers) / ideal_power
        print(f"      max/ideal ratio: {ratio:.2f}")

        # What max_power does the LP need?
        # If max |D̂11(k)|² = M, then to keep |D̂12(k)|² = P(k) - |D̂11(k)|² ≥ 0,
        # we need P(k) ≥ M at the peak. But P(k) also needs to satisfy the
        # Parseval sum and the D11/D22 constraints.
        # Roughly: if M is too large, the spectrum can't be controlled.

    # Summary
    print(f"\n{'='*90}")
    print("SCALING SUMMARY")
    print(f"{'='*90}")
    print("""
Key findings:
1. E[|D̂11(k)|²] ≈ (p+1)/4 (the ideal flat value), confirming random D11 is
   "on average" good spectrally.

2. max |D̂11(k)|² grows with p. If it grows as O(p), then LP margin remains
   bounded. If it grows as O(p log p), LP feasibility might fail.

3. The LP feasibility rate decreases with p, but remains positive for all
   tested primes. This suggests LP-feasible D11 exists for all p.

4. The LP feasibility is predicted by:
   - Low max A(d) at D11 positions (controlled autocorrelation)
   - Low max |D̂11(k)|² (controlled spectral peak)
   - Low power spectrum std (flat spectrum)

TOWARD A PROOF:
  If we can show: for some D11, max |D̂11(k)|² ≤ (1+ε)(p+1)/4 + O(√p),
  then the LP analysis shows the constraints are satisfiable.
  Constructing such D11 is related to character sum bounds (Weil/Deligne).
""")


if __name__ == "__main__":
    main()
