"""
Fourier-analytic reformulation of the Ramsey book graph constraints.

Key insight: Express Delta(D11,D11,d) and Delta(D12,D12,d) in terms of
Fourier coefficients, then characterize the constraint conditions spectrally.

For a set S ⊂ Z_m, define hat{S}(k) = sum_{s in S} omega^{ks} where omega = e^{2pi i/m}.

Then Delta(S, S, d, m) = (1/m) * sum_k |hat{S}(k)|^2 * omega^{-kd}

This is the inverse Fourier transform of |hat{S}(k)|^2.

The V1V1 constraint: A(d) + B(d) <= n-2 for d in D11, where:
  A(d) = Delta(D11, D11, d) = (1/m) sum_k |hat{D11}(k)|^2 omega^{-kd}
  B(d) = Delta(D12, D12, d) = (1/m) sum_k |hat{D12}(k)|^2 omega^{-kd}
"""

import sys, os, json, cmath, math
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import Delta, Sigma

KNOWN = {
    6: {"D11": {3, 5, 6, 8}, "D12": {0, 1, 4, 6, 7}, "D22": {1, 2, 4, 7, 9, 10}},
    8: {"D11": {3, 6, 7, 8, 9, 12}, "D12": {0, 1, 4, 6, 8, 9, 13}, "D22": {1, 2, 4, 5, 10, 11, 13, 14}},
    10: {"D11": {4, 5, 7, 9, 10, 12, 14, 15}, "D12": {0, 1, 2, 6, 7, 10, 11, 13, 17}, "D22": {1, 2, 3, 6, 8, 11, 13, 16, 17, 18}},
    12: {"D11": {5, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17, 18}, "D12": {0, 1, 2, 6, 10, 13, 14, 16, 18, 20, 21}, "D22": {1, 2, 3, 4, 10, 13, 19, 20, 21, 22}},
    14: {"D11": {5, 7, 8, 9, 10, 11, 13, 14, 16, 17, 18, 19, 20, 22}, "D12": {0, 1, 2, 7, 8, 10, 13, 14, 17, 18, 21, 23, 25}, "D22": {1, 2, 3, 4, 6, 12, 15, 21, 23, 24, 25, 26}},
    16: {"D11": {6, 7, 8, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 23, 24, 25}, "D12": {0, 1, 2, 3, 8, 11, 12, 13, 15, 18, 20, 21, 24, 27, 29}, "D22": {1, 2, 3, 4, 5, 9, 13, 18, 22, 26, 27, 28, 29, 30}},
    18: {"D11": {6, 8, 9, 10, 11, 13, 14, 15, 17, 18, 20, 21, 22, 24, 25, 26, 27, 29}, "D12": {0, 1, 2, 3, 8, 9, 10, 12, 14, 15, 19, 20, 22, 24, 25, 28, 32}, "D22": {1, 2, 3, 4, 5, 7, 12, 16, 19, 23, 28, 30, 31, 32, 33, 34}},
    20: {"D11": {7, 8, 9, 10, 12, 13, 14, 16, 18, 19, 20, 21, 23, 25, 26, 27, 29, 30, 31, 32}, "D12": {0, 1, 2, 5, 6, 8, 12, 14, 15, 17, 18, 22, 23, 25, 27, 28, 33, 36, 37}, "D22": {1, 2, 3, 4, 5, 6, 11, 15, 17, 22, 24, 28, 33, 34, 35, 36, 37, 38}},
}

n22_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "solution_n22.json")
if os.path.exists(n22_path):
    with open(n22_path) as f:
        data = json.load(f)
    KNOWN[22] = {"D11": set(data["parameters"]["D11"]), "D12": set(data["parameters"]["D12"]), "D22": set(data["parameters"]["D22"])}


def fourier_coeff(S, k, m):
    """hat{S}(k) = sum_{s in S} omega^{ks}."""
    omega = cmath.exp(2j * cmath.pi / m)
    return sum(omega ** (k * s) for s in S)


def power_spectrum(S, m):
    """Compute |hat{S}(k)|^2 for k = 0, ..., m-1."""
    return [abs(fourier_coeff(S, k, m)) ** 2 for k in range(m)]


print("FOURIER ANALYSIS OF RAMSEY CONSTRAINTS")
print("=" * 80)

for n in sorted(KNOWN.keys()):
    m = 2 * n - 1
    D11 = KNOWN[n]["D11"]
    D12 = KNOWN[n]["D12"]

    ps_d11 = power_spectrum(D11, m)
    ps_d12 = power_spectrum(D12, m)

    print(f"\nn={n}, m={m}")

    # Verify: A(d) = (1/m) sum_k |hat{D11}(k)|^2 omega^{-kd}
    omega = cmath.exp(2j * cmath.pi / m)
    for d in [1, 2, m // 2]:
        A_direct = Delta(D11, D11, d, m)
        A_fourier = sum(ps_d11[k] * omega ** (-k * d) for k in range(m)) / m
        assert abs(A_direct - A_fourier.real) < 1e-6, f"Fourier mismatch at d={d}"

    # Since D11 is symmetric: hat{D11}(k) is REAL for all k.
    # Because hat{D11}(k) = sum omega^{ks} and for every s in D11, -s is also in D11:
    # hat{D11}(k) = sum_{s in D11} (omega^{ks} + omega^{-ks})/2 * 2 = sum 2*cos(2*pi*k*s/m)/2...
    # Actually: hat{D11}(k) = sum omega^{ks}. Since D11 = -D11:
    # hat{D11}(k) = sum_{s in D11} omega^{ks} = sum_{s in D11} omega^{k(-s)} = conj(hat{D11}(k))
    # So hat{D11}(k) is real. ✓

    is_real = all(abs(fourier_coeff(D11, k, m).imag) < 1e-8 for k in range(m))
    print(f"  hat{{D11}}(k) is real: {is_real}")

    # The constraint A(d) + B(d) <= n-2 for d in D11 becomes:
    # (1/m) sum_k (|hat{D11}(k)|^2 + |hat{D12}(k)|^2) omega^{-kd} <= n-2
    #
    # The average of A(d)+B(d) over all d != 0:
    # (1/(m-1)) sum_{d=1}^{m-1} [A(d) + B(d)]
    # = (1/(m-1)) * [(|D11|^2 + |D12|^2 - |D11| - |D12|) + ... hmm]
    # Actually sum_{d=1}^{m-1} A(d) = sum_{d=1}^{m-1} Delta(D11,D11,d)
    # = |{(a,b) in D11 x D11 : a-b != 0}| = |D11|^2 - |D11| = |D11|(|D11|-1)
    # Similarly sum B(d) = |D12|^2 - |D12| for d = 1..m-1?
    # No: sum_{d=0}^{m-1} Delta(S,S,d) = |S|^2 (all pairs). Delta(S,S,0) = |S|.
    # So sum_{d=1}^{m-1} Delta(S,S,d) = |S|^2 - |S|.
    # But B(d) = Delta(D12,D12,d) and d ranges over 1..m-1.
    # sum_{d=1}^{m-1} B(d) = |D12|^2 - |D12|.

    avg_AB = (len(D11) * (len(D11) - 1) + len(D12) * (len(D12) - 1)) / (m - 1)
    print(f"  Average A(d)+B(d) over d=1..m-1: {avg_AB:.4f}")
    print(f"  Red threshold: {n - 2}")
    print(f"  Average excess over threshold: {avg_AB - (n - 2):.4f}")

    # For the constraint to be satisfiable, we need A(d)+B(d) <= n-2 for d in D11,
    # but A(d)+B(d) can exceed n-2 for d not in D11 (blue edges in V1V1).
    #
    # The total "budget" for the D11 sum:
    # sum_{d in D11} (A(d)+B(d)) + sum_{d not in D11, d>0} (A(d)+B(d)) = total
    # We need sum_{d in D11} (A(d)+B(d)) <= |D11| * (n-2)
    #
    # Similarly, for blue: the blue common neighbors at d (d not in D11) are
    # (N-2) - 2*d1 + (A(d)+B(d)), which must be <= n-1.
    # So A(d)+B(d) <= n-1 + 2*d1 - (N-2) = n-1 + 2*(|D11|+|D12|) - (4n-4)
    # For the standard sizes |D11|=n, |D12|=n-1:
    # = n-1 + 2*(2n-1) - 4n + 4 = n-1 + 4n - 2 - 4n + 4 = n + 1
    # So blue constraint: A(d)+B(d) <= n+1 for d not in D11.

    d1 = len(D11) + len(D12)
    N = 4 * n - 2
    blue_max = n - 1 + 2 * d1 - (N - 2)
    print(f"  Blue upper bound on A(d)+B(d): {blue_max}")

    # The Fourier picture: define P(k) = |hat{D11}(k)|^2 + |hat{D12}(k)|^2.
    # Then A(d)+B(d) = (1/m) sum_k P(k) omega^{-kd}.
    # The k=0 term: P(0) = |D11|^2 + |D12|^2.
    # So A(d)+B(d) = (|D11|^2 + |D12|^2)/m + (1/m) sum_{k>0} P(k) omega^{-kd}
    # The first term is the average contribution. The second is the fluctuation.

    P_spectrum = [ps_d11[k] + ps_d12[k] for k in range(m)]
    P0 = P_spectrum[0]
    avg_from_P0 = P0 / m
    print(f"  P(0)/m = (|D11|^2+|D12|^2)/m = {avg_from_P0:.4f}")

    # The fluctuation at d: F(d) = (1/m) sum_{k=1}^{m-1} P(k) omega^{-kd}
    # = A(d)+B(d) - P0/m
    # We need F(d) <= (n-2) - P0/m for d in D11 (red constraint)
    # and F(d) <= blue_max - P0/m for d not in D11 (blue constraint)

    red_fluctuation_bound = (n - 2) - avg_from_P0
    blue_fluctuation_bound = blue_max - avg_from_P0
    print(f"  Red needs fluctuation F(d) <= {red_fluctuation_bound:.4f}")
    print(f"  Blue needs fluctuation F(d) <= {blue_fluctuation_bound:.4f}")

    # Compute actual fluctuations
    flucts = {}
    for d in range(1, m):
        ab = Delta(D11, D11, d, m) + Delta(D12, D12, d, m)
        flucts[d] = ab - avg_from_P0

    red_flucts = [flucts[d] for d in D11]
    blue_flucts = [flucts[d] for d in range(1, m) if d not in D11]
    print(f"  Red fluctuations: min={min(red_flucts):.4f}, max={max(red_flucts):.4f} (need <= {red_fluctuation_bound:.4f})")
    print(f"  Blue fluctuations: min={min(blue_flucts):.4f}, max={max(blue_flucts):.4f} (need <= {blue_fluctuation_bound:.4f})")

    # Key question: what fraction of the Fourier energy is in the "useful" direction?
    # P(k) for k>0: energy that creates fluctuation
    total_energy = sum(P_spectrum[k] for k in range(1, m))
    print(f"  Total P(k>0) energy: {total_energy:.4f}")
    print(f"  Energy per frequency: {total_energy / (m - 1):.4f}")

    # The D11 indicator as a Fourier series: 1_{D11}(d) = (1/m) sum_k hat{D11}(k) omega^{-kd}
    # For the constraint to work, we need 1_{D11}(d) * F(d) to be "correlated negatively"
    # i.e., F(d) is small where 1_{D11}(d) = 1

    # Cross-correlation: sum_{d=1}^{m-1} 1_{D11}(d) * (A(d)+B(d))
    # = sum_{d in D11} (A(d)+B(d))
    cross = sum(Delta(D11, D11, d, m) + Delta(D12, D12, d, m) for d in D11)
    total_ab = sum(Delta(D11, D11, d, m) + Delta(D12, D12, d, m) for d in range(1, m))
    print(f"  sum_{{d in D11}} (A+B) = {cross}, vs uniform prediction {total_ab * len(D11) / (m-1):.1f}")
    print(f"  Negative correlation: {cross < total_ab * len(D11) / (m - 1)}")

print("\n" + "=" * 80)
print("FOURIER CONSTRAINT SUMMARY")
print("=" * 80)
print("""
The V1V1 constraint in Fourier space:

Define P(k) = |hat{D11}(k)|^2 + |hat{D12}(k)|^2 for k = 0,...,m-1.

Then the common neighbor count at difference d is:
  lambda(d) = (1/m) * sum_k P(k) omega^{-kd}

This must be:
  <= n-2 for d in D11 (red edges)
  <= n+1 for d not in D11 (blue edges, after inclusion-exclusion)

The average is P(0)/m = (|D11|^2 + |D12|^2)/m ~ n + (n-1)^2/m ~ 2n-1 as m = 2n-1.
Wait: (n^2 + (n-1)^2)/(2n-1) = (2n^2 - 2n + 1)/(2n-1) = n - 1/(2n-1) + ... ~ n.

Actually for |D11|=n, |D12|=n-1, m=2n-1:
P(0)/m = (n^2 + (n-1)^2)/(2n-1) = (2n^2 - 2n + 1)/(2n-1)

For large n: ~ n - 0.5 + O(1/n). Close to n-1, so avg is just below n.
Red threshold is n-2, so we need F(d) ~ -1 for red edges.
Blue threshold is n+1, so F(d) < +2 for blue edges -- very generous.

The key requirement: fluctuations of ~1 unit below average on a set of size n out of m-1 ~ 2n.
This is about half the elements being below average by ~1 -- very achievable!

For a PROOF: need to show that for any symmetric D11 with |D11| = (m+1)/2 in Z_p,
there exists D12 with |D12| = (m+1)/2 - 1 such that:
  max_{d in D11} [A(d) + B(d)] <= (m+1)/2 - 2
  max_{d not in D11} [A(d) + B(d)] <= (m+1)/2 + 1
""")

# Now let's check: for each known case, what is the Fourier decomposition of
# the "excess" function A(d)+B(d) - threshold, restricted to D11?
print("\nPROOF DIRECTION: Fourier characterization of valid (D11, D12) pairs")
print("=" * 80)

for n in sorted(KNOWN.keys()):
    m = 2 * n - 1
    D11 = KNOWN[n]["D11"]
    D12 = KNOWN[n]["D12"]

    # The constraint sum: what is sum_{d in D11} (A(d)+B(d)) vs |D11|*(n-2)?
    total_red = sum(Delta(D11, D11, d, m) + Delta(D12, D12, d, m) for d in D11)
    budget = len(D11) * (n - 2)
    slack = budget - total_red

    total_blue = sum(Delta(D11, D11, d, m) + Delta(D12, D12, d, m) for d in range(1, m) if d not in D11)
    d1 = len(D11) + len(D12)
    blue_threshold = n - 1 + 2 * d1 - (4 * n - 4)
    blue_budget = len(set(range(1, m)) - D11) * blue_threshold
    blue_slack = blue_budget - total_blue

    print(f"n={n}: red_sum={total_red}, red_budget={budget}, red_slack={slack}")
    print(f"       blue_sum={total_blue}, blue_budget={blue_budget}, blue_slack={blue_slack}")
