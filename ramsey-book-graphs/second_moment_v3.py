#!/usr/bin/env python3
"""
Second moment computation for Paley-Zygmund.

For a FIXED D11, N = #{valid D12} is deterministic.
The second moment method applies to RANDOM D11 (or to the "random variable"
interpretation where we think of choosing D12 at random).

The correct Paley-Zygmund setup:
  Let f(S) = 1[{0} union S is valid for D11], where S is a random k-subset.
  E[f] = p_0 = N/C(p-1,k)
  E[f(S)*f(S')] for independent S, S' would be p_0^2.
  But S, S' are NOT independent -- they are drawn from the same ground set.

  Actually for the second moment method:
  N = sum_{S: |S|=k} f(S)
  E_S[N] = E_S[sum f] is the EXPECTED number of valid S when we count all S.
  Wait, N is just the total count. We want Pr[N > 0].

  The right framing: we want to show N > 0. We know:
    N = sum_{|S|=k} f(S)
    By first moment: N >= 0 (trivially).
    By Paley-Zygmund applied to a random variable X:
      If X is the count in a RANDOM experiment, Pr[X > 0] >= E[X]^2/E[X^2].

  But N is not random -- D11 is fixed.

  The CORRECT approach: make D11 random, and consider
    N(D11) = #{valid D12 for D11}
  Then:
    E[N] = E_{D11}[N(D11)] (average over D11)
    E[N^2] = E_{D11}[N(D11)^2]
    Pr[N > 0] >= E[N]^2 / E[N^2]  (Paley-Zygmund)

  This shows: some D11 has N(D11) > 0.

  Alternatively: fix D11 and use the "random D12" version:
    Choose a random k-subset S. Let X = 1[S valid].
    E[X] = p_0 = N/C(p-1,k).
    We want to show E[X] > 0, i.e., N > 0.
    But E[X] > 0 iff N > 0 -- no need for second moment!

  So the second moment method is only useful when we randomize D11.
  OR: when we can't compute E[X] = p_0 directly (which we can't for large p).

  The actual proof strategy should be:
  1. For a specific A-flat D11, bound E[N(D11)] from below.
  2. This is done by product-of-marginals (which we've done: ~ 2^{(p-1)/2}).
  3. The "correlation loss" c_0 is what makes this fail.
  4. The second moment method can bound Var[N(D11)] to show N(D11) > 0 w.h.p.

  Wait, N(D11) for fixed D11 is deterministic! There's no variance.

  I think the correct interpretation is: we DON'T fix D11 first. Instead:
  - Choose D11 randomly (symmetric, size n).
  - Choose D12 randomly (with 0, size (p-1)/2).
  - Let X(D11, D12) = 1[valid construction].
  - N = sum_{D12} sum_{D11} X(D11, D12)
  - Or N_fixed = sum_{D12} X(D11, D12) for fixed D11.

  Actually, let me re-read the task description more carefully...

  "For fixed A-flat D11, let N = #{valid D12}"
  "E[N] = C(p-1,k) × Pr[valid] grows exponentially"
  "E[N^2] = Σ_{D12,D12'} Pr[both valid], grouped by overlap s"

  OK so they're computing E[N^2] = sum_{D12,D12'} Pr[both valid] where
  the "probability" is over... the randomness in D12?

  I think the setup is: N is a random variable because D12 is random.
  N = sum_{D12 in all k-subsets} 1[D12 valid].
  This is a FIXED number (deterministic).

  BUT: if we think of "drawing D12 at random" and asking "will a valid D12 exist?",
  we can use the second moment method on the COUNT.

  Let X_S = 1[{0} union S is valid]. N = sum_S X_S.
  E[N] = sum_S E[X_S] = C(p-1,k) * Pr[random S is valid].

  Hmm E[N] = N since X_S is deterministic (0 or 1, not random).
  The expectation is trivial.

  UNLESS the randomness is in D11 as well. Let me try:

  Choose D11 uniformly at random among A-flat D11.
  Then N(D11) is random.
  E_{D11}[N(D11)] = average over A-flat D11 of #{valid D12}.
  E_{D11}[N(D11)^2] = average over A-flat D11 of #{valid D12}^2.

  Paley-Zygmund: Pr_{D11}[N(D11) > 0] >= E[N]^2 / E[N^2].

  This is indeed the right approach when D11 is random.

Alternatively, fix D11 and consider the PAIR (D12, D12'):
  P(s) = Pr_{D12, D12'}[D12 valid AND D12' valid | |D12 ∩ D12'| = s+1]
  where D12, D12' are two INDEPENDENT random k-subsets conditioned on overlap s.

  Then: E[N(N-1)] = sum_{S != S'} Pr[both valid]
  and if we can show E[N^2] / E[N]^2 = O(poly(p)), done.

  But N is deterministic so E[N^2] = N^2 and E[N]^2 = N^2. Ratio = 1. Trivial.

OK I think the only meaningful approach is to randomize D11 as well.
Let me compute this for small p.
"""

from itertools import combinations
from collections import defaultdict
from math import comb, log2, sqrt


def compute_A(D11_set, p):
    A = [0]*p
    for a in D11_set:
        for b in D11_set:
            if a != b:
                A[(a - b) % p] += 1
    return A


def compute_B(S_set, p):
    B = [0]*p
    for a in S_set:
        for b in S_set:
            if a != b:
                B[(a - b) % p] += 1
    return B


def check_valid(B, T, p):
    for d in range(1, p):
        if B[d] > T[d]:
            return False
    return True


def get_thresholds(D11_set, A, p):
    T = [0]*p
    for d in range(1, p):
        if d in D11_set:
            T[d] = (p - 3) // 2 - A[d]
        else:
            T[d] = (p + 3) // 2 - A[d]
    return T


def second_moment_random_d11(p):
    """Second moment with D11 randomized."""
    n = (p + 1) // 2
    k = (p - 3) // 2

    # Generate all symmetric D11
    pairs = []
    seen = set()
    for d in range(1, p):
        if d not in seen:
            pairs.append((d, p - d))
            seen.add(d)
            seen.add(p - d)

    all_d11s = []
    for combo in combinations(range(len(pairs)), n // 2):
        D11 = []
        for idx in combo:
            d1, d2 = pairs[idx]
            D11.extend([d1, d2])
        if len(D11) == n:
            all_d11s.append(sorted(D11))

    print(f"p = {p}, n = {n}, k = {k}")
    print(f"Total symmetric D11: {len(all_d11s)}")

    # For each D11, count valid D12
    elements = list(range(1, p))
    all_S = list(combinations(elements, k))
    total_S = len(all_S)

    # Precompute B for all S
    S_B = {}
    for S in all_S:
        S_set = frozenset(S)
        B = compute_B(S_set, p)
        S_B[S_set] = B

    d11_counts = {}
    for D11 in all_d11s:
        D11_set = set(D11)
        A = compute_A(D11_set, p)
        T = get_thresholds(D11_set, A, p)
        count = 0
        for S_set, B in S_B.items():
            if check_valid(B, T, p):
                count += 1
        d11_counts[tuple(D11)] = count
        max_A = max(A[d] for d in range(1, p))

    # Statistics
    counts = list(d11_counts.values())
    E_N = sum(counts) / len(counts)
    E_N2 = sum(c**2 for c in counts) / len(counts)
    ratio = E_N2 / E_N**2 if E_N > 0 else float('inf')

    print(f"\nAll D11 statistics:")
    print(f"  E[N] = {E_N:.2f}")
    print(f"  E[N^2] = {E_N2:.2f}")
    print(f"  E[N^2]/E[N]^2 = {ratio:.4f}")
    print(f"  Paley-Zygmund: Pr[N>0] >= {1/ratio:.6f}" if ratio < float('inf') else "  N=0 everywhere")
    print(f"  Pr[N>0] (exact) = {sum(1 for c in counts if c > 0)/len(counts):.4f}")
    print(f"  N values: min={min(counts)}, max={max(counts)}, median={sorted(counts)[len(counts)//2]}")

    # Also: restrict to A-flat D11 (those with min max_A)
    d11_maxA = {}
    for D11 in all_d11s:
        D11_set = set(D11)
        A = compute_A(D11_set, p)
        max_A = max(A[d] for d in range(1, p))
        d11_maxA[tuple(D11)] = max_A

    min_maxA = min(d11_maxA.values())
    flat_d11s = [D11 for D11 in all_d11s if d11_maxA[tuple(D11)] == min_maxA]
    flat_counts = [d11_counts[tuple(D11)] for D11 in flat_d11s]

    if flat_counts:
        E_N_flat = sum(flat_counts) / len(flat_counts)
        E_N2_flat = sum(c**2 for c in flat_counts) / len(flat_counts)
        ratio_flat = E_N2_flat / E_N_flat**2 if E_N_flat > 0 else float('inf')

        print(f"\nA-flat D11 (max_A = {min_maxA}): {len(flat_d11s)} total")
        print(f"  E[N] = {E_N_flat:.2f}")
        print(f"  E[N^2] = {E_N2_flat:.2f}")
        print(f"  E[N^2]/E[N]^2 = {ratio_flat:.4f}")
        print(f"  Paley-Zygmund: Pr[N>0] >= {1/ratio_flat:.6f}" if ratio_flat < float('inf') else "  N=0 everywhere")
        print(f"  N values: {sorted(flat_counts)}")

    # Distribution of N
    print(f"\nN distribution (all D11):")
    hist = defaultdict(int)
    for c in counts:
        hist[c] += 1
    for n_val in sorted(hist.keys()):
        print(f"  N = {n_val:4d}: {hist[n_val]:4d} D11s ({hist[n_val]/len(counts)*100:.1f}%)")


def main():
    for p in [11, 19]:
        if p % 4 != 3:
            continue
        print("=" * 70)
        second_moment_random_d11(p)
        print()


if __name__ == '__main__':
    main()
