#!/usr/bin/env python3
"""Verify the constant z-sum claim from the second moment proof.

Claim: For every symmetric D11 of size n = (p+1)/2:
  sum_{i=1}^R z_i = (p-7)/(4*sigma)

where z_i = (tau(d_i) - mu_B) / sigma are the constraint z-scores,
tau(d) = T(d) - A(d), and sigma = sqrt(Var[B(d)]).
"""
from itertools import combinations
from math import sqrt

def verify_constant_zsum(p):
    n = (p + 1) // 2
    k = (p - 3) // 2
    mu_B = (p - 3) / 4
    var_B = (p - 3) * (p + 1) / (16 * (p - 2))
    sigma = sqrt(var_B)

    # Generate all complementary pairs
    pairs = []
    for d in range(1, (p + 1) // 2):
        pairs.append((d, p - d))

    # Number of pairs to select for D11
    num_pairs = n // 2  # n is even since p = 3 mod 4 => n = (p+1)/2 is even

    # Enumerate all symmetric D11
    z_sums = []
    s2_values = []

    for chosen_pairs in combinations(range(len(pairs)), num_pairs):
        D11 = set()
        for idx in chosen_pairs:
            d, pd = pairs[idx]
            D11.add(d)
            D11.add(pd)
        D22 = set(range(1, p)) - D11

        # Compute A(d) for all nonzero d
        D11_list = sorted(D11)
        A = {}
        for d in range(1, p):
            count = 0
            for a in D11_list:
                b = (a - d) % p
                if b != 0 and b in D11:
                    count += 1
            A[d] = count

        # Compute z-scores for representatives
        z_scores = []
        reps = [pairs[idx][0] for idx in chosen_pairs]  # D11 reps
        d22_reps = [pairs[idx][0] for idx in range(len(pairs)) if idx not in chosen_pairs]  # D22 reps

        for d in reps:
            tau = (p - 3) / 2 - A[d]
            z = (tau - mu_B) / sigma
            z_scores.append(z)

        for d in d22_reps:
            tau = (p + 3) / 2 - A[d]
            z = (tau - mu_B) / sigma
            z_scores.append(z)

        z_sum = sum(z_scores)
        s2 = sum(z**2 for z in z_scores)
        z_sums.append(z_sum)
        s2_values.append(s2)

    # Theoretical prediction
    predicted = (p - 7) / (4 * sigma)

    print(f"p={p}: sigma={sigma:.4f}")
    print(f"  Predicted sum(z) = (p-7)/(4*sigma) = {predicted:.6f}")
    print(f"  Actual sum(z) range: [{min(z_sums):.6f}, {max(z_sums):.6f}]")
    print(f"  All equal: {abs(max(z_sums) - min(z_sums)) < 1e-10}")
    print(f"  S2 range: [{min(s2_values):.2f}, {max(s2_values):.2f}]")
    print(f"  #{len(z_sums)} symmetric D11 checked")
    print()

for p in [11, 19, 23]:
    verify_constant_zsum(p)
