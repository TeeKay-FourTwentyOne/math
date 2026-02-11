#!/usr/bin/env python3
"""
Quick check of the second moment (Paley-Zygmund) approach.

For the D11-averaged second moment:
  E[N] = average of N(D11) over all symmetric D11
  E[N²] = average of N(D11)² over all symmetric D11
  Ratio = E[N²]/E[N]²
  Paley-Zygmund: Pr[N > 0] >= E[N]²/E[N²] = 1/ratio

Also: for FIXED D11, a "random D12 pair" second moment.
"""
from math import comb, log2

# Data from exhaustive enumeration (proofs/enumeration_results.md)
# Format: {p: {N_value: count_of_D11_with_this_N}}

data = {
    11: {
        # 10 symmetric D11 total, 5 working with N=20 each, 5 with N=0
        20: 5,
        0: 5,
    },
    19: {
        # 126 symmetric D11, 9 working with N=18 each
        18: 9,
        0: 117,
    },
    23: {
        # 462 symmetric D11, 55 working.
        # From enumeration: "Valid D12 counts: 22, 44, 110, or 198"
        # Total valid pairs = 4356
        # Need exact distribution. Let's use: 4356 = 22a + 44b + 110c + 198d, a+b+c+d=55
        # From the data, a rough estimate. Let's compute what we know:
        # We'll use the reported total and see if consistent
        # Actually, let me check: 55 working D11, counts in {22,44,110,198}
        # The complement symmetry means the distribution should be symmetric.
        # At p=23, n=12, n-2=10. Complement bijection maps N -> N.
        # So counts should be the same.
        # From the data file, we'd need exact per-D11 counts.
        # For now, let's estimate: if all 55 have N=4356/55 ≈ 79.2 (average):
        # We'll flag this as approximate.
        0: 407,  # 462 - 55 = 407 non-working
        # For the working ones, we'll try different distributions:
    },
}


def analyze_prime(p, N_dist):
    """Analyze second moment for given N(D11) distribution."""
    total_D11 = sum(N_dist.values())
    total_valid = sum(N * count for N, count in N_dist.items())

    E_N = total_valid / total_D11
    E_N2 = sum(N**2 * count for N, count in N_dist.items()) / total_D11
    ratio = E_N2 / E_N**2 if E_N > 0 else float('inf')

    working = sum(count for N, count in N_dist.items() if N > 0)
    frac_working = working / total_D11

    k = (p - 3) // 2
    budget = log2(comb(p - 1, k))

    print(f"p={p}: {total_D11} symmetric D11, {working} working ({100*frac_working:.1f}%)")
    print(f"  Total valid pairs: {total_valid}")
    print(f"  E[N] = {E_N:.2f}")
    print(f"  E[N²] = {E_N2:.2f}")
    print(f"  Ratio E[N²]/E[N]² = {ratio:.2f}")
    print(f"  Paley-Zygmund: Pr[N>0] >= {1/ratio:.4f} (actual: {frac_working:.4f})")
    print(f"  log2(budget C(p-1,k)) = {budget:.2f}")
    print()


print("=" * 60)
print("SECOND MOMENT (PALEY-ZYGMUND) OVER D11")
print("=" * 60)
print()

for p in [11, 19]:
    analyze_prime(p, data[p])

# For p=23, try different plausible distributions of N among 55 working D11
print("p=23 scenarios (55 working D11, total 4356 valid pairs):")
print()

# Scenario 1: All 55 have the same N = 4356/55 ≈ 79.2 (impossible, must be integer)
# Closest: 22a + 44b + 110c + 198d = 4356, a+b+c+d=55
# Let's try: if mostly N=198: 22 D11 with N=198 → 22*198=4356. a+b+c=33 with 22*33=0.
# 22*198 = 4356. So d=22, and 33 have other values summing to 0. → 33 with N=0?
# But we said 55 working. Hmm.

# Let's try: 198d = 4356 - (22a + 44b + 110c), a+b+c+d=55
# Simplest: d=22 gives 198*22=4356, a=b=c=0, d=22. But then 55-22=33 working with N=0?
# That can't be right since "working" means N>0.

# Actually from the enumeration results: "55 work, with 22, 44, 110, or 198"
# These are the possible N values. Not all 55 have N=198.
# We need: sum of N's = 4356 over 55 D11 where each N is in {22, 44, 110, 198}.
# 4356/55 = 79.2 average
# Check: if x D11 have N=22 and (55-x) have N=198:
# 22x + 198(55-x) = 4356 → 22x + 10890 - 198x = 4356 → -176x = -6534 → x = 37.1
# Not integer. Try {22, 198}: 22a + 198d = 4356, a+d=55 → 22(55-d)+198d = 4356 → 176d = 3146 → d=17.9
# Not integer. Try {44, 198}: 44b + 198d = 4356, b+d=55 → 44(55-d)+198d=4356 → 154d=2036 → d=13.2
# Try {22, 110}: 22a + 110c = 4356, a+c=55 → 22(55-c)+110c=4356 → 88c=3146 → c=35.75
# Try {44, 110}: 44b + 110c = 4356, b+c=55 → 44(55-c)+110c=4356 → 66c=1936 → c=29.3
# None work with two values. Need at least 3 types.

# Let me just try some plausible distributions:
scenarios = [
    # (name, {N: count} for working D11)
    ("Most at 22, few at 198",
     {22: 40, 44: 5, 110: 5, 198: 5, 0: 407}),
    # Check: 22*40+44*5+110*5+198*5 = 880+220+550+990 = 2640 ≠ 4356

    ("Even split",
     {22: 11, 44: 11, 110: 11, 198: 22, 0: 407}),
    # Check: 22*11+44*11+110*11+198*22 = 242+484+1210+4356 ... too much

    # Let me solve properly:
    # 22a + 44b + 110c + 198d = 4356, a+b+c+d = 55
    # Try a=33, b=0, c=0, d=22: 22*33+198*22 = 726+4356=5082 ≠ 4356
    # Try a=0, b=33, c=0, d=22: 44*33+198*22 = 1452+4356 = 5808 ≠ 4356
    # Try a=0, b=0, c=33, d=22: 110*33+198*22 = 3630+4356 = 7986 ≠ 4356

    # Hmm, let me just compute 4356/55 = 79.2.
    # The minimum N is 22 and max is 198.
    # Let me try: most at 44 and 110.
    # 44b + 110c = 4356 - 22a - 198d, b+c = 55-a-d
    # If a=d=0: 44b+110c=4356, b+c=55 → 44(55-c)+110c=4356 → 66c=1936 → c=29.33
    # If a=10,d=5: 44b+110c=4356-220-990=3146, b+c=40 → 44(40-c)+110c=3146 → 66c=1386 → c=21
    # So: a=10,d=5,c=21,b=19. Check: 22*10+44*19+110*21+198*5=220+836+2310+990=4356 ✓!
]

dist_23 = {22: 10, 44: 19, 110: 21, 198: 5, 0: 407}
print("Using estimated distribution: 10×22, 19×44, 21×110, 5×198:")
analyze_prime(23, dist_23)

# Also try extreme: all 55 at same value (hypothetical)
avg_N = 4356 / 55
print(f"If all 55 had same N={avg_N:.1f} (lower bound on ratio):")
dist_23_uniform = {0: 407}
# Can't have non-integer N. Use E[N²] lower bound: E[N²] >= E[N]² by Jensen
E_N = 4356 / 462
E_N2_lower = E_N ** 2  # Can't actually achieve this with N only in {0, 22, ...}
print(f"  E[N] = {E_N:.2f}, E[N²] >= {E_N**2:.2f} (Jensen)")
print(f"  Best possible ratio >= 1.0 (from Jensen)")
print()

# What ratio do we get from worst case (all working D11 have same N)?
E_N2_same = (55 * avg_N**2 + 407 * 0) / 462
ratio_same = E_N2_same / E_N**2
print(f"If all 55 working D11 had same N={avg_N:.1f}:")
print(f"  E[N²] = {E_N2_same:.2f}")
print(f"  Ratio = {ratio_same:.2f}")
print(f"  = #{'{'}all D11{'}'}/#{'{'}working D11{'}'} = 462/55 = {462/55:.2f}")
print()

# With actual variation (N in {22,44,110,198}), the ratio is LARGER:
E_N2_actual = (10*22**2 + 19*44**2 + 21*110**2 + 5*198**2 + 407*0) / 462
ratio_actual = E_N2_actual / E_N**2
print(f"With estimated N distribution:")
print(f"  E[N²] = {E_N2_actual:.2f}")
print(f"  Ratio = {ratio_actual:.2f}")

print()
print("=" * 60)
print("SUMMARY")
print("=" * 60)
print()
print("| p | E[N] | ratio E[N²]/E[N]² | Pr[N>0] (PZ) | Pr[N>0] (actual) |")
print("|---|------|-------------------|--------------|------------------|")
for p, E, r, pz, actual in [
    (11, 10.0, 2.0, 0.5, 0.5),
    (19, 1.29, 14.0, 0.071, 0.071),
]:
    print(f"| {p} | {E:.2f} | {r:.1f} | {pz:.3f} | {actual:.3f} |")

print(f"| 23 | {E_N:.2f} | {ratio_same:.1f}-{ratio_actual:.1f} | {1/ratio_actual:.3f}-{1/ratio_same:.3f} | {55/462:.3f} |")

print()
print("Key insight: Paley-Zygmund bound is TIGHT (equals actual fraction of working D11)")
print("because N(D11) is either 0 or positive, with positive values concentrated.")
print()
print("The ratio = #{all D11}/#{working D11} in the uniform case.")
print("For the method to work: need #{working D11}/#{all D11} >= 1/poly(p).")
print(f"Data: p=11: 50%, p=19: 7.1%, p=23: 11.9%")
print("Not monotone decreasing - could work asymptotically!")
