#!/usr/bin/env python3
"""
Compare Gaussian hyperplane probability to exact discrete probability.
Shows that the Gaussian approximation fails badly for this problem.
"""
from math import comb, log2, sqrt
from scipy.stats import norm

def compare(p):
    n = (p + 1) // 2
    k = (p - 3) // 2  # |S| = |D12| - 1
    r = (p + 1) // 4
    r_prime = (p - 3) // 4
    R = r + r_prime

    E_B = (p - 3) / 4
    E_A = (p + 1) / 4
    sigma_cond = sqrt((p+1)*(p-1)/(16*(p-2)) * (1 - 2/(p-1)))

    T_red = (p - 3) / 2 - E_A  # binding threshold minus E[A]
    T_blue = (p + 3) / 2 - E_A

    z_red = (T_red - E_B) / sigma_cond
    z_blue = (T_blue - E_B) / sigma_cond

    log2_prod_phi = r * log2(norm.cdf(z_red)) + r_prime * log2(norm.cdf(z_blue))
    log2_C = log2(comb(p - 1, k))

    # Known exact data from enumeration
    exact_data = {
        11: (20, 210),     # best D11: 20 valid out of C(10,4)=210
        19: (18, 43758),   # best D11: 18 valid out of C(18,8)=43758
        23: (198, 646646), # best D11: 198 valid out of C(22,10)
    }

    print(f"p={p}: R={R}, r={r}, r'={r_prime}")
    print(f"  Product of Gaussian marginals: -log2 = {-log2_prod_phi:.2f}")
    print(f"  Budget log2(C(p-1,k)): {log2_C:.2f}")
    print(f"  Gaussian E[valid] (if c0=1): 2^{{{log2_C + log2_prod_phi:.2f}}}")

    if p in exact_data:
        n_valid, n_total = exact_data[p]
        pr_exact = n_valid / n_total
        log2_pr = log2(pr_exact)
        log2_E_exact = log2(n_valid)  # = log2(C * Pr) = log2(n_valid)

        # Effective c0 (exact vs product of Gaussian marginals)
        c0_effective = pr_exact / (2 ** log2_prod_phi)

        print(f"  EXACT Pr[valid] = {n_valid}/{n_total} = {pr_exact:.6e}")
        print(f"  EXACT -log2(Pr) = {-log2_pr:.2f}")
        print(f"  EXACT E[valid] = {n_valid}")
        print(f"  Effective c0 (exact/prod_Phi) = {c0_effective:.2f}")
        print(f"  log2(effective c0) = {log2(c0_effective):.2f}")
        print(f"  EXACT is {c0_effective:.1f}x BETTER than Gaussian product")
    print()

print("=" * 60)
print("COMPARISON: GAUSSIAN vs EXACT DISCRETE")
print("=" * 60)
print()
print("Key finding: Gaussian hyperplane MC gives c0 << 1")
print("(~2^{-8} at p=11, ~2^{-11} at p=19)")
print("But EXACT discrete probability is MUCH larger.")
print()

for p in [11, 19, 23]:
    compare(p)

print("CONCLUSION:")
print("The Gaussian approximation FAILS for this problem.")
print("The product of Gaussian marginals is a good lower bound")
print("on the EXACT Pr[valid] (c0_effective > 1), but this is")
print("NOT because Slepian helps -- it's because the discrete")
print("distribution has much higher joint probability than the")
print("Gaussian, overwhelming the negative correlation loss.")
