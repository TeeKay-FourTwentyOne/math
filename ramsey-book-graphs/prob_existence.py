"""
Probabilistic existence proof for R(B_{n-1}, B_n) = 4n-1 when m = 2n-1 ≡ 3 mod 4.

THEOREM (to prove): For every prime p ≡ 3 mod 4 (sufficiently large), there exist
D11, D12 satisfying all constraints.

APPROACH: Lovász Local Lemma (symmetric form).

For fixed D11, choose D12 uniformly at random among (p-1)/2-element subsets of Z_p
containing 0. Define bad events:
  E_d : A(d) + B(d) > n-2   for d ∈ D11    (binding constraint)
  F_d : A(d) + B(d) > n+1   for d ∈ D22    (loose constraint)

If we can show these bad events satisfy LLL conditions, then Pr[no bad event] > 0,
proving existence.

KEY QUANTITIES TO COMPUTE:
  1. Pr[E_d] for each d ∈ D11
  2. Pr[F_d] for each d ∈ D22
  3. Dependency structure between events
  4. LLL condition: e·d_max·p_max < 1

This script computes these exactly for small primes and estimates them for larger ones.
"""

import sys
import os
import json
import time
import math
import numpy as np
from itertools import combinations
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import BlockCirculantGraph, verify_construction, Delta


def autocorrelation_fft(indicator, m):
    fft_val = np.fft.fft(indicator)
    autocorr = np.fft.ifft(np.abs(fft_val) ** 2).real
    return np.round(autocorr).astype(int)


# Known solutions
KNOWN = {
    11: {"D11": {1, 2, 4, 7, 9, 10}, "D12": {0, 4, 5, 7, 10}},
    19: {"D11": {1, 2, 3, 6, 8, 11, 13, 16, 17, 18},
         "D12": {0, 2, 6, 8, 9, 12, 13, 17, 18}},
    23: {"D11": {5, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17, 18},
         "D12": {0, 1, 2, 6, 10, 13, 14, 16, 18, 20, 21}},
    31: {"D11": {6, 7, 8, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 23, 24, 25},
         "D12": {0, 1, 2, 3, 8, 11, 12, 13, 15, 18, 20, 21, 24, 27, 29}},
    43: {"D11": {1, 2, 5, 10, 11, 13, 16, 17, 18, 19, 20, 23, 24, 25, 26, 27,
                 30, 32, 33, 38, 41, 42},
         "D12": {0, 2, 5, 6, 8, 11, 15, 16, 20, 24, 25, 27, 28, 31, 32, 34,
                 35, 36, 37, 39, 41}},
}


def exact_b_distribution(p, D11, d_query, num_samples=500000):
    """
    Sample B(d_query) = Delta(D12, D12, d_query) for random D12 containing 0
    with |D12| = (p-1)/2.

    Returns histogram of B values.
    """
    d12_size = (p - 1) // 2
    rng = np.random.default_rng(42)
    histogram = defaultdict(int)

    for _ in range(num_samples):
        others = rng.choice(range(1, p), size=d12_size - 1, replace=False)
        D12 = set([0]) | set(others)
        b = sum(1 for x in D12 if (x - d_query) % p in D12)
        histogram[b] += 1

    return dict(histogram)


def exact_ab_distribution_mc(p, D11, d_query, A_d, num_samples=500000):
    """
    Sample A(d)+B(d) for random D12 at a given position d.
    A(d) is fixed, B(d) is random.
    Returns: Pr[A+B > threshold] for each threshold.
    """
    d12_size = (p - 1) // 2
    n = (p + 1) // 2
    tight_thresh = n - 2
    loose_thresh = n + 1

    rng = np.random.default_rng(42)

    exceed_tight = 0
    exceed_loose = 0
    b_values = []

    for _ in range(num_samples):
        others = rng.choice(range(1, p), size=d12_size - 1, replace=False)
        D12_set = set([0]) | set(int(x) for x in others)
        b = sum(1 for x in D12_set if (x - d_query) % p in D12_set)
        ab = A_d + b
        b_values.append(b)

        if ab > tight_thresh:
            exceed_tight += 1
        if ab > loose_thresh:
            exceed_loose += 1

    return {
        "A_d": A_d,
        "mean_B": np.mean(b_values),
        "std_B": np.std(b_values),
        "mean_AB": A_d + np.mean(b_values),
        "pr_exceed_tight": exceed_tight / num_samples,
        "pr_exceed_loose": exceed_loose / num_samples,
        "tight_thresh": tight_thresh,
        "loose_thresh": loose_thresh,
    }


def theoretical_b_moments(p, d):
    """
    Compute E[B(d)] and Var[B(d)] for random D12 of size k=(p-1)/2 containing 0.

    B(d) = #{x ∈ D12 : (x-d) mod p ∈ D12}

    Since 0 ∈ D12 always, decompose:
    B(d) = [(-d) ∈ D12] + [d ∈ D12] + #{x ∈ D12\{0} : x≠d and (x-d)∈D12\{0}}

    Wait, more carefully:
    B(d) = #{x ∈ D12 : (x-d) ∈ D12}

    x=0: contributes 1 iff (-d)≡p-d ∈ D12
    x=d: contributes 1 iff 0 ∈ D12 = 1 (always, since 0∈D12) — but only if d∈D12

    Actually, if d≠0: x can be any element of D12.
    If 0∈D12 and d∈D12: x=d gives (d-d)=0∈D12, contributing 1.
    If 0∈D12 and d∉D12: x=d not in D12, doesn't contribute.
    If d∈D12 and x=0: gives (0-d)=-d, contributes 1 iff -d∈D12.

    Let's use the exact combinatorial formula.
    D12 has k = (p-1)/2 elements, with 0 always included. The other k-1 elements
    are chosen from {1,...,p-1}.

    We want E[B(d)] and Var[B(d)] where B(d) = Σ_{x∈Z_p} 1_{D12}(x) · 1_{D12}(x-d).
    """
    k = (p - 1) // 2  # |D12|
    pool = p - 1  # elements to choose k-1 from (i.e., {1,...,p-1})

    # For d ≠ 0:
    # B(d) = Σ_{x∈Z_p} f(x)f(x-d) where f = indicator of D12, f(0)=1 always
    #
    # = f(0)f(-d) + f(d)f(0) + Σ_{x≠0,x≠d} f(x)f(x-d)
    # = f(-d) + f(d) + Σ_{x≠0,x≠d} f(x)f(x-d)
    #
    # But we need to handle x-d mod p. If x ≠ 0 and x ≠ d, then x-d ≠ -d and x-d ≠ 0.
    # Also need x-d ≠ x, which gives d ≠ 0 (true).
    # And x-d can equal d only if x = 2d.

    # For the random model: each element of {1,...,p-1} is independently in D12 with
    # probability (k-1)/(p-1) marginally, but they're drawn without replacement.
    # Actually, D12\{0} is a uniform random (k-1)-subset of {1,...,p-1}.

    # Probability that a specific element a ∈ {1,...,p-1} is in D12:
    q = (k - 1) / pool  # = ((p-1)/2 - 1)/(p-1) = (p-3)/(2(p-1))

    # E[f(a)] = q for a ≠ 0
    # E[f(a)f(b)] for a,b ≠ 0, a ≠ b:
    q2 = (k - 1) * (k - 2) / (pool * (pool - 1))

    # E[B(d)] = E[f(-d)] + E[f(d)] + Σ_{x≠0,d} E[f(x)f(x-d)]
    # If -d ≠ d (i.e., d ≠ p/2, always true for odd p):
    # and -d ≠ 0, d ≠ 0 (true since d>0):

    # Number of terms in the sum Σ_{x≠0,d}: p-2 terms
    # For each x ≠ 0, d: x ∈ {1,...,p-1}\{d}, and x-d ∈ {1,...,p-1}\{-d mod p}
    # Actually x-d mod p. If x ≠ 0 and x ≠ d, then x-d ≠ -d and x-d ≠ 0.
    # So both x and x-d are in {1,...,p-1}, and they're distinct (since d≠0).
    # BUT x-d might equal d (when x=2d) or might equal -d (when x=0, excluded).
    # And x might equal -d (when x = p-d).
    # These are all possible but don't change the expectation (just different elements).

    # E[f(x)f(x-d)] depends on whether x and x-d are the same element:
    # They're always different (d ≠ 0), so E[f(x)f(x-d)] = q2 for each.

    # Wait, but there might be cases where {x, x-d} = {a, -d} or similar overlaps
    # with the f(-d) and f(d) terms. Let me be more careful.

    # Actually, the cleanest way:
    # B(d) = Σ_{x=0}^{p-1} f(x)f(x-d)
    # Group into: x=0, x=d, and x ∈ {1,...,p-1}\{d}

    # x=0: f(0)f(-d) = 1 · f(-d) = f(p-d)
    #   E = q (since p-d ∈ {1,...,p-1})

    # x=d: f(d)f(0) = f(d) · 1 = f(d)
    #   E = q

    # x ∈ {1,...,p-1}\{d}: f(x)f(x-d)
    #   There are p-2 such terms.
    #   For each, x ∈ {1,...,p-1} and (x-d) mod p ∈ {1,...,p-1}\{0}
    #   (since x ≠ d means x-d ≠ 0)
    #   And x ≠ x-d (since d ≠ 0).
    #   Special case: is (x-d) = d when x=2d? Yes but that's fine.
    #   Is x = p-d? Then x-d = p-2d, which is some other element.
    #
    #   For each such pair (x, x-d), both in {1,...,p-1}, distinct:
    #   E[f(x)f(x-d)] = q2
    #
    #   BUT: some of these pairs involve p-d or d.
    #   For x = p-d: pair is (p-d, p-2d). These are both random, correlated with
    #   the f(p-d) from the x=0 term and f(d) from x=d term.
    #   This doesn't change E[f(x)f(x-d)] = q2, but it DOES affect
    #   the covariance structure when computing Var[B(d)].

    EB = 2 * q + (p - 2) * q2

    # Var[B(d)] is more complex. Let me compute E[B(d)^2] - E[B(d)]^2.
    # B(d)^2 = (Σ_x f(x)f(x-d))^2 = Σ_{x,y} f(x)f(x-d)f(y)f(y-d)
    # This involves fourth moments of the hypergeometric distribution.

    return {
        "E_B": EB,
        "q": q,
        "q2": q2,
    }


def monte_carlo_lll(p, D11, num_samples=500000):
    """
    Monte Carlo computation of LLL parameters for the binding constraints.

    For each d ∈ D11, compute Pr[E_d] where E_d = {A(d)+B(d) > n-2}.
    Also estimate the dependency neighborhood size.
    """
    n = (p + 1) // 2
    d12_size = (p - 1) // 2
    D11_set = set(D11)
    D22 = set(range(1, p)) - D11_set

    # Compute A(d)
    ind11 = np.zeros(p, dtype=np.float64)
    for j in D11:
        ind11[j] = 1.0
    A = autocorrelation_fft(ind11, p)

    tight_thresh = n - 2
    loose_thresh = n + 1

    rng = np.random.default_rng(42)

    # For each constraint position, track how often it's violated
    d11_violations = defaultdict(int)
    d22_violations = defaultdict(int)

    # Track joint violations (for dependency estimation)
    joint_violations = defaultdict(int)  # (d1, d2) -> count

    # Also track B(d) statistics
    B_sums = defaultdict(float)
    B_sq_sums = defaultdict(float)

    all_valid = 0

    for trial in range(num_samples):
        others = rng.choice(range(1, p), size=d12_size - 1, replace=False)
        D12_set = set([0]) | set(int(x) for x in others)

        # Compute B(d) for all d
        ind12 = np.zeros(p, dtype=np.float64)
        for j in D12_set:
            ind12[j] = 1.0
        B = autocorrelation_fft(ind12, p)

        violated_d11 = set()
        violated_d22 = set()

        for d in D11:
            ab = int(A[d]) + int(B[d])
            B_sums[d] += int(B[d])
            B_sq_sums[d] += int(B[d]) ** 2
            if ab > tight_thresh:
                d11_violations[d] += 1
                violated_d11.add(d)

        for d in D22:
            ab = int(A[d]) + int(B[d])
            if ab > loose_thresh:
                d22_violations[d] += 1
                violated_d22.add(d)

        if not violated_d11 and not violated_d22:
            all_valid += 1

        # Joint violations (subsample for speed)
        if trial < min(num_samples, 100000):
            for d1 in violated_d11:
                for d2 in violated_d11:
                    if d1 < d2:
                        joint_violations[(d1, d2)] += 1

    # Compute results
    d11_sorted = sorted(D11)
    d22_sorted = sorted(D22)

    # Individual violation probabilities
    pr_d11 = {d: d11_violations[d] / num_samples for d in d11_sorted}
    pr_d22 = {d: d22_violations[d] / num_samples for d in d22_sorted}

    # B(d) statistics
    B_means = {d: B_sums[d] / num_samples for d in d11_sorted}
    B_stds = {d: math.sqrt(B_sq_sums[d] / num_samples - (B_sums[d] / num_samples) ** 2)
              for d in d11_sorted}

    # Max individual probability (for symmetric LLL)
    pr_d11_vals = [pr_d11[d] for d in d11_sorted]
    pr_d22_vals = [pr_d22[d] for d in d22_sorted]

    max_pr_d11 = max(pr_d11_vals) if pr_d11_vals else 0
    max_pr_d22 = max(pr_d22_vals) if pr_d22_vals else 0
    max_pr = max(max_pr_d11, max_pr_d22)

    # Dependency: for LLL, events E_d and E_{d'} are dependent if they share
    # random variables. B(d) depends on ALL elements of D12, so technically
    # all events are dependent. But the EFFECTIVE dependency can be measured
    # by excess joint probability.

    # For LLL to work, we need: e * p_max * (total_events) < 1
    total_events = len(D11) + len(D22)
    lll_check = math.e * max_pr * total_events

    # More refined: symmetric LLL with neighborhood size
    # In the symmetric case, need e * p * d < 1 where d = dependency neighborhood size.
    # Since all events depend on D12, d = total_events - 1 (worst case).
    # This only works if p < 1/(e * total_events).

    # Compute conditional probabilities for dependency analysis
    cond_probs = {}
    joint_sample_size = min(num_samples, 100000)
    for (d1, d2), count in sorted(joint_violations.items())[:20]:
        if d11_violations[d1] > 0:
            pr_d2_given_d1 = count / min(d11_violations[d1], joint_sample_size * d11_violations[d1] / num_samples)
            cond_probs[(d1, d2)] = {
                "pr_d1": d11_violations[d1] / num_samples,
                "pr_d2": d11_violations[d2] / num_samples,
                "pr_joint": count / joint_sample_size,
                "pr_d2_given_d1": min(pr_d2_given_d1, 1.0),
            }

    return {
        "p": p,
        "n": n,
        "num_samples": num_samples,
        "A_values": {d: int(A[d]) for d in d11_sorted},
        "B_means": {d: B_means[d] for d in d11_sorted},
        "B_stds": {d: B_stds[d] for d in d11_sorted},
        "pr_d11": pr_d11,
        "pr_d22": pr_d22,
        "max_pr_d11": max_pr_d11,
        "max_pr_d22": max_pr_d22,
        "all_valid_count": all_valid,
        "all_valid_rate": all_valid / num_samples,
        "total_constraints": total_events,
        "lll_quantity": lll_check,
        "lll_satisfied": lll_check < 1,
        "cond_probs_sample": cond_probs,
    }


def analyze_random_d11_lll(p, num_d11=100, num_mc=50000):
    """
    For random D11, compute LLL quantities.
    This tests whether the proof works for GENERIC D11 (not just known solutions).
    """
    n = (p + 1) // 2
    d12_size = (p - 1) // 2
    tight_thresh = n - 2
    loose_thresh = n + 1

    pairs = []
    seen = set()
    for d in range(1, p):
        if d not in seen:
            pairs.append((d, (-d) % p))
            seen.add(d)
            seen.add((-d) % p)

    num_pairs_needed = (p + 1) // 4

    rng = np.random.default_rng(99)

    results = []

    for trial in range(num_d11):
        # Random symmetric D11
        chosen = rng.choice(len(pairs), size=num_pairs_needed, replace=False)
        D11 = set()
        for i in chosen:
            D11.add(pairs[i][0])
            D11.add(pairs[i][1])

        D22 = set(range(1, p)) - D11

        # Compute A(d)
        ind11 = np.zeros(p, dtype=np.float64)
        for j in D11:
            ind11[j] = 1.0
        A = autocorrelation_fft(ind11, p)

        # MC for violation probabilities
        d11_viol = 0
        d22_viol = 0
        all_valid = 0
        max_d11_pr = 0
        max_d22_pr = 0

        per_d_viol = defaultdict(int)

        for _ in range(num_mc):
            others = rng.choice(range(1, p), size=d12_size - 1, replace=False)
            D12_set = set([0]) | set(int(x) for x in others)

            ind12 = np.zeros(p, dtype=np.float64)
            for j in D12_set:
                ind12[j] = 1.0
            B = autocorrelation_fft(ind12, p)

            any_viol = False
            for d in D11:
                if int(A[d]) + int(B[d]) > tight_thresh:
                    per_d_viol[d] += 1
                    any_viol = True
            for d in D22:
                if int(A[d]) + int(B[d]) > loose_thresh:
                    any_viol = True

            if not any_viol:
                all_valid += 1

        max_pr = max((per_d_viol[d] / num_mc for d in D11), default=0)
        total_events = len(D11) + len(D22)
        lll_q = math.e * max_pr * total_events

        results.append({
            "max_pr": max_pr,
            "all_valid_rate": all_valid / num_mc,
            "lll_quantity": lll_q,
        })

    return results


def main():
    primes = [11, 19, 23, 31, 43]

    print("=" * 90)
    print("PROBABILISTIC EXISTENCE ANALYSIS")
    print("=" * 90)
    print("""
Goal: Apply Lovász Local Lemma to prove existence of valid (D11, D12) pairs.

For random D12 (uniform among (p-1)/2-subsets of Z_p containing 0):
  - Bad event E_d (d ∈ D11): A(d) + B(d) > n-2
  - Bad event F_d (d ∈ D22): A(d) + B(d) > n+1

LLL (symmetric form): If e·p_max·D < 1, where
  p_max = max violation probability
  D = dependency neighborhood size
then Pr[no bad event] > 0.
""")

    all_results = {}

    for p in primes:
        n = (p + 1) // 2
        D11 = KNOWN[p]["D11"]

        print(f"\n{'='*80}")
        print(f"p = {p}, n = {n}")
        print(f"{'='*80}")

        # Theoretical B moments
        theory = theoretical_b_moments(p, 1)
        print(f"\n  Theoretical B(d) moments (for generic d):")
        print(f"    E[B(d)] = {theory['E_B']:.4f}")
        print(f"    q = Pr[a ∈ D12] = {theory['q']:.4f}")
        print(f"    q2 = Pr[a,b ∈ D12] = {theory['q2']:.4f}")
        print(f"    E[A+B] ≈ (p-1)/2 = {(p-1)/2:.1f}")
        print(f"    Tight threshold: {n-2}, gap from mean: {n-2 - (p-1)/2:.1f}")

        # Monte Carlo LLL analysis with known D11
        num_mc = min(500000, max(100000, 20000000 // p))
        print(f"\n  Monte Carlo LLL analysis ({num_mc:,} samples):")
        t0 = time.time()
        mc = monte_carlo_lll(p, D11, num_mc)
        elapsed = time.time() - t0

        print(f"    Time: {elapsed:.1f}s")
        b_means_str = [f"{mc['B_means'][d]:.2f}" for d in sorted(D11)[:6]]
        b_stds_str = [f"{mc['B_stds'][d]:.2f}" for d in sorted(D11)[:6]]
        print(f"    B(d) means: {b_means_str}...")
        print(f"    B(d) stds:  {b_stds_str}...")

        print(f"\n    Pr[E_d] for d ∈ D11 (binding constraint violations):")
        d11_sorted = sorted(D11)
        for d in d11_sorted:
            A_d = mc["A_values"][d]
            pr = mc["pr_d11"][d]
            slack = n - 2 - A_d - mc["B_means"][d]
            print(f"      d={d:3d}: A(d)={A_d}, E[B]={mc['B_means'][d]:.2f}, "
                  f"slack={slack:.2f}, Pr[exceed]={pr:.6f}")

        print(f"\n    Pr[F_d] for d ∈ D22 (loose constraint violations):")
        d22_sorted = sorted(set(range(1, p)) - D11)
        for d in d22_sorted[:8]:
            pr = mc["pr_d22"][d]
            print(f"      d={d:3d}: Pr[exceed]={pr:.6f}")
        if len(d22_sorted) > 8:
            print(f"      ... ({len(d22_sorted)} total D22 positions)")

        print(f"\n    LLL ANALYSIS:")
        print(f"      max Pr[E_d] = {mc['max_pr_d11']:.6f}")
        print(f"      max Pr[F_d] = {mc['max_pr_d22']:.6f}")
        print(f"      Total constraints: {mc['total_constraints']}")
        print(f"      e·p_max·D = {mc['lll_quantity']:.4f}")
        print(f"      LLL condition (e·p·D < 1): {'SATISFIED' if mc['lll_satisfied'] else 'NOT satisfied'}")

        print(f"\n    Direct count: {mc['all_valid_count']}/{num_mc} fully valid "
              f"({mc['all_valid_rate']:.6f})")

        if mc['all_valid_count'] > 0 and mc['total_constraints'] > 0:
            # Union bound comparison
            union_bound = sum(mc["pr_d11"][d] for d in d11_sorted) + sum(mc["pr_d22"][d] for d in d22_sorted)
            print(f"    Union bound on Pr[any violation]: {union_bound:.4f}")
            print(f"    Actual Pr[any violation]: {1 - mc['all_valid_rate']:.6f}")

        # Joint probabilities (dependency structure)
        if mc["cond_probs_sample"]:
            print(f"\n    DEPENDENCY STRUCTURE (sample of joint violations):")
            for (d1, d2), info in list(mc["cond_probs_sample"].items())[:5]:
                pr_indep = info["pr_d1"] * info["pr_d2"]
                ratio = info["pr_joint"] / pr_indep if pr_indep > 0 else float("inf")
                print(f"      ({d1},{d2}): Pr[joint]={info['pr_joint']:.6f}, "
                      f"Pr[indep]={pr_indep:.6f}, ratio={ratio:.2f}")

        all_results[p] = mc

    # Summary and trend analysis
    print(f"\n{'='*90}")
    print("SUMMARY AND SCALING TRENDS")
    print(f"{'='*90}")

    print(f"\n  {'p':>5s} {'n':>4s} {'max_pr':>10s} {'total_C':>8s} {'e·p·D':>10s} "
          f"{'LLL?':>6s} {'valid_rate':>12s}")
    print(f"  {'-'*5} {'-'*4} {'-'*10} {'-'*8} {'-'*10} {'-'*6} {'-'*12}")

    for p in primes:
        mc = all_results[p]
        n = (p + 1) // 2
        print(f"  {p:5d} {n:4d} {mc['max_pr_d11']:10.6f} {mc['total_constraints']:8d} "
              f"{mc['lll_quantity']:10.4f} {'YES' if mc['lll_satisfied'] else 'NO':>6s} "
              f"{mc['all_valid_rate']:12.6f}")

    # Extrapolation: if max_pr scales as C/sqrt(p), LLL needs C·p/sqrt(p) < 1/e
    # i.e., C·sqrt(p) < 1/e
    if len(primes) >= 3:
        print(f"\n  SCALING ANALYSIS:")
        for p in primes:
            mc = all_results[p]
            if mc["max_pr_d11"] > 0:
                scaled = mc["max_pr_d11"] * math.sqrt(p)
                print(f"    p={p}: max_pr * sqrt(p) = {scaled:.4f}")

    # Random D11 analysis for p=11 and p=19
    print(f"\n{'='*90}")
    print("RANDOM D11 LLL ANALYSIS")
    print(f"{'='*90}")

    for p in [11, 19]:
        n = (p + 1) // 2
        print(f"\n  p={p}: Testing 50 random D11 ({10000} MC each):")
        results = analyze_random_d11_lll(p, num_d11=50, num_mc=10000)

        max_prs = [r["max_pr"] for r in results]
        lll_qs = [r["lll_quantity"] for r in results]
        valid_rates = [r["all_valid_rate"] for r in results]

        print(f"    max_pr: min={min(max_prs):.4f}, max={max(max_prs):.4f}, "
              f"mean={sum(max_prs)/len(max_prs):.4f}")
        print(f"    e·p·D: min={min(lll_qs):.4f}, max={max(lll_qs):.4f}, "
              f"mean={sum(lll_qs)/len(lll_qs):.4f}")
        print(f"    LLL satisfied: {sum(1 for q in lll_qs if q < 1)}/{len(lll_qs)}")
        print(f"    Valid rate: min={min(valid_rates):.6f}, max={max(valid_rates):.6f}")

    # Save results
    output = {str(p): {
        "max_pr_d11": all_results[p]["max_pr_d11"],
        "max_pr_d22": all_results[p]["max_pr_d22"],
        "lll_quantity": all_results[p]["lll_quantity"],
        "lll_satisfied": all_results[p]["lll_satisfied"],
        "all_valid_rate": all_results[p]["all_valid_rate"],
    } for p in primes}

    outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "prob_existence_results.json")
    with open(outpath, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to {outpath}")


if __name__ == "__main__":
    main()
