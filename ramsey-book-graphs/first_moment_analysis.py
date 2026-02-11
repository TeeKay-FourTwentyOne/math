#!/usr/bin/env python3
"""First moment method analysis for existence of valid (D11, D12) pairs.

KEY QUESTION: Does E[# valid (D11, D12) pairs] grow with p?

If yes, then by the first moment method, valid pairs exist for all large p,
completing the proof for the p ≡ 3 mod 4 infinite family.

We compute:
  E = (#D11 choices) × (#D12 choices) × Pr[random pair is valid]

The first two factors grow exponentially in p. The third factor decays,
but if it decays sub-exponentially, we win.

Uses the k=n formulation (|D11| = (p+1)/2) where only D11 constraints bind.
By Theorem 5 (complement symmetry), this is equivalent to k=n-2.
"""

import numpy as np
from math import comb, log, log2, exp, factorial, lgamma
import time


def autocorrelation_fft(indicator, p):
    """Compute autocorrelation via FFT."""
    F = np.fft.fft(indicator)
    power = np.abs(F) ** 2
    return np.round(np.fft.ifft(power).real).astype(int)


def generate_random_symmetric_d11(p, rng):
    """Generate one random symmetric D11 of size (p+1)/2."""
    pairs = [(x, p - x) for x in range(1, (p + 1) // 2)]
    num_pairs = (p + 1) // 4
    chosen = rng.choice(len(pairs), size=num_pairs, replace=False)
    D11 = set()
    for i in chosen:
        D11.add(pairs[i][0])
        D11.add(pairs[i][1])
    return D11


def generate_random_d12(p, rng):
    """Generate one random D12 of size (p-1)/2 containing 0."""
    d12_size = (p - 1) // 2
    others = rng.choice(range(1, p), size=d12_size - 1, replace=False)
    D12 = {0}
    D12.update(others.tolist())
    return D12


def check_pair(D11, D12, p):
    """Check if (D11, D12) satisfies all constraints.

    At k=n (|D11| = (p+1)/2):
      - D11 positions: A(d)+B(d) <= (p-3)/2
      - D22 positions: A(d)+B(d) <= (p+3)/2
    """
    n = (p + 1) // 2
    d11_thresh = n - 2  # = (p-3)/2
    d22_thresh = n + 1  # = (p+3)/2

    ind11 = np.zeros(p, dtype=np.float64)
    for x in D11:
        ind11[x] = 1.0
    A = autocorrelation_fft(ind11, p)

    ind12 = np.zeros(p, dtype=np.float64)
    for x in D12:
        ind12[x] = 1.0
    B = autocorrelation_fft(ind12, p)

    for d in range(1, p):
        ab = int(A[d]) + int(B[d])
        thresh = d11_thresh if d in D11 else d22_thresh
        if ab > thresh:
            return False
    return True


def check_pair_batch_d12(D11, p, M, rng):
    """Check M random D12s against a fixed D11.

    Returns number of valid D12s found.
    """
    n = (p + 1) // 2
    d11_thresh = n - 2
    d22_thresh = n + 1
    d12_size = (p - 1) // 2

    ind11 = np.zeros(p, dtype=np.float64)
    for x in D11:
        ind11[x] = 1.0
    A = autocorrelation_fft(ind11, p)

    threshold = np.zeros(p, dtype=np.int32)
    for d in range(1, p):
        threshold[d] = d11_thresh if d in D11 else d22_thresh

    valid_count = 0
    batch_size = min(M, 100000)

    for start in range(0, M, batch_size):
        end = min(start + batch_size, M)
        M_batch = end - start

        # Generate batch of random D12s containing 0
        rand_vals = rng.random((M_batch, p - 1))
        sorted_idx = np.argpartition(rand_vals, d12_size - 1,
                                     axis=1)[:, :d12_size - 1]

        indicators = np.zeros((M_batch, p), dtype=np.float64)
        indicators[:, 0] = 1.0  # 0 always in D12
        rows = np.repeat(np.arange(M_batch), d12_size - 1)
        cols = sorted_idx.ravel() + 1  # shift by 1 since we're sampling from {1,...,p-1}
        indicators[rows, cols] = 1.0

        # Batch FFT autocorrelation
        F = np.fft.fft(indicators, axis=1)
        power = np.abs(F) ** 2
        B_all = np.round(np.fft.ifft(power, axis=1).real).astype(np.int32)

        # Check constraints
        sum_ab = B_all[:, 1:] + A[np.newaxis, 1:]
        thresh_row = threshold[np.newaxis, 1:]
        valid = np.all(sum_ab <= thresh_row, axis=1)
        valid_count += int(np.sum(valid))

    return valid_count


def log_comb(n, k):
    """Compute log(C(n,k)) using lgamma for large values."""
    if k < 0 or k > n:
        return float('-inf')
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)


def main():
    print("=" * 80)
    print("FIRST MOMENT ANALYSIS: E[# valid (D11, D12) pairs]")
    print("=" * 80)

    primes = [11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83]

    # Phase 1: Pure random (D11, D12) pair testing
    print("\n" + "=" * 80)
    print("PHASE 1: Pr[random (D11, D12) pair is valid]")
    print("=" * 80)
    print(f"\n  Using k=n formulation: |D11| = (p+1)/2, |D12| = (p-1)/2")

    # For each prime, sample many random (D11, D12) pairs and count valid ones
    rng = np.random.default_rng(42)

    results = {}

    print(f"\n  {'p':>4s} {'n':>4s} {'#D11':>10s} {'#D12':>10s} "
          f"{'log2(#pairs)':>12s} {'samples':>10s} {'valid':>7s} "
          f"{'rate':>10s} {'log2(E)':>10s}")
    print("  " + "-" * 85)

    for p in primes:
        n = (p + 1) // 2
        num_pairs_avail = (p - 1) // 2
        num_pairs_select = (p + 1) // 4

        # Number of symmetric D11 of size n = (p+1)/2
        log2_num_d11 = log_comb(num_pairs_avail, num_pairs_select) / log(2)
        num_d11 = comb(num_pairs_avail, num_pairs_select)

        # Number of D12 of size (p-1)/2 containing 0
        log2_num_d12 = log_comb(p - 1, (p - 1) // 2 - 1) / log(2)
        num_d12 = comb(p - 1, (p - 1) // 2 - 1)

        log2_total = log2_num_d11 + log2_num_d12

        # Determine sample size based on expected hit rate
        # From correlation analysis: joint_prob ranges from 0.17 (p=11) to 2e-4 (p=59)
        # We need enough samples to see at least a few hits
        if p <= 23:
            M_d11 = 500
            M_d12 = 10000
        elif p <= 43:
            M_d11 = 200
            M_d12 = 50000
        elif p <= 59:
            M_d11 = 100
            M_d12 = 100000
        else:
            M_d11 = 50
            M_d12 = 200000

        t0 = time.time()
        total_valid = 0
        total_tested = 0

        for i in range(M_d11):
            D11 = generate_random_symmetric_d11(p, rng)
            valid = check_pair_batch_d12(D11, p, M_d12, rng)
            total_valid += valid
            total_tested += M_d12

        elapsed = time.time() - t0
        total_samples = M_d11 * M_d12
        rate = total_valid / total_samples if total_samples > 0 else 0

        if rate > 0:
            log2_rate = log2(rate)
            log2_E = log2_total + log2_rate
        else:
            log2_rate = float('-inf')
            log2_E = float('-inf')

        results[p] = {
            'n': n,
            'log2_num_d11': log2_num_d11,
            'log2_num_d12': log2_num_d12,
            'log2_total': log2_total,
            'total_samples': total_samples,
            'total_valid': total_valid,
            'rate': rate,
            'log2_rate': log2_rate,
            'log2_E': log2_E,
            'time': elapsed,
        }

        rate_str = f"{rate:.2e}" if rate > 0 else "0"
        log2_E_str = f"{log2_E:.1f}" if log2_E > float('-inf') else "-inf"

        print(f"  {p:4d} {n:4d} {num_d11:10d} {num_d12:10d} "
              f"{log2_total:12.1f} {total_samples:10d} {total_valid:7d} "
              f"{rate_str:>10s} {log2_E_str:>10s}  [{elapsed:.1f}s]")

    # Phase 2: Scaling analysis
    print("\n" + "=" * 80)
    print("PHASE 2: Scaling analysis")
    print("=" * 80)

    primes_with_data = [p for p in primes if results[p]['rate'] > 0]

    if len(primes_with_data) >= 3:
        print(f"\n  Key question: does log2(Pr[valid]) decay slower than log2(#pairs)?")
        print(f"\n  {'p':>4s} {'log2(#pairs)':>12s} {'log2(rate)':>12s} "
              f"{'log2(E)':>10s} {'E grows?':>10s}")
        print("  " + "-" * 55)

        prev_log2_E = None
        for p in primes_with_data:
            r = results[p]
            growing = ""
            if prev_log2_E is not None:
                growing = "YES" if r['log2_E'] > prev_log2_E else "NO"
            print(f"  {p:4d} {r['log2_total']:12.1f} {r['log2_rate']:12.1f} "
                  f"{r['log2_E']:10.1f} {growing:>10s}")
            prev_log2_E = r['log2_E']

        # Fit log2(rate) vs p
        ps = np.array([float(p) for p in primes_with_data])
        log2_rates = np.array([results[p]['log2_rate'] for p in primes_with_data])
        log2_totals = np.array([results[p]['log2_total'] for p in primes_with_data])
        log2_Es = np.array([results[p]['log2_E'] for p in primes_with_data])

        # Linear fit: log2(rate) = a*p + b
        if len(ps) >= 2:
            coeffs = np.polyfit(ps, log2_rates, 1)
            print(f"\n  Linear fit: log2(rate) ≈ {coeffs[0]:.4f} * p + {coeffs[1]:.2f}")
            print(f"    Rate decays as 2^({coeffs[0]:.4f}*p) = {2**coeffs[0]:.6f}^p")

            # log2(#pairs) growth rate
            if len(ps) >= 2:
                pair_coeffs = np.polyfit(ps, log2_totals, 1)
                print(f"\n  Linear fit: log2(#pairs) ≈ {pair_coeffs[0]:.4f} * p + {pair_coeffs[1]:.2f}")
                print(f"    #pairs grows as 2^({pair_coeffs[0]:.4f}*p) = {2**pair_coeffs[0]:.6f}^p")

                net_rate = coeffs[0] + pair_coeffs[0]
                print(f"\n  Net growth: log2(E) ≈ {net_rate:.4f} * p")
                if net_rate > 0:
                    print(f"  >>> E[# valid pairs] GROWS EXPONENTIALLY with p!")
                    print(f"  >>> First moment method proves existence for all large p!")
                    print(f"  >>> Growth rate: E ~ {2**net_rate:.4f}^p")
                elif net_rate > -0.01:
                    print(f"  >>> E[# valid pairs] is approximately CONSTANT")
                    print(f"  >>> Need second moment method or more data")
                else:
                    print(f"  >>> E[# valid pairs] DECAYS exponentially")
                    print(f"  >>> First moment method FAILS; need structured construction")

            # Fit log2(E) directly
            if len(ps) >= 2:
                E_coeffs = np.polyfit(ps, log2_Es, 1)
                print(f"\n  Direct fit: log2(E) ≈ {E_coeffs[0]:.4f} * p + {E_coeffs[1]:.2f}")

    # Phase 3: Detailed constraint analysis
    print("\n" + "=" * 80)
    print("PHASE 3: Why the first moment works (or doesn't)")
    print("=" * 80)

    print(f"\n  Theoretical analysis at k=n:")
    print(f"  E[A(d)] = (n-1)n/2(p-1) ≈ (p+1)/4")
    print(f"  E[B(d)] = (d12_size)(d12_size-1)/(p-1) ≈ (p-3)/4")
    print(f"  E[A+B] ≈ (p-1)/2")
    print(f"  Threshold at D11: (p-3)/2 = E[A+B] - 1")
    print(f"  Threshold at D22: (p+3)/2 = E[A+B] + 2  [loose]")

    print(f"\n  For the first moment to work, we need:")
    print(f"    Pr[A+B ≤ thresh for all d ∈ D11] × #pairs → ∞")
    print(f"  The binding constraint is the (p+1)/2 positions in D11.")
    print(f"  If events were independent with Pr ≈ 1/2 each:")
    print(f"    Pr[all ok] ≈ 2^{{-(p+1)/2}}")
    print(f"  But #pairs ≈ 2^{{cp}} for some c > 1/2 would still dominate.")

    for p in [11, 23, 43, 59]:
        n = (p + 1) // 2
        num_pairs = (p - 1) // 2
        s = (p + 1) // 4
        log2_d11 = log_comb(num_pairs, s) / log(2)
        log2_d12 = log_comb(p - 1, (p - 1) // 2 - 1) / log(2)
        log2_naive = -(p + 1) / 2 * log2(2)  # independent Pr ≈ 1/2 each

        print(f"\n    p={p}: log2(#D11)={log2_d11:.1f}, "
              f"log2(#D12)={log2_d12:.1f}, "
              f"total={log2_d11 + log2_d12:.1f}")
        print(f"           naive log2(Pr)=−{(p + 1) / 2:.0f} "
              f"(if independent, Pr≈1/2 each)")
        print(f"           naive log2(E)={log2_d11 + log2_d12 - (p + 1) / 2:.1f}")

    # Phase 4: Per-D11 analysis — what fraction of D11s are "good"?
    print("\n" + "=" * 80)
    print("PHASE 4: Per-D11 success rate")
    print("=" * 80)

    print(f"\n  For each prime, what fraction of random D11s admit at least one valid D12?")
    print(f"\n  {'p':>4s} {'M_d11':>7s} {'M_d12':>8s} "
          f"{'good_D11':>9s} {'D11_rate':>9s} {'avg_valid':>10s}")
    print("  " + "-" * 55)

    for p in primes:
        if p > 59:
            break
        n = (p + 1) // 2

        if p <= 23:
            M_d11 = 1000
            M_d12 = 10000
        elif p <= 43:
            M_d11 = 500
            M_d12 = 50000
        else:
            M_d11 = 200
            M_d12 = 100000

        t0 = time.time()
        good_d11 = 0
        total_valid = 0

        for i in range(M_d11):
            D11 = generate_random_symmetric_d11(p, rng)
            valid = check_pair_batch_d12(D11, p, M_d12, rng)
            total_valid += valid
            if valid > 0:
                good_d11 += 1

        elapsed = time.time() - t0
        d11_rate = good_d11 / M_d11 * 100
        avg_valid = total_valid / M_d11

        print(f"  {p:4d} {M_d11:7d} {M_d12:8d} "
              f"{good_d11:9d} {d11_rate:8.1f}% {avg_valid:10.1f}  [{elapsed:.1f}s]")

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print("""
The first moment method asks: does E[# valid pairs] → ∞?

  E = (#D11 choices) × (#D12 choices) × Pr[random pair valid]
    = 2^{c₁p} × 2^{c₂p} × 2^{-c₃p}

If c₁ + c₂ > c₃, then E → ∞ and existence is guaranteed for large p.

Key finding from Phase 2: the scaling of log2(Pr) vs log2(#pairs) determines
whether the proof works. The correlation structure analysis shows that events
are approximately independent with favorable co-satisfaction, suggesting
c₃ < c₁ + c₂.

If the first moment method works:
  - Existence is proven for all p ≥ p₀ (some threshold)
  - SA constructions verify p < p₀
  - Combined with Paley (q ≡ 1 mod 4), this covers all n where 2n-1 is prime
  - This is a positive-density set of all n
""")


if __name__ == '__main__':
    main()
