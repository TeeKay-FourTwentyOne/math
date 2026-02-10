"""
Deep pattern analysis of all known Ramsey book graph constructions.

Computes:
1. Spectral/Fourier analysis of difference sets
2. Cyclotomic class analysis
3. D12 structural relationships
4. Categorization of all n <= 50
5. Composite m=45 (n=23) analysis
6. Delta fluctuation bounds

Outputs comprehensive findings for pattern_analysis.md.
"""

import sys, os, math, json
from collections import defaultdict
import cmath

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import BlockCirculantGraph, verify_construction, Delta, Sigma

# ============================================================================
# All known constructions
# ============================================================================

KNOWN = {
    6: {
        "D11": {3, 5, 6, 8},
        "D12": {0, 1, 4, 6, 7},
        "D22": {1, 2, 4, 7, 9, 10},
    },
    8: {
        "D11": {3, 6, 7, 8, 9, 12},
        "D12": {0, 1, 4, 6, 8, 9, 13},
        "D22": {1, 2, 4, 5, 10, 11, 13, 14},
    },
    10: {
        "D11": {4, 5, 7, 9, 10, 12, 14, 15},
        "D12": {0, 1, 2, 6, 7, 10, 11, 13, 17},
        "D22": {1, 2, 3, 6, 8, 11, 13, 16, 17, 18},
    },
    12: {
        "D11": {5, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17, 18},
        "D12": {0, 1, 2, 6, 10, 13, 14, 16, 18, 20, 21},
        "D22": {1, 2, 3, 4, 10, 13, 19, 20, 21, 22},
    },
    14: {
        "D11": {5, 7, 8, 9, 10, 11, 13, 14, 16, 17, 18, 19, 20, 22},
        "D12": {0, 1, 2, 7, 8, 10, 13, 14, 17, 18, 21, 23, 25},
        "D22": {1, 2, 3, 4, 6, 12, 15, 21, 23, 24, 25, 26},
    },
    16: {
        "D11": {6, 7, 8, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 23, 24, 25},
        "D12": {0, 1, 2, 3, 8, 11, 12, 13, 15, 18, 20, 21, 24, 27, 29},
        "D22": {1, 2, 3, 4, 5, 9, 13, 18, 22, 26, 27, 28, 29, 30},
    },
    18: {
        "D11": {6, 8, 9, 10, 11, 13, 14, 15, 17, 18, 20, 21, 22, 24, 25, 26, 27, 29},
        "D12": {0, 1, 2, 3, 8, 9, 10, 12, 14, 15, 19, 20, 22, 24, 25, 28, 32},
        "D22": {1, 2, 3, 4, 5, 7, 12, 16, 19, 23, 28, 30, 31, 32, 33, 34},
    },
    20: {
        "D11": {7, 8, 9, 10, 12, 13, 14, 16, 18, 19, 20, 21, 23, 25, 26, 27, 29, 30, 31, 32},
        "D12": {0, 1, 2, 5, 6, 8, 12, 14, 15, 17, 18, 22, 23, 25, 27, 28, 33, 36, 37},
        "D22": {1, 2, 3, 4, 5, 6, 11, 15, 17, 22, 24, 28, 33, 34, 35, 36, 37, 38},
    },
}

# Load n=22 solution
n22_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "solution_n22.json")
if os.path.exists(n22_path):
    with open(n22_path) as f:
        data = json.load(f)
    KNOWN[22] = {
        "D11": set(data["parameters"]["D11"]),
        "D12": set(data["parameters"]["D12"]),
        "D22": set(data["parameters"]["D22"]),
    }


# ============================================================================
# Helper functions
# ============================================================================

def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0: return False
        i += 6
    return True

def is_prime_power(n):
    """Check if n is a prime power, return (p, k) if so, else None."""
    if n < 2: return None
    for p in range(2, int(n**0.5) + 2):
        if n % p == 0:
            k = 0
            val = n
            while val % p == 0:
                val //= p
                k += 1
            if val == 1:
                return (p, k)
            else:
                return None
    # n is prime
    return (n, 1)

def factorize(n):
    """Return prime factorization as list of (p, e)."""
    factors = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            e = 0
            while n % d == 0:
                n //= d
                e += 1
            factors.append((d, e))
        d += 1
    if n > 1:
        factors.append((n, 1))
    return factors

def euler_phi(n):
    result = n
    for p, _ in factorize(n):
        result = result * (p - 1) // p
    return result

def primitive_root(p):
    """Find a primitive root mod p (p must be prime)."""
    if p == 2: return 1
    phi = p - 1
    factors = factorize(phi)
    for g in range(2, p):
        ok = True
        for (q, _) in factors:
            if pow(g, phi // q, p) == 1:
                ok = False
                break
        if ok:
            return g
    return None

def quadratic_residues(m):
    """Return set of quadratic residues mod m (nonzero)."""
    qr = set()
    for x in range(1, m):
        qr.add((x * x) % m)
    qr.discard(0)
    return qr

def units_mod(m):
    """Return the set of units in Z_m."""
    return {x for x in range(1, m) if math.gcd(x, m) == 1}


# ============================================================================
# 1. SPECTRAL / FOURIER ANALYSIS
# ============================================================================

def character_sum(S, k, m):
    """Compute chi_k(S) = sum_{s in S} omega^{ks} where omega = e^{2*pi*i/m}."""
    omega = cmath.exp(2j * cmath.pi / m)
    return sum(omega ** (k * s) for s in S)

def spectral_analysis(n_param, D11, m):
    """Compute Fourier spectrum of D11 as a subset of Z_m."""
    results = {}
    for k in range(m):
        chi = character_sum(D11, k, m)
        results[k] = {
            "chi": chi,
            "magnitude": abs(chi),
            "magnitude_sq": abs(chi) ** 2,
        }

    # chi_0 = |D11|
    mags_sq = [results[k]["magnitude_sq"] for k in range(1, m)]
    avg_mag_sq = sum(mags_sq) / len(mags_sq)
    max_mag_sq = max(mags_sq)
    min_mag_sq = min(mags_sq)

    # For a perfect (m, |D11|, lambda)-difference set: all |chi_k|^2 = |D11| - lambda for k != 0
    # Flatness ratio: max/min of |chi_k|^2 for k != 0
    flatness = max_mag_sq / min_mag_sq if min_mag_sq > 0.001 else float('inf')

    return {
        "n": n_param,
        "m": m,
        "|D11|": len(D11),
        "chi_0": len(D11),
        "avg_|chi_k|^2 (k>0)": avg_mag_sq,
        "max_|chi_k|^2": max_mag_sq,
        "min_|chi_k|^2": min_mag_sq,
        "flatness_ratio": flatness,
        "spectrum": results,
    }

def run_spectral_analysis():
    """Run spectral analysis on all known D11 sets."""
    print("\n" + "=" * 80)
    print("SPECTRAL / FOURIER ANALYSIS OF D11 SETS")
    print("=" * 80)

    all_results = {}
    for n_param in sorted(KNOWN.keys()):
        m = 2 * n_param - 1
        D11 = KNOWN[n_param]["D11"]
        res = spectral_analysis(n_param, D11, m)
        all_results[n_param] = res

        print(f"\nn={n_param}, m={m}, |D11|={len(D11)}")
        print(f"  chi_0 = |D11| = {len(D11)}")
        print(f"  avg |chi_k|^2 (k>0) = {res['avg_|chi_k|^2 (k>0)']:.4f}")
        print(f"  max |chi_k|^2 = {res['max_|chi_k|^2']:.4f}")
        print(f"  min |chi_k|^2 = {res['min_|chi_k|^2']:.4f}")
        print(f"  flatness ratio (max/min) = {res['flatness_ratio']:.4f}")

        # Parseval check: sum |chi_k|^2 = m * |D11|
        total = sum(res["spectrum"][k]["magnitude_sq"] for k in range(m))
        print(f"  Parseval check: sum|chi_k|^2 = {total:.2f}, m*|D11| = {m * len(D11)}")

        # If D11 were a perfect difference set with lambda, then |chi_k|^2 = |D11| - lambda for all k>0
        # and |D11|(|D11|-1) = lambda(m-1)
        # So lambda = |D11|(|D11|-1)/(m-1) and |chi_k|^2 = |D11| - lambda
        lambda_diff = len(D11) * (len(D11) - 1) / (m - 1) if m > 1 else 0
        ideal_chi_sq = len(D11) - lambda_diff
        print(f"  If perfect difference set: lambda = {lambda_diff:.4f}, |chi_k|^2 = {ideal_chi_sq:.4f}")

        # Print top-5 largest and smallest |chi_k|^2
        sorted_k = sorted(range(1, m), key=lambda k: res["spectrum"][k]["magnitude_sq"], reverse=True)
        print(f"  Top-5 |chi_k|^2: ", end="")
        for k in sorted_k[:5]:
            print(f"k={k}:{res['spectrum'][k]['magnitude_sq']:.3f} ", end="")
        print()
        print(f"  Bottom-5 |chi_k|^2: ", end="")
        for k in sorted_k[-5:]:
            print(f"k={k}:{res['spectrum'][k]['magnitude_sq']:.3f} ", end="")
        print()

    return all_results


# ============================================================================
# 2. CYCLOTOMIC CLASS ANALYSIS
# ============================================================================

def cyclotomic_analysis(n_param, D11, D12, m):
    """Classify D11 and D12 elements by quadratic residue status and cyclotomic classes."""
    results = {}

    qr = quadratic_residues(m)
    qnr = set(range(1, m)) - qr  # quadratic non-residues

    D11_qr = D11 & qr
    D11_qnr = D11 & qnr
    D12_qr = (D12 - {0}) & qr
    D12_qnr = (D12 - {0}) & qnr

    results["m"] = m
    results["is_prime"] = is_prime(m)
    results["m_mod4"] = m % 4
    results["|QR|"] = len(qr)
    results["|QNR|"] = len(qnr)
    results["|D11 ∩ QR|"] = len(D11_qr)
    results["|D11 ∩ QNR|"] = len(D11_qnr)
    results["|D12\\{0} ∩ QR|"] = len(D12_qr)
    results["|D12\\{0} ∩ QNR|"] = len(D12_qnr)
    results["D11_qr_fraction"] = len(D11_qr) / len(D11) if len(D11) > 0 else 0
    results["D12_qr_fraction"] = len(D12_qr) / (len(D12) - (1 if 0 in D12 else 0)) if len(D12) > (1 if 0 in D12 else 0) else 0

    # For prime m: compute higher-order cyclotomic classes
    if is_prime(m) and m > 2:
        g = primitive_root(m)
        results["primitive_root"] = g

        # Index map: x = g^{index(x)}
        index_map = {}
        val = 1
        for i in range(m - 1):
            index_map[val] = i
            val = (val * g) % m

        # Classify D11 elements by index mod various divisors
        D11_indices = sorted([index_map[x] for x in D11 if x in index_map])
        D12_indices = sorted([index_map[x] for x in D12 if x in index_map and x != 0])

        results["D11_indices_mod2"] = [i % 2 for i in D11_indices]
        results["D11_index_mod2_counts"] = {
            0: sum(1 for i in D11_indices if i % 2 == 0),
            1: sum(1 for i in D11_indices if i % 2 == 1),
        }

        # Check if D11 = union of specific cosets of index subgroup
        for divisor in [2, 3, 4, 6]:
            if (m - 1) % divisor == 0:
                classes = defaultdict(set)
                for x in range(1, m):
                    classes[index_map[x] % divisor].add(x)

                # Which classes does D11 cover?
                d11_coverage = {}
                for c in range(divisor):
                    d11_coverage[c] = len(D11 & classes[c])
                results[f"D11_cyclotomic_{divisor}_coverage"] = d11_coverage

                d12_coverage = {}
                for c in range(divisor):
                    d12_coverage[c] = len((D12 - {0}) & classes[c])
                results[f"D12_cyclotomic_{divisor}_coverage"] = d12_coverage

    return results

def run_cyclotomic_analysis():
    """Run cyclotomic analysis on all known constructions."""
    print("\n" + "=" * 80)
    print("CYCLOTOMIC CLASS ANALYSIS")
    print("=" * 80)

    all_results = {}
    for n_param in sorted(KNOWN.keys()):
        m = 2 * n_param - 1
        D11 = KNOWN[n_param]["D11"]
        D12 = KNOWN[n_param]["D12"]
        res = cyclotomic_analysis(n_param, D11, D12, m)
        all_results[n_param] = res

        print(f"\nn={n_param}, m={m} ({'prime' if res['is_prime'] else 'composite'}), m mod 4 = {res['m_mod4']}")
        print(f"  |QR| = {res['|QR|']}, |QNR| = {res['|QNR|']}")
        print(f"  |D11| = {len(D11)}: {res['|D11 ∩ QR|']} QR + {res['|D11 ∩ QNR|']} QNR (QR fraction: {res['D11_qr_fraction']:.3f})")
        d12_no0 = len(D12) - (1 if 0 in D12 else 0)
        d12_qr_count = res['|D12\\{0} ∩ QR|']
        d12_qnr_count = res['|D12\\{0} ∩ QNR|']
        d12_qr_frac = res['D12_qr_fraction']
        print(f"  |D12 \\ {{0}}| = {d12_no0}: {d12_qr_count} QR + {d12_qnr_count} QNR (QR fraction: {d12_qr_frac:.3f})")

        if res["is_prime"] and m > 2:
            print(f"  Primitive root g = {res.get('primitive_root', '?')}")
            if "D11_index_mod2_counts" in res:
                print(f"  D11 index mod 2: {res['D11_index_mod2_counts']}")
            for divisor in [2, 3, 4, 6]:
                key = f"D11_cyclotomic_{divisor}_coverage"
                if key in res:
                    print(f"  D11 cyclotomic-{divisor} class coverage: {dict(res[key])}")
                key12 = f"D12_cyclotomic_{divisor}_coverage"
                if key12 in res:
                    print(f"  D12 cyclotomic-{divisor} class coverage: {dict(res[key12])}")

    return all_results


# ============================================================================
# 3. D12 STRUCTURE ANALYSIS
# ============================================================================

def d12_structure_analysis(n_param, D11, D12, D22, m):
    """Analyze D12's relationship to D11 and D22."""
    results = {}

    # Basic
    results["|D11 ∩ D12|"] = len(D11 & D12)
    results["|D22 ∩ D12|"] = len(D22 & D12)
    results["0 in D12"] = 0 in D12

    # D12 transposed
    D12T = {(-x) % m for x in D12}
    results["D12 == D12^T"] = (D12 == D12T)
    results["|D12 ∩ D12^T|"] = len(D12 & D12T)
    results["|D12 \\ D12^T|"] = len(D12 - D12T)

    # Check if D12 = D11 + c for some constant
    for c in range(m):
        shifted = {(x + c) % m for x in D11}
        if shifted == D12:
            results["D12 = D11 + c"] = c
            break
    else:
        results["D12 = D11 + c"] = None

    # Check if D12 = D22 + c for some constant
    for c in range(m):
        shifted = {(x + c) % m for x in D22}
        if shifted == D12:
            results["D12 = D22 + c"] = c
            break
    else:
        results["D12 = D22 + c"] = None

    # Check if D12 = multiplier * D11 for some multiplier
    for t in range(1, m):
        if math.gcd(t, m) == 1:
            mult = {(t * x) % m for x in D11}
            if mult == D12:
                results["D12 = t * D11"] = t
                break
    else:
        results["D12 = t * D11"] = None

    # Intersection numbers |D11 ∩ (D12 + k)| for various k
    intersections = {}
    for k in range(m):
        shifted_D12 = {(x + k) % m for x in D12}
        intersections[k] = len(D11 & shifted_D12)
    results["D11_D12_intersection_numbers"] = intersections

    # |D12 ∩ (D12 + k)| (autocorrelation of D12)
    autocorr = {}
    for k in range(m):
        shifted = {(x + k) % m for x in D12}
        autocorr[k] = len(D12 & shifted)
    results["D12_autocorrelation"] = autocorr

    # Check: is D12 determined by some "threshold" on D11's Delta function?
    # I.e., is there a function f(d) computed from D11 that determines membership in D12?

    # Compute Delta(D11, D11, d) for each d
    d11_delta = {}
    for d in range(m):
        d11_delta[d] = Delta(D11, D11, d, m)

    results["Delta_D11_D11"] = d11_delta

    # Check if D12 = {d : Delta(D11,D11,d) >= threshold} for some threshold
    for threshold in range(m):
        candidate = {d for d in range(m) if d11_delta[d] >= threshold}
        if candidate == D12:
            results["D12 = {d: Delta(D11,D11,d) >= threshold}"] = threshold
            break
    else:
        results["D12 = {d: Delta(D11,D11,d) >= threshold}"] = None

    # Check Sigma(D11, D11, d) values for D12 vs complement
    sigma_vals = {}
    for d in range(m):
        sigma_vals[d] = Sigma(D11, D11, d, m)
    results["Sigma_D11_D11"] = sigma_vals

    # sigma + delta combined
    combined = {}
    for d in range(m):
        combined[d] = d11_delta[d] + sigma_vals[d]
    results["Delta+Sigma_D11_D11"] = combined

    # Check if D12 correlates with Delta+Sigma threshold
    for threshold in range(2 * m):
        candidate = {d for d in range(m) if combined[d] >= threshold}
        if candidate == D12:
            results["D12 = {d: Delta+Sigma >= threshold}"] = threshold
            break
    else:
        results["D12 = {d: Delta+Sigma >= threshold}"] = None

    return results

def run_d12_analysis():
    """Run D12 structure analysis on all known constructions."""
    print("\n" + "=" * 80)
    print("D12 STRUCTURE ANALYSIS")
    print("=" * 80)

    all_results = {}
    for n_param in sorted(KNOWN.keys()):
        m = 2 * n_param - 1
        D11 = KNOWN[n_param]["D11"]
        D12 = KNOWN[n_param]["D12"]
        D22 = KNOWN[n_param]["D22"]
        res = d12_structure_analysis(n_param, D11, D12, D22, m)
        all_results[n_param] = res

        print(f"\nn={n_param}, m={m}")
        print(f"  |D11 cap D12| = {res['|D11 ∩ D12|']}")
        print(f"  |D22 cap D12| = {res['|D22 ∩ D12|']}")
        print(f"  0 in D12: {res['0 in D12']}")
        d12_eq_d12t = res['D12 == D12^T']
        d12_cap_d12t = res['|D12 ∩ D12^T|']
        d12_minus_d12t = res['|D12 \\ D12^T|']
        print(f"  D12 == D12^T: {d12_eq_d12t}")
        print(f"  |D12 cap D12^T| = {d12_cap_d12t}")
        print(f"  |D12 minus D12^T| = {d12_minus_d12t}")
        print(f"  D12 = D11 + c: {res['D12 = D11 + c']}")
        print(f"  D12 = D22 + c: {res['D12 = D22 + c']}")
        print(f"  D12 = t * D11: {res['D12 = t * D11']}")
        delta_thresh = res['D12 = {d: Delta(D11,D11,d) >= threshold}']
        ds_thresh = res['D12 = {d: Delta+Sigma >= threshold}']
        print(f"  D12 = {{d: Delta(D11,D11,d) >= threshold}}: {delta_thresh}")
        print(f"  D12 = {{d: Delta+Sigma >= threshold}}: {ds_thresh}")

        # Print Delta(D11,D11,d) for d in D12 vs d not in D12
        d11_d = res["Delta_D11_D11"]
        in_d12 = sorted([d11_d[d] for d in D12])
        not_d12 = sorted([d11_d[d] for d in range(m) if d not in D12])
        print(f"  Delta(D11,D11,d) for d in D12: {in_d12}")
        print(f"  Delta(D11,D11,d) for d not in D12: {not_d12}")

        # Print intersection numbers distribution
        ints = res["D11_D12_intersection_numbers"]
        int_vals = sorted(set(ints.values()))
        print(f"  |D11 ∩ (D12+k)| values: {int_vals}")
        for v in int_vals:
            ks = [k for k in range(m) if ints[k] == v]
            print(f"    value {v}: {len(ks)} occurrences, k={ks[:10]}{'...' if len(ks) > 10 else ''}")

    return all_results


# ============================================================================
# 4. CATEGORIZE ALL n <= 50
# ============================================================================

def categorize_n_values():
    """Categorize all n from 3 to 50 by which proof method applies."""
    print("\n" + "=" * 80)
    print("CATEGORIZATION OF ALL n <= 50 BY PROOF METHOD")
    print("=" * 80)

    categories = {
        "paley_1mod4": [],        # 2n-1 is prime power ≡ 1 (mod 4) -- Paley construction
        "prime_3mod4": [],        # 2n-1 is prime ≡ 3 (mod 4) -- SA search
        "prime_power_3mod4": [],  # 2n-1 is prime power (not prime) ≡ 3 (mod 4)
        "composite": [],          # 2n-1 is composite, not a prime power
        "solved_sa": [],          # Already solved by SA
    }

    print(f"\n{'n':>3} {'m=2n-1':>6} {'factorization':>20} {'m mod 4':>7} {'type':>25} {'status':>15}")
    print("-" * 85)

    for n in range(3, 51):
        m = 2 * n - 1
        factors = factorize(m)
        factor_str = " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in factors)
        pp = is_prime_power(m)

        m_mod4 = m % 4

        if pp is not None:
            p, k = pp
            if m_mod4 == 1:
                cat = "paley_1mod4"
                type_str = f"prime power ≡ 1(4)"
                if k == 1:
                    type_str = f"prime ≡ 1(4)"
            else:  # m_mod4 == 3
                if k == 1:
                    cat = "prime_3mod4"
                    type_str = f"prime ≡ 3(4)"
                else:
                    cat = "prime_power_3mod4"
                    type_str = f"p^{k} ≡ 3(4)"
        else:
            cat = "composite"
            type_str = "composite"

        status = "proven" if n in KNOWN else ("paley" if cat == "paley_1mod4" else "OPEN")
        if cat == "prime_3mod4" and n <= 22:
            status = "proven (SA)"

        categories[cat].append(n)
        print(f"{n:3d} {m:6d} {factor_str:>20} {m_mod4:7d} {type_str:>25} {status:>15}")

    print(f"\n--- Summary ---")
    print(f"Paley (prime power ≡ 1 mod 4): {len(categories['paley_1mod4'])} cases: {categories['paley_1mod4']}")
    print(f"Prime ≡ 3 mod 4: {len(categories['prime_3mod4'])} cases: {categories['prime_3mod4']}")
    print(f"Prime power ≡ 3 mod 4 (not prime): {len(categories['prime_power_3mod4'])} cases: {categories['prime_power_3mod4']}")
    print(f"Composite: {len(categories['composite'])} cases: {categories['composite']}")

    # Identify the hard cases
    hard = [n for n in categories["composite"]]
    print(f"\nHARD CASES (composite m, not prime power): n = {hard}")
    for n in hard:
        m = 2 * n - 1
        factors = factorize(m)
        factor_str = " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in factors)
        print(f"  n={n}, m={m} = {factor_str}")

    return categories


# ============================================================================
# 5. COMPOSITE m=45 ANALYSIS (n=23)
# ============================================================================

def composite_m45_analysis():
    """Analyze Z_45 for constructing D11/D12 for n=23."""
    print("\n" + "=" * 80)
    print("COMPOSITE m=45 ANALYSIS (n=23)")
    print("=" * 80)

    m = 45
    n = 23

    # Z_45 ≅ Z_9 × Z_5 via CRT
    # CRT map: x -> (x mod 9, x mod 5)
    print(f"\nZ_{m} ≅ Z_9 x Z_5")
    print(f"CRT mapping x -> (x mod 9, x mod 5):")

    crt_map = {}
    inv_crt = {}
    for x in range(m):
        pair = (x % 9, x % 5)
        crt_map[x] = pair
        inv_crt[pair] = x

    # Print mapping table
    print(f"\n  x -> (x%9, x%5):")
    for x in range(m):
        if x % 9 == 0:
            print(f"  ", end="")
        print(f"{x:2d}->({crt_map[x][0]},{crt_map[x][1]}) ", end="")
        if (x + 1) % 9 == 0:
            print()

    # Units of Z_45
    units = units_mod(m)
    print(f"\nUnits of Z_{m}: {sorted(units)} ({len(units)} = phi(45) = {euler_phi(m)})")

    # Quadratic residues of Z_45
    qr = quadratic_residues(m)
    print(f"Quadratic residues mod {m}: {sorted(qr)} ({len(qr)} elements)")

    # Quadratic residues in each component
    qr_9 = quadratic_residues(9)
    qr_5 = quadratic_residues(5)
    print(f"QR mod 9: {sorted(qr_9)}")
    print(f"QR mod 5: {sorted(qr_5)}")

    # For n=23: |D11| = n or n+1 (need to figure out the right sizes)
    # Actually from the pattern: for m=2n-1, |D11|=(m+1)/2 or (m-1)/2
    # m=45, so |D11| = 22 or 23, |D12| = 22 or 23, |D22| = 22 or 23
    # With d1 + d2 = N-2+1 = 4n-3 = 89... no
    # d1 = |D11| + |D12|, d2 = |D22| + |D12|
    # For valid construction: typically d1 = m = 2n-1 or d1 = m-2

    # From pattern in known constructions:
    # For n even (n=6,8,10): |D11| = n-2 = 2*(n/2-1), d1 = 2n-2 = N/2-1, etc.
    # For n >= 12 (even): |D11| = n, |D12| = n-1, |D22| = n-2
    # Wait, let me check the actual sizes

    print(f"\n--- Expected sizes for n={n}, m={m} ---")
    # Pattern from known data:
    # n=12: |D11|=12, |D12|=11, |D22|=10 -> d1=23=m, d2=21=m-2
    # n=14: |D11|=14, |D12|=13, |D22|=12 -> d1=27=m, d2=25=m-2
    # n=16: |D11|=16, |D12|=15, |D22|=14 -> d1=31=m, d2=29=m-2
    # n=18: |D11|=18, |D12|=17, |D22|=16 -> d1=35=m, d2=33=m-2
    # n=20: |D11|=20, |D12|=19, |D22|=18 -> d1=39=m, d2=37=m-2
    # n=22: |D11|=21, |D12|=21, |D22|=20 -> d1=42, d2=41
    # Hmm, n=22 breaks the pattern slightly

    for np in sorted(KNOWN.keys()):
        mp = 2 * np - 1
        d11s = len(KNOWN[np]["D11"])
        d12s = len(KNOWN[np]["D12"])
        d22s = len(KNOWN[np]["D22"])
        d1 = d11s + d12s
        d2 = d22s + d12s
        print(f"  n={np}: |D11|={d11s}, |D12|={d12s}, |D22|={d22s}, d1={d1}, d2={d2}, m={mp}")
        print(f"         |D11|/m={d11s/mp:.4f}, |D12|/m={d12s/mp:.4f}, |D22|/m={d22s/mp:.4f}")

    # For n=23, m=45:
    # Expected: |D11| ~ n = 23, |D12| ~ n-1 = 22, |D22| ~ n-2 = 21
    # Or |D11| ~ (m+1)/2 = 23, etc.
    d11_size = 23
    d12_size = 22
    d22_size = 21  # = m - 1 - d11_size + 1? No, D22 = complement of D11 would give |D22| = m-1-|D11| = 21

    print(f"\n  Expected for n=23: |D11|={d11_size}, |D12|={d12_size}, |D22|={d22_size}")
    print(f"  d1={d11_size+d12_size}, d2={d22_size+d12_size}")

    # Product structure analysis
    print(f"\n--- Product structure of potential D11 ---")
    print(f"If D11 = A x B for A ⊂ Z_9, B ⊂ Z_5:")
    print(f"  |D11| = |A| * |B|")
    print(f"  Need |D11| = 22 or 23")
    print(f"  Possible: |A|*|B| ∈ {{22, 23}}")
    print(f"  But 22 = 2*11 or 22*1 (11 > 8, so no), 23 = 23*1 (no product works)")
    print(f"  --> D11 CANNOT be a pure product A x B")

    # Instead, D11 could be a union of products
    print(f"\n  D11 as union: e.g., A1xB1 ∪ A2xB2")
    print(f"  Possible decompositions for |D11| = 23:")
    for a1 in range(1, 9):
        for b1 in range(1, 6):
            rem = 23 - a1 * b1
            if rem > 0:
                for a2 in range(1, 9):
                    for b2 in range(1, 6):
                        if a2 * b2 == rem:
                            # This is just one possible decomposition
                            pass

    # More useful: analyze what QR structure looks like in Z_45
    print(f"\n--- Quadratic residues in CRT coordinates ---")
    print(f"QR mod 45 in (Z_9, Z_5) coordinates:")
    for x in sorted(qr):
        print(f"  {x:2d} -> ({x%9}, {x%5})", end="")
        print(f"  QR9={'Y' if x%9 in qr_9 else 'N'}, QR5={'Y' if x%5 in qr_5 else 'N'}")

    # Multiplier orbits under units
    print(f"\n--- Multiplier orbits in Z_{m} ---")
    # Under multiplication by units, which orbits does each element belong to?
    visited = set()
    orbits = []
    for x in range(m):
        if x not in visited:
            orbit = set()
            for u in units:
                orbit.add((u * x) % m)
            orbits.append(sorted(orbit))
            visited |= orbit

    print(f"Number of orbits under units: {len(orbits)}")
    for i, orb in enumerate(orbits):
        print(f"  Orbit {i}: {orb} (size {len(orb)})")

    # Negation orbits: {x, -x}
    print(f"\n--- Negation orbits (symmetric pairs) ---")
    neg_orbits = []
    neg_visited = set()
    for x in range(m):
        if x not in neg_visited:
            pair = {x, (-x) % m}
            neg_orbits.append(sorted(pair))
            neg_visited |= pair

    print(f"Number of negation orbits: {len(neg_orbits)} (0 is self-paired)")
    sym_pairs = [(a, b) for a, b in [(o[0], o[-1]) for o in neg_orbits] if a != b]
    print(f"Non-trivial symmetric pairs: {len(sym_pairs)}")

    # For D11 symmetric with |D11| = 22:
    # Need to pick 11 pairs from the 22 non-trivial pairs (total 44 nonzero elements, 22 pairs)
    total_pairs = len([o for o in neg_orbits if len(o) == 2])
    print(f"Total symmetric pairs: {total_pairs}")
    print(f"For |D11|=22, need to pick 11 of {total_pairs} pairs")
    print(f"For |D11|=23 (if m odd and 0 not in D11): need to pick 11.5 pairs -- impossible unless...")

    # Wait - m=45 has no element equal to its own negation except 0
    # because 2x ≡ 0 mod 45 means x = 0 or x = 45/gcd(2,45) ... gcd(2,45)=1 so only x=0
    # So all nonzero elements pair up: 44/2 = 22 pairs
    # |D11| must be even (since D11 is symmetric). So |D11| ∈ {22, 24, ...}
    print(f"\n  Since D11 must be symmetric and m=45 is odd, |D11| must be even.")
    print(f"  Most likely |D11| = 22 (= n-1), |D22| = 22 (= complement has 22 elements)")
    print(f"  Then |D12| should be ~ 22 or 23")

    # Alternative: |D11| = 24 gives |D22| = 20
    print(f"\n  Alternative: |D11| = 24, |D22| = 20, and |D12| ~ 21 or 22")

    return {}


# ============================================================================
# 6. DELTA FLUCTUATION BOUNDS
# ============================================================================

def delta_fluctuation_analysis():
    """Analyze Delta fluctuation patterns across all known constructions."""
    print("\n" + "=" * 80)
    print("DELTA FLUCTUATION BOUNDS")
    print("=" * 80)

    results = []

    for n_param in sorted(KNOWN.keys()):
        m = 2 * n_param - 1
        D11 = KNOWN[n_param]["D11"]
        D12 = KNOWN[n_param]["D12"]
        D22 = KNOWN[n_param]["D22"]

        # V1V1 Delta values
        v11_deltas = {}
        for d in range(1, m):
            v11_deltas[d] = Delta(D11, D11, d, m) + Delta(D12, D12, d, m)

        all_v11 = list(v11_deltas.values())
        avg_v11 = sum(all_v11) / len(all_v11)
        max_v11 = max(all_v11)
        min_v11 = min(all_v11)
        range_v11 = max_v11 - min_v11

        # Deviation of red edges from average
        red_v11 = [v11_deltas[d] for d in range(1, m) if d in D11]
        blue_v11 = [v11_deltas[d] for d in range(1, m) if d not in D11]

        # V2V2 Delta values
        D12T = {(-x) % m for x in D12}
        v22_deltas = {}
        for d in range(1, m):
            v22_deltas[d] = Delta(D22, D22, d, m) + Delta(D12T, D12T, d, m)

        all_v22 = list(v22_deltas.values())
        avg_v22 = sum(all_v22) / len(all_v22)

        # V1V2 Delta values
        v12_deltas = {}
        for d in range(m):
            v12_deltas[d] = Sigma(D11, D12, d, m) + Delta(D12, D22, d, m)

        all_v12 = list(v12_deltas.values())
        avg_v12 = sum(all_v12) / len(all_v12)

        # Compute standard deviation
        import statistics
        std_v11 = statistics.stdev(all_v11) if len(all_v11) > 1 else 0
        std_v22 = statistics.stdev(all_v22) if len(all_v22) > 1 else 0
        std_v12 = statistics.stdev(all_v12) if len(all_v12) > 1 else 0

        red_threshold = n_param - 2

        # Key metric: how many standard deviations below avg are the red values?
        max_red_v11 = max(red_v11)
        gap_v11 = avg_v11 - max_red_v11
        gap_in_stds_v11 = gap_v11 / std_v11 if std_v11 > 0 else 0

        # Normalized fluctuation: range / sqrt(m)
        norm_range_v11 = range_v11 / math.sqrt(m)
        norm_std_v11 = std_v11 / math.sqrt(m)

        row = {
            "n": n_param, "m": m,
            "|D11|": len(D11), "|D12|": len(D12), "|D22|": len(D22),
            "avg_V11": avg_v11, "std_V11": std_v11, "range_V11": range_v11,
            "avg_V22": avg_v22, "std_V22": std_v22,
            "avg_V12": avg_v12, "std_V12": std_v12,
            "max_red_V11": max_red_v11, "red_threshold": red_threshold,
            "gap_V11": gap_v11, "gap_in_stds_V11": gap_in_stds_v11,
            "norm_range_V11": norm_range_v11, "norm_std_V11": norm_std_v11,
        }
        results.append(row)

        print(f"\nn={n_param}, m={m}")
        print(f"  V1V1: avg={avg_v11:.3f}, std={std_v11:.3f}, range={range_v11}, norm_range={norm_range_v11:.3f}")
        print(f"    Red max={max_red_v11}, threshold={red_threshold}, gap={gap_v11:.3f}, gap/std={gap_in_stds_v11:.3f}")
        print(f"  V2V2: avg={avg_v22:.3f}, std={std_v22:.3f}")
        print(f"  V1V2: avg={avg_v12:.3f}, std={std_v12:.3f}")

    # Summary table
    print(f"\n--- Summary Table ---")
    print(f"{'n':>3} {'m':>3} {'avg_V11':>8} {'std_V11':>8} {'range_V11':>9} {'norm_range':>10} {'norm_std':>8} {'gap/std':>7}")
    print("-" * 65)
    for row in results:
        print(f"{row['n']:3d} {row['m']:3d} {row['avg_V11']:8.3f} {row['std_V11']:8.3f} "
              f"{row['range_V11']:9d} {row['norm_range_V11']:10.4f} {row['norm_std_V11']:8.4f} "
              f"{row['gap_in_stds_V11']:7.3f}")

    # Check if norm_std is roughly constant
    norm_stds = [r["norm_std_V11"] for r in results]
    print(f"\nNormalized std (std/sqrt(m)) values: {[f'{v:.4f}' for v in norm_stds]}")
    print(f"Mean normalized std: {sum(norm_stds)/len(norm_stds):.4f}")
    print(f"This suggests std ≈ {sum(norm_stds)/len(norm_stds):.3f} * sqrt(m)")

    return results


# ============================================================================
# 7. DEGREE SEQUENCE AND SIZE PATTERN ANALYSIS
# ============================================================================

def size_pattern_analysis():
    """Analyze the pattern of |D11|, |D12|, |D22| sizes across all n."""
    print("\n" + "=" * 80)
    print("SIZE PATTERN ANALYSIS")
    print("=" * 80)

    print(f"\n{'n':>3} {'m':>3} {'|D11|':>5} {'|D12|':>5} {'|D22|':>5} {'d1':>4} {'d2':>4} "
          f"{'|D11|/m':>7} {'|D12|/m':>7} {'d1=m?':>6} {'D22=comp?':>9}")
    print("-" * 75)

    for n_param in sorted(KNOWN.keys()):
        m = 2 * n_param - 1
        D11 = KNOWN[n_param]["D11"]
        D12 = KNOWN[n_param]["D12"]
        D22 = KNOWN[n_param]["D22"]
        d1 = len(D11) + len(D12)
        d2 = len(D22) + len(D12)

        all_nonzero = set(range(1, m))
        is_comp = D22 == all_nonzero - D11

        print(f"{n_param:3d} {m:3d} {len(D11):5d} {len(D12):5d} {len(D22):5d} "
              f"{d1:4d} {d2:4d} {len(D11)/m:7.4f} {len(D12)/m:7.4f} "
              f"{'YES' if d1 == m else 'no':>6} {'YES' if is_comp else 'no':>9}")

    # Verify complement property
    print(f"\nComplement property D22 = {{1,...,m-1}} \\ D11:")
    for n_param in sorted(KNOWN.keys()):
        m = 2 * n_param - 1
        D11 = KNOWN[n_param]["D11"]
        D22 = KNOWN[n_param]["D22"]
        comp = set(range(1, m)) - D11
        print(f"  n={n_param}: D22 == comp(D11): {D22 == comp}")

    # Check symmetry
    print(f"\nSymmetry check:")
    for n_param in sorted(KNOWN.keys()):
        m = 2 * n_param - 1
        D11 = KNOWN[n_param]["D11"]
        D12 = KNOWN[n_param]["D12"]
        D22 = KNOWN[n_param]["D22"]

        D11_sym = all((-x) % m in D11 for x in D11)
        D22_sym = all((-x) % m in D22 for x in D22)
        D12_sym = all((-x) % m in D12 for x in D12)

        print(f"  n={n_param}: D11 symmetric: {D11_sym}, D22 symmetric: {D22_sym}, D12 symmetric: {D12_sym}")


# ============================================================================
# 8. CROSS-BLOCK (V1V2) TIGHTNESS ANALYSIS
# ============================================================================

def v1v2_tightness():
    """Analyze how tight the V1V2 constraints are -- range is always 1!"""
    print("\n" + "=" * 80)
    print("V1V2 TIGHTNESS ANALYSIS")
    print("=" * 80)

    for n_param in sorted(KNOWN.keys()):
        m = 2 * n_param - 1
        D11 = KNOWN[n_param]["D11"]
        D12 = KNOWN[n_param]["D12"]
        D22 = KNOWN[n_param]["D22"]

        v12_deltas = {}
        for d in range(m):
            v12_deltas[d] = Sigma(D11, D12, d, m) + Delta(D12, D22, d, m)

        red_v12 = {d: v12_deltas[d] for d in D12}
        blue_v12 = {d: v12_deltas[d] for d in range(m) if d not in D12}

        print(f"\nn={n_param}, m={m}")
        print(f"  V1V2 red values (d in D12): {sorted(red_v12.values())}")
        print(f"  V1V2 blue values (d not in D12): {sorted(blue_v12.values())}")
        print(f"  All V1V2 values: min={min(v12_deltas.values())}, max={max(v12_deltas.values())}")

        # Check: is every red value = threshold?
        red_thresh = n_param - 2
        all_at_thresh = all(v == red_thresh for v in red_v12.values())
        print(f"  All red V1V2 = threshold ({red_thresh}): {all_at_thresh}")


# ============================================================================
# 9. MULTIPLIER EQUIVALENCE ANALYSIS
# ============================================================================

def multiplier_analysis():
    """Check which multipliers preserve D11 (automorphisms of the construction)."""
    print("\n" + "=" * 80)
    print("MULTIPLIER / AUTOMORPHISM ANALYSIS")
    print("=" * 80)

    for n_param in sorted(KNOWN.keys()):
        m = 2 * n_param - 1
        D11 = KNOWN[n_param]["D11"]
        D12 = KNOWN[n_param]["D12"]

        # Find multipliers t such that t*D11 = D11
        auto_d11 = []
        for t in range(1, m):
            if math.gcd(t, m) == 1:
                mult = {(t * x) % m for x in D11}
                if mult == D11:
                    auto_d11.append(t)

        # Find multipliers t such that t*D12 = D12
        auto_d12 = []
        for t in range(1, m):
            if math.gcd(t, m) == 1:
                mult = {(t * x) % m for x in D12}
                if mult == D12:
                    auto_d12.append(t)

        # Find multipliers t such that t preserves both D11 and D12
        auto_both = [t for t in auto_d11 if t in auto_d12]

        print(f"\nn={n_param}, m={m}")
        print(f"  Multipliers fixing D11: {auto_d11} (count: {len(auto_d11)})")
        print(f"  Multipliers fixing D12: {auto_d12} (count: {len(auto_d12)})")
        print(f"  Multipliers fixing both: {auto_both} (count: {len(auto_both)})")

        # For prime m, check if auto_d11 forms a subgroup of (Z/mZ)*
        if is_prime(m) and len(auto_d11) > 1:
            # Check closure
            closed = True
            for a in auto_d11:
                for b in auto_d11:
                    if (a * b) % m not in auto_d11:
                        closed = False
                        break
                if not closed:
                    break
            print(f"  auto_D11 is subgroup: {closed}")
            if closed:
                # What subgroup is it?
                g = primitive_root(m)
                idx_map = {}
                val = 1
                for i in range(m - 1):
                    idx_map[val] = i
                    val = (val * g) % m
                auto_indices = sorted([idx_map[t] for t in auto_d11])
                print(f"  auto_D11 indices (base g={g}): {auto_indices}")
                # GCD of indices
                from math import gcd
                from functools import reduce
                g_idx = reduce(gcd, auto_indices) if auto_indices else 0
                print(f"  GCD of indices: {g_idx}, subgroup index: {(m-1)//len(auto_d11)}")


# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    print("DEEP PATTERN ANALYSIS OF RAMSEY BOOK GRAPH CONSTRUCTIONS")
    print("=" * 80)

    # 0. Size patterns
    size_pattern_analysis()

    # 1. Spectral analysis
    spectral_results = run_spectral_analysis()

    # 2. Cyclotomic analysis
    cyclotomic_results = run_cyclotomic_analysis()

    # 3. D12 structure
    d12_results = run_d12_analysis()

    # 4. Categorization
    categories = categorize_n_values()

    # 5. Composite m=45
    composite_m45_analysis()

    # 6. Delta fluctuation
    delta_results = delta_fluctuation_analysis()

    # 7. V1V2 tightness
    v1v2_tightness()

    # 8. Multiplier analysis
    multiplier_analysis()

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
