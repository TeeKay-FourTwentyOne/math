#!/usr/bin/env python3
"""
SA orbit scanner: sample random symmetric D11 orbits, search for valid D12.
Goal: estimate p_working = fraction of orbits with N(D11) > 0.
"""

import sys, os, time, random, math, json, argparse
from collections import defaultdict
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import Delta, Sigma


def gen_random_symmetric_d11(p):
    """Generate a random symmetric D11 of size (p+1)/2 from {1,...,p-1}."""
    half = (p - 1) // 2
    choose = (p + 1) // 4
    pairs = list(range(1, half + 1))  # pairs are (d, p-d)
    random.shuffle(pairs)
    selected = pairs[:choose]
    d11 = set()
    for d in selected:
        d11.add(d)
        d11.add(p - d)
    return frozenset(d11)


def multiplicative_canonical(d11, p):
    """Return canonical form under multiplication by quadratic residues."""
    # Find a generator of the multiplicative group mod p
    # Use the first quadratic residue multiplier
    best = d11
    # Multiply by all units and take lexicographic min
    for k in range(1, p):
        mapped = frozenset((d * k) % p for d in d11)
        # Discard if 0 appears (shouldn't for d11 ⊂ {1,...,p-1})
        if 0 in mapped:
            continue
        if sorted(mapped) < sorted(best):
            best = mapped
    return best


def fast_cost(d11_set, d12_set, p):
    """Compute total constraint violation cost."""
    n = (p + 1) // 2
    m = p
    d22_set = set(range(1, p)) - set(d11_set)
    d12t_set = {(-x) % p for x in d12_set}

    cost = 0

    for d in range(1, p):
        # V1V1
        A = Delta(d11_set, d11_set, d, m)
        B = Delta(d12_set, d12_set, d, m)
        lam_v1v1 = A + B
        if d in d11_set:
            if lam_v1v1 > n - 2:
                cost += lam_v1v1 - (n - 2)
        else:
            d1 = len(d11_set) + len(d12_set)
            lam_blue = (4 * n - 4) - 2 * d1 + lam_v1v1
            if lam_blue > n - 1:
                cost += lam_blue - (n - 1)

        # V2V2
        C = Delta(d22_set, d22_set, d, m)
        Bt = Delta(d12t_set, d12t_set, d, m)
        lam_v2v2 = C + Bt
        if d in d22_set:
            if lam_v2v2 > n - 2:
                cost += lam_v2v2 - (n - 2)
        else:
            d2 = len(d22_set) + len(d12_set)
            lam_blue = (4 * n - 4) - 2 * d2 + lam_v2v2
            if lam_blue > n - 1:
                cost += lam_blue - (n - 1)

    for d in range(p):
        # V1V2
        S = Sigma(d11_set, d12_set, d, m)
        Dx = Delta(d12_set, d22_set, d, m)
        lam_v1v2 = S + Dx
        if d in d12_set:
            if lam_v1v2 > n - 2:
                cost += lam_v1v2 - (n - 2)
        else:
            d1 = len(d11_set) + len(d12_set)
            d2 = len(d22_set) + len(d12_set)
            lam_blue = (4 * n - 4) - d1 - d2 + lam_v1v2
            if lam_blue > n - 1:
                cost += lam_blue - (n - 1)

    return cost


def sa_search_d12(d11_set, p, max_iter=500000, n_trials=8,
                  T_init=3.0, cooling=0.99999):
    """SA search for valid D12 with fixed D11. Returns D12 if found, else None."""
    n = (p + 1) // 2
    m = p
    d12_size = (p - 1) // 2
    d22_set = set(range(1, p)) - set(d11_set)

    for trial in range(n_trials):
        # Random D12: always include 0, choose (d12_size-1) from {1,...,p-1}
        rest = list(range(1, p))
        random.shuffle(rest)
        d12 = {0} | set(rest[:d12_size - 1])
        not_d12 = set(range(1, p)) - d12

        cost = fast_cost(d11_set, d12, p)
        if cost == 0:
            return d12

        best_cost = cost
        temp = T_init

        for it in range(max_iter):
            # Swap move: exchange one element in D12\{0} with one not in D12
            d12_movable = list(d12 - {0})
            not_d12_list = list(not_d12)
            a = random.choice(d12_movable)
            b = random.choice(not_d12_list)

            d12.remove(a)
            d12.add(b)
            not_d12.remove(b)
            not_d12.add(a)

            new_cost = fast_cost(d11_set, d12, p)

            if new_cost == 0:
                return d12

            delta = new_cost - cost
            if delta <= 0 or random.random() < math.exp(-delta / temp):
                cost = new_cost
                if cost < best_cost:
                    best_cost = cost
            else:
                # Undo
                d12.remove(b)
                d12.add(a)
                not_d12.remove(a)
                not_d12.add(b)

            temp *= cooling

    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--p', type=int, required=True)
    parser.add_argument('--n-orbits', type=int, default=500,
                        help='Number of random orbits to sample')
    parser.add_argument('--sa-trials', type=int, default=12,
                        help='SA trials per orbit')
    parser.add_argument('--sa-iter', type=int, default=500000,
                        help='SA iterations per trial')
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    p = args.p
    n = (p + 1) // 2
    random.seed(args.seed)

    half = (p - 1) // 2
    choose = (p + 1) // 4
    from math import comb
    total_orbits = comb(half, choose) // half

    print(f"SA orbit scan: p={p}, n={n}")
    print(f"Total orbits: {total_orbits}")
    print(f"Sampling {args.n_orbits} orbits, {args.sa_trials} SA trials each, "
          f"{args.sa_iter} iters/trial")
    print(f"Seed: {args.seed}")
    print()

    seen_orbits = set()
    n_tested = 0
    n_working = 0
    t0 = time.time()

    working_d11s = []

    while n_tested < args.n_orbits:
        # Generate random D11 and canonicalize
        d11 = gen_random_symmetric_d11(p)

        # Simple dedup: use sorted tuple as key
        key = tuple(sorted(d11))
        # For orbit dedup, multiply by smallest QR
        # Just use key directly - collisions within same orbit are fine
        # (they'll have same working status)
        if key in seen_orbits:
            continue

        # Try all multipliers to find canonical form (expensive for large p)
        # For sampling purposes, skip canonicalization and accept some orbit repeats
        seen_orbits.add(key)
        n_tested += 1

        result = sa_search_d12(d11, p, max_iter=args.sa_iter,
                               n_trials=args.sa_trials)
        if result is not None:
            n_working += 1
            working_d11s.append(sorted(d11))

        if n_tested % 20 == 0 or result is not None:
            elapsed = time.time() - t0
            pw = n_working / n_tested
            ppw = p * pw
            print(f"  [{n_tested:>4}/{args.n_orbits}] working={n_working}, "
                  f"p_working={pw:.4f}, p*p_working={ppw:.3f} "
                  f"[{elapsed:.0f}s]")

    elapsed = time.time() - t0
    pw = n_working / n_tested if n_tested > 0 else 0
    ppw = p * pw

    print(f"\n{'='*60}")
    print(f"RESULTS: p={p}")
    print(f"  Orbits tested:  {n_tested}")
    print(f"  Working:        {n_working}")
    print(f"  p_working:      {pw:.4f}")
    print(f"  p × p_working:  {ppw:.4f}")
    print(f"  Total time:     {elapsed:.1f}s")
    print(f"{'='*60}")

    # Save results
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    outfile = os.path.join(base_dir, f'enumeration_data/sa_orbit_scan_p{p}.json')
    with open(outfile, 'w') as f:
        json.dump({
            'p': p, 'n': n,
            'total_orbits': total_orbits,
            'orbits_tested': n_tested,
            'working': n_working,
            'p_working': pw,
            'p_times_p_working': ppw,
            'sa_trials': args.sa_trials,
            'sa_iter': args.sa_iter,
            'seed': args.seed,
            'elapsed': elapsed,
            'working_d11s': working_d11s,
        }, f, indent=2)
    print(f"Saved to {outfile}")


if __name__ == '__main__':
    main()
