#!/usr/bin/env python3
"""
Batch computation of omega (D12-orbit count) for all 124 working D11 orbits at p=43.
Uses count_omega_fast (incremental delta updates) for much faster SA.
Runs in parallel across CPU cores. Saves results incrementally.
"""
import json
import subprocess
import sys
import os
import time
import math
from concurrent.futures import ProcessPoolExecutor, as_completed

DIR = os.path.dirname(os.path.abspath(__file__))
BINARY = os.path.join(DIR, "count_omega_fast")
N_TRIALS = 500  # SA trials per orbit
RESULTS_FILE = os.path.join(DIR, "omega_results_fast_p43.jsonl")
FINAL_FILE = os.path.join(DIR, "p43_omega_all.json")
EXHAUSTIVE_FILE = os.path.join(DIR, "p43_exhaustive_results.jsonl")
OLD_RESULTS_FILE = os.path.join(DIR, "omega_results_p43.jsonl")

# Known omega data from proof sketch Section 14 (23 orbits, already verified)
KNOWN_OMEGA = {
    15015: 63, 14837: 10, 15142: 100, 11409: 4, 9404: 324,
    10590: 3, 6530: 11, 11674: 3, 6123: 13, 14644: 18,
    10082: 7, 15669: 365, 15355: 25, 15972: 6, 10176: 66,
    10668: 2, 8359: 26, 8421: 27, 9549: 6, 8777: 22,
    15815: 2, 6664: 1, 6655: 2,
}


def run_one_orbit(orbit_id, d11, n_trials):
    """Run count_omega_fast for a single D11 orbit and parse the result."""
    args = [BINARY, str(n_trials)] + [str(x) for x in d11]
    try:
        result = subprocess.run(args, capture_output=True, text=True, timeout=7200)
        output = result.stdout
        omega = None
        found = None
        for line in output.split('\n'):
            if 'Distinct D12-orbits (omega):' in line:
                omega = int(line.split(':')[1].strip())
            if 'Found valid:' in line:
                found = int(line.split('Found valid:')[1].split('(')[0].strip())
        return {
            'orbit_id': orbit_id,
            'omega': omega,
            'found': found,
            'n_trials': n_trials,
        }
    except subprocess.TimeoutExpired:
        return {'orbit_id': orbit_id, 'omega': None, 'error': 'timeout'}
    except Exception as e:
        return {'orbit_id': orbit_id, 'omega': None, 'error': str(e)}


def load_working_orbits():
    """Load all 124 working orbits from exhaustive results."""
    working = []
    with open(EXHAUSTIVE_FILE) as f:
        for line in f:
            rec = json.loads(line)
            if rec.get('working'):
                working.append(rec)
    return working


def load_existing_results():
    """Load all existing omega results from both old and new results files."""
    done = {}

    # Load known omega from proof sketch
    for oid, omega in KNOWN_OMEGA.items():
        done[oid] = {'orbit_id': oid, 'omega': omega, 'source': 'proof_sketch'}

    # Load old results (successful ones only)
    if os.path.exists(OLD_RESULTS_FILE):
        with open(OLD_RESULTS_FILE) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                rec = json.loads(line)
                if rec.get('omega') is not None:
                    oid = rec['orbit_id']
                    # Only override if we have a higher omega (more trials found more)
                    if oid not in done or rec['omega'] > done[oid].get('omega', 0):
                        done[oid] = rec

    # Load new fast results (successful ones only)
    if os.path.exists(RESULTS_FILE):
        with open(RESULTS_FILE) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                rec = json.loads(line)
                if rec.get('omega') is not None:
                    oid = rec['orbit_id']
                    if oid not in done or rec['omega'] > done[oid].get('omega', 0):
                        done[oid] = rec

    return done


def compute_stats(done, n_working):
    """Compute E[N], E[N^2], ratio from omega data."""
    total_orbits = 16796  # total D11-orbits at p=43
    orbit_size = 86  # 2*43

    omegas = []
    omega_by_id = {}
    for oid, rec in sorted(done.items()):
        if rec.get('omega') is not None:
            omegas.append(rec['omega'])
            omega_by_id[oid] = rec['omega']

    if len(omegas) < n_working:
        print(f"WARNING: only {len(omegas)}/{n_working} orbits have omega data")

    total_N = sum(w * orbit_size for w in omegas)
    total_N2 = sum((w * orbit_size)**2 for w in omegas)

    E_N = total_N / total_orbits
    E_N2 = total_N2 / total_orbits
    ratio = E_N2 / (E_N**2) if E_N > 0 else float('inf')

    print(f"\n{'='*60}")
    print(f"Second Moment Statistics (p=43)")
    print(f"{'='*60}")
    print(f"Working orbits with omega: {len(omegas)}/{n_working}")
    print(f"Total D11-orbits: {total_orbits}")
    print(f"Orbit size (2p): {orbit_size}")
    print(f"Omega values: min={min(omegas)}, max={max(omegas)}, "
          f"mean={sum(omegas)/len(omegas):.1f}, median={sorted(omegas)[len(omegas)//2]}")
    print(f"sum(omega) = {sum(omegas)}")
    print(f"sum(omega^2) = {sum(w**2 for w in omegas)}")
    print(f"Total N = sum(omega_i * {orbit_size}) = {total_N}")
    print(f"E[N] = {E_N:.6f}")
    print(f"E[N^2] = {E_N2:.6f}")
    print(f"E[N^2] / E[N]^2 = {ratio:.4f}")

    print(f"\nCross-prime comparison:")
    primes = [11, 19, 23, 31, 43]
    ratios = [2.0, 14.0, 14.5, 71.25, ratio]
    for p, r in zip(primes, ratios):
        print(f"  p={p:2d}: ratio = {r:.4f}")

    # Log-log slopes
    print(f"\nLog-log slopes (base p=11):")
    for i in range(1, len(primes)):
        if ratios[i] > 0 and ratios[0] > 0:
            alpha = math.log(ratios[i] / ratios[0]) / math.log(primes[i] / primes[0])
            print(f"  p=11 to p={primes[i]}: alpha = {alpha:.3f}")

    # Omega distribution
    print(f"\nOmega distribution (all {len(omegas)} orbits):")
    for oid, omega in sorted(omega_by_id.items()):
        print(f"  orbit {oid:5d}: omega={omega}")

    return {
        'p': 43,
        'total_orbits': total_orbits,
        'working_orbits': n_working,
        'orbits_with_omega': len(omegas),
        'orbit_size': orbit_size,
        'omega_by_orbit': omega_by_id,
        'sum_omega': sum(omegas),
        'sum_omega2': sum(w**2 for w in omegas),
        'total_N': total_N,
        'E_N': E_N,
        'E_N2': E_N2,
        'ratio': ratio,
        'comparison': {str(p): r for p, r in zip(primes, ratios)},
    }


def main():
    if not os.path.exists(BINARY):
        print(f"ERROR: Binary not found: {BINARY}")
        print("Compile with: gcc -O3 -o count_omega_fast count_omega_fast.c -lm")
        sys.exit(1)

    working = load_working_orbits()
    print(f"Total working orbits: {len(working)}")

    done = load_existing_results()
    print(f"Already have omega for: {len(done)} orbits")
    for oid in sorted(done):
        src = done[oid].get('source', 'computed')
        print(f"  orbit {oid}: omega={done[oid]['omega']} ({src})")

    remaining = [w for w in working if w['orbit_id'] not in done]
    print(f"\nRemaining to compute: {len(remaining)} orbits")

    if not remaining:
        print("All orbits computed! Computing final statistics...")
        stats = compute_stats(done, len(working))
        with open(FINAL_FILE, 'w') as f:
            json.dump(stats, f, indent=2)
        print(f"\nFinal results saved to {FINAL_FILE}")
        return

    n_workers = min(6, len(remaining))  # use 6 cores to leave some headroom
    est_per_orbit = 550  # seconds per orbit at 500 trials
    est_total = len(remaining) * est_per_orbit / n_workers
    print(f"Running with {n_workers} parallel workers, {N_TRIALS} trials each")
    print(f"Estimated ~{est_per_orbit}s per orbit, total ~{est_total/3600:.1f}h")
    print(f"Results saved incrementally to {RESULTS_FILE}\n")

    start = time.time()
    completed = 0
    n_success = 0
    n_fail = 0

    with open(RESULTS_FILE, 'a') as fout:
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {}
            for w in remaining:
                f = executor.submit(run_one_orbit, w['orbit_id'], w['d11'], N_TRIALS)
                futures[f] = w['orbit_id']

            for f in as_completed(futures):
                result = f.result()
                completed += 1
                fout.write(json.dumps(result) + '\n')
                fout.flush()

                if result.get('omega') is not None:
                    done[result['orbit_id']] = result
                    n_success += 1
                else:
                    n_fail += 1

                elapsed = time.time() - start
                rate = completed / elapsed if elapsed > 0 else 0
                eta = (len(remaining) - completed) / rate if rate > 0 else 0

                omega_str = str(result.get('omega', 'FAIL'))
                print(f"[{completed}/{len(remaining)}] orbit {result['orbit_id']:5d}: "
                      f"omega={omega_str:>4s}, "
                      f"ok={n_success} fail={n_fail} "
                      f"ETA={eta/60:.0f}m ({elapsed:.0f}s elapsed)", flush=True)

    elapsed = time.time() - start
    print(f"\nBatch done in {elapsed:.0f}s ({elapsed/3600:.1f}h)")
    print(f"Success: {n_success}, Failed: {n_fail}")

    stats = compute_stats(done, len(working))
    with open(FINAL_FILE, 'w') as f:
        json.dump(stats, f, indent=2)
    print(f"\nFinal results saved to {FINAL_FILE}")


if __name__ == '__main__':
    main()
