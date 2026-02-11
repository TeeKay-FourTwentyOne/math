"""
Feasibility analysis for n=32, m=63.

Key question: Is it even possible to have A(d)+B(d) <= 30 for ALL d in D11,
given that sum(A(d)+B(d)) = 1922 and |D11| constraints?

We use the Welch bound and combinatorial arguments.
"""

import numpy as np
from itertools import combinations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from n32_fast_numpy import compute_delta_fft

m = 63
n = 32

print("=== Feasibility Analysis for n=32 ===\n")

# For |D11|=32:
# sum(A(d), d=1..62) = |D11|^2 - |D11| = 32*31 = 992
# sum(B(d), d=1..62) = |D12|^2 - |D12| = 31*30 = 930
# Total: 1922
# Need A(d)+B(d) <= 30 for all d in D11 (32 differences)
# Need A(d)+B(m-d) <= 30 for all d in D11 (32 differences)
# These 32+32=64 constraints (on pairs (A(d)+B(d)) and (A(d)+B(m-d)))
# must be satisfied simultaneously.

# If d in D11 and (m-d) in D11 (which happens since D11 is symmetric):
# Then both A(d)+B(d) <= 30 and A(m-d)+B(m-d) <= 30 AND A(d)+B(m-d) <= 30 AND A(m-d)+B(d) <= 30
# Since A is symmetric (A(d) = A(m-d)), these reduce to:
# A(d) + B(d) <= 30 and A(d) + B(m-d) <= 30

# Sum of A(d)+B(d) over all d in D11: <= 32*30 = 960
# Sum over d not in D11: 1922 - (sum over D11) >= 1922 - 960 = 962
# For d not in D11 (30 differences): need <= 33 each
# Total capacity: 30*33 = 990 >= 962 ✓

# But we also need A(d)+B(m-d) <= 30 for d in D11
# sum(A(d)+B(m-d), d in D11) = sum(A(d), d in D11) + sum(B(m-d), d in D11)
# Since D11 is symmetric, {m-d : d in D11} = D11, so sum(B(m-d), d in D11) = sum(B(d), d in D11)
# So this is the same as sum(A(d)+B(d), d in D11) <= 960
# Same constraint!

print("For |D11|=32, |D12|=31:")
print("  Total sum = 1922")
print("  Budget for D11 (32 diffs, thresh 30): 960")
print("  Budget for D22 (30 diffs, thresh 33): 990")
print("  Total budget: 960 + 990 = 1950 >= 1922 ✓")
print("  Slack: 1950 - 1922 = 28")
print()

# But the constraint is much tighter than just the sum.
# We need EACH d in D11 to satisfy A(d)+B(d) <= 30,
# and ALSO A(d)+B(m-d) <= 30.
# These are TWO constraints per symmetric pair.

# Let's think in terms of symmetric pairs.
# D11 has 16 symmetric pairs {d, m-d}.
# For each pair (d, m-d) where both are in D11:
#   A(d) + B(d) <= 30   [since d in D11]
#   A(d) + B(m-d) <= 30 [since d in D11, V2V2 constraint]
# But A(d) = A(m-d), so:
#   A(d) + B(d) <= 30
#   A(d) + B(m-d) <= 30
# And similarly from the (m-d) side:
#   A(m-d) + B(m-d) <= 30  =>  A(d) + B(m-d) <= 30  (same)
#   A(m-d) + B(d) <= 30    =>  A(d) + B(d) <= 30     (same)
# So effectively: A(d) + max(B(d), B(m-d)) <= 30 for each pair.

# What about B? Is B(d) = B(m-d)?
# B(d) = Delta(D12, D12, d). D12 is NOT necessarily symmetric.
# So B(d) != B(m-d) in general.
# However, B(d) = number of pairs (x,y) in D12 with x-y=d mod m.
# B(m-d) = number of pairs with x-y=m-d = -d mod m.
# So B(m-d) = |{(x,y): y-x = d}| = B(d) if we're counting unordered,
# but for ordered pairs: B(m-d) = |{(x,y): x-y=-d}| = |{(y,x): y-x=d}| = B(d)!
# Wait, Delta(D12,D12,d) = |{a in D12: (a-d) in D12}| = |{(a,b) in D12xD12: a-b=d}|
# Delta(D12,D12,m-d) = |{a in D12: (a-m+d) in D12}| = |{(a,b): a-b=-d}| = |{(b,a): b-a=d}| = Delta(D12,D12,d)
# So B(d) = B(m-d)!!

# Therefore A(d)+B(d) = A(d)+B(m-d) always! The V1V1 and V2V2 constraints are actually identical!

print("Key observation: B(d) = B(m-d) always (autocorrelation is symmetric).")
print("Therefore A(d)+B(d) = A(d)+B(m-d) for all d.")
print("The V1V1 and V2V2 constraints are redundant!")
print()

# Wait, but the V2V2 uses Delta_12T, not Delta_12.
# V2V2 common = Delta(D11,D11,d) + Delta(D12T,D12T,d) + v22_const + 2[d in D11]
# where D12T = {-x mod m : x in D12}
# Delta(D12T, D12T, d) = |{a in D12T: (a-d) in D12T}|
#   = |{a: -a in D12, -(a-d) in D12}|
#   = |{a: -a in D12, d-a in D12}|
#   = |{b: b in D12, d-(-b) in D12}| (let b = -a)
#   = |{b: b in D12, (d+b) in D12}|
#   = |{b: b in D12, (b-(-d)) in D12}|
#   = Delta(D12, D12, -d)
#   = Delta(D12, D12, m-d)

print("Delta(D12T, D12T, d) = Delta(D12, D12, m-d) = B(m-d)")
print("But we showed B(m-d) = B(d).")
print("So Delta(D12T, D12T, d) = B(d) always.")
print()
print("This means the V2V2 constraint IS identical to V1V1!")
print("The 8 violations always come as 4 pairs precisely because V1V1(d) mirrors V2V2(d).")
print()

# Now let's verify this claim numerically
print("=== Numerical Verification ===")
import random
random.seed(42)
D11_ind = np.zeros(m, dtype=np.int64)
# Pick 16 random pairs
pairs = []
for x in range(1, m):
    neg_x = (-x) % m
    if x <= neg_x:
        pairs.append((x, neg_x))
selected = random.sample(range(len(pairs)), 16)
for i in selected:
    D11_ind[pairs[i][0]] = 1
    D11_ind[pairs[i][1]] = 1

D12_ind = np.zeros(m, dtype=np.int64)
D12_ind[0] = 1
for x in random.sample(range(1, m), 30):
    D12_ind[x] = 1

delta_11 = compute_delta_fft(D11_ind, m)
delta_12 = compute_delta_fft(D12_ind, m)
D12T_ind = np.zeros(m, dtype=np.int64)
for x in range(m):
    if D12_ind[x]:
        D12T_ind[(-x) % m] = 1
delta_12T = compute_delta_fft(D12T_ind, m)

# Verify B(d) = B(m-d)
all_equal = True
for d in range(1, m):
    if delta_12[d] != delta_12[(m - d) % m]:
        all_equal = False
        print(f"  B({d}) = {delta_12[d]} != B({m-d}) = {delta_12[(m-d)%m]}")
print(f"B(d) = B(m-d) for all d: {all_equal}")

# Verify delta_12T[d] = delta_12[m-d]
all_equal2 = True
for d in range(1, m):
    if delta_12T[d] != delta_12[(m - d) % m]:
        all_equal2 = False
print(f"Delta(D12T,D12T,d) = B(m-d) for all d: {all_equal2}")

print()
print("=== Summary ===")
print("The V1V1 and V2V2 constraints are mathematically identical.")
print("The 4 violated differences produce 8 violations (4 V1V1 + 4 V2V2).")
print("So the real question is: can we find D11, D12 where A(d)+B(d) <= 30")
print("for all 32 differences d in D11, given total sum = 1922?")
print()

# Check: is the constraint satisfiable with sum 1922 and 32 diffs at <=30 and 30 diffs at <=33?
# Total capacity = 32*30 + 30*33 = 960 + 990 = 1950
# Need 1922. Excess capacity = 28.
# So we can have at most 28 units distributed across the 30 "slack" positions.
# Each slack position has 3 units of slack (33-30=3). So 30*3 = 90 available slack.
# We need the 32 "tight" positions to average exactly 30 (total 960 vs budget 960).
# And the 30 "loose" positions to total exactly 962 (1922 - 960).
# Average loose = 962/30 = 32.07, with max 33. This is fine.

# BUT: we also need A(d) to satisfy its own constraints.
# A(d) is an autocorrelation of a set of size 32.
# sum(A(d), d=1..62) = 32*31 = 992
# For d in D11: A(d) + B(d) <= 30, so A(d) <= 30 - B(d)
# For d not in D11: A(d) + B(d) <= 33, so no strong constraint on A(d)
# B(d) >= 0, so A(d) <= 30 for d in D11

# Average A(d) = 992/62 = 16.0 for D11 diffs.
# We need A(d) <= 30 for d in D11 -- easy, since average is 16!
# The real constraint is A(d) + B(d) jointly.

# So the joint constraint IS satisfiable from a counting perspective.
# The issue is whether we can find specific sets D11, D12 that achieve it.

print("Counting argument says the constraints are satisfiable (capacity 1950 > need 1922).")
print("But this doesn't guarantee an actual construction exists.")
print("The autocorrelation values are highly correlated (not independent).")
print()

# Try to compute the actual distribution of A(d)+B(d) for many random D11/D12
print("=== Random Sampling: Distribution of max(A(d)+B(d)) over d in D11 ===")
min_max_sum = 999
for trial in range(1000):
    random.seed(trial * 137)
    D11i = np.zeros(m, dtype=np.int64)
    sel = random.sample(range(len(pairs)), 16)
    for i in sel:
        D11i[pairs[i][0]] = 1
        D11i[pairs[i][1]] = 1

    D12i = np.zeros(m, dtype=np.int64)
    D12i[0] = 1
    for x in random.sample(range(1, m), 30):
        D12i[x] = 1

    dA = compute_delta_fft(D11i, m)
    dB = compute_delta_fft(D12i, m)

    max_in_D11 = max(int(dA[d] + dB[d]) for d in range(1, m) if D11i[d])
    if max_in_D11 < min_max_sum:
        min_max_sum = max_in_D11

if min_max_sum <= 30:
    print(f"Found random pair with max(A+B) for d in D11 = {min_max_sum} <= 30!")
else:
    print(f"Best random pair: max(A+B) for d in D11 = {min_max_sum} (need <= 30)")
    print("1000 random trials couldn't achieve threshold. This IS hard.")
