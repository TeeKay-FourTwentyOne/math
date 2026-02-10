"""
Modified Paley construction for m ≡ 1 (mod 4).

The pure Paley graph on GF(q) with q ≡ 1 (mod 4) is a strongly regular graph
with parameters (q, (q-1)/2, (q-5)/4, (q-1)/4).

For q = 49: SRG(49, 24, 11, 12). This means each red edge has 11 common
red neighbors and each blue edge has 12 common blue neighbors.

For the Ramsey book graph R(B_24, B_25):
- Red edges need < 24 common red neighbors (threshold 23)
- Blue edges need < 25 common blue neighbors (threshold 24)

The Paley graph has 11 common red neighbors per red edge (well under 23),
BUT that's only counting within one block. With the 2-block structure,
the cross-block contribution adds more neighbors.

The issue with my attempt: using D11 = D22 = QR and D12 = QR ∪ {0} gives
d1 = d2 = 49 = m, meaning every V1 vertex is adjacent to every V2 vertex.
That makes the common neighbor count huge.

The correct construction for m ≡ 1 (mod 4) must use different sets.
Let me look at this more carefully.

For the book graph B_n = K_2 + K̄_n, we need:
- N = 4n-2 = 2m vertices
- No red B_{n-1}: no red edge with n-1 common red neighbors
- No blue B_n: no blue edge with n common blue neighbors

For m ≡ 1 (mod 4), the Paley construction typically gives:
D11 = QR, |D11| = (m-1)/2
D22 = QR, |D22| = (m-1)/2
D12 depends on the specific construction.

But we showed that for m ≡ 3 (mod 4), the universal pattern is:
D22 = complement(D11), D11 symmetric, |D12| = n-1.

For m ≡ 1 (mod 4), different patterns may apply.
Let me check if the complement approach works here.

For n=25, m=49:
- If D22 = complement(D11): |D22| = 48 - |D11|
- For V1V2 theorem to apply (V1V2 free): need |D12| = n-1 = 24
- Standard sizes: |D11| = n = 25, |D12| = 24, |D22| = 23
- But |D11| = 25 is odd and D11 symmetric in Z_49...
  In Z_49 (odd m), symmetric D11 must have even size. So |D11| = 25 is impossible.
- Try |D11| = 24, |D12| = 24, |D22| = 24.
  Then d1 = 48, d2 = 48. And V1V2 theorem gives common neighbors = 23 or 24.
  Red threshold = 23, blue threshold = 24. ✓

Wait, but Z_49 is fine for D11 symmetric with |D11|=24 since 24 is even (12 pairs).

Actually wait -- Z_49 also has the element 0, and D11 ⊂ {1,...,48}.
48 elements, and negation gives pairs {x, 49-x} for x=1,...,24. That's 24 pairs.
For |D11| = 24: need exactly 12 of 24 pairs. OK, even and fine.

For |D11| = 24, |D12| = 24 (= n-1), D22 = complement(D11) giving |D22| = 24:
- V1V2 theorem: common neighbors = |D12| - [d in D12] = 23 or 24. ✓
- d1 = |D11| + |D12| = 48
- d2 = |D22| + |D12| = 48
- Symmetric! Both blocks have degree 48.

Now the V1V1 constraint:
- Average A(d)+B(d) = (24*23 + 24*23)/48 = 2*24*23/48 = 23.0
- Red threshold: n-2 = 23. Average = threshold!
- So we need A(d)+B(d) <= 23 for ALL d in D11.
- Since average = 23, this means ALL red values must be <= average.
  Any value above 23 violates, so the fluctuation must be <= 0 for red.

This is MUCH tighter than the m ≡ 3 mod 4 case!

For the Paley case, the strongly regular property gives:
Delta(QR, QR, d) = (q-5)/4 for d in QR, and (q-1)/4 for d in QNR.
For q=49: Delta(QR,QR,d) = 11 for d in QR, 12 for d in QNR.

If D11 = QR (using GF(49) Cayley graph, not Z_49 circulant):
A(d) = 11 for d in D11, A(d) = 12 for d not in D11.

Then V1V1(d) = 11 + B(d) for d in D11. Need 11 + B(d) <= 23, so B(d) <= 12.
avg B(d) = 24*23/48 = 11.5. Need B(d) <= 12 for d in D11 -- that's just 0.5 above average.

And V1V1 blue: 12 + B(d) for d not in D11. Need blue_common <= 24.
blue_common = N-2-2*d1 + (12 + B(d)) = 96 - 96 + 12 + B(d) = 12 + B(d).
So need 12 + B(d) <= 24, i.e., B(d) <= 12. Same constraint!

So for Paley on GF(q), the constraint is simply: B(d) <= 12 for all d != 0.
Since avg B(d) = 11.5, we need B(d) to never exceed average by more than 0.5.
This means all B(d) values must be either 11 or 12!

If D12 is also chosen as a "Paley-like" set:
Delta(D12, D12, d) would be either 11 or 12 for all d != 0.
THIS WORKS if D12 is a set with constant-ish autocorrelation.

In fact, if D12 is also a set of 24 elements such that Delta(D12,D12,d) ∈ {11,12}
for all d != 0, then ALL constraints are satisfied.

A perfect (49, 24, 11)-difference set would give Delta = 11 for all d.
Does such a set exist? The Paley difference set: QR in GF(49) IS such a set
when viewed over the correct additive group.

Let me verify: in GF(49), the QR set has |QR| = 24, and for the additive
group of GF(49), Delta(QR, QR, d) should be constant for d != 0 if QR
is a difference set.
"""

import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

m = 49
n = 25

class GF49:
    def __init__(self, a, b):
        self.a = a % 7
        self.b = b % 7

    def __mul__(self, other):
        return GF49(self.a * other.a - self.b * other.b,
                     self.a * other.b + self.b * other.a)

    def __eq__(self, other):
        return self.a == other.a and self.b == other.b

    def __hash__(self):
        return hash((self.a, self.b))

    def __repr__(self):
        return f"({self.a}+{self.b}i)"

    def is_zero(self):
        return self.a == 0 and self.b == 0

    def to_int(self):
        return self.a + 7 * self.b

    @staticmethod
    def from_int(n):
        return GF49(n % 7, n // 7)

def order_gf49(elem):
    if elem.is_zero():
        return 0
    power = elem
    for k in range(1, 49):
        if power.a == 1 and power.b == 0:
            return k
        power = power * elem
    return -1

# Find generator
gen = None
for a in range(7):
    for b in range(7):
        elem = GF49(a, b)
        if not elem.is_zero() and order_gf49(elem) == 48:
            gen = elem
            break
    if gen:
        break

print(f"GF(49) generator: {gen}")

# Compute QR
powers = []
power = GF49(1, 0)
for k in range(48):
    powers.append(power)
    power = power * gen

qr_set = {powers[2*k] for k in range(24)}
qr_ints = {e.to_int() for e in qr_set}

print(f"QR as ints (GF49 encoding): {sorted(qr_ints)}")
print(f"|QR| = {len(qr_ints)}")

# GF(49) subtraction
def gf_sub_int(a, b):
    ea = GF49.from_int(a)
    eb = GF49.from_int(b)
    return GF49((ea.a - eb.a) % 7, (ea.b - eb.b) % 7).to_int()

# Compute Delta(QR, QR, d) over GF(49)+ for each nonzero d
print("\nDelta(QR, QR, d) over GF(49)+:")
deltas = {}
for d in range(1, 49):
    count = 0
    for a in qr_ints:
        diff = gf_sub_int(a, d)
        if diff in qr_ints:
            count += 1
    deltas[d] = count

vals = sorted(set(deltas.values()))
print(f"Unique values: {vals}")
for v in vals:
    count = sum(1 for d in range(1, 49) if deltas[d] == v)
    in_qr = sum(1 for d in range(1, 49) if deltas[d] == v and d in qr_ints)
    print(f"  Delta = {v}: {count} times ({in_qr} in QR, {count - in_qr} in QNR)")

# Perfect! QR is a (49, 24, 11)-difference set in GF(49)+.
# This means Delta(QR, QR, d) = 11 for ALL d != 0.

# Now: the Ramsey construction.
# We work over GF(49)+ (i.e., Z_7 x Z_7 as additive group).
# D11 = D22 = QR, D12 = QR ∪ {0}.
# Common neighbors:
# V1V1 at d: Delta(D11,D11,d) + Delta(D12,D12,d) = 11 + Delta(QR∪{0}, QR∪{0}, d)
# Delta(QR∪{0}, QR∪{0}, d) = |{a in QR∪{0}: (a-d) in QR∪{0}}|
#   = |{a in QR: (a-d) in QR}| + |{a in QR: (a-d) = 0}| + |{0: (-d) in QR}| + |{0: (-d) = 0}|
#   = 11 + [d in QR] + [d in QR] + 0 (since d != 0, -d != 0... wait, -d could be in QR)
# Hmm, let me be more careful.
# QR∪{0} has 25 elements. Delta(QR∪{0}, QR∪{0}, d) for d != 0:
# = |{a in QR∪{0}: (a-d) in QR∪{0}}|
# a ranges over QR and {0}:
# - a in QR: (a-d) in QR∪{0} iff (a-d) in QR or a-d=0 (i.e., a=d)
#   = Delta(QR,QR,d) + [d in QR] = 11 + [d in QR]
# - a = 0: (0-d) = -d in QR∪{0} iff -d in QR or -d=0
#   Since d != 0, -d != 0. And -d in QR iff d in QR (since -1 is QR in GF(49)).
#   So [a=0 contributes 1] iff d in QR.
# Total: 11 + [d in QR] + [d in QR] = 11 + 2*[d in QR]

# So V1V1(d) = 11 + (11 + 2*[d in QR]) = 22 + 2*[d in QR]
# For d in QR (= D11): V1V1 = 24. Red threshold = 23. VIOLATION!
# For d in QNR: V1V1 = 22. Red threshold = 23. OK (blue common = N-2-2*48+22 = -2... hmm)

print("\nV1V1 with D11=D22=QR, D12=QR∪{0}:")
print(f"V1V1(d) = 22 + 2*[d in QR]")
print(f"For d in QR: V1V1 = 24 > {n-2} = red threshold -> VIOLATION")
print(f"This construction FAILS for book graphs.")

# So the simple D11 = D22 = QR, D12 = QR ∪ {0} doesn't work.
# We need a different D12 or different setup.

# Let's try the COMPLEMENT approach: D22 = complement(D11)
# With D11 = QR (24 elements), D22 = QNR (24 elements), |D12| = 24.
# V1V2 is free (theorem applies since D22 = complement).
# V1V1(d) = Delta(QR,QR,d) + Delta(D12,D12,d) = 11 + B(d)
# Need 11 + B(d) <= 23 for d in QR, so B(d) <= 12.
# avg B = 24*23/48 = 11.5. Need B(d) <= 12 = avg + 0.5.

# V2V2: by our theorem, V2V2(d) = A(d) + B(48-d in GF49) + (m-2-2|D11|) + 2*[d in D11]
# Wait, the V2V2 theorem was for Z_m circulants. Here we're on GF(49)+.
# But the algebra is the same since the theorem only uses group properties.
# V2V2(d) = Delta(QNR,QNR,d) + Delta(D12T,D12T,d)
# Delta(QNR,QNR,d) for d != 0: since QR is a (49,24,11)-diff set,
# QNR = complement of QR in GF(49)*.
# Delta(QNR,QNR,d) = Delta({1..48}\QR, {1..48}\QR, d)
# By the complement formula: = A(d) + (m-2-2|D11|) + 2[d in D11]
# = 11 + (49-2-48) + 2[d in QR] = 11 - 1 + 2[d in QR] = 10 + 2[d in QR]

# For d in QNR (red in V2V2): Delta(QNR,QNR,d) = 10 (since d not in QR)
# For d in QR (blue in V2V2): Delta(QNR,QNR,d) = 12

# D12T = {-x : x in D12} over GF(49).
# Delta(D12T,D12T,d) = Delta(D12,D12,-d) = B(-d in GF49) = B(neg_d)

# V2V2 red (d in QNR): 10 + B(neg_d) <= 23, so B(neg_d) <= 13.
# V2V2 blue (d in QR): blue_common = (N-2) - 2*d2 + 12 + B(neg_d)
#   = 96 - 96 + 12 + B(neg_d) = 12 + B(neg_d) <= 24, so B(neg_d) <= 12.

# So: V1V1 red needs B(d) <= 12 for d in QR
#     V2V2 blue needs B(-d) <= 12 for d in QR, i.e., B(d) <= 12 for d in QR (since -QR = QR)
#     V2V2 red needs B(-d) <= 13 for d in QNR, i.e., B(d) <= 13 for d in QNR (since -QNR = QNR)
#     V1V1 blue needs B(d) s.t. blue_common <= 24:
#       blue_common = 96 - 96 + 12 + B(d) = 12 + B(d) for d in QNR
#       Need B(d) <= 12 for d in QNR.

# Summary: need B(d) <= 12 for ALL d != 0.
# B(d) = Delta(D12,D12,d) with |D12| = 24, avg = 23*24/48 = 11.5.
# Need max B(d) <= 12. This is exactly the condition that D12 is an
# "almost-difference set" or has nearly flat autocorrelation.

# QR itself has Delta(QR,QR,d) = 11 for all d != 0 (PERFECT!).
# So D12 = QR would give B(d) = 11 <= 12 for all d. But |D12| = |QR| = 24 ✓ and 0 not in QR.
# We need 0 in D12!

# D12 = QR ∪ {0} \ {some element}: |D12| = 24, 0 in D12.
# Or D12 = QR shifted: {x + c : x in QR} for some c?
# Or some other set of size 24 containing 0 with Delta <= 12.

# Actually, we could use D12 = QR ∪ {0} and then |D12| = 25.
# But we need |D12| = n-1 = 24. Conflict!

# Hmm, wait. Let me recheck. For general n:
# |D12| = n - 1 was the pattern for m ≡ 3 mod 4 cases.
# For m ≡ 1 mod 4, the pattern might be different.
# Let me reconsider.

# Actually, the V1V2 theorem said: if D22 = comp(D11), D11 symmetric, |D12| = n-1,
# then V1V2 is free. For n=25: |D12| = 24. And QR has 24 elements. But 0 not in QR.

# Can we take D12 = QR? Then 0 not in D12. Is that OK?
# The "0 in D12" pattern was universal for m ≡ 3 mod 4 cases.
# For m ≡ 1 mod 4, maybe 0 is NOT in D12.

# Let's test: D11 = QR, D22 = QNR, D12 = QR (without 0), |D12| = 24.
# V1V2 theorem still applies (requires |D12| = n-1 = 24 ✓ and D22 = comp(D11) ✓, D11 sym ✓).
# BUT: the theorem gives V1V2(d) = |D12| - [d in D12].
# For d in D12 = QR: V1V2 = 23 = n-2 ✓
# For d not in D12 = QNR ∪ {0}: V1V2 = 24 = n-1 ✓
# STILL WORKS even without 0 in D12!

# V1V1(d) = 11 + B(d) where B(d) = Delta(QR, QR, d) = 11 for all d.
# So V1V1(d) = 22 for all d. Red (d in QR): 22 <= 23 ✓.
# Blue (d in QNR): blue_common = 96 - 2*48 + 22 = 22. 22 <= 24 ✓.

# V2V2: D12T = {-x : x in QR} = QR (since -1 is QR in GF(49), m ≡ 1 mod 4).
# So D12T = QR. Delta(D12T, D12T, d) = Delta(QR, QR, d) = 11.
# Delta(QNR, QNR, d) = 10 + 2[d in QR] (from complement formula).
# V2V2(d) = (10 + 2[d in QR]) + 11 = 21 + 2[d in QR].
# For d in QNR (= D22, red): V2V2 = 21 <= 23 ✓
# For d in QR (blue in V2V2): blue_common = 96 - 96 + 23 = 23 <= 24 ✓

print("\n" + "=" * 80)
print("PALEY CONSTRUCTION: D11 = QR, D22 = QNR, D12 = QR (over GF(49)+)")
print("=" * 80)

D11 = qr_ints
D22 = {e.to_int() for e in (set(GF49(a,b) for a in range(7) for b in range(7)) - {GF49(0,0)} - qr_set)}
D12 = qr_ints  # D12 = QR (no 0)

print(f"D11 (QR): {sorted(D11)}")
print(f"D22 (QNR): {sorted(D22)}")
print(f"D12 (QR): {sorted(D12)}")
print(f"|D11| = {len(D11)}, |D12| = {len(D12)}, |D22| = {len(D22)}")
print(f"0 in D12: {0 in D12}")

d1 = len(D11) + len(D12)
d2 = len(D22) + len(D12)
print(f"d1 = {d1}, d2 = {d2}")

# Full verification using GF(49) arithmetic
max_red = 0
max_blue = 0
violations = []
red_thresh = n - 2  # 23
blue_thresh = n - 1  # 24

# V1V1
for d in range(1, 49):
    A_d = 0
    for a in D11:
        if gf_sub_int(a, d) in D11:
            A_d += 1
    B_d = 0
    for a in D12:
        if gf_sub_int(a, d) in D12:
            B_d += 1
    common = A_d + B_d

    if d in D11:
        max_red = max(max_red, common)
        if common > red_thresh:
            violations.append(("V1V1_red", d, common))
    else:
        blue_common = (98 - 2) - d1 - d1 + common
        max_blue = max(max_blue, blue_common)
        if blue_common > blue_thresh:
            violations.append(("V1V1_blue", d, blue_common))

# V2V2
D12T_gf = {gf_sub_int(0, x) for x in D12}
for d in range(1, 49):
    C_d = 0
    for a in D22:
        if gf_sub_int(a, d) in D22:
            C_d += 1
    E_d = 0
    for a in D12T_gf:
        if gf_sub_int(a, d) in D12T_gf:
            E_d += 1
    common = C_d + E_d

    if d in D22:
        max_red = max(max_red, common)
        if common > red_thresh:
            violations.append(("V2V2_red", d, common))
    else:
        blue_common = (98 - 2) - d2 - d2 + common
        max_blue = max(max_blue, blue_common)
        if blue_common > blue_thresh:
            violations.append(("V2V2_blue", d, blue_common))

# V1V2
for d in range(49):
    # Sigma(D11, D12, d) using GF(49)
    sig = 0
    for a in D11:
        if gf_sub_int(d, a) in D12:
            sig += 1
    # Delta(D12, D22, d)
    delt = 0
    for a in D12:
        if gf_sub_int(a, d) in D22:
            delt += 1
    common = sig + delt

    if d in D12:
        max_red = max(max_red, common)
        if common > red_thresh:
            violations.append(("V1V2_red", d, common))
    else:
        blue_common = (98 - 2) - d1 - d2 + common
        max_blue = max(max_blue, blue_common)
        if blue_common > blue_thresh:
            violations.append(("V1V2_blue", d, blue_common))

print(f"\nMax red common: {max_red} (threshold {red_thresh})")
print(f"Max blue common: {max_blue} (threshold {blue_thresh})")
print(f"Violations: {len(violations)}")

if violations:
    print("INVALID! First violations:")
    for v in violations[:10]:
        print(f"  {v}")
else:
    print("VALID! The Paley construction works for n=25!")

    result = {
        "n": n,
        "target_vertices": 98,
        "solver": "Paley-GF49",
        "field": "GF(49) = GF(7)[i], i^2=-1",
        "encoding": "a+bi -> a+7*b",
        "construction": "D11 = D12 = QR(GF(49)*), D22 = QNR(GF(49)*)",
        "note": "Uses GF(49) addition (Z_7 x Z_7), NOT Z_49 arithmetic. 0 not in D12.",
        "parameters": {
            "m": m,
            "D11": sorted(D11),
            "D12": sorted(D12),
            "D22": sorted(D22)
        },
        "verification": {
            "max_red_common": max_red,
            "max_blue_common": max_blue,
            "red_threshold": red_thresh,
            "blue_threshold": blue_thresh,
            "valid": True
        }
    }

    outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "solution_n25_paley.json")
    with open(outpath, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"Saved to {outpath}")
