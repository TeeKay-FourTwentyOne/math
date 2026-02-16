# V1V2 Identity and Threshold Equality

**Date**: 2026-02-13
**Status**: PROVEN (unconditional)

---

## 1. V1V2 Identity

**Theorem.** For any symmetric D11 subset of {1,...,p-1} with |D11| = (p+1)/2 and any D12 subset of {0,...,p-1} with 0 in D12, the V1V2 common red neighbor count satisfies:

  X(d) = |D12| - 1_{d in D12}    for all d in {1,...,p-1}

where X(d) = Sigma(D11, D12, d, p) + Delta(D12, D22, d, p).

**Proof.** Define indicator functions on Z_p:
- h(x) = 1_{x in D11}  (note h(0) = 0)
- g(x) = 1_{x in D22}  (note g(0) = 0)
- f(x) = 1_{x in D12}  (note f(0) = 1)

For all x in Z_p: g(x) = 1 - h(x) - delta(x, 0), since:
- x = 0: g(0) = 0 = 1 - 0 - 1
- x in {1,...,p-1}: g(x) = 1 - h(x) (D11 and D22 partition {1,...,p-1})

The V1V2 red common neighbor count is:

  X(d) = Sum_{a in Z_p} h(a) f(d-a) + Sum_{a in Z_p} f(a) g(a-d)

For the second sum, substitute g(a-d) = 1 - h(a-d) - delta(a-d, 0):

  Sum_a f(a) g(a-d) = Sum_a f(a) - Sum_a f(a) h(a-d) - f(d)
                    = |D12| - Sum_a f(a) h(a-d) - f(d)

So:

  X(d) = Sum_a h(a) f(d-a) + |D12| - Sum_a f(a) h(a-d) - f(d)

**Key step** (uses D11 symmetry): In Sum_a f(a) h(a-d), substitute b = d - a (so a = d - b):

  Sum_a f(a) h(a-d) = Sum_b f(d-b) h(-b) = Sum_b f(d-b) h(p-b)

Since D11 is symmetric (d in D11 iff p-d in D11), we have h(p-b) = h(b). Therefore:

  Sum_a f(a) h(a-d) = Sum_b h(b) f(d-b) = Sum_a h(a) f(d-a)

The two sums cancel:

  X(d) = Sum_a h(a) f(d-a) + |D12| - Sum_a h(a) f(d-a) - f(d)
       = |D12| - f(d)
       = |D12| - 1_{d in D12}

QED.

---

## 2. V1V2 Constraints Are Automatically Satisfied

**Corollary.** For the 2-block circulant construction with symmetric D11 and |D12| = (p-1)/2:

- **V1V2 red** (d in D12): X(d) = |D12| - 1 = (p-3)/2 = n - 2, which equals the threshold n - 2.
- **V1V2 blue** (d not in D12): X(d) = |D12| = (p-1)/2 = n - 1, which equals the threshold n - 1.

Both constraints are satisfied with EQUALITY for every D12 and every symmetric D11. The V1V2 constraints are therefore **non-binding** and add **zero additional cost** to the probabilistic analysis.

---

## 3. A(d) - C(d) Identity

**Theorem.** For any symmetric D11 subset of {1,...,p-1} with |D11| = (p+1)/2:

  A(d) - C(d) = 3 - 2 * 1_{d in D11}

That is:
- d in D11: A(d) - C(d) = 1
- d in D22: A(d) - C(d) = 3

where A(d) = #{(a,b) in D11 x D11: a - b = d mod p} and C(d) = #{(a,b) in D22 x D22: a - b = d mod p}.

**Proof.** Define h(x) = 1_{x in D11} for x in Z_p, with h(0) = 0.

  A(d) = Sum_{x=1}^{p-1} h(x) h((x-d) mod p)
  C(d) = Sum_{x=1}^{p-1} (1-h(x)-delta(x,0))(1-h((x-d) mod p)-delta((x-d) mod p, 0))

For x in {1,...,p-1} and d in {1,...,p-1}:
- h(x) h((x-d) mod p) - (case analysis on whether (x-d) mod p = 0)

More directly: define h_ext(x) for all x in Z_p with h_ext(0) = 0.
Then:

  A(d) = Sum_{x in Z_p} h(x) h((x-d) mod p)

(The x = 0 term is 0 since h(0) = 0.)

Similarly:

  C(d) = Sum_{x in Z_p} j(x) j((x-d) mod p)

where j(x) = 1_{x in D22} = 1 - h(x) - delta(x,0).

  A(d) - C(d) = Sum_x [h(x) h(x-d) - j(x) j(x-d)]

**Case x != d and x-d != 0 (i.e., x != d mod p):**
  j(x) = 1 - h(x), j(x-d) = 1 - h(x-d)
  h(x)h(x-d) - (1-h(x))(1-h(x-d)) = h(x) + h(x-d) - 1

**Case x = d:** h(d)h(0) - j(d)j(0) = 0 - 0 = 0.

So:

  A(d) - C(d) = Sum_{x=1, x!=d}^{p-1} [h(x) + h((x-d) mod p) - 1]

Split the sum:

  = Sum_{x!=d} h(x) + Sum_{x!=d} h((x-d) mod p) - (p-2)

First sum: Sum_{x=1, x!=d}^{p-1} h(x) = |D11| - h(d) = n - h(d)
  (since Sum_{x=1}^{p-1} h(x) = |D11| = n, and we remove the x=d term)

Second sum: As x ranges over {1,...,p-1}\{d}, y = (x-d) mod p ranges over {1,...,p-1}
  (since x=d is removed, y=0 is removed, and the map is bijective)
  Sum_{x!=d} h((x-d) mod p) = Sum_{y=1}^{p-1} h(y) = n

Wait, but the map x -> (x-d) mod p sends {1,...,p-1}\{d} to {0,...,p-1}\{0,p-d}... no.
Let y = (x-d) mod p. When x = d: y = 0. Since we exclude x = d, y ranges over
{(1-d) mod p, ..., ((p-1)-d) mod p} \ {0} = {1,...,p-1} \ {(p-d) mod p}... wait, not quite.

Actually: x in {1,...,p-1} maps to y = (x-d) mod p in {(1-d),...,(p-1-d)} mod p = {0,...,p-1}\{(-d) mod p}.
Since d in {1,...,p-1}, (-d) mod p = p-d in {1,...,p-1}.
So y ranges over {0,...,p-1}\{p-d} = {0,1,...,p-1} minus the element p-d.
Removing x=d (which maps to y=0): y ranges over {1,...,p-1}\{p-d}.

Sum_{y in {1,...,p-1}\{p-d}} h(y) = n - h(p-d).

For symmetric D11: h(p-d) = h(d). So the second sum = n - h(d).

Therefore:

  A(d) - C(d) = (n - h(d)) + (n - h(d)) - (p-2) = 2n - 2h(d) - p + 2

With n = (p+1)/2: 2n = p+1. So:

  A(d) - C(d) = p + 1 - 2h(d) - p + 2 = 3 - 2h(d)

- d in D11: h(d) = 1, A-C = 1
- d in D22: h(d) = 0, A-C = 3

QED.

---

## 4. V1V1/V2V2 Threshold Equality

**Theorem.** For symmetric D11, the V1V1 and V2V2 binding thresholds are EQUAL at every distance d.

Specifically:
- **d in D11:** T_{V1V1 red} = n-2 - A(d) and T_{V2V2 blue} = n-3 - C(d). These are equal since:
  T_{V1V1 red} - T_{V2V2 blue} = (n-2-A) - (n-3-C) = 1 - (A-C) = 1 - 1 = 0.

- **d in D22:** T_{V1V1 blue} = n+1 - A(d) and T_{V2V2 red} = n-2 - C(d). These are equal since:
  T_{V1V1 blue} - T_{V2V2 red} = (n+1-A) - (n-2-C) = 3 - (A-C) = 3 - 3 = 0.

**Corollary.** The constraint model has exactly R = (p-1)/2 independent binding constraints (one per complementary pair {d, p-d}), all of the form:

  B(d) <= T(d)

where T(d) = (p-3)/2 - A(d) for d in D11, and T(d) = (p-3)/2 - C(d) for d in D22.

The V1V2 constraints contribute nothing. The V1V1 blue (d in D22) and V2V2 red (d in D22) constraints give the SAME threshold. The V1V1 red (d in D11) and V2V2 blue (d in D11) constraints give the SAME threshold.

---

## 5. Implications for c_0 Computation

The c_0 computation in the first moment analysis is **already correct** using only the R = (p-1)/2 binding constraints B(d) <= T(d):

  c_0(D11) = Pr[all B(d_i) <= T(d_i)] / prod Pr[B(d_i) <= T(d_i)]

No additional constraints from V1V2 need to be included. The values computed at p = 11, 19, 23 are:

| p | c_0 (working D11) | log_2(c_0) |
|---|---|---|
| 11 | 0.700 | -0.51 |
| 19 | 0.033 | -4.94 |
| 23 | 0.006-0.036 | -4.8 to -7.4 |

These are correct and represent the full correlation penalty.

---

## 6. Verified Numerically

All results verified at p = 11, 19, 23:
- V1V2 identity: 600+ (D11, D12) pairs tested, all pass
- A(d)-C(d) identity: all D11 orbits at p = 11, 19, 23 tested
- Threshold equality: all D11 orbits tested, T_{V1V1} = T_{V2V2} in every case
- V1V2 never violated: exhaustive check at p = 11 (all 210 D12), p = 19 (all 43758 D12)
