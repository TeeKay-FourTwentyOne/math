# Proof of Equi-Covariance: Cov[B(d1), B(d2)] is Constant for Non-Complementary Pairs

**Date**: 2026-02-10
**Status**: Complete
**Purpose**: Fill the gap in Fact 3 of proof_L6_multivariate_clt.md by direct indicator computation.

---

## 1. Statement

**Theorem.** Let p be an odd prime, k = (p-3)/2, N = p-1, and S a uniformly random k-subset of {1,...,p-1}. Define Y_a = 1[a in S] and

  B(d) = Y_{p-d} + Y_d + Q(d),    where Q(d) = sum_{a in T(d)} Y_a * Y_{(a-d) mod p}

with T(d) = {1,...,p-1} \ {d}. Then for all pairs (d1, d2) with d1, d2 in {1,...,p-1}, d1 != d2, and d1 + d2 != 0 (mod p):

  Cov[B(d1), B(d2)] = -(p+1) / (8(p-2)).

In particular, the covariance is the **same** for all non-complementary pairs, regardless of the specific values of d1 and d2.

---

## 2. Background: The Gap in the Symmetry Argument

The earlier proof (Fact 3 in proof_L6_multivariate_clt.md) argued that (Z/pZ)* acts on {1,...,p-1} by multiplication, preserving the uniform k-subset measure, and maps B(d) to B(gd). This shows that the joint distribution of (B(d1), B(d2)) depends only on the ratio d2/d1 mod p. The further involution d -> p-d (using B(d) = B(p-d)) reduces the dependence to the unordered pair {d2/d1, -d2/d1} mod p.

However, this only shows that Cov[B(d1), B(d2)] depends on the "orbit" of d2/d1 under g -> -g. For distinct ratios r != r' (with r, r' not in {1, -1}), the symmetry argument alone does not force Cov to be the same. The argument that "all other ratios give the same joint distribution" is the gap.

We fill this gap by direct computation, showing that the expectation E[B(d1)*B(d2)] has the **same value** for all non-complementary pairs without appealing to any group action beyond the structure already implicit in the indicator decomposition.

---

## 3. Setup and Notation

**Random model**: S is a uniformly random k-subset of [N] = {1,...,p-1} where k = (p-3)/2 and N = p-1. For a in [N], define Y_a = 1[a in S].

**Hypergeometric moments**: The expectations of products of distinct indicators are:
- q_j = E[Y_{a_1} ... Y_{a_j}] = k^{(j)} / N^{(j)}

where x^{(j)} = x(x-1)...(x-j+1) is the falling factorial, and the a_i are distinct elements of [N]. By the indicator identity Y_a^2 = Y_a, the expectation of any product of Y's depends only on the **number of distinct indices**.

**Decomposition of B(d)**:
- L_1(d) = Y_{p-d}  (the "alpha" linear term)
- L_2(d) = Y_d  (the "beta" linear term)
- Q(d) = sum_{a in T(d)} Y_a * Y_{(a-d) mod p}  (the quadratic sum, T(d) = {1,...,p-1} \ {d})

So B(d) = L_1(d) + L_2(d) + Q(d).

**Key notation**: For d1 != 0, d2 != 0 with d1 != d2:
- alpha_1 = p - d1,  beta_1 = d1
- alpha_2 = p - d2,  beta_2 = d2
- T_1 = {1,...,p-1} \ {d1},  T_2 = {1,...,p-1} \ {d2}

Note that alpha_1, beta_1, alpha_2, beta_2 are all in {1,...,p-1}, and the four values {alpha_1, beta_1, alpha_2, beta_2} = {p-d1, d1, p-d2, d2}.

---

## 4. Expansion of E[B(d1) * B(d2)]

Expanding the product:

  B(d1) * B(d2) = [L_1(d1) + L_2(d1) + Q(d1)] * [L_1(d2) + L_2(d2) + Q(d2)]

yields nine cross-terms. We analyze each and show the total depends only on p and on whether d1 + d2 = 0 (mod p).

### 4.1 Terms (1)-(4): Y * Y products

These are:
- (1) Y_{alpha_1} * Y_{alpha_2},  E = q_{|{alpha_1, alpha_2}|}
- (2) Y_{alpha_1} * Y_{beta_2},   E = q_{|{alpha_1, beta_2}|}
- (3) Y_{beta_1} * Y_{alpha_2},   E = q_{|{beta_1, alpha_2}|}
- (4) Y_{beta_1} * Y_{beta_2},    E = q_{|{beta_1, beta_2}|}

Each evaluates to q_1 if the two indices coincide, or q_2 if they are distinct.

**Analysis of coincidences**:
- (1): alpha_1 = alpha_2 iff p-d1 = p-d2 iff d1 = d2 (excluded).
- (2): alpha_1 = beta_2 iff p-d1 = d2 iff d1 + d2 = p = 0 mod p.
- (3): beta_1 = alpha_2 iff d1 = p-d2 iff d1 + d2 = 0 mod p.
- (4): beta_1 = beta_2 iff d1 = d2 (excluded).

**Conclusion**: For non-complementary pairs (d1 + d2 != 0 mod p), all four pairs of indices are distinct, so:

  E_1 + E_2 + E_3 + E_4 = 4 * q_2.

This is **independent of d1, d2**.

(For complementary pairs d1 + d2 = 0: terms (2) and (3) contribute q_1 instead of q_2, giving 2*q_1 + 2*q_2.)

### 4.2 Terms (5)-(8): Y * Q cross-terms

**Term (5)**: E[L_1(d1) * Q(d2)] = sum_{b in T_2} E[Y_{alpha_1} * Y_b * Y_{(b-d2) mod p}].

For each b in T_2, the three indices are {alpha_1, b, (b-d2) mod p}. Since b != d2, we have b != (b-d2) (as d2 != 0, and b-d2 != 0 since b != d2 in T_2). So the pair {b, (b-d2)} always contributes 2 distinct indices. The triple has 3 distinct indices unless alpha_1 coincides with b or with (b-d2) mod p.

Collisions:
- alpha_1 = b: occurs for exactly one value b = alpha_1 = p - d1 (if alpha_1 in T_2, i.e., alpha_1 != d2, i.e., p - d1 != d2, i.e., d1 + d2 != p).
- alpha_1 = (b-d2) mod p: occurs for b = alpha_1 + d2 = p - d1 + d2 mod p (if this b is in T_2, i.e., b != d2, i.e., p - d1 != 0 mod p, which holds since d1 != p = 0 mod p).
- Both simultaneously: alpha_1 = b and alpha_1 = (b-d2), so d2 = 0 (excluded).

For non-complementary pairs (d1 + d2 != 0):
- alpha_1 = p - d1 != d2 (since d1 + d2 != 0 mod p), so alpha_1 in T_2: **one collision at b = alpha_1**.
- b = p - d1 + d2: need to check if in T_2. This b = d2 iff p - d1 = 0 iff d1 = 0 (excluded). Also, b = 0 mod p iff d1 = d2 (excluded). So b in T_2: **one collision at b = p - d1 + d2**.
- These two collision values of b are distinct (alpha_1 = p-d1 vs p-d1+d2; equal iff d2 = 0).

So there are exactly **2 values of b** giving 2 distinct indices (contributing q_2 each), and **(p-2) - 2 = p - 4 values of b** giving 3 distinct indices (contributing q_3 each). But we need to verify that b = alpha_1 and b = p-d1+d2 are both in T_2 and are in {1,...,p-1}:

- b = alpha_1 = p - d1: in {1,...,p-1} since d1 in {1,...,p-1}. In T_2 iff p-d1 != d2, i.e., d1+d2 != 0 mod p (our assumption). Check.
- b = (p-d1+d2) mod p: if d1 != d2, this is nonzero (since d1 != d2 means p-d1+d2 != p = 0 mod p). In T_2 iff b != d2, i.e., p-d1 != 0 mod p (true). Check.

However, we also need b in {1,...,p-1}, which means b != 0. For b = (p-d1+d2) mod p: this is 0 iff d1 = d2 (excluded). Check.

Also, we excluded b = d2 from T_2: the collision value b = p-d1+d2 equals d2 iff p-d1 = 0 iff d1 = 0 (excluded). And b = p-d1 equals d2 iff d1+d2 = 0 (excluded). Check.

**Conclusion for Term (5)**:
  E_5 = 2 * q_2 + (p - 4) * q_3.

By an identical argument (swapping roles):

**Term (6)**: E[L_2(d1) * Q(d2)] = sum_{b in T_2} E[Y_{beta_1} * Y_b * Y_{(b-d2)}].
Collisions at b = beta_1 = d1 and b = d1 + d2 mod p.
- b = d1: in T_2 iff d1 != d2 (true). In {1,...,p-1} since d1 >= 1. Check.
- b = d1 + d2 mod p: in T_2 iff d1 + d2 != d2 mod p, i.e., d1 != 0 (true). Nonzero iff d1 + d2 != 0 (our assumption). Check.

So: E_6 = 2 * q_2 + (p - 4) * q_3.

**Term (7)**: E[Q(d1) * L_1(d2)] = sum_{a in T_1} E[Y_a * Y_{(a-d1)} * Y_{alpha_2}].
Collisions of alpha_2 = p-d2 with a or with (a-d1):
- alpha_2 = a: a = p-d2, in T_1 iff p-d2 != d1 iff d1+d2 != 0 (our assumption). Check.
- alpha_2 = (a-d1): a = p-d2+d1, in T_1 iff a != d1, i.e., p-d2 != 0, i.e., d2 != 0 (true). Also a in {1,...,p-1}: a = 0 iff d1 = d2 (excluded). Check.

So: E_7 = 2 * q_2 + (p - 4) * q_3.

**Term (8)**: E[Q(d1) * L_2(d2)] = sum_{a in T_1} E[Y_a * Y_{(a-d1)} * Y_{beta_2}].
Collisions of beta_2 = d2 with a or (a-d1):
- beta_2 = a: a = d2, in T_1 iff d2 != d1 (true). Check.
- beta_2 = (a-d1): a = d2+d1, in T_1 iff d1+d2 != d1 iff d2 != 0 (true). Also a != 0 iff d1+d2 != 0 (our assumption). Check.

So: E_8 = 2 * q_2 + (p - 4) * q_3.

**Summary for Terms (5)-(8)**: Each equals 2*q_2 + (p-4)*q_3, and their sum is:

  E_5 + E_6 + E_7 + E_8 = 4 * [2*q_2 + (p-4)*q_3] = 8*q_2 + 4*(p-4)*q_3.

This is **independent of d1, d2** (for non-complementary pairs).

### 4.3 Term (9): Q(d1) * Q(d2) (the double sum)

  E_9 = sum_{a in T_1} sum_{b in T_2} E[Y_a * Y_{(a-d1)} * Y_b * Y_{(b-d2)}]

where T_1 = {1,...,p-1} \ {d1}, T_2 = {1,...,p-1} \ {d2}.

For each (a,b) pair, the four indices are {a, a' = (a-d1) mod p, b, b' = (b-d2) mod p}. We always have a != a' (since d1 != 0) and b != b' (since d2 != 0). The expectation depends on the number of distinct indices.

**Possible cross-equalities**: The only equalities that can arise among {a, a', b, b'} beyond a != a' and b != b' are:
- a = b
- a = b' (i.e., a = b - d2, so b = a + d2)
- a' = b (i.e., a - d1 = b, so b = a - d1)
- a' = b' (i.e., a - d1 = b - d2, so b = a + d2 - d1)

Note that a = a' and b = b' are impossible (d1, d2 != 0). The four equalities correspond to four values of the offset delta = b - a (mod p):

- delta_1 = 0  (a = b)
- delta_2 = d2  (a = b')
- delta_3 = -d1 mod p = p - d1  (a' = b)
- delta_4 = d2 - d1 mod p  (a' = b')

**When do two of these deltas coincide?**

| delta_i = delta_j | Condition | Status |
|---|---|---|
| delta_1 = delta_2 | d2 = 0 | Excluded |
| delta_1 = delta_3 | d1 = 0 | Excluded |
| delta_1 = delta_4 | d1 = d2 | Excluded |
| delta_2 = delta_3 | d1 + d2 = 0 mod p | **Excluded for non-complementary** |
| delta_2 = delta_4 | d1 = 0 | Excluded |
| delta_3 = delta_4 | d2 = 0 | Excluded |

**Conclusion**: For non-complementary pairs (d1 + d2 != 0, d1 != d2, d1 != 0, d2 != 0), the four collision offsets delta_1, delta_2, delta_3, delta_4 are **all distinct**.

### Counting (a,b) pairs for each collision offset

For a given offset delta, the collision pairs are {(a, b) : a in T_1, b in T_2, b = a + delta mod p}. For each a in T_1, there is exactly one b = (a + delta) mod p, and this b is valid iff b in T_2, i.e., b != d2 and b in {1,...,p-1} (so b != 0).

**Count of valid a's**: Start with |T_1| = p-2 values of a. Remove those where b = 0 or b = d2:
- b = 0: a = -delta mod p. This a is in T_1 iff -delta != d1 mod p, i.e., delta != -d1 = p-d1 mod p.
- b = d2: a = d2 - delta mod p. This a is in T_1 iff d2 - delta != d1 mod p, i.e., delta != d2 - d1 mod p.

We also need a in {1,...,p-1}, which is given since a in T_1 subset {1,...,p-1}.

But we also need to check that b != 0 actually removes an element from T_1: it does iff a = (-delta) mod p is in T_1, i.e., a != d1 (i.e., -delta != d1, i.e., delta != p-d1 = -d1) AND a in {1,...,p-1} (i.e., -delta != 0 mod p, i.e., delta != 0).

Let us carefully count for each of the four collision offsets:

**Offset delta_1 = 0** (collision: a = b):
- b = a, so b = 0 iff a = 0 (but a in T_1 means a >= 1).
- b = d2 iff a = d2. a = d2 is in T_1 iff d2 != d1 (true).
- So we lose 1 value (a = d2), giving **p - 3** valid pairs.
- Number of distinct indices: a = b, a' = a-d1, b' = b-d2 = a-d2. With a = b, the set is {a, a-d1, a-d2}. These are 3 distinct values iff d1 != d2 (true) and a-d1 != a (true, d1 != 0) and a-d2 != a (true, d2 != 0). Actually we also need a-d1 != a-d2, which requires d1 != d2 (true). So: **3 distinct indices -> q_3**.

**Offset delta_2 = d2** (collision: a = b'):
- b = a + d2. b = 0 iff a = -d2 = p-d2 mod p. Then a = p-d2 is in T_1 iff p-d2 != d1, i.e., d1+d2 != 0 mod p (our assumption). So we lose this a.
- b = d2 iff a = 0, but a in T_1 so a >= 1. Not lost.
- So we lose 1 value (a = p-d2), giving **p - 3** valid pairs.
- Distinct indices: {a, a-d1, a+d2, a+d2-d2} = {a, a-d1, a+d2, a} = {a, a-d1, a+d2}. Wait -- b' = b - d2 = (a+d2) - d2 = a. So the collision is a = b'. The four indices become {a, a-d1, b, a} = {a, a-d1, b} = {a, a-d1, a+d2}. These are 3 distinct iff a != a-d1 (true) and a != a+d2 (true, d2 != 0) and a-d1 != a+d2 (i.e., d1+d2 != 0, which is our assumption). So: **3 distinct indices -> q_3**.

**Offset delta_3 = p-d1** (collision: a' = b):
- b = a + (p-d1) = a - d1 mod p. b = 0 iff a = d1, but a = d1 is excluded from T_1. So no loss from b = 0.
- b = d2 iff a - d1 = d2 iff a = d1 + d2 mod p. This a is in T_1 iff d1+d2 != d1, i.e., d2 != 0 (true). Also a = d1+d2 != 0 since d1+d2 != 0 mod p (our assumption). But a must be in {1,...,p-1}: d1+d2 mod p is in {1,...,p-1} since d1+d2 != 0. In T_1 iff d1+d2 != d1 (i.e., d2 != 0). So we lose this a.
- So we lose 1 value (a = d1+d2), giving **p - 3** valid pairs.
- Distinct indices: a' = b, so {a, a-d1, a-d1, b-d2} = {a, a-d1, (a-d1)-d2} = {a, a-d1, a-d1-d2}. Three distinct iff a != a-d1 (true), a != a-d1-d2 (i.e., d1+d2 != 0, true), a-d1 != a-d1-d2 (i.e., d2 != 0, true). So: **3 distinct indices -> q_3**.

**Offset delta_4 = d2-d1 mod p** (collision: a' = b'):
- b = a + (d2-d1). b = 0 iff a = d1-d2 mod p. This a is in T_1 iff d1-d2 != d1 (i.e., d2 != 0, true) and d1-d2 != 0 (i.e., d1 != d2, true). So a = (d1-d2) mod p is in T_1 and we lose it.
- b = d2 iff a = d2-(d2-d1) = d1. But a = d1 is excluded from T_1. So no additional loss.
- So we lose 1 value (a = d1-d2 mod p), giving **p - 3** valid pairs.
- Distinct indices: a' = b', so {a, a-d1, b, b-d2} = {a, a-d1, b, a-d1} since a' = a-d1 and b' = b-d2, and a-d1 = b-d2 means b = a-d1+d2. So the indices are {a, a-d1, a-d1+d2}. Three distinct iff a != a-d1 (true), a != a-d1+d2 (i.e., d1 != d2, true), a-d1 != a-d1+d2 (i.e., d2 != 0, true). So: **3 distinct indices -> q_3**.

**However**, we need to check that the four collision offsets are distinct from each other, which we already proved. We also need to check for **overcounting**: can an (a,b) pair satisfy two collision conditions simultaneously? This happens when b-a equals two different delta values simultaneously, which is impossible since b-a has a unique value mod p. So there is no overcounting.

**Generic pairs** (no collision): The total number of (a,b) pairs in T_1 x T_2 is (p-2)^2. The number of collision pairs is at most 4 * (p-2) (at most p-2 pairs per offset, before exclusions). More precisely, the total collision count is 4 * (p-3) = 4p - 12 (each offset contributes p-3 pairs). The generic count is:

  (p-2)^2 - 4(p-3) = p^2 - 4p + 4 - 4p + 12 = p^2 - 8p + 16 = (p-4)^2.

Wait, let me recount. But actually, we need to be more careful: some of the "lost" pairs in the collision offsets may involve a values that are excluded due to being collision values of other offsets. But since the offsets are distinct and b is determined by a, different offsets cannot share the same (a,b) pair. So the total collision pair count is indeed:

  4 * (p - 3) = 4p - 12.

But we also need to subtract the cases where a collision offset hits a **second** collision condition. Since each (a,b) pair has a unique offset b-a, it belongs to at most one collision class. So the total collision count is the sum of individual counts: 4(p-3).

Actually, I need to be more careful. For each offset delta, I counted "p-3 valid pairs." But I should also check: is a = (-delta) mod p always in T_1, and is a = (d2-delta) mod p always in T_1? Let me recheck each offset.

For **delta_1 = 0**: Lost a = d2 (from b = d2 condition). What about b = 0? b = a = 0 requires a = 0, but a >= 1. No loss. Total exclusions from T_1: just a = d1. From validity: a = d2. Since d1 != d2, these are different. So valid a: p-1 (total {1,...,p-1}) minus 1 (a = d1, not in T_1) minus 1 (a = d2, b = d2 excluded) = p-3. Check.

For **delta_2 = d2**: Lost a = p-d2 (from b = 0 condition). b = d2 requires a = 0, excluded. So valid a: p-1 minus 1 (a = d1) minus 1 (a = p-d2) = p-3. But: is a = p-d2 = d1 possible? That would mean d1 + d2 = p = 0 mod p, excluded. So these are distinct removals. Count: p-3. Check.

For **delta_3 = p-d1**: b = 0 requires a = d1, already excluded from T_1. b = d2 requires a = d1+d2. Is a = d1+d2 in T_1? We need d1+d2 != d1 (true since d2 != 0) and d1+d2 in {1,...,p-1} (true since d1+d2 != 0 mod p). So valid a: p-1 minus 1 (a = d1) minus 1 (a = d1+d2) = p-3. But is d1+d2 = d1 possible? Only if d2 = 0, excluded. Check.

For **delta_4 = d2-d1**: b = 0 requires a = d1-d2. Is a = d1-d2 in T_1? Need d1-d2 != d1 (i.e., d2 != 0, true) and d1-d2 != 0 (i.e., d1 != d2, true). So a = d1-d2 is a valid member of T_1, and we remove it. b = d2 requires a = d1, already excluded. So valid a: p-1 minus 1 (a = d1) minus 1 (a = d1-d2) = p-3. Check.

**However**, we also need to verify that each lost a-value is distinct from a = d1 (which is not in T_1). We checked this above in each case. So each offset contributes exactly **p - 3** collision pairs, each with 3 distinct indices.

**Total collision pair count**: 4(p-3).

**Generic pair count**: (p-2)^2 - 4(p-3) = p^2 - 4p + 4 - 4p + 12 = p^2 - 8p + 16.

Generic pairs have 4 distinct indices -> contribute q_4.

**Summary for Term (9)**:

  E_9 = 4(p-3) * q_3 + (p^2 - 8p + 16) * q_4.

This expression depends **only on p**, not on d1 or d2.

---

## 5. Combining All Nine Terms

For non-complementary pairs (d1 + d2 != 0, d1 != d2, d1, d2 in {1,...,p-1}):

  E[B(d1) * B(d2)] = [E_1 + E_2 + E_3 + E_4] + [E_5 + E_6 + E_7 + E_8] + E_9

where:
- Terms (1)-(4): 4 * q_2
- Terms (5)-(8): 8 * q_2 + 4(p-4) * q_3
- Term (9): 4(p-3) * q_3 + (p^2 - 8p + 16) * q_4

Total:

  E[B(d1)*B(d2)] = (4 + 8) q_2 + [4(p-4) + 4(p-3)] q_3 + (p^2 - 8p + 16) q_4
                  = 12 q_2 + (8p - 28) q_3 + (p^2 - 8p + 16) q_4.

**This depends only on p.** Since E[B(d)] = (p-3)/4 is also independent of d, we have:

  Cov[B(d1), B(d2)] = E[B(d1)*B(d2)] - E[B(d1)] * E[B(d2)]
                     = [12 q_2 + (8p-28) q_3 + (p^2-8p+16) q_4] - [(p-3)/4]^2

which is the **same for all non-complementary pairs**. QED (equi-covariance)

---

## 6. Closed-Form Derivation

Now we compute the explicit value. Substituting k = (p-3)/2, N = p-1:

  q_1 = k/N = (p-3) / (2(p-1))

  q_2 = k(k-1) / (N(N-1)) = [(p-3)(p-5)] / [4(p-1)(p-2)]

  q_3 = k(k-1)(k-2) / (N(N-1)(N-2)) = [(p-3)(p-5)(p-7)] / [8(p-1)(p-2)(p-3)]
       = (p-5)(p-7) / [8(p-1)(p-2)]

  q_4 = k(k-1)(k-2)(k-3) / (N(N-1)(N-2)(N-3))
       = [(p-3)(p-5)(p-7)(p-9)] / [16(p-1)(p-2)(p-3)(p-4)]
       = (p-5)(p-7)(p-9) / [16(p-1)(p-2)(p-4)]

**E[B]^2**:

  E[B]^2 = [(p-3)/4]^2 = (p-3)^2 / 16

**Computing E[B(d1)*B(d2)]**:

Denote the three contributions:
- A = 12 q_2 = 12(p-3)(p-5) / [4(p-1)(p-2)] = 3(p-3)(p-5) / [(p-1)(p-2)]
- B_coeff = (8p-28) q_3 = (8p-28)(p-5)(p-7) / [8(p-1)(p-2)] = (p-3.5)(p-5)(p-7) / [(p-1)(p-2)]

Wait, let me be more careful with Fraction arithmetic. It is cleaner to work with a common denominator.

Using N = p-1, the common denominator for q_2, q_3, q_4 is N(N-1)(N-2)(N-3) = (p-1)(p-2)(p-3)(p-4).

  q_2 = k(k-1)(N-2)(N-3) / [N(N-1)(N-2)(N-3)]
  q_3 = k(k-1)(k-2)(N-3) / [N(N-1)(N-2)(N-3)]
  q_4 = k(k-1)(k-2)(k-3) / [N(N-1)(N-2)(N-3)]

Let D = N^{(4)} = (p-1)(p-2)(p-3)(p-4) be the common denominator.

Numerator of the E[B1*B2] sum (times D):

  12 * k(k-1) * (N-2)(N-3) + (8p-28) * k(k-1)(k-2) * (N-3) + (p^2-8p+16) * k(k-1)(k-2)(k-3)

Factor out k(k-1):

  k(k-1) * [12(N-2)(N-3) + (8p-28)(k-2)(N-3) + (p^2-8p+16)(k-2)(k-3)]

Substituting k = (p-3)/2, N = p-1, N-2 = p-3, N-3 = p-4, k-2 = (p-7)/2, k-3 = (p-9)/2:

  k(k-1) = (p-3)(p-5)/4

  12(p-3)(p-4) + (8p-28) * (p-7)/2 * (p-4) + (p-4)^2 * (p-7)(p-9)/4

  = 12(p-3)(p-4) + 4(p-3.5)(p-7)(p-4) + ...

This is getting unwieldy. Instead, let us use the approach of computing Cov from the sum constraint.

### Alternative derivation via sum constraint

Since all non-complementary covariances are equal (proven in Sections 4-5), we call this common value c. The sum constraint and complementary-pair structure give c directly.

**From Fact 2 (Parseval/sum constraint)**: sum_{d=1}^{p-1} B(d) = s(s-1) is constant, where s = |D12| = (p-1)/2. Therefore Var(sum) = 0.

Expanding:

  0 = sum_{d1, d2 = 1}^{p-1} Cov[B(d1), B(d2)]
    = sum_d Var[B(d)] + sum_{d1 != d2} Cov[B(d1), B(d2)]

For each d, there are:
- 1 complementary partner p-d (with Cov = Var since B(d) = B(p-d))
- p-3 non-complementary partners (with Cov = c)

So:

  0 = (p-1) Var + (p-1) Var + (p-1)(p-3) c

(The first term: sum of diagonal = (p-1) * Var. The second: sum over complementary off-diagonal pairs, each d has one complementary partner, giving (p-1) terms of value Var. The third: (p-1)(p-3) non-complementary off-diagonal pairs, each with value c.)

Solving:

  c = -2 Var / (p-3)

Using Var[B(d)] = (p-3)(p+1) / (16(p-2)) (proven in exact_moments.py):

  **c = -2 * (p-3)(p+1) / (16(p-2)(p-3)) = -(p+1) / (8(p-2))**

The non-complementary correlation is:

  **rho = c / Var = -2/(p-3)**

---

## 7. Verification

The formula Cov = -(p+1)/(8(p-2)) has been verified by exact computation (Python Fraction arithmetic) against the brute-force indicator expansion for all primes p in {7, 11, 19, 23, 31, 43, 47, 59, 67, 83, 199, 997}, and for ALL pairs (d1, d2) with d1+d2 != 0 mod p for p = 7, 11, 19, 23. See `verify_equi_covariance.py`.

Additionally, the collision counts derived in Section 4 (specifically: the E_9 collision count of 4(p-3) and the term-by-term structure of E_5 through E_8) have been verified against brute-force enumeration for p = 11 and p = 23.

---

## 8. Summary

The equi-covariance property Cov[B(d1), B(d2)] = const for all non-complementary pairs does **not** require any group-theoretic transitivity argument. It follows from the elementary observation that:

1. E[B(d1)*B(d2)] decomposes into a sum of hypergeometric moments q_j, weighted by the number of terms with j distinct indices.

2. The number of terms at each collision level depends only on whether certain index offsets are zero mod p. For non-complementary pairs (d1+d2 != 0), the conditions d1 != 0, d2 != 0, d1 != d2, and d1+d2 != 0 ensure that all potential collision offsets remain distinct and all collision counts are the same.

3. Therefore E[B(d1)*B(d2)] is the same for all non-complementary pairs.

4. Combined with the sum constraint sum B(d) = const, this yields the closed form Cov = -(p+1)/(8(p-2)).
