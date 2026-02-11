# The V₁V₂ Auto-Satisfaction Theorem for 2-Block Circulant Book Ramsey Constructions

## Setup and Notation

Let m ≥ 3 be odd. We consider a 2-coloring (red/blue) of the complete graph K_{2m} using a **2-block circulant** structure:

- **Vertex set**: V = V₁ ∪ V₂, where V₁ = {0, 1, …, m−1} and V₂ = {0′, 1′, …, (m−1)′}.
- **Difference sets**: D₁₁, D₁₂, D₂₂ ⊆ {0, 1, …, m−1}, with arithmetic in ℤ_m.
- **Edge coloring**:
  - V₁–V₁: edge {i, j} is **red** iff (i − j) mod m ∈ D₁₁
  - V₁–V₂: edge {i, j′} is **red** iff (i − j) mod m ∈ D₁₂
  - V₂–V₂: edge {i′, j′} is **red** iff (i − j) mod m ∈ D₂₂
  - All other edges are **blue**.

The **book graph** B_n = K₂ + K̄_n consists of n triangles sharing a common edge (2n + 2 vertices). To show R(B_{n−1}, B_n) ≥ 4n − 1, we must exhibit a 2-coloring of K_{4n−2} containing no red B_{n−1} and no blue B_n. The validity conditions are:

> **(R)** Every red edge has at most n − 2 red common neighbors.
>
> **(B)** Every blue edge has at most n − 1 blue common neighbors.

Here a **red common neighbor** of edge {u, v} is a vertex w such that {w, u} and {w, v} are both red (and similarly for blue).

We impose the following **structural constraints** on the difference sets, with n = (m + 1)/2:

| Label | Constraint | Consequence |
|-------|-----------|-------------|
| (S1) | D₁₁ ⊆ {1, …, m−1} is symmetric: x ∈ D₁₁ ⟺ −x ∈ D₁₁ | Within-block adjacency is undirected |
| (S2) | D₂₂ = {1, …, m−1} ∖ D₁₁ | D₂₂ is the complement of D₁₁ among nonzero elements |
| (S3) | \|D₁₂\| = n − 1 | Cross-block red degree is controlled |
| (S4) | 0 ∈ D₁₂ | Convention (applicable when m ≡ 3 mod 4) |

Note that 0 ∉ D₁₁ and 0 ∉ D₂₂ since vertices are not self-adjacent.

---

## Theorem (V₁V₂ Auto-Satisfaction)

**Under constraints (S1)–(S3), every V₁–V₂ edge automatically satisfies both validity conditions (R) and (B). Specifically:**

**(a)** For any d ∈ ℤ_m, the V₁–V₂ edge at difference d has exactly |D₁₂| − 𝟙[d ∈ D₁₂] red common neighbors.

**(b)** Consequently:
  - Every **red** V₁–V₂ edge (d ∈ D₁₂) has exactly n − 2 red common neighbors. ✓ (R)
  - Every **blue** V₁–V₂ edge (d ∉ D₁₂) has exactly n − 1 blue common neighbors. ✓ (B)

**Therefore, the problem of constructing a valid 2-block circulant R(B_{n−1}, B_n) lower bound reduces entirely to satisfying constraints (R) and (B) on V₁–V₁ and V₂–V₂ edges.**

---

## Proof

### Part (a): Red common neighbor count

Fix vertices u = a ∈ V₁ and v = b′ ∈ V₂ with difference d = (a − b) mod m. We count red common neighbors from each block.

**Red common neighbors in V₁.** A vertex w ∈ V₁ (w ≠ u) is a red common neighbor iff:
- {w, u} is red: (w − a) mod m ∈ D₁₁
- {w, v} is red: (w − b) mod m ∈ D₁₂

Setting s = w − a, the conditions become s ∈ D₁₁ and s + d ∈ D₁₂. Define:

$$\sigma_1(d) \;=\; \bigl|\{s \in D_{11} : s + d \in D_{12}\}\bigr|$$

**Red common neighbors in V₂.** A vertex w′ ∈ V₂ (w′ ≠ v) is a red common neighbor iff:
- {u, w′} is red: (a − w) mod m ∈ D₁₂
- {w′, v} is red: (w − b) mod m ∈ D₂₂

Setting t = a − w, the conditions become t ∈ D₁₂ and d − t ∈ D₂₂. Define:

$$\sigma_2(d) \;=\; \bigl|\{t \in D_{12} : d - t \in D_{22}\}\bigr|$$

The total red common neighbor count is σ₁(d) + σ₂(d). We now show this equals |D₁₂| − 𝟙[d ∈ D₁₂].

**Step 1: Expand σ₂ using the complement structure (S2).**

Since D₂₂ = {1, …, m−1} ∖ D₁₁, for any nonzero element x we have x ∈ D₂₂ ⟺ x ∉ D₁₁. Note that d − t ∈ D₂₂ forces d − t ≠ 0, hence t ≠ d. Conversely, for t ∈ D₁₂ with t ≠ d, the nonzero element d − t lies in exactly one of D₁₁ or D₂₂. Therefore:

$$\sigma_2(d) \;=\; \bigl|\{t \in D_{12} : t \neq d \text{ and } d - t \notin D_{11}\}\bigr|$$

$$= \bigl|\{t \in D_{12} : t \neq d\}\bigr| \;-\; \bigl|\{t \in D_{12} : t \neq d \text{ and } d - t \in D_{11}\}\bigr|$$

The first term equals |D₁₂| − 𝟙[d ∈ D₁₂].

For the second term, since 0 ∉ D₁₁, the condition d − t ∈ D₁₁ already implies d − t ≠ 0, hence t ≠ d. So the constraint t ≠ d is redundant:

$$\bigl|\{t \in D_{12} : t \neq d,\; d - t \in D_{11}\}\bigr| \;=\; \bigl|\{t \in D_{12} : d - t \in D_{11}\}\bigr|$$

**Step 2: Apply symmetry of D₁₁ (S1).**

By (S1), d − t ∈ D₁₁ if and only if t − d ∈ D₁₁. So:

$$\bigl|\{t \in D_{12} : d - t \in D_{11}\}\bigr| \;=\; \bigl|\{t \in D_{12} : t - d \in D_{11}\}\bigr|$$

**Step 3: Identify with σ₁(d).**

Substituting s = t − d (so t = s + d), the set {t ∈ D₁₂ : t − d ∈ D₁₁} becomes {s + d : s ∈ D₁₁, s + d ∈ D₁₂}, which bijects with {s ∈ D₁₁ : s + d ∈ D₁₂}. This is exactly σ₁(d).

**Step 4: Combine.**

$$\sigma_2(d) = \bigl(|D_{12}| - \mathbb{1}[d \in D_{12}]\bigr) - \sigma_1(d)$$

Therefore:

$$\boxed{\sigma_1(d) + \sigma_2(d) = |D_{12}| - \mathbb{1}[d \in D_{12}]}$$

This completes the proof of part (a). ∎

---

### Part (b): Validity of red V₁V₂ edges

If d ∈ D₁₂ (red edge), then by part (a):

$$\text{red common neighbors} = |D_{12}| - 1 = (n - 1) - 1 = n - 2 \;\leq\; n - 2 \quad\checkmark\text{ (R)}$$

---

### Part (c): Validity of blue V₁V₂ edges

For a blue V₁V₂ edge (d ∉ D₁₂), we must show the **blue** common neighbor count equals n − 1. We count directly.

**Blue common neighbors in V₁.** Vertex w ∈ V₁ (w ≠ u) is a blue common neighbor iff:
- {w, u} is blue: s = w − a ∈ D₂₂ (i.e., s ≠ 0 and s ∉ D₁₁)
- {w, v} is blue: s + d ∉ D₁₂

So the count is:

$$\beta_1(d) = |D_{22}| - \bigl|\{s \in D_{22} : s + d \in D_{12}\}\bigr|$$

**Blue common neighbors in V₂.** Vertex w′ ∈ V₂ (w′ ≠ v) is a blue common neighbor iff:
- {u, w′} is blue: t = a − w ∉ D₁₂
- {w′, v} is blue: d − t ∉ D₂₂ and d − t ≠ 0

Since d ∉ D₁₂ and we need t ∉ D₁₂, setting t = d would give d − t = 0, which is neither in D₁₁ nor D₂₂. For t ≠ d with d − t ≠ 0: the condition d − t ∉ D₂₂ is equivalent to d − t ∈ D₁₁ (for nonzero elements). By symmetry (S1), d − t ∈ D₁₁ ⟺ t − d ∈ D₁₁. Substituting s = t − d:

$$\beta_2(d) = \bigl|\{s \in D_{11} : s + d \notin D_{12}\}\bigr| = |D_{11}| - \sigma_1(d)$$

(Here the condition t ∉ D₁₂ translates to s + d ∉ D₁₂ being replaced by the complementary condition; see below for the clean resolution.)

**Note on the t = d term:** When t = d, we need d ∉ D₁₂ (already given) and d − d = 0, which means {w′, v} = {b′, b′} is not an edge. So t = d does not contribute. When d − t ∈ D₁₁ (which forces d − t ≠ 0, hence t ≠ d), the constraint t ≠ d is automatic.

**Combining via partition.**

We compute β₁(d) + β₂(d). The key observation is that the "lost" terms partition cleanly:

$$\sigma_1(d) + \bigl|\{s \in D_{22} : s + d \in D_{12}\}\bigr| = \bigl|\{s \in D_{11} : s + d \in D_{12}\}\bigr| + \bigl|\{s \in D_{22} : s + d \in D_{12}\}\bigr|$$

$$= \bigl|\{s \in \{1,\ldots,m{-}1\} : s + d \in D_{12}\}\bigr|$$

Since d ∉ D₁₂, the element s = 0 satisfies 0 + d = d ∉ D₁₂, so it does not contribute. Therefore:

$$= \bigl|\{s \in \mathbb{Z}_m : s + d \in D_{12}\}\bigr| = |D_{12}|$$

So |{s ∈ D₂₂ : s + d ∈ D₁₂}| = |D₁₂| − σ₁(d), and:

$$\beta_1(d) = |D_{22}| - |D_{12}| + \sigma_1(d)$$

$$\beta_2(d) = |D_{11}| - \sigma_1(d)$$

$$\beta_1(d) + \beta_2(d) = |D_{22}| + |D_{11}| - |D_{12}| = (m - 1) - (n - 1) = m - n$$

Since m = 2n − 1:

$$\boxed{\beta_1(d) + \beta_2(d) = (2n - 1) - n = n - 1 \;\leq\; n - 1 \quad\checkmark\text{ (B)}}$$

This completes the proof. ∎

---

## Structural Consequence

The theorem shows that for any 2-block circulant satisfying (S1)–(S3), the V₁–V₂ constraints are **identically satisfied at threshold** — there is no freedom or slack. This has two important implications:

1. **Dimensional reduction**: The search space for valid constructions is restricted entirely to satisfying (R) and (B) on V₁–V₁ and V₂–V₂ edges. By Theorem 2 of the proof outline (algebraic determination of V₂₂ constraints from V₁₁), this further reduces to controlling the autocorrelation Δ(D₁₁, D₁₁, d) and the cross-correlation Δ(D₁₂, D₁₂, d).

2. **Extremal tightness at the boundary**: The V₁–V₂ edges achieve *exactly* the maximum allowable common neighbor counts. This is consistent with the empirical observation that all known constructions are extremally tight — operating with zero slack across all edge types.

---

## Algebraic Core

The entire proof rests on a single identity. Define the **cross-correlation** of sets A, B ⊆ ℤ_m at shift d:

$$C(A, B; d) = \bigl|\{a \in A : a + d \in B\}\bigr|$$

Then:

- **Complement partition**: For 0 ∉ A and D₂₂ = \{1,…,m−1\} ∖ D₁₁, and any d, the nonzero elements s with s + d ∈ B are partitioned by membership in D₁₁ vs. D₂₂.

- **Symmetry–bijection lemma**: If A is symmetric (A = −A), then C(B, A; −d) = C(A, B; d).

  *Proof*: C(B, A; −d) = |{b ∈ B : b − d ∈ A}|. Setting a = b − d: b = a + d, and b − d = a ∈ A ⟺ −a ∈ A (symmetry) ⟺ a ∈ A. So C(B, A; −d) = |{a ∈ A : a + d ∈ B}| = C(A, B; d). ∎

The cancellation σ₁(d) + σ₂(d) = |D₁₂| − 𝟙[d ∈ D₁₂] is then a direct consequence: σ₂ decomposes into a "total minus overlap" where the overlap is exactly σ₁, thanks to the symmetry–bijection lemma identifying the two cross-correlation terms.
