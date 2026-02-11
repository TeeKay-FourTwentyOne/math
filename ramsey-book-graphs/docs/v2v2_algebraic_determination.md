# Theorem 2: V₂V₂ Algebraic Determination and the Complete Structural Reduction

## Continuation from the V₁V₂ Auto-Satisfaction Theorem

We retain all notation and constraints (S1)–(S4) from the V₁V₂ theorem. We now analyze the remaining constraint classes — V₁V₁ and V₂V₂ edges — and show that V₂V₂ constraints are algebraically determined by V₁V₁ data.

---

## Notation

Throughout, arithmetic is in ℤ_m with m = 2n − 1. Define the **cross-correlation** of sets A, B ⊆ ℤ_m at shift d:

$$C(A, B;\, d) \;=\; \bigl|\{a \in A : a + d \in B\}\bigr|$$

We write:

$$A(d) \;=\; C(D_{11}, D_{11};\, d), \qquad B(d) \;=\; C(D_{12}, D_{12};\, d)$$

Note that by symmetry (S1), A(d) = A(−d) for all d (proof: C(D₁₁, D₁₁; d) = |{a ∈ D₁₁ : a + d ∈ D₁₁}|; substituting a′ = −(a + d) ∈ D₁₁ when a + d ∈ D₁₁, using −D₁₁ = D₁₁, gives a bijection with {a′ ∈ D₁₁ : a′ − d ∈ D₁₁} = C(D₁₁, D₁₁; −d)).

However, B(d) ≠ B(−d) in general since D₁₂ need not be symmetric.

---

## Recollection of Edge Types

| Edge type | Red condition | Blue condition |
|-----------|--------------|----------------|
| V₁–V₁ at difference d | d ∈ D₁₁ | d ∈ D₂₂ = {1,…,m−1} ∖ D₁₁ |
| V₂–V₂ at difference d | d ∈ D₂₂ | d ∈ D₁₁ |
| V₁–V₂ at difference d | d ∈ D₁₂ | d ∉ D₁₂ |

The complement structure (S2) creates a **red–blue duality** between V₁V₁ and V₂V₂: V₁V₁ edges at difference d ∈ D₁₁ are red while V₂V₂ edges at the same difference are blue, and vice versa.

---

## Part I: V₁V₁ Common Neighbor Counts

### Red V₁V₁ edges (d ∈ D₁₁)

Fix vertices u = a, v = b ∈ V₁ with d = (a − b) mod m, where d ∈ D₁₁.

**Red common neighbors from V₁.** Vertex w ∈ V₁ (w ≠ u, v) contributes iff (w − a) ∈ D₁₁ and (w − b) ∈ D₁₁. Setting s = w − a: require s ∈ D₁₁ and s + d ∈ D₁₁, with s ≠ 0 and s ≠ −d.

Since 0 ∉ D₁₁, s ∈ D₁₁ implies s ≠ 0. Similarly s + d ∈ D₁₁ implies s + d ≠ 0, hence s ≠ −d. Both exclusions are automatic.

$$\text{Count from } V_1 = A(d)$$

**Red common neighbors from V₂.** Vertex w′ ∈ V₂ contributes iff edge {u, w′} and edge {v, w′} are both red. By our convention, {a, w′} is red iff (a − w) ∈ D₁₂, and {b, w′} is red iff (b − w) ∈ D₁₂. Setting t = a − w: require t ∈ D₁₂ and t − d ∈ D₁₂.

$$\text{Count from } V_2 = C(D_{12}, D_{12};\, {-d}) = B(-d)$$

**Total:**

$$\boxed{R_{11}(d) \;=\; A(d) + B(-d)}$$

**Constraint (R):** For red V₁V₁ edges (d ∈ D₁₁): $A(d) + B(-d) \leq n - 2$.

---

### Blue V₁V₁ edges (d ∈ D₂₂, i.e., d ∉ D₁₁, d ≠ 0)

**Blue common neighbors from V₁.** Vertex w ∈ V₁ (w ≠ u, v) contributes iff s = w − a satisfies: s ∉ D₁₁ (blue within V₁), s + d ∉ D₁₁ (blue), s ≠ 0, s + d ≠ 0.

For s ≠ 0 and s ∉ D₁₁: this means s ∈ D₂₂. Similarly s + d ≠ 0 and s + d ∉ D₁₁ means s + d ∈ D₂₂.

$$\text{Count from } V_1 = C(D_{22}, D_{22};\, d)$$

**Blue common neighbors from V₂.** Vertex w′ ∈ V₂ contributes iff (a − w) ∉ D₁₂ and (b − w) ∉ D₁₂. Setting t = a − w: t ∉ D₁₂ and t − d ∉ D₁₂.

$$\text{Count from } V_2 = m - 2|D_{12}| + B(-d)$$

(by inclusion–exclusion: total m elements minus those with t ∈ D₁₂ minus those with t − d ∈ D₁₂ plus those with both)

**Constraint (B):** Blue common neighbors ≤ n − 1.

We will not need to expand these further because the key result is that V₂V₂ is *determined* by V₁V₁ data — the actual computational constraints that must be checked live in V₁V₁.

---

## Part II: V₂V₂ Common Neighbor Counts

### Lemma (Complement Autocorrelation)

For D₂₂ = {1, …, m−1} ∖ D₁₁ and any nonzero d:

$$C(D_{22}, D_{22};\, d) \;=\; A(d) + (m - 2) - 2|D_{11}| + 2\cdot\mathbb{1}[d \in D_{11}]$$

**Proof.** Define indicator functions on ℤ_m: let f(x) = 𝟙[x ∈ D₁₁] and g(x) = 𝟙[x ∈ D₂₂]. For all x: g(x) = 𝟙[x ≠ 0] − f(x), since D₁₁ and D₂₂ partition {1, …, m−1} and both exclude 0.

Expand:

$$C(D_{22}, D_{22};\, d) = \sum_{s \in \mathbb{Z}_m} g(s)\,g(s+d)$$

$$= \sum_s \bigl(\mathbb{1}[s \neq 0] - f(s)\bigr)\bigl(\mathbb{1}[s{+}d \neq 0] - f(s{+}d)\bigr)$$

$$= T_1 - T_2 - T_3 + T_4$$

where:

**$T_1$** $= \sum_s \mathbb{1}[s \neq 0]\cdot\mathbb{1}[s+d \neq 0]$. For d ≠ 0, we exclude s = 0 and s = −d (distinct), giving $T_1 = m - 2$.

**$T_2$** $= \sum_s \mathbb{1}[s \neq 0] \cdot f(s+d) = \sum_{s \neq 0} f(s+d)$. Substituting t = s + d: as s ranges over ℤ_m ∖ {0}, t ranges over ℤ_m ∖ {d}. So $T_2 = |D_{11}| - \mathbb{1}[d \in D_{11}]$.

**$T_3$** $= \sum_s f(s) \cdot \mathbb{1}[s+d \neq 0] = \sum_{s \neq -d} f(s) = |D_{11}| - f(-d) = |D_{11}| - \mathbb{1}[d \in D_{11}]$,

where the last step uses symmetry (S1): $f(-d) = \mathbb{1}[-d \in D_{11}] = \mathbb{1}[d \in D_{11}]$.

**$T_4$** $= \sum_s f(s)\,f(s+d) = A(d)$.

Combining: $C(D_{22}, D_{22};\, d) = (m-2) - 2\bigl(|D_{11}| - \mathbb{1}[d \in D_{11}]\bigr) + A(d)$. ∎

---

### Theorem 2 (V₂V₂ Algebraic Determination)

For any nonzero d ∈ ℤ_m, the red common neighbor count for the V₂V₂ edge at difference d is:

$$\boxed{R_{22}(d) \;=\; A(d) + B(d) + (m - 2 - 2|D_{11}|) + 2\cdot\mathbb{1}[d \in D_{11}]}$$

**Proof.** Fix vertices u = a′, v = b′ ∈ V₂ with d = (a − b) mod m.

**Red common neighbors from V₁.** Vertex w ∈ V₁ contributes iff (w − a) ∈ D₁₂ and (w − b) ∈ D₁₂. Setting s = w − a: require s ∈ D₁₂ and s + d ∈ D₁₂.

$$\text{Count from } V_1 = C(D_{12}, D_{12};\, d) = B(d)$$

**Red common neighbors from V₂.** Vertex w′ ∈ V₂ (w′ ≠ u, v) contributes iff (w − a) ∈ D₂₂ and (w − b) ∈ D₂₂. Setting s = w − a: require s ∈ D₂₂ and s + d ∈ D₂₂, with s ≠ 0 and s + d ≠ 0. As in Part I, both exclusions are automatic since 0 ∉ D₂₂.

$$\text{Count from } V_2 = C(D_{22}, D_{22};\, d) = A(d) + (m-2-2|D_{11}|) + 2\cdot\mathbb{1}[d \in D_{11}]$$

by the Complement Autocorrelation Lemma. Adding yields the result. ∎

---

## Part III: The Duality and Complete Reduction

### Corollary (Red–Blue Duality)

By the complement structure, the four constraint classes reduce to two independent conditions involving $A(d)$ and $B(d)$:

| Constraint class | Difference d in | Condition |
|-----------------|-----------------|-----------|
| Red V₁V₁ | D₁₁ | $A(d) + B(-d) \leq n - 2$ |
| Blue V₂V₂ | D₁₁ | $1 + A(d) + B(d) \leq n - 1$ |
| Blue V₁V₁ | D₂₂ | (see below) |
| Red V₂V₂ | D₂₂ | $A(d) + B(d) + (m - 2 - 2|D_{11}|) \leq n - 2$ |

**Derivation of Blue V₂V₂ condition.** For d ∈ D₁₁ (blue V₂V₂ edge), we count blue common neighbors:

From V₁: vertex w contributes iff (w − a) ∉ D₁₂ and (w − b) ∉ D₁₂. By inclusion–exclusion:

$$\text{Count from } V_1 = m - 2|D_{12}| + B(d)$$

From V₂: vertex w′ (w′ ≠ u, v) contributes iff (w − a) ∉ D₂₂ and (w − b) ∉ D₂₂. For w ≠ a, b: this requires (w − a) ∈ D₁₁ and (w − b) ∈ D₁₁, giving count A(d).

Total blue CN = $(m - 2|D_{12}| + B(d)) + A(d) = m - 2(n-1) + A(d) + B(d) = 1 + A(d) + B(d)$

since $m - 2(n-1) = (2n-1) - 2n + 2 = 1$. Constraint: $\leq n - 1$, i.e., $A(d) + B(d) \leq n - 2$.

---

### Theorem 3 (Complete Structural Reduction)

**Under constraints (S1)–(S3), the problem of constructing a valid 2-block circulant for R(B_{n−1}, B_n) ≥ 4n − 1 reduces to finding D₁₁, D₁₂ satisfying:**

**(I)** For all $d \in D_{11}$:

$$A(d) + B(-d) \leq n - 2 \qquad\text{and}\qquad A(d) + B(d) \leq n - 2$$

**(II)** For all $d \in D_{22}$:

$$A(d) + B(d) + (m - 2 - 2|D_{11}|) \leq n - 2$$

together with the blue V₁V₁ constraint (which we derive below is also a function of A, B only).

**Proof.** By the V₁V₂ Auto-Satisfaction Theorem, all V₁V₂ constraints are identically satisfied. The V₁V₁ and V₂V₂ constraints are enumerated in the corollary above. Every constraint involves only the functions $A(d) = C(D_{11}, D_{11};\, d)$ and $B(d) = C(D_{12}, D_{12};\, d)$, which are the autocorrelations of the two design parameters $D_{11}$ and $D_{12}$. ∎

---

### Remark: Blue V₁V₁ Constraint

For completeness, the blue V₁V₁ constraint (d ∈ D₂₂) is:

**Blue CN from V₁** = $C(D_{22}, D_{22};\, d) = A(d) + (m - 2 - 2|D_{11}|)$ (the indicator $\mathbb{1}[d \in D_{11}] = 0$)

**Blue CN from V₂** = $m - 2|D_{12}| + B(-d) = 1 + B(-d)$

**Total** = $A(d) + B(-d) + (m - 1 - 2|D_{11}|)$

**Constraint**: $\leq n - 1$.

---

## Summary: The Four Constraints

For a valid 2-block circulant under (S1)–(S3), the necessary and sufficient conditions (beyond V₁V₂, which is free) are:

| # | For d ∈ | Constraint | Origin |
|---|---------|-----------|--------|
| (C1) | $D_{11}$ | $A(d) + B(-d) \leq n - 2$ | Red V₁V₁ |
| (C2) | $D_{11}$ | $A(d) + B(d) \leq n - 2$ | Blue V₂V₂ |
| (C3) | $D_{22}$ | $A(d) + B(d) \leq n - 2 - (m - 2 - 2|D_{11}|)$ | Red V₂V₂ |
| (C4) | $D_{22}$ | $A(d) + B(-d) \leq n - 1 - (m - 1 - 2|D_{11}|)$ | Blue V₁V₁ |

Substituting $m = 2n - 1$ and writing $k = |D_{11}|$:

| # | For d ∈ | Simplified constraint |
|---|---------|----------------------|
| (C1) | $D_{11}$ | $A(d) + B(-d) \leq n - 2$ |
| (C2) | $D_{11}$ | $A(d) + B(d) \leq n - 2$ |
| (C3) | $D_{22}$ | $A(d) + B(d) \leq 2k - n + 2$ |
| (C4) | $D_{22}$ | $A(d) + B(-d) \leq 2k - n + 1$ |

Observe the structure: **(C1)** and **(C2)** are "tight" constraints on $D_{11}$ differences, while **(C3)** and **(C4)** on $D_{22}$ differences have thresholds that depend on $|D_{11}|$. When $|D_{11}| = n - 1$ (the typical case, matching the Paley construction where $|QR| = (q-1)/2 = n - 1$), the thresholds become $n$ and $n - 1$ respectively — providing more room than (C1)/(C2).

**In the Paley case** ($k = n - 1$): (C3) becomes $A(d) + B(d) \leq n$, and (C4) becomes $A(d) + B(-d) \leq n - 1$. The binding constraints are always (C1) and (C2), which are the hardest to satisfy.

---

## Algebraic Core (Self-Contained Summary)

The entire structural theory rests on three algebraic facts:

1. **Complement Partition**: $D_{11} \sqcup D_{22} = \{1, \ldots, m{-}1\}$ converts any summation over $D_{22}$ into a summation over all nonzero elements minus $D_{11}$.

2. **Symmetry of $D_{11}$**: The condition $D_{11} = -D_{11}$ provides $A(d) = A(-d)$ and, crucially, enables the cross-correlation identity $C(D_{12}, D_{11};\, d) = C(D_{11}, D_{12};\, -d)$ that drives the V₁V₂ cancellation.

3. **Dimension count**: $|D_{12}| = n - 1 = (m - 1)/2$ ensures that the V₁V₂ red and blue counts land exactly on their respective thresholds, with no slack.

The V₁V₂ theorem (exact cancellation) and the V₂V₂ determination (complement autocorrelation) are both consequences of (1)–(3), involving no additional assumptions.
