# Golay Complementary Pair Literature: Connection to Our Spectral Complementarity

## Summary of Findings

Our valid (D11, D12) pairs exhibit **spectral complementarity**: the combined power spectrum |hat(1_{D11})|^2 + |hat(1_{D12})|^2 is approximately constant across all nonzero frequencies. This property is the defining characteristic of **Golay complementary pairs** (exact case) or **near-complementary pairs** (approximate case). This document surveys the relevant literature and assesses its applicability to our proof of R(B_{n-1}, B_n) = 4n-1.

---

## 1. Golay Complementary Pairs: Definitions

### 1.1 Aperiodic Golay Pairs (Golay, 1961)
Two sequences A = (a_0, ..., a_{N-1}) and B = (b_0, ..., b_{N-1}) of length N (with entries in {+1, -1}) form a **Golay complementary pair** if their aperiodic autocorrelation functions sum to zero for all nonzero shifts:

  C_A(s) + C_B(s) = 0 for all s != 0

where C_A(s) = sum_{k=0}^{N-1-s} a_k a_{k+s}.

Equivalently, in the spectral domain: |hat(A)(k)|^2 + |hat(B)(k)|^2 = 2N for all k.

This is exactly our spectral complementarity condition.

### 1.2 Periodic Golay Pairs
Two {+/-1}-sequences of even length v form a **periodic Golay pair** if their periodic autocorrelation functions sum to zero for all nonzero shifts:

  PAF(A,s) + PAF(B,s) = 0 for s = 1, ..., v-1

where PAF(A,s) = sum_{k=0}^{v-1} a_k a_{k+s mod v}.

**Critical equivalence**: Periodic Golay pairs of length v are equivalent to **2-{v; k1, k2; lambda} supplementary difference sets** in the cyclic group Z_v.

### 1.3 Our Setting
We work with indicator functions of subsets D11, D12 of Z_p (entries in {0,1}, not {+/-1}). Our spectral complementarity is:

  |hat(1_{D11})(k)|^2 + |hat(1_{D12})(k)|^2 ≈ constant for all k != 0

with |D11| = (p+1)/2 and |D12| = (p-1)/2.

This is related to but NOT identical to the standard Golay pair definition, which uses bipolar ({+/-1}) sequences of EQUAL length. Our sequences have UNEQUAL lengths and are {0,1}-valued.

---

## 2. Existence Results: What Is Known

### 2.1 Binary Aperiodic Golay Pairs

**Eliahou-Kervaire-Saffari Theorem (1990)**: If a binary Golay pair of length N exists, then N has NO prime factor congruent to 3 mod 4.

Known to exist: lengths N = 2^a * 10^b * 26^c (Turyn, 1974).

**Conjecture**: These are the ONLY lengths for which binary Golay pairs exist.

**CRITICAL NEGATIVE RESULT FOR US**: Since our p ≡ 3 mod 4 is prime, and (p+1)/2 or (p-1)/2 will often have prime factors ≡ 3 mod 4, binary aperiodic Golay pairs of the needed lengths generally do NOT exist. This rules out a direct application of aperiodic Golay pair theory.

### 2.2 Polyphase Golay Pairs

4-phase (entries in {1, -1, i, -i}) Golay pairs exist for arbitrary lengths (recent result, 2024). But these use complex entries, not {0,1} indicator functions.

### 2.3 Periodic Golay Pairs

More permissive than aperiodic. Known to exist for many lengths where aperiodic pairs do not.

**Necessary condition**: If a periodic Golay pair of length v exists, then a^2 + b^2 = 2v for some integers a, b (the row sums).

Known periodic Golay numbers (up to v ≤ 100): 2, 4, 8, 10, 16, 18, 20, 26, 32, 34, 36, 40, 50, 52, 58, 64, 68, 72, 74, 80, 82, 90(?), 98, 100.

**Exhaustive enumeration** completed for v ≤ 72 (2024 result).

Seven known periodic Golay numbers that are NOT aperiodic Golay numbers: 34, 50, 58, 68, 72, 74, 82.

**NOTE**: These are for {+/-1} sequences. Our {0,1} indicator functions live in a different setting.

### 2.4 Supplementary Difference Sets

Two subsets S1, S2 of Z_v form a 2-{v; k1, k2; lambda} supplementary difference set if the number of representations of each nonzero element as a difference within S1 plus the number within S2 equals lambda.

**Connection to Hadamard matrices**: If a 2-{v; k1, k2; lambda} SDS exists with appropriate parameters, it yields a Hadamard matrix of order 2v.

**Existence for prime power v**: Supplementary difference sets using cyclotomic classes (unions of cosets of 2m-th power residues) exist when v is an appropriate prime power. Specifically, m-(v, k, lambda) SDS exist with v = 2mf + 1 prime, f odd, k = mf, lambda = m(mf-1)/2.

### 2.5 Almost Supplementary Difference Sets (Armario-Flannery, 2020)

ASDS relax the exact condition. For odd m, certain ASDS in Z_m with amicable incidence matrices are equivalent to quaternary sequences of length m with optimal autocorrelation.

**Existence**: ASDS exist when 2m-1 is a prime power, or m ≡ 1 mod 4 is prime.

---

## 3. Connection to Our Problem

### 3.1 Translation Between Settings

Our (D11, D12) with |D11| = (p+1)/2, |D12| = (p-1)/2 over Z_p can be viewed through the SDS lens:

- The indicator functions 1_{D11}, 1_{D12} are {0,1}-valued sequences of length p
- Converting to bipolar: f_i(x) = 2*1_{D_i}(x) - 1 gives {+/-1} sequences
- The bipolar versions satisfy: |hat(f_1)(k)|^2 + |hat(f_2)(k)|^2 = |2*hat(1_{D11})(k) - p*delta(k)|^2 + similar
- This is NOT a clean periodic Golay pair because:
  (a) The sizes |D11| != |D12| (asymmetric)
  (b) The conversion from {0,1} to {+/-1} introduces DC offsets

### 3.2 What Spectral Complementarity Actually Means for SDSs

Our spectral complementarity |hat(1_{D11})(k)|^2 + |hat(1_{D12})(k)|^2 ≈ C for nonzero k translates to a condition on the combined difference function:

  delta_{D11}(d) + delta_{D12}(d) ≈ constant for all d != 0

where delta_S(d) = |{(x,y) in S x S : x - y = d}| is the difference function of S.

This is EXACTLY the supplementary difference set condition (with lambda = constant).

For our parameters (p; (p+1)/2, (p-1)/2; lambda), the expected lambda is:
  lambda = (|D11|^2 + |D12|^2 - |D11| - |D12|) / (p-1) = ((p+1)^2/4 + (p-1)^2/4 - p) / (p-1) = (p-1)/2

### 3.3 Why Standard Golay Theory Does NOT Directly Apply

1. **Length restriction**: Binary Golay pairs cannot have length p for p ≡ 3 mod 4 (EKS theorem). Our sequences live on Z_p.

2. **Unequal sizes**: Standard Golay pairs have equal-length sequences. Our |D11| = (p+1)/2 != |D12| = (p-1)/2.

3. **Alphabet mismatch**: Our sequences are {0,1}-valued (indicator functions), not {+/-1}-valued.

4. **Periodic vs aperiodic**: Our setting is inherently periodic (cyclic group Z_p).

5. **Approximate, not exact**: Our spectral complementarity is approximate (improving with p), not exact.

### 3.4 The Closest Match: Supplementary Difference Sets in Z_p

The closest framework is 2-{p; (p+1)/2, (p-1)/2; (p-1)/2} supplementary difference sets in Z_p.

**Key question**: Do such SDSs exist for all primes p ≡ 3 mod 4?

**What we know**:
- For p ≡ 3 mod 4, the quadratic residues and non-residues provide natural candidates
- Cyclotomic constructions yield SDSs for specific parameter families
- The Paley construction uses QRs to build conference matrices and thence Hadamard matrices for p ≡ 3 mod 4, but with DIFFERENT parameters than what we need
- The Dokovic-Kotsireas work constructs cyclic difference families with v ≡ 3 mod 4 prime, but with THREE blocks (not two)

---

## 4. Assessment for Our Proof

### 4.1 What the Golay Literature DOES Provide

1. **Framework**: The SDS framework precisely captures our spectral complementarity condition as a combinatorial object. This validates that our observed phenomenon has deep combinatorial structure.

2. **Cyclotomic methods**: The construction of SDSs via cyclotomic classes is the same algebraic toolkit used in our Paley-based constructions. This confirms the connection between quadratic residue structure and spectral complementarity.

3. **Asymptotic results**: The fact that polyphase Golay pairs exist for ALL lengths (in the 4-phase alphabet) suggests that approximate spectral complementarity in our {0,1} setting should also be achievable for all sufficiently large p.

4. **Near-complementary theory**: The existence of near-complementary pairs with bounded PMEPR (peak-to-mean envelope power ratio) provides a model for "approximate" spectral flatness.

### 4.2 What the Literature Does NOT Provide

1. **No direct existence theorem** for 2-{p; (p+1)/2, (p-1)/2; (p-1)/2} SDSs over Z_p for all p ≡ 3 mod 4.

2. **No quantitative bound** on spectral flatness achievable for {0,1}-valued sequences of these specific sizes over Z_p.

3. **No proof that approximate SDSs suffice** for our Ramsey graph application (we need the constraint violations to be bounded, not zero).

4. **No connection between SDS parameters and book graph constraints** (our A(d), B(d), C(d) thresholds).

### 4.3 Bottom Line

The Golay/SDS literature **validates our spectral complementarity observation** as a real phenomenon with deep algebraic roots, but it does **NOT provide a direct path to proving L6**. The gap between:

(a) "An approximate SDS exists" (which the literature partially addresses), and
(b) "The approximate SDS gives a valid Ramsey coloring" (which requires all constraint maxima to be bounded)

remains the core difficulty. The literature confirms this is a hard combinatorial optimization problem even in the structured (cyclotomic/Paley) setting.

---

## 5. Key References

1. **Golay, M.J.E.** (1961). Complementary series. IRE Trans. Inform. Theory, 7, 82-87.
   - Original definition of complementary pairs.

2. **Eliahou, S., Kervaire, M., Saffari, B.** (1990). A new restriction on the lengths of Golay complementary sequences. J. Combin. Theory Ser. A, 55, 49-59.
   - Binary Golay pairs exist only for lengths with no prime factor ≡ 3 mod 4.

3. **Turyn, R.J.** (1974). Hadamard matrices, Baumert-Hall units, four-symbol sequences, pulse compression, and surface wave encodings. J. Combin. Theory Ser. A, 16, 313-333.
   - Construction of Golay pairs of length 2^a * 10^b * 26^c.

4. **Fiedler, F., Jedwab, J., Parker, M.G.** (2008). A multi-dimensional approach to the construction and enumeration of Golay complementary sequences. J. Combin. Theory Ser. A, 115, 753-776.
   - Framework for all known binary Golay pairs via Boolean functions.

5. **Djokovic, D.Z., Kotsireas, I.S.** (2016). A class of cyclic (v; k1, k2, k3; lambda) difference families with v ≡ 3 (mod 4) a prime. Special Matrices, 4(1).
   - Three-block difference families using Paley-Todd sets.

6. **Armario, J.A., Flannery, D.L.** (2020). Almost supplementary difference sets and quaternary sequences with optimal autocorrelation. Cryptography and Communications.
   - ASDS existence when 2m-1 is prime power or m ≡ 1 mod 4 is prime.

7. **Crnkovic, D., Egan, R., Svob, A.** (2024). New results on periodic Golay pairs. arXiv:2408.15611.
   - Exhaustive enumeration of periodic Golay pairs for v ≤ 72.

8. **Shen, B. et al.** (2024). Golay complementary sequences of arbitrary length. arXiv:2401.15381.
   - 4-phase Golay pairs exist for all lengths.

9. **Wallis, J.** (1972). On supplementary difference sets. Aequationes Mathematicae.
   - Foundational work on SDS families.

10. **Lehmer, E.** (1974). A family of supplementary difference sets. Bull. Austral. Math. Soc., 11, 1-4.
    - Cyclotomic construction of SDS families.

---

## 6. Implications for Proof Strategy

### What to pursue:
- The SDS framework gives us a LANGUAGE for our spectral complementarity. The condition we need is essentially: "there exists a near-SDS in Z_p with parameters close to (p; (p+1)/2, (p-1)/2; (p-1)/2)."
- Cyclotomic constructions (union of cosets of higher-order residues) are the natural candidates. Our QR-based D11 is the simplest cyclotomic construction.
- The key metric is how far from exact SDS our pairs are -- measured by max_d |delta_{D11}(d) + delta_{D12}(d) - (p-1)/2|.

### What to avoid:
- Binary aperiodic Golay pair theory is a dead end (EKS theorem blocks it).
- Polyphase Golay pairs use complex entries, not {0,1} indicators.
- Periodic Golay pair enumeration is computational, not constructive.

### Most promising connection:
- The **cyclotomic SDS construction** (Lehmer, 1974; Wallis, 1972) directly builds SDSs from power residue classes in Z_p, which is exactly our algebraic setting. If we can show that the Paley/cyclotomic choice of D11 (quadratic residues) achieves near-SDS with deviation O(sqrt(p)) in the difference function, that would be the key spectral input for L6.
- The known deviation bound for a single difference set from QRs is O(sqrt(p)) by the Weil bound. The question is whether the COMBINED deviation delta_{D11}(d) + delta_{D12}(d) cancels better than individual terms -- our spectral complementarity data suggests it does.
