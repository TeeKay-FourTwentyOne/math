# Double-Averaging Computation for E[Z]

**Date**: 2026-02-11
**Purpose**: Compute E[Z] = Pr[random (D11, D12) is valid] and compare to theoretical headroom.

## 1. Setup

For each prime p = 3 mod 4, we consider:
- **D11**: symmetric subset of {1,...,p-1} with |D11| = n = (p+1)/2
- **D12**: subset of {0,...,p-1} containing 0, with |D12| = (p-1)/2
  - Equivalently, D12 = {0} union S where S is a k-subset of {1,...,p-1}, k = (p-3)/2
- **Z(D11, D12)** = 1 if (D11, D12) satisfies all Ramsey constraints

The double average is:

  E[Z] = (1 / #D11) * (1 / #D12) * sum_{D11} sum_{D12} Z(D11, D12)
       = (1 / #D11) * (1 / #D12) * (total valid pairs)

## 2. Counting

| Quantity | Formula | p=11 | p=19 | p=23 | p=31 |
|----------|---------|------|------|------|------|
| n = (p+1)/2 | | 6 | 10 | 12 | 16 |
| k = (p-3)/2 | | 4 | 8 | 10 | 14 |
| #D11 = C(R, n/2), R=(p-1)/2 | | 10 | 126 | 462 | 6435 |
| #D12 = C(p-1, k) | | 210 | 43,758 | 646,646 | 145,422,675 |
| Total valid pairs | sum N(D11) | 100 | 162 | 4,356 | 364,560 |

## 3. Double-Averaging Results

| Quantity | p=11 | p=19 | p=23 | p=31 |
|----------|------|------|------|------|
| E[N(D11)] | 10.0000 | 1.2857 | 9.4286 | 56.6527 |
| E[N(D11)^2] | 200.00 | 23.14 | 1290.67 | 228668.72 |
| E[N^2]/E[N]^2 | 2.0000 | 14.0000 | 14.5185 | 71.2469 |
| E[Z] = Pr[valid] | 4.761905e-02 | 2.938238e-05 | 1.458073e-05 | 3.895725e-07 |
| log2(E[Z]) | -4.3923 | -15.0547 | -16.0656 | -21.2916 |
| log2(#D12) (headroom) | 7.7142 | 15.4173 | 19.3026 | 27.1157 |
| log2(E[N]) = log2(E[Z]) + log2(#D12) | 3.3219 | 0.3626 | 3.2370 | 5.8241 |
| Pr[N(D11) > 0] | 0.5000 | 0.0714 | 0.1190 | 0.0420 |
| Paley-Zygmund bound | 0.5000 | 0.0714 | 0.0689 | 0.0140 |

## 4. Interpretation

### 4.1 Headroom Analysis

The key quantity is **log2(E[Z]) - log2(#D12)**. This measures how many "bits" of headroom
remain after accounting for the constraint cost:

- log2(#D12) = log2(C(p-1, k)) is the entropy of the D12 choice (the "budget")
- log2(E[Z]) = log2(Pr[valid]) combines the budget with the constraint cost
- The difference = log2(E[N]) = log2(expected number of valid D12 per D11)

Note: log2(E[Z]) = log2(E[N]) - log2(#D12), so log2(E[N]) = log2(E[Z]) + log2(#D12).

**p = 11**: log2(E[N]) = -4.39 + 7.71 = 3.32. E[N] = 10.0, so on average each D11 has ~10 valid D12.

**p = 19**: log2(E[N]) = -15.05 + 15.42 = 0.36. E[N] = 1.29, so on average each D11 has ~1.3 valid D12.

**p = 23**: log2(E[N]) = -16.07 + 19.30 = 3.24. E[N] = 9.43, so on average each D11 has ~9.4 valid D12.

**p = 31**: log2(E[N]) = -21.29 + 27.12 = 5.82. E[N] = 56.65, so on average each D11 has ~56.7 valid D12.

### 4.2 Second Moment Ratio

The Paley-Zygmund inequality gives:

  Pr[N(D11) > 0] >= E[N]^2 / E[N^2] = 1 / (E[N^2]/E[N]^2)

| p | E[N^2]/E[N]^2 | PZ bound | Actual Pr[N>0] | Ratio actual/PZ |
|---|---------------|----------|----------------|-----------------|
| 11 | 2.00 | 0.5000 | 0.5000 | 1.00 |
| 19 | 14.00 | 0.0714 | 0.0714 | 1.00 |
| 23 | 14.52 | 0.0689 | 0.1190 | 1.73 |
| 31 | 71.25 | 0.0140 | 0.0420 | 2.99 |

The PZ bound is tight at p=11 and p=19 (ratio = 1.0), meaning the second moment
method is essentially optimal there. At p=23 and p=31, the actual probability exceeds
the PZ bound by a constant factor.

### 4.3 N-value Distributions

**p = 11**: N values = {0: 5, 20: 5}

**p = 19**: N values = {0: 117, 18: 9}

**p = 23**: N values = {0: 407, 22: 22, 44: 11, 110: 11, 198: 11}

**p = 31**: N values = {0: 6165, 62: 105, 124: 15, 186: 30, 558: 15, 620: 15, 1426: 15, 2790: 15, 3472: 15, 3720: 15, 4092: 15, 6696: 15}

### 4.4 Growth of E[N] and the Ratio

| p | E[N] | log2(E[N]) | E[N^2]/E[N]^2 | log2(ratio) | (p-1)/2 (bits budget) |
|---|------|------------|---------------|-------------|----------------------|
| 11 | 10.00 | 3.32 | 2.00 | 1.00 | 5 |
| 19 | 1.29 | 0.36 | 14.00 | 3.81 | 9 |
| 23 | 9.43 | 3.24 | 14.52 | 3.86 | 11 |
| 31 | 56.65 | 5.82 | 71.25 | 6.15 | 15 |

The budget grows as (p-1)/2 bits. For the second moment method to work, we need:
- E[N] >> 1 (equivalently log2(E[N]) >> 0)
- E[N^2]/E[N]^2 = O(poly(p))

Both conditions appear to be satisfied empirically.

## 5. Divisibility Structure of N(D11)

Nonzero N(D11) values exhibit divisibility patterns arising from symmetry:

- **Multiplicative orbit symmetry**: N(D11) = N(gD11) for all g in Z_p^*, so D11
  within the same orbit have the same N. Orbit size = (p-1)/2.
- **D12 orbits under Z_p^*/{+-1}**: each valid D12 orbit has size (p-1)/2.
- **Pair orbits**: each valid (D11,D12) pair orbit has size p-1 = 2 * (p-1)/2.

| p | Nonzero N values | GCD | p-1 | 2p |
|---|-----------------|-----|-----|-----|
| 11 | {20} | 20 | 10 | 22 |
| 19 | {18} | 18 | 18 | 38 |
| 23 | {22, 44, 110, 198} | 22 | 22 | 46 |
| 31 | {62, 124, 186, 558, 620, 1426, 2790, 3472, 3720, 4092, 6696} | 62 | 30 | 62 |

- p=11: GCD=20=2*10=2*(p-1). All N divisible by 2(p-1).
- p=19: GCD=18=p-1. All N divisible by p-1.
- p=23: GCD=22=p-1. All N divisible by p-1.
- p=31: GCD=62=2p. All N divisible by 2p.

The divisibility by p-1 at all primes is explained by the pair orbit structure
(each valid pair orbit contributes p-1 to the count). At p=31, the additional
factor of 2 (giving 2p = 62 instead of p-1 = 30) indicates extra structure.
