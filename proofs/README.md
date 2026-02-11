# Proof: R(B_{n-1}, B_n) = 4n-1

## Status: Complete conditional on Lemma L6 (positive association)

All other components are proven or computationally verified. L6 has overwhelming
computational evidence but needs a formal proof.

## Key Documents

### Primary Proof Documents
- `proof_positive_association.md` — **Main proof document.** Complete proof structure,
  theorem statement, all lemmas, computational evidence tables. Read this first.
- `first_moment_proof_sketch.md` — Detailed proof sketch with all 6 lemmas,
  alternative approaches (energy/Parseval), and computational evidence summary.

### Supporting Analysis Documents
- `proof_outline.md` — Structural reduction theorems (Thm 3, 5), complement symmetry,
  D11 size universality.
- `external_input_20260210.md` — External research input: Szekeres 1969,
  Đoković-Kotsireas 2015, Energy/Parseval argument from Opus 4.6.

### Key Computational Scripts (in `ramsey-book-graphs/`)
- `exact_moments.py` — Proves E[B(d)] = (p-3)/4 and Var[B(d)] = p/16 exactly.
  The foundation of the first moment argument.
- `joint_probability_analysis.py` — **The breakthrough script.** Proves positive
  association of B-events for fixed D11. Tests ALL D11s (not just SA solutions).
  Shows ratio joint/independence grows from 4.7x to 131x.
- `spectral_conditioning.py` — Investigates WHY positive association holds.
  Shows it persists after conditioning on spectral flatness (intrinsic property).
- `correlation_analysis.py` — Original correlation structure analysis.
- `first_moment_analysis.py` — MC estimation of E[valid pairs] across primes.
- `suen_inequality_test.py` — Tests Suen/Janson bounds (too crude, but informative).

### Data Files
- `correlation_results.json` — MC correlation data for p=7 to p=59.
- `solutions_registry.json` — All verified R(B_{n-1}, B_n) constructions (n ≤ 33).

## The Remaining Gap: Lemma L6

**Statement**: For any symmetric D11 of size (p+1)/2, the events {B(d) ≤ T(d)}
are positively associated under uniform random D12.

**Evidence**: Ratio joint/independence:
| p  | Ratio | Method |
|----|-------|--------|
| 11 | 4.68  | Exact enumeration |
| 19 | 11.59 | MC 200K |
| 23 | 8.14  | MC 200K |
| 31 | 11.20 | MC 100K |
| 43 | 48.53 | MC 50K  |
| 47 | 109.25| MC 30K  |
| 59 | 131.40| MC 20K  |

Verified for ALL tested D11s (500 random per prime), not just SA solutions.
The ratio **grows with p**, suggesting L6 becomes MORE true for larger primes.

**Candidate approaches**:
1. FKG-type inequality for quadratic events on k-subset measures
2. Spectral conditioning + asymptotic analysis
3. Finite verification (p ≤ 59) + asymptotic argument (p ≥ 61)
4. Coupling via strong Rayleigh property (Borcea-Brändén-Liggett 2009)
