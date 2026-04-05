# Iteration 8: Antonym Mining — When Hardware Noise Inverts Meaning

**Direction**: Find cases where rank-1 and rank-2 tokens at fragile positions are semantic opposites — the most consequential type of hardware-induced divergence.

---

## Summary

**5.47% of all fragile positions involve some form of opposition**, but the picture has layers:
- The overwhelming majority (83%) are `'a' vs 'the'` — grammatically meaningful but not meaning-inverting
- True semantic inversions (man/woman, black/white, in/out, +/-) account for ~2% of fragile positions
- At realistic noise (ε=1e-2), ~4.6% of diverging prompts have an opposition at their first flip
- This means ~1 in 100 divergent outputs would have a genuinely meaning-inverting flip

---

## Opposition Classification (59,744 fragile positions, margin < 1.0)

| Type | Count | % of fragile |
|------|-------|-------------|
| Curated antonyms | 1,700 | 2.85% |
| Numeric opposition (+/-) | 1,354 | 2.27% |
| Boolean (true/false) | 201 | 0.34% |
| Prefix (un-/dis-/non-) | 12 | 0.02% |
| **Any opposition** | **3,267** | **5.47%** |

## The 'a' vs 'the' Dominance

The most common antonym pair is `'a' vs 'the'` with **1,404 occurrences** — 82.6% of all curated antonyms. This isn't a meaning inversion (nobody gets hurt if "a man" becomes "the man"), but it IS a real difference:
- "The man walked into the room" → definite, previously introduced character
- "A man walked into the room" → indefinite, new character

If we EXCLUDE 'a'/'the', the remaining genuinely meaning-inverting pairs:

| Pair | Count | Impact |
|------|-------|--------|
| man ↔ woman | 100 | Gender inversion |
| black ↔ white | 60 | Color inversion |
| large ↔ small | 41 | Size inversion |
| he ↔ she | 22 | Pronoun inversion |
| north ↔ south | 14 | Direction inversion |
| complex ↔ simple | 9 | Complexity inversion |
| new ↔ old | 6 | Temporal inversion |
| father ↔ mother | 5 | Kinship inversion |
| first ↔ last | 5 | Order inversion |
| big ↔ small | 4 | Size inversion |

Plus 1,354 numeric sign flips (` +` ↔ ` -`) in reasoning prompts.

## Most Spectacular Concrete Examples

### 1. Color inversion (margin = 0.001)
- **Prompt**: "Reports of missing pets in the neighborhood had been increasing..."
- **Position 48**: ` black` ↔ ` white` (margin 0.001087)
- **Context**: "...with a **black** and white coat..."
- A different GPU could describe the pet with reversed colors.

### 2. Direction inversion (margin = 0.003)
- **Prompt**: "The isolated hospital had been locked..."
- **Position 9**: ` out` ↔ ` in` (margin 0.002869)
- **Context**: "...door and get **out**"
- A different GPU: character gets IN instead of OUT.

### 3. Size inversion (margin = 0.008)
- **Prompt**: "It was twilight in a lakeside resort..."
- **Position 18**: ` small` ↔ ` large` (margin 0.007744)
- **Context**: "...carrying a **small** bag"
- A different GPU: character carries a LARGE bag.

### 4. Gender inversion (margin = 0.015)
- **Prompt**: "Three days after the incident at the schoolhouse..."
- **Position 75**: ` man` ↔ ` woman` (margin 0.014572)
- **Context**: '...said a **man** who had been...'
- A different GPU: the witness becomes a woman.

### 5. Arithmetic sign flip (margin = 0.0004)
- **Prompt**: "What is 65 times 50?"
- **Position 50**: ` -` ↔ ` +` (margin 0.000391)
- A different GPU would compute a different sign for a math result.

## Opposition Rates by Category

| Category | Fragile positions | Oppositions | Rate |
|----------|-------------------|-------------|------|
| reasoning | 18,771 | 1,759 | **9.37%** |
| factual | 9,534 | 545 | 5.72% |
| open_ended | 8,517 | 322 | 3.78% |
| continuation | 13,926 | 485 | 3.48% |
| conversational | 8,996 | 156 | 1.73% |

Reasoning dominates due to +/- sign flips. After factual, narrative categories have moderate rates driven by man/woman, a/the, and color pairs.

## Opposition Rate Among First Divergences

| ε | Diverging prompts | First-flip oppositions | Rate |
|---|-------------------|----------------------|------|
| 1e-3 | 87 | 3 | 3.4% |
| 1e-2 | 634 | 29 | **4.6%** |
| 1e-1 | 1,807 | 82 | 4.5% |

The rate is stable at ~4.5% regardless of noise level. At ε=1e-2:
- 634 prompts diverge
- ~29 have an opposition at the first flip point
- Of those 29: ~15 are +/- flips (reasoning), ~8 are a/the, ~6 are true meaning inversions

## Opposition Rate by Margin Size

| Margin range | Fragile positions | Oppositions | Rate |
|-------------|-------------------|-------------|------|
| [0.00, 0.01) | 1,062 | 51 | 4.80% |
| [0.01, 0.05) | 4,354 | 254 | 5.83% |
| [0.05, 0.10) | 5,129 | 214 | 4.17% |
| [0.10, 0.50) | 25,475 | 1,685 | 6.61% |
| [0.50, 1.00) | 23,724 | 1,063 | 4.48% |

The rate is remarkably flat — oppositions are NOT concentrated at tiny margins. They occur uniformly across the fragility spectrum.

## Practical Implications

1. **~1 in 20** hardware-divergent outputs has some form of opposition at the first flip point
2. **~1 in 100** has a genuinely meaning-inverting flip (excluding a/the)
3. The most dangerous flips are in narrative text: **man↔woman, black↔white, in↔out, small↔large**
4. In reasoning, **+/- sign flips** are common and would produce wrong answers with opposite sign
5. These rates are per-divergence, not per-prompt. At ε=1e-2, ~6 out of 1999 prompts (~0.3%) would have a truly meaning-inverting first token on a different GPU.

## Scripts
- `research_loop/iter8_antonym_mining.py`
