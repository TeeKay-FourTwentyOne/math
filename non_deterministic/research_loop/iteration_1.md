# Iteration 1: Semantic Taxonomy of Competing Tokens at Fragile Positions

**Date**: 2026-04-04
**Direction investigated**: What kind of tokens compete at fragile positions? Are hardware-induced divergences mostly cosmetic or meaning-changing?

---

## Summary of Key Findings

### Finding 1: 83.7% of fragile positions are meaning-changing, not cosmetic

At positions where the top-1 and top-2 logits are within 0.1 of each other (the threshold used to define "fragile"), the competing tokens were classified as:

| Category | Count | Percentage |
|----------|-------|------------|
| Meaning-changing (word-vs-word, digit-vs-digit, word-vs-digit) | 8,829 | 83.7% |
| Cosmetic (case variant, whitespace variant, punctuation swap) | 321 | 3.0% |
| Ambiguous (cross-type, newline-involved) | 1,395 | 13.2% |

**Implication**: Hardware-induced bit flips in floating-point arithmetic overwhelmingly cause *semantically meaningful* divergences, not just surface formatting differences. A different GPU executing the same model on the same prompt will produce text with different factual claims, different word choices, and different numerical outputs.

### Finding 2: Two dominant competition types -- word-vs-word (39.7%) and digit-vs-digit (35.5%)

The fragile positions split cleanly into two regimes:

- **word-vs-word** (4,191 / 39.7%): The model is uncertain about which word to use. Examples include `' aunt' vs ' uncle'`, `' American' vs ' South'`, `' was' vs ' is'`, `' philosopher' vs ' mathematic'`.
- **digit-vs-digit** (3,741 / 35.5%): The model is uncertain about which number to produce. Heavily concentrated in reasoning prompts (98.6% of all digit-vs-digit fragility).

Other types: word-vs-punctuation (5.7%), word-vs-digit (2.8%), punctuation-vs-punctuation (2.7%).

### Finding 3: Dramatic position-dependent semantic shift

The character of fragility shifts completely over the generation window:

| Position range | Total fragile | word-vs-word | digit-vs-digit |
|---------------|---------------|--------------|----------------|
| 0-31 | 4,497 | 63.1% | 13.1% |
| 32-63 | 2,000 | 43.8% | 24.2% |
| 64-127 | 1,674 | 20.3% | 58.4% |
| 128-255 | 2,374 | 11.5% | 71.2% |

The crossover point where digit fragility exceeds word fragility occurs at approximately **position 52-67**. This is NOT because the model starts computing later -- it's because:

1. **Early positions (0-31)**: The model is choosing discourse direction. "The" vs "I", "was" vs "is", "aunt" vs "uncle" -- these are narrative-level decisions.
2. **Late positions (64+)**: Reasoning prompts have entered repetitive arithmetic practice loops (see Finding 5), and the fragility is within the digit tokens of those loops.

The shift is category-dependent. For **non-reasoning** categories (factual, continuation, conversational, open_ended), word-vs-word dominates at ALL positions (60-80% throughout). The digit crossover is entirely driven by **reasoning** prompts, where digits dominate from position 0 onward (65.5% digit-vs-digit even in positions 0-15).

### Finding 4: The position-16 claim in STATE.md is INCORRECT

STATE.md claimed: "49.6% of ALL prompts are fragile at position 16, regardless of category."

**Actual finding**: Only 7.3% of prompts are fragile at position 16 (using the standard margin < 0.1 threshold). The 49.6% figure matches the fraction of prompts with margin < 1.0 at position 16 -- a 10x different threshold. This appears to be a threshold confusion in the original analysis.

The position-by-position fragility curve shows no special peak at position 16:
- Position 2: 11.8% (highest)
- Position 3: 11.3%
- Position 13: 12.1%
- Position 16: 7.3% (unremarkable)

The actual pattern is: fragility is highest in the first ~15 positions (when the model commits to a response strategy), then gradually declines. There is no cliff or step function at position 16.

**First fragile position distribution**: 47.1% of prompts have their first fragile position in positions 0-7. By position 31, 87.5% of prompts have encountered at least one fragile position.

### Finding 5: The 1-vs-3 digit confusion is a repetitive loop artifact

The most common competing digit pair is `'1' vs '3'` (861 occurrences), far exceeding adjacent-digit pairs like `'1' vs '2'` or `'2' vs '3'`. Investigation reveals this is NOT a genuine "confusion" between the numbers 1 and 3. Instead:

The model, when given arithmetic prompts like "What is 28 times 87?", generates repetitive practice-problem sequences:
```
"What is the product of -0.1 and -0.[1|3]..."
"What is the product of -0.1 and -0.[1|3]..."
(repeated every ~20 tokens for the entire 256-token window)
```

The 1-vs-3 competition happens at the same structural position in each repetition of this loop. All 861 instances come from just 93 unique reasoning prompts. The top prompt contributes 13 instances alone (one every ~20 tokens for 256 tokens of generation).

Similarly, the 2-vs-10 confusion (449 instances) comes from base-5 arithmetic loops like "In base 5, what is -10 + -[2|10]?".

**This means the raw count of 3,741 digit-vs-digit fragile positions is heavily inflated by repetitive loop patterns.** A more accurate measure would count unique prompt-level digit confusions (~93 prompts with repetitive digit loops, plus a smaller number with genuine single-instance digit fragility).

### Finding 6: Fine-grained word pair taxonomy

Among the 4,357 word-vs-word fragile positions, the breakdown is:

| Sub-category | Count | Percentage | Example |
|-------------|-------|------------|---------|
| Content word swap | 1,830 | 42.0% | 'aunt' vs 'uncle', 'ready' vs 'able' |
| Mixed function/content | 1,088 | 25.0% | 'the' vs 'Mexico', 'a' vs 'used' |
| Function word swap | 838 | 19.2% | 'The' vs 'I', 'of' vs 'is' |
| Article swap | 198 | 4.5% | 'a' vs 'the' |
| Prefix/substring | 166 | 3.8% | 'north' vs 'northern' |
| Pronoun swap | 102 | 2.3% | 'you' vs 'I' |
| Preposition swap | 87 | 2.0% | 'in' vs 'with' |
| Tense swap | 30 | 0.7% | 'was' vs 'is' |
| Conjunction swap | 18 | 0.4% | 'and' vs 'but' |

The content word swaps are particularly interesting because they represent genuine factual uncertainty:
- `' philosopher' vs ' mathematic'` (describing Pythagoras)
- `' wolf' vs ' human'` (discussing lifespans)
- `' fiction' vs ' novels'` (describing a writer's work)
- `' north' vs ' south'` (describing where a language is spoken)
- `' smallest' vs ' largest'` (answering a factual question -- outright contradiction!)

### Finding 7: Most common specific competing pairs

| Pair | Count | % of all fragile |
|------|-------|-----------------|
| '1' vs '3' | 861 | 8.2% |
| '1' vs '2' | 596 | 5.7% |
| '2' vs '10' | 449 | 4.3% |
| ' -' vs ' +' | 311 | 2.9% |
| ' 4' vs ' 5' | 294 | 2.8% |
| '3' vs '10' | 234 | 2.2% |
| ',' vs '.' | 184 | 1.7% |
| ' a' vs ' the' | 109 | 1.0% |
| '\n' vs ' The' | 92 | 0.9% |
| '\n' vs ' I' | 90 | 0.9% |
| 'I' vs 'The' | 67 | 0.6% |

The `' -' vs ' +'` pair (311 instances, 2.9%) is notable: this is the plus/minus operator confusion in arithmetic, which would flip the sign of an entire calculation.

### Finding 8: Cascade runs are short

Fragile positions are overwhelmingly isolated (not consecutive):
- Length 1 (isolated): 93.6%
- Length 2: 6.0%
- Length 3: 0.4%
- Length 4: 0.04%
- Maximum cascade: 4

Mean cascade length: 1.07. This means fragile positions are scattered, not clustered. A flip at one position does not systematically make the next position fragile (at least not when measured by logit margin alone -- the downstream effect via changed conditioning is a different question).

---

## Methodology

All analyses used the 1,999 logit trace files in `experiment_0/logit_traces/`, each containing the top-5 logit values and token IDs at every generation step (256 steps per prompt). "Fragile" = margin between rank-1 and rank-2 logits < 0.1. Token decoding used the Pythia-410M tokenizer (`EleutherAI/pythia-410m`).

Scripts:
- `analyze_competing_tokens.py`: Main 9-phase analysis (token pair classification, position patterns, cascade analysis, position-16 curve)
- `verify_position16.py`: Detailed position-16 claim verification
- `deep_semantic_analysis.py`: Fine-grained word pair taxonomy with contextual examples
- `digit_analysis.py`: Digit confusion matrix, distance analysis, 1-vs-3 investigation
- `position_shift_analysis.py`: Position-dependent word-to-digit shift with per-category breakdown
- `verify_1vs3.py`: Deep investigation of the 1-vs-3 confusion pattern

---

## Implications

1. **Hardware non-determinism in LLM inference is NOT a cosmetic problem.** 83.7% of fragile positions involve meaning-changing token alternatives. Running the same model on different hardware will produce factually different outputs.

2. **The semantic character of divergence depends on where in generation it occurs.** Early divergences change the narrative frame (which word, which perspective), while late divergences change numerical details.

3. **The position-16 "decision point" does not exist** at the margin < 0.1 threshold. The pattern is a smooth decline from high early fragility to low late fragility, with the first ~8 tokens being the most fragile.

4. **Digit-vs-digit fragility statistics are inflated by repetitive generation loops.** Future analyses should account for this by de-duplicating within-prompt repetitions or by counting unique prompt-level events.

5. **Antonym confusions exist** (e.g., `' smallest' vs ' largest'`, `' north' vs ' south'`, `' +' vs ' -'`). These are the most consequential type of fragile position -- a hardware bit flip could invert the meaning of a factual answer.
