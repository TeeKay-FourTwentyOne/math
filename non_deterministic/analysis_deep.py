"""
Deep-dive analysis: additional patterns and surprising findings.
"""

import pandas as pd
import numpy as np
import torch
import json
from collections import Counter, defaultdict
from pathlib import Path
from transformers import AutoTokenizer

BASE = Path("C:/Projects/math/non_deterministic")
LOGIT_DIR = BASE / "experiment_0" / "logit_traces"
tokenizer = AutoTokenizer.from_pretrained("EleutherAI/pythia-410m")

df = pd.read_csv(BASE / "prompts" / "fragility_scores.csv")

###############################################################################
# DEEP DIVE: Reasoning arithmetic sub-analysis
###############################################################################
print("=" * 80)
print("REASONING ARITHMETIC DEEP-DIVE")
print("=" * 80)

reasoning = df[df['category'] == 'reasoning']

# Classify arithmetic operations
def classify_arithmetic(text):
    text_lower = text.lower()
    if '+' in text or 'plus' in text_lower:
        return 'addition'
    if '-' in text and 'miles' not in text_lower:
        return 'subtraction'
    if 'times' in text_lower or '*' in text or 'multiply' in text_lower:
        return 'multiplication'
    if '/' in text or 'divide' in text_lower:
        return 'division'
    if '%' in text or 'percent' in text_lower:
        return 'percentage'
    if 'travel' in text_lower or 'car ' in text_lower:
        return 'word_problem'
    if 'sequence' in text_lower or 'pattern' in text_lower or 'next' in text_lower:
        return 'sequence'
    return 'other'

reasoning = reasoning.copy()
reasoning['op_type'] = reasoning['text'].apply(classify_arithmetic)
print("\nBy operation type:")
op_stats = reasoning.groupby('op_type')['fragility_score'].agg(['mean', 'median', 'std', 'count'])
op_stats = op_stats.sort_values('mean', ascending=False)
print(op_stats.to_string())

# For arithmetic: does number magnitude matter?
import re
def extract_numbers(text):
    return [int(x) for x in re.findall(r'\d+', text)]

reasoning['numbers'] = reasoning['text'].apply(extract_numbers)
reasoning['max_number'] = reasoning['numbers'].apply(lambda ns: max(ns) if ns else 0)
reasoning['num_digits'] = reasoning['max_number'].apply(lambda n: len(str(n)) if n > 0 else 0)

arith = reasoning[reasoning['op_type'].isin(['addition', 'subtraction', 'multiplication', 'percentage'])]
print(f"\n--- Number magnitude vs fragility (arithmetic only, n={len(arith)}) ---")
corr = arith['fragility_score'].corr(arith['max_number'])
corr_digits = arith['fragility_score'].corr(arith['num_digits'])
print(f"Correlation (fragility, max_number): {corr:.4f}")
print(f"Correlation (fragility, num_digits): {corr_digits:.4f}")

# By digit count
for nd in sorted(arith['num_digits'].unique()):
    sub = arith[arith['num_digits'] == nd]
    if len(sub) > 1:
        print(f"  {nd}-digit numbers: mean_frag={sub['fragility_score'].mean():.4f} "
              f"median={sub['fragility_score'].median():.4f} (n={len(sub)})")

###############################################################################
# DEEP DIVE: What makes word problems stable?
###############################################################################
print("\n" + "=" * 80)
print("WORD PROBLEMS vs RAW ARITHMETIC")
print("=" * 80)

word_probs = reasoning[reasoning['op_type'] == 'word_problem']
raw_arith = reasoning[reasoning['op_type'].isin(['addition', 'subtraction', 'multiplication', 'percentage'])]
print(f"Word problems: mean_frag={word_probs['fragility_score'].mean():.4f}, median={word_probs['fragility_score'].median():.4f}, n={len(word_probs)}")
print(f"Raw arithmetic: mean_frag={raw_arith['fragility_score'].mean():.4f}, median={raw_arith['fragility_score'].median():.4f}, n={len(raw_arith)}")

# Look at what the model generates for highest-fragility arithmetic
print("\n--- Sample outputs from most fragile arithmetic prompts ---")
for _, row in raw_arith.nlargest(5, 'fragility_score').iterrows():
    idx = int(row['prompt_idx'])
    pt_file = LOGIT_DIR / f"prompt_{idx:04d}.pt"
    if not pt_file.exists():
        continue
    d = torch.load(pt_file, map_location='cpu', weights_only=True)
    gen_ids = d['generated_ids']
    text = tokenizer.decode(gen_ids[:80].tolist())
    print(f"\n  Prompt {idx} (frag={row['fragility_score']:.4f}): {row['text']}")
    print(f"  Output (first 80 tokens): {text[:200]}")

# And least fragile word problems
print("\n--- Sample outputs from least fragile word problems ---")
for _, row in word_probs.nsmallest(5, 'fragility_score').iterrows():
    idx = int(row['prompt_idx'])
    pt_file = LOGIT_DIR / f"prompt_{idx:04d}.pt"
    if not pt_file.exists():
        continue
    d = torch.load(pt_file, map_location='cpu', weights_only=True)
    gen_ids = d['generated_ids']
    text = tokenizer.decode(gen_ids[:80].tolist())
    print(f"\n  Prompt {idx} (frag={row['fragility_score']:.4f}): {row['text']}")
    print(f"  Output (first 80 tokens): {text[:200]}")

###############################################################################
# DEEP DIVE: Fragile token pairs (chosen vs runner-up)
###############################################################################
print("\n" + "=" * 80)
print("FRAGILE TOKEN PAIRS: WHAT COMPETES?")
print("=" * 80)

pair_counter = Counter()
pair_ultra_counter = Counter()  # margin < 0.01

for i in range(len(df)):
    pt_file = LOGIT_DIR / f"prompt_{i:04d}.pt"
    if not pt_file.exists():
        continue
    d = torch.load(pt_file, map_location='cpu', weights_only=True)
    vals = d['top5_values']
    indices_t = d['top5_indices']
    gen_ids = d['generated_ids']
    input_len = d['input_ids'].shape[0]

    margins = (vals[:, 0] - vals[:, 1]).numpy()

    for pos in range(256):
        m = margins[pos]
        if m < 1.0:
            chosen = tokenizer.decode([indices_t[pos, 0].item()])
            runner = tokenizer.decode([indices_t[pos, 1].item()])
            pair_counter[(chosen, runner)] += 1
            if m < 0.01:
                pair_ultra_counter[(chosen, runner)] += 1

print("\n--- Top 30 Most Common Fragile Token Pairs (margin < 1.0) ---")
for (c, r), cnt in pair_counter.most_common(30):
    print(f"  {cnt:5d}x  {repr(c):<20s} vs {repr(r):<20s}")

print("\n--- Top 20 Ultra-Fragile Token Pairs (margin < 0.01) ---")
for (c, r), cnt in pair_ultra_counter.most_common(20):
    print(f"  {cnt:4d}x  {repr(c):<20s} vs {repr(r):<20s}")

###############################################################################
# DEEP DIVE: Entropy-like analysis of top-5
###############################################################################
print("\n" + "=" * 80)
print("TOP-5 LOGIT ENTROPY AT FRAGILE vs NON-FRAGILE POSITIONS")
print("=" * 80)

# Compute softmax entropy over top-5 for fragile vs non-fragile
fragile_entropies = []
nonfragile_entropies = []

for i in range(0, len(df), 10):  # sample every 10th
    pt_file = LOGIT_DIR / f"prompt_{i:04d}.pt"
    if not pt_file.exists():
        continue
    d = torch.load(pt_file, map_location='cpu', weights_only=True)
    vals = d['top5_values']
    margins = (vals[:, 0] - vals[:, 1]).numpy()

    # Softmax over top-5
    logits = vals.numpy()
    logits_shifted = logits - logits.max(axis=1, keepdims=True)
    probs = np.exp(logits_shifted) / np.exp(logits_shifted).sum(axis=1, keepdims=True)
    entropy = -(probs * np.log(probs + 1e-12)).sum(axis=1)

    for pos in range(256):
        if margins[pos] < 1.0:
            fragile_entropies.append(entropy[pos])
        else:
            nonfragile_entropies.append(entropy[pos])

fragile_entropies = np.array(fragile_entropies)
nonfragile_entropies = np.array(nonfragile_entropies)

print(f"Fragile positions: mean entropy={fragile_entropies.mean():.4f}, median={np.median(fragile_entropies):.4f}")
print(f"Non-fragile positions: mean entropy={nonfragile_entropies.mean():.4f}, median={np.median(nonfragile_entropies):.4f}")
print(f"Ratio: {fragile_entropies.mean() / nonfragile_entropies.mean():.2f}x")

###############################################################################
# DEEP DIVE: Continuation - narrative patterns
###############################################################################
print("\n" + "=" * 80)
print("CONTINUATION CATEGORY: NARRATIVE PATTERNS")
print("=" * 80)

continuation = df[df['category'] == 'continuation'].copy()

# Check for genre/style keywords
def classify_continuation(text):
    text_lower = text.lower()
    if any(w in text_lower for w in ['laboratory', 'scientist', 'experiment', 'research']):
        return 'science'
    if any(w in text_lower for w in ['detective', 'murder', 'crime', 'suspect', 'investigation']):
        return 'mystery'
    if any(w in text_lower for w in ['forest', 'mountain', 'river', 'ocean', 'garden', 'storm', 'rain']):
        return 'nature'
    if any(w in text_lower for w in ['castle', 'king', 'queen', 'sword', 'dragon', 'knight']):
        return 'fantasy'
    if any(w in text_lower for w in ['spaceship', 'planet', 'galaxy', 'alien', 'space']):
        return 'scifi'
    if any(w in text_lower for w in ['village', 'town', 'city', 'street', 'shop', 'market']):
        return 'urban'
    return 'general'

continuation['genre'] = continuation['text'].apply(classify_continuation)
print("By genre:")
genre_stats = continuation.groupby('genre')['fragility_score'].agg(['mean', 'median', 'count'])
genre_stats = genre_stats.sort_values('mean', ascending=False)
print(genre_stats.to_string())

###############################################################################
# DEEP DIVE: Position 16 spike investigation
###############################################################################
print("\n" + "=" * 80)
print("POSITION 16 SPIKE INVESTIGATION")
print("=" * 80)

# Position 16 had the highest fragility rate (0.4957). Why?
pos16_fragile_tokens = Counter()
pos16_context = []

for i in range(len(df)):
    pt_file = LOGIT_DIR / f"prompt_{i:04d}.pt"
    if not pt_file.exists():
        continue
    d = torch.load(pt_file, map_location='cpu', weights_only=True)
    vals = d['top5_values']
    gen_ids = d['generated_ids']
    indices_t = d['top5_indices']
    input_len = d['input_ids'].shape[0]

    margin_16 = (vals[16, 0] - vals[16, 1]).item()
    if margin_16 < 1.0:
        chosen = tokenizer.decode([indices_t[16, 0].item()])
        runner = tokenizer.decode([indices_t[16, 1].item()])
        pos16_fragile_tokens[(chosen, runner)] += 1
        if len(pos16_context) < 10:
            context = tokenizer.decode(gen_ids[max(0, input_len+13):input_len+16].tolist())
            full = tokenizer.decode(gen_ids[max(0, input_len+13):input_len+19].tolist())
            pos16_context.append((i, margin_16, chosen, runner, context, full, df.iloc[i]['category']))

print("Top token pairs at position 16 (when fragile):")
for (c, r), cnt in pos16_fragile_tokens.most_common(15):
    print(f"  {cnt:4d}x  {repr(c):<20s} vs {repr(r):<20s}")

print("\nSample contexts at position 16:")
for idx, m, c, r, ctx, full, cat in pos16_context:
    print(f"  Prompt {idx} ({cat}): margin={m:.4f}  ...{repr(ctx)} -> [{repr(c)} vs {repr(r)}] -> {repr(full)}")

###############################################################################
# DEEP DIVE: Cascading fragility (consecutive fragile positions)
###############################################################################
print("\n" + "=" * 80)
print("CASCADING FRAGILITY: CONSECUTIVE FRAGILE POSITIONS")
print("=" * 80)

max_run_lengths = []
run_length_counts = Counter()

for i in range(len(df)):
    pt_file = LOGIT_DIR / f"prompt_{i:04d}.pt"
    if not pt_file.exists():
        continue
    d = torch.load(pt_file, map_location='cpu', weights_only=True)
    vals = d['top5_values']
    margins = (vals[:, 0] - vals[:, 1]).numpy()
    fragile = margins < 1.0

    # Find consecutive runs
    max_run = 0
    current_run = 0
    for pos in range(256):
        if fragile[pos]:
            current_run += 1
            max_run = max(max_run, current_run)
        else:
            if current_run > 0:
                run_length_counts[current_run] += 1
            current_run = 0
    if current_run > 0:
        run_length_counts[current_run] += 1
    max_run_lengths.append(max_run)

max_run_lengths = np.array(max_run_lengths)
print(f"Max consecutive fragile run per prompt:")
print(f"  Mean: {max_run_lengths.mean():.2f}")
print(f"  Median: {np.median(max_run_lengths):.2f}")
print(f"  Max: {max_run_lengths.max()}")
print(f"  % prompts with run >= 5: {(max_run_lengths >= 5).mean()*100:.1f}%")
print(f"  % prompts with run >= 10: {(max_run_lengths >= 10).mean()*100:.1f}%")
print(f"  % prompts with run >= 20: {(max_run_lengths >= 20).mean()*100:.1f}%")

print(f"\nRun length distribution:")
for rl in sorted(run_length_counts.keys()):
    if rl <= 20 or run_length_counts[rl] >= 5:
        print(f"  Run length {rl:3d}: {run_length_counts[rl]:5d} occurrences")

# Correlation between fragility and max run
corr_run = np.corrcoef(df['fragility_score'].values, max_run_lengths)[0, 1]
print(f"\nCorrelation (fragility_score, max_consecutive_run): {corr_run:.4f}")

###############################################################################
# SUMMARY: Key findings
###############################################################################
print("\n" + "=" * 80)
print("SUMMARY OF KEY FINDINGS")
print("=" * 80)

print("""
1. FRAGILITY LANDSCAPE:
   - Right-skewed distribution (skew=1.90), mode at 0.05-0.06, long tail to 0.625
   - No discrete clusters; continuous distribution with single mode
   - Weak negative correlation with prompt length (r=-0.15)
   - STRONG negative correlation within reasoning (r=-0.82): shorter arithmetic prompts
     are far more fragile than longer word problems

2. MARGIN ANALYSIS:
   - 11.7% of all token positions are fragile (margin < 1.0)
   - 2.1% have margin < 0.1 (extremely fragile)
   - 0.2% have margin < 0.01 (hair-trigger fragile)
   - 93 positions across all prompts have margin < 0.001
   - DRAMATIC position dependence: 46.6% fragile at positions 0-31, dropping
     to 2.9% at positions 192-223. Position 16 peaks at 49.6%!
   - The model "stabilizes" its predictions as it generates more context

3. TOKEN PATTERNS:
   - Digits are 5.1x OVERREPRESENTED at fragile positions (15% vs 3% baseline)
   - Whitespace is 0.38x underrepresented (model is confident about newlines)
   - ' +' has 97.4% fragile rate; '10' has 94.8%; '1' has 71.7%
   - Content words and function words are slightly underrepresented at fragile positions

4. CATEGORY DEEP-DIVE:
   - Reasoning (mean=0.183) >> continuation (0.136) > factual (0.093) ~
     conversational (0.088) ~ open-ended (0.083)
   - Within reasoning: percentage/addition of 4-digit numbers are most fragile
   - Word problems (distance/speed/time) are highly stable (mean ~0.01)
   - Model gets stuck in arithmetic loops, generating digits with near-equal logits

5. DIVERGENCE PREDICTION:
   - eps=1e-7: effectively zero divergence
   - eps=1e-5: ~1 prompt diverges (0.05%)
   - eps=1e-3: ~87 prompts (4.4%) show deterministic flip; ~102 expected probabilistically
   - eps=1e-2: ~634 prompts (31.7%) deterministic; ~692 expected probabilistically
   - eps=0.1:  ~1808 prompts (90.4%) diverge
""")

print("DEEP ANALYSIS COMPLETE.")
