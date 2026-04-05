"""
Comprehensive analysis of LLM inference fragility experiment.
Pythia-410M, 1999 prompts, greedy decoding, 256 generation steps.
"""

import pandas as pd
import numpy as np
import torch
import json
import os
from collections import Counter, defaultdict
from pathlib import Path

BASE = Path("C:/Projects/math/non_deterministic")
LOGIT_DIR = BASE / "experiment_0" / "logit_traces"

# Load CSV
df = pd.read_csv(BASE / "prompts" / "fragility_scores.csv")

###############################################################################
# 1. FRAGILITY LANDSCAPE
###############################################################################
print("=" * 80)
print("SECTION 1: FRAGILITY LANDSCAPE")
print("=" * 80)

# Full distribution statistics
print("\n--- Full Distribution Statistics ---")
print(f"N prompts: {len(df)}")
print(f"Mean fragility: {df['fragility_score'].mean():.6f}")
print(f"Median fragility: {df['fragility_score'].median():.6f}")
print(f"Std dev: {df['fragility_score'].std():.6f}")
print(f"Min: {df['fragility_score'].min():.6f}")
print(f"Max: {df['fragility_score'].max():.6f}")
print(f"Skewness: {df['fragility_score'].skew():.4f}")
print(f"Kurtosis: {df['fragility_score'].kurtosis():.4f}")

# Percentiles
for p in [1, 5, 10, 25, 50, 75, 90, 95, 99]:
    print(f"  P{p}: {df['fragility_score'].quantile(p/100):.6f}")

# Histogram (text-based)
print("\n--- Histogram of Fragility Scores ---")
bins = np.arange(0, 0.70, 0.05)
hist, edges = np.histogram(df['fragility_score'], bins=bins)
for i in range(len(hist)):
    bar = "#" * (hist[i] // 5)
    print(f"  [{edges[i]:.2f}, {edges[i+1]:.2f}): {hist[i]:4d} {bar}")

# Check for clusters via mode detection
print("\n--- Cluster Analysis ---")
fine_bins = np.arange(0, 0.70, 0.01)
fine_hist, fine_edges = np.histogram(df['fragility_score'], bins=fine_bins)
# Find local maxima
local_max = []
for i in range(1, len(fine_hist)-1):
    if fine_hist[i] > fine_hist[i-1] and fine_hist[i] > fine_hist[i+1]:
        local_max.append((fine_edges[i], fine_edges[i+1], fine_hist[i]))
print("Local maxima in fine-grained histogram (0.01 bins):")
for lo, hi, cnt in sorted(local_max, key=lambda x: -x[2])[:10]:
    print(f"  [{lo:.2f}, {hi:.2f}): count={cnt}")

# By category
print("\n--- Fragility by Category ---")
cat_stats = df.groupby('category')['fragility_score'].agg(['mean', 'median', 'std', 'min', 'max'])
cat_stats = cat_stats.sort_values('mean', ascending=False)
print(cat_stats.to_string())

# Prompt length analysis
print("\n--- Fragility vs Prompt Length ---")
df['prompt_len'] = df['text'].str.len()
df['prompt_word_count'] = df['text'].str.split().str.len()

# Correlation
corr_len = df['fragility_score'].corr(df['prompt_len'])
corr_words = df['fragility_score'].corr(df['prompt_word_count'])
print(f"Correlation (fragility, char length): {corr_len:.4f}")
print(f"Correlation (fragility, word count): {corr_words:.4f}")

# Binned by length quartiles
df['len_quartile'] = pd.qcut(df['prompt_len'], 4, labels=['Q1(short)', 'Q2', 'Q3', 'Q4(long)'])
print("\nFragility by prompt length quartile:")
print(df.groupby('len_quartile')['fragility_score'].agg(['mean', 'median', 'count']).to_string())

# By category AND length
print("\nCorrelation by category:")
for cat in df['category'].unique():
    sub = df[df['category'] == cat]
    c = sub['fragility_score'].corr(sub['prompt_len'])
    print(f"  {cat}: corr={c:.4f} (n={len(sub)})")

# Most/least fragile prompts
print("\n--- Top 10 Most Fragile Prompts ---")
top10 = df.nlargest(10, 'fragility_score')
for _, row in top10.iterrows():
    print(f"  idx={row['prompt_idx']:4d}  frag={row['fragility_score']:.4f}  "
          f"cat={row['category']:<15s}  text={row['text'][:80]}")

print("\n--- Top 10 Least Fragile Prompts ---")
bot10 = df.nsmallest(10, 'fragility_score')
for _, row in bot10.iterrows():
    print(f"  idx={row['prompt_idx']:4d}  frag={row['fragility_score']:.4f}  "
          f"cat={row['category']:<15s}  text={row['text'][:80]}")


###############################################################################
# 2. MARGIN ANALYSIS
###############################################################################
print("\n" + "=" * 80)
print("SECTION 2: MARGIN ANALYSIS")
print("=" * 80)

# Collect all margins from all prompts
all_margins = []
fragile_margins = []  # margin < 1.0
margins_by_position = defaultdict(list)

print("Loading logit traces for margin analysis...")
for i in range(len(df)):
    pt_file = LOGIT_DIR / f"prompt_{i:04d}.pt"
    if not pt_file.exists():
        continue
    d = torch.load(pt_file, map_location='cpu', weights_only=True)
    vals = d['top5_values']  # (256, 5)
    margins = (vals[:, 0] - vals[:, 1]).numpy()
    all_margins.extend(margins.tolist())
    for pos, m in enumerate(margins):
        margins_by_position[pos].append(m)
        if m < 1.0:
            fragile_margins.append(m)

all_margins = np.array(all_margins)
fragile_margins = np.array(fragile_margins)

print(f"\nTotal token positions analyzed: {len(all_margins)}")
print(f"Total fragile positions (margin < 1.0): {len(fragile_margins)}")
print(f"Fraction fragile: {len(fragile_margins)/len(all_margins):.6f}")

print(f"\n--- All Margin Statistics ---")
print(f"Mean: {all_margins.mean():.4f}")
print(f"Median: {np.median(all_margins):.4f}")
print(f"Std: {all_margins.std():.4f}")
print(f"Min: {all_margins.min():.6f}")
print(f"Max: {all_margins.max():.4f}")

print(f"\n--- Fragile Margin Statistics (margin < 1.0) ---")
print(f"Mean: {fragile_margins.mean():.6f}")
print(f"Median: {np.median(fragile_margins):.6f}")
print(f"Std: {fragile_margins.std():.6f}")

# Extremely fragile
for threshold in [1.0, 0.5, 0.1, 0.05, 0.01, 0.001]:
    count = (all_margins < threshold).sum()
    frac = count / len(all_margins)
    print(f"  Margin < {threshold}: {count} positions ({frac:.6f} = {frac*100:.4f}%)")

# Distribution of fragile margins
print("\n--- Histogram of Fragile Margins (margin < 1.0) ---")
frag_bins = np.arange(0, 1.05, 0.1)
frag_hist, frag_edges = np.histogram(fragile_margins, bins=frag_bins)
for i in range(len(frag_hist)):
    bar = "#" * (frag_hist[i] // 50)
    print(f"  [{frag_edges[i]:.1f}, {frag_edges[i+1]:.1f}): {frag_hist[i]:5d} {bar}")

# Position-in-sequence analysis
print("\n--- Fragility by Position in Sequence ---")
# Divide 256 positions into 8 blocks of 32
for block_start in range(0, 256, 32):
    block_end = block_start + 32
    block_margins = []
    for pos in range(block_start, block_end):
        block_margins.extend(margins_by_position[pos])
    block_margins = np.array(block_margins)
    frac_fragile = (block_margins < 1.0).mean()
    mean_margin = block_margins.mean()
    print(f"  Positions [{block_start:3d}-{block_end-1:3d}]: "
          f"fragile_frac={frac_fragile:.4f}  mean_margin={mean_margin:.4f}")

# Fine-grained: per-position fragility rate
print("\n--- Per-Position Fragility Rate (sampled every 16 positions) ---")
for pos in range(0, 256, 16):
    m = np.array(margins_by_position[pos])
    frac = (m < 1.0).mean()
    mean_m = m.mean()
    print(f"  Position {pos:3d}: fragile_frac={frac:.4f}  mean_margin={mean_m:.4f}")


###############################################################################
# 3. TOKEN-LEVEL PATTERNS
###############################################################################
print("\n" + "=" * 80)
print("SECTION 3: TOKEN-LEVEL PATTERNS")
print("=" * 80)

# Load tokenizer
from transformers import AutoTokenizer
tokenizer = AutoTokenizer.from_pretrained("EleutherAI/pythia-410m")

# Sample 20 prompts spanning the fragility range
fragility_sorted = df.sort_values('fragility_score')
indices = np.linspace(0, len(df)-1, 20, dtype=int)
sample_prompts = fragility_sorted.iloc[indices]

# Categorize tokens at fragile positions
token_categories = Counter()
fragile_token_examples = []

# Function to categorize a token
def categorize_token(token_str):
    token_str_stripped = token_str.strip()
    if not token_str_stripped:
        return "whitespace"
    if token_str_stripped in '.,;:!?-()[]{}"\'/':
        return "punctuation"
    if token_str_stripped.isdigit() or all(c.isdigit() or c == '.' for c in token_str_stripped):
        return "digit"
    if token_str_stripped == '\n' or token_str_stripped == '\r\n':
        return "newline"
    # Function words (common English)
    func_words = {'the','a','an','of','in','to','for','is','was','are','were',
                  'be','been','being','have','has','had','do','does','did',
                  'will','would','could','should','may','might','shall','can',
                  'and','or','but','if','then','that','this','these','those',
                  'it','its','he','she','they','we','you','I','me','my','his',
                  'her','their','our','your','not','no','on','at','by','with',
                  'from','as','into','about','up','out','so','than','too'}
    if token_str_stripped.lower() in func_words:
        return "function_word"
    if token_str_stripped.isalpha():
        return "content_word"
    if token_str_stripped.startswith('\n'):
        return "newline"
    return "other"

all_fragile_tokens_all_prompts = Counter()
all_nonfragile_tokens_all_prompts = Counter()
fragile_token_cat_all = Counter()
nonfragile_token_cat_all = Counter()

# Analyze all prompts for token categorization (broader sample)
print("Analyzing token patterns across all prompts...")
for i in range(len(df)):
    pt_file = LOGIT_DIR / f"prompt_{i:04d}.pt"
    if not pt_file.exists():
        continue
    d = torch.load(pt_file, map_location='cpu', weights_only=True)
    vals = d['top5_values']
    gen_ids = d['generated_ids']
    input_len = d['input_ids'].shape[0]

    margins = (vals[:, 0] - vals[:, 1]).numpy()

    for pos in range(256):
        token_id = gen_ids[input_len + pos].item()
        token_str = tokenizer.decode([token_id])
        cat = categorize_token(token_str)

        if margins[pos] < 1.0:
            all_fragile_tokens_all_prompts[token_str] += 1
            fragile_token_cat_all[cat] += 1
        else:
            all_nonfragile_tokens_all_prompts[token_str] += 1
            nonfragile_token_cat_all[cat] += 1

total_fragile = sum(fragile_token_cat_all.values())
total_nonfragile = sum(nonfragile_token_cat_all.values())

print(f"\n--- Token Category Distribution ---")
print(f"{'Category':<20s} {'Fragile':>8s} {'%':>7s} {'Non-Fragile':>12s} {'%':>7s} {'Ratio':>7s}")
all_cats = set(fragile_token_cat_all.keys()) | set(nonfragile_token_cat_all.keys())
for cat in sorted(all_cats):
    fc = fragile_token_cat_all[cat]
    nc = nonfragile_token_cat_all[cat]
    fp = fc / total_fragile * 100 if total_fragile else 0
    np_ = nc / total_nonfragile * 100 if total_nonfragile else 0
    ratio = fp / np_ if np_ > 0 else float('inf')
    print(f"  {cat:<20s} {fc:8d} {fp:6.2f}% {nc:12d} {np_:6.2f}% {ratio:6.2f}x")

print(f"\n--- Top 30 Most Common Fragile Tokens ---")
for tok, cnt in all_fragile_tokens_all_prompts.most_common(30):
    nf_cnt = all_nonfragile_tokens_all_prompts.get(tok, 0)
    total_tok = cnt + nf_cnt
    fragile_rate = cnt / total_tok if total_tok > 0 else 0
    print(f"  {repr(tok):<25s} fragile={cnt:5d}  non-fragile={nf_cnt:5d}  fragile_rate={fragile_rate:.4f}")

# Show specific examples from sample prompts
print(f"\n--- Detailed Examples from 20 Sample Prompts ---")
for _, row in sample_prompts.iterrows():
    idx = int(row['prompt_idx'])
    pt_file = LOGIT_DIR / f"prompt_{idx:04d}.pt"
    if not pt_file.exists():
        continue
    d = torch.load(pt_file, map_location='cpu', weights_only=True)
    vals = d['top5_values']
    gen_ids = d['generated_ids']
    indices_t = d['top5_indices']
    input_len = d['input_ids'].shape[0]

    margins = (vals[:, 0] - vals[:, 1]).numpy()
    fragile_positions = np.where(margins < 1.0)[0]

    # Show the 3 most fragile positions
    if len(fragile_positions) == 0:
        continue
    most_fragile = fragile_positions[np.argsort(margins[fragile_positions])[:3]]

    print(f"\n  Prompt {idx} (frag={row['fragility_score']:.4f}, cat={row['category']}):")
    print(f"    Text: {row['text'][:80]}")
    for pos in most_fragile:
        chosen_id = gen_ids[input_len + pos].item()
        chosen_tok = tokenizer.decode([chosen_id])
        runner_up_id = indices_t[pos, 1].item()
        runner_up_tok = tokenizer.decode([runner_up_id])
        margin = margins[pos]
        # Context: 3 tokens before
        context_ids = gen_ids[max(0, input_len+pos-3):input_len+pos]
        context = tokenizer.decode(context_ids.tolist())
        print(f"    Pos {pos:3d}: margin={margin:.6f} "
              f"chosen={repr(chosen_tok)} vs runner_up={repr(runner_up_tok)} "
              f"context=...{repr(context)}")


###############################################################################
# 4. CATEGORY DEEP-DIVE
###############################################################################
print("\n" + "=" * 80)
print("SECTION 4: CATEGORY DEEP-DIVE")
print("=" * 80)

# Reasoning category analysis
print("\n--- Reasoning Category Deep-Dive ---")
reasoning = df[df['category'] == 'reasoning'].copy()
print(f"N reasoning prompts: {len(reasoning)}")
print(f"Mean fragility: {reasoning['fragility_score'].mean():.4f}")
print(f"Median fragility: {reasoning['fragility_score'].median():.4f}")

# Classify math problem types by keywords
def classify_reasoning(text):
    text_lower = text.lower()
    if any(w in text_lower for w in ['calculate', 'compute', 'what is', 'how much', 'how many', 'sum', 'product']):
        return 'arithmetic'
    if any(w in text_lower for w in ['solve', 'equation', 'variable', 'x =', 'find x']):
        return 'algebra'
    if any(w in text_lower for w in ['probability', 'chance', 'likely', 'odds']):
        return 'probability'
    if any(w in text_lower for w in ['sequence', 'pattern', 'series', 'next']):
        return 'sequence'
    if any(w in text_lower for w in ['area', 'perimeter', 'volume', 'triangle', 'circle', 'rectangle']):
        return 'geometry'
    if any(w in text_lower for w in ['logic', 'if', 'then', 'deduce', 'conclude', 'implies']):
        return 'logic'
    if any(w in text_lower for w in ['percent', 'ratio', 'proportion', 'fraction']):
        return 'ratio_percent'
    return 'other'

reasoning['subtype'] = reasoning['text'].apply(classify_reasoning)
print("\nReasoning subtypes:")
subtype_stats = reasoning.groupby('subtype')['fragility_score'].agg(['mean', 'median', 'count'])
subtype_stats = subtype_stats.sort_values('mean', ascending=False)
print(subtype_stats.to_string())

# Top/bottom 5 in reasoning
print("\nMost fragile reasoning prompts:")
for _, row in reasoning.nlargest(5, 'fragility_score').iterrows():
    print(f"  frag={row['fragility_score']:.4f}  {row['text'][:90]}")

print("\nLeast fragile reasoning prompts:")
for _, row in reasoning.nsmallest(5, 'fragility_score').iterrows():
    print(f"  frag={row['fragility_score']:.4f}  {row['text'][:90]}")

# Continuation category analysis
print("\n--- Continuation Category Deep-Dive ---")
continuation = df[df['category'] == 'continuation'].copy()
print(f"N continuation prompts: {len(continuation)}")

# Prompt length vs fragility
corr_cont = continuation['fragility_score'].corr(continuation['prompt_len'])
print(f"Correlation (fragility, char length) in continuation: {corr_cont:.4f}")
corr_cont_w = continuation['fragility_score'].corr(continuation['prompt_word_count'])
print(f"Correlation (fragility, word count) in continuation: {corr_cont_w:.4f}")

# Quartile breakdown
continuation['len_q'] = pd.qcut(continuation['prompt_len'], 4, labels=['Q1', 'Q2', 'Q3', 'Q4'])
print("\nContinuation fragility by length quartile:")
print(continuation.groupby('len_q')['fragility_score'].agg(['mean', 'median', 'count']).to_string())

# All categories comparison with more detail
print("\n--- All Categories: Detailed Statistics ---")
for cat in df['category'].unique():
    sub = df[df['category'] == cat]
    print(f"\n{cat} (n={len(sub)}):")
    print(f"  Mean frag: {sub['fragility_score'].mean():.4f}")
    print(f"  Median frag: {sub['fragility_score'].median():.4f}")
    print(f"  Std: {sub['fragility_score'].std():.4f}")
    print(f"  Mean min_margin: {sub['min_margin'].mean():.6f}")
    print(f"  Mean mean_margin: {sub['mean_margin'].mean():.4f}")
    print(f"  Mean prompt_len: {sub['prompt_len'].mean():.1f} chars")


###############################################################################
# 5. DIVERGENCE PREDICTION
###############################################################################
print("\n" + "=" * 80)
print("SECTION 5: DIVERGENCE PREDICTION")
print("=" * 80)

# For each prompt, find the minimum margin. If FP noise > min_margin, at least
# one token diverges. But the noise acts on logits, so we need the margin in
# logit space (which is what we have).
#
# Model: noise perturbation ~ epsilon applied to logits. A token flips if
# noise magnitude exceeds margin/2 (since both logits can be perturbed).
# More precisely, if |delta_logit1 - delta_logit2| > margin.
# For Gaussian noise with std epsilon, P(flip at one position) = P(|N(0, sqrt(2)*eps)| > margin)
# = 2 * Phi(-margin / (sqrt(2)*eps))
# For at least one flip in 256 positions: 1 - prod(1 - P(flip_i))

from scipy import stats

print("\n--- Method: Direct margin comparison ---")
print("If FP noise magnitude epsilon means any margin < epsilon will flip:")

# Collect per-prompt minimum margins
prompt_min_margins = df['min_margin'].values

for eps in [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.5, 1.0]:
    frac = (prompt_min_margins < eps).mean()
    count = (prompt_min_margins < eps).sum()
    print(f"  eps={eps:.0e}: {count:5d} prompts diverge ({frac*100:.4f}%)")

print("\n--- Method: Probabilistic (Gaussian FP noise on logits) ---")
print("Model: each logit gets independent N(0, eps) noise.")
print("P(flip at position) = P(|N(0, sqrt(2)*eps)| > margin) = 2*Phi(-margin/(sqrt(2)*eps))")
print("P(at least one flip in prompt) = 1 - prod(1-P_i) over 256 positions")

for eps in [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1]:
    # For each prompt, compute P(at least one flip)
    divergent_count = 0
    prompt_probs = []

    for i in range(len(df)):
        pt_file = LOGIT_DIR / f"prompt_{i:04d}.pt"
        if not pt_file.exists():
            continue
        d = torch.load(pt_file, map_location='cpu', weights_only=True)
        vals = d['top5_values']
        margins_arr = (vals[:, 0] - vals[:, 1]).numpy()

        # P(flip) at each position
        p_flip = 2 * stats.norm.cdf(-margins_arr / (np.sqrt(2) * eps))
        # P(at least one flip)
        p_no_flip = np.prod(1 - p_flip)
        p_any_flip = 1 - p_no_flip
        prompt_probs.append(p_any_flip)

    prompt_probs = np.array(prompt_probs)
    mean_prob = prompt_probs.mean()
    frac_likely = (prompt_probs > 0.5).mean()
    frac_nonzero = (prompt_probs > 1e-10).mean()

    print(f"\n  eps={eps:.0e}:")
    print(f"    Mean P(diverge): {mean_prob:.6e}")
    print(f"    Prompts with P(diverge)>0.5: {frac_likely*100:.4f}%")
    print(f"    Prompts with P(diverge)>1e-10: {frac_nonzero*100:.2f}%")
    print(f"    Expected divergent prompts (of 1999): {mean_prob*1999:.2f}")

# Also: expected number of flipped tokens across all prompts
print("\n--- Expected Total Token Flips ---")
for eps in [1e-7, 1e-5, 1e-3, 0.1]:
    total_flips = 0
    for i in range(len(df)):
        pt_file = LOGIT_DIR / f"prompt_{i:04d}.pt"
        if not pt_file.exists():
            continue
        d = torch.load(pt_file, map_location='cpu', weights_only=True)
        vals = d['top5_values']
        margins_arr = (vals[:, 0] - vals[:, 1]).numpy()
        p_flip = 2 * stats.norm.cdf(-margins_arr / (np.sqrt(2) * eps))
        total_flips += p_flip.sum()

    print(f"  eps={eps:.0e}: expected total flips = {total_flips:.4f} "
          f"(across {len(df)*256} positions)")

# Cumulative margin distribution
print("\n--- Cumulative Margin Distribution (All Positions) ---")
sorted_margins = np.sort(all_margins)
for threshold in [0.001, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
    frac_below = (sorted_margins < threshold).mean()
    print(f"  Fraction with margin < {threshold}: {frac_below:.6f}")

print("\n\nANALYSIS COMPLETE.")
