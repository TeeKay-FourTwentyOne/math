"""
Iteration 1: Semantic Analysis of Competing Tokens at Fragile Positions

Investigates what kinds of tokens compete at fragile positions (where margin < 0.1
between rank-1 and rank-2 logits). Key question: are hardware-induced divergences
mostly cosmetic (synonym-like swaps) or do they cause meaningful semantic shifts?

Analysis plan:
1. Collect all fragile positions across all 1999 prompts
2. Classify the token pair types (digit-digit, punctuation-punctuation, word-word, etc.)
3. Measure "semantic distance" between competing tokens using:
   a. Token type classification (same-type vs cross-type)
   b. Character-level similarity (edit distance)
   c. Contextual analysis: does the competing token change the meaning?
4. Check if certain token pair types cluster at specific generation positions
5. Look at cascade effects: when token 1 vs token 2 is chosen, what happens next?
"""

import torch
import numpy as np
import json
import os
from collections import Counter, defaultdict
from pathlib import Path

# We'll use Pythia tokenizer concepts - the model is Pythia-410M based on the 1024-dim hidden size
# and 24-layer architecture.  Let's load a tokenizer to decode tokens.
try:
    from transformers import AutoTokenizer
    tokenizer = AutoTokenizer.from_pretrained("EleutherAI/pythia-410m")
    HAS_TOKENIZER = True
    print("Loaded Pythia-410M tokenizer")
except Exception as e:
    print(f"Cannot load tokenizer: {e}")
    HAS_TOKENIZER = False

LOGIT_DIR = "C:/Projects/math/non_deterministic/experiment_0/logit_traces"
FRAGILITY_CSV = "C:/Projects/math/non_deterministic/prompts/fragility_scores.csv"
PROMPTS_FILE = "C:/Projects/math/non_deterministic/prompts/candidates_2000.jsonl"
MARGIN_THRESHOLD = 0.1  # Same threshold used in the scoring

def decode_token(token_id):
    """Decode a single token ID to its string representation."""
    if HAS_TOKENIZER:
        return tokenizer.decode([token_id])
    return f"<{token_id}>"

def classify_token(text):
    """Classify a token into a category."""
    text_stripped = text.strip()
    if not text_stripped:
        return "whitespace"
    if text_stripped.isdigit():
        return "digit"
    if all(c in '.,;:!?-\'"()[]{}' for c in text_stripped):
        return "punctuation"
    if text_stripped.startswith('\n'):
        return "newline"
    if text_stripped.isalpha():
        if text_stripped[0].isupper():
            return "word_capitalized"
        return "word_lower"
    if any(c.isdigit() for c in text_stripped) and any(c.isalpha() for c in text_stripped):
        return "alphanumeric"
    return "other"

def semantic_similarity_category(tok1_text, tok2_text):
    """Categorize the semantic relationship between two competing tokens."""
    t1 = tok1_text.strip()
    t2 = tok2_text.strip()

    # Same token (shouldn't happen, but just in case)
    if t1 == t2:
        return "identical"

    # Case variation
    if t1.lower() == t2.lower():
        return "case_variant"

    # Whitespace variants
    if tok1_text.replace(' ', '') == tok2_text.replace(' ', ''):
        return "whitespace_variant"

    c1 = classify_token(tok1_text)
    c2 = classify_token(tok2_text)

    # Same type comparisons
    if c1 == "digit" and c2 == "digit":
        return "digit_vs_digit"
    if c1 == "punctuation" and c2 == "punctuation":
        return "punctuation_vs_punctuation"
    if "word" in c1 and "word" in c2:
        # Check if one is a prefix/suffix of the other
        if t1.lower().startswith(t2.lower()) or t2.lower().startswith(t1.lower()):
            return "word_prefix"
        return "word_vs_word"
    if c1 == "newline" or c2 == "newline":
        return "newline_involved"

    # Cross-type
    if (c1 == "digit" or c2 == "digit") and (c1 == "punctuation" or c2 == "punctuation"):
        return "digit_vs_punctuation"
    if ("word" in c1 or "word" in c2) and (c1 == "digit" or c2 == "digit"):
        return "word_vs_digit"
    if ("word" in c1 or "word" in c2) and (c1 == "punctuation" or c2 == "punctuation"):
        return "word_vs_punctuation"

    return f"cross_type:{c1}_vs_{c2}"

print("=" * 70)
print("PHASE 1: Collecting all fragile positions across 1999 prompts")
print("=" * 70)

all_fragile_pairs = []  # (prompt_idx, position, tok1_id, tok2_id, margin, tok1_logit, tok2_logit)
total_positions = 0
total_fragile = 0
per_prompt_stats = []

# Load prompt metadata
prompts = {}
with open(PROMPTS_FILE) as f:
    for line in f:
        d = json.loads(line)
        prompts[d['idx']] = d

for i in range(1999):
    path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
    if not os.path.exists(path):
        continue
    data = torch.load(path, map_location='cpu', weights_only=False)
    top5_vals = data['top5_values']  # (256, 5)
    top5_idx = data['top5_indices']  # (256, 5)

    margins = top5_vals[:, 0] - top5_vals[:, 1]
    n_pos = len(margins)
    total_positions += n_pos

    fragile_mask = margins < MARGIN_THRESHOLD
    n_fragile = fragile_mask.sum().item()
    total_fragile += n_fragile

    per_prompt_stats.append({
        'idx': i,
        'category': prompts[i]['category'],
        'n_fragile': n_fragile,
        'n_total': n_pos,
        'fragility': n_fragile / n_pos
    })

    for pos in range(n_pos):
        if fragile_mask[pos]:
            all_fragile_pairs.append((
                i, pos,
                top5_idx[pos, 0].item(), top5_idx[pos, 1].item(),
                margins[pos].item(),
                top5_vals[pos, 0].item(), top5_vals[pos, 1].item()
            ))

print(f"Total positions scanned: {total_positions:,}")
print(f"Total fragile positions: {total_fragile:,} ({100*total_fragile/total_positions:.1f}%)")
print(f"Prompts processed: {len(per_prompt_stats)}")

print("\n" + "=" * 70)
print("PHASE 2: Token pair classification")
print("=" * 70)

pair_categories = Counter()
pair_examples = defaultdict(list)  # category -> [(prompt, pos, tok1, tok2, margin), ...]
tok1_types = Counter()
tok2_types = Counter()
type_pairs = Counter()

for prompt_idx, pos, tok1_id, tok2_id, margin, logit1, logit2 in all_fragile_pairs:
    tok1_text = decode_token(tok1_id)
    tok2_text = decode_token(tok2_id)

    cat = semantic_similarity_category(tok1_text, tok2_text)
    pair_categories[cat] += 1

    if len(pair_examples[cat]) < 10:
        pair_examples[cat].append((prompt_idx, pos, tok1_text, tok2_text, margin))

    t1 = classify_token(tok1_text)
    t2 = classify_token(tok2_text)
    tok1_types[t1] += 1
    tok2_types[t2] += 1
    type_pairs[f"{t1} vs {t2}"] += 1

print("\nSemantic relationship between competing tokens:")
print("-" * 50)
for cat, count in pair_categories.most_common(20):
    pct = 100 * count / len(all_fragile_pairs)
    print(f"  {cat:35s}: {count:6d} ({pct:5.1f}%)")
    # Show 3 examples
    for ex in pair_examples[cat][:3]:
        prompt_idx, pos, t1, t2, margin = ex
        print(f"      prompt {prompt_idx:4d} pos {pos:3d}: {repr(t1):20s} vs {repr(t2):20s} (margin={margin:.4f})")

print("\nToken type for rank-1 at fragile positions:")
for t, c in tok1_types.most_common(10):
    print(f"  {t:20s}: {c:6d} ({100*c/len(all_fragile_pairs):5.1f}%)")

print("\nToken type for rank-2 at fragile positions:")
for t, c in tok2_types.most_common(10):
    print(f"  {t:20s}: {c:6d} ({100*c/len(all_fragile_pairs):5.1f}%)")

print("\nTop type-pair combinations:")
for tp, c in type_pairs.most_common(15):
    print(f"  {tp:40s}: {c:6d} ({100*c/len(all_fragile_pairs):5.1f}%)")

print("\n" + "=" * 70)
print("PHASE 3: Position-dependent patterns")
print("=" * 70)

# Group fragile pairs by generation position
pos_to_cats = defaultdict(lambda: Counter())
for prompt_idx, pos, tok1_id, tok2_id, margin, logit1, logit2 in all_fragile_pairs:
    tok1_text = decode_token(tok1_id)
    tok2_text = decode_token(tok2_id)
    cat = semantic_similarity_category(tok1_text, tok2_text)
    pos_to_cats[pos][cat] += 1

# Bin positions
pos_bins = [(0, 15), (16, 31), (32, 63), (64, 127), (128, 255)]
print("\nSemantic category distribution by position bin:")
for lo, hi in pos_bins:
    merged = Counter()
    total_in_bin = 0
    for p in range(lo, hi + 1):
        merged += pos_to_cats[p]
        total_in_bin += sum(pos_to_cats[p].values())
    if total_in_bin == 0:
        continue
    print(f"\n  Positions {lo}-{hi} ({total_in_bin} fragile positions):")
    for cat, count in merged.most_common(5):
        print(f"    {cat:35s}: {count:5d} ({100*count/total_in_bin:5.1f}%)")

print("\n" + "=" * 70)
print("PHASE 4: Per-category breakdown of competing token types")
print("=" * 70)

# Load categories
cat_to_pairs = defaultdict(list)
for prompt_idx, pos, tok1_id, tok2_id, margin, logit1, logit2 in all_fragile_pairs:
    cat = prompts[prompt_idx]['category']
    tok1_text = decode_token(tok1_id)
    tok2_text = decode_token(tok2_id)
    sem_cat = semantic_similarity_category(tok1_text, tok2_text)
    cat_to_pairs[cat].append(sem_cat)

for cat in sorted(cat_to_pairs.keys()):
    sem_counts = Counter(cat_to_pairs[cat])
    total = len(cat_to_pairs[cat])
    print(f"\n{cat} ({total} fragile positions):")
    for sc, count in sem_counts.most_common(5):
        print(f"  {sc:35s}: {count:5d} ({100*count/total:5.1f}%)")

print("\n" + "=" * 70)
print("PHASE 5: Cascade analysis - what happens after a fragile position?")
print("=" * 70)

# For each fragile position, check if the NEXT position is also fragile
# and if the semantic category changes
consecutive_fragile = 0
single_fragile = 0
cascade_lengths = []

# Build per-prompt fragile position sets
prompt_fragile_positions = defaultdict(set)
for prompt_idx, pos, tok1_id, tok2_id, margin, logit1, logit2 in all_fragile_pairs:
    prompt_fragile_positions[prompt_idx].add(pos)

for prompt_idx, positions in prompt_fragile_positions.items():
    sorted_pos = sorted(positions)
    # Find consecutive runs
    i = 0
    while i < len(sorted_pos):
        run_start = sorted_pos[i]
        run_len = 1
        while i + run_len < len(sorted_pos) and sorted_pos[i + run_len] == run_start + run_len:
            run_len += 1
        cascade_lengths.append(run_len)
        if run_len > 1:
            consecutive_fragile += run_len
        else:
            single_fragile += 1
        i += run_len

cascade_arr = np.array(cascade_lengths)
print(f"Total fragile runs: {len(cascade_lengths)}")
print(f"Run length distribution:")
for length in range(1, min(11, cascade_arr.max() + 1)):
    count = (cascade_arr == length).sum()
    print(f"  Length {length:2d}: {count:5d} ({100*count/len(cascade_lengths):5.1f}%)")
if cascade_arr.max() >= 11:
    count = (cascade_arr >= 11).sum()
    print(f"  Length 11+: {count:5d} ({100*count/len(cascade_lengths):5.1f}%)")

print(f"\nMax cascade length: {cascade_arr.max()}")
print(f"Mean cascade length: {cascade_arr.mean():.2f}")
print(f"Median cascade length: {np.median(cascade_arr):.1f}")

# Now look at what types of tokens appear in long cascades
print("\n\nLong cascade examples (length >= 5):")
long_cascade_count = 0
for prompt_idx, positions in prompt_fragile_positions.items():
    sorted_pos = sorted(positions)
    i = 0
    while i < len(sorted_pos):
        run_start = sorted_pos[i]
        run_len = 1
        while i + run_len < len(sorted_pos) and sorted_pos[i + run_len] == run_start + run_len:
            run_len += 1
        if run_len >= 5 and long_cascade_count < 5:
            long_cascade_count += 1
            # Load the prompt data
            path = os.path.join(LOGIT_DIR, f"prompt_{prompt_idx:04d}.pt")
            data = torch.load(path, map_location='cpu', weights_only=False)
            top5_idx = data['top5_indices']
            top5_vals = data['top5_values']

            print(f"\n  Prompt {prompt_idx} ({prompts[prompt_idx]['category']}): "
                  f"'{prompts[prompt_idx]['text'][:60]}...'")
            print(f"  Cascade at positions {run_start}-{run_start + run_len - 1} (length {run_len}):")
            for p in range(run_start, min(run_start + run_len, run_start + 8)):
                t1 = decode_token(top5_idx[p, 0].item())
                t2 = decode_token(top5_idx[p, 1].item())
                margin = top5_vals[p, 0].item() - top5_vals[p, 1].item()
                print(f"    pos {p:3d}: {repr(t1):15s} vs {repr(t2):15s} (margin={margin:.4f})")
        i += run_len

print("\n" + "=" * 70)
print("PHASE 6: Meaning-changing vs cosmetic divergences")
print("=" * 70)

# Classify each fragile position as potentially meaning-changing or cosmetic
meaning_changing = 0
cosmetic = 0
ambiguous = 0

COSMETIC_CATS = {'case_variant', 'whitespace_variant', 'identical', 'punctuation_vs_punctuation'}
MEANING_CATS = {'word_vs_word', 'digit_vs_digit', 'word_vs_digit', 'word_vs_punctuation'}

for prompt_idx, pos, tok1_id, tok2_id, margin, logit1, logit2 in all_fragile_pairs:
    tok1_text = decode_token(tok1_id)
    tok2_text = decode_token(tok2_id)
    cat = semantic_similarity_category(tok1_text, tok2_text)

    if cat in COSMETIC_CATS:
        cosmetic += 1
    elif cat in MEANING_CATS:
        meaning_changing += 1
    else:
        ambiguous += 1

total = len(all_fragile_pairs)
print(f"Meaning-changing: {meaning_changing:6d} ({100*meaning_changing/total:.1f}%)")
print(f"Cosmetic:         {cosmetic:6d} ({100*cosmetic/total:.1f}%)")
print(f"Ambiguous:        {ambiguous:6d} ({100*ambiguous/total:.1f}%)")

print("\n" + "=" * 70)
print("PHASE 7: The specific tokens most often involved in fragile pairs")
print("=" * 70)

# What specific token IDs appear most at rank-1 and rank-2 in fragile positions?
rank1_tokens = Counter()
rank2_tokens = Counter()
pair_counter = Counter()

for prompt_idx, pos, tok1_id, tok2_id, margin, logit1, logit2 in all_fragile_pairs:
    rank1_tokens[tok1_id] += 1
    rank2_tokens[tok2_id] += 1
    pair_key = (min(tok1_id, tok2_id), max(tok1_id, tok2_id))
    pair_counter[pair_key] += 1

print("\nMost common rank-1 tokens at fragile positions:")
for tok_id, count in rank1_tokens.most_common(15):
    tok_text = decode_token(tok_id)
    print(f"  {tok_id:6d} {repr(tok_text):20s}: {count:5d} ({100*count/total:.1f}%)")

print("\nMost common rank-2 tokens at fragile positions:")
for tok_id, count in rank2_tokens.most_common(15):
    tok_text = decode_token(tok_id)
    print(f"  {tok_id:6d} {repr(tok_text):20s}: {count:5d} ({100*count/total:.1f}%)")

print("\nMost common competing PAIRS at fragile positions:")
for (t1, t2), count in pair_counter.most_common(20):
    t1_text = decode_token(t1)
    t2_text = decode_token(t2)
    print(f"  {repr(t1_text):15s} vs {repr(t2_text):15s}: {count:5d} ({100*count/total:.1f}%)")

print("\n" + "=" * 70)
print("PHASE 8: Position-16 deep dive with semantic categories")
print("=" * 70)

# Check what's special about position 16
pos16_cats = Counter()
pos16_examples = []
non16_cats = Counter()

for prompt_idx, pos, tok1_id, tok2_id, margin, logit1, logit2 in all_fragile_pairs:
    tok1_text = decode_token(tok1_id)
    tok2_text = decode_token(tok2_id)
    cat = semantic_similarity_category(tok1_text, tok2_text)

    if pos == 16:
        pos16_cats[cat] += 1
        if len(pos16_examples) < 20:
            pos16_examples.append((prompt_idx, tok1_text, tok2_text, margin, prompts[prompt_idx]['category']))
    else:
        non16_cats[cat] += 1

total_16 = sum(pos16_cats.values())
total_non16 = sum(non16_cats.values())

if total_16 > 0:
    print(f"\nPosition 16 ({total_16} fragile instances):")
    for cat, count in pos16_cats.most_common(10):
        pct = 100 * count / total_16
        non16_pct = 100 * non16_cats.get(cat, 0) / total_non16 if total_non16 > 0 else 0
        print(f"  {cat:35s}: {pct:5.1f}% (vs {non16_pct:.1f}% elsewhere)")

    print(f"\nPosition 16 examples:")
    for prompt_idx, t1, t2, margin, cat in pos16_examples[:10]:
        print(f"  [{cat}] prompt {prompt_idx}: {repr(t1):15s} vs {repr(t2):15s} (margin={margin:.4f})")

print("\n" + "=" * 70)
print("PHASE 9: Fine-grained position fragility curve (positions 0-31)")
print("=" * 70)

pos_fragile_count = Counter()
pos_total_count = Counter()

for i in range(1999):
    path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
    if not os.path.exists(path):
        continue
    data = torch.load(path, map_location='cpu', weights_only=False)
    top5_vals = data['top5_values']
    margins = top5_vals[:, 0] - top5_vals[:, 1]

    for pos in range(min(32, len(margins))):
        pos_total_count[pos] += 1
        if margins[pos] < MARGIN_THRESHOLD:
            pos_fragile_count[pos] += 1

print("Position-by-position fragility rate (positions 0-31):")
for pos in range(32):
    rate = pos_fragile_count[pos] / pos_total_count[pos] if pos_total_count[pos] > 0 else 0
    bar = '#' * int(rate * 100)
    print(f"  pos {pos:2d}: {rate*100:5.1f}% ({pos_fragile_count[pos]:4d}/{pos_total_count[pos]:4d}) {bar}")

print("\n\nDone! Analysis complete.")
