"""
Deep analysis of digit-vs-digit fragility.

Digit-vs-digit is 35.5% of all fragile positions. This analysis examines:
1. Which specific digit pairs compete most
2. Whether there's a pattern to which digits get confused
3. Whether adjacent digits (e.g., 4 vs 5) compete more than distant ones
4. The relationship to arithmetic problem types
"""

import torch
import numpy as np
import json
import os
from collections import Counter, defaultdict

from transformers import AutoTokenizer
tokenizer = AutoTokenizer.from_pretrained("EleutherAI/pythia-410m")

LOGIT_DIR = "C:/Projects/math/non_deterministic/experiment_0/logit_traces"
PROMPTS_FILE = "C:/Projects/math/non_deterministic/prompts/candidates_2000.jsonl"
MARGIN_THRESHOLD = 0.1

prompts = {}
with open(PROMPTS_FILE) as f:
    for line in f:
        d = json.loads(line)
        prompts[d['idx']] = d

def decode(tid):
    return tokenizer.decode([tid])

# Collect all digit-vs-digit pairs
digit_pairs = []
# Map token id to its numeric value
digit_token_values = {}

print("Phase 1: Collecting all digit-vs-digit fragile positions...")
for i in range(1999):
    path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
    if not os.path.exists(path):
        continue
    data = torch.load(path, map_location='cpu', weights_only=False)
    top5_vals = data['top5_values']
    top5_idx = data['top5_indices']
    margins = top5_vals[:, 0] - top5_vals[:, 1]

    for pos in range(len(margins)):
        if margins[pos] < MARGIN_THRESHOLD:
            t1 = decode(top5_idx[pos, 0].item())
            t2 = decode(top5_idx[pos, 1].item())
            t1s = t1.strip()
            t2s = t2.strip()
            if t1s.isdigit() and t2s.isdigit():
                v1 = int(t1s)
                v2 = int(t2s)
                digit_pairs.append((i, pos, v1, v2, margins[pos].item(),
                                     top5_idx[pos, 0].item(), top5_idx[pos, 1].item(),
                                     prompts[i]['category']))

print(f"Total digit-vs-digit fragile positions: {len(digit_pairs)}")

# Analysis 1: Digit value distribution
print("\n" + "=" * 70)
print("Digit value distribution at fragile positions")
print("=" * 70)

val1_counter = Counter()
val2_counter = Counter()
for _, _, v1, v2, _, _, _, _ in digit_pairs:
    val1_counter[v1] += 1
    val2_counter[v2] += 1

print("\nRank-1 digit values:")
for v, c in sorted(val1_counter.items()):
    bar = '#' * (c // 20)
    print(f"  {v:4d}: {c:5d} {bar}")

print("\nRank-2 digit values:")
for v, c in sorted(val2_counter.items()):
    bar = '#' * (c // 20)
    print(f"  {v:4d}: {c:5d} {bar}")

# Analysis 2: Digit distance (|v1 - v2|)
print("\n" + "=" * 70)
print("Digit distance analysis")
print("=" * 70)

distances = Counter()
for _, _, v1, v2, _, _, _, _ in digit_pairs:
    d = abs(v1 - v2)
    distances[d] += 1

total = len(digit_pairs)
print("\nAbsolute distance between competing digits:")
for d, c in sorted(distances.items()):
    pct = 100 * c / total
    bar = '#' * int(pct * 2)
    print(f"  |d|={d:3d}: {c:5d} ({pct:5.1f}%) {bar}")

# Analysis 3: Specific pair heatmap (for single digits 0-9)
print("\n" + "=" * 70)
print("Single-digit pair confusion matrix (row=rank1, col=rank2)")
print("=" * 70)

pair_matrix = np.zeros((10, 10), dtype=int)
for _, _, v1, v2, _, _, _, _ in digit_pairs:
    if 0 <= v1 <= 9 and 0 <= v2 <= 9:
        pair_matrix[v1][v2] += 1

# Print matrix
header = "     " + " ".join(f"{j:5d}" for j in range(10))
print(header)
for i in range(10):
    row = f"  {i}: " + " ".join(f"{pair_matrix[i][j]:5d}" for j in range(10))
    print(row)

# Analysis 4: Multi-digit tokens
print("\n" + "=" * 70)
print("Multi-digit token competition")
print("=" * 70)

multi_digit_pairs = Counter()
for _, _, v1, v2, m, _, _, _ in digit_pairs:
    if v1 >= 10 or v2 >= 10:
        multi_digit_pairs[(v1, v2)] += 1

print("Most common multi-digit competing pairs:")
for (v1, v2), c in multi_digit_pairs.most_common(15):
    print(f"  {v1:4d} vs {v2:4d}: {c:5d}")

# Analysis 5: By prompt category
print("\n" + "=" * 70)
print("Digit-vs-digit by prompt category")
print("=" * 70)

cat_digit_counts = Counter()
cat_total = Counter()
for _, _, v1, v2, _, _, _, cat in digit_pairs:
    cat_digit_counts[cat] += 1

# Get total fragile positions per category (from all types, not just digits)
import csv
with open("C:/Projects/math/non_deterministic/prompts/fragility_scores.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        cat = row['category']
        cat_total[cat] += int(row['fragile_positions'])

for cat in sorted(cat_digit_counts.keys()):
    d = cat_digit_counts[cat]
    t = cat_total.get(cat, 1)
    print(f"  {cat:20s}: {d:5d} digit-vs-digit of {t:6d} total fragile ({100*d/t:.1f}%)")

# Analysis 6: Position dependence of digit fragility
print("\n" + "=" * 70)
print("Position dependence of digit-vs-digit fragility")
print("=" * 70)

pos_digit = Counter()
for _, pos, _, _, _, _, _, _ in digit_pairs:
    pos_digit[pos] += 1

pos_bins = [(0, 15), (16, 31), (32, 63), (64, 127), (128, 255)]
for lo, hi in pos_bins:
    count = sum(pos_digit.get(p, 0) for p in range(lo, hi + 1))
    pct = 100 * count / total
    print(f"  pos {lo:3d}-{hi:3d}: {count:5d} ({pct:.1f}%)")

# Analysis 7: Are digit confusions in arithmetic correlated with wrong answers?
# Look for patterns: is the model getting arithmetic wrong at fragile positions?
print("\n" + "=" * 70)
print("Arithmetic context examples (digit confusions in reasoning prompts)")
print("=" * 70)

reasoning_digit_pairs = [(pi, pos, v1, v2, m) for pi, pos, v1, v2, m, _, _, cat in digit_pairs
                          if cat == 'reasoning']

shown = 0
for pi, pos, v1, v2, m in reasoning_digit_pairs[:30]:
    path = os.path.join(LOGIT_DIR, f"prompt_{pi:04d}.pt")
    data = torch.load(path, map_location='cpu', weights_only=False)
    gen_ids = data['generated_ids']
    input_len = len(data['input_ids'])
    actual_pos = input_len + pos

    ctx_start = max(0, actual_pos - 5)
    ctx_end = min(len(gen_ids), actual_pos + 3)
    context = tokenizer.decode(gen_ids[ctx_start:actual_pos])
    after = tokenizer.decode(gen_ids[actual_pos:ctx_end])

    if shown < 15:
        shown += 1
        print(f"\n  p{pi} @{pos}: {v1} vs {v2} (m={m:.4f})")
        print(f"    Prompt: '{prompts[pi]['text'][:80]}'")
        print(f"    Context: '...{context}[{v1}|{v2}]{after}...'")

# Analysis 8: Are adjacent digits (d=1) preferentially confused?
print("\n" + "=" * 70)
print("Adjacent vs non-adjacent digit confusion rates")
print("=" * 70)

adjacent = sum(1 for _, _, v1, v2, _, _, _, _ in digit_pairs if abs(v1 - v2) == 1)
non_adj = total - adjacent
print(f"Adjacent (|d|=1): {adjacent} ({100*adjacent/total:.1f}%)")
print(f"Non-adjacent (|d|>1): {non_adj} ({100*non_adj/total:.1f}%)")

# Expected rates if random: for single digits 0-9, there are 18 adjacent pairs
# and 72 non-adjacent pairs. So random expected adjacent rate = 18/90 = 20%.
print(f"\nExpected if random (single-digit): ~20% adjacent")
single_digit_only = sum(1 for _, _, v1, v2, _, _, _, _ in digit_pairs if 0 <= v1 <= 9 and 0 <= v2 <= 9)
adj_single = sum(1 for _, _, v1, v2, _, _, _, _ in digit_pairs
                 if 0 <= v1 <= 9 and 0 <= v2 <= 9 and abs(v1 - v2) == 1)
if single_digit_only > 0:
    print(f"Actual adjacent rate (single-digit only): {100*adj_single/single_digit_only:.1f}%")

print("\n\nDone!")
