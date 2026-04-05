"""
Position-dependent semantic shift analysis.

Key finding from Phase 3: Early positions (0-31) are dominated by word-vs-word
fragility, but later positions (128+) are dominated by digit-vs-digit.

This analysis investigates:
1. Is the shift gradual or abrupt?
2. Is it category-dependent? (reasoning prompts = digits later, etc.)
3. What's the crossover point?
4. What happens in the transition zone?
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

def classify_broad(t1_text, t2_text):
    """Broad classification: word, digit, punctuation, other."""
    t1s = t1_text.strip()
    t2s = t2_text.strip()

    def tok_type(s):
        if not s:
            return 'whitespace'
        if s.isdigit():
            return 'digit'
        if s.isalpha():
            return 'word'
        if all(c in '.,;:!?-\'"()[]{}' for c in s):
            return 'punct'
        return 'other'

    c1 = tok_type(t1s)
    c2 = tok_type(t2s)
    if c1 == c2:
        return c1 + '_vs_' + c2
    return 'cross_type'

# Collect position-by-position fragility with type
pos_type_counts = defaultdict(lambda: Counter())  # pos -> Counter of types
pos_cat_type_counts = defaultdict(lambda: Counter())  # (cat, pos) -> Counter of types

print("Collecting fragile positions by type and position...")
for i in range(1999):
    path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
    if not os.path.exists(path):
        continue
    data = torch.load(path, map_location='cpu', weights_only=False)
    top5_vals = data['top5_values']
    top5_idx = data['top5_indices']
    margins = top5_vals[:, 0] - top5_vals[:, 1]
    cat = prompts[i]['category']

    for pos in range(len(margins)):
        if margins[pos] < MARGIN_THRESHOLD:
            t1 = decode(top5_idx[pos, 0].item())
            t2 = decode(top5_idx[pos, 1].item())
            broad = classify_broad(t1, t2)
            pos_type_counts[pos][broad] += 1
            pos_cat_type_counts[(cat, pos)][broad] += 1

# Print per-position type fractions (binned by 8-position windows)
print("\n" + "=" * 70)
print("Fragile position types by 8-position windows")
print("=" * 70)

window_size = 8
print(f"\n{'Window':>10s} {'Total':>6s} {'word%':>6s} {'digit%':>7s} {'punct%':>7s} {'cross%':>7s} {'other%':>7s}")
for w_start in range(0, 256, window_size):
    w_end = w_start + window_size - 1
    merged = Counter()
    for p in range(w_start, w_end + 1):
        merged += pos_type_counts[p]
    total = sum(merged.values())
    if total < 10:
        continue
    word_pct = 100 * merged.get('word_vs_word', 0) / total
    digit_pct = 100 * merged.get('digit_vs_digit', 0) / total
    punct_pct = 100 * merged.get('punct_vs_punct', 0) / total
    cross_pct = 100 * merged.get('cross_type', 0) / total
    other_pct = 100 * (merged.get('other_vs_other', 0) + merged.get('whitespace_vs_whitespace', 0)) / total
    label = f"{w_start:3d}-{w_end:3d}"
    # Visual bars
    w_bar = '#' * int(word_pct / 2)
    d_bar = '.' * int(digit_pct / 2)
    print(f"{label:>10s} {total:6d} {word_pct:5.1f}% {digit_pct:6.1f}% {punct_pct:6.1f}% {cross_pct:6.1f}% {other_pct:6.1f}%  W:{w_bar} D:{d_bar}")

# Per category
print("\n" + "=" * 70)
print("Per-category position dependence")
print("=" * 70)

for cat in ['factual', 'reasoning', 'continuation', 'conversational', 'open_ended']:
    print(f"\n--- {cat} ---")
    print(f"{'Window':>10s} {'Total':>6s} {'word%':>6s} {'digit%':>7s}")
    for w_start in range(0, 256, 16):
        w_end = w_start + 15
        merged = Counter()
        for p in range(w_start, w_end + 1):
            merged += pos_cat_type_counts.get((cat, p), Counter())
        total = sum(merged.values())
        if total < 5:
            continue
        word_pct = 100 * merged.get('word_vs_word', 0) / total
        digit_pct = 100 * merged.get('digit_vs_digit', 0) / total
        print(f"{w_start:3d}-{w_end:3d}   {total:6d} {word_pct:5.1f}% {digit_pct:6.1f}%")

# Find the crossover point precisely
print("\n" + "=" * 70)
print("Crossover point: where does digit fragility exceed word fragility?")
print("=" * 70)

# Use a sliding window of 16 positions
for w_start in range(0, 240, 4):
    w_end = w_start + 15
    merged = Counter()
    for p in range(w_start, w_end + 1):
        merged += pos_type_counts[p]
    total = sum(merged.values())
    if total < 10:
        continue
    word_frac = merged.get('word_vs_word', 0) / total
    digit_frac = merged.get('digit_vs_digit', 0) / total
    if digit_frac > word_frac and w_start >= 16:
        print(f"Crossover at positions {w_start}-{w_end}: "
              f"word={100*word_frac:.1f}%, digit={100*digit_frac:.1f}%")
        break

# Also compute the overall summary
print("\n" + "=" * 70)
print("Summary: Semantic character of divergences by position")
print("=" * 70)

early_word = sum(pos_type_counts[p].get('word_vs_word', 0) for p in range(32))
early_digit = sum(pos_type_counts[p].get('digit_vs_digit', 0) for p in range(32))
early_total = sum(sum(pos_type_counts[p].values()) for p in range(32))

late_word = sum(pos_type_counts[p].get('word_vs_word', 0) for p in range(128, 256))
late_digit = sum(pos_type_counts[p].get('digit_vs_digit', 0) for p in range(128, 256))
late_total = sum(sum(pos_type_counts[p].values()) for p in range(128, 256))

print(f"\nEarly (pos 0-31): {early_total} fragile positions")
print(f"  word-vs-word: {early_word} ({100*early_word/early_total:.1f}%)")
print(f"  digit-vs-digit: {early_digit} ({100*early_digit/early_total:.1f}%)")

print(f"\nLate (pos 128-255): {late_total} fragile positions")
print(f"  word-vs-word: {late_word} ({100*late_word/late_total:.1f}%)")
print(f"  digit-vs-digit: {late_digit} ({100*late_digit/late_total:.1f}%)")

print(f"\nword dominance shift: {100*early_word/early_total:.1f}% -> {100*late_word/late_total:.1f}%")
print(f"digit dominance shift: {100*early_digit/early_total:.1f}% -> {100*late_digit/late_total:.1f}%")

print("\nDone!")
