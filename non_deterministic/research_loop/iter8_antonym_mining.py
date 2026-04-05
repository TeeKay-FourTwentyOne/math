"""
Iteration 8: Antonym Mining at Fragile Positions

Find cases where the rank-1 and rank-2 tokens at fragile positions are
semantic opposites. These are the highest-impact divergences — hardware
noise could literally invert the meaning of the output.

Approach:
1. Extract all fragile positions (margin < 1.0) across all 1999 prompts
2. Decode the rank-1 and rank-2 tokens
3. Apply multiple antonym detection methods:
   a. Curated antonym list (manually built)
   b. Prefix-based opposition (un-, in-, non-, dis-)
   c. Numerical opposition (positive/negative numbers, +/-)
   d. Boolean/logical opposition (true/false, yes/no)
4. Estimate: what fraction of hardware divergences would be meaning-inversions?
5. Find the most spectacular concrete examples
"""

import torch
import numpy as np
import json
import os
from collections import Counter, defaultdict
from pathlib import Path

LOGIT_DIR = "C:/Projects/math/non_deterministic/experiment_0/logit_traces"
PROMPTS_PATH = "C:/Projects/math/non_deterministic/prompts/candidates_2000.jsonl"

from transformers import AutoTokenizer
tokenizer = AutoTokenizer.from_pretrained("EleutherAI/pythia-410m-deduped")

prompts = {}
with open(PROMPTS_PATH, encoding="utf-8") as f:
    for line in f:
        p = json.loads(line)
        prompts[p["idx"]] = p

# ============================================================
# Curated antonym/opposition lists
# ============================================================

# Word-level antonyms (lowercase, stripped)
ANTONYM_PAIRS = set()
antonym_raw = [
    ("good", "bad"), ("big", "small"), ("large", "small"), ("hot", "cold"),
    ("up", "down"), ("left", "right"), ("in", "out"), ("yes", "no"),
    ("true", "false"), ("open", "close"), ("open", "closed"),
    ("high", "low"), ("fast", "slow"), ("old", "new"), ("old", "young"),
    ("light", "dark"), ("black", "white"), ("long", "short"),
    ("first", "last"), ("best", "worst"), ("more", "less"), ("more", "fewer"),
    ("most", "least"), ("many", "few"), ("all", "none"),
    ("always", "never"), ("before", "after"), ("above", "below"),
    ("north", "south"), ("east", "west"), ("start", "end"),
    ("begin", "end"), ("beginning", "end"), ("top", "bottom"),
    ("positive", "negative"), ("increase", "decrease"),
    ("win", "lose"), ("success", "failure"), ("right", "wrong"),
    ("correct", "incorrect"), ("possible", "impossible"),
    ("visible", "invisible"), ("happy", "sad"), ("love", "hate"),
    ("rich", "poor"), ("strong", "weak"), ("hard", "soft"),
    ("easy", "hard"), ("safe", "dangerous"), ("war", "peace"),
    ("life", "death"), ("alive", "dead"), ("accept", "reject"),
    ("agree", "disagree"), ("include", "exclude"), ("inside", "outside"),
    ("maximum", "minimum"), ("major", "minor"), ("male", "female"),
    ("man", "woman"), ("boy", "girl"), ("he", "she"),
    ("his", "her"), ("him", "her"), ("father", "mother"),
    ("largest", "smallest"), ("tallest", "shortest"), ("longest", "shortest"),
    ("highest", "lowest"), ("best", "worst"), ("most", "least"),
    ("better", "worse"), ("bigger", "smaller"), ("older", "younger"),
    ("newer", "older"), ("faster", "slower"), ("stronger", "weaker"),
    ("increase", "decrease"), ("rise", "fall"), ("grow", "shrink"),
    ("expand", "contract"), ("add", "subtract"), ("plus", "minus"),
    ("can", "cannot"), ("will", "won't"), ("do", "don't"),
    ("is", "isn't"), ("was", "wasn't"), ("has", "hasn't"),
    ("important", "unimportant"), ("known", "unknown"),
    ("common", "rare"), ("public", "private"), ("legal", "illegal"),
    ("formal", "informal"), ("internal", "external"),
    ("ancient", "modern"), ("rural", "urban"), ("domestic", "foreign"),
    ("simple", "complex"), ("natural", "artificial"),
    ("the", "a"),  # Not antonyms but maximally different determiners
]
for a, b in antonym_raw:
    ANTONYM_PAIRS.add((a, b))
    ANTONYM_PAIRS.add((b, a))

def is_antonym_pair(s1, s2):
    """Check if two strings are antonym pairs."""
    w1 = s1.strip().lower()
    w2 = s2.strip().lower()
    return (w1, w2) in ANTONYM_PAIRS

def is_prefix_opposition(s1, s2):
    """Check if one token is the negation-prefix version of the other."""
    w1 = s1.strip().lower()
    w2 = s2.strip().lower()
    prefixes = ["un", "in", "im", "ir", "il", "non", "dis", "anti"]
    for p in prefixes:
        if w1 == p + w2 or w2 == p + w1:
            return True
        if w1.startswith(p) and w2 == w1[len(p):]:
            return True
        if w2.startswith(p) and w1 == w2[len(p):]:
            return True
    return False

def is_numeric_opposition(s1, s2):
    """Check for +/-, positive/negative numbers."""
    w1, w2 = s1.strip(), s2.strip()
    # +/-
    if (w1 == '+' and w2 == '-') or (w1 == '-' and w2 == '+'):
        return True
    # Positive/negative numbers
    try:
        n1, n2 = float(w1), float(w2)
        if n1 * n2 < 0:  # opposite signs
            return True
    except ValueError:
        pass
    return False

def is_boolean_opposition(s1, s2):
    """Check for true/false, yes/no, etc."""
    w1, w2 = s1.strip().lower(), s2.strip().lower()
    bool_pairs = {("true", "false"), ("yes", "no"), ("True", "False"),
                  ("1", "0"), ("correct", "incorrect"), ("right", "wrong")}
    return (w1, w2) in bool_pairs or (w2, w1) in bool_pairs

# ============================================================
# Extract all fragile position pairs
# ============================================================
print("Extracting fragile position pairs...", flush=True)

all_pairs = []  # (prompt_idx, pos, margin, tok1_str, tok2_str, category)
total_fragile = 0
total_positions = 0

for i in range(1999):
    path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
    if not os.path.exists(path):
        continue
    data = torch.load(path, map_location="cpu", weights_only=True)
    top5_vals = data["top5_values"]
    top5_idxs = data["top5_indices"]
    margins = (top5_vals[:, 0] - top5_vals[:, 1]).numpy()
    cat = prompts[i]["category"]
    total_positions += len(margins)

    for pos in range(len(margins)):
        if margins[pos] < 1.0:
            total_fragile += 1
            t1 = tokenizer.decode([top5_idxs[pos, 0].item()])
            t2 = tokenizer.decode([top5_idxs[pos, 1].item()])
            all_pairs.append({
                "idx": i, "pos": pos, "margin": float(margins[pos]),
                "t1": t1, "t2": t2, "cat": cat
            })

print(f"Total fragile positions (margin < 1.0): {total_fragile}")
print(f"Total positions: {total_positions}")
print(f"Fragile rate: {total_fragile/total_positions:.4f}")

# ============================================================
# PART 1: Classify all fragile pairs
# ============================================================
print("\n" + "=" * 70)
print("PART 1: Opposition classification of fragile pairs")
print("=" * 70)

n_antonym = 0
n_prefix = 0
n_numeric = 0
n_boolean = 0
antonym_examples = []
prefix_examples = []
numeric_examples = []

for p in all_pairs:
    t1, t2 = p["t1"], p["t2"]
    if is_antonym_pair(t1, t2):
        n_antonym += 1
        antonym_examples.append(p)
    if is_prefix_opposition(t1, t2):
        n_prefix += 1
        prefix_examples.append(p)
    if is_numeric_opposition(t1, t2):
        n_numeric += 1
        numeric_examples.append(p)
    if is_boolean_opposition(t1, t2):
        n_boolean += 1

n_any_opposition = 0
for p in all_pairs:
    t1, t2 = p["t1"], p["t2"]
    if (is_antonym_pair(t1, t2) or is_prefix_opposition(t1, t2) or
        is_numeric_opposition(t1, t2) or is_boolean_opposition(t1, t2)):
        n_any_opposition += 1

N = len(all_pairs)
print(f"\nTotal fragile positions analyzed: {N}")
print(f"  Curated antonyms:     {n_antonym:6d} ({n_antonym/N*100:.2f}%)")
print(f"  Prefix opposition:    {n_prefix:6d} ({n_prefix/N*100:.2f}%)")
print(f"  Numeric opposition:   {n_numeric:6d} ({n_numeric/N*100:.2f}%)")
print(f"  Boolean opposition:   {n_boolean:6d} ({n_boolean/N*100:.2f}%)")
print(f"  ANY opposition:       {n_any_opposition:6d} ({n_any_opposition/N*100:.2f}%)")

# ============================================================
# PART 2: Most common antonym pairs
# ============================================================
print("\n" + "=" * 70)
print("PART 2: Most common antonym pairs at fragile positions")
print("=" * 70)

antonym_counter = Counter()
for p in antonym_examples:
    pair = tuple(sorted([p["t1"].strip().lower(), p["t2"].strip().lower()]))
    antonym_counter[pair] += 1

print(f"\nTop 30 antonym pairs:")
for (w1, w2), count in antonym_counter.most_common(30):
    print(f"  '{w1}' vs '{w2}': {count}")

# ============================================================
# PART 3: By margin bin — are oppositions more common at tiny margins?
# ============================================================
print("\n" + "=" * 70)
print("PART 3: Opposition rate by margin size")
print("=" * 70)

bins = [(0, 0.01), (0.01, 0.05), (0.05, 0.1), (0.1, 0.5), (0.5, 1.0)]
for lo, hi in bins:
    group = [p for p in all_pairs if lo <= p["margin"] < hi]
    if not group:
        continue
    n_opp = sum(1 for p in group if
                is_antonym_pair(p["t1"], p["t2"]) or is_prefix_opposition(p["t1"], p["t2"]) or
                is_numeric_opposition(p["t1"], p["t2"]) or is_boolean_opposition(p["t1"], p["t2"]))
    print(f"  margin [{lo:.2f}, {hi:.2f}): {len(group):6d} pairs, "
          f"{n_opp:4d} oppositions ({n_opp/len(group)*100:.2f}%)")

# ============================================================
# PART 4: Concrete antonym examples (with context)
# ============================================================
print("\n" + "=" * 70)
print("PART 4: Spectacular antonym examples (sorted by margin)")
print("=" * 70)

# Sort antonym examples by margin (smallest = most likely to flip)
antonym_examples.sort(key=lambda x: x["margin"])

print(f"\nTop 30 most fragile antonym pairs (smallest margin):")
seen_prompts = set()
for p in antonym_examples:
    if p["idx"] in seen_prompts:
        continue  # one example per prompt
    seen_prompts.add(p["idx"])
    prompt_text = prompts[p["idx"]]["text"][:60]
    # Get surrounding context from generated tokens
    trace = torch.load(os.path.join(LOGIT_DIR, f"prompt_{p['idx']:04d}.pt"),
                       map_location="cpu", weights_only=True)
    gen_ids = trace["generated_ids"].tolist()
    input_len = len(trace["input_ids"].tolist())
    # Tokens around the fragile position
    ctx_start = max(0, p["pos"] - 3)
    ctx_end = min(len(gen_ids) - input_len, p["pos"] + 4)
    context_ids = gen_ids[input_len + ctx_start:input_len + ctx_end]
    context = tokenizer.decode(context_ids).replace("\n", "\\n")

    print(f"\n  idx={p['idx']:4d} pos={p['pos']:3d} margin={p['margin']:.6f} [{p['cat']}]")
    print(f"    FLIP: '{p['t1']}' <-> '{p['t2']}'")
    print(f"    Context: ...{context}...")
    print(f"    Prompt: {prompt_text}")

    if len(seen_prompts) >= 30:
        break

# ============================================================
# PART 5: Per-category and per-threshold opposition rates
# ============================================================
print("\n" + "=" * 70)
print("PART 5: Opposition rate by category")
print("=" * 70)

for cat in sorted(set(p["cat"] for p in all_pairs)):
    cat_pairs = [p for p in all_pairs if p["cat"] == cat]
    n_opp = sum(1 for p in cat_pairs if
                is_antonym_pair(p["t1"], p["t2"]) or is_prefix_opposition(p["t1"], p["t2"]) or
                is_numeric_opposition(p["t1"], p["t2"]) or is_boolean_opposition(p["t1"], p["t2"]))
    print(f"  {cat:15s}: {len(cat_pairs):6d} fragile, {n_opp:4d} oppositions "
          f"({n_opp/max(len(cat_pairs),1)*100:.2f}%)")

# ============================================================
# PART 6: Impact estimate — what fraction of FIRST divergences are oppositions?
# ============================================================
print("\n" + "=" * 70)
print("PART 6: Opposition rate among FIRST divergence points")
print("=" * 70)

for eps in [1e-3, 1e-2, 1e-1]:
    first_flips = []
    for i in range(1999):
        path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
        if not os.path.exists(path):
            continue
        data = torch.load(path, map_location="cpu", weights_only=True)
        top5_vals = data["top5_values"]
        top5_idxs = data["top5_indices"]
        margins = (top5_vals[:, 0] - top5_vals[:, 1]).numpy()
        flip_pos = np.where(margins < eps)[0]
        if len(flip_pos) > 0:
            pos = flip_pos[0]
            t1 = tokenizer.decode([top5_idxs[pos, 0].item()])
            t2 = tokenizer.decode([top5_idxs[pos, 1].item()])
            is_opp = (is_antonym_pair(t1, t2) or is_prefix_opposition(t1, t2) or
                      is_numeric_opposition(t1, t2) or is_boolean_opposition(t1, t2))
            first_flips.append({"idx": i, "t1": t1, "t2": t2, "is_opp": is_opp,
                               "margin": float(margins[pos]), "pos": pos})

    n_opp = sum(1 for f in first_flips if f["is_opp"])
    print(f"\n  eps={eps:.0e}: {len(first_flips)} diverging prompts, "
          f"{n_opp} ({n_opp/max(len(first_flips),1)*100:.1f}%) have opposition at first flip")

    # Show some opposition first-flips
    opp_flips = [f for f in first_flips if f["is_opp"]]
    if opp_flips:
        opp_flips.sort(key=lambda x: x["margin"])
        print(f"  Examples:")
        for f in opp_flips[:5]:
            print(f"    idx={f['idx']:4d} pos={f['pos']:3d} '{f['t1']}'->{f['t2']}' "
                  f"margin={f['margin']:.6f} [{prompts[f['idx']]['category']}]")
