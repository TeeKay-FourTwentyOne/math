"""
Iteration 2 follow-up: Why is position 2 at 82.9% fragility?

Position 0 is a newline 83% of the time (low fragility: 17%).
Position 2 spikes to 82.9%. What's happening?

Hypothesis: Position 0=newline, position 1=first significant token,
position 2=second token of first phrase. The model's uncertainty peaks
at the second content token where it must commit to a clause structure.

Let's investigate:
1. What tokens appear at positions 0, 1, 2, 3?
2. What's the competing pair at position 2?
3. Does the spike depend on what position 0 is?
4. Is position 2 fragility caused by a specific category?
"""

import torch
import numpy as np
import json
import os
from collections import Counter, defaultdict

LOGIT_DIR = "C:/Projects/math/non_deterministic/experiment_0/logit_traces"
PROMPTS_PATH = "C:/Projects/math/non_deterministic/prompts/candidates_2000.jsonl"

from transformers import AutoTokenizer
tokenizer = AutoTokenizer.from_pretrained("EleutherAI/pythia-410m-deduped")

prompts = {}
with open(PROMPTS_PATH) as f:
    for line in f:
        p = json.loads(line)
        prompts[p["idx"]] = p

# Collect data at positions 0, 1, 2, 3
pos_data = {pos: {"tokens": [], "margins": [], "top2": [], "cats": []} for pos in range(4)}

for i in range(1999):
    path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
    if not os.path.exists(path):
        continue
    data = torch.load(path, map_location="cpu", weights_only=True)
    top5_vals = data["top5_values"]
    top5_idxs = data["top5_indices"]
    margins = (top5_vals[:, 0] - top5_vals[:, 1]).numpy()
    cat = prompts[i]["category"]

    for pos in range(min(4, len(margins))):
        pos_data[pos]["tokens"].append(top5_idxs[pos, 0].item())
        pos_data[pos]["margins"].append(margins[pos])
        pos_data[pos]["top2"].append((top5_idxs[pos, 0].item(), top5_idxs[pos, 1].item()))
        pos_data[pos]["cats"].append(cat)

# ============================================================
# Token distribution at positions 0-3
# ============================================================
print("=" * 70)
print("TOKEN DISTRIBUTION AT POSITIONS 0-3")
print("=" * 70)

for pos in range(4):
    d = pos_data[pos]
    m = np.array(d["margins"])
    frag_rate = (m < 1.0).mean()
    print(f"\nPosition {pos}: fragility={frag_rate:.3f}, mean_margin={m.mean():.3f}")
    tc = Counter(d["tokens"])
    print(f"  Top 10 tokens:")
    for tid, count in tc.most_common(10):
        ts = tokenizer.decode([tid]).replace("\n", "\\n")
        f = count / len(d["tokens"])
        print(f"    '{ts}' ({tid}): {count} ({f:.3f})")

# ============================================================
# What COMPETES at position 2?
# ============================================================
print("\n" + "=" * 70)
print("COMPETING PAIRS AT POSITION 2 (fragile only, margin < 1.0)")
print("=" * 70)

d2 = pos_data[2]
m2 = np.array(d2["margins"])
fragile_pairs = [(d2["top2"][i], d2["margins"][i], d2["cats"][i])
                 for i in range(len(d2["margins"])) if d2["margins"][i] < 1.0]

pair_counter = Counter()
for (t1, t2), margin, cat in fragile_pairs:
    s1 = tokenizer.decode([t1]).replace("\n", "\\n")
    s2 = tokenizer.decode([t2]).replace("\n", "\\n")
    pair_counter[(s1, s2)] += 1

print(f"Fragile at pos 2: {len(fragile_pairs)}/{len(d2['margins'])} = {len(fragile_pairs)/len(d2['margins']):.3f}")
print(f"\nTop 20 competing pairs:")
for (s1, s2), count in pair_counter.most_common(20):
    print(f"  '{s1}' vs '{s2}': {count}")

# ============================================================
# Per-category fragility at position 2
# ============================================================
print("\n" + "=" * 70)
print("PER-CATEGORY FRAGILITY AT POSITION 2")
print("=" * 70)

cat_frags = defaultdict(list)
for i in range(len(d2["margins"])):
    cat_frags[d2["cats"][i]].append(d2["margins"][i])

for cat in sorted(cat_frags.keys()):
    m = np.array(cat_frags[cat])
    rate = (m < 1.0).mean()
    print(f"  {cat:15s}: {rate:.4f} ({int((m < 1.0).sum())}/{len(m)})")

# ============================================================
# Does pos 2 fragility depend on what pos 0 was?
# ============================================================
print("\n" + "=" * 70)
print("POS 2 FRAGILITY CONDITIONED ON POS 0 TOKEN")
print("=" * 70)

pos0_tok = pos_data[0]["tokens"]
pos2_m = np.array(pos_data[2]["margins"])

# Group by pos 0 token
p0_groups = defaultdict(list)
for i in range(min(len(pos0_tok), len(pos2_m))):
    p0_groups[pos0_tok[i]].append(pos2_m[i])

for tid, ms in sorted(p0_groups.items(), key=lambda x: -len(x[1]))[:8]:
    ms_arr = np.array(ms)
    rate = (ms_arr < 1.0).mean()
    ts = tokenizer.decode([tid]).replace("\n", "\\n")
    print(f"  pos0='{ts}' ({len(ms)} prompts): pos2_fragility={rate:.4f}")

# ============================================================
# Concrete examples: what does the model generate at positions 0-5?
# ============================================================
print("\n" + "=" * 70)
print("CONCRETE EXAMPLES: First 6 generated tokens (5 high-fragility, 5 low-fragility at pos 2)")
print("=" * 70)

# Sort prompts by margin at position 2
indexed = [(i, pos_data[2]["margins"][i]) for i in range(len(pos_data[2]["margins"]))]
indexed.sort(key=lambda x: x[1])

print("\n--- Most fragile at position 2 (smallest margins) ---")
for rank, (i, margin) in enumerate(indexed[:5]):
    path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
    data = torch.load(path, map_location="cpu", weights_only=True)
    gen_ids = data["generated_ids"].tolist()
    first6 = gen_ids[:6]
    decoded = [tokenizer.decode([t]).replace("\n", "\\n") for t in first6]
    t2_1 = tokenizer.decode([data["top5_indices"][2, 0].item()])
    t2_2 = tokenizer.decode([data["top5_indices"][2, 1].item()])
    print(f"  #{rank+1} idx={i} margin={margin:.6f} [{prompts[i]['category']}]")
    print(f"    Prompt: {prompts[i]['text'][:70]}")
    print(f"    First 6 tokens: {decoded}")
    print(f"    Pos 2 contest: '{t2_1}' vs '{t2_2}'")

print("\n--- Least fragile at position 2 (largest margins) ---")
for rank, (i, margin) in enumerate(indexed[-5:]):
    path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
    data = torch.load(path, map_location="cpu", weights_only=True)
    gen_ids = data["generated_ids"].tolist()
    first6 = gen_ids[:6]
    decoded = [tokenizer.decode([t]).replace("\n", "\\n") for t in first6]
    print(f"  #{rank+1} idx={i} margin={margin:.4f} [{prompts[i]['category']}]")
    print(f"    Prompt: {prompts[i]['text'][:70]}")
    print(f"    First 6 tokens: {decoded}")
