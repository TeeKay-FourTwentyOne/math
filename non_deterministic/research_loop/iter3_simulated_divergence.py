"""
Iteration 3: Simulated Cross-Platform Divergence

We have top-5 logits at every generation step. We can simulate what happens
when another platform introduces small FP noise:
- If noise > margin at a position, the token flips from rank-1 to rank-2
- The FIRST flip point is reliable (original logits are valid up to that point)
- After the first flip, all subsequent logits would be different (autoregressive)

Analysis:
1. First-divergence-point distribution at various noise levels
2. Cascade analysis at threshold 1.0 (open question from STATE.md)
3. What token flips at the first divergence point?
4. Semantic impact: is the first flip cosmetic or meaning-changing?
5. Per-category and experimental/control breakdown
"""

import torch
import numpy as np
import json
import os
from collections import Counter, defaultdict
from pathlib import Path

LOGIT_DIR = "C:/Projects/math/non_deterministic/experiment_0/logit_traces"
PROMPTS_PATH = "C:/Projects/math/non_deterministic/prompts/candidates_2000.jsonl"
EXP_PATH = "C:/Projects/math/non_deterministic/prompts/experimental_200.jsonl"
CTRL_PATH = "C:/Projects/math/non_deterministic/prompts/control_200.jsonl"

from transformers import AutoTokenizer
tokenizer = AutoTokenizer.from_pretrained("EleutherAI/pythia-410m-deduped")

# Load prompt metadata
prompts = {}
with open(PROMPTS_PATH) as f:
    for line in f:
        p = json.loads(line)
        prompts[p["idx"]] = p

# Build experimental/control index
exp_set = set()
ctrl_set = set()
for path, s in [(EXP_PATH, exp_set), (CTRL_PATH, ctrl_set)]:
    with open(path) as f:
        for line in f:
            d = json.loads(line)
            s.add(d["prompt_idx"])

# ============================================================
# Load all logit traces (margins and top-2 token ids)
# ============================================================
print("Loading logit traces...", flush=True)
all_margins = []  # list of arrays, one per prompt
all_top2_ids = []  # list of (n_steps, 2) arrays
all_idxs = []

for i in range(1999):
    path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
    if not os.path.exists(path):
        continue
    data = torch.load(path, map_location="cpu", weights_only=True)
    top5_vals = data["top5_values"]
    top5_idxs = data["top5_indices"]
    margins = (top5_vals[:, 0] - top5_vals[:, 1]).numpy()
    top2 = top5_idxs[:, :2].numpy()  # (n_steps, 2)
    all_margins.append(margins)
    all_top2_ids.append(top2)
    all_idxs.append(i)

N = len(all_margins)
print(f"Loaded {N} prompts", flush=True)

# ============================================================
# PART 1: Cascade analysis at threshold 1.0
# ============================================================
print("\n" + "=" * 70)
print("PART 1: Cascade (consecutive fragile run) analysis at threshold 1.0")
print("=" * 70)

run_lengths_all = []
max_runs = []
prompts_with_run_ge = {k: 0 for k in [1, 3, 5, 10, 20, 50]}

for margins in all_margins:
    fragile = margins < 1.0
    # Find consecutive runs
    runs = []
    current_run = 0
    for f in fragile:
        if f:
            current_run += 1
        else:
            if current_run > 0:
                runs.append(current_run)
            current_run = 0
    if current_run > 0:
        runs.append(current_run)

    if runs:
        run_lengths_all.extend(runs)
        max_run = max(runs)
        max_runs.append(max_run)
        for k in prompts_with_run_ge:
            if max_run >= k:
                prompts_with_run_ge[k] += 1
    else:
        max_runs.append(0)

rl = np.array(run_lengths_all)
print(f"Total runs: {len(rl)}")
print(f"Mean run length: {rl.mean():.2f}")
print(f"Median run length: {np.median(rl):.1f}")
print(f"Max run length: {rl.max()}")
print(f"Std: {rl.std():.2f}")

print(f"\nRun length distribution:")
for length in range(1, 21):
    count = (rl == length).sum()
    pct = count / len(rl) * 100
    bar = "#" * int(pct)
    print(f"  {length:3d}: {count:5d} ({pct:5.1f}%) {bar}")
longer = (rl > 20).sum()
print(f"  >20: {longer:5d} ({longer/len(rl)*100:5.1f}%)")

print(f"\nPrompts with max run >= K:")
for k, count in sorted(prompts_with_run_ge.items()):
    print(f"  >= {k:3d}: {count:5d} ({count/N*100:.1f}%)")

# ============================================================
# PART 2: First divergence point at various noise levels
# ============================================================
print("\n" + "=" * 70)
print("PART 2: Simulated first-divergence-point distribution")
print("=" * 70)

noise_levels = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1, 1.0]

for eps in noise_levels:
    first_div_points = []
    n_diverged = 0
    n_exp_div = 0
    n_ctrl_div = 0

    for j, (margins, idx) in enumerate(zip(all_margins, all_idxs)):
        # Find first position where margin < eps
        flip_positions = np.where(margins < eps)[0]
        if len(flip_positions) > 0:
            first_pos = flip_positions[0]
            first_div_points.append(first_pos)
            n_diverged += 1
            if idx in exp_set:
                n_exp_div += 1
            if idx in ctrl_set:
                n_ctrl_div += 1

    if n_diverged > 0:
        fdp = np.array(first_div_points)
        print(f"\n  eps={eps:.0e}: {n_diverged}/{N} diverge ({n_diverged/N*100:.1f}%)"
              f"  [exp: {n_exp_div}/200, ctrl: {n_ctrl_div}/200]")
        print(f"    First div point: mean={fdp.mean():.1f}  median={np.median(fdp):.1f}"
              f"  min={fdp.min()}  max={fdp.max()}")
        # Distribution in bins
        for lo, hi in [(0, 2), (3, 7), (8, 15), (16, 31), (32, 63), (64, 127), (128, 255)]:
            count = ((fdp >= lo) & (fdp <= hi)).sum()
            if count > 0:
                print(f"    pos {lo:3d}-{hi:3d}: {count:4d} ({count/len(fdp)*100:.1f}%)")
    else:
        print(f"\n  eps={eps:.0e}: 0/{N} diverge (0.0%)")

# ============================================================
# PART 3: What happens at the first divergence point (eps=1e-2)?
# ============================================================
print("\n" + "=" * 70)
print("PART 3: Token flip analysis at first divergence (eps=1e-2)")
print("=" * 70)

eps = 1e-2
flip_data = []
for j, (margins, top2, idx) in enumerate(zip(all_margins, all_top2_ids, all_idxs)):
    flip_positions = np.where(margins < eps)[0]
    if len(flip_positions) > 0:
        pos = flip_positions[0]
        old_tok = top2[pos, 0]
        new_tok = top2[pos, 1]
        old_str = tokenizer.decode([old_tok])
        new_str = tokenizer.decode([new_tok])
        margin_val = margins[pos]
        cat = prompts[idx]["category"]
        flip_data.append({
            "idx": idx, "pos": pos, "margin": margin_val,
            "old": old_str, "new": new_str,
            "old_id": int(old_tok), "new_id": int(new_tok),
            "cat": cat,
            "in_exp": idx in exp_set,
            "in_ctrl": idx in ctrl_set,
        })

print(f"Total flips at eps=1e-2: {len(flip_data)}")

# Classify flip types
def classify_flip(old, new):
    old_s, new_s = old.strip(), new.strip()
    if old_s.isdigit() and new_s.isdigit():
        return "digit-digit"
    elif old_s in "+-*/=" and new_s in "+-*/=":
        return "operator"
    elif old in ["\n", "\r", "\r\n", " "] and new in ["\n", "\r", "\r\n", " "]:
        return "whitespace"
    elif old.strip() and new.strip() and old_s[0].isupper() and new_s[0].isupper():
        return "sentence-start"
    elif old in [".", ",", "!", "?", ";", ":"] or new in [".", ",", "!", "?", ";", ":"]:
        return "punctuation"
    else:
        return "word"

type_counter = Counter()
for fd in flip_data:
    ft = classify_flip(fd["old"], fd["new"])
    type_counter[ft] += 1

print(f"\nFlip type distribution:")
for ft, count in type_counter.most_common():
    print(f"  {ft:15s}: {count:4d} ({count/len(flip_data)*100:.1f}%)")

# Position distribution of first flips
print(f"\nFirst flip position distribution:")
positions = [fd["pos"] for fd in flip_data]
for lo, hi in [(0, 2), (3, 7), (8, 15), (16, 31), (32, 63), (64, 127), (128, 255)]:
    count = sum(1 for p in positions if lo <= p <= hi)
    if count > 0:
        print(f"  pos {lo:3d}-{hi:3d}: {count:4d} ({count/len(positions)*100:.1f}%)")

# Per-category
print(f"\nFirst flip by category:")
cat_flips = defaultdict(int)
cat_total = Counter(prompts[idx]["category"] for idx in all_idxs)
for fd in flip_data:
    cat_flips[fd["cat"]] += 1
for cat in sorted(cat_total.keys()):
    n_cat = cat_total[cat]
    n_flip = cat_flips[cat]
    print(f"  {cat:15s}: {n_flip}/{n_cat} ({n_flip/n_cat*100:.1f}%)")

# Show concrete examples
print(f"\nTop 20 concrete first-flips (smallest margins):")
flip_data.sort(key=lambda x: x["margin"])
for i, fd in enumerate(flip_data[:20]):
    text = prompts[fd["idx"]]["text"][:50]
    old_repr = repr(fd["old"])
    new_repr = repr(fd["new"])
    set_tag = "EXP" if fd["in_exp"] else ("CTRL" if fd["in_ctrl"] else "---")
    print(f"  #{i+1:2d} idx={fd['idx']:4d} pos={fd['pos']:3d} margin={fd['margin']:.6f} "
          f"[{fd['cat'][:6]:6s}|{set_tag}] {old_repr} -> {new_repr}")
    print(f"      {text}")

# ============================================================
# PART 4: What fraction of EXPERIMENTAL vs CONTROL prompts diverge?
# ============================================================
print("\n" + "=" * 70)
print("PART 4: Experimental vs Control divergence rates")
print("=" * 70)

print(f"{'eps':>8s} {'Exp div':>10s} {'Ctrl div':>10s} {'Ratio':>8s}")
print("-" * 40)
for eps in [1e-5, 1e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1]:
    n_exp, n_ctrl = 0, 0
    for j, (margins, idx) in enumerate(zip(all_margins, all_idxs)):
        if (margins < eps).any():
            if idx in exp_set:
                n_exp += 1
            if idx in ctrl_set:
                n_ctrl += 1
    ratio = n_exp / max(n_ctrl, 1)
    print(f"{eps:8.0e} {n_exp:5d}/200  {n_ctrl:5d}/200  {ratio:8.1f}x")

# ============================================================
# PART 5: "Semantic severity" of first divergence
# ============================================================
print("\n" + "=" * 70)
print("PART 5: Semantic severity at various noise levels")
print("=" * 70)

# For each noise level, classify the first flip as:
# cosmetic (whitespace, punctuation), structural (sentence-start), or content (word, digit)
for eps in [1e-3, 1e-2, 1e-1]:
    severity = Counter()
    for j, (margins, top2, idx) in enumerate(zip(all_margins, all_top2_ids, all_idxs)):
        flip_positions = np.where(margins < eps)[0]
        if len(flip_positions) > 0:
            pos = flip_positions[0]
            old_str = tokenizer.decode([top2[pos, 0]])
            new_str = tokenizer.decode([top2[pos, 1]])
            ft = classify_flip(old_str, new_str)
            if ft in ("whitespace",):
                severity["cosmetic"] += 1
            elif ft in ("punctuation",):
                severity["punctuation"] += 1
            elif ft in ("sentence-start",):
                severity["structural"] += 1
            else:
                severity["content"] += 1
    total = sum(severity.values())
    print(f"\n  eps={eps:.0e}: {total} diverging prompts")
    for sev in ["content", "structural", "punctuation", "cosmetic"]:
        n = severity[sev]
        print(f"    {sev:12s}: {n:4d} ({n/max(total,1)*100:.1f}%)")
