"""
Iteration 2: Definitive position-by-position fragility analysis.

The seed STATE.md claims "49.6% of ALL prompts are fragile at position 16."
A previous script checked this with threshold 0.1 and got 7.3%.
But the experiment's fragility threshold is 1.0. Let's check BOTH thresholds
and map the full position curve.

Key questions:
1. Per-position fragility rate at threshold=1.0 — does pos 16 stand out?
2. Per-position fragility rate at threshold=0.1 — same question
3. Does the pattern hold per-category?
4. What TOKENS appear at position 16? Is this a tokenization artifact (newline, etc.)?
5. What does the fine-grained ramp look like in positions 0-63?
"""

import torch
import numpy as np
import json
import os
from collections import Counter, defaultdict

LOGIT_DIR = "C:/Projects/math/non_deterministic/experiment_0/logit_traces"
PROMPTS_PATH = "C:/Projects/math/non_deterministic/prompts/candidates_2000.jsonl"

# Load prompt metadata
prompts = {}
with open(PROMPTS_PATH) as f:
    for line in f:
        p = json.loads(line)
        prompts[p["idx"]] = p

# ============================================================
# PHASE 1: Per-position fragility at BOTH thresholds
# ============================================================
print("=" * 70)
print("PHASE 1: Per-position fragility rates")
print("=" * 70)

N_STEPS = 256
# Accumulators: per position, count of fragile and total
frag_count_1_0 = np.zeros(N_STEPS)  # threshold 1.0
frag_count_0_1 = np.zeros(N_STEPS)  # threshold 0.1
total_count = np.zeros(N_STEPS)

# Also per-category
categories = ["factual", "open_ended", "reasoning", "continuation", "conversational"]
cat_frag_1_0 = {c: np.zeros(N_STEPS) for c in categories}
cat_total = {c: np.zeros(N_STEPS) for c in categories}

# Also collect the actual margin values at position 16 for histogram
margins_at_16 = []
tokens_at_16 = []

for i in range(1999):
    path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
    if not os.path.exists(path):
        continue
    data = torch.load(path, map_location="cpu", weights_only=True)
    top5_vals = data["top5_values"]
    top5_idxs = data["top5_indices"]
    margins = (top5_vals[:, 0] - top5_vals[:, 1]).numpy()
    n = len(margins)
    cat = prompts[i]["category"]

    for pos in range(min(n, N_STEPS)):
        total_count[pos] += 1
        cat_total[cat][pos] += 1
        if margins[pos] < 1.0:
            frag_count_1_0[pos] += 1
            cat_frag_1_0[cat][pos] += 1
        if margins[pos] < 0.1:
            frag_count_0_1[pos] += 1

    if n > 16:
        margins_at_16.append(margins[16])
        tokens_at_16.append(top5_idxs[16, 0].item())

# Compute rates
rate_1_0 = np.where(total_count > 0, frag_count_1_0 / total_count, 0)
rate_0_1 = np.where(total_count > 0, frag_count_0_1 / total_count, 0)

print("\nPer-position fragility (threshold=1.0), positions 0-63:")
print(f"{'Pos':>4s} {'Rate':>7s}  {'Pos':>4s} {'Rate':>7s}  {'Pos':>4s} {'Rate':>7s}  {'Pos':>4s} {'Rate':>7s}")
for row in range(16):
    parts = []
    for col in range(4):
        pos = row + col * 16
        parts.append(f"{pos:4d} {rate_1_0[pos]:7.3f}")
    print("  ".join(parts))

# Find the peak position
peak_pos_1_0 = np.argmax(rate_1_0)
peak_rate_1_0 = rate_1_0[peak_pos_1_0]
print(f"\nPeak position (threshold=1.0): pos {peak_pos_1_0} at {peak_rate_1_0:.3f}")

peak_pos_0_1 = np.argmax(rate_0_1)
peak_rate_0_1 = rate_0_1[peak_pos_0_1]
print(f"Peak position (threshold=0.1): pos {peak_pos_0_1} at {peak_rate_0_1:.3f}")

# Position 16 specifically
print(f"\nPosition 16 specifically:")
print(f"  Threshold 1.0: {rate_1_0[16]:.4f} ({int(frag_count_1_0[16])}/{int(total_count[16])} prompts)")
print(f"  Threshold 0.1: {rate_0_1[16]:.4f} ({int(frag_count_0_1[16])}/{int(total_count[16])} prompts)")

# Binned summary
print(f"\nBinned fragility rates (threshold=1.0):")
bins = [(0, 7), (8, 15), (16, 23), (24, 31), (32, 63), (64, 127), (128, 191), (192, 255)]
for lo, hi in bins:
    rate = rate_1_0[lo:hi + 1].mean()
    print(f"  pos {lo:3d}-{hi:3d}: {rate:.4f}")

# ============================================================
# PHASE 2: Per-category breakdown at position 16
# ============================================================
print("\n" + "=" * 70)
print("PHASE 2: Per-category fragility at position 16 (threshold=1.0)")
print("=" * 70)

for cat in categories:
    if cat_total[cat][16] > 0:
        rate = cat_frag_1_0[cat][16] / cat_total[cat][16]
        print(f"  {cat:15s}: {rate:.4f} ({int(cat_frag_1_0[cat][16])}/{int(cat_total[cat][16])})")

# ============================================================
# PHASE 3: What tokens appear at position 16?
# ============================================================
print("\n" + "=" * 70)
print("PHASE 3: Token analysis at position 16")
print("=" * 70)

from transformers import AutoTokenizer
tokenizer = AutoTokenizer.from_pretrained("EleutherAI/pythia-410m-deduped")

# Decode the most common tokens at position 16
tok_counter = Counter(tokens_at_16)
print(f"Total prompts with pos 16: {len(tokens_at_16)}")
print(f"\nTop 20 tokens at position 16:")
for tok_id, count in tok_counter.most_common(20):
    tok_str = tokenizer.decode([tok_id])
    frac = count / len(tokens_at_16)
    print(f"  id={tok_id:5d}  '{tok_str}'  count={count:4d} ({frac:.3f})")

# Check if newline is dominant
newline_ids = [tokenizer.encode("\n")[-1]]
newline_count = sum(tok_counter.get(nid, 0) for nid in newline_ids)
print(f"\nNewline at pos 16: {newline_count}/{len(tokens_at_16)} ({newline_count/len(tokens_at_16):.3f})")

# ============================================================
# PHASE 4: Is position 0 special? (the very first generated token)
# ============================================================
print("\n" + "=" * 70)
print("PHASE 4: Position 0 (first generated token) analysis")
print("=" * 70)

# Margins at position 0
margins_at_0 = []
tokens_at_0 = []
for i in range(1999):
    path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
    if not os.path.exists(path):
        continue
    data = torch.load(path, map_location="cpu", weights_only=True)
    top5_vals = data["top5_values"]
    top5_idxs = data["top5_indices"]
    margins = (top5_vals[:, 0] - top5_vals[:, 1]).numpy()
    margins_at_0.append(margins[0])
    tokens_at_0.append(top5_idxs[0, 0].item())

m0 = np.array(margins_at_0)
print(f"Position 0 fragility (thresh 1.0): {(m0 < 1.0).mean():.4f}")
print(f"Position 0 fragility (thresh 0.1): {(m0 < 0.1).mean():.4f}")
print(f"Position 0 mean margin: {m0.mean():.4f}")
print(f"Position 0 median margin: {np.median(m0):.4f}")

tok0_counter = Counter(tokens_at_0)
print(f"\nTop 10 tokens at position 0:")
for tok_id, count in tok0_counter.most_common(10):
    tok_str = tokenizer.decode([tok_id])
    frac = count / len(tokens_at_0)
    print(f"  id={tok_id:5d}  '{tok_str}'  count={count:4d} ({frac:.3f})")

# ============================================================
# PHASE 5: Monotonicity check — does fragility always decrease?
# ============================================================
print("\n" + "=" * 70)
print("PHASE 5: Monotonicity — do any prompts have INCREASING fragility?")
print("=" * 70)

# For each prompt, compute fragility in windows of 32 positions
increasing_prompts = []
for i in range(1999):
    path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
    if not os.path.exists(path):
        continue
    data = torch.load(path, map_location="cpu", weights_only=True)
    top5_vals = data["top5_values"]
    margins = (top5_vals[:, 0] - top5_vals[:, 1]).numpy()

    # Window fragility: 8 windows of 32 positions each
    windows = []
    for w in range(8):
        start = w * 32
        end = min(start + 32, len(margins))
        if end <= start:
            break
        wfrag = (margins[start:end] < 1.0).mean()
        windows.append(wfrag)

    # Check if any late window has HIGHER fragility than the first window
    if len(windows) >= 4:
        early = windows[0]
        for w_idx in range(4, len(windows)):
            if windows[w_idx] > early + 0.1 and windows[w_idx] > 0.15:
                increasing_prompts.append({
                    "idx": i,
                    "category": prompts[i]["category"],
                    "text": prompts[i]["text"][:80],
                    "windows": [round(w, 3) for w in windows],
                    "early": round(early, 3),
                    "late": round(windows[w_idx], 3),
                })
                break

print(f"Prompts with increasing fragility (late > early + 0.1): {len(increasing_prompts)}")
for p in increasing_prompts[:10]:
    print(f"  idx={p['idx']:4d} [{p['category']:12s}] early={p['early']} late={p['late']}")
    print(f"    windows={p['windows']}")
    print(f"    text: {p['text']}")

# ============================================================
# PHASE 6: Fine-grained curve — single-position fragility 0-63
# ============================================================
print("\n" + "=" * 70)
print("PHASE 6: Fine-grained position curve (threshold=1.0)")
print("=" * 70)

# Check for sharp transitions
print("Position-by-position fragility rate, positions 0-63:")
for pos in range(64):
    r = rate_1_0[pos]
    bar = "#" * int(r * 100)
    print(f"  pos {pos:2d}: {r:.4f}  {bar}")

# Is there a specific cliff?
max_drop = 0
max_drop_pos = 0
for pos in range(1, 64):
    drop = rate_1_0[pos - 1] - rate_1_0[pos]
    if drop > max_drop:
        max_drop = drop
        max_drop_pos = pos
print(f"\nLargest single-position drop: pos {max_drop_pos-1} -> {max_drop_pos} "
      f"({rate_1_0[max_drop_pos-1]:.4f} -> {rate_1_0[max_drop_pos]:.4f}, delta={max_drop:.4f})")

# ============================================================
# PHASE 7: Margin histogram at position 16
# ============================================================
print("\n" + "=" * 70)
print("PHASE 7: Margin distribution at position 16")
print("=" * 70)

m16 = np.array(margins_at_16)
print(f"N = {len(m16)}")
print(f"Mean margin: {m16.mean():.4f}")
print(f"Median margin: {np.median(m16):.4f}")
print(f"Std: {m16.std():.4f}")

thresholds = [0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0]
for t in thresholds:
    frac = (m16 < t).mean()
    print(f"  < {t:5.2f}: {frac:.4f} ({int((m16 < t).sum())} prompts)")
