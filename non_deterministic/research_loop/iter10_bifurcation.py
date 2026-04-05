"""
Iteration 10: Verify the Bifurcation Claim

Unverified claim: "by position ~32, the trajectory is determined."
Questions:
1. How well does early fragility (pos 0-31) predict late fragility (pos 32-255)?
2. Is there a single critical position where the prediction becomes sharp?
3. Can we find prompts that switch attractor state mid-generation?
4. What does the transition look like — is it a cliff or a ramp?
"""

import torch
import numpy as np
import json
import os
from collections import Counter
from pathlib import Path

LOGIT_DIR = "C:/Projects/math/non_deterministic/experiment_0/logit_traces"
PROMPTS_PATH = "C:/Projects/math/non_deterministic/prompts/candidates_2000.jsonl"

prompts = {}
with open(PROMPTS_PATH, encoding="utf-8") as f:
    for line in f:
        p = json.loads(line)
        prompts[p["idx"]] = p

# ============================================================
# Load all margin data in windows
# ============================================================
print("Loading margin data...", flush=True)

data = []
for i in range(1999):
    path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
    if not os.path.exists(path):
        continue
    d = torch.load(path, map_location="cpu", weights_only=True)
    margins = (d["top5_values"][:, 0] - d["top5_values"][:, 1]).numpy()
    n = len(margins)

    # Compute fragility in various windows
    windows = {}
    for (lo, hi, label) in [
        (0, 8, "w0_8"), (8, 16, "w8_16"), (16, 24, "w16_24"),
        (24, 32, "w24_32"), (32, 48, "w32_48"), (48, 64, "w48_64"),
        (64, 96, "w64_96"), (96, 128, "w96_128"), (128, 192, "w128_192"),
        (192, 256, "w192_256"),
        # Coarse windows
        (0, 16, "early16"), (16, 32, "mid16_32"),
        (0, 32, "early32"), (32, 256, "late32"),
        (0, 64, "early64"), (64, 256, "late64"),
    ]:
        seg = margins[lo:min(hi, n)]
        if len(seg) > 0:
            windows[label] = float((seg < 1.0).mean())
        else:
            windows[label] = 0.0

    # Also track margin means (not just fragility)
    windows["mean_margin_0_32"] = float(margins[:min(32, n)].mean())
    windows["mean_margin_32_256"] = float(margins[min(32, n):].mean()) if n > 32 else 0.0
    windows["min_margin"] = float(margins.min())
    windows["overall_frag"] = float((margins < 1.0).mean())

    data.append({
        "idx": i,
        "category": prompts[i]["category"],
        **windows,
    })

N = len(data)
print(f"Loaded {N} prompts", flush=True)

# ============================================================
# PART 1: How well does early fragility predict late fragility?
# ============================================================
print("\n" + "=" * 70)
print("PART 1: Prediction of late fragility from early fragility")
print("=" * 70)

# Try various cutpoints
cutpoints = [8, 16, 24, 32, 48, 64, 96]
for cut in cutpoints:
    early_key = f"early{cut}" if f"early{cut}" in data[0] else None
    late_key = f"late{cut}" if f"late{cut}" in data[0] else None

    if not early_key or not late_key:
        # Compute manually
        early_vals = []
        late_vals = []
        for d in data:
            path = os.path.join(LOGIT_DIR, f"prompt_{d['idx']:04d}.pt")
            dd = torch.load(path, map_location="cpu", weights_only=True)
            margins = (dd["top5_values"][:, 0] - dd["top5_values"][:, 1]).numpy()
            n = len(margins)
            early_vals.append(float((margins[:min(cut, n)] < 1.0).mean()))
            late_vals.append(float((margins[min(cut, n):] < 1.0).mean()) if n > cut else 0.0)
    else:
        early_vals = [d[early_key] for d in data]
        late_vals = [d[late_key] for d in data]

    early_arr = np.array(early_vals)
    late_arr = np.array(late_vals)
    r = np.corrcoef(early_arr, late_arr)[0, 1]
    print(f"  Cutpoint {cut:3d}: r(early, late) = {r:.4f}")

# ============================================================
# PART 2: Information gain at each position
# ============================================================
print("\n" + "=" * 70)
print("PART 2: Cumulative predictive power — at what position can we predict overall fragility?")
print("=" * 70)

overall_frags = np.array([d["overall_frag"] for d in data])

# For each prefix length, compute correlation of prefix fragility with overall
print(f"\n  Prefix length -> r(prefix_frag, overall_frag):")
for prefix in [1, 2, 3, 4, 8, 12, 16, 20, 24, 28, 32, 48, 64, 96, 128]:
    prefix_frags = []
    for d in data:
        path = os.path.join(LOGIT_DIR, f"prompt_{d['idx']:04d}.pt")
        dd = torch.load(path, map_location="cpu", weights_only=True)
        margins = (dd["top5_values"][:, 0] - dd["top5_values"][:, 1]).numpy()
        pf = float((margins[:min(prefix, len(margins))] < 1.0).mean())
        prefix_frags.append(pf)
    prefix_arr = np.array(prefix_frags)
    r = np.corrcoef(prefix_arr, overall_frags)[0, 1]
    # Also: can the prefix median-split classify top/bottom fragility quartiles?
    q75 = np.percentile(overall_frags, 75)
    q25 = np.percentile(overall_frags, 25)
    high_mask = overall_frags >= q75
    low_mask = overall_frags <= q25
    if high_mask.sum() > 0 and low_mask.sum() > 0:
        threshold = np.median(prefix_arr)
        pred_high = prefix_arr >= threshold
        tp = (pred_high & high_mask).sum()
        tn = (~pred_high & low_mask).sum()
        acc = (tp + tn) / (high_mask.sum() + low_mask.sum())
    else:
        acc = 0
    print(f"  prefix={prefix:3d}: r={r:.4f}  quartile_acc={acc:.3f}")

# ============================================================
# PART 3: Find prompts that SWITCH attractor states
# ============================================================
print("\n" + "=" * 70)
print("PART 3: State-switching prompts")
print("=" * 70)

# Define: "state switch" = early fragility and late fragility are in different quartiles
early32 = np.array([d["early32"] for d in data])
late32 = np.array([d["late32"] for d in data])

# Category 1: Calm start, fragile finish (early < 0.3, late > 0.15)
calm_to_fragile = [(i, d) for i, d in enumerate(data)
                   if d["early32"] < 0.30 and d["late32"] > 0.15]
# Category 2: Fragile start, calm finish (early > 0.5, late < 0.05)
fragile_to_calm = [(i, d) for i, d in enumerate(data)
                   if d["early32"] > 0.50 and d["late32"] < 0.05]

print(f"Calm->Fragile (early32 < 0.3, late32 > 0.15): {len(calm_to_fragile)}")
print(f"Fragile->Calm (early32 > 0.5, late32 < 0.05): {len(fragile_to_calm)}")

# Show examples of calm->fragile (the interesting state switches)
if calm_to_fragile:
    calm_to_fragile.sort(key=lambda x: x[1]["late32"] - x[1]["early32"], reverse=True)
    print(f"\nTop 10 Calm->Fragile switches (biggest late-early gap):")
    for rank, (i, d) in enumerate(calm_to_fragile[:10]):
        print(f"  #{rank+1} idx={d['idx']:4d} [{d['category']:12s}] "
              f"early32={d['early32']:.3f} late32={d['late32']:.3f} "
              f"gap={d['late32']-d['early32']:+.3f}")
        print(f"    Prompt: {prompts[d['idx']]['text'][:70]}")

if fragile_to_calm:
    fragile_to_calm.sort(key=lambda x: x[1]["early32"] - x[1]["late32"], reverse=True)
    print(f"\nTop 10 Fragile->Calm switches (biggest early-late gap):")
    for rank, (i, d) in enumerate(fragile_to_calm[:10]):
        print(f"  #{rank+1} idx={d['idx']:4d} [{d['category']:12s}] "
              f"early32={d['early32']:.3f} late32={d['late32']:.3f} "
              f"gap={d['early32']-d['late32']:+.3f}")
        print(f"    Prompt: {prompts[d['idx']]['text'][:70]}")

# ============================================================
# PART 4: Fine-grained window trajectory shapes
# ============================================================
print("\n" + "=" * 70)
print("PART 4: Trajectory shape classification")
print("=" * 70)

# Classify each prompt's trajectory as:
# A) Monotonic decrease (pos 2 is peak, then all downhill)
# B) U-shape (early high, mid low, late high)
# C) Monotonic increase (rare — getting MORE fragile)
# D) Flat (roughly constant)

win_keys = ["w0_8", "w8_16", "w16_24", "w24_32", "w32_48", "w48_64",
            "w64_96", "w96_128", "w128_192", "w192_256"]

shapes = Counter()
for d in data:
    trajectory = [d[k] for k in win_keys]

    # Classify
    # Get early (avg of first 2), mid (avg of middle), late (avg of last 2)
    early = np.mean(trajectory[:2])
    mid = np.mean(trajectory[3:6])
    late = np.mean(trajectory[-2:])

    if early > mid + 0.05 and late > mid + 0.05:
        shape = "U-shape"
    elif early > mid > late and early - late > 0.05:
        shape = "Monotonic decrease"
    elif late > mid > early and late - early > 0.05:
        shape = "Monotonic increase"
    elif abs(early - late) < 0.03 and abs(early - mid) < 0.05:
        shape = "Flat"
    else:
        shape = "Irregular"
    shapes[shape] += 1

print(f"Trajectory shape distribution:")
for shape, count in shapes.most_common():
    print(f"  {shape:22s}: {count:5d} ({count/N*100:.1f}%)")

# ============================================================
# PART 5: The bifurcation sharpness — scatter plot data
# ============================================================
print("\n" + "=" * 70)
print("PART 5: Bifurcation sharpness at position 32")
print("=" * 70)

# How well does fragility in pos 0-31 separate into clusters?
# Look at the joint distribution of early32 and late32

# Bin early32 into deciles and show mean late32
early_sorted = sorted(data, key=lambda d: d["early32"])
decile_size = N // 10
print(f"\nMean late32 by early32 decile:")
print(f"{'Decile':>8s} {'Early32 range':>20s} {'Mean late32':>12s} {'Std':>8s}")
for dec in range(10):
    group = early_sorted[dec * decile_size:(dec + 1) * decile_size]
    early_range = f"{group[0]['early32']:.3f}-{group[-1]['early32']:.3f}"
    late_vals = [d["late32"] for d in group]
    print(f"{dec+1:8d} {early_range:>20s} {np.mean(late_vals):12.4f} {np.std(late_vals):8.4f}")

# Key test: is there a GAP in the early32 distribution?
early32_arr = np.array([d["early32"] for d in data])
print(f"\nEarly32 distribution:")
print(f"  Mean: {early32_arr.mean():.4f}")
print(f"  Std: {early32_arr.std():.4f}")
print(f"  Quantiles: P10={np.percentile(early32_arr, 10):.3f} "
      f"P25={np.percentile(early32_arr, 25):.3f} "
      f"P50={np.percentile(early32_arr, 50):.3f} "
      f"P75={np.percentile(early32_arr, 75):.3f} "
      f"P90={np.percentile(early32_arr, 90):.3f}")

# Histogram
bins = np.arange(0, 1.01, 0.05)
hist, _ = np.histogram(early32_arr, bins=bins)
print(f"\nEarly32 histogram:")
for i in range(len(hist)):
    bar = "#" * (hist[i] * 60 // max(hist))
    print(f"  {bins[i]:.2f}-{bins[i+1]:.2f}: {hist[i]:4d} {bar}")
