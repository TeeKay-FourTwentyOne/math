"""
Verify the "position 16" claim from STATE.md.

STATE.md claims: "49.6% of ALL prompts are fragile at position 16"
Our PHASE 9 found: only 7.3% fragility rate at position 16.

Possible interpretations:
1. "49.6% of prompts have at least one fragile position in 0-31" (cumulative)
2. "Position 16" is the most common first-fragile-position
3. Something else

Let's figure out what the original claim means.
"""

import torch
import numpy as np
import os
from collections import Counter

LOGIT_DIR = "C:/Projects/math/non_deterministic/experiment_0/logit_traces"
MARGIN_THRESHOLD = 0.1

# For each prompt, find:
# 1. Whether it has any fragile position in 0-31
# 2. Whether position 16 specifically is fragile
# 3. The first fragile position
# 4. Whether it has a fragile position at position 15 or 16 (±1)

has_fragile_0_31 = 0
has_fragile_at_16 = 0
has_fragile_at_15_or_16 = 0
has_fragile_at_14_to_17 = 0
first_fragile_pos = Counter()
any_fragile_0_15 = 0
any_fragile_16_31 = 0

# More detailed: for each prompt check if its FIRST fragile position is around 16
first_fragile_positions = []

for i in range(1999):
    path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
    if not os.path.exists(path):
        continue
    data = torch.load(path, map_location='cpu', weights_only=False)
    top5_vals = data['top5_values']
    margins = top5_vals[:, 0] - top5_vals[:, 1]

    fragile_positions = (margins < MARGIN_THRESHOLD).nonzero().flatten().tolist()

    if fragile_positions:
        first_pos = fragile_positions[0]
        first_fragile_positions.append(first_pos)
        first_fragile_pos[first_pos] += 1

    if any(p <= 31 for p in fragile_positions):
        has_fragile_0_31 += 1
    if 16 in fragile_positions:
        has_fragile_at_16 += 1
    if 15 in fragile_positions or 16 in fragile_positions:
        has_fragile_at_15_or_16 += 1
    if any(14 <= p <= 17 for p in fragile_positions):
        has_fragile_at_14_to_17 += 1
    if any(p <= 15 for p in fragile_positions):
        any_fragile_0_15 += 1
    if any(16 <= p <= 31 for p in fragile_positions):
        any_fragile_16_31 += 1

N = 1999
print("Position-16 claim verification:")
print(f"  Prompts with fragile at exactly pos 16: {has_fragile_at_16}/{N} = {100*has_fragile_at_16/N:.1f}%")
print(f"  Prompts with fragile at pos 15 or 16:   {has_fragile_at_15_or_16}/{N} = {100*has_fragile_at_15_or_16/N:.1f}%")
print(f"  Prompts with fragile at pos 14-17:      {has_fragile_at_14_to_17}/{N} = {100*has_fragile_at_14_to_17/N:.1f}%")
print(f"  Prompts with ANY fragile in 0-15:       {any_fragile_0_15}/{N} = {100*any_fragile_0_15/N:.1f}%")
print(f"  Prompts with ANY fragile in 16-31:      {any_fragile_16_31}/{N} = {100*any_fragile_16_31/N:.1f}%")
print(f"  Prompts with ANY fragile in 0-31:       {has_fragile_0_31}/{N} = {100*has_fragile_0_31/N:.1f}%")

# Most common first fragile position
print(f"\nFirst fragile position distribution (top 20):")
for pos, count in first_fragile_pos.most_common(20):
    print(f"  pos {pos:3d}: {count:4d} ({100*count/N:.1f}%)")

# What percentage have their first fragile in bins?
bins = [(0, 7), (8, 15), (16, 23), (24, 31), (32, 63), (64, 127), (128, 255)]
print(f"\nFirst fragile position by bin:")
for lo, hi in bins:
    count = sum(1 for p in first_fragile_positions if lo <= p <= hi)
    print(f"  {lo:3d}-{hi:3d}: {count:4d} ({100*count/N:.1f}%)")

# Now check: maybe the 49.6% figure is about the position WHERE most prompts
# "make their first decision" (including non-fragile)
# Let's look at margin distribution at position 16 specifically
print(f"\nMargin distribution at position 16:")
margins_at_16 = []
for i in range(1999):
    path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
    if not os.path.exists(path):
        continue
    data = torch.load(path, map_location='cpu', weights_only=False)
    top5_vals = data['top5_values']
    margins = top5_vals[:, 0] - top5_vals[:, 1]
    if len(margins) > 16:
        margins_at_16.append(margins[16].item())

m16 = np.array(margins_at_16)
print(f"  N = {len(m16)}")
print(f"  Mean margin: {m16.mean():.4f}")
print(f"  Median margin: {np.median(m16):.4f}")
print(f"  Fraction < 0.1: {(m16 < 0.1).mean():.3f}")
print(f"  Fraction < 0.5: {(m16 < 0.5).mean():.3f}")
print(f"  Fraction < 1.0: {(m16 < 1.0).mean():.3f}")
print(f"  Fraction < 2.0: {(m16 < 2.0).mean():.3f}")

# Compare with other positions
print(f"\n\nMean margin by position (0-31):")
for pos in range(32):
    margins_at_pos = []
    for i in range(1999):
        path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
        if not os.path.exists(path):
            continue
        data = torch.load(path, map_location='cpu', weights_only=False)
        top5_vals = data['top5_values']
        margins = top5_vals[:, 0] - top5_vals[:, 1]
        if len(margins) > pos:
            margins_at_pos.append(margins[pos].item())
    arr = np.array(margins_at_pos)
    fragile_rate = (arr < 0.1).mean()
    print(f"  pos {pos:2d}: mean={arr.mean():.3f} median={np.median(arr):.3f} fragile_rate={fragile_rate:.3f}")
