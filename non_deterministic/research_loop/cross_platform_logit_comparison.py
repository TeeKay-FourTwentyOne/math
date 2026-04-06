"""
Definitive Cross-Platform Logit Comparison: H_nvidia vs H_mac
Both files now have aligned format with top5_values and top5_indices.
"""

import json
import numpy as np
from collections import Counter, defaultdict

def load_run(path):
    """Load a run, strip input tokens from generated_ids, return indexed dict."""
    data = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            d = json.loads(line)
            offset = d["n_input_tokens"]
            d["gen_only"] = d["generated_ids"][offset:]
            data[d["prompt_idx"]] = d
    return data

nvidia = load_run("C:/Projects/math/non_deterministic/experiment_1/outputs/h_nvidia_run1.jsonl")
mac = load_run("C:/Projects/math/non_deterministic/experiment_1/outputs/h_mac_run1.jsonl")

# Load exp/ctrl sets
exp_idxs = set()
ctrl_idxs = set()
with open("C:/Projects/math/non_deterministic/prompts/experimental_200.jsonl", encoding="utf-8") as f:
    for line in f:
        exp_idxs.add(json.loads(line)["prompt_idx"])
with open("C:/Projects/math/non_deterministic/prompts/control_200.jsonl", encoding="utf-8") as f:
    for line in f:
        ctrl_idxs.add(json.loads(line)["prompt_idx"])

idxs = sorted(set(nvidia) & set(mac))

# ============================================================
print("=" * 70)
print("CROSS-PLATFORM COMPARISON: H_nvidia vs H_mac (aligned format)")
print("=" * 70)

# PART 1: Token-level comparison
print("\n--- PART 1: Token-level ---")
token_match = 0
token_div = []
for idx in idxs:
    n_gen = nvidia[idx]["gen_only"]
    m_gen = mac[idx]["gen_only"]
    if n_gen == m_gen:
        token_match += 1
    else:
        min_len = min(len(n_gen), len(m_gen))
        fd = next((j for j in range(min_len) if n_gen[j] != m_gen[j]), min_len)
        token_div.append({"idx": idx, "first_div": fd})

print(f"Prompts: {len(idxs)}")
print(f"Token exact match: {token_match}/{len(idxs)} ({token_match/len(idxs)*100:.1f}%)")
print(f"Token diverged: {len(token_div)}/{len(idxs)}")

# PART 2: Logit-level comparison (the main event)
print("\n--- PART 2: Logit-level (per-step top-1 value) ---")

all_diffs = []          # |nvidia_logit - mac_logit| at rank-1
all_margins_n = []      # nvidia margin (rank1 - rank2)
all_margins_m = []      # mac margin
all_margin_diffs = []   # |nvidia_margin - mac_margin|
all_ratios = []         # noise / min(margin)
closest_calls = []
per_prompt_max_diff = []
per_prompt_max_ratio = []

# Also track rank-agreement: do nvidia and mac agree on which token is rank-1?
rank1_mismatches = 0
rank1_total = 0

for idx in idxs:
    n_vals = np.array(nvidia[idx]["top5_values"])  # (256, 5)
    m_vals = np.array(mac[idx]["top5_values"])      # (256, 5)
    n_idxs_arr = np.array(nvidia[idx]["top5_indices"])  # (256, 5)
    m_idxs_arr = np.array(mac[idx]["top5_indices"])

    n_steps = min(len(n_vals), len(m_vals))
    prompt_max_diff = 0
    prompt_max_ratio = 0

    for s in range(n_steps):
        # Rank-1 logit difference
        diff = abs(n_vals[s, 0] - m_vals[s, 0])
        all_diffs.append(diff)

        # Margins
        n_margin = n_vals[s, 0] - n_vals[s, 1]
        m_margin = m_vals[s, 0] - m_vals[s, 1]
        all_margins_n.append(n_margin)
        all_margins_m.append(m_margin)
        all_margin_diffs.append(abs(n_margin - m_margin))

        # Noise-to-margin ratio
        min_margin = min(n_margin, m_margin)
        if min_margin > 0:
            ratio = diff / min_margin
            all_ratios.append(ratio)
            if ratio > prompt_max_ratio:
                prompt_max_ratio = ratio

            if ratio > 0.1:
                closest_calls.append({
                    "idx": idx, "pos": s, "noise": diff,
                    "n_margin": n_margin, "m_margin": m_margin,
                    "ratio": ratio,
                    "cat": nvidia[idx]["category"],
                    "text": nvidia[idx]["text"][:50],
                    "in_exp": idx in exp_idxs,
                })

        if diff > prompt_max_diff:
            prompt_max_diff = diff

        # Rank-1 token agreement
        rank1_total += 1
        if n_idxs_arr[s, 0] != m_idxs_arr[s, 0]:
            rank1_mismatches += 1

    per_prompt_max_diff.append(prompt_max_diff)
    per_prompt_max_ratio.append(prompt_max_ratio)

diffs = np.array(all_diffs)
margins_n = np.array(all_margins_n)
margins_m = np.array(all_margins_m)
margin_diffs = np.array(all_margin_diffs)
ratios = np.array(all_ratios)

print(f"\nPositions analyzed: {len(diffs):,}")
print(f"\nLogit noise (|nvidia_rank1 - mac_rank1|):")
print(f"  Mean:   {diffs.mean():.6f}")
print(f"  Median: {np.median(diffs):.6f}")
print(f"  P95:    {np.percentile(diffs, 95):.6f}")
print(f"  P99:    {np.percentile(diffs, 99):.6f}")
print(f"  Max:    {diffs.max():.6f}")
print(f"  Zeros:  {(diffs == 0).sum():,} ({(diffs == 0).mean()*100:.2f}%)")

print(f"\nMargin difference (|nvidia_margin - mac_margin|):")
print(f"  Mean:   {margin_diffs.mean():.6f}")
print(f"  Median: {np.median(margin_diffs):.6f}")
print(f"  Max:    {margin_diffs.max():.6f}")

print(f"\nNoise / margin ratio:")
print(f"  Mean:   {ratios.mean():.6f}")
print(f"  Median: {np.median(ratios):.6f}")
print(f"  P99:    {np.percentile(ratios, 99):.6f}")
print(f"  Max:    {ratios.max():.6f}")
print(f"  > 0.5:  {(ratios > 0.5).sum()}")
print(f"  > 0.1:  {(ratios > 0.1).sum()}")
print(f"  > 0.01: {(ratios > 0.01).sum()}")

print(f"\nRank-1 token agreement: {rank1_total - rank1_mismatches}/{rank1_total} "
      f"({(1 - rank1_mismatches/rank1_total)*100:.4f}%)")
print(f"Rank-1 mismatches: {rank1_mismatches}")

# PART 3: Per-prompt summary
print("\n--- PART 3: Per-prompt noise summary ---")
pmd = np.array(per_prompt_max_diff)
pmr = np.array(per_prompt_max_ratio)
print(f"Max logit noise per prompt:")
print(f"  Mean:   {pmd.mean():.6f}")
print(f"  Median: {np.median(pmd):.6f}")
print(f"  Max:    {pmd.max():.6f}")
print(f"Max noise/margin ratio per prompt:")
print(f"  Mean:   {pmr.mean():.6f}")
print(f"  Median: {np.median(pmr):.6f}")
print(f"  Max:    {pmr.max():.6f}")
print(f"  Prompts with ratio > 0.5: {(pmr > 0.5).sum()}")
print(f"  Prompts with ratio > 0.1: {(pmr > 0.1).sum()}")

# PART 4: Exp vs Ctrl noise levels
print("\n--- PART 4: Experimental vs Control noise ---")
for label, idx_set in [("Experimental", exp_idxs), ("Control", ctrl_idxs)]:
    group_diffs = []
    group_ratios = []
    for idx in idxs:
        if idx not in idx_set:
            continue
        n_vals = np.array(nvidia[idx]["top5_values"])
        m_vals = np.array(mac[idx]["top5_values"])
        for s in range(min(len(n_vals), len(m_vals))):
            diff = abs(n_vals[s, 0] - m_vals[s, 0])
            group_diffs.append(diff)
            margin = min(n_vals[s, 0] - n_vals[s, 1], m_vals[s, 0] - m_vals[s, 1])
            if margin > 0:
                group_ratios.append(diff / margin)
    gd = np.array(group_diffs)
    gr = np.array(group_ratios)
    print(f"  {label}:")
    print(f"    Logit noise mean: {gd.mean():.6f}, max: {gd.max():.6f}")
    print(f"    Ratio mean: {gr.mean():.6f}, max: {gr.max():.6f}, >0.1: {(gr > 0.1).sum()}")

# PART 5: Closest calls
print("\n--- PART 5: Closest calls (noise > 10% of margin) ---")
closest_calls.sort(key=lambda c: -c["ratio"])
for i, c in enumerate(closest_calls[:20]):
    tag = "EXP" if c["in_exp"] else "CTRL"
    print(f"  #{i+1} idx={c['idx']:4d} pos={c['pos']:3d} noise={c['noise']:.6f} "
          f"n_margin={c['n_margin']:.6f} m_margin={c['m_margin']:.6f} "
          f"ratio={c['ratio']:.4f} [{c['cat'][:6]}|{tag}]")
    print(f"       {c['text']}")

# PART 6: Position-dependent noise
print("\n--- PART 6: Noise by generation position ---")
pos_diffs = defaultdict(list)
for idx in idxs:
    n_vals = np.array(nvidia[idx]["top5_values"])
    m_vals = np.array(mac[idx]["top5_values"])
    for s in range(min(len(n_vals), len(m_vals))):
        pos_diffs[s].append(abs(n_vals[s, 0] - m_vals[s, 0]))

print(f"{'Pos range':>12s} {'Mean noise':>12s} {'Max noise':>12s}")
for lo, hi in [(0, 7), (8, 15), (16, 31), (32, 63), (64, 127), (128, 191), (192, 255)]:
    vals = []
    for p in range(lo, hi + 1):
        if p in pos_diffs:
            vals.extend(pos_diffs[p])
    if vals:
        arr = np.array(vals)
        print(f"  {lo:3d}-{hi:3d}:    {arr.mean():.6f}    {arr.max():.6f}")

# PART 7: Do margins ever DISAGREE on which token is better?
print("\n--- PART 7: Margin sign check ---")
# Cases where nvidia_rank1 != mac_rank1 in top5_indices
# (This would mean the platforms choose different tokens despite both using greedy)
print(f"Positions where rank-1 token ID differs: {rank1_mismatches}/{rank1_total}")
if rank1_mismatches > 0:
    print("WARNING: Some positions have different greedy-selected tokens!")
else:
    print("All positions agree on the greedy-selected token (rank-1 ID matches).")
