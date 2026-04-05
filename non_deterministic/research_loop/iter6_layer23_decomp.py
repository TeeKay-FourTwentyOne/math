"""
Iteration 6: Layer 23 Residual Decomposition

Layer 23 has r=+0.97 correlation with fragility in its norms.
Its residual (output - layer_22_output) exceeds its own output norm.
But WHICH dimensions drive this? Is there a small "uncertainty subspace"?

Plan:
1. Compute layer-23 residual for all 50 activation prompts
2. Find dimensions with largest exp/ctrl difference
3. Test if a small set of dims can predict fragility
4. Check if this is a low-rank phenomenon (like dim 487 for BOS)
5. Look at how these dimensions relate to generated text positions
"""

import torch
import numpy as np
import json
import os
from pathlib import Path
from collections import defaultdict

ACT_DIR = "C:/Projects/math/non_deterministic/experiment_1/activations/h_nvidia"
SUBSET_PATH = "C:/Projects/math/non_deterministic/experiment_1/activations/activation_subset.jsonl"
LOGIT_DIR = "C:/Projects/math/non_deterministic/experiment_0/logit_traces"

# Load subset metadata
subset = []
with open(SUBSET_PATH, encoding="utf-8") as f:
    for line in f:
        subset.append(json.loads(line))

print(f"Loading {len(subset)} activation profiles...", flush=True)

# Load all activation data
exp_data = []  # experimental prompts
ctrl_data = []  # control prompts

for meta in subset:
    idx = meta["prompt_idx"]
    path = os.path.join(ACT_DIR, f"prompt_{idx:04d}.pt")
    data = torch.load(path, map_location="cpu", weights_only=False)

    # Compute layer-23 residual
    l23 = data["layer_23"]  # (seq_len, 1024)
    l22 = data["layer_22"]  # (seq_len, 1024)
    residual = l23 - l22     # what layer 23 adds

    input_len = data["input_length"]
    gen_start = input_len  # generated tokens start here

    # Also load fragility from logit traces
    logit_path = os.path.join(LOGIT_DIR, f"prompt_{idx:04d}.pt")
    logit_data = torch.load(logit_path, map_location="cpu", weights_only=True)
    margins = (logit_data["top5_values"][:, 0] - logit_data["top5_values"][:, 1]).numpy()
    fragility = float((margins < 1.0).mean())

    entry = {
        "idx": idx,
        "set": meta["set"],
        "category": meta["category"],
        "fragility": fragility,
        "input_len": input_len,
        "l23": l23,
        "l22": l22,
        "residual": residual,
        "gen_start": gen_start,
    }

    if meta["set"] == "experimental":
        exp_data.append(entry)
    else:
        ctrl_data.append(entry)

print(f"Experimental: {len(exp_data)}, Control: {len(ctrl_data)}", flush=True)

# ============================================================
# PART 1: Per-dimension residual norms (averaged over gen positions)
# ============================================================
print("\n" + "=" * 70)
print("PART 1: Per-dimension analysis of layer-23 residual")
print("=" * 70)

# For each prompt, compute per-dimension mean absolute residual over generated positions
def get_gen_dim_profile(entry):
    """Mean absolute value of each dimension over generated positions."""
    res = entry["residual"][entry["gen_start"]:, :]  # (gen_len, 1024)
    return res.abs().mean(dim=0).numpy()  # (1024,)

exp_profiles = np.stack([get_gen_dim_profile(e) for e in exp_data])  # (25, 1024)
ctrl_profiles = np.stack([get_gen_dim_profile(e) for e in ctrl_data])  # (25, 1024)

exp_mean = exp_profiles.mean(axis=0)  # (1024,)
ctrl_mean = ctrl_profiles.mean(axis=0)  # (1024,)
diff = exp_mean - ctrl_mean  # positive = larger in experimental

# Find dimensions with largest absolute difference
top_dims = np.argsort(-np.abs(diff))

print(f"\nTop 20 dimensions by |exp_mean - ctrl_mean| of layer-23 residual:")
print(f"{'Dim':>5s} {'ExpMean':>9s} {'CtrlMean':>9s} {'Diff':>9s} {'Ratio':>7s}")
print("-" * 45)
for rank, dim in enumerate(top_dims[:20]):
    e, c = exp_mean[dim], ctrl_mean[dim]
    ratio = e / c if c > 0.001 else float("inf")
    print(f"{dim:5d} {e:9.4f} {c:9.4f} {diff[dim]:+9.4f} {ratio:7.2f}x")

# ============================================================
# PART 2: Is there a dominant dimension (like dim 487 for BOS)?
# ============================================================
print("\n" + "=" * 70)
print("PART 2: Concentration of exp/ctrl difference")
print("=" * 70)

abs_diff = np.abs(diff)
total_diff = abs_diff.sum()
sorted_diffs = np.sort(abs_diff)[::-1]
cumulative = np.cumsum(sorted_diffs) / total_diff

for k in [1, 3, 5, 10, 20, 50, 100]:
    print(f"  Top {k:3d} dims explain {cumulative[k-1]*100:.1f}% of total |diff|")

# ============================================================
# PART 3: Per-dimension correlation with fragility
# ============================================================
print("\n" + "=" * 70)
print("PART 3: Dimension-fragility correlations")
print("=" * 70)

all_profiles = np.vstack([exp_profiles, ctrl_profiles])  # (50, 1024)
all_frags = np.array([e["fragility"] for e in exp_data + ctrl_data])

dim_corrs = np.array([np.corrcoef(all_profiles[:, d], all_frags)[0, 1]
                       for d in range(1024)])

top_pos_corr = np.argsort(-dim_corrs)
top_neg_corr = np.argsort(dim_corrs)

print(f"\nTop 10 positively correlated dimensions (more active when fragile):")
for rank, dim in enumerate(top_pos_corr[:10]):
    print(f"  dim {dim:4d}: r = {dim_corrs[dim]:+.4f}  exp_mean={exp_mean[dim]:.4f}  ctrl_mean={ctrl_mean[dim]:.4f}")

print(f"\nTop 10 negatively correlated dimensions (more active when confident):")
for rank, dim in enumerate(top_neg_corr[:10]):
    print(f"  dim {dim:4d}: r = {dim_corrs[dim]:+.4f}  exp_mean={exp_mean[dim]:.4f}  ctrl_mean={ctrl_mean[dim]:.4f}")

# How many dims have |r| > 0.5?
strong_corr = (np.abs(dim_corrs) > 0.5).sum()
very_strong = (np.abs(dim_corrs) > 0.7).sum()
extreme = (np.abs(dim_corrs) > 0.9).sum()
print(f"\nDimensions with |r| > 0.5: {strong_corr}/1024 ({strong_corr/1024*100:.1f}%)")
print(f"Dimensions with |r| > 0.7: {very_strong}/1024 ({very_strong/1024*100:.1f}%)")
print(f"Dimensions with |r| > 0.9: {extreme}/1024 ({extreme/1024*100:.1f}%)")

# ============================================================
# PART 4: Can a small set of dims predict exp vs ctrl?
# ============================================================
print("\n" + "=" * 70)
print("PART 4: Classification accuracy using top-k dimensions")
print("=" * 70)

labels = np.array([1] * len(exp_data) + [0] * len(ctrl_data))

for k in [1, 3, 5, 10, 20, 50]:
    # Use top-k dimensions by abs correlation
    top_k_dims = top_pos_corr[:k]
    features = all_profiles[:, top_k_dims].mean(axis=1)  # mean of top-k dims
    threshold = np.median(features)
    predictions = (features > threshold).astype(int)
    accuracy = (predictions == labels).mean()
    print(f"  Top {k:3d} dims (mean, median-split): accuracy = {accuracy:.1%}")

# Also try: sum of signed contributions
print(f"\n  Using signed correlation weighting:")
for k in [5, 10, 20, 50]:
    top_k = np.argsort(-np.abs(dim_corrs))[:k]
    weights = dim_corrs[top_k]
    score = (all_profiles[:, top_k] * weights[None, :]).sum(axis=1)
    threshold = np.median(score)
    predictions = (score > threshold).astype(int)
    accuracy = (predictions == labels).mean()
    print(f"  Top {k:3d} dims (correlation-weighted): accuracy = {accuracy:.1%}")

# ============================================================
# PART 5: Position-dependent residual analysis
# ============================================================
print("\n" + "=" * 70)
print("PART 5: Layer-23 residual magnitude over generation positions")
print("=" * 70)

# For a few representative prompts, show how the residual evolves
# across generated positions

# Compare: one very fragile, one very confident
exp_sorted = sorted(exp_data, key=lambda e: -e["fragility"])
ctrl_sorted = sorted(ctrl_data, key=lambda e: e["fragility"])

for label, entries in [("Most fragile exp", exp_sorted[:3]),
                        ("Most confident ctrl", ctrl_sorted[:3])]:
    print(f"\n{label}:")
    for entry in entries:
        res = entry["residual"]
        gs = entry["gen_start"]
        gen_res = res[gs:, :]  # generated portion
        n_gen = gen_res.shape[0]

        # Residual L2 norm at each generated position
        pos_norms = gen_res.norm(dim=1).numpy()

        # Show in windows
        windows = [(0, min(8, n_gen)),
                    (max(0, n_gen//4), min(n_gen//4+8, n_gen)),
                    (max(0, n_gen//2), min(n_gen//2+8, n_gen)),
                    (max(0, n_gen-8), n_gen)]
        win_means = []
        for lo, hi in windows:
            if hi > lo:
                win_means.append(f"{pos_norms[lo:hi].mean():.2f}")
            else:
                win_means.append("N/A")
        print(f"  idx={entry['idx']:4d} frag={entry['fragility']:.3f} "
              f"res_norms=[start:{win_means[0]}, Q1:{win_means[1]}, "
              f"mid:{win_means[2]}, end:{win_means[3]}]  "
              f"total_mean={pos_norms.mean():.2f}")

# ============================================================
# PART 6: The "uncertainty subspace" — PCA of the difference
# ============================================================
print("\n" + "=" * 70)
print("PART 6: PCA of layer-23 residual difference")
print("=" * 70)

# Center each prompt's profile, then do PCA
all_centered = all_profiles - all_profiles.mean(axis=0, keepdims=True)

# SVD
U, S, Vt = np.linalg.svd(all_centered, full_matrices=False)

# How much variance do top components explain?
total_var = (S ** 2).sum()
for k in [1, 2, 3, 5, 10]:
    explained = (S[:k] ** 2).sum() / total_var
    print(f"  Top {k:2d} PCs explain {explained*100:.1f}% of variance")

# Does PC1 separate exp/ctrl?
pc1_scores = U[:, 0] * S[0]
exp_pc1 = pc1_scores[:len(exp_data)]
ctrl_pc1 = pc1_scores[len(exp_data):]
print(f"\n  PC1 scores: exp mean={exp_pc1.mean():.4f} ctrl mean={ctrl_pc1.mean():.4f}")
d = (exp_pc1.mean() - ctrl_pc1.mean()) / np.sqrt((exp_pc1.var() + ctrl_pc1.var()) / 2)
print(f"  PC1 Cohen's d = {d:.2f}")

# Which dimensions load most on PC1?
pc1_loadings = Vt[0, :]  # (1024,)
top_pc1_dims = np.argsort(-np.abs(pc1_loadings))
print(f"\n  Top 10 dimensions loading on PC1:")
for dim in top_pc1_dims[:10]:
    print(f"    dim {dim:4d}: loading = {pc1_loadings[dim]:+.4f}")

# ============================================================
# PART 7: Overlap between dim-487 (BOS sink) and uncertainty dims
# ============================================================
print("\n" + "=" * 70)
print("PART 7: Is dim 487 involved in the uncertainty signal?")
print("=" * 70)

print(f"  Dim 487 correlation with fragility: r = {dim_corrs[487]:+.4f}")
print(f"  Dim 487 exp mean: {exp_mean[487]:.4f}, ctrl mean: {ctrl_mean[487]:.4f}")
print(f"  Dim 487 rank in |correlation|: {np.where(np.argsort(-np.abs(dim_corrs)) == 487)[0][0] + 1}")
print(f"  Dim 487 PC1 loading: {pc1_loadings[487]:+.6f}")

# Check if the top uncertainty dims overlap with BOS-sink dims
# (BOS sink was dominated by dim 487 at layer 5; check if same dims appear at layer 23)
print(f"\n  Do uncertainty dims overlap with BOS-active dims?")
# Load one activation to check BOS pattern at layer 5
sample = torch.load(os.path.join(ACT_DIR, f"prompt_{subset[0]['prompt_idx']:04d}.pt"),
                     map_location="cpu", weights_only=False)
l5_bos = sample["layer_5"][0, :].numpy()  # BOS token at layer 5
l5_top = np.argsort(-np.abs(l5_bos))[:10]
l23_top_corr = top_pos_corr[:10]
overlap = set(l5_top) & set(l23_top_corr)
print(f"  Layer-5 BOS top-10 dims: {list(l5_top)}")
print(f"  Layer-23 uncertainty top-10 dims: {list(l23_top_corr)}")
print(f"  Overlap: {overlap if overlap else 'NONE'}")
