"""
Iteration 7: Tracing the Uncertainty Axis Through All Layers

Iter 6 found a single PC explaining 83.5% of variance in layer-23 residuals,
with perfect fragility classification. Layer 17 has r=-0.97 with fragility
(inverted sign!). Key questions:

1. At each layer, how much variance does PC1 explain?
2. How well does each layer's PC1 correlate with fragility?
3. Are the PC1 DIRECTIONS consistent across layers (cosine similarity)?
4. When does the axis emerge, and when does it flip?
5. Is layer 23's axis just the negation of an earlier layer's axis?
"""

import torch
import numpy as np
import json
import os
from pathlib import Path

ACT_DIR = "C:/Projects/math/non_deterministic/experiment_1/activations/h_nvidia"
SUBSET_PATH = "C:/Projects/math/non_deterministic/experiment_1/activations/activation_subset.jsonl"
LOGIT_DIR = "C:/Projects/math/non_deterministic/experiment_0/logit_traces"

# Load subset metadata + fragility
subset = []
with open(SUBSET_PATH, encoding="utf-8") as f:
    for line in f:
        subset.append(json.loads(line))

frags = []
for meta in subset:
    idx = meta["prompt_idx"]
    logit_path = os.path.join(LOGIT_DIR, f"prompt_{idx:04d}.pt")
    data = torch.load(logit_path, map_location="cpu", weights_only=True)
    margins = (data["top5_values"][:, 0] - data["top5_values"][:, 1]).numpy()
    frags.append(float((margins < 1.0).mean()))
frags = np.array(frags)

# Layer names in order
layer_names = ["embed"] + [f"layer_{i}" for i in range(24)] + ["final_norm"]

# ============================================================
# Load activations and compute per-layer profiles
# ============================================================
print("Loading activations...", flush=True)

all_acts = []  # list of dicts, one per prompt
for meta in subset:
    idx = meta["prompt_idx"]
    path = os.path.join(ACT_DIR, f"prompt_{idx:04d}.pt")
    data = torch.load(path, map_location="cpu", weights_only=False)
    all_acts.append(data)

N = len(all_acts)
print(f"Loaded {N} prompts", flush=True)

def get_gen_mean_profile(acts_list, layer_name):
    """For each prompt, mean activation over generated positions at given layer."""
    profiles = []
    for data in acts_list:
        act = data[layer_name]  # (seq_len, 1024)
        gs = data["input_length"]
        gen = act[gs:, :]  # generated portion
        profiles.append(gen.mean(dim=0).numpy())  # (1024,)
    return np.stack(profiles)  # (N, 1024)

def get_residual_profile(acts_list, layer_name, prev_layer_name):
    """Layer residual = layer_k - layer_{k-1}, mean over generated positions."""
    profiles = []
    for data in acts_list:
        gs = data["input_length"]
        cur = data[layer_name][gs:, :]
        prev = data[prev_layer_name][gs:, :]
        res = cur - prev
        profiles.append(res.mean(dim=0).numpy())
    return np.stack(profiles)

# ============================================================
# PART 1: PCA at each layer (using layer OUTPUTS, not residuals)
# ============================================================
print("\n" + "=" * 70)
print("PART 1: PCA of layer outputs across all layers")
print("=" * 70)

layer_pca_results = {}

for lname in layer_names:
    profiles = get_gen_mean_profile(all_acts, lname)  # (50, 1024)

    # Center
    centered = profiles - profiles.mean(axis=0, keepdims=True)

    # SVD
    U, S, Vt = np.linalg.svd(centered, full_matrices=False)

    # Variance explained by PC1
    total_var = (S ** 2).sum()
    pc1_var = (S[0] ** 2) / total_var if total_var > 0 else 0

    # PC1 score correlation with fragility
    pc1_scores = U[:, 0] * S[0]
    r_frag = np.corrcoef(pc1_scores, frags)[0, 1]

    # Store
    layer_pca_results[lname] = {
        "pc1_var": pc1_var,
        "r_frag": r_frag,
        "pc1_direction": Vt[0, :],  # (1024,) — the direction
        "pc1_scores": pc1_scores,
        "S": S,
    }

print(f"\n{'Layer':12s} {'PC1 VarExp':>10s} {'r(PC1,frag)':>12s} {'|r|':>6s}")
print("-" * 44)
for lname in layer_names:
    res = layer_pca_results[lname]
    print(f"{lname:12s} {res['pc1_var']:10.3f} {res['r_frag']:+12.4f} {abs(res['r_frag']):6.4f}")

# ============================================================
# PART 2: PCA of RESIDUALS (layer_k - layer_{k-1})
# ============================================================
print("\n" + "=" * 70)
print("PART 2: PCA of layer residuals")
print("=" * 70)

residual_pca = {}
prev_name = "embed"
for i in range(24):
    cur_name = f"layer_{i}"
    profiles = get_residual_profile(all_acts, cur_name, prev_name)
    centered = profiles - profiles.mean(axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(centered, full_matrices=False)
    total_var = (S ** 2).sum()
    pc1_var = (S[0] ** 2) / total_var if total_var > 0 else 0
    pc1_scores = U[:, 0] * S[0]
    r_frag = np.corrcoef(pc1_scores, frags)[0, 1]
    residual_pca[cur_name] = {
        "pc1_var": pc1_var,
        "r_frag": r_frag,
        "pc1_direction": Vt[0, :],
    }
    prev_name = cur_name

# Also final_norm residual
profiles = get_residual_profile(all_acts, "final_norm", "layer_23")
centered = profiles - profiles.mean(axis=0, keepdims=True)
U, S, Vt = np.linalg.svd(centered, full_matrices=False)
total_var = (S ** 2).sum()
pc1_var = (S[0] ** 2) / total_var if total_var > 0 else 0
pc1_scores = U[:, 0] * S[0]
r_frag = np.corrcoef(pc1_scores, frags)[0, 1]
residual_pca["final_norm"] = {"pc1_var": pc1_var, "r_frag": r_frag, "pc1_direction": Vt[0, :]}

print(f"\n{'Layer residual':15s} {'PC1 VarExp':>10s} {'r(PC1,frag)':>12s}")
print("-" * 40)
for lname in [f"layer_{i}" for i in range(24)] + ["final_norm"]:
    res = residual_pca[lname]
    print(f"{lname:15s} {res['pc1_var']:10.3f} {res['r_frag']:+12.4f}")

# ============================================================
# PART 3: Cosine similarity of PC1 directions across layers
# ============================================================
print("\n" + "=" * 70)
print("PART 3: Cosine similarity of PC1 directions between layers")
print("=" * 70)

# Use the OUTPUT-based PCA directions
def cosine_sim(a, b):
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-10))

# Compare each layer's PC1 to layer 23's PC1
ref_dir = layer_pca_results["layer_23"]["pc1_direction"]
print(f"\nCosine similarity of each layer's PC1 with LAYER 23's PC1:")
print(f"{'Layer':12s} {'cos(PC1, L23)':>14s}")
print("-" * 28)
for lname in layer_names:
    d = layer_pca_results[lname]["pc1_direction"]
    cs = cosine_sim(d, ref_dir)
    bar = "+" * int(abs(cs) * 30) if cs > 0 else "-" * int(abs(cs) * 30)
    print(f"{lname:12s} {cs:+14.4f}  {bar}")

# Pairwise between a few key layers
print(f"\nPairwise cosine similarities between key layers:")
key_layers = ["layer_5", "layer_10", "layer_15", "layer_17", "layer_20", "layer_22", "layer_23"]
print(f"{'':12s}", end="")
for l2 in key_layers:
    print(f" {l2[:7]:>7s}", end="")
print()
for l1 in key_layers:
    print(f"{l1:12s}", end="")
    for l2 in key_layers:
        d1 = layer_pca_results[l1]["pc1_direction"]
        d2 = layer_pca_results[l2]["pc1_direction"]
        cs = cosine_sim(d1, d2)
        print(f" {cs:+7.3f}", end="")
    print()

# ============================================================
# PART 4: Where does the sign flip happen?
# ============================================================
print("\n" + "=" * 70)
print("PART 4: The sign flip — tracking r(PC1, fragility) across layers")
print("=" * 70)

print(f"\nLayer-by-layer correlation of PC1 with fragility:")
for lname in layer_names:
    r = layer_pca_results[lname]["r_frag"]
    bar_len = int(abs(r) * 40)
    if r > 0:
        bar = " " * 40 + "|" + "+" * bar_len
    else:
        bar = " " * (40 - bar_len) + "-" * bar_len + "|"
    print(f"{lname:12s} {r:+.4f} {bar}")

# ============================================================
# PART 5: When does the uncertainty axis EMERGE?
# ============================================================
print("\n" + "=" * 70)
print("PART 5: Emergence of the uncertainty axis")
print("=" * 70)

# Classification accuracy at each layer (using PC1 median-split)
print(f"\n{'Layer':12s} {'|r|':>6s} {'Accuracy':>9s} {'VarExp':>8s}")
print("-" * 38)
labels = np.array([1 if m["set"] == "experimental" else 0 for m in subset])
for lname in layer_names:
    res = layer_pca_results[lname]
    scores = res["pc1_scores"]
    threshold = np.median(scores)
    preds = (scores > threshold).astype(int)
    # If r is negative, flip predictions
    if res["r_frag"] < 0:
        preds = 1 - preds
    acc = (preds == labels).mean()
    print(f"{lname:12s} {abs(res['r_frag']):6.4f} {acc:9.1%} {res['pc1_var']:8.3f}")
