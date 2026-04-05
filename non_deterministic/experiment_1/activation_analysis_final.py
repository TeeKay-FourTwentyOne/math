"""
Final targeted analyses: input boundary effects, fragility correlations, effect sizes.
"""

import torch
import numpy as np
import json
from pathlib import Path
from scipy import stats

DATA_DIR = Path("C:/Projects/math/non_deterministic/experiment_1/activations/h_nvidia")
META_FILE = Path("C:/Projects/math/non_deterministic/experiment_1/activations/activation_subset.jsonl")

LAYER_NAMES = ["embed"] + [f"layer_{i}" for i in range(24)] + ["final_norm"]

meta = {}
with open(META_FILE) as f:
    for line in f:
        if line.strip():
            d = json.loads(line)
            meta[d["prompt_idx"]] = d

data = {"experimental": [], "control": []}
for pt_file in sorted(DATA_DIR.glob("prompt_*.pt")):
    d = torch.load(pt_file, map_location="cpu", weights_only=False)
    entry = {
        "prompt_idx": d["prompt_idx"],
        "input_length": d["input_length"],
        "total_length": d["total_length"],
        "set": d["set"],
        "category": d["category"],
        "fragility": meta[d["prompt_idx"]]["fragility_score"],
    }
    for ln in LAYER_NAMES:
        entry[ln] = d[ln].numpy()
    data[entry["set"]].append(entry)


print("="*80)
print("FINAL ANALYSIS 1: INPUT BOUNDARY TRANSITION")
print("="*80)
print("What happens at the exact boundary between input and generated tokens?")

# For each prompt, look at positions around input_length
for ln in ["layer_6", "layer_12", "layer_23"]:
    print(f"\n  {ln}:")
    for s in ["experimental", "control"]:
        norms_before = []  # last 3 input tokens
        norms_after = []   # first 3 generated tokens
        for e in data[s]:
            il = e["input_length"]
            act = e[ln]
            before = np.linalg.norm(act[max(0,il-3):il], axis=1)
            after = np.linalg.norm(act[il:il+3], axis=1)
            norms_before.append(np.mean(before))
            norms_after.append(np.mean(after))
        ratio = np.mean(norms_after) / np.mean(norms_before)
        print(f"    {s:<15} last 3 input: {np.mean(norms_before):.2f}  first 3 gen: {np.mean(norms_after):.2f}  ratio: {ratio:.4f}")


print("\n" + "="*80)
print("FINAL ANALYSIS 2: FRAGILITY SUBGROUPS WITHIN EXPERIMENTAL")
print("="*80)

# Split experimental into higher vs lower fragility
exp_sorted = sorted(data["experimental"], key=lambda x: x["fragility"])
low_frag = exp_sorted[:12]  # lower half of experimental
high_frag = exp_sorted[13:]  # upper half of experimental

print(f"Low-frag experimental: n={len(low_frag)}, mean frag={np.mean([e['fragility'] for e in low_frag]):.4f}")
print(f"High-frag experimental: n={len(high_frag)}, mean frag={np.mean([e['fragility'] for e in high_frag]):.4f}")

print(f"\n{'Layer':<15} {'Low-frag norm':>15} {'High-frag norm':>15} {'Ratio (H/L)':>12}")
for ln in LAYER_NAMES:
    lf_norms = [np.mean(np.linalg.norm(e[ln], axis=1)) for e in low_frag]
    hf_norms = [np.mean(np.linalg.norm(e[ln], axis=1)) for e in high_frag]
    ratio = np.mean(hf_norms) / np.mean(lf_norms) if np.mean(lf_norms) > 0 else 0
    print(f"  {ln:<13} {np.mean(lf_norms):>13.2f} {np.mean(hf_norms):>13.2f} {ratio:>10.4f}")


print("\n" + "="*80)
print("FINAL ANALYSIS 3: COHEN'S d EFFECT SIZES")
print("="*80)

def cohens_d(group1, group2):
    n1, n2 = len(group1), len(group2)
    m1, m2 = np.mean(group1), np.mean(group2)
    s1, s2 = np.std(group1, ddof=1), np.std(group2, ddof=1)
    sp = np.sqrt(((n1-1)*s1**2 + (n2-1)*s2**2) / (n1+n2-2))
    return (m1 - m2) / (sp + 1e-10)

print(f"\n{'Layer':<15} {'Cohen d (norm)':>15} {'Effect size':>15}")
for ln in LAYER_NAMES:
    exp_vals = [np.mean(np.linalg.norm(e[ln], axis=1)) for e in data["experimental"]]
    ctrl_vals = [np.mean(np.linalg.norm(e[ln], axis=1)) for e in data["control"]]
    d = cohens_d(exp_vals, ctrl_vals)
    if abs(d) < 0.2:
        mag = "negligible"
    elif abs(d) < 0.5:
        mag = "small"
    elif abs(d) < 0.8:
        mag = "medium"
    else:
        mag = "LARGE"
    print(f"  {ln:<13} {d:>+13.4f} {mag:>15}")


print("\n" + "="*80)
print("FINAL ANALYSIS 4: NORM PROFILE OF GENERATED TOKENS OVER TIME")
print("="*80)
print("How do norms evolve across generated positions (first 50 gen tokens)?")

for ln in ["layer_12", "layer_23", "final_norm"]:
    print(f"\n  {ln}:")
    for s in ["experimental", "control"]:
        bucket_norms = {i: [] for i in range(10)}  # 10 buckets of 5 tokens each
        for e in data[s]:
            il = e["input_length"]
            gen_act = e[ln][il:]  # generated portion
            for bucket in range(10):
                start = bucket * 5
                end = start + 5
                if start < gen_act.shape[0]:
                    chunk = gen_act[start:min(end, gen_act.shape[0])]
                    bucket_norms[bucket].append(np.mean(np.linalg.norm(chunk, axis=1)))
        print(f"    {s:<13}", end="")
        for bucket in range(10):
            if bucket_norms[bucket]:
                print(f" [{bucket*5}-{bucket*5+4}]:{np.mean(bucket_norms[bucket]):.1f}", end="")
        print()


print("\n" + "="*80)
print("FINAL ANALYSIS 5: DIMENSION 487 — THE DOMINANT FEATURE")
print("="*80)
print("Dim 487 carries ~85% of layer_5 variance. What is its activation profile?")

for ln in ["layer_4", "layer_5", "layer_6", "layer_12", "layer_23"]:
    exp_d487 = [np.mean(e[ln][:, 487]) for e in data["experimental"]]
    ctrl_d487 = [np.mean(e[ln][:, 487]) for e in data["control"]]
    exp_d487_pos0 = [e[ln][0, 487] for e in data["experimental"]]
    ctrl_d487_pos0 = [e[ln][0, 487] for e in data["control"]]
    print(f"  {ln:<15} dim487 mean: exp={np.mean(exp_d487):.2f} ctrl={np.mean(ctrl_d487):.2f}"
          f"  |  pos0: exp={np.mean(exp_d487_pos0):.2f} ctrl={np.mean(ctrl_d487_pos0):.2f}")

# Check if dim 487 is responsible for the BOS anomaly
print("\nBOS token norm with vs without dim 487:")
for ln in ["layer_5", "layer_12"]:
    for s in ["experimental", "control"]:
        full_norms = []
        no487_norms = []
        for e in data[s]:
            bos = e[ln][0]  # (1024,)
            full_norms.append(np.linalg.norm(bos))
            bos_no487 = np.delete(bos, 487)
            no487_norms.append(np.linalg.norm(bos_no487))
        print(f"  {ln} {s:<13}: with_487={np.mean(full_norms):.2f}  without_487={np.mean(no487_norms):.2f}  "
              f"dim487_contribution={np.mean(full_norms)**2 - np.mean(no487_norms)**2:.2f}")


print("\n" + "="*80)
print("FINAL ANALYSIS 6: LAYER 23 LAST-TOKEN vs FRAGILITY (direct regression)")
print("="*80)

all_entries = data["experimental"] + data["control"]
frags = np.array([e["fragility"] for e in all_entries])
l23_last = np.array([np.linalg.norm(e["layer_23"][-1]) for e in all_entries])
l23_mean = np.array([np.mean(np.linalg.norm(e["layer_23"], axis=1)) for e in all_entries])

slope_last, intercept_last, r_last, p_last, se_last = stats.linregress(frags, l23_last)
slope_mean, intercept_mean, r_mean, p_mean, se_mean = stats.linregress(frags, l23_mean)

print(f"  Last-token norm vs fragility: r={r_last:.4f}, p={p_last:.2e}, slope={slope_last:.2f}")
print(f"  Mean norm vs fragility:       r={r_mean:.4f}, p={p_mean:.2e}, slope={slope_mean:.2f}")

# Also check final_norm
fn_last = np.array([np.linalg.norm(e["final_norm"][-1]) for e in all_entries])
fn_mean = np.array([np.mean(np.linalg.norm(e["final_norm"], axis=1)) for e in all_entries])
slope_fn, intercept_fn, r_fn, p_fn, se_fn = stats.linregress(frags, fn_last)
print(f"  Final_norm last-token:        r={r_fn:.4f}, p={p_fn:.2e}")
print(f"  Final_norm mean:              r={np.corrcoef(frags, fn_mean)[0,1]:.4f}")


print("\n" + "="*80)
print("SUMMARY TABLE: KEY FINDINGS")
print("="*80)

print("""
+----------------------------------+--------------------------------------------+
| Finding                          | Detail                                     |
+----------------------------------+--------------------------------------------+
| BOS token anomaly                | Pos-0 norm 40-50x larger at layers 5-21    |
|                                  | (single dim 487 accounts for ~85%)         |
+----------------------------------+--------------------------------------------+
| Layer 5 phase transition         | Variance jumps 30x from layer 4 to 5       |
|                                  | Dim 487 activates explosively at BOS        |
+----------------------------------+--------------------------------------------+
| Activation growth                | ~144-154x from embed to final_norm         |
|                                  | Neither linear (R2=.63) nor exponential    |
|                                  | (.52); multi-phase growth                  |
+----------------------------------+--------------------------------------------+
| Control > Experimental (L10-21)  | Control has LARGER norms in middle/late     |
|                                  | layers; Cohen's d up to -8.6 at layer_11   |
+----------------------------------+--------------------------------------------+
| Layer 23 reversal                | Experimental norms 60% LARGER than control |
|                                  | (d=+3.35); only deep layer where exp>ctrl  |
+----------------------------------+--------------------------------------------+
| Fragility-norm correlation       | r=+0.97 at layer_23, r=-0.97 at layer_17   |
|                                  | Perfect sign flip between L17 and L23      |
+----------------------------------+--------------------------------------------+
| Input vs generated tokens        | Generated/input ratio = 0.12-0.16 (exp)    |
|                                  | at layers 5-16, vs 0.25-0.34 (ctrl)        |
|                                  | Experimental input tokens dominate more     |
+----------------------------------+--------------------------------------------+
| Effective dimensionality         | Control uses 10-30x MORE dimensions in     |
|                                  | layers 13-21 (PR: 114 vs 22 at L20)       |
|                                  | Experimental is more concentrated           |
+----------------------------------+--------------------------------------------+
| Intra-group coherence            | Experimental more similar to each other     |
|                                  | at embed/L0/L23; control more similar at   |
|                                  | middle layers                              |
+----------------------------------+--------------------------------------------+
| Last-token layer_23 norm         | exp=53.6 vs ctrl=26.3 (2.04x ratio)       |
|                                  | Last token cos(t,t-1): 0.61 vs 0.35       |
|                                  | Uncertain model repeats more at end        |
+----------------------------------+--------------------------------------------+
| Layer 23 residual > output       | Residual/output ratio = 1.08; layer 23     |
|                                  | partially CANCELS the residual stream      |
+----------------------------------+--------------------------------------------+
""")

print("All analyses complete.")
