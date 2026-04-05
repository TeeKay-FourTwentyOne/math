"""
Deep-dive analysis on surprising findings from initial pass.
"""

import torch
import numpy as np
import json
from pathlib import Path
from collections import defaultdict

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
        "text": meta[d["prompt_idx"]]["text"],
    }
    for ln in LAYER_NAMES:
        entry[ln] = d[ln].numpy()
    data[entry["set"]].append(entry)

print("="*80)
print("DEEP DIVE 1: POSITION-0 ANOMALY (BOS TOKEN)")
print("="*80)

# Position 0 had massive norms in middle layers. Investigate.
print("\nL2 norm at position 0 vs mean of positions 1+ at every layer:")
print(f"{'Layer':<15} {'Pos0 (exp)':>12} {'Rest (exp)':>12} {'Ratio':>8} {'Pos0 (ctrl)':>12} {'Rest (ctrl)':>12} {'Ratio':>8}")
for ln in LAYER_NAMES:
    exp_pos0 = np.mean([np.linalg.norm(e[ln][0]) for e in data["experimental"]])
    exp_rest = np.mean([np.mean(np.linalg.norm(e[ln][1:], axis=1)) for e in data["experimental"]])
    ctrl_pos0 = np.mean([np.linalg.norm(e[ln][0]) for e in data["control"]])
    ctrl_rest = np.mean([np.mean(np.linalg.norm(e[ln][1:], axis=1)) for e in data["control"]])
    er = exp_pos0 / exp_rest if exp_rest > 0 else 0
    cr = ctrl_pos0 / ctrl_rest if ctrl_rest > 0 else 0
    print(f"  {ln:<13} {exp_pos0:>12.2f} {exp_rest:>12.2f} {er:>7.1f}x {ctrl_pos0:>12.2f} {ctrl_rest:>12.2f} {cr:>7.1f}x")

print("\n" + "="*80)
print("DEEP DIVE 2: THE LAYER_23 REVERSAL")
print("="*80)
print("\nLayer 23 is the only deep layer where experimental > control.")
print("Individual prompt norms at layer_23:")

exp_l23_norms = []
ctrl_l23_norms = []
for e in data["experimental"]:
    n = np.mean(np.linalg.norm(e["layer_23"], axis=1))
    exp_l23_norms.append((n, e["fragility"], e["text"][:60]))
for e in data["control"]:
    n = np.mean(np.linalg.norm(e["layer_23"], axis=1))
    ctrl_l23_norms.append((n, e["fragility"], e["text"][:60]))

print("\nExperimental (sorted by norm):")
for n, f, t in sorted(exp_l23_norms, reverse=True)[:10]:
    print(f"  norm={n:.2f}  frag={f:.4f}  {t}")
print("\nControl (sorted by norm):")
for n, f, t in sorted(ctrl_l23_norms, reverse=True)[:10]:
    print(f"  norm={n:.2f}  frag={f:.4f}  {t}")

# Correlation between fragility and layer_23 norm
all_frags = [e["fragility"] for e in data["experimental"] + data["control"]]
all_l23 = [np.mean(np.linalg.norm(e["layer_23"], axis=1)) for e in data["experimental"] + data["control"]]
corr = np.corrcoef(all_frags, all_l23)[0, 1]
print(f"\nCorrelation(fragility, layer_23 mean norm) = {corr:.4f}")

# Check per-layer correlation with fragility
print("\nCorrelation of fragility with mean L2 norm at each layer:")
for ln in LAYER_NAMES:
    all_norms = [np.mean(np.linalg.norm(e[ln], axis=1)) for e in data["experimental"] + data["control"]]
    c = np.corrcoef(all_frags, all_norms)[0, 1]
    print(f"  {ln:<15} r = {c:+.4f}")


print("\n" + "="*80)
print("DEEP DIVE 3: LAYER 5 PHASE TRANSITION")
print("="*80)
print("\nLayers 4->5 show a massive jump in variance (0.26 -> 6.07 for experimental).")
print("Examining what changes:")

for s_name, entries in [("Experimental", data["experimental"][:3]), ("Control", data["control"][:3])]:
    print(f"\n{s_name} (first 3 prompts):")
    for e in entries:
        l4 = e["layer_4"]
        l5 = e["layer_5"]
        diff = l5 - l4
        # Check if specific dimensions light up
        diff_mean = np.mean(np.abs(diff), axis=0)  # mean across positions per dim
        top_dims = np.argsort(diff_mean)[::-1][:10]
        print(f"  prompt {e['prompt_idx']}: top changed dims = {top_dims.tolist()}")
        print(f"    layer_4 norm: mean={np.mean(np.linalg.norm(l4, axis=1)):.2f}, max_dim_var={np.max(np.var(l4, axis=0)):.4f}")
        print(f"    layer_5 norm: mean={np.mean(np.linalg.norm(l5, axis=1)):.2f}, max_dim_var={np.max(np.var(l5, axis=0)):.4f}")

# Check which dimensions have the highest variance in layer 5
print("\nDimension analysis at layer_5 (averaged across all prompts):")
all_vars_l5 = []
for e in data["experimental"] + data["control"]:
    v = np.var(e["layer_5"], axis=0)  # (1024,)
    all_vars_l5.append(v)
avg_var_l5 = np.mean(all_vars_l5, axis=0)
top20 = np.argsort(avg_var_l5)[::-1][:20]
print(f"  Top 20 high-variance dimensions: {top20.tolist()}")
print(f"  Their variances: {avg_var_l5[top20].tolist()[:10]}")
print(f"  Total variance sum: {np.sum(avg_var_l5):.2f}, top20 share: {np.sum(avg_var_l5[top20])/np.sum(avg_var_l5)*100:.1f}%")

# Same for layer 4
all_vars_l4 = []
for e in data["experimental"] + data["control"]:
    v = np.var(e["layer_4"], axis=0)
    all_vars_l4.append(v)
avg_var_l4 = np.mean(all_vars_l4, axis=0)
top20_l4 = np.argsort(avg_var_l4)[::-1][:20]
print(f"\n  Layer 4 top 20 dims: {top20_l4.tolist()}")
print(f"  Their variances: {avg_var_l4[top20_l4].tolist()[:10]}")
print(f"  Total variance sum: {np.sum(avg_var_l4):.2f}")


print("\n" + "="*80)
print("DEEP DIVE 4: COSINE SIMILARITY WITHIN GROUPS")
print("="*80)
print("How similar are prompts within the same group at each layer?")

def mean_pairwise_cosine(entries, layer_name, n_samples=None):
    """Mean pairwise cosine similarity of position-averaged activations."""
    vecs = []
    for e in entries:
        v = np.mean(e[layer_name], axis=0)
        vecs.append(v / (np.linalg.norm(v) + 1e-10))
    sims = []
    for i in range(len(vecs)):
        for j in range(i+1, len(vecs)):
            sims.append(np.dot(vecs[i], vecs[j]))
    return np.mean(sims), np.std(sims)

print(f"\n{'Layer':<15} {'Exp intra-cos':>15} {'Ctrl intra-cos':>15} {'Diff':>10}")
for ln in LAYER_NAMES:
    em, es = mean_pairwise_cosine(data["experimental"], ln)
    cm, cs = mean_pairwise_cosine(data["control"], ln)
    print(f"  {ln:<13} {em:>13.4f}±{es:.4f} {cm:>13.4f}±{cs:.4f} {em-cm:>+10.4f}")


print("\n" + "="*80)
print("DEEP DIVE 5: LAST-TOKEN ANALYSIS (CRITICAL FOR NEXT-TOKEN PREDICTION)")
print("="*80)
print("The last token's activation drives the next-token prediction.")

for ln in ["layer_20", "layer_21", "layer_22", "layer_23", "final_norm"]:
    exp_last = [np.linalg.norm(e[ln][-1]) for e in data["experimental"]]
    ctrl_last = [np.linalg.norm(e[ln][-1]) for e in data["control"]]
    print(f"  {ln:<15} last-token norm: exp={np.mean(exp_last):.2f}±{np.std(exp_last):.2f}  ctrl={np.mean(ctrl_last):.2f}±{np.std(ctrl_last):.2f}  ratio={np.mean(exp_last)/np.mean(ctrl_last):.4f}")

# Check if the last token has unusual activation patterns for experimental
print("\nLast token cosine similarity to second-to-last token:")
for ln in ["layer_22", "layer_23", "final_norm"]:
    exp_cos = []
    ctrl_cos = []
    for e in data["experimental"]:
        v1 = e[ln][-1]
        v2 = e[ln][-2]
        cos = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-10)
        exp_cos.append(cos)
    for e in data["control"]:
        v1 = e[ln][-1]
        v2 = e[ln][-2]
        cos = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-10)
        ctrl_cos.append(cos)
    print(f"  {ln:<15} exp={np.mean(exp_cos):.4f}  ctrl={np.mean(ctrl_cos):.4f}")


print("\n" + "="*80)
print("DEEP DIVE 6: CATEGORY BREAKDOWN")
print("="*80)

# Are the effects driven by category differences?
cats = defaultdict(lambda: defaultdict(list))
for s in ["experimental", "control"]:
    for e in data[s]:
        for ln in ["embed", "layer_6", "layer_12", "layer_23", "final_norm"]:
            norm = np.mean(np.linalg.norm(e[ln], axis=1))
            cats[(e["category"], e["set"])][ln].append(norm)

print("\nMean L2 norms by category and set:")
for (cat, s), layer_dict in sorted(cats.items()):
    n = len(list(layer_dict.values())[0])
    print(f"\n  {cat} ({s}, n={n}):")
    for ln in ["embed", "layer_6", "layer_12", "layer_23", "final_norm"]:
        vals = layer_dict[ln]
        print(f"    {ln:<15} {np.mean(vals):.2f} ± {np.std(vals):.2f}")


print("\n" + "="*80)
print("DEEP DIVE 7: RESIDUAL CONTRIBUTION RATIO AT LAYER 23")
print("="*80)

# Layer 23 has residual/output ratio > 1 for experimental. This means
# layer 23's contribution EXCEEDS its output — it's partially canceling previous layers.
for e in data["experimental"][:5]:
    l22 = e["layer_22"]
    l23 = e["layer_23"]
    residual = l23 - l22
    res_norm = np.mean(np.linalg.norm(residual, axis=1))
    out_norm = np.mean(np.linalg.norm(l23, axis=1))
    l22_norm = np.mean(np.linalg.norm(l22, axis=1))
    # How much does layer 23 change the direction?
    cos_sims = []
    for pos in range(l22.shape[0]):
        c = np.dot(l22[pos], l23[pos]) / (np.linalg.norm(l22[pos]) * np.linalg.norm(l23[pos]) + 1e-10)
        cos_sims.append(c)
    avg_cos = np.mean(cos_sims)
    print(f"  prompt {e['prompt_idx']}: l22_norm={l22_norm:.2f}, l23_norm={out_norm:.2f}, "
          f"residual={res_norm:.2f}, cos(l22,l23)={avg_cos:.4f}")

for e in data["control"][:5]:
    l22 = e["layer_22"]
    l23 = e["layer_23"]
    residual = l23 - l22
    res_norm = np.mean(np.linalg.norm(residual, axis=1))
    out_norm = np.mean(np.linalg.norm(l23, axis=1))
    l22_norm = np.mean(np.linalg.norm(l22, axis=1))
    cos_sims = []
    for pos in range(l22.shape[0]):
        c = np.dot(l22[pos], l23[pos]) / (np.linalg.norm(l22[pos]) * np.linalg.norm(l23[pos]) + 1e-10)
        cos_sims.append(c)
    avg_cos = np.mean(cos_sims)
    print(f"  prompt {e['prompt_idx']}: l22_norm={l22_norm:.2f}, l23_norm={out_norm:.2f}, "
          f"residual={res_norm:.2f}, cos(l22,l23)={avg_cos:.4f}")


print("\n" + "="*80)
print("DEEP DIVE 8: EFFECTIVE DIMENSIONALITY (PARTICIPATION RATIO)")
print("="*80)
print("PR = (sum |x_i|^2)^2 / (sum |x_i|^4) — measures effective #dimensions used.")

def participation_ratio(entries, layer_name):
    prs = []
    for e in entries:
        mean_act = np.mean(e[layer_name], axis=0)
        sq = mean_act ** 2
        s2 = np.sum(sq)
        s4 = np.sum(sq ** 2)
        if s4 > 0:
            pr = s2 ** 2 / s4
        else:
            pr = 0
        prs.append(pr)
    return np.mean(prs), np.std(prs)

print(f"\n{'Layer':<15} {'Exp PR':>12} {'Ctrl PR':>12} {'Ratio':>10} {'(out of 1024)'}")
for ln in LAYER_NAMES:
    em, es = participation_ratio(data["experimental"], ln)
    cm, cs = participation_ratio(data["control"], ln)
    ratio = em / cm if cm > 0 else float('inf')
    print(f"  {ln:<13} {em:>10.1f}±{es:.1f} {cm:>10.1f}±{cs:.1f} {ratio:>10.4f}")

print("\nAll deep-dive analyses complete.")
