"""
Comprehensive activation analysis for Experiment 1B.
Pythia-410M: 24 layers, hidden_size=1024, FP32.
50 prompts: 25 experimental (high fragility), 25 control (low fragility).
"""

import torch
import numpy as np
import json
import os
from pathlib import Path
from collections import defaultdict

DATA_DIR = Path("C:/Projects/math/non_deterministic/experiment_1/activations/h_nvidia")
META_FILE = Path("C:/Projects/math/non_deterministic/experiment_1/activations/activation_subset.jsonl")
OUT_DIR = Path("C:/Projects/math/non_deterministic/experiment_1/analysis_output")
OUT_DIR.mkdir(exist_ok=True)

LAYER_NAMES = ["embed"] + [f"layer_{i}" for i in range(24)] + ["final_norm"]
N_LAYERS = len(LAYER_NAMES)  # 26 total: embed + 24 layers + final_norm

# ─── Load metadata ───
meta = {}
with open(META_FILE) as f:
    for line in f:
        if line.strip():
            d = json.loads(line)
            meta[d["prompt_idx"]] = d

# ─── Load all activations ───
print("Loading activation files...")
data = {"experimental": [], "control": []}
for pt_file in sorted(DATA_DIR.glob("prompt_*.pt")):
    d = torch.load(pt_file, map_location="cpu", weights_only=False)
    prompt_idx = d["prompt_idx"]
    set_label = d["set"]
    entry = {
        "prompt_idx": prompt_idx,
        "input_length": d["input_length"],
        "total_length": d["total_length"],
        "category": d["category"],
        "set": set_label,
        "fragility": meta[prompt_idx]["fragility_score"],
        "text": meta[prompt_idx]["text"],
    }
    # Store per-layer activations as numpy
    for ln in LAYER_NAMES:
        entry[ln] = d[ln].numpy()
    data[set_label].append(entry)

print(f"Loaded: {len(data['experimental'])} experimental, {len(data['control'])} control")

# Quick summary
for s in ["experimental", "control"]:
    frags = [e["fragility"] for e in data[s]]
    lengths = [e["total_length"] for e in data[s]]
    input_lengths = [e["input_length"] for e in data[s]]
    cats = defaultdict(int)
    for e in data[s]:
        cats[e["category"]] += 1
    print(f"\n{s.upper()}:")
    print(f"  Fragility: mean={np.mean(frags):.4f}, range=[{min(frags):.4f}, {max(frags):.4f}]")
    print(f"  Total length: mean={np.mean(lengths):.1f}, range=[{min(lengths)}, {max(lengths)}]")
    print(f"  Input length: mean={np.mean(input_lengths):.1f}, range=[{min(input_lengths)}, {max(input_lengths)}]")
    print(f"  Categories: {dict(cats)}")


# ═══════════════════════════════════════════════════════════════════
# ANALYSIS 1: Activation Growth Profile
# ═══════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("ANALYSIS 1: ACTIVATION GROWTH PROFILE")
print("="*80)

def compute_mean_l2_norm_per_layer(entries):
    """For each layer, compute L2 norm per position, then average across positions and prompts."""
    layer_norms = {ln: [] for ln in LAYER_NAMES}
    for e in entries:
        for ln in LAYER_NAMES:
            act = e[ln]  # (seq_len, 1024)
            l2_per_pos = np.linalg.norm(act, axis=1)  # (seq_len,)
            mean_l2 = np.mean(l2_per_pos)
            layer_norms[ln].append(mean_l2)
    return {ln: (np.mean(layer_norms[ln]), np.std(layer_norms[ln])) for ln in LAYER_NAMES}

exp_norms = compute_mean_l2_norm_per_layer(data["experimental"])
ctrl_norms = compute_mean_l2_norm_per_layer(data["control"])

print(f"\n{'Layer':<15} {'Experimental':>15} {'Control':>15} {'Ratio (E/C)':>12} {'Jump from prev':>15}")
prev_exp, prev_ctrl = None, None
for ln in LAYER_NAMES:
    em, es = exp_norms[ln]
    cm, cs = ctrl_norms[ln]
    ratio = em / cm if cm > 0 else float('inf')
    if prev_exp is not None:
        jump_exp = em - prev_exp
        jump_ctrl = cm - prev_ctrl
        print(f"  {ln:<13} {em:>12.2f} ±{es:<5.2f} {cm:>10.2f} ±{cs:<5.2f} {ratio:>10.4f}   exp:{jump_exp:+.2f} ctrl:{jump_ctrl:+.2f}")
    else:
        print(f"  {ln:<13} {em:>12.2f} ±{es:<5.2f} {cm:>10.2f} ±{cs:<5.2f} {ratio:>10.4f}")
    prev_exp, prev_ctrl = em, cm

# Compute growth ratios
print("\nGrowth pattern analysis:")
exp_vals = [exp_norms[ln][0] for ln in LAYER_NAMES]
ctrl_vals = [ctrl_norms[ln][0] for ln in LAYER_NAMES]
for name, vals in [("Experimental", exp_vals), ("Control", ctrl_vals)]:
    total_growth = vals[-1] / vals[0] if vals[0] > 0 else 0
    mid_growth = vals[13] / vals[0] if vals[0] > 0 else 0
    print(f"  {name}: embed-to-final ratio = {total_growth:.2f}x, embed-to-mid ratio = {mid_growth:.2f}x")
    # Check linearity: compute R^2 of linear fit
    x = np.arange(len(vals))
    coeffs = np.polyfit(x, vals, 1)
    predicted = np.polyval(coeffs, x)
    ss_res = np.sum((vals - predicted)**2)
    ss_tot = np.sum((vals - np.mean(vals))**2)
    r2_linear = 1 - ss_res / ss_tot
    # Check exponential: fit log
    log_vals = np.log(np.array(vals) + 1e-10)
    coeffs_exp = np.polyfit(x, log_vals, 1)
    predicted_exp = np.polyval(coeffs_exp, x)
    ss_res_exp = np.sum((log_vals - predicted_exp)**2)
    ss_tot_exp = np.sum((log_vals - np.mean(log_vals))**2)
    r2_exp = 1 - ss_res_exp / ss_tot_exp
    print(f"    Linear fit R² = {r2_linear:.4f}, Exponential fit R² = {r2_exp:.4f}")


# ═══════════════════════════════════════════════════════════════════
# ANALYSIS 2: Activation Variance
# ═══════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("ANALYSIS 2: ACTIVATION VARIANCE ACROSS POSITIONS")
print("="*80)

def compute_variance_per_layer(entries):
    """Variance of activations across positions (average over dimensions, then over prompts)."""
    layer_vars = {ln: [] for ln in LAYER_NAMES}
    for e in entries:
        for ln in LAYER_NAMES:
            act = e[ln]  # (seq_len, 1024)
            # Variance across positions for each dimension, then mean across dimensions
            var_per_dim = np.var(act, axis=0)  # (1024,)
            mean_var = np.mean(var_per_dim)
            layer_vars[ln].append(mean_var)
    return {ln: (np.mean(layer_vars[ln]), np.std(layer_vars[ln])) for ln in LAYER_NAMES}

exp_vars = compute_variance_per_layer(data["experimental"])
ctrl_vars = compute_variance_per_layer(data["control"])

print(f"\n{'Layer':<15} {'Exp Variance':>15} {'Ctrl Variance':>15} {'Ratio (E/C)':>12}")
for ln in LAYER_NAMES:
    em, es = exp_vars[ln]
    cm, cs = ctrl_vars[ln]
    ratio = em / cm if cm > 0 else float('inf')
    print(f"  {ln:<13} {em:>13.4f} ±{es:<8.4f} {cm:>11.4f} ±{cs:<8.4f} {ratio:>10.4f}")


# ═══════════════════════════════════════════════════════════════════
# ANALYSIS 3: Layer-by-Layer Divergence
# ═══════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("ANALYSIS 3: DIVERGENCE BETWEEN EXPERIMENTAL AND CONTROL")
print("="*80)

def compute_mean_activation_vector(entries):
    """Compute mean activation vector per layer (averaged across positions, then across prompts)."""
    layer_means = {ln: [] for ln in LAYER_NAMES}
    for e in entries:
        for ln in LAYER_NAMES:
            act = e[ln]  # (seq_len, 1024)
            mean_vec = np.mean(act, axis=0)  # (1024,)
            layer_means[ln].append(mean_vec)
    return {ln: np.mean(layer_means[ln], axis=0) for ln in LAYER_NAMES}

exp_means = compute_mean_activation_vector(data["experimental"])
ctrl_means = compute_mean_activation_vector(data["control"])

print(f"\n{'Layer':<15} {'L2 Distance':>12} {'Cosine Sim':>12} {'Rel. Distance':>14}")
for ln in LAYER_NAMES:
    diff = exp_means[ln] - ctrl_means[ln]
    l2_dist = np.linalg.norm(diff)
    # Cosine similarity
    dot = np.dot(exp_means[ln], ctrl_means[ln])
    norm_e = np.linalg.norm(exp_means[ln])
    norm_c = np.linalg.norm(ctrl_means[ln])
    cos_sim = dot / (norm_e * norm_c + 1e-10)
    rel_dist = l2_dist / (0.5 * (norm_e + norm_c) + 1e-10)
    print(f"  {ln:<13} {l2_dist:>12.4f} {cos_sim:>12.6f} {rel_dist:>14.6f}")

# Also compute per-prompt pairwise distances at each layer
print("\nPer-prompt mean L2 distance between group centroids and individual prompts:")
print("(How scattered are prompts within each group?)")
for s, entries, means in [("Experimental", data["experimental"], exp_means),
                           ("Control", data["control"], ctrl_means)]:
    print(f"\n  {s}:")
    for ln in ["embed", "layer_0", "layer_6", "layer_12", "layer_18", "layer_23", "final_norm"]:
        dists = []
        for e in entries:
            mean_vec = np.mean(e[ln], axis=0)
            dists.append(np.linalg.norm(mean_vec - means[ln]))
        print(f"    {ln:<13} mean intra-group distance: {np.mean(dists):.4f} ± {np.std(dists):.4f}")


# ═══════════════════════════════════════════════════════════════════
# ANALYSIS 4: Position-Dependent Patterns
# ═══════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("ANALYSIS 4: POSITION-DEPENDENT PATTERNS (INPUT vs GENERATED)")
print("="*80)

def compute_positional_norms(entries, part):
    """Compute L2 norms for input or generated portions."""
    layer_norms = {ln: [] for ln in LAYER_NAMES}
    for e in entries:
        inp_len = e["input_length"]
        tot_len = e["total_length"]
        for ln in LAYER_NAMES:
            act = e[ln]
            if part == "input":
                portion = act[:inp_len]
            elif part == "generated":
                portion = act[inp_len:]
            else:  # "first5" or "last5"
                if part == "first5":
                    portion = act[:min(5, act.shape[0])]
                else:
                    portion = act[max(0, act.shape[0]-5):]
            if portion.shape[0] > 0:
                l2_per_pos = np.linalg.norm(portion, axis=1)
                layer_norms[ln].append(np.mean(l2_per_pos))
    return {ln: (np.mean(layer_norms[ln]), np.std(layer_norms[ln])) for ln in LAYER_NAMES}

for part_name in ["input", "generated"]:
    print(f"\n--- {part_name.upper()} portion ---")
    print(f"{'Layer':<15} {'Experimental':>15} {'Control':>15} {'Ratio':>10}")
    exp_pn = compute_positional_norms(data["experimental"], part_name)
    ctrl_pn = compute_positional_norms(data["control"], part_name)
    for ln in ["embed", "layer_0", "layer_6", "layer_12", "layer_18", "layer_23", "final_norm"]:
        em, _ = exp_pn[ln]
        cm, _ = ctrl_pn[ln]
        ratio = em / cm if cm > 0 else float('inf')
        print(f"  {ln:<13} {em:>13.2f} {cm:>13.2f} {ratio:>10.4f}")

# Ratio of generated-to-input norms
print("\n--- Ratio of GENERATED / INPUT norms (amplification of generated tokens) ---")
print(f"{'Layer':<15} {'Experimental':>15} {'Control':>15}")
exp_inp = compute_positional_norms(data["experimental"], "input")
exp_gen = compute_positional_norms(data["experimental"], "generated")
ctrl_inp = compute_positional_norms(data["control"], "input")
ctrl_gen = compute_positional_norms(data["control"], "generated")
for ln in LAYER_NAMES:
    er = exp_gen[ln][0] / exp_inp[ln][0] if exp_inp[ln][0] > 0 else 0
    cr = ctrl_gen[ln][0] / ctrl_inp[ln][0] if ctrl_inp[ln][0] > 0 else 0
    print(f"  {ln:<13} {er:>13.4f} {cr:>13.4f}")


# ═══════════════════════════════════════════════════════════════════
# ANALYSIS 5: Residual Stream Analysis
# ═══════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("ANALYSIS 5: RESIDUAL STREAM — LAYER CONTRIBUTIONS")
print("="*80)

def compute_residual_contributions(entries):
    """Compute ||layer_k - layer_(k-1)|| for each layer."""
    ordered = LAYER_NAMES  # embed, layer_0, ..., layer_23, final_norm
    contribs = {ordered[i]: [] for i in range(1, len(ordered))}
    for e in entries:
        for i in range(1, len(ordered)):
            curr = e[ordered[i]]
            prev = e[ordered[i-1]]
            diff = curr - prev
            # Mean L2 norm of residual across positions
            l2_per_pos = np.linalg.norm(diff, axis=1)
            contribs[ordered[i]].append(np.mean(l2_per_pos))
    return {k: (np.mean(v), np.std(v)) for k, v in contribs.items()}

exp_contribs = compute_residual_contributions(data["experimental"])
ctrl_contribs = compute_residual_contributions(data["control"])

print(f"\n{'Layer':<15} {'Exp Contribution':>18} {'Ctrl Contribution':>18} {'Ratio (E/C)':>12}")
for ln in LAYER_NAMES[1:]:
    em, es = exp_contribs[ln]
    cm, cs = ctrl_contribs[ln]
    ratio = em / cm if cm > 0 else float('inf')
    print(f"  {ln:<13} {em:>14.4f} ±{es:<6.4f} {cm:>14.4f} ±{cs:<6.4f} {ratio:>10.4f}")

# Relative contribution (contribution / layer norm)
print(f"\n--- Relative contribution (residual magnitude / layer output magnitude) ---")
print(f"{'Layer':<15} {'Experimental':>15} {'Control':>15}")
for ln in LAYER_NAMES[1:]:
    exp_rel = exp_contribs[ln][0] / exp_norms[ln][0] if exp_norms[ln][0] > 0 else 0
    ctrl_rel = ctrl_contribs[ln][0] / ctrl_norms[ln][0] if ctrl_norms[ln][0] > 0 else 0
    print(f"  {ln:<13} {exp_rel:>13.4f} {ctrl_rel:>13.4f}")

# Top 5 contributing layers
print("\nTop 5 contributing layers by absolute residual magnitude:")
for name, contribs in [("Experimental", exp_contribs), ("Control", ctrl_contribs)]:
    sorted_layers = sorted(contribs.items(), key=lambda x: x[1][0], reverse=True)
    print(f"  {name}: ", end="")
    for ln, (m, s) in sorted_layers[:5]:
        print(f"{ln}({m:.2f})", end="  ")
    print()


# ═══════════════════════════════════════════════════════════════════
# ANALYSIS 6: Entropy / Concentration
# ═══════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("ANALYSIS 6: ACTIVATION CONCENTRATION")
print("="*80)

def compute_concentration(entries, top_frac=0.1):
    """Fraction of total L2 norm² in top 10% of dimensions."""
    layer_conc = {ln: [] for ln in LAYER_NAMES}
    layer_gini = {ln: [] for ln in LAYER_NAMES}
    for e in entries:
        for ln in LAYER_NAMES:
            act = e[ln]  # (seq_len, 1024)
            # Mean activation vector across positions
            mean_act = np.mean(act, axis=0)  # (1024,)
            sq = mean_act ** 2
            total = np.sum(sq)
            if total < 1e-12:
                layer_conc[ln].append(top_frac)
                layer_gini[ln].append(0.0)
                continue
            sorted_sq = np.sort(sq)[::-1]
            k = int(np.ceil(top_frac * len(sq)))
            top_mass = np.sum(sorted_sq[:k]) / total
            layer_conc[ln].append(top_mass)

            # Gini coefficient of |activations|
            abs_vals = np.abs(mean_act)
            sorted_abs = np.sort(abs_vals)
            n = len(sorted_abs)
            cumsum = np.cumsum(sorted_abs)
            gini = (2 * np.sum((np.arange(1, n+1) * sorted_abs)) / (n * np.sum(sorted_abs))) - (n + 1) / n
            layer_gini[ln].append(gini)

    return ({ln: (np.mean(layer_conc[ln]), np.std(layer_conc[ln])) for ln in LAYER_NAMES},
            {ln: (np.mean(layer_gini[ln]), np.std(layer_gini[ln])) for ln in LAYER_NAMES})

exp_conc, exp_gini = compute_concentration(data["experimental"])
ctrl_conc, ctrl_gini = compute_concentration(data["control"])

print(f"\nFraction of L2 norm² in top 10% of dimensions (uniform = 0.10):")
print(f"{'Layer':<15} {'Experimental':>15} {'Control':>15} {'Ratio':>10}")
for ln in LAYER_NAMES:
    em, _ = exp_conc[ln]
    cm, _ = ctrl_conc[ln]
    ratio = em / cm if cm > 0 else float('inf')
    print(f"  {ln:<13} {em:>13.4f} {cm:>13.4f} {ratio:>10.4f}")

print(f"\nGini coefficient of |activations| (0=uniform, 1=concentrated):")
print(f"{'Layer':<15} {'Experimental':>15} {'Control':>15} {'Diff (E-C)':>12}")
for ln in LAYER_NAMES:
    em, _ = exp_gini[ln]
    cm, _ = ctrl_gini[ln]
    diff = em - cm
    print(f"  {ln:<13} {em:>13.4f} {cm:>13.4f} {diff:>+12.4f}")


# ═══════════════════════════════════════════════════════════════════
# ADDITIONAL: Per-position activation norm profiles
# ═══════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("ADDITIONAL: PER-POSITION NORM PROFILE AT KEY LAYERS")
print("="*80)

def position_norm_profile(entries, layer_name, max_pos=30):
    """L2 norm at each position (first max_pos positions), averaged across prompts."""
    profiles = []
    for e in entries:
        act = e[layer_name]
        norms = np.linalg.norm(act[:max_pos], axis=1)
        profiles.append(norms)
    # Pad to same length
    min_len = min(len(p) for p in profiles)
    profiles = [p[:min_len] for p in profiles]
    return np.mean(profiles, axis=0), np.std(profiles, axis=0)

for ln in ["embed", "layer_0", "layer_12", "layer_23", "final_norm"]:
    print(f"\n  {ln} — first 20 positions:")
    exp_prof, exp_std = position_norm_profile(data["experimental"], ln, 20)
    ctrl_prof, ctrl_std = position_norm_profile(data["control"], ln, 20)
    for pos in range(min(20, len(exp_prof))):
        print(f"    pos {pos:>2}: exp={exp_prof[pos]:.3f}  ctrl={ctrl_prof[pos]:.3f}  diff={exp_prof[pos]-ctrl_prof[pos]:+.3f}")


# ═══════════════════════════════════════════════════════════════════
# STATISTICAL TESTS
# ═══════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("STATISTICAL SIGNIFICANCE: Mann-Whitney U tests")
print("="*80)

from scipy import stats

# Test: do experimental prompts have higher activation norms at each layer?
print(f"\n{'Layer':<15} {'U statistic':>12} {'p-value':>12} {'Significant?':>14}")
for ln in LAYER_NAMES:
    exp_vals = []
    ctrl_vals = []
    for e in data["experimental"]:
        l2_per_pos = np.linalg.norm(e[ln], axis=1)
        exp_vals.append(np.mean(l2_per_pos))
    for e in data["control"]:
        l2_per_pos = np.linalg.norm(e[ln], axis=1)
        ctrl_vals.append(np.mean(l2_per_pos))
    u_stat, p_val = stats.mannwhitneyu(exp_vals, ctrl_vals, alternative='two-sided')
    sig = "YES" if p_val < 0.05 else "no"
    print(f"  {ln:<13} {u_stat:>12.1f} {p_val:>12.6f} {sig:>14}")

print("\nDone. All analyses complete.")
