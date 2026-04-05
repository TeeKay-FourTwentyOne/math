"""
Iteration 4: Characterize the 70 "unflippable" control prompts.

These prompts have ALL 256 margins > 0.1, meaning even the model's most
uncertain moment is still fairly confident. What makes them special?

Also characterize the opposite extreme: the most flippable prompts.
"""

import torch
import numpy as np
import json
import os
from collections import Counter, defaultdict
from pathlib import Path

LOGIT_DIR = "C:/Projects/math/non_deterministic/experiment_0/logit_traces"
PROMPTS_PATH = "C:/Projects/math/non_deterministic/prompts/candidates_2000.jsonl"
CTRL_PATH = "C:/Projects/math/non_deterministic/prompts/control_200.jsonl"
EXP_PATH = "C:/Projects/math/non_deterministic/prompts/experimental_200.jsonl"
RUN1_PATH = "C:/Projects/math/non_deterministic/experiment_1/outputs/h_nvidia_run1.jsonl"

from transformers import AutoTokenizer
tokenizer = AutoTokenizer.from_pretrained("EleutherAI/pythia-410m-deduped")

# Load all prompts
prompts = {}
with open(PROMPTS_PATH) as f:
    for line in f:
        p = json.loads(line)
        prompts[p["idx"]] = p

# Load control set
ctrl_idxs = set()
ctrl_prompts = {}
with open(CTRL_PATH, encoding="utf-8") as f:
    for line in f:
        d = json.loads(line)
        ctrl_idxs.add(d["prompt_idx"])
        ctrl_prompts[d["prompt_idx"]] = d

# Load experimental set
exp_idxs = set()
with open(EXP_PATH, encoding="utf-8") as f:
    for line in f:
        d = json.loads(line)
        exp_idxs.add(d["prompt_idx"])

# Load generated text from experiment 1A
run1_data = {}
with open(RUN1_PATH, encoding="utf-8") as f:
    for line in f:
        d = json.loads(line)
        run1_data[d["prompt_idx"]] = d

# ============================================================
# Identify unflippable and most-flippable prompts
# ============================================================
print("Analyzing all prompts...", flush=True)

prompt_stats = []
for i in range(1999):
    path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
    if not os.path.exists(path):
        continue
    data = torch.load(path, map_location="cpu", weights_only=True)
    top5_vals = data["top5_values"]
    margins = (top5_vals[:, 0] - top5_vals[:, 1]).numpy()

    prompt_stats.append({
        "idx": i,
        "category": prompts[i]["category"],
        "text": prompts[i]["text"],
        "min_margin": float(margins.min()),
        "mean_margin": float(margins.mean()),
        "median_margin": float(np.median(margins)),
        "max_margin": float(margins.max()),
        "std_margin": float(margins.std()),
        "fragility_1_0": float((margins < 1.0).mean()),
        "fragility_0_1": float((margins < 0.1).mean()),
        "n_tokens": len(margins),
        "in_ctrl": i in ctrl_idxs,
        "in_exp": i in exp_idxs,
        "text_len": len(prompts[i]["text"]),
        "word_count": len(prompts[i]["text"].split()),
    })

# Identify unflippable (min_margin > 0.1 among control prompts)
unflippable = [p for p in prompt_stats if p["in_ctrl"] and p["min_margin"] > 0.1]
flippable_ctrl = [p for p in prompt_stats if p["in_ctrl"] and p["min_margin"] <= 0.1]

# Most flippable (experimental with highest fragility)
most_flippable = sorted([p for p in prompt_stats if p["in_exp"]],
                        key=lambda x: x["min_margin"])[:70]

print(f"Unflippable control prompts: {len(unflippable)}")
print(f"Flippable control prompts: {len(flippable_ctrl)}")

# ============================================================
# PART 1: Category distribution
# ============================================================
print("\n" + "=" * 70)
print("PART 1: Category distribution")
print("=" * 70)

# Control set baseline
ctrl_cats = Counter(p["category"] for p in prompt_stats if p["in_ctrl"])
unflip_cats = Counter(p["category"] for p in unflippable)
flip_ctrl_cats = Counter(p["category"] for p in flippable_ctrl)
most_flip_cats = Counter(p["category"] for p in most_flippable)

print(f"\n{'Category':15s} {'Control':>8s} {'Unflip':>8s} {'%Unflip':>8s} {'FlipCtrl':>8s} {'MostFlip':>8s}")
print("-" * 55)
for cat in sorted(ctrl_cats.keys()):
    n_ctrl = ctrl_cats[cat]
    n_unflip = unflip_cats.get(cat, 0)
    n_fctrl = flip_ctrl_cats.get(cat, 0)
    n_mflip = most_flip_cats.get(cat, 0)
    pct = n_unflip / n_ctrl * 100 if n_ctrl > 0 else 0
    print(f"{cat:15s} {n_ctrl:8d} {n_unflip:8d} {pct:7.1f}% {n_fctrl:8d} {n_mflip:8d}")

# ============================================================
# PART 2: Margin profile comparison
# ============================================================
print("\n" + "=" * 70)
print("PART 2: Margin statistics comparison")
print("=" * 70)

groups = {
    "Unflippable (70)": unflippable,
    "Flippable ctrl (130)": flippable_ctrl,
    "Most flippable exp (70)": most_flippable,
}

print(f"\n{'Group':25s} {'MinMarg':>8s} {'MeanMarg':>9s} {'MedMarg':>8s} {'MaxMarg':>8s} {'StdMarg':>8s}")
print("-" * 70)
for name, group in groups.items():
    min_m = np.mean([p["min_margin"] for p in group])
    mean_m = np.mean([p["mean_margin"] for p in group])
    med_m = np.mean([p["median_margin"] for p in group])
    max_m = np.mean([p["max_margin"] for p in group])
    std_m = np.mean([p["std_margin"] for p in group])
    print(f"{name:25s} {min_m:8.3f} {mean_m:9.3f} {med_m:8.3f} {max_m:8.3f} {std_m:8.3f}")

# ============================================================
# PART 3: Prompt text characteristics
# ============================================================
print("\n" + "=" * 70)
print("PART 3: Prompt characteristics")
print("=" * 70)

for name, group in groups.items():
    text_lens = [p["text_len"] for p in group]
    word_counts = [p["word_count"] for p in group]
    print(f"\n{name}:")
    print(f"  Text length: mean={np.mean(text_lens):.0f} median={np.median(text_lens):.0f} "
          f"min={min(text_lens)} max={max(text_lens)}")
    print(f"  Word count:  mean={np.mean(word_counts):.1f} median={np.median(word_counts):.0f}")

# ============================================================
# PART 4: What do unflippable prompts look like?
# ============================================================
print("\n" + "=" * 70)
print("PART 4: Unflippable prompt examples (sorted by min_margin, highest first)")
print("=" * 70)

unflippable.sort(key=lambda x: -x["min_margin"])
for i, p in enumerate(unflippable[:15]):
    text = p["text"][:75]
    print(f"  #{i+1:2d} [{p['category']:12s}] min_margin={p['min_margin']:.3f} "
          f"mean={p['mean_margin']:.2f}")
    print(f"      {text}")

print(f"\n--- Bottom of unflippable (smallest min_margin, barely unflippable) ---")
unflippable.sort(key=lambda x: x["min_margin"])
for i, p in enumerate(unflippable[:10]):
    text = p["text"][:75]
    print(f"  #{i+1:2d} [{p['category']:12s}] min_margin={p['min_margin']:.3f} "
          f"mean={p['mean_margin']:.2f}")
    print(f"      {text}")

# ============================================================
# PART 5: Generated text analysis
# ============================================================
print("\n" + "=" * 70)
print("PART 5: Generated text analysis")
print("=" * 70)

# For prompts in the experiment 1A run, compare generated text
def get_gen_text(idx):
    if idx in run1_data:
        return run1_data[idx].get("generated_text", "")
    return ""

# Measure generated text properties
for name, group in groups.items():
    gen_lens = []
    has_newlines = 0
    has_repetition = 0
    starts_with = Counter()

    for p in group:
        gt = get_gen_text(p["idx"])
        if gt:
            gen_lens.append(len(gt))
            if "\n\n" in gt[10:]:
                has_newlines += 1
            # Check for repetition (same 20-char substring appearing 3+ times)
            if len(gt) > 60:
                for start in range(0, len(gt) - 20, 10):
                    chunk = gt[start:start+20]
                    if gt.count(chunk) >= 3:
                        has_repetition += 1
                        break
            # First meaningful word
            stripped = gt.lstrip("\n ").split()
            if stripped:
                starts_with[stripped[0]] += 1

    n_with_gen = sum(1 for p in group if get_gen_text(p["idx"]))
    print(f"\n{name} (with generated text: {n_with_gen}):")
    if gen_lens:
        print(f"  Gen text length: mean={np.mean(gen_lens):.0f}")
    if n_with_gen > 0:
        print(f"  Has paragraph breaks: {has_newlines}/{n_with_gen} ({has_newlines/n_with_gen*100:.0f}%)")
        print(f"  Has repetition: {has_repetition}/{n_with_gen} ({has_repetition/n_with_gen*100:.0f}%)")
        print(f"  Top starting words: {starts_with.most_common(5)}")

# ============================================================
# PART 6: Margin CURVE shape — do unflippable prompts start confident and stay confident?
# ============================================================
print("\n" + "=" * 70)
print("PART 6: Margin curve shape across generation")
print("=" * 70)

# For each group, compute average margin at each position window
windows = [(0, 7), (8, 15), (16, 31), (32, 63), (64, 127), (128, 191), (192, 255)]

for name, group in groups.items():
    print(f"\n{name} — average margin by window:")
    for lo, hi in windows:
        win_margins = []
        for p in group:
            path = os.path.join(LOGIT_DIR, f"prompt_{p['idx']:04d}.pt")
            data = torch.load(path, map_location="cpu", weights_only=True)
            margins = (data["top5_values"][:, 0] - data["top5_values"][:, 1]).numpy()
            seg = margins[lo:hi+1]
            if len(seg) > 0:
                win_margins.append(seg.mean())
        if win_margins:
            m = np.mean(win_margins)
            s = np.std(win_margins)
            print(f"  pos {lo:3d}-{hi:3d}: mean={m:6.2f} std={s:5.2f}")

# ============================================================
# PART 7: Is there a simple unflippability predictor?
# ============================================================
print("\n" + "=" * 70)
print("PART 7: Simple predictors of unflippability")
print("=" * 70)

# Among control prompts, can we predict unflippability from:
# a) category
# b) prompt length
# c) margin at position 2 (the sentence-start)
# d) margin at position 0

ctrl_data = [p for p in prompt_stats if p["in_ctrl"]]

# Logistic-style analysis: what separates unflippable from flippable?
for feature_name, feature_fn in [
    ("prompt_text_length", lambda p: p["text_len"]),
    ("prompt_word_count", lambda p: p["word_count"]),
    ("mean_margin", lambda p: p["mean_margin"]),
    ("median_margin", lambda p: p["median_margin"]),
]:
    vals_unflip = [feature_fn(p) for p in unflippable]
    vals_flip = [feature_fn(p) for p in flippable_ctrl]

    if vals_unflip and vals_flip:
        u_mean = np.mean(vals_unflip)
        f_mean = np.mean(vals_flip)
        # Simple effect size
        pooled_std = np.sqrt((np.var(vals_unflip) + np.var(vals_flip)) / 2)
        d = (u_mean - f_mean) / pooled_std if pooled_std > 0 else 0
        print(f"  {feature_name:25s}: unflip={u_mean:.2f}  flip={f_mean:.2f}  Cohen's d={d:.2f}")

# ============================================================
# PART 8: The unflippable margin at position 2
# ============================================================
print("\n" + "=" * 70)
print("PART 8: Position-2 margin for unflippable vs flippable control")
print("=" * 70)

for name, group in [("Unflippable", unflippable), ("Flippable ctrl", flippable_ctrl)]:
    pos2_margins = []
    for p in group:
        path = os.path.join(LOGIT_DIR, f"prompt_{p['idx']:04d}.pt")
        data = torch.load(path, map_location="cpu", weights_only=True)
        margins = (data["top5_values"][:, 0] - data["top5_values"][:, 1]).numpy()
        if len(margins) > 2:
            pos2_margins.append(margins[2])
    m = np.array(pos2_margins)
    print(f"  {name:20s}: mean={m.mean():.3f} median={np.median(m):.3f} "
          f"min={m.min():.3f} max={m.max():.3f} frac<1.0={float((m < 1.0).mean()):.3f}")

# ============================================================
# PART 9: What specific prompts make up the most unflippable?
# ============================================================
print("\n" + "=" * 70)
print("PART 9: The 10 MOST unflippable prompts (highest min_margin)")
print("=" * 70)

all_sorted = sorted(prompt_stats, key=lambda x: -x["min_margin"])
for i, p in enumerate(all_sorted[:10]):
    gt = get_gen_text(p["idx"])
    gen_preview = gt[:100].replace("\n", "\\n") if gt else "(no gen text)"
    set_tag = "CTRL" if p["in_ctrl"] else ("EXP" if p["in_exp"] else "---")
    print(f"\n  #{i+1} idx={p['idx']} [{p['category']}|{set_tag}] "
          f"min_margin={p['min_margin']:.3f} mean_margin={p['mean_margin']:.2f}")
    print(f"    Prompt: {p['text'][:80]}")
    print(f"    Gen: {gen_preview}")
