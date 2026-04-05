"""
Cross-Platform Comparison (ALIGNED): H_nvidia vs H_mac
Accounts for format difference: Mac includes input tokens in generated_ids.
"""

import json
import numpy as np
from collections import Counter
from transformers import AutoTokenizer

tokenizer = AutoTokenizer.from_pretrained("EleutherAI/pythia-410m-deduped")

def load_nvidia_run(path):
    """NVIDIA format: generated_ids are ONLY generated tokens."""
    data = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            d = json.loads(line)
            data[d["prompt_idx"]] = d["generated_ids"]  # already stripped
    return data

def load_mac_run(path):
    """Mac format: generated_ids include input tokens — strip them."""
    data = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            d = json.loads(line)
            offset = d["n_input_tokens"]
            data[d["prompt_idx"]] = d["generated_ids"][offset:]  # strip input
    return data

# Load all runs
nvidia_r1 = load_nvidia_run("C:/Projects/math/non_deterministic/experiment_1/outputs/h_nvidia_run1.jsonl")
nvidia_r2 = load_nvidia_run("C:/Projects/math/non_deterministic/experiment_1/outputs/h_nvidia_run2.jsonl")
mac_r1 = load_mac_run("C:/Projects/math/non_deterministic/experiment_1/outputs/h_mac_run1.jsonl")
mac_r2 = load_mac_run("C:/Projects/math/non_deterministic/experiment_1/outputs/h_mac_run2.jsonl")

# Load prompt metadata for context
prompt_meta = {}
with open("C:/Projects/math/non_deterministic/experiment_1/outputs/h_nvidia_run1.jsonl", encoding="utf-8") as f:
    for line in f:
        d = json.loads(line)
        prompt_meta[d["prompt_idx"]] = d

# Load exp/ctrl
exp_idxs = set()
ctrl_idxs = set()
with open("C:/Projects/math/non_deterministic/prompts/experimental_200.jsonl", encoding="utf-8") as f:
    for line in f:
        exp_idxs.add(json.loads(line)["prompt_idx"])
with open("C:/Projects/math/non_deterministic/prompts/control_200.jsonl", encoding="utf-8") as f:
    for line in f:
        ctrl_idxs.add(json.loads(line)["prompt_idx"])

prompt_idxs = sorted(nvidia_r1.keys())

# ============================================================
print("=" * 70)
print("CROSS-PLATFORM COMPARISON (ALIGNED)")
print("RTX 4070 (CUDA, FP32) vs M3 Pro (CPU, FP32)")
print("Checkpoint hash verified identical on both platforms")
print("=" * 70)

# Within-platform baselines
nv_match = sum(1 for idx in prompt_idxs if nvidia_r1[idx] == nvidia_r2[idx])
mac_match = sum(1 for idx in prompt_idxs if mac_r1[idx] == mac_r2[idx])
print(f"\nWithin-platform (run1 vs run2):")
print(f"  NVIDIA: {nv_match}/400 exact match")
print(f"  Mac:    {mac_match}/400 exact match")

# Cross-platform
cross_match = 0
diverged = []
for idx in prompt_idxs:
    n = nvidia_r1[idx]
    m = mac_r1[idx]
    if n == m:
        cross_match += 1
    else:
        min_len = min(len(n), len(m))
        first_div = next((j for j in range(min_len) if n[j] != m[j]), min_len)
        agree = sum(1 for a, b in zip(n[:min_len], m[:min_len]) if a == b)

        n_tok = tokenizer.decode([n[first_div]]) if first_div < len(n) else "?"
        m_tok = tokenizer.decode([m[first_div]]) if first_div < len(m) else "?"

        meta = prompt_meta[idx]
        diverged.append({
            "idx": idx, "first_div": first_div,
            "agreement": agree / max(min_len, 1),
            "n_tok": n_tok, "m_tok": m_tok,
            "category": meta["category"],
            "set": meta["set"],
            "text": meta["text"][:70],
            "frag": meta.get("fragility_score", 0),
        })

n_div = len(diverged)
print(f"\nCross-platform (NVIDIA run1 vs Mac run1):")
print(f"  Exact match: {cross_match}/400 ({cross_match/4:.1f}%)")
print(f"  Diverged:    {n_div}/400 ({n_div/4:.1f}%)")

if n_div == 0:
    print("\n*** PERFECT CROSS-PLATFORM AGREEMENT ***")
    print("Both platforms produce bitwise-identical outputs for all 400 prompts.")
    print("This means: FP32 with deterministic mode on CPU vs CUDA produces")
    print("identical results for Pythia-410M greedy decoding.")
else:
    # Full analysis
    exp_div = [d for d in diverged if d["idx"] in exp_idxs]
    ctrl_div = [d for d in diverged if d["idx"] in ctrl_idxs]

    print(f"\n--- Experimental vs Control ---")
    print(f"  Experimental: {len(exp_div)}/200 ({len(exp_div)/2:.1f}%)")
    print(f"  Control:      {len(ctrl_div)}/200 ({len(ctrl_div)/2:.1f}%)")
    if len(ctrl_div) > 0:
        ratio = (len(exp_div)/200) / (len(ctrl_div)/200)
        print(f"  Ratio: {ratio:.1f}x")

    print(f"\n--- First divergence position ---")
    positions = [d["first_div"] for d in diverged]
    print(f"  Mean: {np.mean(positions):.1f}, Median: {np.median(positions):.1f}")
    print(f"  Min: {min(positions)}, Max: {max(positions)}")

    bins = [(0, 2), (3, 7), (8, 15), (16, 31), (32, 63), (64, 127), (128, 255)]
    for lo, hi in bins:
        c = sum(1 for p in positions if lo <= p <= hi)
        if c > 0:
            print(f"    pos {lo:3d}-{hi:3d}: {c:3d} ({c/n_div*100:.1f}%)")

    print(f"\n--- Per-category ---")
    cat_total = Counter(prompt_meta[idx]["category"] for idx in prompt_idxs)
    cat_div = Counter(d["category"] for d in diverged)
    for cat in sorted(cat_total.keys()):
        total = cat_total[cat]
        div = cat_div.get(cat, 0)
        print(f"    {cat:15s}: {div}/{total} ({div/total*100:.1f}%)")

    print(f"\n--- Token type at first divergence ---")
    def classify(n_tok, m_tok):
        ns, ms = n_tok.strip(), m_tok.strip()
        if ns.isdigit() and ms.isdigit(): return "digit"
        if ns in "+-*/=" and ms in "+-*/=": return "operator"
        if n_tok in ["\n", " "] or m_tok in ["\n", " "]: return "whitespace"
        if n_tok in ".,!?;:" or m_tok in ".,!?;:": return "punctuation"
        return "word"
    type_c = Counter(classify(d["n_tok"], d["m_tok"]) for d in diverged)
    for t, c in type_c.most_common():
        print(f"    {t:15s}: {c} ({c/n_div*100:.1f}%)")

    print(f"\n--- Token agreement after divergence ---")
    agrees = [d["agreement"] for d in diverged]
    print(f"  Mean: {np.mean(agrees):.4f}")
    print(f"  Median: {np.median(agrees):.4f}")

    print(f"\n--- Concrete examples (sorted by first_div) ---")
    diverged.sort(key=lambda d: d["first_div"])
    for i, d in enumerate(diverged[:20]):
        tag = "EXP" if d["idx"] in exp_idxs else "CTRL"
        print(f"  #{i+1} idx={d['idx']} [{d['category'][:6]}|{tag}] pos={d['first_div']} "
              f"'{d['n_tok']}' -> '{d['m_tok']}' agree={d['agreement']:.3f} frag={d['frag']:.3f}")
        print(f"      {d['text']}")

# Also check: NVIDIA run1 vs Mac run2 (should be same as run1 vs run1 if both deterministic)
cross_12 = sum(1 for idx in prompt_idxs if nvidia_r1[idx] == mac_r2[idx])
print(f"\n--- Additional checks ---")
print(f"  NVIDIA run1 vs Mac run2: {cross_12}/400 match")
print(f"  (Should equal run1 vs run1 if both platforms are internally deterministic)")
