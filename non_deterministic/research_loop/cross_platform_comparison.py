"""
Cross-Platform Comparison: H_nvidia (RTX 4070, CUDA) vs H_mac (M3 Pro, CPU)

This is the moment of truth — comparing actual outputs from two different
hardware platforms running the same model with the same weights.
"""

import json
import numpy as np
from collections import Counter, defaultdict
from transformers import AutoTokenizer

tokenizer = AutoTokenizer.from_pretrained("EleutherAI/pythia-410m-deduped")

# Load both platforms
def load_run(path):
    data = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            data.append(json.loads(line))
    return {d["prompt_idx"]: d for d in data}

nvidia_r1 = load_run("C:/Projects/math/non_deterministic/experiment_1/outputs/h_nvidia_run1.jsonl")
nvidia_r2 = load_run("C:/Projects/math/non_deterministic/experiment_1/outputs/h_nvidia_run2.jsonl")
mac_r1 = load_run("C:/Projects/math/non_deterministic/experiment_1/outputs/h_mac_run1.jsonl")
mac_r2 = load_run("C:/Projects/math/non_deterministic/experiment_1/outputs/h_mac_run2.jsonl")

# Load experimental/control membership
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
# PART 1: Overall divergence statistics
# ============================================================
print("=" * 70)
print("CROSS-PLATFORM COMPARISON: H_nvidia vs H_mac")
print("RTX 4070 (CUDA, FP32) vs M3 Pro (CPU, FP32)")
print("=" * 70)

# First: within-platform baselines
print("\n--- Within-platform baselines ---")
nvidia_within = sum(1 for idx in prompt_idxs
                    if nvidia_r1[idx]["generated_ids"] == nvidia_r2[idx]["generated_ids"])
mac_within = sum(1 for idx in prompt_idxs
                 if mac_r1[idx]["generated_ids"] == mac_r2[idx]["generated_ids"])
print(f"NVIDIA run1 vs run2: {nvidia_within}/400 exact match ({nvidia_within/4:.1f}%)")
print(f"Mac run1 vs run2:    {mac_within}/400 exact match ({mac_within/4:.1f}%)")

# Cross-platform
print("\n--- Cross-platform: NVIDIA run1 vs Mac run1 ---")
cross_exact = 0
cross_diverged = []
for idx in prompt_idxs:
    n_gen = nvidia_r1[idx]["generated_ids"]
    m_gen = mac_r1[idx]["generated_ids"]
    if n_gen == m_gen:
        cross_exact += 1
    else:
        # Find first divergence point
        min_len = min(len(n_gen), len(m_gen))
        first_div = None
        for j in range(min_len):
            if n_gen[j] != m_gen[j]:
                first_div = j
                break
        if first_div is None:
            first_div = min_len  # length difference

        # Token agreement
        agree = sum(1 for a, b in zip(n_gen[:min_len], m_gen[:min_len]) if a == b)
        agreement_rate = agree / max(min_len, 1)

        cross_diverged.append({
            "idx": idx,
            "first_div": first_div,
            "agreement_rate": agreement_rate,
            "n_tok": tokenizer.decode([n_gen[first_div]]) if first_div < len(n_gen) else "?",
            "m_tok": tokenizer.decode([m_gen[first_div]]) if first_div < len(m_gen) else "?",
            "category": nvidia_r1[idx]["category"],
            "set": nvidia_r1[idx]["set"],
            "text": nvidia_r1[idx]["text"][:70],
        })

n_diverged = len(cross_diverged)
print(f"Exact match:  {cross_exact}/400 ({cross_exact/4:.1f}%)")
print(f"Diverged:     {n_diverged}/400 ({n_diverged/4:.1f}%)")

# ============================================================
# PART 2: Experimental vs Control breakdown
# ============================================================
print("\n" + "=" * 70)
print("PART 2: Experimental vs Control")
print("=" * 70)

exp_div = [d for d in cross_diverged if d["idx"] in exp_idxs]
ctrl_div = [d for d in cross_diverged if d["idx"] in ctrl_idxs]
exp_total = sum(1 for idx in prompt_idxs if idx in exp_idxs)
ctrl_total = sum(1 for idx in prompt_idxs if idx in ctrl_idxs)

print(f"Experimental: {len(exp_div)}/{exp_total} diverged ({len(exp_div)/exp_total*100:.1f}%)")
print(f"Control:      {len(ctrl_div)}/{ctrl_total} diverged ({len(ctrl_div)/ctrl_total*100:.1f}%)")
if len(ctrl_div) > 0:
    ratio = len(exp_div)/exp_total / (len(ctrl_div)/ctrl_total)
    print(f"Ratio:        {ratio:.1f}x")
else:
    print(f"Ratio:        inf (no control divergences)")

# ============================================================
# PART 3: First divergence position distribution
# ============================================================
print("\n" + "=" * 70)
print("PART 3: First divergence position")
print("=" * 70)

if cross_diverged:
    div_positions = [d["first_div"] for d in cross_diverged]
    print(f"First divergence position:")
    print(f"  Mean:   {np.mean(div_positions):.1f}")
    print(f"  Median: {np.median(div_positions):.1f}")
    print(f"  Min:    {min(div_positions)}")
    print(f"  Max:    {max(div_positions)}")

    bins = [(0, 2), (3, 7), (8, 15), (16, 31), (32, 63), (64, 127), (128, 255)]
    print(f"\n  Distribution:")
    for lo, hi in bins:
        count = sum(1 for p in div_positions if lo <= p <= hi)
        pct = count / len(div_positions) * 100
        bar = "#" * int(pct)
        print(f"    pos {lo:3d}-{hi:3d}: {count:3d} ({pct:5.1f}%) {bar}")

# ============================================================
# PART 4: Per-category breakdown
# ============================================================
print("\n" + "=" * 70)
print("PART 4: Per-category divergence")
print("=" * 70)

cat_counts = Counter(nvidia_r1[idx]["category"] for idx in prompt_idxs)
cat_div = Counter(d["category"] for d in cross_diverged)

print(f"{'Category':15s} {'Total':>6s} {'Diverged':>9s} {'Rate':>7s}")
print("-" * 40)
for cat in sorted(cat_counts.keys()):
    total = cat_counts[cat]
    div = cat_div.get(cat, 0)
    rate = div / total * 100
    print(f"{cat:15s} {total:6d} {div:9d} {rate:6.1f}%")

# ============================================================
# PART 5: What tokens diverge?
# ============================================================
print("\n" + "=" * 70)
print("PART 5: Token-level analysis of divergences")
print("=" * 70)

if cross_diverged:
    # Classify divergence types
    def classify_div(n_tok, m_tok):
        ns, ms = n_tok.strip(), m_tok.strip()
        if ns.isdigit() and ms.isdigit():
            return "digit-digit"
        elif ns in "+-*/=" and ms in "+-*/=":
            return "operator"
        elif n_tok in ["\n", " "] or m_tok in ["\n", " "]:
            return "whitespace"
        elif ns and ms and ns[0].isupper() and ms[0].isupper():
            return "sentence-start"
        elif n_tok in [".", ",", "!", "?", ";", ":"] or m_tok in [".", ",", "!", "?", ";", ":"]:
            return "punctuation"
        else:
            return "word"

    type_counter = Counter()
    for d in cross_diverged:
        type_counter[classify_div(d["n_tok"], d["m_tok"])] += 1

    print(f"First-divergence token type:")
    for t, count in type_counter.most_common():
        print(f"  {t:15s}: {count:4d} ({count/len(cross_diverged)*100:.1f}%)")

# ============================================================
# PART 6: Token agreement rate distribution
# ============================================================
print("\n" + "=" * 70)
print("PART 6: Token agreement rates for diverged prompts")
print("=" * 70)

if cross_diverged:
    agree_rates = [d["agreement_rate"] for d in cross_diverged]
    print(f"Token agreement rate (among diverged prompts):")
    print(f"  Mean:   {np.mean(agree_rates):.4f}")
    print(f"  Median: {np.median(agree_rates):.4f}")
    print(f"  Min:    {min(agree_rates):.4f}")
    print(f"  Max:    {max(agree_rates):.4f}")

    # How many tokens typically agree after divergence?
    print(f"\n  After first divergence, do they reconverge?")
    reconverge_count = 0
    for d in cross_diverged:
        idx = d["idx"]
        n_gen = nvidia_r1[idx]["generated_ids"]
        m_gen = mac_r1[idx]["generated_ids"]
        fd = d["first_div"]
        # Check if any tokens match AFTER the first divergence
        post_div_matches = sum(1 for a, b in zip(n_gen[fd+1:], m_gen[fd+1:]) if a == b)
        post_div_total = min(len(n_gen), len(m_gen)) - fd - 1
        if post_div_total > 0 and post_div_matches / post_div_total > 0.5:
            reconverge_count += 1
    print(f"  Prompts with >50% post-divergence agreement: {reconverge_count}/{len(cross_diverged)}")

# ============================================================
# PART 7: Concrete examples
# ============================================================
print("\n" + "=" * 70)
print("PART 7: Concrete divergence examples (first 15)")
print("=" * 70)

if cross_diverged:
    cross_diverged.sort(key=lambda d: d["first_div"])
    for i, d in enumerate(cross_diverged[:15]):
        idx = d["idx"]
        set_tag = "EXP" if idx in exp_idxs else "CTRL"
        n_gen_text = nvidia_r1[idx].get("generated_text", "")[:80]
        m_gen_text = mac_r1[idx].get("generated_text", "")[:80]

        print(f"\n  #{i+1} idx={idx} [{d['category'][:6]}|{set_tag}] first_div=pos {d['first_div']}")
        print(f"    Prompt: {d['text']}")
        print(f"    NVIDIA tok: '{d['n_tok']}' -> Mac tok: '{d['m_tok']}'")
        print(f"    NVIDIA: {n_gen_text}...")
        print(f"    Mac:    {m_gen_text}...")
        print(f"    Agreement: {d['agreement_rate']:.3f}")

# ============================================================
# PART 8: Compare with simulated predictions
# ============================================================
print("\n" + "=" * 70)
print("PART 8: Actual vs Predicted divergence")
print("=" * 70)

print(f"\nActual cross-platform divergence: {n_diverged}/400 = {n_diverged/4:.1f}%")
print(f"Exp: {len(exp_div)}/{exp_total} = {len(exp_div)/exp_total*100:.1f}%")
print(f"Ctrl: {len(ctrl_div)}/{ctrl_total} = {len(ctrl_div)/ctrl_total*100:.1f}%")
if len(ctrl_div) > 0:
    actual_ratio = (len(exp_div)/exp_total) / (len(ctrl_div)/ctrl_total)
else:
    actual_ratio = float('inf')

print(f"\nSimulated predictions (from iteration 3/11):")
print(f"  eps=1e-3: 4.4% diverge, ratio 9.0x  (NVIDIA vs NVIDIA FP32)")
print(f"  eps=5e-3: 19.7% diverge, ratio 4.2x  (NVIDIA vs Apple Silicon FP32)")
print(f"  eps=1e-2: 31.7% diverge, ratio 3.9x  (BF16 cross-platform)")
print(f"\nActual: {n_diverged/4:.1f}% diverge, ratio {actual_ratio:.1f}x")

# Which epsilon best matches?
pred_table = [(1e-4, 0.3), (1e-3, 4.4), (5e-3, 19.7), (1e-2, 31.7), (5e-2, 74.5)]
actual_pct = n_diverged / 4
best_eps = min(pred_table, key=lambda x: abs(x[1] - actual_pct))
print(f"Best matching epsilon: ~{best_eps[0]:.0e} (predicted {best_eps[1]}%, actual {actual_pct:.1f}%)")
