"""
Iteration 5: Quality vs Stability Tradeoff

Key question: Is there a fragility "sweet spot" where the model produces
both coherent text AND hardware-stable output?

We know:
- Low fragility (unflippable) → repetitive, low-quality text
- High fragility → uncertain arithmetic, digit soup
- What about the middle?

Approach: decode generated text for all 1999 prompts from logit traces,
compute text quality heuristics, and bin by fragility.
"""

import torch
import numpy as np
import json
import os
from collections import Counter, defaultdict
from pathlib import Path

LOGIT_DIR = "C:/Projects/math/non_deterministic/experiment_0/logit_traces"
PROMPTS_PATH = "C:/Projects/math/non_deterministic/prompts/candidates_2000.jsonl"

from transformers import AutoTokenizer
tokenizer = AutoTokenizer.from_pretrained("EleutherAI/pythia-410m-deduped")

# Load prompt metadata
prompts = {}
with open(PROMPTS_PATH, encoding="utf-8") as f:
    for line in f:
        p = json.loads(line)
        prompts[p["idx"]] = p

# ============================================================
# Load all prompts: margins + generated text
# ============================================================
print("Loading and decoding all 1999 prompts...", flush=True)

all_data = []
for i in range(1999):
    path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
    if not os.path.exists(path):
        continue
    data = torch.load(path, map_location="cpu", weights_only=True)
    top5_vals = data["top5_values"]
    margins = (top5_vals[:, 0] - top5_vals[:, 1]).numpy()
    gen_ids = data["generated_ids"].tolist()
    input_ids = data["input_ids"].tolist()

    # Decode generated text
    gen_text = tokenizer.decode(gen_ids[len(input_ids):], skip_special_tokens=True)

    all_data.append({
        "idx": i,
        "category": prompts[i]["category"],
        "text": prompts[i]["text"],
        "gen_text": gen_text,
        "fragility": float((margins < 1.0).mean()),
        "min_margin": float(margins.min()),
        "mean_margin": float(margins.mean()),
    })

print(f"Loaded {len(all_data)} prompts", flush=True)

# ============================================================
# Text quality metrics
# ============================================================
def compute_quality_metrics(text):
    """Compute heuristic text quality metrics."""
    if not text or len(text.strip()) < 10:
        return {"empty": True}

    words = text.split()
    n_words = len(words)
    if n_words == 0:
        return {"empty": True}

    # 1. Vocabulary diversity (type-token ratio)
    unique_words = len(set(w.lower() for w in words))
    ttr = unique_words / n_words

    # 2. Repetition: fraction of 4-grams that appear more than once
    if n_words >= 4:
        ngrams = [" ".join(words[j:j+4]) for j in range(n_words - 3)]
        ngram_counts = Counter(ngrams)
        repeated_ngrams = sum(1 for c in ngram_counts.values() if c > 1)
        rep_rate = repeated_ngrams / len(ngram_counts) if ngram_counts else 0
    else:
        rep_rate = 0

    # 3. Sentence count (crude)
    sentences = [s.strip() for s in text.replace("!", ".").replace("?", ".").split(".") if s.strip()]
    n_sentences = len(sentences)

    # 4. Mean sentence length
    mean_sent_len = n_words / max(n_sentences, 1)

    # 5. Character-level repetition (longest repeated substring >= 20 chars)
    has_long_repeat = False
    for start in range(0, min(len(text) - 20, 500), 10):
        chunk = text[start:start+20]
        if text.count(chunk) >= 3:
            has_long_repeat = True
            break

    # 6. Newline density (paragraph structure)
    newline_count = text.count("\n")
    newline_density = newline_count / max(len(text), 1) * 100

    # 7. Digit density
    digit_count = sum(1 for c in text if c.isdigit())
    digit_density = digit_count / max(len(text), 1)

    # 8. "Composite quality" — heuristic combining diversity and anti-repetition
    # Higher = better text. TTR rewards diversity, (1-rep_rate) penalizes repetition.
    composite = ttr * (1 - rep_rate)

    return {
        "empty": False,
        "n_words": n_words,
        "n_unique_words": unique_words,
        "ttr": ttr,  # type-token ratio (higher = more diverse)
        "rep_rate": rep_rate,  # 4-gram repetition rate (lower = less repetitive)
        "n_sentences": n_sentences,
        "mean_sent_len": mean_sent_len,
        "has_long_repeat": has_long_repeat,
        "newline_density": newline_density,
        "digit_density": digit_density,
        "composite": composite,
    }

print("Computing quality metrics...", flush=True)
for d in all_data:
    d.update(compute_quality_metrics(d["gen_text"]))

# Filter out empty
valid = [d for d in all_data if not d.get("empty", True)]
print(f"Valid (non-empty) prompts: {len(valid)}")

# ============================================================
# PART 1: Bin by fragility and compute average quality
# ============================================================
print("\n" + "=" * 70)
print("PART 1: Quality metrics by fragility bin")
print("=" * 70)

bins = [
    (0.00, 0.03, "Very low (0-3%)"),
    (0.03, 0.06, "Low (3-6%)"),
    (0.06, 0.10, "Medium-low (6-10%)"),
    (0.10, 0.15, "Medium (10-15%)"),
    (0.15, 0.25, "Medium-high (15-25%)"),
    (0.25, 0.40, "High (25-40%)"),
    (0.40, 1.00, "Very high (40%+)"),
]

print(f"\n{'Bin':25s} {'N':>5s} {'TTR':>6s} {'RepRate':>8s} {'Composite':>10s} "
      f"{'LongRep%':>8s} {'Digits%':>8s} {'MinMarg':>8s}")
print("-" * 85)

bin_data = {}
for lo, hi, label in bins:
    group = [d for d in valid if lo <= d["fragility"] < hi]
    if not group:
        continue
    ttrs = [d["ttr"] for d in group]
    reps = [d["rep_rate"] for d in group]
    comps = [d["composite"] for d in group]
    long_reps = [d["has_long_repeat"] for d in group]
    digits = [d["digit_density"] for d in group]
    min_margins = [d["min_margin"] for d in group]

    bin_data[label] = {
        "n": len(group),
        "ttr": np.mean(ttrs),
        "rep": np.mean(reps),
        "comp": np.mean(comps),
        "long_rep_pct": np.mean(long_reps) * 100,
        "digit_pct": np.mean(digits) * 100,
        "min_margin": np.mean(min_margins),
    }

    print(f"{label:25s} {len(group):5d} {np.mean(ttrs):6.3f} {np.mean(reps):8.3f} "
          f"{np.mean(comps):10.4f} {np.mean(long_reps)*100:7.1f}% "
          f"{np.mean(digits)*100:7.2f}% {np.mean(min_margins):8.4f}")

# ============================================================
# PART 2: Per-category quality vs fragility
# ============================================================
print("\n" + "=" * 70)
print("PART 2: Quality by category (composite score)")
print("=" * 70)

categories = sorted(set(d["category"] for d in valid))
for cat in categories:
    cat_data = [d for d in valid if d["category"] == cat]
    # Split into tertiles by fragility within category
    cat_data.sort(key=lambda x: x["fragility"])
    n = len(cat_data)
    t1 = cat_data[:n//3]
    t2 = cat_data[n//3:2*n//3]
    t3 = cat_data[2*n//3:]

    print(f"\n  {cat} ({n} prompts):")
    for label, group in [("Bottom third (stable)", t1),
                          ("Middle third", t2),
                          ("Top third (fragile)", t3)]:
        comp = np.mean([d["composite"] for d in group])
        ttr = np.mean([d["ttr"] for d in group])
        rep = np.mean([d["rep_rate"] for d in group])
        frag = np.mean([d["fragility"] for d in group])
        print(f"    {label:25s}: composite={comp:.4f} ttr={ttr:.3f} "
              f"rep={rep:.3f} frag={frag:.3f}")

# ============================================================
# PART 3: Find the sweet spot — highest composite score
# ============================================================
print("\n" + "=" * 70)
print("PART 3: Sweet spot analysis — sliding window over fragility")
print("=" * 70)

# Sort all by fragility
valid.sort(key=lambda x: x["fragility"])

# Sliding window of 100 prompts
window = 100
best_comp = 0
best_center = 0
results = []
for start in range(0, len(valid) - window, 10):
    group = valid[start:start+window]
    mean_frag = np.mean([d["fragility"] for d in group])
    mean_comp = np.mean([d["composite"] for d in group])
    mean_ttr = np.mean([d["ttr"] for d in group])
    mean_rep = np.mean([d["rep_rate"] for d in group])
    results.append((mean_frag, mean_comp, mean_ttr, mean_rep))
    if mean_comp > best_comp:
        best_comp = mean_comp
        best_center = mean_frag

print(f"Fragility vs Composite Quality (window={window}):")
print(f"{'Fragility':>10s} {'Composite':>10s} {'TTR':>6s} {'RepRate':>8s}")
print("-" * 40)
for frag, comp, ttr, rep in results[::3]:  # every 3rd for readability
    bar = "#" * int(comp * 100)
    print(f"{frag:10.3f} {comp:10.4f} {ttr:6.3f} {rep:8.3f}  {bar}")

print(f"\nBest composite quality at fragility ~ {best_center:.3f} (composite={best_comp:.4f})")

# ============================================================
# PART 4: Concrete examples from the sweet spot
# ============================================================
print("\n" + "=" * 70)
print("PART 4: Concrete examples from the sweet spot")
print("=" * 70)

# Find prompts near the sweet spot (within 0.02 of best_center)
sweet_spot = [d for d in valid if abs(d["fragility"] - best_center) < 0.02]
sweet_spot.sort(key=lambda x: -x["composite"])

print(f"Prompts near sweet spot (fragility ~ {best_center:.3f}, ±0.02): {len(sweet_spot)}")
for i, d in enumerate(sweet_spot[:8]):
    gen_preview = d["gen_text"][:200].replace("\n", "\\n")
    print(f"\n  #{i+1} [{d['category']}] frag={d['fragility']:.3f} comp={d['composite']:.4f} "
          f"ttr={d['ttr']:.3f} rep={d['rep_rate']:.3f}")
    print(f"    Prompt: {d['text'][:70]}")
    print(f"    Gen: {gen_preview}")

# ============================================================
# PART 5: Examples from extremes for comparison
# ============================================================
print("\n" + "=" * 70)
print("PART 5: Extreme examples for comparison")
print("=" * 70)

# Most stable, high quality (low fragility, high composite)
low_frag_sorted = sorted([d for d in valid if d["fragility"] < 0.05],
                          key=lambda x: -x["composite"])
print("\n--- Most stable + highest quality (frag < 5%) ---")
for i, d in enumerate(low_frag_sorted[:5]):
    gen_preview = d["gen_text"][:150].replace("\n", "\\n")
    print(f"  #{i+1} [{d['category']}] frag={d['fragility']:.3f} comp={d['composite']:.4f}")
    print(f"    Prompt: {d['text'][:70]}")
    print(f"    Gen: {gen_preview}")

# Most fragile, any quality
high_frag_sorted = sorted([d for d in valid if d["fragility"] > 0.4],
                           key=lambda x: -x["composite"])
print("\n--- Most fragile + best quality among them (frag > 40%) ---")
for i, d in enumerate(high_frag_sorted[:5]):
    gen_preview = d["gen_text"][:150].replace("\n", "\\n")
    print(f"  #{i+1} [{d['category']}] frag={d['fragility']:.3f} comp={d['composite']:.4f}")
    print(f"    Prompt: {d['text'][:70]}")
    print(f"    Gen: {gen_preview}")

# ============================================================
# PART 6: Correlation between fragility and quality metrics
# ============================================================
print("\n" + "=" * 70)
print("PART 6: Correlations")
print("=" * 70)

frags = np.array([d["fragility"] for d in valid])
ttrs = np.array([d["ttr"] for d in valid])
reps = np.array([d["rep_rate"] for d in valid])
comps = np.array([d["composite"] for d in valid])
digits = np.array([d["digit_density"] for d in valid])

for name, arr in [("TTR", ttrs), ("RepRate", reps), ("Composite", comps), ("DigitDensity", digits)]:
    r = np.corrcoef(frags, arr)[0, 1]
    print(f"  Fragility vs {name:15s}: r = {r:+.4f}")

# Non-reasoning only
non_reason = [d for d in valid if d["category"] != "reasoning"]
frags_nr = np.array([d["fragility"] for d in non_reason])
comps_nr = np.array([d["composite"] for d in non_reason])
ttrs_nr = np.array([d["ttr"] for d in non_reason])
reps_nr = np.array([d["rep_rate"] for d in non_reason])
print(f"\n  Non-reasoning only ({len(non_reason)} prompts):")
for name, arr in [("TTR", ttrs_nr), ("RepRate", reps_nr), ("Composite", comps_nr)]:
    r = np.corrcoef(frags_nr, arr)[0, 1]
    print(f"    Fragility vs {name:15s}: r = {r:+.4f}")
