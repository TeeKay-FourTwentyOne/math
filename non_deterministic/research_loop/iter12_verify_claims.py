"""
Iteration 12: Verify two longstanding unverified claims.

Claim 1: r=-0.82 length-fragility within reasoning is confounded by problem type.
Test: within each problem type (addition, multiplication, etc.), does length
still predict fragility? If yes, length matters independently. If no, it's confounded.

Claim 2: 100% repetition rate in unflippable is a detection artifact.
Test: use multiple repetition detection methods (20-char, 10-char, 4-gram,
exact sentence repeat) and compare rates.
"""

import torch
import numpy as np
import json
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

LOGIT_DIR = "C:/Projects/math/non_deterministic/experiment_0/logit_traces"
PROMPTS_PATH = "C:/Projects/math/non_deterministic/prompts/candidates_2000.jsonl"
CTRL_PATH = "C:/Projects/math/non_deterministic/prompts/control_200.jsonl"

from transformers import AutoTokenizer
tokenizer = AutoTokenizer.from_pretrained("EleutherAI/pythia-410m-deduped")

prompts = {}
with open(PROMPTS_PATH, encoding="utf-8") as f:
    for line in f:
        p = json.loads(line)
        prompts[p["idx"]] = p

ctrl_idxs = set()
with open(CTRL_PATH, encoding="utf-8") as f:
    for line in f:
        ctrl_idxs.add(json.loads(line)["prompt_idx"])

# ============================================================
# CLAIM 1: Length-fragility confounded by problem type?
# ============================================================
print("=" * 70)
print("CLAIM 1: Is r=-0.82 length-fragility within reasoning confounded?")
print("=" * 70)

# Load reasoning prompts with fragility and classify by problem type
reasoning_data = []
for i in range(1999):
    if prompts[i]["category"] != "reasoning":
        continue
    path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
    if not os.path.exists(path):
        continue
    data = torch.load(path, map_location="cpu", weights_only=True)
    margins = (data["top5_values"][:, 0] - data["top5_values"][:, 1]).numpy()
    frag = float((margins < 1.0).mean())
    text = prompts[i]["text"]

    # Classify problem type
    text_lower = text.lower()
    if "+" in text or "plus" in text_lower:
        ptype = "addition"
    elif "times" in text_lower or "*" in text:
        ptype = "multiplication"
    elif "-" in text and ("what is" in text_lower):
        ptype = "subtraction"
    elif "divided" in text_lower or "/" in text:
        ptype = "division"
    elif "%" in text or "percent" in text_lower:
        ptype = "percentage"
    elif "sequence" in text_lower or "comes next" in text_lower:
        ptype = "sequence"
    elif "recipe" in text_lower or "cups" in text_lower or "flour" in text_lower:
        ptype = "scaling"
    elif "workers" in text_lower or "days" in text_lower:
        ptype = "workers"
    elif "travel" in text_lower or "miles" in text_lower or "speed" in text_lower or "mph" in text_lower:
        ptype = "distance"
    elif "buy" in text_lower or "$" in text or "cost" in text_lower or "price" in text_lower:
        ptype = "money"
    elif "shirt" in text_lower or "sale" in text_lower or "off" in text_lower:
        ptype = "discount"
    else:
        ptype = "other"

    reasoning_data.append({
        "idx": i, "text": text, "frag": frag,
        "char_len": len(text), "word_count": len(text.split()),
        "ptype": ptype,
    })

# Overall correlation
all_lens = np.array([d["char_len"] for d in reasoning_data])
all_frags = np.array([d["frag"] for d in reasoning_data])
r_overall = np.corrcoef(all_lens, all_frags)[0, 1]
print(f"\nOverall (all reasoning): r(char_len, frag) = {r_overall:.4f} (N={len(reasoning_data)})")

# Per problem type
print(f"\nPer problem type:")
print(f"{'Type':15s} {'N':>5s} {'MeanLen':>8s} {'MeanFrag':>9s} {'r(len,frag)':>12s}")
print("-" * 55)
type_stats = defaultdict(list)
for d in reasoning_data:
    type_stats[d["ptype"]].append(d)

for ptype in sorted(type_stats.keys(), key=lambda k: -len(type_stats[k])):
    group = type_stats[ptype]
    if len(group) < 5:
        continue
    lens = np.array([d["char_len"] for d in group])
    frags = np.array([d["frag"] for d in group])
    r = np.corrcoef(lens, frags)[0, 1] if len(group) > 2 else float("nan")
    print(f"{ptype:15s} {len(group):5d} {lens.mean():8.1f} {frags.mean():9.4f} {r:+12.4f}")

# Key test: after controlling for problem type, does length still matter?
# Residualize fragility against problem type, then correlate with length
ptype_means = {pt: np.mean([d["frag"] for d in group]) for pt, group in type_stats.items()}
residual_frags = np.array([d["frag"] - ptype_means[d["ptype"]] for d in reasoning_data])
r_residual = np.corrcoef(all_lens, residual_frags)[0, 1]
print(f"\nAfter removing problem-type means: r(char_len, residual_frag) = {r_residual:.4f}")
print(f"Original r = {r_overall:.4f}")
print(f"Reduction: {abs(r_overall) - abs(r_residual):.4f} ({(1 - abs(r_residual)/abs(r_overall))*100:.0f}% of correlation explained by type)")

# Even more rigorous: within JUST addition prompts (same type, varying length)
addition = type_stats.get("addition", [])
if len(addition) > 10:
    add_lens = np.array([d["char_len"] for d in addition])
    add_frags = np.array([d["frag"] for d in addition])
    r_add = np.corrcoef(add_lens, add_frags)[0, 1]
    print(f"\nWithin addition only (N={len(addition)}): r = {r_add:.4f}")
    # Show examples
    addition.sort(key=lambda d: d["char_len"])
    print(f"  Shortest: len={addition[0]['char_len']:3d} frag={addition[0]['frag']:.3f} '{addition[0]['text'][:50]}'")
    print(f"  Longest:  len={addition[-1]['char_len']:3d} frag={addition[-1]['frag']:.3f} '{addition[-1]['text'][:50]}'")

mult = type_stats.get("multiplication", [])
if len(mult) > 10:
    mult_lens = np.array([d["char_len"] for d in mult])
    mult_frags = np.array([d["frag"] for d in mult])
    r_mult = np.corrcoef(mult_lens, mult_frags)[0, 1]
    print(f"Within multiplication only (N={len(mult)}): r = {r_mult:.4f}")

distance = type_stats.get("distance", [])
if len(distance) > 10:
    dist_lens = np.array([d["char_len"] for d in distance])
    dist_frags = np.array([d["frag"] for d in distance])
    r_dist = np.corrcoef(dist_lens, dist_frags)[0, 1]
    print(f"Within distance only (N={len(distance)}): r = {r_dist:.4f}")

# ============================================================
# CLAIM 2: 100% repetition rate in unflippable — artifact?
# ============================================================
print("\n\n" + "=" * 70)
print("CLAIM 2: Is 100% repetition in unflippable a detection artifact?")
print("=" * 70)

# Load generated text for control prompts
ctrl_data = []
for i in range(1999):
    if i not in ctrl_idxs:
        continue
    path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
    if not os.path.exists(path):
        continue
    data = torch.load(path, map_location="cpu", weights_only=True)
    margins = (data["top5_values"][:, 0] - data["top5_values"][:, 1]).numpy()
    gen_ids = data["generated_ids"].tolist()
    input_ids = data["input_ids"].tolist()
    gen_text = tokenizer.decode(gen_ids[len(input_ids):], skip_special_tokens=True)
    min_margin = float(margins.min())

    ctrl_data.append({
        "idx": i, "gen_text": gen_text, "min_margin": min_margin,
        "unflippable": min_margin > 0.1,
    })

unflippable = [d for d in ctrl_data if d["unflippable"]]
flippable = [d for d in ctrl_data if not d["unflippable"]]

print(f"\nControl prompts: {len(ctrl_data)} total, {len(unflippable)} unflippable, {len(flippable)} flippable")

# Multiple repetition detection methods
def detect_rep_20char(text):
    """Original method: 20-char substring appearing 3+ times."""
    for start in range(0, min(len(text) - 20, 500), 10):
        chunk = text[start:start+20]
        if text.count(chunk) >= 3:
            return True
    return False

def detect_rep_sentence(text):
    """Exact sentence repeated 2+ times."""
    sentences = [s.strip() for s in re.split(r'[.!?\n]', text) if len(s.strip()) > 10]
    if not sentences:
        return False
    counts = Counter(sentences)
    return any(c >= 2 for c in counts.values())

def detect_rep_4gram(text):
    """4-word n-gram repeated 3+ times (from iter 5)."""
    words = text.split()
    if len(words) < 4:
        return False
    ngrams = [" ".join(words[j:j+4]) for j in range(len(words) - 3)]
    counts = Counter(ngrams)
    return any(c >= 3 for c in counts.values())

def detect_rep_paragraph(text):
    """Full paragraph (split by double newline) repeated."""
    paras = [p.strip() for p in text.split("\n\n") if len(p.strip()) > 20]
    if len(paras) < 2:
        return False
    counts = Counter(paras)
    return any(c >= 2 for c in counts.values())

def detect_rep_strict(text):
    """Strict: 50-char substring appearing 2+ times, non-overlapping."""
    for start in range(0, min(len(text) - 50, 300), 25):
        chunk = text[start:start+50]
        # Count non-overlapping occurrences
        count = 0
        pos = 0
        while True:
            idx = text.find(chunk, pos)
            if idx == -1:
                break
            count += 1
            pos = idx + 50
        if count >= 2:
            return True
    return False

methods = {
    "20-char x3 (original)": detect_rep_20char,
    "Sentence repeat x2": detect_rep_sentence,
    "4-gram x3": detect_rep_4gram,
    "Paragraph repeat": detect_rep_paragraph,
    "50-char x2 strict": detect_rep_strict,
}

print(f"\nRepetition rates by detection method:")
print(f"{'Method':30s} {'Unflippable':>12s} {'Flippable':>12s} {'Difference':>12s}")
print("-" * 70)

for method_name, method_fn in methods.items():
    unflip_rate = np.mean([method_fn(d["gen_text"]) for d in unflippable])
    flip_rate = np.mean([method_fn(d["gen_text"]) for d in flippable])
    diff = unflip_rate - flip_rate
    print(f"{method_name:30s} {unflip_rate:11.1%} {flip_rate:11.1%} {diff:+11.1%}")

# Show the unflippable prompts that DON'T show repetition under strict method
strict_non_rep = [d for d in unflippable if not detect_rep_strict(d["gen_text"])]
print(f"\nUnflippable WITHOUT 50-char strict repetition: {len(strict_non_rep)}/{len(unflippable)}")
for d in strict_non_rep[:5]:
    preview = d["gen_text"][:150].replace("\n", "\\n")
    print(f"  idx={d['idx']:4d} min_margin={d['min_margin']:.3f}")
    print(f"    {preview}")

# Also: how repetitive are ALL 1999 prompts by the strict method?
print(f"\n\nRepetition rate (50-char strict) by fragility quartile:")
all_1999 = []
for i in range(1999):
    path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
    if not os.path.exists(path):
        continue
    data = torch.load(path, map_location="cpu", weights_only=True)
    margins = (data["top5_values"][:, 0] - data["top5_values"][:, 1]).numpy()
    gen_ids = data["generated_ids"].tolist()
    input_ids = data["input_ids"].tolist()
    gen_text = tokenizer.decode(gen_ids[len(input_ids):], skip_special_tokens=True)
    frag = float((margins < 1.0).mean())
    all_1999.append({"frag": frag, "gen_text": gen_text})

all_1999.sort(key=lambda d: d["frag"])
q_size = len(all_1999) // 4
for q in range(4):
    group = all_1999[q * q_size:(q + 1) * q_size]
    rate = np.mean([detect_rep_strict(d["gen_text"]) for d in group])
    frag_range = f"{group[0]['frag']:.3f}-{group[-1]['frag']:.3f}"
    print(f"  Q{q+1} (frag {frag_range}): {rate:.1%} repetitive (strict)")
