"""
Iteration 13: Predicting fragility from the prompt alone.

Can we predict which prompts will be hardware-fragile without running inference?

Features derived only from prompt TEXT (no model required):
- Category
- Character length, word count
- Contains digits, operators, question marks
- Number of numeric values in prompt
- First word / question word
- Presence of specific keywords (arithmetic, story, describe, etc.)

Target: overall fragility score (proportion of 256 positions with margin < 1.0)
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

prompts = {}
with open(PROMPTS_PATH, encoding="utf-8") as f:
    for line in f:
        p = json.loads(line)
        prompts[p["idx"]] = p

# ============================================================
# Load fragility scores
# ============================================================
print("Loading fragility data...", flush=True)
data = []
for i in range(1999):
    path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
    if not os.path.exists(path):
        continue
    d = torch.load(path, map_location="cpu", weights_only=True)
    margins = (d["top5_values"][:, 0] - d["top5_values"][:, 1]).numpy()
    frag = float((margins < 1.0).mean())
    text = prompts[i]["text"]
    cat = prompts[i]["category"]
    data.append({"idx": i, "text": text, "category": cat, "frag": frag})

N = len(data)
print(f"Loaded {N} prompts", flush=True)

# ============================================================
# Feature extraction (from prompt text only)
# ============================================================
def extract_features(text, category):
    """Extract predictive features from prompt text only."""
    features = {}

    # Basic text stats
    features["char_len"] = len(text)
    features["word_count"] = len(text.split())

    # Category one-hot
    for cat in ["factual", "open_ended", "reasoning", "continuation", "conversational"]:
        features[f"cat_{cat}"] = 1.0 if category == cat else 0.0

    # Numeric content
    numbers = re.findall(r'\d+', text)
    features["num_count"] = len(numbers)
    features["has_digits"] = 1.0 if numbers else 0.0
    features["max_number_digits"] = max((len(n) for n in numbers), default=0)
    features["sum_number_magnitudes"] = sum(len(n) for n in numbers)

    # Operators
    features["has_plus"] = 1.0 if "+" in text else 0.0
    features["has_minus"] = 1.0 if " - " in text or text.startswith("What is ") and "-" in text else 0.0
    features["has_times"] = 1.0 if "times" in text.lower() or "*" in text else 0.0
    features["has_divided"] = 1.0 if "divided" in text.lower() or "/" in text else 0.0
    features["has_percent"] = 1.0 if "%" in text or "percent" in text.lower() else 0.0

    # Question structure
    features["has_question_mark"] = 1.0 if "?" in text else 0.0
    first_word = text.split()[0].lower() if text.split() else ""
    features["starts_what"] = 1.0 if first_word == "what" else 0.0
    features["starts_who"] = 1.0 if first_word == "who" else 0.0
    features["starts_how"] = 1.0 if first_word == "how" else 0.0
    features["starts_can"] = 1.0 if first_word == "can" else 0.0
    features["starts_write"] = 1.0 if first_word == "write" else 0.0
    features["starts_describe"] = 1.0 if first_word == "describe" else 0.0
    features["starts_imagine"] = 1.0 if first_word == "imagine" else 0.0
    features["starts_the"] = 1.0 if first_word == "the" or first_word == '"the' else 0.0

    # Content keywords
    text_lower = text.lower()
    features["kw_sequence"] = 1.0 if "sequence" in text_lower or "comes next" in text_lower else 0.0
    features["kw_recipe"] = 1.0 if "recipe" in text_lower or "flour" in text_lower or "cups" in text_lower else 0.0
    features["kw_workers"] = 1.0 if "workers" in text_lower or "days" in text_lower else 0.0
    features["kw_travel"] = 1.0 if "travel" in text_lower or "miles" in text_lower or "speed" in text_lower else 0.0
    features["kw_buy"] = 1.0 if "buy" in text_lower or "cost" in text_lower or "$" in text else 0.0
    features["kw_capital"] = 1.0 if "capital" in text_lower else 0.0
    features["kw_wrote"] = 1.0 if "wrote" in text_lower else 0.0
    features["kw_symbol"] = 1.0 if "symbol" in text_lower else 0.0

    # Arithmetic detection (strong predictor from iter 12)
    features["is_arithmetic"] = 1.0 if (features["has_plus"] or features["has_minus"] or
                                         features["has_times"] or features["has_divided"] or
                                         features["has_percent"]) and features["cat_reasoning"] else 0.0

    return features

# Extract all features
print("Extracting features...", flush=True)
for d in data:
    d["features"] = extract_features(d["text"], d["category"])

feature_names = sorted(data[0]["features"].keys())
X = np.array([[d["features"][f] for f in feature_names] for d in data])
y = np.array([d["frag"] for d in data])

print(f"Features: {len(feature_names)}")
print(f"Feature names: {feature_names}")

# ============================================================
# PART 1: Individual feature correlations
# ============================================================
print("\n" + "=" * 70)
print("PART 1: Individual feature correlations with fragility")
print("=" * 70)

correlations = {}
for i, fname in enumerate(feature_names):
    r = np.corrcoef(X[:, i], y)[0, 1]
    correlations[fname] = r

print(f"\n{'Feature':35s} {'r':>8s} {'|r|':>8s}")
print("-" * 55)
for fname in sorted(correlations.keys(), key=lambda k: -abs(correlations[k])):
    r = correlations[fname]
    bar = "#" * int(abs(r) * 40)
    print(f"{fname:35s} {r:+8.4f} {abs(r):8.4f}  {bar}")

# ============================================================
# PART 2: Simple rule-based predictor
# ============================================================
print("\n" + "=" * 70)
print("PART 2: Simple rule-based fragility prediction")
print("=" * 70)

# Rule: if is_arithmetic, predict high fragility; if continuation, medium; else low
# Based on established findings

def predict_fragility_simple(features):
    """Simple rule-based predictor."""
    if features["is_arithmetic"]:
        return 0.28  # mean for arithmetic
    elif features["cat_reasoning"] and features["kw_sequence"]:
        return 0.054
    elif features["cat_reasoning"] and (features["kw_travel"] or features["kw_workers"] or
                                         features["kw_recipe"] or features["kw_buy"]):
        return 0.05
    elif features["cat_continuation"]:
        return 0.136
    elif features["cat_reasoning"]:
        return 0.10  # other reasoning
    elif features["cat_factual"]:
        return 0.093
    elif features["cat_conversational"]:
        return 0.088
    else:  # open_ended
        return 0.083

y_pred_simple = np.array([predict_fragility_simple(d["features"]) for d in data])
residuals = y - y_pred_simple
mae = np.abs(residuals).mean()
rmse = np.sqrt((residuals ** 2).mean())
r2 = 1 - (residuals ** 2).sum() / ((y - y.mean()) ** 2).sum()
r_pred = np.corrcoef(y_pred_simple, y)[0, 1]

print(f"Simple rule-based predictor:")
print(f"  R-squared: {r2:.4f}")
print(f"  r(pred, actual): {r_pred:.4f}")
print(f"  MAE: {mae:.4f}")
print(f"  RMSE: {rmse:.4f}")

# ============================================================
# PART 3: Linear regression with all features
# ============================================================
print("\n" + "=" * 70)
print("PART 3: Linear regression (all features)")
print("=" * 70)

# Add intercept
X_aug = np.column_stack([np.ones(N), X])
# OLS solution
try:
    beta = np.linalg.lstsq(X_aug, y, rcond=None)[0]
    y_pred_lr = X_aug @ beta
    residuals_lr = y - y_pred_lr
    r2_lr = 1 - (residuals_lr ** 2).sum() / ((y - y.mean()) ** 2).sum()
    mae_lr = np.abs(residuals_lr).mean()
    r_lr = np.corrcoef(y_pred_lr, y)[0, 1]
    print(f"Linear regression:")
    print(f"  R-squared: {r2_lr:.4f}")
    print(f"  r(pred, actual): {r_lr:.4f}")
    print(f"  MAE: {mae_lr:.4f}")

    # Top coefficients
    print(f"\n  Top 10 coefficients:")
    coef_names = ["intercept"] + feature_names
    coef_ranked = sorted(range(len(coef_names)), key=lambda i: -abs(beta[i]))
    for rank, idx in enumerate(coef_ranked[:10]):
        print(f"    {coef_names[idx]:35s}: {beta[idx]:+.6f}")
except Exception as e:
    print(f"  Error: {e}")

# ============================================================
# PART 4: Can we predict which QUARTILE a prompt falls in?
# ============================================================
print("\n" + "=" * 70)
print("PART 4: Quartile prediction accuracy")
print("=" * 70)

q25, q50, q75 = np.percentile(y, [25, 50, 75])
actual_quartile = np.digitize(y, [q25, q50, q75])

# Using linear regression predictions
pred_quartile_lr = np.digitize(y_pred_lr, [q25, q50, q75])
exact_match = (pred_quartile_lr == actual_quartile).mean()
off_by_one = (np.abs(pred_quartile_lr - actual_quartile) <= 1).mean()
print(f"Linear regression quartile prediction:")
print(f"  Exact quartile match: {exact_match:.1%}")
print(f"  Within 1 quartile: {off_by_one:.1%}")

# Using simple rules
pred_quartile_simple = np.digitize(y_pred_simple, [q25, q50, q75])
exact_simple = (pred_quartile_simple == actual_quartile).mean()
off_simple = (np.abs(pred_quartile_simple - actual_quartile) <= 1).mean()
print(f"\nSimple rule-based quartile prediction:")
print(f"  Exact quartile match: {exact_simple:.1%}")
print(f"  Within 1 quartile: {off_simple:.1%}")

# ============================================================
# PART 5: Can we predict experimental vs control membership?
# ============================================================
print("\n" + "=" * 70)
print("PART 5: Can prompt text predict experimental vs control set?")
print("=" * 70)

from pathlib import Path
exp_idxs = set()
ctrl_idxs = set()
with open("C:/Projects/math/non_deterministic/prompts/experimental_200.jsonl", encoding="utf-8") as f:
    for line in f:
        exp_idxs.add(json.loads(line)["prompt_idx"])
with open("C:/Projects/math/non_deterministic/prompts/control_200.jsonl", encoding="utf-8") as f:
    for line in f:
        ctrl_idxs.add(json.loads(line)["prompt_idx"])

# For each prompt, check if linear regression would have selected it
exp_ctrl_data = [(d, d["idx"] in exp_idxs, d["idx"] in ctrl_idxs) for d in data
                 if d["idx"] in exp_idxs or d["idx"] in ctrl_idxs]

# Get predicted fragility for exp/ctrl prompts
exp_preds = [y_pred_lr[i] for i, d in enumerate(data) if d["idx"] in exp_idxs]
ctrl_preds = [y_pred_lr[i] for i, d in enumerate(data) if d["idx"] in ctrl_idxs]

print(f"Predicted fragility — Experimental: mean={np.mean(exp_preds):.4f}, Control: mean={np.mean(ctrl_preds):.4f}")
print(f"Actual fragility — Experimental: mean={np.mean([d['frag'] for d in data if d['idx'] in exp_idxs]):.4f}, "
      f"Control: mean={np.mean([d['frag'] for d in data if d['idx'] in ctrl_idxs]):.4f}")

# Threshold: predict experimental if predicted fragility > some threshold
threshold = np.median(y_pred_lr)
tp = sum(1 for p in exp_preds if p > threshold)
tn = sum(1 for p in ctrl_preds if p <= threshold)
accuracy = (tp + tn) / 400
print(f"\nMedian-split classification of exp vs ctrl:")
print(f"  TP (exp correctly predicted): {tp}/200")
print(f"  TN (ctrl correctly predicted): {tn}/200")
print(f"  Accuracy: {accuracy:.1%}")

# ============================================================
# PART 6: The "is_arithmetic" feature alone
# ============================================================
print("\n" + "=" * 70)
print("PART 6: How much does 'is_arithmetic' alone explain?")
print("=" * 70)

is_arith = np.array([d["features"]["is_arithmetic"] for d in data])
arith_frags = y[is_arith == 1]
non_arith_frags = y[is_arith == 0]
print(f"Arithmetic prompts (N={int(is_arith.sum())}): mean frag = {arith_frags.mean():.4f}")
print(f"Non-arithmetic (N={int((1-is_arith).sum())}): mean frag = {non_arith_frags.mean():.4f}")
print(f"Ratio: {arith_frags.mean() / non_arith_frags.mean():.1f}x")

# Just category alone
print(f"\nCategory alone as predictor:")
for cat in sorted(set(d["category"] for d in data)):
    cat_frags = [d["frag"] for d in data if d["category"] == cat]
    print(f"  {cat:15s}: mean={np.mean(cat_frags):.4f} std={np.std(cat_frags):.4f} N={len(cat_frags)}")

# R-squared from category alone
cat_means = {cat: np.mean([d["frag"] for d in data if d["category"] == cat])
             for cat in set(d["category"] for d in data)}
y_pred_cat = np.array([cat_means[d["category"]] for d in data])
r2_cat = 1 - ((y - y_pred_cat) ** 2).sum() / ((y - y.mean()) ** 2).sum()
print(f"\nR-squared from category alone: {r2_cat:.4f}")
print(f"R-squared from is_arithmetic alone: {1 - ((y - is_arith * arith_frags.mean() - (1-is_arith) * non_arith_frags.mean()) ** 2).sum() / ((y - y.mean()) ** 2).sum():.4f}")
print(f"R-squared from linear regression (all features): {r2_lr:.4f}")
