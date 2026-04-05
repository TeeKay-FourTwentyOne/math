"""
Verify the 1-vs-3 digit confusion dominance.

The confusion matrix shows 1 vs 3 appearing 603 times as the most common pair.
This is surprising - why would 1 and 3 be confused more than adjacent digits?

Hypotheses:
1. Artifact of BPE tokenization (the tokens '1' and '3' are used in contexts
   where other representations like '10', '13', etc. are also likely)
2. Specific to arithmetic exponents (x^2 vs x^3)
3. Related to the model's internal number representation
"""

import torch
import numpy as np
import json
import os
from collections import Counter, defaultdict

from transformers import AutoTokenizer
tokenizer = AutoTokenizer.from_pretrained("EleutherAI/pythia-410m")

LOGIT_DIR = "C:/Projects/math/non_deterministic/experiment_0/logit_traces"
PROMPTS_FILE = "C:/Projects/math/non_deterministic/prompts/candidates_2000.jsonl"
MARGIN_THRESHOLD = 0.1

prompts = {}
with open(PROMPTS_FILE) as f:
    for line in f:
        d = json.loads(line)
        prompts[d['idx']] = d

def decode(tid):
    return tokenizer.decode([tid])

# Collect all 1-vs-3 and 3-vs-1 instances
one_three_cases = []

for i in range(1999):
    path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
    if not os.path.exists(path):
        continue
    data = torch.load(path, map_location='cpu', weights_only=False)
    top5_vals = data['top5_values']
    top5_idx = data['top5_indices']
    margins = top5_vals[:, 0] - top5_vals[:, 1]
    gen_ids = data['generated_ids']
    input_len = len(data['input_ids'])

    for pos in range(len(margins)):
        if margins[pos] < MARGIN_THRESHOLD:
            t1 = decode(top5_idx[pos, 0].item())
            t2 = decode(top5_idx[pos, 1].item())
            if (t1.strip() == '1' and t2.strip() == '3') or (t1.strip() == '3' and t2.strip() == '1'):
                actual_pos = input_len + pos
                ctx_start = max(0, actual_pos - 8)
                ctx_end = min(len(gen_ids), actual_pos + 5)
                context = tokenizer.decode(gen_ids[ctx_start:ctx_end])
                one_three_cases.append({
                    'prompt_idx': i,
                    'pos': pos,
                    'tok1': t1.strip(),
                    'tok2': t2.strip(),
                    'margin': margins[pos].item(),
                    'category': prompts[i]['category'],
                    'context': context,
                    'prompt': prompts[i]['text'][:80]
                })

print(f"Total 1-vs-3 fragile positions: {len(one_three_cases)}")

# How many unique prompts?
unique_prompts = len(set(c['prompt_idx'] for c in one_three_cases))
print(f"Unique prompts with 1-vs-3 confusion: {unique_prompts}")

# Category breakdown
cat_counts = Counter(c['category'] for c in one_three_cases)
print(f"\nBy category:")
for cat, count in cat_counts.most_common():
    print(f"  {cat}: {count}")

# Show examples
print(f"\nFirst 30 examples:")
for i, case in enumerate(one_three_cases[:30]):
    print(f"\n{i+1}. [{case['category']}] p{case['prompt_idx']} @{case['pos']}: "
          f"'{case['tok1']}' vs '{case['tok2']}' (m={case['margin']:.4f})")
    print(f"   Prompt: '{case['prompt']}'")
    print(f"   Context: '{case['context']}'")

# Is this exponent-related? Check if context contains ** or ^
exponent_count = 0
for case in one_three_cases:
    ctx = case['context']
    if '**' in ctx or '^' in ctx or '\\times' in ctx:
        exponent_count += 1

print(f"\n\nExponent/math context: {exponent_count}/{len(one_three_cases)} ({100*exponent_count/len(one_three_cases):.1f}%)")

# Check: is this the same few prompts repeated?
per_prompt = Counter(c['prompt_idx'] for c in one_three_cases)
print(f"\nTop 10 prompts by 1-vs-3 count:")
for pi, count in per_prompt.most_common(10):
    print(f"  p{pi} ({count}x): '{prompts[pi]['text'][:60]}...' [{prompts[pi]['category']}]")

# Now check: 2-vs-10 confusion (449 times) - what's going on there?
print("\n" + "=" * 70)
print("2-vs-10 and 10-vs-2 confusion analysis")
print("=" * 70)

two_ten_cases = []
for i in range(1999):
    path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
    if not os.path.exists(path):
        continue
    data = torch.load(path, map_location='cpu', weights_only=False)
    top5_vals = data['top5_values']
    top5_idx = data['top5_indices']
    margins = top5_vals[:, 0] - top5_vals[:, 1]
    gen_ids = data['generated_ids']
    input_len = len(data['input_ids'])

    for pos in range(len(margins)):
        if margins[pos] < MARGIN_THRESHOLD:
            t1 = decode(top5_idx[pos, 0].item())
            t2 = decode(top5_idx[pos, 1].item())
            pair = (t1.strip(), t2.strip())
            if pair == ('10', '2') or pair == ('2', '10') or pair == ('10', '3') or pair == ('3', '10'):
                actual_pos = input_len + pos
                ctx_start = max(0, actual_pos - 8)
                ctx_end = min(len(gen_ids), actual_pos + 5)
                context = tokenizer.decode(gen_ids[ctx_start:ctx_end])
                two_ten_cases.append({
                    'prompt_idx': i,
                    'pos': pos,
                    'pair': pair,
                    'margin': margins[pos].item(),
                    'context': context,
                    'prompt': prompts[i]['text'][:80]
                })

print(f"Total 2/3-vs-10 cases: {len(two_ten_cases)}")

# Show first 15
for i, case in enumerate(two_ten_cases[:15]):
    print(f"\n{i+1}. p{case['prompt_idx']} @{case['pos']}: {case['pair']} (m={case['margin']:.4f})")
    print(f"   Prompt: '{case['prompt']}'")
    print(f"   Context: '{case['context']}'")

# Check if 10 appears as the first token of "10^" or similar
ten_exp_count = 0
for case in two_ten_cases:
    if '10^' in case['context'] or '10 ^' in case['context'] or '10**' in case['context'] or '\\times 10' in case['context'] or '10\\' in case['context']:
        ten_exp_count += 1
print(f"\n10 in scientific notation context: {ten_exp_count}/{len(two_ten_cases)} ({100*ten_exp_count/len(two_ten_cases) if two_ten_cases else 0:.1f}%)")

print("\n\nDone!")
