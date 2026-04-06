"""
Counterfactual flip: What happens if the closest-call token had flipped?

Prompt idx=334: "In what year did the Great Fire of London occur?"
Position 61: noise was 84.4% of margin (the closest call in the entire dataset).

We take the original generation up to position 61, swap rank-1 for rank-2,
then let the model continue generating. This shows the ACTUAL downstream
effect of a single hardware-induced token flip.

Also do this for the top 5 closest calls to see the pattern.
"""

import os
import sys
import json
import torch
import numpy as np
import random
import logging

logging.getLogger("transformers").setLevel(logging.ERROR)

SEED = 42
MODEL_NAME = "EleutherAI/pythia-410m-deduped"
MAX_NEW_TOKENS = 256

def setup_determinism():
    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(SEED)
        os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.use_deterministic_algorithms(True)

setup_determinism()
device = torch.device("cuda")

from transformers import AutoModelForCausalLM, AutoTokenizer

print("Loading model...", flush=True)
tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
if tokenizer.pad_token is None:
    tokenizer.pad_token = tokenizer.eos_token
model = AutoModelForCausalLM.from_pretrained(MODEL_NAME, torch_dtype=torch.float32).to(device)
model.eval()
print("Model loaded.", flush=True)

# Load the NVIDIA run data
nvidia_data = {}
with open("C:/Projects/math/non_deterministic/experiment_1/outputs/h_nvidia_run1.jsonl",
          encoding="utf-8") as f:
    for line in f:
        d = json.loads(line)
        nvidia_data[d["prompt_idx"]] = d

mac_data = {}
with open("C:/Projects/math/non_deterministic/experiment_1/outputs/h_mac_run1.jsonl",
          encoding="utf-8") as f:
    for line in f:
        d = json.loads(line)
        mac_data[d["prompt_idx"]] = d

# The closest calls from the comparison
closest_calls = [
    {"idx": 334, "pos": 61, "ratio": 0.844,
     "text": "In what year did the Great Fire of London occur?"},
    {"idx": 801, "pos": 164, "ratio": 0.769,
     "text": "What is 22 times 58?"},
    {"idx": 357, "pos": 23, "ratio": 0.525,
     "text": "In what year did the signing of the Treaty of Westphalia occur?"},
    {"idx": 1050, "pos": 80, "ratio": 0.488,
     "text": "What is 15% of 1751?"},
    {"idx": 1831, "pos": 23, "ratio": 0.466,
     "text": "Can you recommend a good recipe for a beginner cook?"},
]

print("=" * 70)
print("COUNTERFACTUAL FLIP EXPERIMENT")
print("What happens when the closest-call token actually flips?")
print("=" * 70)

for case in closest_calls:
    idx = case["idx"]
    flip_pos = case["pos"]  # position in GENERATED tokens (0-indexed)

    nd = nvidia_data[idx]
    md = mac_data[idx]
    n_input = nd["n_input_tokens"]

    # Full original sequence (input + generated)
    original_full = nd["generated_ids"]  # includes input tokens
    original_gen = original_full[n_input:]  # generated only

    # What are the competing tokens at the flip position?
    rank1_id = nd["top5_indices"][flip_pos][0]
    rank2_id = nd["top5_indices"][flip_pos][1]
    rank1_val = nd["top5_values"][flip_pos][0]
    rank2_val = nd["top5_values"][flip_pos][1]
    margin = rank1_val - rank2_val

    rank1_tok = tokenizer.decode([rank1_id])
    rank2_tok = tokenizer.decode([rank2_id])

    print(f"\n{'=' * 70}")
    print(f"PROMPT idx={idx}: \"{case['text']}\"")
    print(f"Flip position: {flip_pos} (in generated tokens)")
    print(f"Noise/margin ratio: {case['ratio']:.3f}")
    print(f"Token contest: '{rank1_tok}' (logit {rank1_val:.4f}) vs "
          f"'{rank2_tok}' (logit {rank2_val:.4f}), margin={margin:.6f}")
    print(f"{'=' * 70}")

    # Build the counterfactual prefix: original up to flip_pos, then rank-2 token
    # The prefix is: input_tokens + generated_tokens[0:flip_pos] + rank2_token
    prefix = original_full[:n_input + flip_pos] + [rank2_id]

    # Generate the rest from the flipped prefix
    setup_determinism()
    prefix_tensor = torch.tensor([prefix], dtype=torch.long, device=device)
    remaining_tokens = MAX_NEW_TOKENS - flip_pos - 1  # how many more to generate

    with torch.no_grad():
        counterfactual_out = model.generate(
            prefix_tensor,
            max_new_tokens=remaining_tokens,
            do_sample=False,
        )

    counterfactual_full = counterfactual_out[0].cpu().tolist()
    counterfactual_gen = counterfactual_full[n_input:]  # generated only

    # Decode both
    original_text = tokenizer.decode(original_gen, skip_special_tokens=True)
    counterfactual_text = tokenizer.decode(counterfactual_gen, skip_special_tokens=True)

    # Compare token-by-token
    min_len = min(len(original_gen), len(counterfactual_gen))
    token_matches = sum(1 for a, b in zip(original_gen[:min_len], counterfactual_gen[:min_len])
                        if a == b)
    agreement_rate = token_matches / max(min_len, 1)

    # Find reconvergence: after the flip, how many tokens match?
    post_flip_matches = 0
    post_flip_total = min_len - flip_pos - 1
    reconverge_pos = None
    consecutive_match = 0
    for j in range(flip_pos + 1, min_len):
        if original_gen[j] == counterfactual_gen[j]:
            post_flip_matches += 1
            consecutive_match += 1
            if consecutive_match >= 5 and reconverge_pos is None:
                reconverge_pos = j - 4  # start of the matching run
        else:
            consecutive_match = 0

    post_agreement = post_flip_matches / max(post_flip_total, 1)

    # Show context around the flip
    ctx_start = max(0, flip_pos - 5)
    ctx_end = min(len(original_gen), flip_pos + 20)
    orig_context = tokenizer.decode(original_gen[ctx_start:ctx_end], skip_special_tokens=True)
    cf_context = tokenizer.decode(counterfactual_gen[ctx_start:min(len(counterfactual_gen), ctx_end)],
                                   skip_special_tokens=True)

    print(f"\n--- Token statistics ---")
    print(f"Original gen length:       {len(original_gen)}")
    print(f"Counterfactual gen length: {len(counterfactual_gen)}")
    print(f"Overall token agreement:   {agreement_rate:.4f} ({token_matches}/{min_len})")
    print(f"Pre-flip agreement:        1.000 ({flip_pos}/{flip_pos}) [identical by construction]")
    print(f"Post-flip agreement:       {post_agreement:.4f} ({post_flip_matches}/{post_flip_total})")
    if reconverge_pos is not None:
        print(f"Reconvergence at:          position {reconverge_pos} ({reconverge_pos - flip_pos} tokens after flip)")
    else:
        print(f"Reconvergence:             NEVER (no 5+ consecutive matching tokens after flip)")

    print(f"\n--- Context around flip (pos {flip_pos}) ---")
    print(f"  Original:       ...{orig_context}...")
    print(f"  Counterfactual: ...{cf_context}...")

    print(f"\n--- Full generated text comparison ---")
    # Show first 300 chars of each, highlighting the divergence region
    print(f"  ORIGINAL (first 300 chars):")
    print(f"    {original_text[:300]}")
    print(f"  COUNTERFACTUAL (first 300 chars):")
    print(f"    {counterfactual_text[:300]}")

    # Levenshtein-like: count character edits
    from difflib import SequenceMatcher
    sm = SequenceMatcher(None, original_text[:500], counterfactual_text[:500])
    text_similarity = sm.ratio()
    print(f"\n  Text similarity (first 500 chars): {text_similarity:.4f}")
