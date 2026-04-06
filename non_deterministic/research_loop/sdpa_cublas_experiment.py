"""
Experiments 1 & 2: SDPA backend and cuBLAS isolation.

Runs a subset of 50 prompts under different configurations to measure
the contribution of each component to logit divergence.

Configurations:
  A) math SDPA + deterministic cuBLAS     (baseline — our existing run)
  B) flash SDPA + deterministic cuBLAS     (isolates attention tiling)
  C) mem_efficient SDPA + det cuBLAS       (different attention tiling)
  D) math SDPA + non-deterministic cuBLAS  (isolates cuBLAS algorithm)
  E) flash SDPA + non-deterministic cuBLAS (both changed)

We keep torch.use_deterministic_algorithms(True, warn_only=True) so
non-deterministic ops warn but don't error. For configs A, we use strict.
"""

import os
import sys
import json
import time
import random
import logging
import numpy as np
import torch
from pathlib import Path
from contextlib import contextmanager

logging.getLogger("transformers").setLevel(logging.ERROR)

MODEL_NAME = "EleutherAI/pythia-410m-deduped"
MAX_NEW_TOKENS = 256
SEED = 42
N_PROMPTS = 50  # subset for speed

ROOT = Path(__file__).parent.parent
PROMPTS_DIR = ROOT / "prompts"

device = torch.device("cuda")

# Load prompts (first 50 from experimental + control combined)
prompts = []
for name in ["experimental_200.jsonl", "control_200.jsonl"]:
    with open(PROMPTS_DIR / name, "r", encoding="utf-8") as f:
        for line in f:
            if line.strip():
                prompts.append(json.loads(line))
prompts.sort(key=lambda p: p["prompt_idx"])
prompts = prompts[:N_PROMPTS]

print(f"Using {len(prompts)} prompts for comparison", flush=True)

# Load model
from transformers import AutoModelForCausalLM, AutoTokenizer
print("Loading model...", flush=True)
tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
if tokenizer.pad_token is None:
    tokenizer.pad_token = tokenizer.eos_token
model = AutoModelForCausalLM.from_pretrained(MODEL_NAME, torch_dtype=torch.float32).to(device)
model.eval()
print("Model loaded.", flush=True)


def full_determinism():
    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    torch.cuda.manual_seed_all(SEED)
    os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.use_deterministic_algorithms(True)


def warn_only_determinism():
    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    torch.cuda.manual_seed_all(SEED)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.use_deterministic_algorithms(True, warn_only=True)


def non_det_cublas():
    """Deterministic everything EXCEPT cuBLAS."""
    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    torch.cuda.manual_seed_all(SEED)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    # Remove the workspace config to allow non-deterministic cuBLAS
    if "CUBLAS_WORKSPACE_CONFIG" in os.environ:
        del os.environ["CUBLAS_WORKSPACE_CONFIG"]
    torch.use_deterministic_algorithms(True, warn_only=True)


@contextmanager
def sdpa_backend(backend_name):
    """Force a specific SDPA backend."""
    if backend_name == "math":
        with torch.backends.cuda.sdp_kernel(
            enable_flash=False, enable_math=True, enable_mem_efficient=False
        ):
            yield
    elif backend_name == "flash":
        with torch.backends.cuda.sdp_kernel(
            enable_flash=True, enable_math=False, enable_mem_efficient=False
        ):
            yield
    elif backend_name == "mem_efficient":
        with torch.backends.cuda.sdp_kernel(
            enable_flash=False, enable_math=False, enable_mem_efficient=True
        ):
            yield
    elif backend_name == "default":
        yield  # no override
    else:
        raise ValueError(f"Unknown backend: {backend_name}")


def run_config(prompts, config_name, setup_fn, sdpa_name):
    """Run prompts under a specific configuration, return per-step logit data."""
    setup_fn()

    results = {}
    t0 = time.time()

    for i, prompt in enumerate(prompts):
        idx = prompt["prompt_idx"]
        inputs = tokenizer(prompt["text"], return_tensors="pt").to(device)

        try:
            with torch.no_grad():
                with sdpa_backend(sdpa_name):
                    out = model.generate(
                        **inputs,
                        max_new_tokens=MAX_NEW_TOKENS,
                        do_sample=False,
                        output_scores=True,
                        return_dict_in_generate=True,
                    )
        except Exception as e:
            print(f"  Config {config_name}: FAILED on prompt {idx}: {e}", flush=True)
            return None

        gen_ids = out.sequences[0].cpu().tolist()
        n_input = inputs.input_ids.shape[1]
        top5_vals = []
        top5_idxs = []
        for step_scores in out.scores:
            v, ix = torch.topk(step_scores[0].float().cpu(), k=5)
            top5_vals.append(v.tolist())
            top5_idxs.append(ix.tolist())

        results[idx] = {
            "gen_ids": gen_ids[n_input:],  # generated only
            "top5_values": top5_vals,
            "top5_indices": top5_idxs,
        }

        if (i + 1) % 25 == 0:
            elapsed = time.time() - t0
            print(f"  [{config_name}] {i+1}/{len(prompts)} in {elapsed:.0f}s", flush=True)

    elapsed = time.time() - t0
    print(f"  [{config_name}] Done: {len(results)} prompts in {elapsed:.0f}s", flush=True)
    return results


def compare_configs(baseline, other, name_b, name_o):
    """Compare two config results at the logit and token level."""
    common = sorted(set(baseline) & set(other))
    token_match = 0
    all_diffs = []
    all_ratios = []
    rank1_mismatches = 0
    total_positions = 0
    max_ratio = 0
    max_ratio_info = None

    for idx in common:
        b = baseline[idx]
        o = other[idx]

        # Token match
        if b["gen_ids"] == o["gen_ids"]:
            token_match += 1

        # Logit comparison
        b_vals = np.array(b["top5_values"])
        o_vals = np.array(o["top5_values"])
        b_idxs = np.array(b["top5_indices"])
        o_idxs = np.array(o["top5_indices"])
        n_steps = min(len(b_vals), len(o_vals))

        for s in range(n_steps):
            total_positions += 1
            diff = abs(b_vals[s, 0] - o_vals[s, 0])
            all_diffs.append(diff)

            if b_idxs[s, 0] != o_idxs[s, 0]:
                rank1_mismatches += 1

            margin = min(b_vals[s, 0] - b_vals[s, 1], o_vals[s, 0] - o_vals[s, 1])
            if margin > 0:
                ratio = diff / margin
                all_ratios.append(ratio)
                if ratio > max_ratio:
                    max_ratio = ratio
                    max_ratio_info = {"idx": idx, "pos": s, "diff": diff, "margin": margin}

    diffs = np.array(all_diffs)
    ratios = np.array(all_ratios) if all_ratios else np.array([0])

    print(f"\n  {name_b} vs {name_o}:")
    print(f"    Token match:      {token_match}/{len(common)}")
    print(f"    Rank-1 mismatches:{rank1_mismatches}/{total_positions}")
    print(f"    Logit noise: mean={diffs.mean():.2e}  median={np.median(diffs):.2e}  "
          f"max={diffs.max():.2e}  zeros={int((diffs==0).sum())}/{len(diffs)}")
    print(f"    Noise/margin: mean={ratios.mean():.2e}  max={ratios.max():.4f}  "
          f">0.1={int((ratios>0.1).sum())}  >0.5={int((ratios>0.5).sum())}")
    if max_ratio_info:
        print(f"    Closest call: idx={max_ratio_info['idx']} pos={max_ratio_info['pos']} "
              f"ratio={max_ratio:.4f}")


# ============================================================
# Run all configurations
# ============================================================
configs = [
    ("A_math_det", full_determinism, "math",
     "Baseline: math SDPA + deterministic cuBLAS"),
    ("B_flash_det", warn_only_determinism, "flash",
     "Flash SDPA + deterministic cuBLAS"),
    ("C_memeff_det", warn_only_determinism, "mem_efficient",
     "MemEfficient SDPA + deterministic cuBLAS"),
    ("D_math_nondet", non_det_cublas, "math",
     "Math SDPA + non-deterministic cuBLAS"),
    ("E_flash_nondet", non_det_cublas, "flash",
     "Flash SDPA + non-deterministic cuBLAS"),
]

print("=" * 70)
print("SDPA & cuBLAS ISOLATION EXPERIMENT")
print("=" * 70)

all_results = {}
for config_name, setup_fn, sdpa_name, description in configs:
    print(f"\n--- Config {config_name}: {description} ---", flush=True)
    result = run_config(prompts, config_name, setup_fn, sdpa_name)
    if result is not None:
        all_results[config_name] = result
    else:
        print(f"  SKIPPED (failed)")

# ============================================================
# Compare all pairs against baseline
# ============================================================
print("\n" + "=" * 70)
print("COMPARISON RESULTS")
print("=" * 70)

baseline_name = "A_math_det"
if baseline_name in all_results:
    for config_name in all_results:
        if config_name == baseline_name:
            continue
        compare_configs(all_results[baseline_name], all_results[config_name],
                       baseline_name, config_name)

# Also: run baseline twice to confirm it's deterministic
print(f"\n--- Determinism check: running baseline twice ---", flush=True)
baseline2 = run_config(prompts, "A_math_det_v2", full_determinism, "math")
if baseline2:
    compare_configs(all_results[baseline_name], baseline2, "A_math_det", "A_math_det_v2")

# ============================================================
# Summary table
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY: What each component contributes")
print("=" * 70)
print("""
Config changes vs baseline:
  B: SDPA math -> flash       (attention tiling only)
  C: SDPA math -> mem_eff     (different attention tiling)
  D: cuBLAS det -> non-det    (GEMM algorithm only)
  E: Both SDPA + cuBLAS       (everything changed)

If B and D are independent, E's noise should be ~sqrt(B^2 + D^2).
""")
