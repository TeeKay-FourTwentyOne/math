"""
Experiment 1A: Inference Divergence — Token-level divergence measurement.

Runs all 400 prompts (200 experimental + 200 control) through Pythia-410M
with greedy decoding at FP32. Performs 2 runs per platform to establish a
within-platform determinism baseline, then outputs are compared cross-platform.

Each run produces a JSONL file with one record per prompt containing:
  - prompt_idx, category, text
  - generated_token_ids (full sequence including input)
  - top5_values, top5_indices per generation step
  - timing metadata

Usage:
    python experiment_1.py              # auto-detect platform, run both runs
    python experiment_1.py --run 1      # run only run 1
    python experiment_1.py --run 2      # run only run 2
    python experiment_1.py --device cpu  # force CPU (default: auto)
    python experiment_1.py --device mps  # force MPS (if available)
"""
import os
import sys
import json
import time
import argparse
import random
import hashlib
import logging
import numpy as np
import torch
from pathlib import Path

# Suppress noisy warnings
logging.getLogger("transformers").setLevel(logging.ERROR)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
MODEL_NAME = "EleutherAI/pythia-410m-deduped"
MAX_NEW_TOKENS = 256
SEED = 42

ROOT = Path(__file__).parent
PROMPTS_DIR = ROOT / "prompts"
EXP1_DIR = ROOT / "experiment_1"
OUTPUTS_DIR = EXP1_DIR / "outputs"


def setup_determinism(device_type):
    """Set all determinism flags. Returns mode string."""
    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(SEED)
        os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    try:
        torch.use_deterministic_algorithms(True)
        return "strict"
    except Exception:
        try:
            torch.use_deterministic_algorithms(True, warn_only=True)
            return "warn_only"
        except Exception:
            return "off"


def select_device(requested):
    """Select compute device. Returns torch.device and a label string."""
    if requested == "mps":
        if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            return torch.device("mps"), "mps"
        print("WARNING: MPS requested but not available, falling back to CPU")
        return torch.device("cpu"), "cpu"
    elif requested == "cuda":
        if torch.cuda.is_available():
            return torch.device("cuda"), "cuda"
        print("WARNING: CUDA requested but not available, falling back to CPU")
        return torch.device("cpu"), "cpu"
    elif requested == "cpu":
        return torch.device("cpu"), "cpu"
    else:
        # Auto-detect: prefer CUDA > CPU (skip MPS for determinism safety)
        if torch.cuda.is_available():
            return torch.device("cuda"), "cuda"
        return torch.device("cpu"), "cpu"


def detect_platform_label():
    """Return platform label matching hardware_manifest convention."""
    if torch.cuda.is_available():
        return "h_nvidia"
    elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return "h_mac"
    return "h_cpu"


def load_prompts():
    """Load the 400-prompt set (200 experimental + 200 control)."""
    prompts = []
    for name in ["experimental_200.jsonl", "control_200.jsonl"]:
        path = PROMPTS_DIR / name
        if not path.exists():
            print(f"ERROR: {path} not found. Run experiment_0.py first.")
            sys.exit(1)
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                if line.strip():
                    prompts.append(json.loads(line))
    # Sort by prompt_idx for deterministic ordering
    prompts.sort(key=lambda p: p["prompt_idx"])
    return prompts


def hash_model_weights(model):
    """SHA-256 hash of all model parameter tensors (for cross-platform verification)."""
    h = hashlib.sha256()
    for name, param in sorted(model.named_parameters()):
        h.update(name.encode())
        h.update(param.detach().cpu().numpy().tobytes())
    return h.hexdigest()


def run_inference(model, tokenizer, prompts, device, run_id, platform_label):
    """Run greedy inference on all prompts, save results to JSONL."""
    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUTPUTS_DIR / f"{platform_label}_run{run_id}.jsonl"

    # Check for existing partial results (resume support)
    existing = set()
    if out_path.exists():
        with open(out_path, "r", encoding="utf-8") as f:
            for line in f:
                if line.strip():
                    rec = json.loads(line)
                    existing.add(rec["prompt_idx"])
        if len(existing) == len(prompts):
            print(f"  Run {run_id} already complete ({len(existing)} prompts). Skipping.")
            return out_path
        print(f"  Resuming run {run_id}: {len(existing)}/{len(prompts)} already done")

    n_new = 0
    t_start = time.time()

    with open(out_path, "a", encoding="utf-8") as fout:
        for i, prompt in enumerate(prompts):
            idx = prompt["prompt_idx"]
            if idx in existing:
                continue

            # Greedy generation with logit recording
            inputs = tokenizer(prompt["text"], return_tensors="pt").to(device)
            with torch.no_grad():
                out = model.generate(
                    **inputs,
                    max_new_tokens=MAX_NEW_TOKENS,
                    do_sample=False,
                    output_scores=True,
                    return_dict_in_generate=True,
                )

            generated_ids = out.sequences[0].cpu().tolist()
            n_steps = len(out.scores)
            top5_values = []
            top5_indices = []
            for step_scores in out.scores:
                v, ix = torch.topk(step_scores[0].float().cpu(), k=5)
                top5_values.append(v.tolist())
                top5_indices.append(ix.tolist())

            record = {
                "prompt_idx": idx,
                "category": prompt["category"],
                "text": prompt["text"],
                "generated_ids": generated_ids,
                "n_input_tokens": inputs.input_ids.shape[1],
                "n_generated_tokens": n_steps,
                "top5_values": top5_values,
                "top5_indices": top5_indices,
            }
            fout.write(json.dumps(record, ensure_ascii=False) + "\n")
            fout.flush()
            n_new += 1

            # Progress
            total_done = len(existing) + n_new
            if n_new % 25 == 0 or total_done == len(prompts):
                elapsed = time.time() - t_start
                rate = n_new / elapsed if elapsed > 0 else 0
                remaining = (len(prompts) - total_done) / rate if rate > 0 else 0
                print(
                    f"  Run {run_id}: [{total_done:3d}/{len(prompts)}] "
                    f"{rate:.2f} prompts/s | "
                    f"ETA {remaining/60:.0f} min"
                )

    elapsed_total = time.time() - t_start
    print(f"  Run {run_id} complete: {n_new} new prompts in {elapsed_total/60:.1f} min")
    return out_path


def verify_within_platform(path1, path2):
    """Compare two runs on the same platform. They should be identical."""
    records1, records2 = {}, {}
    with open(path1, "r") as f:
        for line in f:
            r = json.loads(line)
            records1[r["prompt_idx"]] = r
    with open(path2, "r") as f:
        for line in f:
            r = json.loads(line)
            records2[r["prompt_idx"]] = r

    common = sorted(set(records1) & set(records2))
    n_match = 0
    n_diverged = 0
    first_divergences = []

    for idx in common:
        ids1 = records1[idx]["generated_ids"]
        ids2 = records2[idx]["generated_ids"]
        if ids1 == ids2:
            n_match += 1
        else:
            n_diverged += 1
            # Find first divergence point
            for pos in range(min(len(ids1), len(ids2))):
                if ids1[pos] != ids2[pos]:
                    first_divergences.append((idx, pos, records1[idx]["category"]))
                    break

    print(f"\n{'=' * 60}")
    print(f"WITHIN-PLATFORM VERIFICATION (Run 1 vs Run 2)")
    print(f"{'=' * 60}")
    print(f"Prompts compared: {len(common)}")
    print(f"Exact match:      {n_match}/{len(common)} ({100*n_match/len(common):.1f}%)")
    print(f"Diverged:         {n_diverged}/{len(common)}")
    if first_divergences:
        print(f"\nDivergence points:")
        for idx, pos, cat in first_divergences[:10]:
            print(f"  prompt {idx} [{cat}]: first diff at token {pos}")
    else:
        print("Perfect within-platform determinism confirmed.")


def main():
    parser = argparse.ArgumentParser(description="Experiment 1A: Inference Divergence")
    parser.add_argument("--run", type=int, choices=[1, 2], default=None,
                        help="Run only this run number (default: both)")
    parser.add_argument("--device", choices=["cpu", "mps", "cuda", "auto"], default="auto",
                        help="Compute device (default: auto)")
    args = parser.parse_args()

    print("=" * 70)
    print("EXPERIMENT 1A: Inference Divergence — Token-level Measurement")
    print("=" * 70)

    # Device and determinism
    device, device_label = select_device(args.device)
    det_mode = setup_determinism(device_label)
    platform_label = detect_platform_label()
    print(f"Platform:     {platform_label}")
    print(f"Device:       {device} ({device_label})")
    print(f"Determinism:  {det_mode}")

    # Load prompts
    prompts = load_prompts()
    print(f"Prompts:      {len(prompts)} (200 experimental + 200 control)")

    # Load model
    from transformers import AutoModelForCausalLM, AutoTokenizer
    print(f"\nLoading {MODEL_NAME}...")
    t0 = time.time()
    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token
    model = AutoModelForCausalLM.from_pretrained(
        MODEL_NAME, torch_dtype=torch.float32
    ).to(device)
    model.eval()
    t_load = time.time() - t0
    print(f"Loaded in {t_load:.1f}s")

    # Verify checkpoint integrity
    print("Hashing model weights (for cross-platform verification)...")
    weight_hash = hash_model_weights(model)
    print(f"Weight SHA-256: {weight_hash[:16]}...")

    # Handle deterministic mode fallback
    try:
        test_ids = tokenizer("test", return_tensors="pt").to(device)
        with torch.no_grad():
            model.generate(**test_ids, max_new_tokens=2, do_sample=False)
    except RuntimeError as e:
        if "deterministic" in str(e).lower():
            print(f"WARNING: strict determinism failed during generation, falling back to warn_only")
            torch.use_deterministic_algorithms(True, warn_only=True)
            det_mode = "warn_only"
        else:
            raise

    # Save run metadata
    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    meta = {
        "platform_label": platform_label,
        "device": device_label,
        "determinism_mode": det_mode,
        "model": MODEL_NAME,
        "weight_hash_sha256": weight_hash,
        "max_new_tokens": MAX_NEW_TOKENS,
        "n_prompts": len(prompts),
        "precision": "fp32",
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
    meta_path = OUTPUTS_DIR / f"{platform_label}_metadata.json"
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")
    print(f"Metadata saved to {meta_path}")

    # Run inference
    runs_to_do = [1, 2] if args.run is None else [args.run]
    paths = {}

    for run_id in runs_to_do:
        print(f"\n{'=' * 70}")
        print(f"RUN {run_id}")
        print(f"{'=' * 70}")
        # Re-seed before each run for reproducibility
        setup_determinism(device_label)
        if det_mode == "warn_only":
            torch.use_deterministic_algorithms(True, warn_only=True)
        paths[run_id] = run_inference(model, tokenizer, prompts, device, run_id, platform_label)

    # Within-platform verification
    if len(paths) == 2:
        verify_within_platform(paths[1], paths[2])

    print(f"\n{'=' * 70}")
    print("EXPERIMENT 1A COMPLETE")
    print(f"{'=' * 70}")
    print(f"Outputs in: {OUTPUTS_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
