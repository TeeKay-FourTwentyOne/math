"""
Experiment 1A: Token-level inference divergence measurement (H_nvidia).
Runs all 400 selected prompts twice with greedy decoding,
records full token sequences, computes within-platform nondeterminism baseline.
Also computes and saves a checkpoint hash for cross-platform verification.
"""
import os
import sys
import json
import time
import hashlib
import random
import logging
import numpy as np
import torch
from pathlib import Path

logging.getLogger("transformers").setLevel(logging.ERROR)

MODEL_NAME = "EleutherAI/pythia-410m-deduped"
MAX_NEW_TOKENS = 256
SEED = 42

ROOT = Path(__file__).parent
PROMPTS_DIR = ROOT / "prompts"
EXP1_DIR = ROOT / "experiment_1"
OUTPUTS_DIR = EXP1_DIR / "outputs"


def setup_determinism():
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
        torch.use_deterministic_algorithms(True, warn_only=True)
        return "warn_only"


def compute_checkpoint_hash(model):
    """SHA-256 of all model parameters (sorted by name) for cross-platform verification."""
    h = hashlib.sha256()
    for name in sorted(model.state_dict().keys()):
        h.update(name.encode("utf-8"))
        h.update(model.state_dict()[name].cpu().numpy().tobytes())
    return h.hexdigest()


def load_prompt_set():
    """Load the 400 selected prompts (200 experimental + 200 control)."""
    prompts = []
    for fname, set_name in [("experimental_200.jsonl", "experimental"),
                            ("control_200.jsonl", "control")]:
        with open(PROMPTS_DIR / fname, "r", encoding="utf-8") as f:
            for line in f:
                d = json.loads(line)
                d["set"] = set_name
                prompts.append(d)
    prompts.sort(key=lambda p: p["prompt_idx"])
    return prompts


def run_inference(model, tokenizer, prompts, device, label):
    """Run greedy inference on all prompts, return results list."""
    results = []
    t0 = time.time()
    for i, p in enumerate(prompts):
        inputs = tokenizer(p["text"], return_tensors="pt").to(device)
        with torch.no_grad():
            out = model.generate(
                **inputs,
                max_new_tokens=MAX_NEW_TOKENS,
                do_sample=False,
            )
        gen_ids = out[0].cpu().tolist()
        input_len = inputs.input_ids.shape[1]
        gen_text = tokenizer.decode(gen_ids[input_len:], skip_special_tokens=True)

        results.append({
            "prompt_idx": p["prompt_idx"],
            "category": p["category"],
            "set": p["set"],
            "fragility_score": p["fragility_score"],
            "text": p["text"],
            "input_ids": gen_ids[:input_len],
            "generated_ids": gen_ids[input_len:],
            "generated_text": gen_text,
        })

        if (i + 1) % 50 == 0 or (i + 1) == len(prompts):
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed
            eta = (len(prompts) - i - 1) / rate if rate > 0 else 0
            print(f"  [{label}] {i+1:4d}/{len(prompts)} | "
                  f"{rate:.1f} prompts/s | ETA {eta:.0f}s",
                  flush=True)

    elapsed = time.time() - t0
    print(f"  [{label}] Complete in {elapsed:.1f}s ({elapsed/len(prompts):.2f}s/prompt)",
          flush=True)
    return results


def save_results(results, path):
    with open(path, "w", encoding="utf-8") as f:
        for r in results:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")


def compare_runs(run1, run2, prompts):
    """Compare two runs for within-platform nondeterminism."""
    n = len(run1)
    exact = 0
    agreement_rates = []
    divergence_points = []
    per_set = {"experimental": {"total": 0, "match": 0},
               "control": {"total": 0, "match": 0}}

    for r1, r2 in zip(run1, run2):
        g1 = r1["generated_ids"]
        g2 = r2["generated_ids"]
        min_len = min(len(g1), len(g2))
        s = r1["set"]
        per_set[s]["total"] += 1

        if g1 == g2:
            exact += 1
            agreement_rates.append(1.0)
            divergence_points.append(None)
            per_set[s]["match"] += 1
        else:
            agree = sum(1 for a, b in zip(g1[:min_len], g2[:min_len]) if a == b)
            agreement_rates.append(agree / max(min_len, 1))
            div_pt = next((j for j in range(min_len) if g1[j] != g2[j]), min_len)
            divergence_points.append(div_pt)

    return {
        "total": n,
        "exact_match": exact,
        "match_rate": exact / n,
        "mean_agreement": float(np.mean(agreement_rates)),
        "per_set": per_set,
        "divergence_points": [p for p in divergence_points if p is not None],
    }


def main():
    print("=" * 70)
    print("EXPERIMENT 1A: Token-Level Inference (H_nvidia)")
    print("=" * 70)

    det_mode = setup_determinism()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}  Determinism: {det_mode}", flush=True)

    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    (EXP1_DIR / "analysis").mkdir(parents=True, exist_ok=True)

    # Load model
    from transformers import AutoModelForCausalLM, AutoTokenizer
    print(f"\nLoading {MODEL_NAME}...", flush=True)
    t0 = time.time()
    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token
    model = AutoModelForCausalLM.from_pretrained(
        MODEL_NAME, torch_dtype=torch.float32
    ).to(device)
    model.eval()
    print(f"Loaded in {time.time()-t0:.1f}s", flush=True)

    # Checkpoint hash
    print("Computing checkpoint hash...", flush=True)
    ckpt_hash = compute_checkpoint_hash(model)
    print(f"Checkpoint SHA-256: {ckpt_hash}", flush=True)
    with open(EXP1_DIR / "checkpoint_hash.json", "w") as f:
        json.dump({
            "model": MODEL_NAME,
            "sha256": ckpt_hash,
            "platform": "H_nvidia",
            "pytorch_version": torch.__version__,
            "dtype": "float32",
        }, f, indent=2)

    # Load prompts
    prompts = load_prompt_set()
    n_exp = sum(1 for p in prompts if p["set"] == "experimental")
    n_ctrl = sum(1 for p in prompts if p["set"] == "control")
    print(f"Loaded {len(prompts)} prompts: {n_exp} experimental + {n_ctrl} control",
          flush=True)

    # Verify deterministic generation works
    try:
        test_in = tokenizer("test", return_tensors="pt").to(device)
        with torch.no_grad():
            model.generate(**test_in, max_new_tokens=2, do_sample=False)
    except RuntimeError as e:
        if "deterministic" in str(e).lower():
            print("WARNING: strict determinism failed, using warn_only")
            torch.use_deterministic_algorithms(True, warn_only=True)
            det_mode = "warn_only"
        else:
            raise

    # --- Run 1 ---
    print(f"\n--- Run 1 ---", flush=True)
    setup_determinism()
    run1 = run_inference(model, tokenizer, prompts, device, "Run1")
    save_results(run1, OUTPUTS_DIR / "h_nvidia_run1.jsonl")

    # --- Run 2 ---
    print(f"\n--- Run 2 ---", flush=True)
    setup_determinism()
    run2 = run_inference(model, tokenizer, prompts, device, "Run2")
    save_results(run2, OUTPUTS_DIR / "h_nvidia_run2.jsonl")

    # --- Within-platform comparison ---
    print(f"\n{'='*70}")
    print("WITHIN-PLATFORM COMPARISON (Run 1 vs Run 2)")
    print(f"{'='*70}")
    stats = compare_runs(run1, run2, prompts)

    print(f"Exact match:          {stats['exact_match']}/{stats['total']} "
          f"({stats['match_rate']:.1%})")
    print(f"Mean token agreement: {stats['mean_agreement']:.6f}")

    for s in ["experimental", "control"]:
        ps = stats["per_set"][s]
        print(f"  {s:13s}: {ps['match']}/{ps['total']} exact match")

    if stats["exact_match"] < stats["total"]:
        dp = stats["divergence_points"]
        print(f"\nDiverged prompts: {len(dp)}")
        print(f"First divergence point: mean={np.mean(dp):.1f} "
              f"median={np.median(dp):.1f} min={min(dp)} max={max(dp)}")
    else:
        print("\nAll prompts identical across both runs (as expected with strict determinism).")

    with open(EXP1_DIR / "analysis" / "within_platform_baseline.json", "w") as f:
        json.dump({
            "platform": "H_nvidia",
            "determinism_mode": det_mode,
            "total_prompts": stats["total"],
            "exact_match": stats["exact_match"],
            "match_rate": stats["match_rate"],
            "mean_token_agreement": stats["mean_agreement"],
            "per_set": {s: stats["per_set"][s] for s in stats["per_set"]},
        }, f, indent=2)

    print(f"\n{'='*70}")
    print("EXPERIMENT 1A COMPLETE")
    print(f"{'='*70}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
