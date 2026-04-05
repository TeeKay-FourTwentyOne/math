"""
Experiment 0: Prompt Set Construction.
Runs all 2000 candidate prompts through Pythia-410M with greedy decoding,
records top-5 logits at each generation step, computes fragility scores,
and selects 200 experimental (high-fragility) + 200 control (low-fragility) prompts.

Supports resuming: skips prompts whose logit traces already exist on disk.
"""
import os
import sys
import json
import time
import csv
import random
import numpy as np
import torch
from pathlib import Path
from collections import Counter

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
MODEL_NAME = "EleutherAI/pythia-410m-deduped"
MAX_NEW_TOKENS = 256
FRAGILITY_THRESHOLD = 1.0   # logit margin between rank-1 and rank-2
NUM_EXPERIMENTAL = 200
NUM_CONTROL = 200
SEED = 42

ROOT = Path(__file__).parent
PROMPTS_DIR = ROOT / "prompts"
EXP0_DIR = ROOT / "experiment_0"
TRACES_DIR = EXP0_DIR / "logit_traces"


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
        try:
            torch.use_deterministic_algorithms(True, warn_only=True)
            return "warn_only"
        except Exception:
            return "off"


def load_prompts(path):
    prompts = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.strip():
                prompts.append(json.loads(line))
    return prompts


def score_one_prompt(model, tokenizer, text, device):
    """Run greedy generation, return token ids and top-5 logits per step."""
    inputs = tokenizer(text, return_tensors="pt").to(device)
    with torch.no_grad():
        out = model.generate(
            **inputs,
            max_new_tokens=MAX_NEW_TOKENS,
            do_sample=False,
            output_scores=True,
            return_dict_in_generate=True,
        )
    generated_ids = out.sequences[0].cpu()
    input_ids = inputs.input_ids[0].cpu()
    n_steps = len(out.scores)
    top5_values = torch.zeros(n_steps, 5)
    top5_indices = torch.zeros(n_steps, 5, dtype=torch.long)
    for step, step_scores in enumerate(out.scores):
        v, ix = torch.topk(step_scores[0].float().cpu(), k=5)
        top5_values[step] = v
        top5_indices[step] = ix
    return {
        "input_ids": input_ids,
        "generated_ids": generated_ids,
        "top5_values": top5_values,
        "top5_indices": top5_indices,
        "n_steps": n_steps,
    }


def compute_fragility(top5_values):
    """Compute fragility metrics from top-5 logit values (n_steps x 5)."""
    margins = (top5_values[:, 0] - top5_values[:, 1]).numpy()
    n = len(margins)
    fragile_mask = margins < FRAGILITY_THRESHOLD
    return {
        "fragility_score": float(fragile_mask.sum()) / n,
        "min_margin": float(margins.min()),
        "mean_margin": float(margins.mean()),
        "median_margin": float(np.median(margins)),
        "fragile_positions": int(fragile_mask.sum()),
        "total_positions": n,
    }


def select_prompts(results, n_exp=NUM_EXPERIMENTAL, n_ctrl=NUM_CONTROL):
    """Select top-fragile and least-fragile prompts, balanced across categories."""
    categories = sorted(set(r["category"] for r in results))
    n_cats = len(categories)
    per_cat_exp = n_exp // n_cats
    per_cat_ctrl = n_ctrl // n_cats

    experimental, control = [], []
    for cat in categories:
        cat_results = sorted(
            [r for r in results if r["category"] == cat],
            key=lambda r: r["fragility_score"],
            reverse=True,
        )
        experimental.extend(cat_results[:per_cat_exp])
        control.extend(cat_results[-per_cat_ctrl:])

    return experimental, control


def main():
    print("=" * 70)
    print("EXPERIMENT 0: Prompt Set Construction")
    print("=" * 70)

    # Setup
    det_mode = setup_determinism()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}  Determinism: {det_mode}")

    # Create output dirs
    TRACES_DIR.mkdir(parents=True, exist_ok=True)

    # Load prompts
    prompt_path = PROMPTS_DIR / "candidates_2000.jsonl"
    if not prompt_path.exists():
        print(f"ERROR: {prompt_path} not found. Run generate_prompts.py first.")
        return 1
    prompts = load_prompts(prompt_path)
    print(f"Loaded {len(prompts)} prompts from {prompt_path}")

    # Check for existing traces (resume support)
    existing_traces = set()
    for f in TRACES_DIR.glob("prompt_*.pt"):
        try:
            idx = int(f.stem.split("_")[1])
            existing_traces.add(idx)
        except (ValueError, IndexError):
            pass
    if existing_traces:
        print(f"Found {len(existing_traces)} existing traces -- will resume from where we left off")

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
    mem_gb = torch.cuda.memory_allocated() / 1e9 if device.type == "cuda" else 0
    print(f"Loaded in {t_load:.1f}s | VRAM: {mem_gb:.2f} GB")

    # Handle deterministic mode fallback during generation
    try:
        # Quick test generation
        test_ids = tokenizer("test", return_tensors="pt").to(device)
        with torch.no_grad():
            model.generate(**test_ids, max_new_tokens=2, do_sample=False)
    except RuntimeError as e:
        if "deterministic" in str(e).lower():
            print(f"WARNING: strict determinism failed: {e}")
            print("Falling back to warn_only mode")
            torch.use_deterministic_algorithms(True, warn_only=True)
            det_mode = "warn_only"
        else:
            raise

    # --- Main scoring loop ---
    print(f"\nScoring {len(prompts)} prompts (max {MAX_NEW_TOKENS} tokens each)...")
    results = []
    t_start = time.time()
    n_done = 0
    n_skipped = 0

    for prompt in prompts:
        idx = prompt["idx"]

        # Skip if trace already exists
        trace_path = TRACES_DIR / f"prompt_{idx:04d}.pt"
        if idx in existing_traces:
            # Load existing trace to get fragility metrics
            trace = torch.load(trace_path, weights_only=True)
            frag = compute_fragility(trace["top5_values"])
            results.append({
                "prompt_idx": idx,
                "category": prompt["category"],
                "text": prompt["text"],
                **frag,
            })
            n_skipped += 1
            continue

        # Score this prompt
        out = score_one_prompt(model, tokenizer, prompt["text"], device)

        # Save trace
        torch.save({
            "prompt_idx": idx,
            "input_ids": out["input_ids"],
            "generated_ids": out["generated_ids"],
            "top5_values": out["top5_values"],
            "top5_indices": out["top5_indices"],
        }, trace_path)

        # Compute fragility
        frag = compute_fragility(out["top5_values"])
        results.append({
            "prompt_idx": idx,
            "category": prompt["category"],
            "text": prompt["text"],
            **frag,
        })
        n_done += 1

        # Progress reporting
        total_done = n_done + n_skipped
        if n_done > 0 and (n_done % 50 == 0 or total_done == len(prompts)):
            elapsed = time.time() - t_start
            rate = n_done / elapsed
            remaining = (len(prompts) - total_done) / rate if rate > 0 else 0
            avg_frag = np.mean([r["fragility_score"] for r in results[-50:]])
            print(
                f"  [{total_done:4d}/{len(prompts)}] "
                f"{rate:.1f} prompts/s | "
                f"ETA {remaining/60:.0f} min | "
                f"avg_fragility={avg_frag:.3f}"
            )

        # Periodic cache clearing
        if n_done % 200 == 0 and device.type == "cuda":
            torch.cuda.empty_cache()

    elapsed_total = time.time() - t_start
    print(f"\nScoring complete: {n_done} new + {n_skipped} resumed in {elapsed_total/60:.1f} min")

    # --- Save fragility scores ---
    scores_path = PROMPTS_DIR / "fragility_scores.csv"
    results_sorted = sorted(results, key=lambda r: r["prompt_idx"])
    with open(scores_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "prompt_idx", "category", "fragility_score", "min_margin",
            "mean_margin", "median_margin", "fragile_positions", "total_positions", "text",
        ])
        writer.writeheader()
        for r in results_sorted:
            writer.writerow(r)
    print(f"Saved fragility scores to {scores_path}")

    # --- Summary statistics ---
    print(f"\n{'=' * 70}")
    print("FRAGILITY SUMMARY")
    print(f"{'=' * 70}")
    all_frag = [r["fragility_score"] for r in results]
    all_min = [r["min_margin"] for r in results]
    print(f"Fragility score: mean={np.mean(all_frag):.4f}  median={np.median(all_frag):.4f}  "
          f"max={max(all_frag):.4f}  min={min(all_frag):.4f}")
    print(f"Min margin:      mean={np.mean(all_min):.4f}  median={np.median(all_min):.4f}")

    cats = sorted(set(r["category"] for r in results))
    for cat in cats:
        cat_frag = [r["fragility_score"] for r in results if r["category"] == cat]
        print(f"  {cat:15s}: mean_fragility={np.mean(cat_frag):.4f}  max={max(cat_frag):.4f}")

    # --- Select experimental and control sets ---
    experimental, control = select_prompts(results)

    exp_path = PROMPTS_DIR / "experimental_200.jsonl"
    ctrl_path = PROMPTS_DIR / "control_200.jsonl"

    for prompts_subset, path in [(experimental, exp_path), (control, ctrl_path)]:
        # Re-attach full prompt text from the original candidates
        with open(path, "w", encoding="utf-8") as f:
            for r in sorted(prompts_subset, key=lambda x: x["prompt_idx"]):
                f.write(json.dumps(r, ensure_ascii=False) + "\n")

    print(f"\nSelected {len(experimental)} experimental prompts -> {exp_path}")
    print(f"Selected {len(control)} control prompts -> {ctrl_path}")

    # Show top-10 most fragile
    top_fragile = sorted(results, key=lambda r: r["fragility_score"], reverse=True)[:10]
    print(f"\nTop 10 most fragile prompts:")
    for r in top_fragile:
        text = r["text"][:70] + ("..." if len(r["text"]) > 70 else "")
        print(f"  [{r['category']:12s}] frag={r['fragility_score']:.3f}  min_margin={r['min_margin']:.3f}  {text}")

    print(f"\n{'=' * 70}")
    print("EXPERIMENT 0 COMPLETE")
    print(f"{'=' * 70}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
