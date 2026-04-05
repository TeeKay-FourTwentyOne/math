"""
Experiment 1B: Activation-level divergence profiling (H_nvidia).
Records all transformer layer outputs for 50 prompts (25 experimental + 25 control)
via forward hooks. Uses the token sequences from Experiment 1A Run 1.

Each prompt produces a .pt file containing:
  - embed: (seq_len, hidden_size) — embedding output
  - layer_0..layer_23: (seq_len, hidden_size) — transformer block outputs
  - final_norm: (seq_len, hidden_size) — final layer norm output
  - metadata: prompt_idx, input_length, total_length, set, category
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

logging.getLogger("transformers").setLevel(logging.ERROR)

MODEL_NAME = "EleutherAI/pythia-410m-deduped"
SEED = 42
N_PER_SET = 25  # 25 experimental + 25 control = 50 total

ROOT = Path(__file__).parent
PROMPTS_DIR = ROOT / "prompts"
EXP1_DIR = ROOT / "experiment_1"
OUTPUTS_DIR = EXP1_DIR / "outputs"
ACTIVATIONS_DIR = EXP1_DIR / "activations" / "h_nvidia"


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


def select_activation_subset():
    """Top 25 experimental (most fragile) + bottom 25 control (least fragile)."""
    selected = []
    for fname, set_name, reverse in [
        ("experimental_200.jsonl", "experimental", True),
        ("control_200.jsonl", "control", False),
    ]:
        prompts = []
        with open(PROMPTS_DIR / fname, "r", encoding="utf-8") as f:
            for line in f:
                d = json.loads(line)
                d["set"] = set_name
                prompts.append(d)
        prompts.sort(key=lambda p: p["fragility_score"], reverse=reverse)
        selected.extend(prompts[:N_PER_SET])
    selected.sort(key=lambda p: p["prompt_idx"])
    return selected


def load_run1_results():
    """Load Experiment 1A Run 1 results, indexed by prompt_idx."""
    results = {}
    with open(OUTPUTS_DIR / "h_nvidia_run1.jsonl", "r", encoding="utf-8") as f:
        for line in f:
            d = json.loads(line)
            results[d["prompt_idx"]] = d
    return results


def record_activations(model, token_ids, device):
    """Single forward pass with hooks to capture all layer outputs."""
    activations = {}
    hooks = []

    # Embedding
    def embed_hook(module, inp, out):
        activations["embed"] = out.detach().cpu().float()
    hooks.append(model.gpt_neox.embed_in.register_forward_hook(embed_hook))

    # Transformer layers
    for i, layer in enumerate(model.gpt_neox.layers):
        name = f"layer_{i}"
        def make_hook(n):
            def hook_fn(module, inp, out):
                activations[n] = out[0].detach().cpu().float()
            return hook_fn
        hooks.append(layer.register_forward_hook(make_hook(name)))

    # Final layer norm
    def norm_hook(module, inp, out):
        activations["final_norm"] = out.detach().cpu().float()
    hooks.append(model.gpt_neox.final_layer_norm.register_forward_hook(norm_hook))

    # Forward pass
    input_tensor = torch.tensor([token_ids], dtype=torch.long, device=device)
    with torch.no_grad():
        model(input_tensor)

    # Cleanup
    for h in hooks:
        h.remove()

    # Squeeze batch dim: (1, seq_len, hidden) -> (seq_len, hidden)
    return {k: v.squeeze(0) for k, v in activations.items()}


def main():
    print("=" * 70)
    print("EXPERIMENT 1B: Activation-Level Profiling (H_nvidia)")
    print("=" * 70)

    det_mode = setup_determinism()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}  Determinism: {det_mode}", flush=True)

    # Check prerequisites
    run1_path = OUTPUTS_DIR / "h_nvidia_run1.jsonl"
    if not run1_path.exists():
        print(f"ERROR: Run experiment_1a.py first (need {run1_path.name})")
        return 1

    ACTIVATIONS_DIR.mkdir(parents=True, exist_ok=True)

    # Select prompt subset
    subset = select_activation_subset()
    n_exp = sum(1 for p in subset if p["set"] == "experimental")
    n_ctrl = sum(1 for p in subset if p["set"] == "control")
    print(f"Selected {len(subset)} prompts: {n_exp} experimental + {n_ctrl} control")

    # Show fragility range
    exp_frags = [p["fragility_score"] for p in subset if p["set"] == "experimental"]
    ctrl_frags = [p["fragility_score"] for p in subset if p["set"] == "control"]
    print(f"  Experimental fragility: {min(exp_frags):.3f} - {max(exp_frags):.3f}")
    print(f"  Control fragility:      {min(ctrl_frags):.3f} - {max(ctrl_frags):.3f}")

    # Save subset metadata
    meta_path = EXP1_DIR / "activations" / "activation_subset.jsonl"
    with open(meta_path, "w", encoding="utf-8") as f:
        for p in subset:
            f.write(json.dumps({
                "prompt_idx": p["prompt_idx"],
                "category": p["category"],
                "set": p["set"],
                "fragility_score": p["fragility_score"],
                "text": p["text"],
            }, ensure_ascii=False) + "\n")

    # Load model
    from transformers import AutoModelForCausalLM, AutoTokenizer
    print(f"\nLoading {MODEL_NAME}...", flush=True)
    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token
    model = AutoModelForCausalLM.from_pretrained(
        MODEL_NAME, torch_dtype=torch.float32
    ).to(device)
    model.eval()

    n_layers = model.config.num_hidden_layers
    hidden_size = model.config.hidden_size
    print(f"Architecture: {n_layers} layers, hidden_size={hidden_size}", flush=True)

    # Test forward + deterministic mode
    try:
        with torch.no_grad():
            model(torch.tensor([[0, 1, 2]], device=device))
    except RuntimeError as e:
        if "deterministic" in str(e).lower():
            torch.use_deterministic_algorithms(True, warn_only=True)

    # Load Run 1 token sequences
    print("Loading Experiment 1A Run 1 results...", flush=True)
    run1 = load_run1_results()

    # Record activations
    print(f"\nRecording activations ({n_layers} layers + embed + norm)...", flush=True)
    t_start = time.time()

    # Accumulators for summary stats
    layer_norms_exp = {k: [] for k in ["embed"] + [f"layer_{i}" for i in range(n_layers)] + ["final_norm"]}
    layer_norms_ctrl = {k: [] for k in layer_norms_exp}

    for i, p in enumerate(subset):
        idx = p["prompt_idx"]
        r1 = run1[idx]
        full_ids = r1["input_ids"] + r1["generated_ids"]

        acts = record_activations(model, full_ids, device)

        # Save to disk
        save_data = {
            "prompt_idx": idx,
            "input_length": len(r1["input_ids"]),
            "total_length": len(full_ids),
            "set": p["set"],
            "category": p["category"],
        }
        save_data.update(acts)
        torch.save(save_data, ACTIVATIONS_DIR / f"prompt_{idx:04d}.pt")

        # Collect per-layer L2 norms for summary
        target = layer_norms_exp if p["set"] == "experimental" else layer_norms_ctrl
        for key in acts:
            target[key].append(acts[key].norm(dim=-1).mean().item())

        if (i + 1) % 10 == 0 or (i + 1) == len(subset):
            elapsed = time.time() - t_start
            print(f"  [{i+1:3d}/{len(subset)}] {elapsed:.1f}s | "
                  f"prompt_{idx:04d} ({p['set'][:4]}, {p['category'][:6]}) | "
                  f"seq_len={len(full_ids)}",
                  flush=True)

        # Free GPU cache periodically
        if (i + 1) % 20 == 0 and device.type == "cuda":
            torch.cuda.empty_cache()

    total_time = time.time() - t_start
    total_bytes = sum(f.stat().st_size for f in ACTIVATIONS_DIR.glob("*.pt"))
    print(f"\nComplete in {total_time:.1f}s | {total_bytes/1e9:.2f} GB on disk", flush=True)

    # Summary: per-layer activation norms
    print(f"\n{'='*70}")
    print("PER-LAYER ACTIVATION NORMS (mean L2 across positions)")
    print(f"{'='*70}")
    layer_keys = ["embed"] + [f"layer_{i}" for i in range(n_layers)] + ["final_norm"]
    print(f"  {'Layer':12s} {'Experimental':>14s} {'Control':>14s} {'Ratio':>8s}")
    print(f"  {'-'*52}")
    for key in layer_keys:
        e = np.mean(layer_norms_exp[key]) if layer_norms_exp[key] else 0
        c = np.mean(layer_norms_ctrl[key]) if layer_norms_ctrl[key] else 0
        ratio = e / c if c > 0 else float("inf")
        print(f"  {key:12s} {e:14.2f} {c:14.2f} {ratio:8.3f}")

    print(f"\n{'='*70}")
    print("EXPERIMENT 1B COMPLETE")
    print(f"{'='*70}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
