"""
Smoke test: validates Experiment 0 pipeline before full run.
Tests model loading, greedy generation, logit recording, determinism.
Expected runtime: ~2 min (excluding model download).
"""
import os
import sys
import time
import torch
import numpy as np
import random

SEED = 42
MODEL_NAME = "EleutherAI/pythia-410m-deduped"
MAX_NEW_TOKENS = 256


def setup_determinism():
    """Set all determinism flags. Returns (strict_ok, warn_mode)."""
    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(SEED)
        os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    # Try strict first, fall back to warn_only
    try:
        torch.use_deterministic_algorithms(True)
        return "strict"
    except Exception:
        pass
    try:
        torch.use_deterministic_algorithms(True, warn_only=True)
        return "warn_only"
    except Exception:
        return "off"


def generate_one(model, tokenizer, text, device):
    """Greedy-generate from a prompt; return token ids + top-5 logits per step."""
    inputs = tokenizer(text, return_tensors="pt").to(device)
    with torch.no_grad():
        out = model.generate(
            **inputs,
            max_new_tokens=MAX_NEW_TOKENS,
            do_sample=False,
            output_scores=True,
            return_dict_in_generate=True,
        )
    ids = out.sequences[0].cpu().tolist()
    top5v, top5i = [], []
    for s in out.scores:
        v, i = torch.topk(s[0].float().cpu(), k=5)
        top5v.append(v.tolist())
        top5i.append(i.tolist())
    return ids, top5v, top5i


def main():
    sep = "=" * 70
    print(sep)
    print("SMOKE TEST -- Experiment 0 Pipeline Validation")
    print(sep)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    gpu_name = torch.cuda.get_device_name(0) if device.type == "cuda" else "N/A"
    print(f"Device: {device} ({gpu_name})")

    det_mode = setup_determinism()
    print(f"Deterministic mode: {det_mode}")

    # ---- Load model ----
    from transformers import AutoModelForCausalLM, AutoTokenizer

    print(f"\nLoading {MODEL_NAME} (this downloads ~800 MB on first run)...")
    t0 = time.time()
    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token
    model = AutoModelForCausalLM.from_pretrained(
        MODEL_NAME, torch_dtype=torch.float32
    ).to(device)
    model.eval()
    t_load = time.time() - t0

    n_params = sum(p.numel() for p in model.parameters())
    n_layers = model.config.num_hidden_layers
    vocab = model.config.vocab_size
    mem_gb = torch.cuda.memory_allocated() / 1e9 if device.type == "cuda" else 0
    print(f"Loaded in {t_load:.1f}s | {n_params:,} params | {n_layers} layers | vocab {vocab}")
    print(f"VRAM used: {mem_gb:.2f} GB")

    # Architecture sanity
    if n_layers != 24:
        print(f"WARNING: expected 24 layers, got {n_layers}")

    # ---- Test prompts (one per category) ----
    test_prompts = [
        ("factual",       "What is the capital of France?"),
        ("open_ended",    "Write a paragraph about the feeling of rain on a tin roof."),
        ("reasoning",     "What is 247 + 389?"),
        ("continuation",  "The old house at the end of the road had been abandoned for years. "
                          "No one could remember who had lived there. One autumn evening,"),
        ("conversational","Hey, how's it going? I was wondering if you could help me with something."),
    ]

    # ---- Run 1 ----
    print(f"\n{sep}\nRUN 1 -- Baseline generation\n{sep}")
    run1_results = []
    t_run1 = time.time()
    deterministic_gen_ok = True
    for cat, text in test_prompts:
        try:
            ids, v, i = generate_one(model, tokenizer, text, device)
        except RuntimeError as e:
            if "deterministic" in str(e).lower():
                print(f"  Deterministic op failed: {e}")
                print("  Retrying with warn_only mode...")
                torch.use_deterministic_algorithms(True, warn_only=True)
                det_mode = "warn_only"
                deterministic_gen_ok = False
                ids, v, i = generate_one(model, tokenizer, text, device)
            else:
                raise

        margins = [row[0] - row[1] for row in v]
        min_m = min(margins) if margins else float("inf")
        frag = sum(1 for m in margins if m < 1.0)
        n_gen = len(v)
        input_len = len(tokenizer.encode(text))
        gen_text = tokenizer.decode(ids[input_len:], skip_special_tokens=True)
        run1_results.append({"ids": ids, "top5v": v, "n_gen": n_gen})
        print(f"  [{cat}] {n_gen} tokens | min_margin={min_m:.3f} | fragile={frag}/{n_gen}")
        print(f"    In:  {text[:70]}{'...' if len(text) > 70 else ''}")
        print(f"    Out: {gen_text[:100]}{'...' if len(gen_text) > 100 else ''}")
    t_run1 = time.time() - t_run1

    # ---- Run 2 (determinism check) ----
    print(f"\n{sep}\nRUN 2 -- Determinism verification (re-seeded)\n{sep}")
    setup_determinism()
    if det_mode == "warn_only":
        torch.use_deterministic_algorithms(True, warn_only=True)
    run2_results = []
    t_run2 = time.time()
    for cat, text in test_prompts:
        ids, v, i = generate_one(model, tokenizer, text, device)
        run2_results.append({"ids": ids, "top5v": v})
    t_run2 = time.time() - t_run2

    # ---- Compare ----
    print(f"\n{sep}\nDETERMINISM RESULTS\n{sep}")
    tok_pass = True
    logit_pass = True
    for idx in range(len(test_prompts)):
        r1, r2 = run1_results[idx], run2_results[idx]
        cat = test_prompts[idx][0]
        ids_match = r1["ids"] == r2["ids"]
        max_diff = 0.0
        for row1, row2 in zip(r1["top5v"], r2["top5v"]):
            for a, b in zip(row1, row2):
                max_diff = max(max_diff, abs(a - b))

        if not ids_match:
            tok_pass = False
            div_pt = next(
                j for j in range(min(len(r1["ids"]), len(r2["ids"])))
                if r1["ids"][j] != r2["ids"][j]
            )
            print(f"  [{cat}] TOKENS DIVERGED at position {div_pt}")
        else:
            print(f"  [{cat}] tokens=MATCH  max_logit_diff={max_diff:.2e}")
        if max_diff > 1e-5:
            logit_pass = False

    # ---- Summary ----
    per_prompt = (t_run1 + t_run2) / (2 * len(test_prompts))
    est_hours = per_prompt * 2000 / 3600
    peak_gb = torch.cuda.max_memory_allocated() / 1e9 if device.type == "cuda" else 0

    print(f"\n{sep}\nSUMMARY\n{sep}")
    print(f"Token determinism:   {'PASS' if tok_pass else 'FAIL'}")
    print(f"Logit determinism:   {'PASS' if logit_pass else 'WARN (>1e-5 diff)'}")
    print(f"Deterministic mode:  {det_mode}")
    print(f"Time per prompt:     {per_prompt:.2f}s  (gen only: {t_run1/len(test_prompts):.2f}s)")
    print(f"Est. 2000 prompts:   {est_hours:.1f} hours")
    print(f"Peak VRAM:           {peak_gb:.2f} GB / {torch.cuda.get_device_properties(0).total_memory/1e9:.1f} GB")
    if not deterministic_gen_ok:
        print(f"NOTE: Fell back to warn_only after deterministic op error during generation")

    overall = tok_pass
    print(f"\n>>> {'SMOKE TEST PASSED -- pipeline validated, ready for full experiment'          if overall else 'SMOKE TEST FAILED -- investigate before proceeding'} <<<")
    return 0 if overall else 1


if __name__ == "__main__":
    sys.exit(main())
